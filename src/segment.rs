use arr_macro::arr;
use std::cmp;
use std::mem;

use crate::constants::{ceil_div, SPOKE, SPOKE_GAPS, SPOKE_SIZE, WHEEL, WHEEL_SIZE};
use crate::spoke::Spoke;

/// A Segment of the sieve consists of several spokes within the range [segment_start, segment_end).
/// Each spoke is a residue class of the multiplicative group of integers modulo n: (Z / n Z)^x,
/// where n is a primorial. For instance, for n = 2 * 3 and range [12, 30) we have:
///     # 1 mod 6   # 5 mod 6
///     spokes[0]   spokes[1]
///    [13,        [17,
///     19,         23,
///     25]         29]
///
/// Sieving only within (Z / n Z)^x lets us skip a guaranteed SPOKE_SIZE / WHEEL_SIZE fraction of
/// candidate primes. This trades off with the increased cache/memory use of storing residues, which
/// gives a practical limit of n = 2 * 3 * 5 * 7 for modern CPU + RAM combinations.
///
/// When striking a prime p, we need to find the multiples of p that occur exactly in spokes,
/// otherwise we'd waste cycles missing spokes. Finding the next spoke multiple is done by
/// exploiting the fact (Z / n Z)^x is a group. For our example where n = 2 * 3 with WHEEL = [1, 5],
/// we know that the next multiples after p * 1 are p * 5, p * (6 + 1), p * (6 + 5), etc.
/// Multiplying by any other residue would map us outside the group because the multiplicand is not
/// a group element.
///
/// This iteration can be done efficiently by adding p multiplied by the gap _between_ spokes
/// (stored in the SPOKE_GAP array) to the current multiple.
///
/// Iterating like this with a variable SPOKE_GAPS increment involves a multiplication and an
/// addition to find each multiple. If we instead slice by spoke, then each iteration only requires
/// adding p * WHEEL_SIZE, or equivalently, just the index p _within_ the spoke vector. Compilers
/// can further optimize constant p-increment iteration with loop unrolling and SIMD instructions.

pub struct WheelSegment {
    next_prime: KWayMerge,
}

impl Iterator for WheelSegment {
    type Item = usize;

    #[inline]
    fn next(&mut self) -> Option<Self::Item> {
        self.next_prime.next()
    }
}

impl WheelSegment {
    pub fn new(origin_primes: &[usize], segment_start: usize, segment_end: usize) -> WheelSegment {
        let mut segment = Segment::new(segment_start, segment_end);
        WheelSegment::strike_primes(origin_primes, segment_start, &mut segment);
        let next_prime = KWayMerge::new(segment);

        WheelSegment { next_prime }
    }

    /// Strike primes for all spokes in this segment.
    fn strike_primes(origin_primes: &[usize], segment_start: usize, segment: &mut Segment) {
        let mut multiples = arr![Vec::with_capacity(origin_primes.len()); 48];
        for &p in origin_primes {
            // Start at p^2, or skip ahead to the first spoke in this segment.
            let factor = cmp::max(p, WheelSegment::first_wheel_factor(p, segment_start));
            let mut multiple = p * factor;
            for spoke_gap in Segment::wheel_iter(factor) {
                multiples[Segment::spoke(multiple)].push(multiple);
                multiple += p * spoke_gap;
            }
        }

        // Iterate through all primes in all spokes in spoke-major order for cache locality.
        for (spoke, spoke_multiples) in segment.spokes.iter_mut().zip(multiples.iter()) {
            for (&p, &multiple) in origin_primes.iter().zip(spoke_multiples.iter()) {
                spoke.strike_prime(p, multiple);
            }
        }
    }

    fn first_wheel_factor(p: usize, segment_start: usize) -> usize {
        let n = ceil_div(segment_start, p);
        (n / WHEEL_SIZE) * WHEEL_SIZE + WHEEL[Segment::spoke(n)]
    }
}

struct KWayMerge {
    winner: Node,
    losers: Vec<Node>,
    level_offset: usize,
    segment: Segment,
}

impl Iterator for KWayMerge {
    type Item = usize;

    fn next(&mut self) -> Option<Self::Item> {
        let (value, index) = self.winner.0.take()?;
        let next_value = self.segment.spokes[index].next();
        self.push(next_value, index);
        Some(value)
    }
}

impl KWayMerge {
    fn new(mut segment: Segment) -> KWayMerge {
        let level_offset = KWayMerge::level_offset(segment.spokes.len());
        let (winner, losers) = KWayMerge::initialize_tree(level_offset, &mut segment);

        KWayMerge {
            winner,
            losers,
            level_offset,
            segment,
        }
    }

    fn push(&mut self, option_value: Option<usize>, index: usize) {
        let mut current_index = self.leaf_index(index);
        let mut current_value = match option_value {
            Some(value) => Node(Some((value, index))),
            None => Node(None),
        };
        while let Some(parent_index) = KWayMerge::parent(current_index) {
            let parent_value = &mut self.losers[parent_index];
            if parent_value < &mut current_value {
                mem::swap(parent_value, &mut current_value);
            }
            current_index = parent_index;
        }
        // If winner had any existing value, it would be overwritten. self.push() should only
        // be called after taking self.winner.
        self.winner = current_value;
    }

    fn initialize_tree(level_offset: usize, segment: &mut Segment) -> (Node, Vec<Node>) {
        let winners_size = KWayMerge::left_child(level_offset);
        let mut winners = vec![Node(None); winners_size];
        for (index, spoke) in segment.spokes.iter_mut().enumerate() {
            winners[index + level_offset] = Node(match spoke.next() {
                Some(value) => Some((value, index)),
                None => None,
            });
        }

        let mut losers = vec![Node(None); level_offset];
        let mut curr_offset = level_offset;
        while let Some(next_offset) = KWayMerge::parent(curr_offset) {
            for next_index in next_offset..curr_offset {
                let left_node = winners[KWayMerge::left_child(next_index)];
                let right_node = winners[KWayMerge::right_child(next_index)];
                let (winner, loser) = if left_node < right_node {
                    (left_node, right_node)
                } else {
                    (right_node, left_node)
                };
                if winner.0.is_some() {
                    winners[next_index] = winner;
                }
                if loser.0.is_some() {
                    losers[next_index] = loser;
                }
            }
            curr_offset = next_offset;
        }

        let winner = *winners.get(0).unwrap_or(&Node(None));

        (winner, losers)
    }

    fn level_offset(k: usize) -> usize {
        let mut power_of_two = 1;
        while power_of_two < k {
            power_of_two *= 2;
        }
        power_of_two - 1
    }
    fn leaf_index(&self, index: usize) -> usize {
        index + self.level_offset
    }
    fn parent(index: usize) -> Option<usize> {
        if index == 0 {
            return None;
        }
        Some((index - 1) / 2)
    }
    fn left_child(index: usize) -> usize {
        2 * index + 1
    }
    fn right_child(index: usize) -> usize {
        2 * index + 2
    }
}

#[derive(PartialEq, Eq, Ord, Copy, Clone)]
struct Node(Option<(usize, usize)>);

impl PartialOrd for Node {
    fn partial_cmp(&self, other: &Node) -> Option<cmp::Ordering> {
        Some(match (self.0, other.0) {
            (Some(self_tuple), Some(other_tuple)) => self_tuple.cmp(&other_tuple),
            (Some(_), None) => cmp::Ordering::Less,
            (None, Some(_)) => cmp::Ordering::Greater,
            (None, None) => cmp::Ordering::Equal,
        })
    }
}

struct Segment {
    spokes: [Spoke; SPOKE_SIZE],
}

impl Segment {
    /// Create an unsieved Segment in [start, end).
    fn new(segment_start: usize, segment_end: usize) -> Segment {
        let spoke_builder = |spoke| Spoke::new(WHEEL[spoke], segment_start, segment_end);
        let mut spoke = 0;
        let spokes = arr![spoke_builder({spoke += 1; spoke - 1}); 48];

        Segment { spokes }
    }

    fn wheel_iter(factor: usize) -> impl Iterator<Item = &'static usize> {
        SPOKE_GAPS
            .iter()
            .cycle()
            .skip(Segment::spoke(factor))
            .take(SPOKE_SIZE)
    }
    fn spoke(n: usize) -> usize {
        SPOKE[n % WHEEL_SIZE]
    }
}
