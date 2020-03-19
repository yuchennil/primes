use arr_macro::arr;
use std::cmp;
use std::collections;

use crate::spoke::Spoke;
use crate::constants::{SPOKE, SPOKE_GAPS, SPOKE_SIZE, WHEEL, WHEEL_SIZE, ceil_div};

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

pub struct OriginSegment {
    segment: Segment,
}

impl OriginSegment {
    pub fn new(origin_end: usize) -> OriginSegment {
        let segment = Segment::new(0, origin_end);

        OriginSegment { segment }
    }

    /// Strike p for all spokes in the origin segment.
    pub fn strike_prime(&mut self, p: usize) {
        // Optimize by striking multiples from p^2. Smaller multiples should already have been
        // struck by previous primes.
        let factor = p;
        let mut multiple = p * factor;
        for spoke_gap in Segment::wheel_iter(factor) {
            self.segment.spokes[Segment::spoke(multiple)].strike_prime(p, multiple);
            multiple += p * spoke_gap;
        }
    }

    /// Find the next prime at or after n in the segment.
    #[inline]
    pub fn find_prime(&self, n: usize) -> Option<usize> {
        let mut min_p = None;
        for spoke in self.segment.spokes.iter() {
            if let Some(p) = spoke.find_prime(n) {
                min_p = match min_p {
                    Some(min_p) if min_p <= p => Some(min_p),
                    _ => Some(p),
                }
            }
        }
        min_p
    }
}

pub struct WheelSegment {
    segment: Segment,
    next_prime: collections::BinaryHeap<(cmp::Reverse<usize>, usize)>,
}

impl Iterator for WheelSegment {
    type Item = usize;

    #[inline]
    fn next(&mut self) -> Option<Self::Item> {
        let (cmp::Reverse(p), spoke_index) = self.next_prime.pop()?;
        if let Some(next_p) = self.segment.spokes[spoke_index].next() {
            self.next_prime.push((cmp::Reverse(next_p), spoke_index));
        }
        Some(p)
    }
}

impl WheelSegment {
    pub fn new(origin_primes: &[usize], segment_start: usize, segment_end: usize) -> WheelSegment {
        let mut segment = Segment::new(segment_start, segment_end);
        WheelSegment::strike_primes(origin_primes, segment_start, &mut segment);
        let next_prime = WheelSegment::initialize_next_prime(&mut segment);

        WheelSegment {
            segment,
            next_prime,
        }
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

    // Populate the next_primes heap with the first primes of each spoke
    fn initialize_next_prime(
        segment: &mut Segment,
    ) -> collections::BinaryHeap<(cmp::Reverse<usize>, usize)> {
        let mut next_prime = collections::BinaryHeap::new();
        for (spoke_index, spoke) in segment.spokes.iter_mut().enumerate() {
            if let Some(next_spoke_prime) = spoke.next() {
                next_prime.push((cmp::Reverse(next_spoke_prime), spoke_index));
            }
        }
        next_prime
    }

    fn first_wheel_factor(p: usize, segment_start: usize) -> usize {
        let n = ceil_div(segment_start, p);
        (n / WHEEL_SIZE) * WHEEL_SIZE + WHEEL[Segment::spoke(n)]
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

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn origin_segment_correct() {
        let mut origin_segment = OriginSegment::new(240);

        // find_prime() when all true
        assert_eq!(Some(143), origin_segment.find_prime(142));
        assert_eq!(Some(143), origin_segment.find_prime(143));
        assert_eq!(Some(149), origin_segment.find_prime(144));
        assert_eq!(Some(209), origin_segment.find_prime(208));
        assert_eq!(Some(209), origin_segment.find_prime(209));
        assert_eq!(Some(211), origin_segment.find_prime(210));
        assert_eq!(Some(211), origin_segment.find_prime(211));
        assert_eq!(Some(221), origin_segment.find_prime(212));
        assert_eq!(Some(239), origin_segment.find_prime(238));
        assert_eq!(Some(239), origin_segment.find_prime(239));
        assert_eq!(None, origin_segment.find_prime(240));
        assert_eq!(None, origin_segment.find_prime(241));

        // strike_prime()
        origin_segment.strike_prime(11);
        origin_segment.strike_prime(13);

        // find_prime() after sieving
        assert_eq!(Some(149), origin_segment.find_prime(142));
        assert_eq!(Some(149), origin_segment.find_prime(143));
        assert_eq!(Some(149), origin_segment.find_prime(144));
        assert_eq!(Some(211), origin_segment.find_prime(208));
        assert_eq!(Some(211), origin_segment.find_prime(209));
        assert_eq!(Some(211), origin_segment.find_prime(210));
        assert_eq!(Some(211), origin_segment.find_prime(211));
        assert_eq!(Some(223), origin_segment.find_prime(212));
        assert_eq!(Some(239), origin_segment.find_prime(238));
        assert_eq!(Some(239), origin_segment.find_prime(239));
        assert_eq!(None, origin_segment.find_prime(240));
        assert_eq!(None, origin_segment.find_prime(241));
    }
}