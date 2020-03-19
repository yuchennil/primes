use arr_macro::arr;
use std::cmp;
use std::collections;

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

pub struct Segment {
    spokes: [Spoke; SPOKE_SIZE],
    next_prime: collections::BinaryHeap<(cmp::Reverse<usize>, usize)>,
}

impl Iterator for Segment {
    type Item = usize;

    #[inline]
    fn next(&mut self) -> Option<Self::Item> {
        let (cmp::Reverse(p), spoke_index) = self.next_prime.pop()?;
        if let Some(next_p) = self.spokes[spoke_index].next() {
            self.next_prime.push((cmp::Reverse(next_p), spoke_index));
        }
        Some(p)
    }
}

impl Segment {
    pub fn new(origin_primes: &[usize], segment_start: usize, segment_end: usize) -> Segment {
        let spoke_builder = |spoke| Spoke::new(WHEEL[spoke], segment_start, segment_end);
        let mut spoke = 0;
        let spokes = arr![spoke_builder({spoke += 1; spoke - 1}); 48];
        let next_prime = collections::BinaryHeap::new();

        let mut segment = Segment { spokes, next_prime };
        segment.strike_primes(origin_primes, segment_start);
        segment.initialize_next_prime();

        segment
    }

    /// Strike primes for all spokes in this segment.
    fn strike_primes(&mut self, origin_primes: &[usize], segment_start: usize) {
        let mut multiples = arr![Vec::with_capacity(origin_primes.len()); 48];
        for &p in origin_primes {
            // Start at p^2, or skip ahead to the first spoke in this segment.
            let factor = cmp::max(p, Segment::first_wheel_factor(p, segment_start));
            let mut multiple = p * factor;
            for spoke_gap in Segment::wheel_iter(factor) {
                multiples[Segment::spoke(multiple)].push(multiple);
                multiple += p * spoke_gap;
            }
        }

        // Iterate through all primes in all spokes in spoke-major order for cache locality.
        for (spoke, spoke_multiples) in self.spokes.iter_mut().zip(multiples.iter()) {
            for (&p, &multiple) in origin_primes.iter().zip(spoke_multiples.iter()) {
                spoke.strike_prime(p, multiple);
            }
        }
    }

    // Populate the next_primes heap with the first primes of each spoke
    fn initialize_next_prime(&mut self) {
        for (spoke_index, spoke) in self.spokes.iter_mut().enumerate() {
            if let Some(next_spoke_prime) = spoke.next() {
                self.next_prime
                    .push((cmp::Reverse(next_spoke_prime), spoke_index));
            }
        }
    }

    fn first_wheel_factor(p: usize, segment_start: usize) -> usize {
        let n = ceil_div(segment_start, p);
        (n / WHEEL_SIZE) * WHEEL_SIZE + WHEEL[Segment::spoke(n)]
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
