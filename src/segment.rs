use std::cmp;

use crate::constants::{ceil_div, SPOKE, SPOKE_GAPS, SPOKE_SIZE, WHEEL, WHEEL_SIZE};
use crate::k_way_merge::KWayMerge;
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
    next_prime: KWayMerge<usize, Spoke>,
}

impl Iterator for WheelSegment {
    type Item = usize;

    fn next(&mut self) -> Option<Self::Item> {
        self.next_prime.next()
    }
}

impl WheelSegment {
    pub fn new(origin_primes: &[usize], segment_start: usize, segment_end: usize) -> WheelSegment {
        let mut segment = Segment::new(segment_start, segment_end);
        WheelSegment::strike_primes(origin_primes, segment_start, &mut segment);
        let next_prime = KWayMerge::new(segment.spokes);

        WheelSegment { next_prime }
    }

    /// Strike primes for all spokes in this segment.
    fn strike_primes(origin_primes: &[usize], segment_start: usize, segment: &mut Segment) {
        let mut multiples = vec![Vec::with_capacity(origin_primes.len()); SPOKE_SIZE];
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

struct Segment {
    spokes: Vec<Spoke>,
}

impl Segment {
    /// Create an unsieved Segment in [start, end).
    fn new(segment_start: usize, segment_end: usize) -> Segment {
        let spoke_builder = |spoke| Spoke::new(WHEEL[spoke], segment_start, segment_end);
        let spokes = (0..SPOKE_SIZE).map(spoke_builder).collect();

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
