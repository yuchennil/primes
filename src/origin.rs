use std::cmp;
use std::ops;

use crate::basis::Basis;
use crate::constants::FIRST_NON_BASIS_PRIME;

/// Iterate through origin_primes within [start, end)
///
/// This requires first sieving through all primes below sqrt(end). We use a custom OriginSegment
/// to sieve, as the WheelSegment's many Spokes carry extra initialization and sorting overhead
/// that isn't relevant for smaller ranges of primes.
pub struct Origin {
    pub start: usize,
    pub end: usize,
    pub origin_primes: Vec<usize>,
    origin_primes_index_iter: ops::Range<usize>,
}

impl Iterator for Origin {
    type Item = usize;

    #[inline]
    fn next(&mut self) -> Option<Self::Item> {
        let origin_primes_index = self.origin_primes_index_iter.next()?;
        Some(self.origin_primes[origin_primes_index])
    }
}

impl From<Basis> for Origin {
    fn from(state: Basis) -> Origin {
        let start = state.start;
        let end = state.end;

        let origin_end = Origin::end(end);
        let origin_primes = Origin::primes(origin_end);
        let origin_primes_index_iter =
            Origin::primes_index(start, origin_end, &origin_primes)..origin_primes.len();

        // Advance start to origin_end so it will be correct when passed to Wheel.
        let start = cmp::max(start, origin_end);

        Origin {
            start,
            end,
            origin_primes,
            origin_primes_index_iter,
        }
    }
}

impl Origin {
    // origin_end is just above sqrt(end) so that the origin primes suffice to sieve
    // all remaining segments (which hence don't need to be kept in memory after we've finished
    // sieving through them).
    fn end(end: usize) -> usize {
        (end as f64).sqrt().ceil() as usize
    }

    /// Sieve an origin segment [0, origin_end) using Eratosthenes, skipping non-wheel numbers.
    pub fn primes(origin_end: usize) -> Vec<usize> {
        let mut origin_primes = Vec::new();

        let mut origin_segment = OriginSegment::new(origin_end);
        let mut n = FIRST_NON_BASIS_PRIME;
        while let Some(p) = origin_segment.find_prime(n) {
            origin_primes.push(p);
            if p * p < origin_end {
                origin_segment.strike_prime(p);
            }
            n = p + 1;
        }

        origin_primes
    }

    /// Find the index of the first prime in origin_primes greater than n.
    fn primes_index(n: usize, origin_end: usize, origin_primes: &[usize]) -> usize {
        if n < origin_end {
            origin_primes
                .iter()
                .position(|&p| p >= n)
                .unwrap_or_else(|| origin_primes.len())
        } else {
            origin_primes.len()
        }
    }
}

struct OriginSegment {
    spoke: Vec<bool>,
}

impl OriginSegment {
    /// Create an unsieved OriginSegment in [start, end).
    fn new(origin_end: usize) -> OriginSegment {
        let spoke_length = OriginSegment::n_to_spoke(origin_end);
        let spoke = vec![true; spoke_length];

        let mut origin_segment = OriginSegment { spoke };
        origin_segment.strike_prime(3);
        origin_segment.strike_prime(5);
        origin_segment.strike_prime(7);
        origin_segment
    }

    /// Strike multiples of prime in this spoke.
    ///
    /// Note that a step size of p in the spoke corresponds to a step of WHEEL_SIZE * p in u64s.
    fn strike_prime(&mut self, p: usize) {
        // This while loop is equivalent to a for loop that steps by p, except it's 30-40% more
        // efficient, according to benchmarks.
        let mut spoke_multiple = OriginSegment::n_to_spoke(p * p);
        while spoke_multiple < self.spoke.len() {
            self.spoke[spoke_multiple] = false;
            spoke_multiple += p;
        }
    }

    /// Find the next prime at or after n in the spoke.
    #[inline]
    fn find_prime(&self, n: usize) -> Option<usize> {
        for spoke_n in OriginSegment::n_to_spoke(n)..self.spoke.len() {
            if self.spoke[spoke_n] {
                let p = OriginSegment::spoke_to_n(spoke_n);
                return Some(p);
            }
        }
        None
    }

    /// Convert between number space and spoke space.
    fn n_to_spoke(n: usize) -> usize {
        n / 2
    }
    fn spoke_to_n(spoke_n: usize) -> usize {
        2 * spoke_n + 1
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
