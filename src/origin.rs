use std::cmp;
use std::ops;

use crate::basis::Basis;
use crate::segment::OriginSegment;

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
        while let Some(p) = origin_segment.next() {
            origin_primes.push(p);
            if p * p < origin_end {
                origin_segment.strike_prime(p);
            }
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
