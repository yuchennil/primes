use std::cmp;
use std::ops;

use crate::constants::{BASIS_PRIMES, FIRST_NON_BASIS_PRIME};

/// Iterate through BASIS_PRIMES within [start, end)
pub struct Basis {
    pub start: usize,
    pub end: usize,
    basis_primes_index_iter: ops::Range<usize>,
}

impl Iterator for Basis {
    type Item = usize;

    fn next(&mut self) -> Option<Self::Item> {
        let basis_primes_index = self.basis_primes_index_iter.next()?;
        Some(BASIS_PRIMES[basis_primes_index])
    }
}

impl Basis {
    pub fn new(start: usize, end: usize) -> Basis {
        let basis_primes_index_iter = Basis::primes_index(start)..Basis::primes_index(end);

        // Advance start to FIRST_NON_BASIS_PRIME so it will be correct when passed to Origin.
        let start = cmp::max(start, FIRST_NON_BASIS_PRIME);

        Basis {
            start,
            end,
            basis_primes_index_iter,
        }
    }

    /// Find the index of the first prime in BASIS_PRIMES greater than n.
    fn primes_index(n: usize) -> usize {
        if n < FIRST_NON_BASIS_PRIME {
            BASIS_PRIMES
                .iter()
                .position(|&p| p >= n)
                .unwrap_or_else(|| BASIS_PRIMES.len())
        } else {
            BASIS_PRIMES.len()
        }
    }
}
