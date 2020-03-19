use crate::bit_vec::BitVec;
use crate::constants::WHEEL_SIZE;

pub struct Spoke {
    sieve: BitVec,
    residue: usize,
    spoke_start: usize,
    spoke_length: usize,
}

impl Iterator for Spoke {
    type Item = usize;

    fn next(&mut self) -> Option<Self::Item> {
        let spoke_p = self.sieve.next()?;
        Some(self.spoke_to_n(spoke_p))
    }
}

impl Spoke {
    /// Create an unsieved Spoke in [start, end).
    pub fn new(residue: usize, segment_start: usize, segment_end: usize) -> Spoke {
        let spoke_start = Spoke::n_to_spoke_start(segment_start, residue);
        let spoke_length = Spoke::n_to_spoke_start(segment_end, residue) - spoke_start;
        let sieve = BitVec::new(spoke_length);

        Spoke {
            sieve,
            residue,
            spoke_start,
            spoke_length,
        }
    }

    /// Strike multiples of prime in this spoke.
    ///
    /// Note that a step size of p in the spoke corresponds to a step of WHEEL_SIZE * p in u64s.
    pub fn strike_prime(&mut self, p: usize, multiple: usize) {
        // This while loop is equivalent to a for loop that steps by p, except it's 30-40% more
        // efficient, according to benchmarks.
        let mut spoke_multiple = self.n_to_spoke(multiple);
        while spoke_multiple < self.spoke_length {
            self.sieve.unset(spoke_multiple);
            spoke_multiple += p;
        }
    }

    /// Skip the spoke iterator forward to n.
    pub fn skip_to(&mut self, n: usize) {
        let spoke_n = self.n_to_spoke(n);
        self.sieve.skip_to(spoke_n);
    }

    /// Convert between number space and spoke space.
    fn n_to_spoke_start(n: usize, residue: usize) -> usize {
        n / WHEEL_SIZE + (n % WHEEL_SIZE > residue) as usize
    }
    fn spoke_start_to_n(spoke_start_n: usize, residue: usize) -> usize {
        spoke_start_n * WHEEL_SIZE + residue
    }
    fn n_to_spoke(&self, n: usize) -> usize {
        Spoke::n_to_spoke_start(n, self.residue) - self.spoke_start
    }
    fn spoke_to_n(&self, spoke_n: usize) -> usize {
        Spoke::spoke_start_to_n(spoke_n + self.spoke_start, self.residue)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn spoke_to_n_to_spoke_correct() {
        let spoke = Spoke::new(13, 20, 50);
        assert_eq!(5, spoke.n_to_spoke(spoke.spoke_to_n(5)));
        assert_eq!(6, spoke.n_to_spoke(spoke.spoke_to_n(6)));
        assert_eq!(7, spoke.n_to_spoke(spoke.spoke_to_n(7)));
        assert_eq!(8, spoke.n_to_spoke(spoke.spoke_to_n(8)));
        assert_eq!(9, spoke.n_to_spoke(spoke.spoke_to_n(9)));
    }
}