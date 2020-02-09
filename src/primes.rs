/// Library for utilities related to primes
use std::cmp;

/// Euclidean algorithm
pub fn gcd(a: u64, b: u64) -> u64 {
    let mut a = a;
    let mut b = b;
    while b > 0 {
        let t = b;
        b = a % b;
        a = t;
    }
    a
}

pub struct Sieve {
    sieve: SieveSegment,
    origin_primes: Vec<usize>,
    multiples: Vec<usize>,
    limit: usize,
    segment_length: usize,
    segment_start: usize,
    segment_end: usize,
    n: usize,
}

impl Sieve {
    /// Run a segmented Sieve of Eratosthenes to save memory while generating primes.
    ///
    /// segment_length is just above sqrt(limit) so that the origin segment suffices to sieve
    /// all remaining segments (which hence don't need to be kept in memory after we've finished
    /// sieving through them).
    pub fn segmented(limit: u64) -> Sieve {
        let limit = limit as usize;
        let segment_length = (limit as f64).sqrt().ceil() as usize;
        let segment_start = 0;
        let segment_end = segment_length;
        let n = 2;

        let sieve = SieveSegment::new(segment_start, segment_end);
        let origin_primes = Vec::new();
        let multiples = Vec::new();

        let mut sieve = Sieve {
            sieve,
            origin_primes,
            multiples,
            limit,
            segment_length,
            segment_start,
            segment_end,
            n,
        };
        sieve.sieve_origin();
        sieve.sieve_segment();
        sieve
    }

    /// Sieve an origin segment [0, limit) using Eratosthenes, skipping non-wheel numbers.
    ///
    /// Optimize by starting the multiples search at p^2 (smaller multiples should already have been
    /// struck by previous primes)
    fn sieve_origin(&mut self) {
        let mut n = 3;
        while n < self.segment_end {
            match self.sieve.find_prime(n) {
                (next_n, Some(p)) => {
                    n = next_n;
                    self.origin_primes.push(p);
                    self.multiples.push(p * p);
                    self.sieve.strike_prime(p, p * p);
                }
                (_, None) => break,
            };
        }
    }

    /// Sieve a segment [start, end) based on an origin segment.
    ///
    /// Optimize by starting the multiples search at the first wheel multiple of p after start,
    /// which should already be set in multiples
    fn sieve_segment(&mut self) {
        self.sieve.reset(self.segment_start, self.segment_end);
        for (&p, multiple) in self.origin_primes.iter().zip(self.multiples.iter_mut()) {
            *multiple = self.sieve.strike_prime(p, *multiple);
        }
    }
}

impl Iterator for Sieve {
    type Item = u64;

    fn next(&mut self) -> Option<Self::Item> {
        // Wheel basis primes aren't included in the sieve output
        if self.n == 2 && self.limit > 2 {
            let result = self.n as u64;
            self.n = 3;
            return Some(result);
        }
        // Now iterate through wheel primes
        loop {
            match self.sieve.find_prime(self.n) {
                (next_n, Some(p)) => {
                    self.n = next_n;
                    return Some(p as u64);
                }
                (next_n, None) => self.n = next_n,
            }
            if self.segment_end == self.limit {
                return None;
            }
            self.segment_start += self.segment_length;
            self.segment_end = cmp::min(self.segment_start + self.segment_length, self.limit);
            self.sieve_segment();
        }
    }
}

struct SieveSegment {
    sieve: Vec<bool>,
    sieve_segment_start: usize,
    sieve_segment_length: usize,
}

impl SieveSegment {
    /// Create an unsieved SieveSegment in [start, end).
    fn new(segment_start: usize, segment_end: usize) -> SieveSegment {
        let mut sieve_segment = SieveSegment {
            sieve: Vec::new(),
            sieve_segment_start: 0,
            sieve_segment_length: 0,
        };
        sieve_segment.reset(segment_start, segment_end);
        sieve_segment
    }

    /// Reset this SieveSegment in-place to avoid constructor/deconstructor costs.
    fn reset(&mut self, segment_start: usize, segment_end: usize) {
        self.sieve_segment_start = SieveSegment::n_to_sieve(segment_start);
        self.sieve_segment_length =
            SieveSegment::n_to_sieve(segment_end) - self.sieve_segment_start;
        self.sieve = vec![true; self.sieve_segment_length];
    }

    /// Strike multiples of prime in sieve. Return the next multiple past the end of this
    /// SieveSegment.
    ///
    /// Note that a step size of p in the sieve corresponds to a step of 2 * p in u64s.
    fn strike_prime(&mut self, p: usize, multiple: usize) -> usize {
        let mut sieve_segment_multiple =
            SieveSegment::n_to_sieve(multiple) - self.sieve_segment_start;
        while sieve_segment_multiple < self.sieve_segment_length {
            self.sieve[sieve_segment_multiple] = false;
            sieve_segment_multiple += p;
        }
        SieveSegment::sieve_to_n(sieve_segment_multiple + self.sieve_segment_start)
    }

    /// Find the next prime at or after n in the sieve. Also return one n past the next prime.
    fn find_prime(&self, n: usize) -> (usize, Option<usize>) {
        let mut sieve_segment_n = SieveSegment::n_to_sieve(n) - self.sieve_segment_start;
        while sieve_segment_n < self.sieve_segment_length {
            if self.sieve[sieve_segment_n] {
                let p = SieveSegment::sieve_to_n(sieve_segment_n + self.sieve_segment_start);
                let next_n = p + 2;
                return (next_n, Some(p));
            }
            sieve_segment_n += 1;
        }
        let next_n = SieveSegment::sieve_to_n(sieve_segment_n + self.sieve_segment_start);
        (next_n, None)
    }

    /// Convert between number space and sieve space.
    ///
    /// If n = 2k + 1 is odd, then sieve_to_n(n_to_sieve(n)) = 2k + 1 is an identity, as desired.
    ///
    /// But if n = 2k is even, then sieve_to_n(n_to_sieve(n)) = 2k + 1 = n + 1 is not an identity.
    /// This error is acceptable because primes past 2 are all odd, so:
    /// - a lower bound of 2k + 1 won't skip any prime candidate.
    /// - an upper bound of 2k + 1 means 2k - 1 is the last prime candidate. 2k isn't possible.
    fn n_to_sieve(n: usize) -> usize {
        n / 2
    }
    fn sieve_to_n(sieve_n: usize) -> usize {
        2 * sieve_n + 1
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn gcd_correct() {
        assert_eq!(4, gcd(8, 12));
        assert_eq!(4, gcd(12, 8));
        assert_eq!(15, gcd(15, 15));
        assert_eq!(1, gcd(47, 23));
    }

    #[test]
    fn sieve_segmented_correct() {
        assert_eq!(vec![0; 0], Sieve::segmented(0).collect::<Vec<_>>());
        assert_eq!(vec![0; 0], Sieve::segmented(1).collect::<Vec<_>>());
        assert_eq!(vec![0; 0], Sieve::segmented(2).collect::<Vec<_>>());
        assert_eq!(vec![2], Sieve::segmented(3).collect::<Vec<_>>());
        assert_eq!(vec![2, 3], Sieve::segmented(4).collect::<Vec<_>>());
        assert_eq!(vec![2, 3, 5, 7], Sieve::segmented(9).collect::<Vec<_>>());
        assert_eq!(
            vec![
                2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79,
                83, 89, 97, 101, 103, 107, 109, 113, 127, 131, 137, 139, 149, 151, 157, 163, 167,
                173, 179, 181, 191, 193, 197, 199, 211, 223, 227, 229, 233, 239, 241, 251, 257,
                263, 269, 271, 277, 281, 283, 293, 307, 311, 313, 317, 331, 337, 347, 349, 353,
                359, 367, 373, 379, 383, 389, 397, 401, 409, 419, 421, 431, 433, 439, 443, 449,
                457, 461, 463, 467, 479, 487, 491, 499, 503, 509, 521, 523, 541, 547, 557, 563,
                569, 571, 577, 587, 593, 599, 601, 607, 613, 617, 619, 631, 641, 643, 647, 653,
                659, 661, 673, 677, 683, 691, 701, 709, 719, 727, 733, 739, 743, 751, 757, 761,
                769, 773, 787, 797, 809, 811, 821, 823, 827, 829, 839, 853, 857, 859, 863, 877,
                881, 883, 887, 907, 911, 919, 929, 937, 941, 947, 953, 967, 971, 977, 983, 991,
                997
            ],
            Sieve::segmented(1000).collect::<Vec<_>>()
        );
    }

    #[test]
    fn sieve_segment_correct() {
        let mut sieve_segment = SieveSegment::new(5, 16);

        // find_prime() when all true
        assert_eq!((7, Some(5)), sieve_segment.find_prime(5));
        assert_eq!((9, Some(7)), sieve_segment.find_prime(7));
        assert_eq!((11, Some(9)), sieve_segment.find_prime(9));
        assert_eq!((13, Some(11)), sieve_segment.find_prime(11));
        assert_eq!((15, Some(13)), sieve_segment.find_prime(13));
        assert_eq!((17, Some(15)), sieve_segment.find_prime(15));
        assert_eq!((17, None), sieve_segment.find_prime(17));

        // strike_prime()
        assert_eq!(21, sieve_segment.strike_prime(3, 9));
        assert_eq!(25, sieve_segment.strike_prime(5, 25));

        // find_prime() after sieving
        assert_eq!((7, Some(5)), sieve_segment.find_prime(5));
        assert_eq!((9, Some(7)), sieve_segment.find_prime(7));
        assert_eq!((13, Some(11)), sieve_segment.find_prime(9));
        assert_eq!((15, Some(13)), sieve_segment.find_prime(13));
        assert_eq!((17, None), sieve_segment.find_prime(15));

        // reset()
        sieve_segment.reset(16, 22);

        // find_prime() after reset
        assert_eq!((19, Some(17)), sieve_segment.find_prime(17));
        assert_eq!((21, Some(19)), sieve_segment.find_prime(19));
        assert_eq!((23, Some(21)), sieve_segment.find_prime(21));
        assert_eq!((23, None), sieve_segment.find_prime(23));
    }

    #[test]
    fn sieve_to_n_to_sieve_correct() {
        assert_eq!(0, SieveSegment::n_to_sieve(SieveSegment::sieve_to_n(0)));
        assert_eq!(1, SieveSegment::n_to_sieve(SieveSegment::sieve_to_n(1)));
        assert_eq!(2, SieveSegment::n_to_sieve(SieveSegment::sieve_to_n(2)));
        assert_eq!(3, SieveSegment::n_to_sieve(SieveSegment::sieve_to_n(3)));
        assert_eq!(4, SieveSegment::n_to_sieve(SieveSegment::sieve_to_n(4)));
        assert_eq!(5, SieveSegment::n_to_sieve(SieveSegment::sieve_to_n(5)));
        assert_eq!(6, SieveSegment::n_to_sieve(SieveSegment::sieve_to_n(6)));
        assert_eq!(7, SieveSegment::n_to_sieve(SieveSegment::sieve_to_n(7)));
        assert_eq!(8, SieveSegment::n_to_sieve(SieveSegment::sieve_to_n(8)));
        assert_eq!(9, SieveSegment::n_to_sieve(SieveSegment::sieve_to_n(9)));
    }
}
