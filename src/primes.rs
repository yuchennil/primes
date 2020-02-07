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
    sieve: Vec<bool>,
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

        let sieve = vec![true; segment_length];
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

    fn sieve_origin(&mut self) {
        for p in (3..self.segment_length).step_by(2) {
            if self.sieve[p] {
                self.origin_primes.push(p);
                self.multiples.push(p * p);
                for multiple in (p * p..self.segment_length).step_by(2 * p) {
                    self.sieve[multiple] = false;
                }
            }
        }
    }

    fn sieve_segment(&mut self) {
        self.sieve = vec![true; self.segment_length];
        for (&prime, multiple) in self.origin_primes.iter().zip(self.multiples.iter_mut()) {
            let mut segment_multiple = *multiple - self.segment_start;
            while segment_multiple < self.segment_length {
                self.sieve[segment_multiple] = false;
                segment_multiple += 2 * prime;
            }
            *multiple = segment_multiple + self.segment_start;
        }
    }
}

impl Iterator for Sieve {
    type Item = u64;

    fn next(&mut self) -> Option<Self::Item> {
        if self.n == 2 && self.limit > 2 {
            let result = self.n as u64;
            self.n = 3;
            return Some(result);
        }
        loop {
            while self.n < self.segment_end {
                if self.sieve[self.n - self.segment_start] {
                    let result = self.n as u64;
                    self.n += 2;
                    return Some(result);
                }
                self.n += 2;
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
    fn sieve_segmented() {
        assert_eq!(vec![0; 0], Sieve::segmented(0).collect::<Vec<_>>());
        assert_eq!(vec![0; 0], Sieve::segmented(1).collect::<Vec<_>>());
        assert_eq!(vec![0; 0], Sieve::segmented(2).collect::<Vec<_>>());
        assert_eq!(vec![2], Sieve::segmented(3).collect::<Vec<_>>());
        assert_eq!(vec![2, 3], Sieve::segmented(4).collect::<Vec<_>>());
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
}
