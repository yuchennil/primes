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

        let sieve = vec![true; Sieve::n_to_sieve(segment_end)];
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
        let sieve_segment_end = Sieve::n_to_sieve(self.segment_end);
        for sieve_n in Sieve::n_to_sieve(3)..sieve_segment_end {
            if self.sieve[sieve_n] {
                let p = Sieve::sieve_to_n(sieve_n);
                self.origin_primes.push(p);
                self.multiples.push(p * p);
                for sieve_multiple in (Sieve::n_to_sieve(p * p)..sieve_segment_end).step_by(p) {
                    self.sieve[sieve_multiple] = false;
                }
            }
        }
    }

    fn sieve_segment(&mut self) {
        let sieve_segment_start = Sieve::n_to_sieve(self.segment_start);
        let sieve_segment_end = Sieve::n_to_sieve(self.segment_end);
        let sieve_segment_length = sieve_segment_end - sieve_segment_start;
        self.sieve = vec![true; sieve_segment_length];
        for (&p, multiple) in self.origin_primes.iter().zip(self.multiples.iter_mut()) {
            let mut sieve_segment_multiple = Sieve::n_to_sieve(*multiple) - sieve_segment_start;
            while sieve_segment_multiple < sieve_segment_length {
                self.sieve[sieve_segment_multiple] = false;
                sieve_segment_multiple += p;
            }
            *multiple = Sieve::sieve_to_n(sieve_segment_multiple + sieve_segment_start);
        }
    }

    fn n_to_sieve(n: usize) -> usize {
        n / 2
    }
    fn sieve_to_n(sieve_n: usize) -> usize {
        2 * sieve_n + 1
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
            let sieve_segment_start = Sieve::n_to_sieve(self.segment_start);
            let sieve_segment_length = Sieve::n_to_sieve(self.segment_end) - sieve_segment_start;
            let mut sieve_segment_n = Sieve::n_to_sieve(self.n) - sieve_segment_start;
            while sieve_segment_n < sieve_segment_length {
                if self.sieve[sieve_segment_n] {
                    let p = Sieve::sieve_to_n(sieve_segment_n + sieve_segment_start) as u64;
                    self.n = Sieve::sieve_to_n(sieve_segment_n + sieve_segment_start + 1);
                    return Some(p);
                }
                sieve_segment_n += 1;
            }
            self.n = Sieve::sieve_to_n(sieve_segment_n + sieve_segment_start);
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
}
