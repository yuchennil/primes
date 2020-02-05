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

#[derive(Debug)]
pub struct Sieve {
    origin: SieveSegment,
    current_iter: std::vec::IntoIter<u64>,
    current_end: u64,
    segment_length: u64,
    limit: u64,
}

impl Sieve {
    /// Run a segmented Sieve of Eratosthenes to save memory while generating primes.
    ///
    /// segment_length is just above sqrt(limit) so that the origin segment suffices to sieve
    /// all remaining segments (which hence don't need to be kept in memory after we've finished
    /// sieving through them).
    pub fn segmented(limit: u64) -> Sieve {
        let segment_length = (limit as f64).sqrt().ceil() as u64;
        let origin = SieveSegment::eratosthenes(segment_length);
        // TODO avoid cloning the origin segment
        let current_iter = origin.clone().into_iter();
        let current_end = segment_length;

        Sieve {
            origin,
            current_iter,
            current_end,
            segment_length,
            limit,
        }
    }
}

impl Iterator for Sieve {
    type Item = u64;

    fn next(&mut self) -> Option<Self::Item> {
        loop {
            if let Some(item) = self.current_iter.next() {
                return Some(item);
            }
            let start = self.current_end;
            if start == self.limit {
                return None;
            }
            self.current_end = cmp::min(start + self.segment_length, self.limit);
            self.current_iter =
                SieveSegment::from_origin(start, self.current_end, &self.origin).into_iter();
        }
    }
}

#[derive(Clone, Debug)]
struct SieveSegment {
    primes: Vec<u64>,
    start: u64,
    end: u64,
}

impl SieveSegment {
    /// Sieve an origin segment [0, limit) using Eratosthenes, skipping non-wheel numbers.
    fn eratosthenes(limit: u64) -> SieveSegment {
        let mut sieve = WheelSieveSegment::create(0, limit);

        let mut p = FIRST_NON_BASIS_PRIME;
        while p * p <= limit {
            sieve.strike_prime(p);
            if let Some(next_p) = sieve.next(p) {
                p = next_p;
            } else {
                break;
            }
        }

        SieveSegment {
            primes: sieve.unpack_primes(),
            start: 0,
            end: limit,
        }
    }

    /// Sieve a segment [start, end) based on an origin segment.
    fn from_origin(start: u64, end: u64, origin: &SieveSegment) -> SieveSegment {
        assert_eq!(0, origin.start);
        assert!(origin.end.pow(2) >= end);
        assert!(start <= end, "start greater than end");

        let mut sieve = WheelSieveSegment::create(start, end);
        // Don't sieve BASIS_PRIMES if they exist.
        for &p in origin.primes.iter().skip(BASIS_PRIMES.len()) {
            sieve.strike_prime(p);
        }

        SieveSegment {
            primes: sieve.unpack_primes(),
            start,
            end,
        }
    }
}

impl IntoIterator for SieveSegment {
    type Item = u64;
    type IntoIter = std::vec::IntoIter<Self::Item>;

    fn into_iter(self) -> Self::IntoIter {
        self.primes.into_iter()
    }
}

#[derive(Debug)]
struct WheelSieveSegment {
    data: Vec<bool>,
    start: usize,
    end: usize,
}

const BASIS_PRIMES: [u64; 4] = [2, 3, 5, 7];
const FIRST_NON_BASIS_PRIME: u64 = 11;
const WHEEL_SIZE: usize = 210;
const WHEEL_SPOKE_DIFFS: [usize; 210] = [
    0, 10, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 0, 4, 0, 0, 0, 2, 0, 4, 0, 0, 0, 6, 0, 0, 0, 0, 0, 2, 0,
    6, 0, 0, 0, 0, 0, 4, 0, 0, 0, 2, 0, 4, 0, 0, 0, 6, 0, 0, 0, 0, 0, 6, 0, 0, 0, 0, 0, 2, 0, 6, 0,
    0, 0, 0, 0, 4, 0, 0, 0, 2, 0, 6, 0, 0, 0, 0, 0, 4, 0, 0, 0, 6, 0, 0, 0, 0, 0, 8, 0, 0, 0, 0, 0,
    0, 0, 4, 0, 0, 0, 2, 0, 4, 0, 0, 0, 2, 0, 4, 0, 0, 0, 8, 0, 0, 0, 0, 0, 0, 0, 6, 0, 0, 0, 0, 0,
    4, 0, 0, 0, 6, 0, 0, 0, 0, 0, 2, 0, 4, 0, 0, 0, 6, 0, 0, 0, 0, 0, 2, 0, 6, 0, 0, 0, 0, 0, 6, 0,
    0, 0, 0, 0, 4, 0, 0, 0, 2, 0, 4, 0, 0, 0, 6, 0, 0, 0, 0, 0, 2, 0, 6, 0, 0, 0, 0, 0, 4, 0, 0, 0,
    2, 0, 4, 0, 0, 0, 2, 0, 10, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2,
];
const CEILING_SPOKE: [usize; 210] = [
    1, 1, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 13, 13, 17, 17, 17, 17, 19, 19, 23, 23, 23, 23,
    29, 29, 29, 29, 29, 29, 31, 31, 37, 37, 37, 37, 37, 37, 41, 41, 41, 41, 43, 43, 47, 47, 47, 47,
    53, 53, 53, 53, 53, 53, 59, 59, 59, 59, 59, 59, 61, 61, 67, 67, 67, 67, 67, 67, 71, 71, 71, 71,
    73, 73, 79, 79, 79, 79, 79, 79, 83, 83, 83, 83, 89, 89, 89, 89, 89, 89, 97, 97, 97, 97, 97, 97,
    97, 97, 101, 101, 101, 101, 103, 103, 107, 107, 107, 107, 109, 109, 113, 113, 113, 113, 121,
    121, 121, 121, 121, 121, 121, 121, 127, 127, 127, 127, 127, 127, 131, 131, 131, 131, 137, 137,
    137, 137, 137, 137, 139, 139, 143, 143, 143, 143, 149, 149, 149, 149, 149, 149, 151, 151, 157,
    157, 157, 157, 157, 157, 163, 163, 163, 163, 163, 163, 167, 167, 167, 167, 169, 169, 173, 173,
    173, 173, 179, 179, 179, 179, 179, 179, 181, 181, 187, 187, 187, 187, 187, 187, 191, 191, 191,
    191, 193, 193, 197, 197, 197, 197, 199, 199, 209, 209, 209, 209, 209, 209, 209, 209, 209, 209,
];

impl WheelSieveSegment {
    // Create an unsieved WheelSieveSegment in [start, end).
    fn create(start: u64, end: u64) -> WheelSieveSegment {
        let start = start as usize;
        let end = end as usize;

        let mut data = vec![false; end - start];

        let mut spoke = WheelSieveSegment::first_spoke(start);
        while spoke < end {
            data[spoke - start] = true;
            spoke = WheelSieveSegment::next_spoke(spoke);
        }

        if start < 2 {
            // Strike 1 from the sieve if it's within bounds, as sieving wouldn't remove it.
            if let Some(e) = data.get_mut(1 - start) {
                *e = false;
            }
        }
        WheelSieveSegment { data, start, end }
    }

    // Strike multiples of prime in sieve.
    //
    // Optimize by starting the multiples search at the greater of
    // - p^2 (smaller multiples should already have been struck by previous primes)
    // - the first wheel multiple of p after start
    //
    // Note that a step size of p in the sieve corresponds to a step of 2 * p in u64s.
    fn strike_prime(&mut self, p: u64) {
        let p = p as usize;
        fn ceil_div(a: usize, b: usize) -> usize {
            a / b + (a % b != 0) as usize
        }
        let mut factor = cmp::max(p, WheelSieveSegment::first_spoke(ceil_div(self.start, p)));
        while p * factor < self.end {
            self.data[p * factor - self.start] = false;
            factor = WheelSieveSegment::next_spoke(factor);
        }
    }

    /// Find the next prime after p in the sieve, or None
    fn next(&self, p: u64) -> Option<u64> {
        let mut p = p as usize;
        p = WheelSieveSegment::next_spoke(p);
        while p < self.end {
            if self.data[p - self.start] {
                return Some(p as u64);
            }
            p = WheelSieveSegment::next_spoke(p);
        }
        None
    }

    /// Consume this WheelSieveSegment to unpack primes
    fn unpack_primes(self) -> Vec<u64> {
        let mut primes = Vec::new();
        // We've only sieved primes after BASIS_PRIMES, so these will need to be manually prepended,
        // if necessary.
        primes.extend(
            BASIS_PRIMES
                .iter()
                .filter(|&&p| p >= self.start as u64 && p < self.end as u64),
        );
        let start = self.start;
        primes.extend(
            self.data
                .into_iter()
                .enumerate()
                .filter_map(|(segment_p, x)| {
                    if x {
                        Some((segment_p + start) as u64)
                    } else {
                        None
                    }
                }),
        );
        primes
    }

    fn first_spoke(start: usize) -> usize {
        (start / WHEEL_SIZE) * WHEEL_SIZE + CEILING_SPOKE[start % WHEEL_SIZE]
    }
    fn next_spoke(p: usize) -> usize {
        p + WHEEL_SPOKE_DIFFS[p % WHEEL_SIZE]
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

    #[test]
    fn sieve_segment_eratosthenes() {
        let sieve_segment = SieveSegment::eratosthenes(0);
        assert_eq!(0, sieve_segment.start);
        assert_eq!(0, sieve_segment.end);
        assert_eq!(vec![0; 0], sieve_segment.primes.to_vec());

        let sieve_segment = SieveSegment::eratosthenes(1);
        assert_eq!(0, sieve_segment.start);
        assert_eq!(1, sieve_segment.end);
        assert_eq!(vec![0; 0], sieve_segment.primes.to_vec());

        let sieve_segment = SieveSegment::eratosthenes(2);
        assert_eq!(0, sieve_segment.start);
        assert_eq!(2, sieve_segment.end);
        assert_eq!(vec![0; 0], sieve_segment.primes.to_vec());

        let sieve_segment = SieveSegment::eratosthenes(3);
        assert_eq!(0, sieve_segment.start);
        assert_eq!(3, sieve_segment.end);
        assert_eq!(vec![2], sieve_segment.primes.to_vec());

        let sieve_segment = SieveSegment::eratosthenes(4);
        assert_eq!(0, sieve_segment.start);
        assert_eq!(4, sieve_segment.end);
        assert_eq!(vec![2, 3], sieve_segment.primes.to_vec());

        let sieve_segment = SieveSegment::eratosthenes(1000);
        assert_eq!(0, sieve_segment.start);
        assert_eq!(1000, sieve_segment.end);
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
            sieve_segment.primes.to_vec()
        );
    }

    #[test]
    fn sieve_segment_from_origin() {
        let origin = SieveSegment::eratosthenes(200);

        let sieve_segment = SieveSegment::from_origin(0, 200, &origin);
        assert_eq!(0, sieve_segment.start);
        assert_eq!(200, sieve_segment.end);
        assert_eq!(
            vec![
                2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79,
                83, 89, 97, 101, 103, 107, 109, 113, 127, 131, 137, 139, 149, 151, 157, 163, 167,
                173, 179, 181, 191, 193, 197, 199,
            ],
            sieve_segment.primes.to_vec()
        );

        let sieve_segment = SieveSegment::from_origin(200, 400, &origin);
        assert_eq!(200, sieve_segment.start);
        assert_eq!(400, sieve_segment.end);
        assert_eq!(
            vec![
                211, 223, 227, 229, 233, 239, 241, 251, 257, 263, 269, 271, 277, 281, 283, 293,
                307, 311, 313, 317, 331, 337, 347, 349, 353, 359, 367, 373, 379, 383, 389, 397
            ],
            sieve_segment.primes.to_vec()
        );

        let sieve_segment = SieveSegment::from_origin(211, 400, &origin);
        assert_eq!(211, sieve_segment.start);
        assert_eq!(400, sieve_segment.end);
        assert_eq!(
            vec![
                211, 223, 227, 229, 233, 239, 241, 251, 257, 263, 269, 271, 277, 281, 283, 293,
                307, 311, 313, 317, 331, 337, 347, 349, 353, 359, 367, 373, 379, 383, 389, 397
            ],
            sieve_segment.primes.to_vec()
        );

        let sieve_segment = SieveSegment::from_origin(212, 400, &origin);
        assert_eq!(212, sieve_segment.start);
        assert_eq!(400, sieve_segment.end);
        assert_eq!(
            vec![
                223, 227, 229, 233, 239, 241, 251, 257, 263, 269, 271, 277, 281, 283, 293, 307,
                311, 313, 317, 331, 337, 347, 349, 353, 359, 367, 373, 379, 383, 389, 397
            ],
            sieve_segment.primes.to_vec()
        );

        let sieve_segment = SieveSegment::from_origin(200, 397, &origin);
        assert_eq!(200, sieve_segment.start);
        assert_eq!(397, sieve_segment.end);
        assert_eq!(
            vec![
                211, 223, 227, 229, 233, 239, 241, 251, 257, 263, 269, 271, 277, 281, 283, 293,
                307, 311, 313, 317, 331, 337, 347, 349, 353, 359, 367, 373, 379, 383, 389
            ],
            sieve_segment.primes.to_vec()
        );

        let sieve_segment = SieveSegment::from_origin(200, 398, &origin);
        assert_eq!(200, sieve_segment.start);
        assert_eq!(398, sieve_segment.end);
        assert_eq!(
            vec![
                211, 223, 227, 229, 233, 239, 241, 251, 257, 263, 269, 271, 277, 281, 283, 293,
                307, 311, 313, 317, 331, 337, 347, 349, 353, 359, 367, 373, 379, 383, 389, 397
            ],
            sieve_segment.primes.to_vec()
        );
    }
}
