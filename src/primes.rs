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
    segments: Vec<SieveSegment>,
    limit: u64,
}

impl Sieve {
    /// Run the Sieve of Eratosthenes to generate primes at or below limit
    pub fn eratosthenes(limit: u64) -> Sieve {
        Sieve {
            segments: vec![SieveSegment::eratosthenes(limit)],
            limit,
        }
    }

    /// Run a segmented Sieve of Eratosthenes to save memory while generating primes.
    ///
    /// segment_length is just above sqrt(limit) so that the origin segment suffices to sieve
    /// all remaining segments (which hence don't need to be kept in memory after we've finished
    /// sieving through them).
    pub fn segmented(limit: u64) -> Sieve {
        let segment_length = cmp::max(1, (limit as f64).sqrt().ceil() as u64);

        let mut segments = vec![SieveSegment::eratosthenes(segment_length)];
        for start in (segment_length..limit).step_by(segment_length as usize) {
            let end = cmp::min(start + segment_length, limit);
            segments.push(SieveSegment::from_origin(start, end, &segments[0]));
        }

        Sieve { segments, limit }
    }
}

impl IntoIterator for Sieve {
    type Item = u64;
    type IntoIter = std::vec::IntoIter<Self::Item>;

    fn into_iter(self) -> Self::IntoIter {
        let result = self
            .segments
            .into_iter()
            .map(|segment| segment.primes.into_iter())
            .flatten()
            .collect::<Vec<_>>();
        result.into_iter()
    }
}

#[derive(Debug)]
struct SieveSegment {
    primes: Vec<u64>,
    start: u64,
    end: u64,
}

impl SieveSegment {
    /// Sieve an origin segment [0, limit) using Eratosthenes optimized by skipping evens after 2.
    fn eratosthenes(limit: u64) -> SieveSegment {
        let mut sieve = WheelSieveSegment::create(0, limit);

        let mut p = 3;
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

        let mut sieve = WheelSieveSegment::create(start, end);
        // Don't sieve 2 if it exists.
        for &p in origin.primes.iter().skip(1) {
            sieve.strike_prime(p);
        }

        SieveSegment {
            primes: sieve.unpack_primes(),
            start,
            end,
        }
    }
}

#[derive(Debug)]
struct WheelSieveSegment {
    data: Vec<bool>,
    offset: usize,
    start: u64,
    end: u64,
}

impl WheelSieveSegment {
    // Create an unsieved WheelSieveSegment in [start, end).
    fn create(start: u64, end: u64) -> WheelSieveSegment {
        // If start is even, then set sieve_offset to point at the first odd after it: start + 1.
        // If start is odd, then adding 1 won't make a difference: u64_to_sieve(start + 1) will
        // still point at start. Likewise with end.
        let offset = WheelSieveSegment::u64_to_sieve(start + 1);
        let size = WheelSieveSegment::u64_to_sieve(end + 1) - offset;
        let mut data = vec![true; size];
        if offset == 0 {
            // Strike 1 from the sieve if it's within bounds, as sieving wouldn't remove it.
            if let Some(e) = data.get_mut(0) {
                *e = false;
            }
        }
        WheelSieveSegment {
            data,
            offset,
            start,
            end,
        }
    }

    // Strike multiples of prime in sieve.
    //
    // Optimize by starting the multiples search at the greater of
    // - p^2 (smaller multiples should already have been struck by previous primes)
    // - the first odd multiple of p after start
    //
    // Note that a step size of p in the sieve corresponds to a step of 2 * p in u64s.
    fn strike_prime(&mut self, p: u64) {
        /// Convenience functions to find u64 ceilings
        fn ceil_div(numerator: u64, divisor: u64) -> u64 {
            numerator / divisor + ((numerator % divisor) != 0) as u64
        }
        fn ceil_odd(n: u64) -> u64 {
            n + (n + 1) % 2
        }

        let p_sieve_start =
            WheelSieveSegment::u64_to_sieve(cmp::max(p * p, p * ceil_odd(ceil_div(self.start, p))))
                - self.offset;
        for multiple in (p_sieve_start..self.data.len()).step_by(p as usize) {
            self.data[multiple] = false;
        }
    }

    /// Find the next prime after p in the sieve, or None
    fn next(&self, p: u64) -> Option<u64> {
        Some(WheelSieveSegment::sieve_to_u64(
            (WheelSieveSegment::u64_to_sieve(p) + 1..self.data.len()).find(|&n| self.data[n])?,
        ))
    }

    /// Consume this WheelSieveSegment to unpack primes
    fn unpack_primes(self) -> Vec<u64> {
        let mut primes = Vec::new();
        // We've only sieved odd primes, so 2 will need to be manually prepended, if necessary.
        if self.start <= 2 && self.end > 2 {
            primes.push(2);
        }
        let offset = self.offset;
        primes.extend(self.data.into_iter().enumerate().filter_map(|(p, x)| {
            if x {
                Some(WheelSieveSegment::sieve_to_u64(p + offset))
            } else {
                None
            }
        }));
        primes
    }

    /// Convenience functions to convert between prime space (u64 numbers) and sieve space
    /// (usizes corresponding to encoded odd u64s).
    fn u64_to_sieve(prime: u64) -> usize {
        ((prime - 1) / 2) as usize
    }
    fn sieve_to_u64(sieve: usize) -> u64 {
        2 * sieve as u64 + 1
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
    fn sieve_eratosthenes() {
        assert_eq!(
            vec![0; 0],
            Sieve::eratosthenes(0).into_iter().collect::<Vec<_>>()
        );
        assert_eq!(
            vec![0; 0],
            Sieve::eratosthenes(1).into_iter().collect::<Vec<_>>()
        );
        assert_eq!(
            vec![0; 0],
            Sieve::eratosthenes(2).into_iter().collect::<Vec<_>>()
        );
        assert_eq!(
            vec![2],
            Sieve::eratosthenes(3).into_iter().collect::<Vec<_>>()
        );
        assert_eq!(
            vec![2, 3],
            Sieve::eratosthenes(4).into_iter().collect::<Vec<_>>()
        );
        assert_eq!(
            vec![2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47],
            Sieve::eratosthenes(50).into_iter().collect::<Vec<_>>()
        );
    }

    #[test]
    fn sieve_segmented() {
        assert_eq!(
            vec![0; 0],
            Sieve::segmented(0).into_iter().collect::<Vec<_>>()
        );
        assert_eq!(
            vec![0; 0],
            Sieve::segmented(1).into_iter().collect::<Vec<_>>()
        );
        assert_eq!(
            vec![0; 0],
            Sieve::segmented(2).into_iter().collect::<Vec<_>>()
        );
        assert_eq!(vec![2], Sieve::segmented(3).into_iter().collect::<Vec<_>>());
        assert_eq!(
            vec![2, 3],
            Sieve::segmented(4).into_iter().collect::<Vec<_>>()
        );
        assert_eq!(
            vec![2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47],
            Sieve::segmented(50).into_iter().collect::<Vec<_>>()
        );
    }

    #[test]
    fn sieve_segment_eratosthenes() {
        let sieve_segment = SieveSegment::eratosthenes(50);
        assert_eq!(0, sieve_segment.start);
        assert_eq!(50, sieve_segment.end);
        assert_eq!(
            vec![2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47],
            sieve_segment.primes.to_vec()
        );
    }

    #[test]
    fn sieve_segment_from_origin() {
        let origin = SieveSegment::eratosthenes(50);

        let sieve_segment = SieveSegment::from_origin(0, 50, &origin);
        assert_eq!(0, sieve_segment.start);
        assert_eq!(50, sieve_segment.end);
        assert_eq!(
            vec![2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47],
            sieve_segment.primes.to_vec()
        );

        let sieve_segment = SieveSegment::from_origin(50, 100, &origin);
        assert_eq!(50, sieve_segment.start);
        assert_eq!(100, sieve_segment.end);
        assert_eq!(
            vec![53, 59, 61, 67, 71, 73, 79, 83, 89, 97],
            sieve_segment.primes.to_vec()
        );

        let sieve_segment = SieveSegment::from_origin(53, 100, &origin);
        assert_eq!(53, sieve_segment.start);
        assert_eq!(100, sieve_segment.end);
        assert_eq!(
            vec![53, 59, 61, 67, 71, 73, 79, 83, 89, 97],
            sieve_segment.primes.to_vec()
        );

        let sieve_segment = SieveSegment::from_origin(54, 100, &origin);
        assert_eq!(54, sieve_segment.start);
        assert_eq!(100, sieve_segment.end);
        assert_eq!(
            vec![59, 61, 67, 71, 73, 79, 83, 89, 97],
            sieve_segment.primes.to_vec()
        );

        let sieve_segment = SieveSegment::from_origin(50, 97, &origin);
        assert_eq!(50, sieve_segment.start);
        assert_eq!(97, sieve_segment.end);
        assert_eq!(
            vec![53, 59, 61, 67, 71, 73, 79, 83, 89],
            sieve_segment.primes.to_vec()
        );

        let sieve_segment = SieveSegment::from_origin(50, 98, &origin);
        assert_eq!(50, sieve_segment.start);
        assert_eq!(98, sieve_segment.end);
        assert_eq!(
            vec![53, 59, 61, 67, 71, 73, 79, 83, 89, 97],
            sieve_segment.primes.to_vec()
        );
    }
}
