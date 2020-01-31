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
    /// The segment length is just above sqrt(limit) so that the origin segment suffices to sieve
    /// all remaining segments (which hence don't need to be kept in memory after we've finished
    /// sieving through them).
    pub fn segmented(limit: u64) -> Sieve {
        let segment_length = cmp::max(1, (limit as f64).sqrt().ceil() as u64);
        let origin = SieveSegment::eratosthenes(segment_length);

        let mut segments = Vec::new();
        for start in (segment_length..limit).step_by(segment_length as usize) {
            let end = cmp::min(start + segment_length, limit);
            segments.push(SieveSegment::from_origin(start, end, &origin));
        }

        Sieve {
            segments: vec![origin].into_iter().chain(segments).collect(),
            limit,
        }
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
        if limit <= 2 {
            return SieveSegment {
                primes: Vec::new(),
                start: 0,
                end: limit,
            };
        }

        let mut sieve = vec![true; to_sieve(limit + 1)];
        // to_sieve(1) == 0
        sieve[0] = false;

        let mut p = 3;
        while p * p <= limit {
            // Optimize by starting the multiples search at p^2, p^2 + 2p, ...
            // instead of 2p, 3p, ...
            //
            // Note that a step size of p in the sieve corresponds to a step of 2 * p in u64s.
            for multiple in (to_sieve(p * p)..sieve.len()).step_by(p as usize) {
                sieve[multiple] = false;
            }
            match (to_sieve(p) + 1..sieve.len()).find(|&n| sieve[n]) {
                Some(next) => p = to_prime(next),
                None => break,
            }
        }

        SieveSegment {
            primes: vec![2]
                .into_iter()
                .chain(sieve.iter().enumerate().filter_map(|(p, &x)| {
                    if x {
                        Some(to_prime(p))
                    } else {
                        None
                    }
                }))
                .collect(),
            start: 0,
            end: limit,
        }
    }

    /// Sieve a segment [start, end) based on an origin segment.
    fn from_origin(start: u64, end: u64, origin: &SieveSegment) -> SieveSegment {
        assert_eq!(0, origin.start);
        assert!(origin.end.pow(2) >= end);

        // If start is even, then set sieve_start to point at the first odd after it: start + 1.
        // If start is odd, then adding 1 won't make a difference: to_sieve(start + 1) will still
        // point at start. Same with sieve_end.
        let sieve_start = to_sieve(start + 1);
        let sieve_end = to_sieve(end + 1);
        let sieve_size = sieve_end - sieve_start;
        let mut sieve = vec![true; sieve_size];
        // Don't sieve 2 if it exists.
        for &p in origin.primes.iter().skip(1) {
            // Optimize by starting the multiples search at p^2, or the first odd multiple of p
            // after start, whichever is greater.
            let p_sieve_start =
                to_sieve(cmp::max(p * p, p * ceil_odd(ceil_div(start, p)))) - sieve_start;
            for multiple in (p_sieve_start..sieve_size).step_by(p as usize) {
                sieve[multiple] = false;
            }
        }

        let mut primes = Vec::new();
        // We've only sieved odd primes, so 2 will need to be manually prepended, if necessary.
        if start <= 2 && end > 2 {
            primes.push(2);
        }
        primes.extend(sieve.into_iter().enumerate().filter_map(|(p, x)| {
            if x {
                Some(to_prime(p + sieve_start))
            } else {
                None
            }
        }));
        SieveSegment { primes, start, end }
    }
}

/// Convenience functions to convert between prime space (u64 numbers) and sieve space
/// (usizes corresponding to encoded odd u64s).
fn to_sieve(prime: u64) -> usize {
    ((prime - 1) / 2) as usize
}
fn to_prime(sieve: usize) -> u64 {
    2 * sieve as u64 + 1
}
/// Convenience functions to find u64 ceilings
fn ceil_div(numerator: u64, divisor: u64) -> u64 {
    numerator / divisor + ((numerator % divisor) != 0) as u64
}
fn ceil_odd(n: u64) -> u64 {
    n + (n + 1) % 2
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

    #[test]
    fn segmented_bench() {
        let sieve = Sieve::segmented(10_u64.pow(9));
    }
}
