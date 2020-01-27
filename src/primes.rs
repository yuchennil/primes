/// Library for utilities related to primes

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
    primes: Vec<u64>,
    limit: u64,
}

impl Sieve {
    /// Run the Sieve of Eratosthenes to generate primes at or below limit
    pub fn eratosthenes(limit: u64) -> Sieve {
        if limit < 2 {
            return Sieve {
                primes: Vec::new(),
                limit,
            };
        }

        fn to_sieve(prime: u64) -> usize {
            ((prime - 1) / 2) as usize
        }
        fn to_prime(sieve: usize) -> u64 {
            2 * sieve as u64 + 1
        }

        let mut sieve = vec![true; to_sieve(limit) + 1];
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

        Sieve {
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
            limit,
        }
    }
}

impl IntoIterator for Sieve {
    type Item = u64;
    type IntoIter = std::vec::IntoIter<Self::Item>;

    fn into_iter(self) -> Self::IntoIter {
        self.primes.into_iter()
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
            vec![2],
            Sieve::eratosthenes(2).into_iter().collect::<Vec<_>>()
        );
        assert_eq!(
            vec![2, 3],
            Sieve::eratosthenes(3).into_iter().collect::<Vec<_>>()
        );
        assert_eq!(
            vec![2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47],
            Sieve::eratosthenes(50).into_iter().collect::<Vec<_>>()
        );
    }
}
