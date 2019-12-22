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

        let size_limit = limit as usize + 1;
        let mut sieve = vec![true; size_limit];
        sieve[0] = false;
        sieve[1] = false;

        let mut p = 2;
        while p * p < size_limit {
            // Optimize by starting the multiples search at p^2, p^2 + p, ...
            // instead of 2p, 3p, ...
            for multiple_p in (p * p..size_limit).step_by(p) {
                sieve[multiple_p] = false;
            }
            match (p + 1..size_limit).find(|&x| sieve[x]) {
                Some(next_p) => p = next_p,
                None => break,
            }
        }

        Sieve {
            primes: sieve
                .iter()
                .enumerate()
                .filter_map(|(i, &x)| if x { Some(i as u64) } else { None })
                .collect(),
            limit,
        }
    }

    /// Filter from a prime sieve to get Pythagorean primes satisfying p % 4 == 1
    pub fn pythagorean(limit: u64) -> Sieve {
        Sieve {
            primes: Sieve::eratosthenes(limit)
                .into_iter()
                .filter(|p| p % 4 == 1)
                .collect(),
            limit,
        }
    }

    /// Return true iff n is coprime to this sieve's primes. If this sieve's primes are complete
    /// then this should always return false except when n == 1.
    pub fn is_coprime(&self, n: u64) -> bool {
        assert!(
            self.limit >= n,
            "Sieve must contain all primes at or below n"
        );
        for &p in self.primes.iter() {
            if n % p == 0 {
                return false;
            }
            if n < p {
                break;
            }
        }
        true
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
            vec![2, 3, 5, 7, 11, 13, 17, 19],
            Sieve::eratosthenes(20).into_iter().collect::<Vec<_>>()
        );
    }

    #[test]
    fn sieve_pythagorean() {
        assert_eq!(
            vec![0; 0],
            Sieve::pythagorean(0).into_iter().collect::<Vec<_>>()
        );
        assert_eq!(
            vec![0; 0],
            Sieve::pythagorean(1).into_iter().collect::<Vec<_>>()
        );
        assert_eq!(
            vec![5],
            Sieve::pythagorean(5).into_iter().collect::<Vec<_>>()
        );
        assert_eq!(
            vec![5, 13],
            Sieve::pythagorean(13).into_iter().collect::<Vec<_>>()
        );
        assert_eq!(
            vec![5, 13, 17, 29, 37, 41],
            Sieve::pythagorean(50).into_iter().collect::<Vec<_>>()
        );
    }

    #[test]
    fn is_coprime_eratosthenes() {
        let sieve = Sieve::eratosthenes(9);
        assert_eq!(false, sieve.is_coprime(0));
        assert_eq!(true, sieve.is_coprime(1));
        assert_eq!(false, sieve.is_coprime(2));
        assert_eq!(false, sieve.is_coprime(3));
        assert_eq!(false, sieve.is_coprime(4));
        assert_eq!(false, sieve.is_coprime(5));
    }

    #[test]
    fn is_coprime_pythagorean() {
        let sieve = Sieve::pythagorean(9);
        assert_eq!(false, sieve.is_coprime(0));
        assert_eq!(true, sieve.is_coprime(1));
        assert_eq!(true, sieve.is_coprime(2));
        assert_eq!(true, sieve.is_coprime(3));
        assert_eq!(true, sieve.is_coprime(4));
        assert_eq!(false, sieve.is_coprime(5));
    }
}
