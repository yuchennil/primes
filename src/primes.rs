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
}
