/// Library for utilities related to primes

use std::collections;
use std::slice;

use crate::combinatorics;

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
                limit: limit,
            }
        }

        let size_limit = limit as usize + 1;
        let mut sieve = vec![true; size_limit];
        sieve[0] = false;
        sieve[1] = false;

        let mut p = 2;
        while p * p < size_limit {
            // Optimize by starting the multiples search at p^2, p^2 + p, ...
            // instead of 2p, 3p, ...
            for multiple_p in (p*p..size_limit).step_by(p) {
                sieve[multiple_p] = false;
            }
            match (p+1..size_limit).find(|&x| sieve[x]) {
                Some(next_p) => p = next_p,
                None => break,
            }
        }

        Sieve {
            primes: sieve.iter().enumerate().filter_map(
                |(i, &x)| match x {
                    true => Some(i as u64),
                    false => None,
                }).collect(),
            limit: limit,
        }
    }

    /// Filter from a prime sieve to get Pythagorean primes satisfying p % 4 == 1
    pub fn pythagorean(limit: u64) -> Sieve {
        Sieve {
            primes: Sieve::eratosthenes(limit).into_iter().filter(|p| p % 4 == 1).collect(),
            limit: limit,
        }
    }

    /// Return true iff n is coprime to this sieve's primes. If this sieve's primes are complete
    /// then this should always return false except when n == 1.
    pub fn is_coprime(&self, n: u64) -> bool {
        assert!(self.limit >= n, "Sieve must contain all primes at or below n");
        for &p in self.iter() {
            if n % p == 0 {
                return false
            }
            if n < p {
                break;
            }
        }
        true
    }

    // Using a prime sieve to factor n, calculate phi with the product formula
    // phi(n) = n product_{distinct p|n} (1 - 1/p)
    pub fn euler_totient(&self, n: u64) -> u64 {
        assert_ne!(0, n, "euler_totient(n) undefined for n == 0");
        assert!(self.limit * self.limit >= n,
                "Sieve must contain all primes at or below sqrt(n)");

        let mut totient = n;
        for p in self.iter() {
            if totient % p == 0 {
                totient -= totient / p;
            }
            // Shortcut out of the loop when primes get too large
            if p >= &totient || p * p >= n {
                break;
            }
        }
        // The above factorization only checks p <= sqrt(n)
        // If n itself is prime then phi(n) = n - 1
        if n > 1 && totient == n {
            totient -= 1;
        }
        totient
    }

    fn prime_factor(&self, n: u64) -> collections::BTreeMap<u64, u32> {
        assert_ne!(0, n, "prime factors of n undefined for n == 0");
        assert!(n <= self.limit, "Sieve must contain all primes at or below n");
        let mut n = n;
        let mut prime_factors = collections::BTreeMap::new();
        for &p in self.iter() {
            while n % p == 0 {
                *prime_factors.entry(p).or_insert(0) += 1;
                n /= p;
            }
            if n == 1 {
                break;
            }
        }
        prime_factors
    }

    pub fn factor(&self, n: u64) -> Vec<u64> {
        let prime_factors = self.prime_factor(n).into_iter().collect::<Vec<_>>();
        let mut factors = Vec::new();
        for factor_combination in combinatorics::CombinationsWithRepeats::new(&prime_factors) {
            factors.push(factor_combination.into_iter().map(
                |(prime, power)| prime.pow(power)).fold(1, |a, b| a * b));
        }
        factors.sort();
        factors
    }

    pub fn iter(&self) -> slice::Iter<u64> {
        self.primes.iter()
    }

    pub fn into_iter(self) -> std::vec::IntoIter<u64> {
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
        assert_eq!(vec![&0; 0], Sieve::eratosthenes(0).iter().collect::<Vec<_>>());
        assert_eq!(vec![&0; 0], Sieve::eratosthenes(1).iter().collect::<Vec<_>>());
        assert_eq!(vec![&2], Sieve::eratosthenes(2).iter().collect::<Vec<_>>());
        assert_eq!(vec![&2, &3], Sieve::eratosthenes(3).iter().collect::<Vec<_>>());
        assert_eq!(vec![&2, &3, &5, &7, &11, &13, &17, &19],
                   Sieve::eratosthenes(20).iter().collect::<Vec<_>>());
    }

    #[test]
    fn sieve_pythagorean() {
        assert_eq!(vec![&0; 0], Sieve::pythagorean(0).iter().collect::<Vec<_>>());
        assert_eq!(vec![&0; 0], Sieve::pythagorean(1).iter().collect::<Vec<_>>());
        assert_eq!(vec![&5], Sieve::pythagorean(5).iter().collect::<Vec<_>>());
        assert_eq!(vec![&5, &13], Sieve::pythagorean(13).iter().collect::<Vec<_>>());
        assert_eq!(vec![&5, &13, &17, &29, &37, &41],
                   Sieve::pythagorean(50).iter().collect::<Vec<_>>());
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

    #[test]
    fn euler_totient_correct() {
        let sieve = Sieve::eratosthenes(9);
        assert_eq!(1, sieve.euler_totient(1));
        assert_eq!(1, sieve.euler_totient(2));
        assert_eq!(2, sieve.euler_totient(3));
        assert_eq!(2, sieve.euler_totient(4));
        assert_eq!(16, sieve.euler_totient(17));
        assert_eq!(8, sieve.euler_totient(24));
        assert_eq!(54, sieve.euler_totient(81));
    }

    #[test]
    #[should_panic]
    fn euler_totient_zero() {
        let sieve = Sieve::eratosthenes(7);
        sieve.euler_totient(0);
    }

    #[test]
    #[should_panic]
    fn euler_totient_sieve_too_small() {
        let sieve = Sieve::eratosthenes(7);
        sieve.euler_totient(50);
    }

    #[test]
    fn prime_factor_correct() {
        let sieve = Sieve::eratosthenes(12);
        assert_eq!(vec![(&0, &0); 0], sieve.prime_factor(1).iter().collect::<Vec<_>>());
        assert_eq!(vec![(&2, &1)], sieve.prime_factor(2).iter().collect::<Vec<_>>());
        assert_eq!(vec![(&3, &1)], sieve.prime_factor(3).iter().collect::<Vec<_>>());
        assert_eq!(vec![(&2, &2)], sieve.prime_factor(4).iter().collect::<Vec<_>>());
        assert_eq!(vec![(&2, &2), (&3, &1)], sieve.prime_factor(12).iter().collect::<Vec<_>>());
    }

    #[test]
    #[should_panic]
    fn prime_factor_zero() {
        let sieve = Sieve::eratosthenes(7);
        sieve.prime_factor(0);
    }

    #[test]
    #[should_panic]
    fn prime_factor_sieve_too_small() {
        let sieve = Sieve::eratosthenes(7);
        sieve.prime_factor(50);
    }

    #[test]
    fn factor_correct() {
        let sieve = Sieve::eratosthenes(12);
        assert_eq!(vec![&1], sieve.factor(1).iter().collect::<Vec<_>>());
        assert_eq!(vec![&1, &2], sieve.factor(2).iter().collect::<Vec<_>>());
        assert_eq!(vec![&1, &3], sieve.factor(3).iter().collect::<Vec<_>>());
        assert_eq!(vec![&1, &2, &4], sieve.factor(4).iter().collect::<Vec<_>>());
        assert_eq!(vec![&1, &2, &3, &4, &6, &12], sieve.factor(12).iter().collect::<Vec<_>>());
    }
}