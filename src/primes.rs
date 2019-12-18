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

    pub fn iter(&self) -> slice::Iter<u64> {
        self.primes.iter()
    }

    pub fn into_iter(self) -> std::vec::IntoIter<u64> {
        self.primes.into_iter()
    }
}

pub struct Factory {
    sieve: Sieve,
}

impl Factory {

    pub fn new(limit: u64) -> Factory {
        Factory {
            sieve: Sieve::eratosthenes(limit),
        }
    }

    pub fn prime_factor(&self, n: u64) -> Vec<(u64, u32)> {
        let mut prime_factors = collections::BTreeMap::new();
        for p in self.trial_division(n) {
            *prime_factors.entry(p).or_insert(0) += 1;
        }
        prime_factors.into_iter().collect()
    }

    pub fn all_factors(&self, n: u64) -> Vec<u64> {
        let prime_factors = self.prime_factor(n);
        let mut factors = Vec::new();
        for factor_combination in combinatorics::CombinationsWithRepeats::new(&prime_factors) {
            factors.push(factor_combination.into_iter().map(
                |(prime, power)| prime.pow(power)).fold(1, |a, b| a * b));
        }
        factors.sort();
        factors
    }

    /// Repeatedly divide by primes smaller than n until n == 1.
    fn trial_division(&self, n: u64) -> Vec<u64> {
        assert_ne!(0, n, "prime factors of n undefined for n == 0");
        assert!(n <= self.sieve.limit, "Sieve must contain all primes at or below n");
        let mut n = n;
        let mut prime_factors = Vec::new();
        for &p in self.sieve.iter() {
            while n % p == 0 {
                prime_factors.push(p);
                n /= p;
            }
            if n == 1 {
                break;
            }
        }
        prime_factors
    }

    /// Pollard's rho Monte Carlo factor-finding with Floyd's cycle detection
    /// Return a nontrivial factor, or failure.
    fn pollards_rho(n: u64) -> Option<u64> {
        fn f(n: u64) -> u64 {
            n * n + 1
        }
        fn difference(a: u64, b: u64) -> u64 {
            if a > b {
                a - b
            } else {
                b - a
            }
        }
        if n < 2 {
            return None
        }
        if n == 2 {
            return Some(2)
        }
        let mut x = 2;
        let mut y = 2;
        let mut d = 1;
        while d == 1 {
            x = f(x) % n;
            y = f(f(y)) % n;
            d = gcd(difference(x, y), n);
        }
        match d == n {
            true => None,
            false => Some(d),
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
        let factory = Factory::new(12);
        assert_eq!(vec![(0, 0); 0], factory.prime_factor(1));
        assert_eq!(vec![(2, 1)], factory.prime_factor(2));
        assert_eq!(vec![(3, 1)], factory.prime_factor(3));
        assert_eq!(vec![(2, 2)], factory.prime_factor(4));
        assert_eq!(vec![(2, 2), (3, 1)], factory.prime_factor(12));
    }

    #[test]
    #[should_panic]
    fn prime_factor_zero() {
        let factory = Factory::new(7);
        factory.prime_factor(0);
    }

    #[test]
    #[should_panic]
    fn prime_factor_sieve_too_small() {
        let factory = Factory::new(7);
        factory.prime_factor(50);
    }

    #[test]
    fn all_factors_correct() {
        let factory = Factory::new(12);
        assert_eq!(vec![1], factory.all_factors(1));
        assert_eq!(vec![1, 2], factory.all_factors(2));
        assert_eq!(vec![1, 3], factory.all_factors(3));
        assert_eq!(vec![1, 2, 4], factory.all_factors(4));
        assert_eq!(vec![1, 2, 3, 4, 6, 12], factory.all_factors(12));
    }

    #[test]
    fn pollards_rho_correct() {
        assert_eq!(None, Factory::pollards_rho(0));
        assert_eq!(None, Factory::pollards_rho(1));
        assert_eq!(Some(2), Factory::pollards_rho(2));
        assert_eq!(Some(3), Factory::pollards_rho(2 * 3));
        assert_eq!(Some(5), Factory::pollards_rho(5 * 17));
        assert_eq!(Some(271), Factory::pollards_rho(271 * 241));
    }
}