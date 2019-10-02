/// Library for utilities related to primes

use std::slice;

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

/// Unit tests

#[test]
fn test_gcd() {
    assert_eq!(4, gcd(8, 12));
    assert_eq!(4, gcd(12, 8));
    assert_eq!(15, gcd(15, 15));
    assert_eq!(1, gcd(47, 23));
}

#[test]
fn test_sieve_eratosthenes() {
    assert_eq!(vec![&0; 0], Sieve::eratosthenes(0).iter().collect::<Vec<_>>());
    assert_eq!(vec![&0; 0], Sieve::eratosthenes(1).iter().collect::<Vec<_>>());
    assert_eq!(vec![&2], Sieve::eratosthenes(2).iter().collect::<Vec<_>>());
    assert_eq!(vec![&2, &3], Sieve::eratosthenes(3).iter().collect::<Vec<_>>());
    assert_eq!(vec![&2, &3, &5, &7, &11, &13, &17, &19],
               Sieve::eratosthenes(20).iter().collect::<Vec<_>>());
}

#[test]
fn test_sieve_pythagorean() {
    assert_eq!(vec![&0; 0], Sieve::pythagorean(0).iter().collect::<Vec<_>>());
    assert_eq!(vec![&0; 0], Sieve::pythagorean(1).iter().collect::<Vec<_>>());
    assert_eq!(vec![&5], Sieve::pythagorean(5).iter().collect::<Vec<_>>());
    assert_eq!(vec![&5, &13], Sieve::pythagorean(13).iter().collect::<Vec<_>>());
    assert_eq!(vec![&5, &13, &17, &29, &37, &41],
               Sieve::pythagorean(50).iter().collect::<Vec<_>>());
}

#[test]
fn test_is_coprime_eratosthenes() {
    let sieve = Sieve::eratosthenes(9);
    assert_eq!(false, sieve.is_coprime(0));
    assert_eq!(true, sieve.is_coprime(1));
    assert_eq!(false, sieve.is_coprime(2));
    assert_eq!(false, sieve.is_coprime(3));
    assert_eq!(false, sieve.is_coprime(4));
    assert_eq!(false, sieve.is_coprime(5));
}

#[test]
fn test_is_coprime_pythagorean() {
    let sieve = Sieve::pythagorean(9);
    assert_eq!(false, sieve.is_coprime(0));
    assert_eq!(true, sieve.is_coprime(1));
    assert_eq!(true, sieve.is_coprime(2));
    assert_eq!(true, sieve.is_coprime(3));
    assert_eq!(true, sieve.is_coprime(4));
    assert_eq!(false, sieve.is_coprime(5));
}

#[test]
fn test_euler_totient() {
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
fn test_zero_euler_totient() {
    let sieve = Sieve::eratosthenes(7);
    sieve.euler_totient(0);
}

#[test]
#[should_panic]
fn test_sieve_too_small_euler_totient() {
    let sieve = Sieve::eratosthenes(7);
    sieve.euler_totient(50);
}