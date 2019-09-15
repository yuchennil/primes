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

/// Run the Sieve of Eratosthenes to generate primes at or below limit.
///
/// Optimize by starting the multiples search at p^2, p^2 + p, ...
/// instead of 2p, 3p, ...
pub fn sieve_of_eratosthenes(limit: u64) -> Vec<u64> {
    let limit = limit as usize;
    let mut sieve = vec![true; limit];
    let mut p = 2;
    while p * p < limit {
        for multiple_p in (p*p..limit).step_by(p) {
            sieve[multiple_p] = false;
        }
        match (p+1..limit).find(|&x| sieve[x]) {
            Some(next_p) => p = next_p,
            None => break,
        }
    }
    (2..limit).filter(|&x| sieve[x]).map(|x| x as u64).collect()
}

// Use a prime sieve to factor n, then calculate phi with the product formula
// phi(n) = n product_{distinct p|n} (1 - 1/p)
//
// sieve must contain all primes less than n.
pub fn euler_totient(n: u64, sieve: &Vec<u64>) -> u64 {
    assert_ne!(0, n, "euler_totient(n) undefined for n = 0");
    let mut totient = n;
    for p in sieve {
        if totient % p == 0 {
            totient -= totient / p;
        }
        if p >= &totient {
            break;
        }
    }
    totient
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
fn test_sieve_of_eratosthenes() {
    assert_eq!(vec![0; 0], sieve_of_eratosthenes(0));
    assert_eq!(vec![0; 0], sieve_of_eratosthenes(1));
    assert_eq!(vec![0; 0], sieve_of_eratosthenes(2));
    assert_eq!(vec![2], sieve_of_eratosthenes(3));
    assert_eq!(vec![2, 3, 5, 7, 11, 13, 17, 19], sieve_of_eratosthenes(20));
}

#[test]
fn test_euler_totient() {
    let sieve = sieve_of_eratosthenes(100);
    assert_eq!(1, euler_totient(1, &sieve));
    assert_eq!(1, euler_totient(2, &sieve));
    assert_eq!(2, euler_totient(3, &sieve));
    assert_eq!(2, euler_totient(4, &sieve));
    assert_eq!(16, euler_totient(17, &sieve));
    assert_eq!(8, euler_totient(24, &sieve));
    assert_eq!(54, euler_totient(81, &sieve));
}

#[test]
#[should_panic]
fn test_zero_euler_totient() {
    let sieve = sieve_of_eratosthenes(10);
    euler_totient(0, &sieve);
}