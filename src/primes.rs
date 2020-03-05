/// Library for utilities related to primes
use std::cmp;
use std::mem;

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

/// {2}-wheel segmented sieve of Eratosthenes to generate all primes below a given limit
///
/// The naive sieve of Eratosthenes strikes multiples of a given prime from a fixed array. The
/// first unstruck number in the array is then the next smallest prime, and we repeat this process
/// until we've consumed all primes in the array.
///
/// The naive sieve requires keeping the entire array in memory, with frequent cache misses due to
/// striding through the whole array once for each prime. By contrast, a segmented Sieve reduces
/// memory use and improves locality by partitioning numbers below limit into contiguous segments.
/// Once we've processed a segment, its memory can be reset and reused. Moreover, if a segment fits
/// in a CPU's L1/L2/L3 cache, then we can strike all primes within it without extra memory loads.
///
/// A {2, 3, ..., p_k}-wheel sieve further optimizes memory by skipping multiples of the
/// small primes 2, 3, ..., p_k. The only prime candidates we need to consider will be coprime
/// to these small primes, saving us 1/2 strikes for a {2}-wheel, 4/6 strikes for a {2, 3}-wheel,
/// etc.
///
/// Usage:
///
///     use primes::primes::Sieve;
///
///     assert_eq!(vec![2, 3, 5, 7, 11, 13, 17, 19], Sieve::segmented(20).collect::<Vec<_>>());
///     assert_eq!(vec![83, 89, 97], Sieve::range(80, 100).collect::<Vec<_>>());
///
/// Algorithm primarily based on Jonathan Sorenson's 1990 "An Introduction to Prime Number Sieves":
/// - https://minds.wisconsin.edu/handle/1793/59248
/// Modern CPU + RAM optimizations due to Kim Walisch's primesieve:
/// - https://github.com/kimwalisch/primesieve/wiki/Segmented-sieve-of-Eratosthenes
pub struct Sieve {
    state_machine: SieveStateMachine,
}

impl Sieve {
    const WHEEL_BASIS_PRIME: usize = 2;
    const FIRST_NON_BASIS_PRIME: usize = 3;

    pub fn segmented(limit: u64) -> Sieve {
        Sieve::range(0, limit)
    }

    pub fn range(start: u64, end: u64) -> Sieve {
        let limit = end as usize;
        // origin_limit is just above sqrt(limit) so that the origin primes suffice to sieve
        // all remaining segments (which hence don't need to be kept in memory after we've finished
        // sieving through them).
        let origin_limit = (limit as f64).sqrt().ceil() as usize;
        let origin_primes = Sieve::sieve_origin(origin_limit);

        let state_machine =
            SieveStateMachine::new(start as usize, limit, origin_limit, origin_primes);

        Sieve { state_machine }
    }

    /// Sieve an origin segment [0, origin_limit) using Eratosthenes, skipping non-wheel numbers.
    fn sieve_origin(origin_limit: usize) -> Vec<usize> {
        let mut origin_primes = Vec::new();

        let mut segment = SieveSegment::new(0, origin_limit);
        let mut n = Sieve::FIRST_NON_BASIS_PRIME;
        while let Some(p) = segment.find_prime(n) {
            origin_primes.push(p);
            // Optimize by striking multiples from p^2. Smaller multiples should already
            // have been struck by previous primes.
            segment.strike_prime_and_get_next_multiple(p, p * p);
            n = p + 2;
        }
        origin_primes
    }
}

impl Iterator for Sieve {
    type Item = u64;

    fn next(&mut self) -> Option<Self::Item> {
        loop {
            // If the current state has a next value, return it.
            use SieveStateMachine::*;
            match self.state_machine {
                Basis(ref mut sieve_state) => {
                    if let Some(p) = sieve_state.next() {
                        return Some(p);
                    }
                }
                Origin(ref mut sieve_state) => {
                    if let Some(p) = sieve_state.next() {
                        return Some(p);
                    }
                }
                Wheel(ref mut sieve_state) => {
                    if let Some(p) = sieve_state.next() {
                        return Some(p);
                    }
                }
                Done => return None,
            };
            // Otherwise step to the next state of the state machine and try again.
            self.state_machine.step();
        }
    }
}

enum SieveStateMachine {
    Basis(SieveState<Basis>),
    Origin(SieveState<Origin>),
    Wheel(SieveState<Wheel>),
    Done,
}

impl SieveStateMachine {
    fn new(
        n: usize,
        limit: usize,
        origin_limit: usize,
        origin_primes: Vec<usize>,
    ) -> SieveStateMachine {
        SieveStateMachine::Basis(SieveState {
            limit,
            origin_limit,
            origin_primes,
            state: Basis { n },
        })
    }

    /// Step the state machine in place (consuming the previous state without moving it).
    ///
    /// We have to mem swap a temporary state machine because callers to step only have a
    /// mutable reference to self.
    fn step(&mut self) {
        use SieveStateMachine::*;

        let mut next = Done;
        mem::swap(self, &mut next);
        next = match next {
            Basis(sieve_state) => Origin(sieve_state.into()),
            Origin(sieve_state) => Wheel(sieve_state.into()),
            Wheel(sieve_state) if !sieve_state.done() => Wheel(sieve_state.advance_segment()),
            Wheel(_) | Done => Done,
        };
        mem::swap(self, &mut next);
    }
}

struct SieveState<S> {
    limit: usize,
    origin_limit: usize,
    origin_primes: Vec<usize>,
    state: S,
}

struct Basis {
    n: usize,
}

impl Iterator for SieveState<Basis> {
    type Item = u64;

    fn next(&mut self) -> Option<Self::Item> {
        if self.state.n <= Sieve::WHEEL_BASIS_PRIME {
            self.state.n = Sieve::FIRST_NON_BASIS_PRIME;
            if Sieve::WHEEL_BASIS_PRIME < self.limit {
                return Some(Sieve::WHEEL_BASIS_PRIME as u64);
            }
        }
        None
    }
}

struct Origin {
    origin_prime_index: usize,
    n: usize,
}

impl Iterator for SieveState<Origin> {
    type Item = u64;

    fn next(&mut self) -> Option<Self::Item> {
        if let Some(&p) = self.origin_primes.get(self.state.origin_prime_index) {
            self.state.origin_prime_index += 1;
            return Some(p as u64);
        }
        self.state.n = cmp::max(self.state.n, self.origin_limit);
        None
    }
}

impl From<SieveState<Basis>> for SieveState<Origin> {
    fn from(sieve_state: SieveState<Basis>) -> SieveState<Origin> {
        let limit = sieve_state.limit;
        let origin_limit = sieve_state.origin_limit;
        let origin_primes = sieve_state.origin_primes;

        let n = sieve_state.state.n;
        let origin_prime_index = if n < origin_limit {
            origin_primes
                .iter()
                .position(|&p| p >= n)
                .unwrap_or_else(|| origin_primes.len())
        } else {
            origin_primes.len()
        };

        SieveState {
            limit,
            origin_limit,
            origin_primes,
            state: Origin {
                n,
                origin_prime_index,
            },
        }
    }
}

struct Wheel {
    segment: SieveSegment,
    multiples: Vec<usize>,
    segment_start: usize,
    segment_end: usize,
    n: usize,
}

impl From<SieveState<Origin>> for SieveState<Wheel> {
    fn from(sieve_state: SieveState<Origin>) -> SieveState<Wheel> {
        let limit = sieve_state.limit;
        let origin_limit = sieve_state.origin_limit;
        let origin_primes = sieve_state.origin_primes;

        let n = Wheel::ceil_wheel(sieve_state.state.n);
        let segment_start = cmp::min(n, limit);
        let segment_end = cmp::min(segment_start + origin_limit, limit);

        let segment = SieveSegment::new(segment_start, segment_end);
        let multiples = Wheel::create_multiples(&origin_primes, segment_start);

        let mut state = Wheel {
            segment,
            multiples,
            segment_start,
            segment_end,
            n,
        };
        state.sieve_segment(&origin_primes);

        SieveState {
            limit,
            origin_limit,
            origin_primes,
            state,
        }
    }
}

impl Iterator for SieveState<Wheel> {
    type Item = u64;

    fn next(&mut self) -> Option<Self::Item> {
        if let Some(p) = self.state.segment.find_prime(self.state.n) {
            self.state.n = p + 2;
            return Some(p as u64);
        }
        None
    }
}

impl SieveState<Wheel> {
    // Was this the last segment to sieve?
    fn done(&self) -> bool {
        self.state.segment_end == self.limit
    }

    // Sieve the next segment, consuming and returning self
    fn advance_segment(mut self) -> Self {
        let segment_length = self.state.segment_end - self.state.segment_start;

        self.state.segment_start = self.state.segment_end;
        self.state.segment_end = cmp::min(self.state.segment_start + segment_length, self.limit);
        self.state.n = Wheel::ceil_wheel(self.state.segment_start);

        self.state.sieve_segment(&self.origin_primes);

        self
    }
}

impl Wheel {
    // Initialize a vector of prime multiples starting at segment_start
    fn create_multiples(origin_primes: &[usize], segment_start: usize) -> Vec<usize> {
        let mut multiples = Vec::new();
        for &p in origin_primes {
            let multiple = p * cmp::max(p, Wheel::ceil_wheel(Wheel::ceil_div(segment_start, p)));
            multiples.push(multiple);
        }
        multiples
    }

    /// Sieve a segment [start, end) based on an origin segment.
    fn sieve_segment(&mut self, origin_primes: &[usize]) {
        self.segment.reset(self.segment_start, self.segment_end);
        // Optimize by starting the multiples search at the first wheel multiple of p after start.
        // This should already be set in self.multiples
        for (&p, multiple) in origin_primes.iter().zip(self.multiples.iter_mut()) {
            *multiple = self
                .segment
                .strike_prime_and_get_next_multiple(p, *multiple);
        }
    }

    // Return the first wheel number at least as large as n
    fn ceil_wheel(n: usize) -> usize {
        SieveSegment::sieve_segment_start_to_n(SieveSegment::n_to_sieve_segment_start(n))
    }
    // Return the dividend a / b, rounded up
    fn ceil_div(a: usize, b: usize) -> usize {
        a / b + (a % b != 0) as usize
    }
}

struct SieveSegment {
    sieve: Vec<bool>,
    sieve_segment_start: usize,
    sieve_segment_length: usize,
}

impl SieveSegment {
    /// Create an unsieved SieveSegment in [start, end).
    fn new(segment_start: usize, segment_end: usize) -> SieveSegment {
        let mut sieve_segment = SieveSegment {
            sieve: Vec::new(),
            sieve_segment_start: 0,
            sieve_segment_length: 0,
        };
        sieve_segment.reset(segment_start, segment_end);
        sieve_segment
    }

    /// Reset this SieveSegment in-place to avoid constructor/deconstructor costs.
    fn reset(&mut self, segment_start: usize, segment_end: usize) {
        self.sieve_segment_start = SieveSegment::n_to_sieve_segment_start(segment_start);
        self.sieve_segment_length = self.n_to_sieve(segment_end);
        self.sieve = vec![true; self.sieve_segment_length];
    }

    /// Strike multiples of prime in sieve.
    ///
    /// Return the next multiple past the end of this SieveSegment.
    ///
    /// Note that a step size of p in the sieve corresponds to a step of 2 * p in u64s.
    fn strike_prime_and_get_next_multiple(&mut self, p: usize, multiple: usize) -> usize {
        let mut sieve_segment_multiple = self.n_to_sieve(multiple);
        while sieve_segment_multiple < self.sieve_segment_length {
            self.sieve[sieve_segment_multiple] = false;
            sieve_segment_multiple += p;
        }
        self.sieve_to_n(sieve_segment_multiple)
    }

    /// Find the next prime at or after n in the sieve.
    fn find_prime(&self, n: usize) -> Option<usize> {
        for sieve_segment_n in self.n_to_sieve(n)..self.sieve_segment_length {
            if self.sieve[sieve_segment_n] {
                let p = self.sieve_to_n(sieve_segment_n);
                return Some(p);
            }
        }
        None
    }

    /// Convert between number space and sieve space.
    ///
    /// If n = 2k + 1 is odd, then sieve_to_n(n_to_sieve(n)) = 2k + 1 is an identity, as desired.
    ///
    /// But if n = 2k is even, then sieve_to_n(n_to_sieve(n)) = 2k + 1 = n + 1 is not an identity.
    /// This error is acceptable because primes past 2 are all odd, so:
    /// - a lower bound of 2k + 1 won't skip any prime candidate.
    /// - an upper bound of 2k + 1 means 2k - 1 is the last prime candidate. 2k isn't possible.
    fn n_to_sieve_segment_start(n: usize) -> usize {
        n / 2
    }
    fn sieve_segment_start_to_n(sieve_segment_start_n: usize) -> usize {
        2 * sieve_segment_start_n + 1
    }
    fn n_to_sieve(&self, n: usize) -> usize {
        SieveSegment::n_to_sieve_segment_start(n) - self.sieve_segment_start
    }
    fn sieve_to_n(&self, sieve_n: usize) -> usize {
        SieveSegment::sieve_segment_start_to_n(sieve_n + self.sieve_segment_start)
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
    fn sieve_segmented_correct() {
        assert_eq!(vec![0; 0], Sieve::segmented(0).collect::<Vec<_>>());
        assert_eq!(vec![0; 0], Sieve::segmented(1).collect::<Vec<_>>());
        assert_eq!(vec![0; 0], Sieve::segmented(2).collect::<Vec<_>>());
        assert_eq!(vec![2], Sieve::segmented(3).collect::<Vec<_>>());
        assert_eq!(vec![2, 3], Sieve::segmented(4).collect::<Vec<_>>());
        assert_eq!(vec![2, 3, 5, 7], Sieve::segmented(9).collect::<Vec<_>>());
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
    fn sieve_range_correct() {
        assert_eq!(vec![2, 3, 5, 7], Sieve::range(0, 10).collect::<Vec<_>>());
        assert_eq!(vec![2, 3, 5, 7], Sieve::range(1, 11).collect::<Vec<_>>());
        assert_eq!(
            vec![2, 3, 5, 7, 11],
            Sieve::range(2, 12).collect::<Vec<_>>()
        );
        assert_eq!(vec![3, 5, 7, 11], Sieve::range(3, 13).collect::<Vec<_>>());
        assert_eq!(vec![5, 7, 11, 13], Sieve::range(4, 14).collect::<Vec<_>>());
        assert_eq!(vec![83, 89, 97], Sieve::range(80, 100).collect::<Vec<_>>());
    }

    #[test]
    fn sieve_segment_correct() {
        let mut sieve_segment = SieveSegment::new(5, 16);

        // find_prime() when all true
        assert_eq!(Some(5), sieve_segment.find_prime(5));
        assert_eq!(Some(7), sieve_segment.find_prime(7));
        assert_eq!(Some(9), sieve_segment.find_prime(9));
        assert_eq!(Some(11), sieve_segment.find_prime(11));
        assert_eq!(Some(13), sieve_segment.find_prime(13));
        assert_eq!(Some(15), sieve_segment.find_prime(15));
        assert_eq!(None, sieve_segment.find_prime(17));

        // strike_prime_and_get_next_multiple()
        assert_eq!(21, sieve_segment.strike_prime_and_get_next_multiple(3, 9));
        assert_eq!(25, sieve_segment.strike_prime_and_get_next_multiple(5, 25));

        // find_prime() after sieving
        assert_eq!(Some(5), sieve_segment.find_prime(5));
        assert_eq!(Some(7), sieve_segment.find_prime(7));
        assert_eq!(Some(11), sieve_segment.find_prime(9));
        assert_eq!(Some(13), sieve_segment.find_prime(13));
        assert_eq!(None, sieve_segment.find_prime(15));

        // reset()
        sieve_segment.reset(16, 22);

        // find_prime() after reset
        assert_eq!(Some(17), sieve_segment.find_prime(17));
        assert_eq!(Some(19), sieve_segment.find_prime(19));
        assert_eq!(Some(21), sieve_segment.find_prime(21));
        assert_eq!(None, sieve_segment.find_prime(23));
    }

    #[test]
    fn sieve_to_n_to_sieve_correct() {
        let sieve_segment = SieveSegment::new(5, 10);
        assert_eq!(5, sieve_segment.n_to_sieve(sieve_segment.sieve_to_n(5)));
        assert_eq!(6, sieve_segment.n_to_sieve(sieve_segment.sieve_to_n(6)));
        assert_eq!(7, sieve_segment.n_to_sieve(sieve_segment.sieve_to_n(7)));
        assert_eq!(8, sieve_segment.n_to_sieve(sieve_segment.sieve_to_n(8)));
        assert_eq!(9, sieve_segment.n_to_sieve(sieve_segment.sieve_to_n(9)));
    }
}
