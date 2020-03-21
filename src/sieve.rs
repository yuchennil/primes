use crate::sieve_state_machine::SieveStateMachine;

/// {2, 3, 5, 7}-wheel segmented sieve of Eratosthenes to generate all primes below a given end
///
/// The naive sieve of Eratosthenes strikes multiples of a given prime from a fixed array. The
/// first unstruck number in the array is then the next smallest prime, and we repeat this process
/// until we've consumed all primes in the array.
///
/// The naive sieve requires keeping the entire array in memory, with frequent cache misses due to
/// striding through the whole array once for each prime. By contrast, a segmented Sieve reduces
/// memory use and improves locality by partitioning numbers below end into contiguous segments.
/// Once we've processed a segment, its memory can be discarded. Moreover, if a segment fits in
/// a CPU's L1/L2/L3 cache, then we can strike all primes within it without extra memory loads.
///
/// In Rust a bool is represented with a single byte. Keeping a vector of bits saves eight times
/// the memory compared to a vector of bools representing the same segment, although access costs
/// a couple extra CPU cycles due to bit operations.
///
/// A {2, 3, ..., p_k}-wheel sieve further optimizes memory by skipping multiples of the
/// small primes 2, 3, ..., p_k. The only prime candidates we need to consider will be coprime
/// to these small primes, saving us 1/2 strikes for a {2}-wheel, 4/6 strikes for a {2, 3}-wheel,
/// etc.
///
/// Usage:
///
///     use primes::Sieve;
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
    pub fn segmented(end: u64) -> Sieve {
        Sieve::range(0, end)
    }

    pub fn range(start: u64, end: u64) -> Sieve {
        Sieve {
            state_machine: SieveStateMachine::new(start as usize, end as usize),
        }
    }
}

impl Iterator for Sieve {
    type Item = u64;

    #[inline]
    fn next(&mut self) -> Option<Self::Item> {
        loop {
            if let Some(p) = self.state_machine.next() {
                return Some(p as u64);
            }
            self.state_machine.step();
            if let SieveStateMachine::Done = self.state_machine {
                return None;
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use crate::constants::SEGMENT_LENGTH;
    use crate::origin::Origin;
    use crate::segment::WheelSegment;

    use super::*;

    #[bench]
    fn enumerate_primes(b: &mut test::Bencher) {
        b.iter(|| {
            for _ in Sieve::segmented(1_000_000) {
                continue;
            }
        })
    }

    #[bench]
    fn only_sieve(b: &mut test::Bencher) {
        b.iter(|| {
            let mut state = SieveStateMachine::new(0, 1_000_000);
            loop {
                if let SieveStateMachine::Done = state {
                    break;
                }
                state.step();
            }
        })
    }

    #[bench]
    fn origin_primes(b: &mut test::Bencher) {
        b.iter(|| Origin::primes(10_000))
    }

    #[bench]
    fn wheel_segment(b: &mut test::Bencher) {
        let origin_primes = Origin::primes(10_000);
        let segment_start = 9_876_543_210;
        let segment_end = segment_start + SEGMENT_LENGTH;
        b.iter(|| WheelSegment::new(&origin_primes, segment_start, segment_end))
    }
}
