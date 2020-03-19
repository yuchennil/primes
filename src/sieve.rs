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
