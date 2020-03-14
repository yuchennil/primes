use arr_macro::arr;
use std::cmp;
use std::mem;

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
    const BASIS_PRIMES: [usize; 4] = [2, 3, 5, 7];
    const FIRST_NON_BASIS_PRIME: usize = 11;
    const SPOKE_SIZE: usize = 48;
    const WHEEL_SIZE: usize = 210;
    const L1_CACHE_BYTES: usize = 32_768;
    const BITS_PER_BOOL: usize = 8;
    const SEGMENT_LENGTH: usize = Sieve::WHEEL_SIZE * Sieve::L1_CACHE_BYTES * Sieve::BITS_PER_BOOL;

    pub fn segmented(end: u64) -> Sieve {
        Sieve::range(0, end)
    }

    pub fn range(start: u64, end: u64) -> Sieve {
        let start = start as usize;
        let end = end as usize;

        Sieve {
            state_machine: SieveStateMachine::new(start, end),
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

/// There are three kinds of primes produced by this sieve as it iterates:
/// - basis primes
///   - smallest primes below origin
///   - memorized
///   - used to determine the sieving wheel
/// - origin primes
///   - above basis and below sqrt(end)
///   - calculated and permanently stored during initialization
///   - used to strike wheel segments above sqrt(end)
/// - wheel primes
///   - above sqrt(end) and below end
///   - generated by sieving segments with origin primes
///   - discarded after iteration
///
/// The sieve exhibits different behavior while iterating through each of these kinds of primes,
/// so we enforce correct computation with a finite state machine:
///     Basis -> Origin -> Wheel -> Done
///                        ^___|
///
/// As a result, we can transparently get a single iterator through all primes, while ensuring
/// mutable state is managed safely.
///
/// Note that we decided *not* to directly make Sieve an enum because that would expose private
/// implementation details (i.e., the sieve states). Rust considers pub enum variants to be pub.
enum SieveStateMachine {
    Basis(Basis),
    Origin(Origin),
    Wheel(Wheel),
    Done,
}

impl Iterator for SieveStateMachine {
    type Item = usize;

    /// Return the next prime in the sieve, if possible.
    fn next(&mut self) -> Option<Self::Item> {
        use SieveStateMachine::*;

        match self {
            Basis(ref mut state) => state.next(),
            Origin(ref mut state) => state.next(),
            Wheel(ref mut state) => state.next(),
            Done => None,
        }
    }
}

impl SieveStateMachine {
    /// All new SieveStateMachines must start from a Basis state.
    fn new(start: usize, end: usize) -> SieveStateMachine {
        SieveStateMachine::Basis(Basis::new(start, end))
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
            Basis(state) => Origin(state.into()),
            Origin(state) => Wheel(state.into()),
            Wheel(state) if !state.done() => Wheel(state.advance_segment()),
            Wheel(_) | Done => Done,
        };
        mem::swap(self, &mut next);
    }
}

struct Basis {
    start: usize,
    end: usize,
    basis_primes_index: usize,
}

impl Iterator for Basis {
    type Item = usize;

    fn next(&mut self) -> Option<Self::Item> {
        match Sieve::BASIS_PRIMES.get(self.basis_primes_index) {
            Some(&p) if p < self.end => {
                self.basis_primes_index += 1;
                Some(p)
            }
            _ => None,
        }
    }
}

impl Basis {
    fn new(start: usize, end: usize) -> Basis {
        let basis_primes_index = Basis::primes_index(start);

        // Advance start to FIRST_NON_BASIS_PRIME so it will be correct when passed to Origin.
        let start = cmp::max(start, Sieve::FIRST_NON_BASIS_PRIME);

        Basis {
            start,
            end,
            basis_primes_index,
        }
    }

    /// Find the index of the first prime in BASIS_PRIMES greater than n.
    fn primes_index(n: usize) -> usize {
        if n < Sieve::FIRST_NON_BASIS_PRIME {
            Sieve::BASIS_PRIMES
                .iter()
                .position(|&p| p >= n)
                .unwrap_or_else(|| Sieve::BASIS_PRIMES.len())
        } else {
            Sieve::BASIS_PRIMES.len()
        }
    }
}

struct Origin {
    start: usize,
    end: usize,
    origin_primes: Vec<usize>,
    origin_primes_index: usize,
}

impl Iterator for Origin {
    type Item = usize;

    fn next(&mut self) -> Option<Self::Item> {
        let p = *self.origin_primes.get(self.origin_primes_index)?;
        self.origin_primes_index += 1;
        Some(p)
    }
}

impl From<Basis> for Origin {
    fn from(state: Basis) -> Origin {
        let start = state.start;
        let end = state.end;

        let origin_end = Origin::end(end);
        let origin_primes = Origin::primes(origin_end);
        let origin_primes_index = Origin::primes_index(start, origin_end, &origin_primes);

        // Advance start to origin_end so it will be correct when passed to Wheel.
        let start = cmp::max(start, origin_end);

        Origin {
            start,
            end,
            origin_primes,
            origin_primes_index,
        }
    }
}

impl Origin {
    // origin_end is just above sqrt(end) so that the origin primes suffice to sieve
    // all remaining segments (which hence don't need to be kept in memory after we've finished
    // sieving through them).
    fn end(end: usize) -> usize {
        (end as f64).sqrt().ceil() as usize
    }
    /// Sieve an origin segment [0, origin_end) using Eratosthenes, skipping non-wheel numbers.
    fn primes(origin_end: usize) -> Vec<usize> {
        let mut origin_primes = Vec::new();

        let mut segment = Segment::new(0, origin_end);
        let mut n = Sieve::FIRST_NON_BASIS_PRIME;
        while let Some(p) = segment.find_prime(n) {
            origin_primes.push(p);
            segment.strike_primes(&[p]);
            n = p + 1;
        }
        origin_primes
    }

    /// Find the index of the first prime in origin_primes greater than n.
    fn primes_index(n: usize, origin_end: usize, origin_primes: &[usize]) -> usize {
        if n < origin_end {
            origin_primes
                .iter()
                .position(|&p| p >= n)
                .unwrap_or_else(|| origin_primes.len())
        } else {
            origin_primes.len()
        }
    }
}

struct Wheel {
    start: usize,
    end: usize,
    origin_primes: Vec<usize>,
    segment: Segment,
    segment_start: usize,
    segment_end: usize,
}

impl From<Origin> for Wheel {
    fn from(state: Origin) -> Wheel {
        let start = state.start;
        let end = state.end;
        let origin_primes = state.origin_primes;

        let segment_start = cmp::min(start, end);
        let segment_end = cmp::min(segment_start + Sieve::SEGMENT_LENGTH, end);
        let segment = Segment::new(segment_start, segment_end);

        let mut state = Wheel {
            start,
            end,
            origin_primes,
            segment,
            segment_start,
            segment_end,
        };
        state.sieve_segment();
        state
    }
}

impl Iterator for Wheel {
    type Item = usize;

    fn next(&mut self) -> Option<Self::Item> {
        let p = self.segment.find_prime(self.start)?;
        self.start = p + 1;
        Some(p)
    }
}

impl Wheel {
    /// Sieve a segment [start, end) based on an origin segment.
    fn sieve_segment(&mut self) {
        self.segment = Segment::new(self.segment_start, self.segment_end);
        self.segment.strike_primes(&self.origin_primes);
    }

    // Was this the last segment to sieve?
    fn done(&self) -> bool {
        self.segment_end == self.end
    }

    // Sieve the next segment, consuming and returning self
    fn advance_segment(mut self) -> Self {
        self.segment_start = self.segment_end;
        self.segment_end = cmp::min(self.segment_start + Sieve::SEGMENT_LENGTH, self.end);
        self.start = self.segment_start;

        self.sieve_segment();

        self
    }
}

struct Segment {
    spokes: [Spoke; Sieve::SPOKE_SIZE],
    segment_start: usize,
}

impl Segment {
    const WHEEL: [usize; Sieve::SPOKE_SIZE] = [
        1, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97, 101,
        103, 107, 109, 113, 121, 127, 131, 137, 139, 143, 149, 151, 157, 163, 167, 169, 173, 179,
        181, 187, 191, 193, 197, 199, 209,
    ];
    const SPOKE_GAPS: [usize; Sieve::SPOKE_SIZE] = [
        10, 2, 4, 2, 4, 6, 2, 6, 4, 2, 4, 6, 6, 2, 6, 4, 2, 6, 4, 6, 8, 4, 2, 4, 2, 4, 8, 6, 4, 6,
        2, 4, 6, 2, 6, 6, 4, 2, 4, 6, 2, 6, 4, 2, 4, 2, 10, 2,
    ];
    const SPOKE: [usize; Sieve::WHEEL_SIZE] = [
        0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 3, 3, 3, 3, 4, 4, 5, 5, 5, 5, 6, 6, 6, 6, 6, 6,
        7, 7, 8, 8, 8, 8, 8, 8, 9, 9, 9, 9, 10, 10, 11, 11, 11, 11, 12, 12, 12, 12, 12, 12, 13, 13,
        13, 13, 13, 13, 14, 14, 15, 15, 15, 15, 15, 15, 16, 16, 16, 16, 17, 17, 18, 18, 18, 18, 18,
        18, 19, 19, 19, 19, 20, 20, 20, 20, 20, 20, 21, 21, 21, 21, 21, 21, 21, 21, 22, 22, 22, 22,
        23, 23, 24, 24, 24, 24, 25, 25, 26, 26, 26, 26, 27, 27, 27, 27, 27, 27, 27, 27, 28, 28, 28,
        28, 28, 28, 29, 29, 29, 29, 30, 30, 30, 30, 30, 30, 31, 31, 32, 32, 32, 32, 33, 33, 33, 33,
        33, 33, 34, 34, 35, 35, 35, 35, 35, 35, 36, 36, 36, 36, 36, 36, 37, 37, 37, 37, 38, 38, 39,
        39, 39, 39, 40, 40, 40, 40, 40, 40, 41, 41, 42, 42, 42, 42, 42, 42, 43, 43, 43, 43, 44, 44,
        45, 45, 45, 45, 46, 46, 47, 47, 47, 47, 47, 47, 47, 47, 47, 47,
    ];

    /// Create an unsieved Segment in [start, end).
    fn new(segment_start: usize, segment_end: usize) -> Segment {
        let spoke_builder = |spoke| Spoke::new(Segment::WHEEL[spoke], segment_start, segment_end);
        let mut spoke = 0;
        let spokes = arr![spoke_builder({spoke += 1; spoke - 1}); 48];

        Segment {
            spokes,
            segment_start,
        }
    }

    /// Strike primes for all spokes in this segment.
    fn strike_primes(&mut self, primes: &[usize]) {
        let mut multiples = arr![Vec::new(); 48];
        for &p in primes {
            // Optimize by striking multiples from p^2. Smaller multiples should already have been
            // struck by previous primes. Also skip ahead to the first spoke in this segment.
            let factor = cmp::max(p, self.first_wheel_factor(p));
            let mut multiple = p * factor;
            let wheel_iter = Segment::SPOKE_GAPS
                .iter()
                .cycle()
                .skip(Segment::spoke(factor))
                .take(Sieve::SPOKE_SIZE);
            for spoke_gap in wheel_iter {
                multiples[Segment::spoke(multiple)].push(multiple);
                multiple += p * spoke_gap;
            }
        }

        // Iterate through all primes in all spokes in spoke-major order for cache locality.
        for (spoke, spoke_multiples) in self.spokes.iter_mut().zip(multiples.iter()) {
            for (&p, &multiple) in primes.iter().zip(spoke_multiples.iter()) {
                spoke.strike_prime(p, multiple);
            }
        }
    }

    /// Find the next prime at or after n in the segment.
    ///
    /// TODO optimize this to O(log k) instead of O(k)
    fn find_prime(&self, n: usize) -> Option<usize> {
        let mut min_p = None;
        for spoke in self.spokes.iter() {
            if let Some(p) = spoke.find_prime(n) {
                min_p = match min_p {
                    Some(min_p) if min_p <= p => Some(min_p),
                    _ => Some(p),
                }
            }
        }
        min_p
    }

    fn first_wheel_factor(&self, p: usize) -> usize {
        let n = ceil_div(self.segment_start, p);
        (n / Sieve::WHEEL_SIZE) * Sieve::WHEEL_SIZE + Segment::WHEEL[Segment::spoke(n)]
    }
    fn spoke(n: usize) -> usize {
        Segment::SPOKE[n % Sieve::WHEEL_SIZE]
    }
}

struct Spoke {
    sieve: BitVec,
    residue: usize,
    spoke_start: usize,
    spoke_length: usize,
}

impl Spoke {
    /// Create an unsieved Spoke in [start, end).
    fn new(residue: usize, segment_start: usize, segment_end: usize) -> Spoke {
        let spoke_start = Spoke::n_to_spoke_start(segment_start, residue);
        let spoke_length = Spoke::n_to_spoke_start(segment_end, residue) - spoke_start;
        let sieve = BitVec::new(spoke_length);

        Spoke {
            sieve,
            residue,
            spoke_start,
            spoke_length,
        }
    }

    /// Strike multiples of prime in this spoke.
    ///
    /// Note that a step size of p in the spoke corresponds to a step of WHEEL_SIZE * p in u64s.
    fn strike_prime(&mut self, p: usize, multiple: usize) {
        for spoke_multiple in (self.n_to_spoke(multiple)..self.spoke_length).step_by(p) {
            self.sieve.unset(spoke_multiple);
        }
    }

    /// Find the next prime at or after n in the spoke.
    fn find_prime(&self, n: usize) -> Option<usize> {
        let spoke_n = self.n_to_spoke(n);
        let spoke_p = self.sieve.find(spoke_n)?;
        Some(self.spoke_to_n(spoke_p))
    }

    /// Convert between number space and spoke space.
    fn n_to_spoke_start(n: usize, residue: usize) -> usize {
        n / Sieve::WHEEL_SIZE + (n % Sieve::WHEEL_SIZE > residue) as usize
    }
    fn spoke_start_to_n(spoke_start_n: usize, residue: usize) -> usize {
        spoke_start_n * Sieve::WHEEL_SIZE + residue
    }
    fn n_to_spoke(&self, n: usize) -> usize {
        Spoke::n_to_spoke_start(n, self.residue) - self.spoke_start
    }
    fn spoke_to_n(&self, spoke_n: usize) -> usize {
        Spoke::spoke_start_to_n(spoke_n + self.spoke_start, self.residue)
    }
}

struct BitVec(Vec<u64>);

impl BitVec {
    const WORD_BITS: usize = 64;
    const SHIFT: usize = 6;
    const MASK: usize = 0b11_1111;
    const ONES: u64 = std::u64::MAX;
    const UNSET_BIT: [u64; 64] = [
        !(1 << 0),
        !(1 << 1),
        !(1 << 2),
        !(1 << 3),
        !(1 << 4),
        !(1 << 5),
        !(1 << 6),
        !(1 << 7),
        !(1 << 8),
        !(1 << 9),
        !(1 << 10),
        !(1 << 11),
        !(1 << 12),
        !(1 << 13),
        !(1 << 14),
        !(1 << 15),
        !(1 << 16),
        !(1 << 17),
        !(1 << 18),
        !(1 << 19),
        !(1 << 20),
        !(1 << 21),
        !(1 << 22),
        !(1 << 23),
        !(1 << 24),
        !(1 << 25),
        !(1 << 26),
        !(1 << 27),
        !(1 << 28),
        !(1 << 29),
        !(1 << 30),
        !(1 << 31),
        !(1 << 32),
        !(1 << 33),
        !(1 << 34),
        !(1 << 35),
        !(1 << 36),
        !(1 << 37),
        !(1 << 38),
        !(1 << 39),
        !(1 << 40),
        !(1 << 41),
        !(1 << 42),
        !(1 << 43),
        !(1 << 44),
        !(1 << 45),
        !(1 << 46),
        !(1 << 47),
        !(1 << 48),
        !(1 << 49),
        !(1 << 50),
        !(1 << 51),
        !(1 << 52),
        !(1 << 53),
        !(1 << 54),
        !(1 << 55),
        !(1 << 56),
        !(1 << 57),
        !(1 << 58),
        !(1 << 59),
        !(1 << 60),
        !(1 << 61),
        !(1 << 62),
        !(1 << 63),
    ];
    const GREATER_OR_EQUAL_BITS: [u64; 64] = [
        BitVec::ONES << 0,
        BitVec::ONES << 1,
        BitVec::ONES << 2,
        BitVec::ONES << 3,
        BitVec::ONES << 4,
        BitVec::ONES << 5,
        BitVec::ONES << 6,
        BitVec::ONES << 7,
        BitVec::ONES << 8,
        BitVec::ONES << 9,
        BitVec::ONES << 10,
        BitVec::ONES << 11,
        BitVec::ONES << 12,
        BitVec::ONES << 13,
        BitVec::ONES << 14,
        BitVec::ONES << 15,
        BitVec::ONES << 16,
        BitVec::ONES << 17,
        BitVec::ONES << 18,
        BitVec::ONES << 19,
        BitVec::ONES << 20,
        BitVec::ONES << 21,
        BitVec::ONES << 22,
        BitVec::ONES << 23,
        BitVec::ONES << 24,
        BitVec::ONES << 25,
        BitVec::ONES << 26,
        BitVec::ONES << 27,
        BitVec::ONES << 28,
        BitVec::ONES << 29,
        BitVec::ONES << 30,
        BitVec::ONES << 31,
        BitVec::ONES << 32,
        BitVec::ONES << 33,
        BitVec::ONES << 34,
        BitVec::ONES << 35,
        BitVec::ONES << 36,
        BitVec::ONES << 37,
        BitVec::ONES << 38,
        BitVec::ONES << 39,
        BitVec::ONES << 40,
        BitVec::ONES << 41,
        BitVec::ONES << 42,
        BitVec::ONES << 43,
        BitVec::ONES << 44,
        BitVec::ONES << 45,
        BitVec::ONES << 46,
        BitVec::ONES << 47,
        BitVec::ONES << 48,
        BitVec::ONES << 49,
        BitVec::ONES << 50,
        BitVec::ONES << 51,
        BitVec::ONES << 52,
        BitVec::ONES << 53,
        BitVec::ONES << 54,
        BitVec::ONES << 55,
        BitVec::ONES << 56,
        BitVec::ONES << 57,
        BitVec::ONES << 58,
        BitVec::ONES << 59,
        BitVec::ONES << 60,
        BitVec::ONES << 61,
        BitVec::ONES << 62,
        BitVec::ONES << 63,
    ];

    fn new(len: usize) -> BitVec {
        let mut bit_vec = vec![BitVec::ONES; ceil_div(len, BitVec::WORD_BITS)];
        if let Some(end) = bit_vec.get_mut(len >> BitVec::SHIFT) {
            *end &= !BitVec::GREATER_OR_EQUAL_BITS[len & BitVec::MASK];
        }
        BitVec(bit_vec)
    }

    fn find(&self, index: usize) -> Option<usize> {
        let first_word_index = index >> BitVec::SHIFT;
        for (word_index, &word) in self.0[first_word_index..].iter().enumerate() {
            let masked_word = if word_index == 0 {
                word & BitVec::GREATER_OR_EQUAL_BITS[index & BitVec::MASK]
            } else {
                word
            };
            if let Some(bit_index) = BitVec::find_first_set(masked_word) {
                return Some(((first_word_index + word_index) << BitVec::SHIFT) + bit_index);
            }
        }
        None
    }

    /// Find the first set bit in word.
    fn find_first_set(word: u64) -> Option<usize> {
        if word == 0 {
            return None;
        }
        Some(word.trailing_zeros() as usize)
    }

    fn unset(&mut self, index: usize) {
        self.0[index >> BitVec::SHIFT] &= BitVec::UNSET_BIT[index & BitVec::MASK]
    }
}

// Return the dividend a / b, rounded up
fn ceil_div(a: usize, b: usize) -> usize {
    a / b + (a % b != 0) as usize
}

#[cfg(test)]
mod tests {
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

    #[test]
    fn segment_correct() {
        let mut segment = Segment::new(140, 240);

        // find_prime() when all true
        assert_eq!(Some(143), segment.find_prime(142));
        assert_eq!(Some(143), segment.find_prime(143));
        assert_eq!(Some(149), segment.find_prime(144));
        assert_eq!(Some(209), segment.find_prime(208));
        assert_eq!(Some(209), segment.find_prime(209));
        assert_eq!(Some(211), segment.find_prime(210));
        assert_eq!(Some(211), segment.find_prime(211));
        assert_eq!(Some(221), segment.find_prime(212));
        assert_eq!(Some(239), segment.find_prime(238));
        assert_eq!(Some(239), segment.find_prime(239));
        assert_eq!(None, segment.find_prime(240));
        assert_eq!(None, segment.find_prime(241));

        // strike_prime()
        segment.strike_primes(&[11, 13]);

        // find_prime() after sieving
        assert_eq!(Some(149), segment.find_prime(142));
        assert_eq!(Some(149), segment.find_prime(143));
        assert_eq!(Some(149), segment.find_prime(144));
        assert_eq!(Some(211), segment.find_prime(208));
        assert_eq!(Some(211), segment.find_prime(209));
        assert_eq!(Some(211), segment.find_prime(210));
        assert_eq!(Some(211), segment.find_prime(211));
        assert_eq!(Some(223), segment.find_prime(212));
        assert_eq!(Some(239), segment.find_prime(238));
        assert_eq!(Some(239), segment.find_prime(239));
        assert_eq!(None, segment.find_prime(240));
        assert_eq!(None, segment.find_prime(241));
    }

    #[test]
    fn spoke_to_n_to_spoke_correct() {
        let spoke = Spoke::new(13, 20, 50);
        assert_eq!(5, spoke.n_to_spoke(spoke.spoke_to_n(5)));
        assert_eq!(6, spoke.n_to_spoke(spoke.spoke_to_n(6)));
        assert_eq!(7, spoke.n_to_spoke(spoke.spoke_to_n(7)));
        assert_eq!(8, spoke.n_to_spoke(spoke.spoke_to_n(8)));
        assert_eq!(9, spoke.n_to_spoke(spoke.spoke_to_n(9)));
    }

    #[test]
    fn bit_vec_correct() {
        let mut bit_vec = BitVec::new(12);

        bit_vec.unset(2);
        bit_vec.unset(3);
        bit_vec.unset(5);
        bit_vec.unset(7);
        bit_vec.unset(11);
        assert_eq!(Some(0), bit_vec.find(0));
        assert_eq!(Some(1), bit_vec.find(1));
        assert_eq!(Some(4), bit_vec.find(2));
        assert_eq!(Some(4), bit_vec.find(3));
        assert_eq!(Some(4), bit_vec.find(4));
        assert_eq!(Some(6), bit_vec.find(5));
        assert_eq!(Some(6), bit_vec.find(6));
        assert_eq!(Some(8), bit_vec.find(7));
        assert_eq!(Some(8), bit_vec.find(8));
        assert_eq!(Some(9), bit_vec.find(9));
        assert_eq!(Some(10), bit_vec.find(10));
        assert_eq!(None, bit_vec.find(11));
    }

    #[test]
    fn bit_vec_bounds() {
        let bit_vec = BitVec::new(12);

        // Past the end. BitVec does no bounds checking so rustc can inline its methods
        assert_eq!(None, bit_vec.find(12));
        assert_eq!(None, bit_vec.find(13));
        assert_eq!(None, bit_vec.find(14));
        assert_eq!(None, bit_vec.find(15));

        let bit_vec = BitVec::new(64);

        // Assure the last byte is correct even when it's a multiple of BitVec::WORD_BITS
        assert_eq!(Some(56), bit_vec.find(56));
        assert_eq!(Some(57), bit_vec.find(57));
        assert_eq!(Some(58), bit_vec.find(58));
        assert_eq!(Some(59), bit_vec.find(59));
        assert_eq!(Some(60), bit_vec.find(60));
        assert_eq!(Some(61), bit_vec.find(61));
        assert_eq!(Some(62), bit_vec.find(62));
        assert_eq!(Some(63), bit_vec.find(63));
    }

    #[test]
    fn bit_vec_find_first_set() {
        assert_eq!(None, BitVec::find_first_set(0b00000000));
        assert_eq!(Some(0), BitVec::find_first_set(0b00000001));
        assert_eq!(Some(1), BitVec::find_first_set(0b00000010));
        assert_eq!(Some(2), BitVec::find_first_set(0b00000100));
        assert_eq!(Some(3), BitVec::find_first_set(0b00001000));
        assert_eq!(Some(4), BitVec::find_first_set(0b00010000));
        assert_eq!(Some(5), BitVec::find_first_set(0b00100000));
        assert_eq!(Some(6), BitVec::find_first_set(0b01000000));
        assert_eq!(Some(7), BitVec::find_first_set(0b10000000));

        assert_eq!(Some(3), BitVec::find_first_set(0b10101000));
        assert_eq!(Some(4), BitVec::find_first_set(0b01010000));
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
            let mut state = SieveStateMachine::new(1_000_000, 0);
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
        let mut state = SieveStateMachine::new(10_usize.pow(10), 9_876_543_210);
        state.step();
        state.step();
        if let SieveStateMachine::Wheel(mut wheel) = state {
            b.iter(|| wheel.sieve_segment())
        }
    }
}
