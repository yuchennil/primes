use crate::constants::ceil_div;

/// Vector of boolean values optimized for prime sieving
///
/// Each bool in Rust is represented by a full 8-bit byte. For memory-constrained applications like
/// prime sieving, it is eight times as memory-efficient to compress these values into a single byte
/// and provide wrapper methods for access.
///
/// As it turns out, the only methods needed by a sieve are
/// - create a new BitVec of trues
/// - set a single index's bit to false
/// - iterate over true bits in order
///
/// We implement a number of optimizations to improve performance:
/// - use u64s as the backing data type. A modern 64-bit processor is just as fast at accessing and
///   writing to a single 64-bit word as it is with a 32-bit or 8-bit subword
/// - check u64 equality with 0 when iterating through words. Again, it's just as fast to check a
///   whole word as it is to check a subword part (or even a bit).
/// - use the builtin u64::trailing_zeros() method to calculate the first set bit in a u64. This
///   is much faster than a naive search with bit masks, or even a custom bit hack using de Bruijn
///   sequences.
/// - use lookup tables to replace bitshift operations. If the UNSET_BITS and GREATER_OR_EQUAL_BITS
///   tables are already in L1 cache, then they're faster to query than a variable bitshift.
pub struct BitVec {
    bit_vec: Vec<u64>,
    word_index: usize,
    bit_index: usize,
}

impl Iterator for BitVec {
    type Item = usize;

    fn next(&mut self) -> Option<Self::Item> {
        loop {
            let word = self.bit_vec.get(self.word_index)?;
            let masked_word = word & BitVec::GREATER_OR_EQUAL_BITS[self.bit_index];
            match BitVec::find_first_set(masked_word) {
                Some(bit_index) => {
                    let result = (self.word_index << BitVec::SHIFT) + bit_index;
                    if bit_index + 1 == BitVec::WORD_BITS {
                        self.word_index += 1;
                        self.bit_index = 0;
                    } else {
                        self.bit_index = bit_index + 1;
                    }
                    return Some(result);
                }
                None => {
                    self.word_index += 1;
                    self.bit_index = 0;
                }
            }
        }
    }
}

impl BitVec {
    const WORD_BITS: usize = 64;
    const SHIFT: usize = 6;
    const MASK: usize = 0b11_1111;
    const ONES: u64 = std::u64::MAX;
    const UNSET_BIT: [u64; 64] = unset_bit();
    const GREATER_OR_EQUAL_BITS: [u64; 64] = greater_or_equal_bits();

    pub fn new(len: usize) -> BitVec {
        let mut bit_vec = vec![BitVec::ONES; ceil_div(len, BitVec::WORD_BITS)];
        if let Some(end) = bit_vec.get_mut(len >> BitVec::SHIFT) {
            *end &= !BitVec::GREATER_OR_EQUAL_BITS[len & BitVec::MASK];
        }
        let word_index = 0;
        let bit_index = 0;

        BitVec {
            bit_vec,
            word_index,
            bit_index,
        }
    }

    pub fn unset(&mut self, index: usize) {
        self.bit_vec[index >> BitVec::SHIFT] &= BitVec::UNSET_BIT[index & BitVec::MASK]
    }

    /// Find the first set bit in word. This index is equal to the number of word's trailing zeros.
    fn find_first_set(word: u64) -> Option<usize> {
        if word == 0 {
            return None;
        }
        Some(word.trailing_zeros() as usize)
    }
}

const fn unset_bit() -> [u64; 64] {
    let mut unset_bit = [0; 64];
    let mut i = 0;
    while i < 64 {
        unset_bit[i] = !(1 << i);
        i += 1;
    }
    unset_bit
}

const fn greater_or_equal_bits() -> [u64; 64] {
    let mut greater_or_equal_bits = [0; 64];
    let mut i = 0;
    while i < 64 {
        greater_or_equal_bits[i] = BitVec::ONES << i;
        i += 1;
    }
    greater_or_equal_bits
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn bit_vec_correct() {
        let mut bit_vec = BitVec::new(12);

        bit_vec.unset(2);
        bit_vec.unset(3);
        bit_vec.unset(5);
        bit_vec.unset(7);
        bit_vec.unset(11);

        assert_eq!(vec![0, 1, 4, 6, 8, 9, 10], bit_vec.collect::<Vec<_>>());
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
}
