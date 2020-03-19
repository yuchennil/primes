use crate::constants::ceil_div;

pub struct BitVec {
    bit_vec: Vec<u64>,
    word_index: usize,
    bit_index: usize,
}

impl Iterator for BitVec {
    type Item = usize;

    #[inline]
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

    pub fn find(&self, index: usize) -> Option<usize> {
        let first_word_index = index >> BitVec::SHIFT;
        for (word_index, &word) in self.bit_vec[first_word_index..].iter().enumerate() {
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

    // #[test]
    // fn bit_vec_bounds() {
    //     let mut bit_vec = BitVec::new(12);

    //     // Past the end. BitVec does no bounds checking so rustc can inline its methods
    //     bit_vec.skip_to(12);
    //     assert_eq!(vec![0; 0], bit_vec.collect::<Vec<_>>());
    // }

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
