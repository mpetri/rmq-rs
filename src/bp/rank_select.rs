use super::BitVec64;
use super::BpBitVec;

pub const BLOCK_SIZE: usize = 8; // in 64bit words
const SELECT_ONES_PER_HINT: usize = 64 * BLOCK_SIZE * 2; // must be > block_size * 64
const SELECT_ZEROS_PER_HINT: usize = SELECT_ONES_PER_HINT;

impl BpBitVec {
    pub fn num_ones(&self) -> u64 {
        self.block_rank_pairs.last().unwrap() - 2
    }

    pub fn rank1(&self, offset: u64) -> u64 {
        if offset as usize == self.len() {
            return self.num_ones();
        }
        let sub_block = offset / 64;
        let mut r = self.sub_block_rank(sub_block);
        let sub_left = offset % 64;
        if sub_left != 0 {
            r += (self.bit_word(sub_block) << (64 - sub_left)).count_ones() as u64;
        }
        r
    }

    pub fn select0(&self, offset: usize) -> usize {
        let offset = offset as u64;
        let mut a = 0;
        let chunk = offset as usize / SELECT_ZEROS_PER_HINT;
        if chunk != 0 {
            a = self.select0_hints[chunk - 1];
        }
        let mut b = self.select0_hints[chunk] as u64 + 1;

        while b - a > 1 {
            let mid = a + (b - a) / 2;
            let x = self.block_rank0(mid);
            if x <= offset {
                a = mid;
            } else {
                b = mid;
            }
        }
        let block = a;

        let block_offset = block * BLOCK_SIZE as u64;
        let mut cur_rank0 = self.block_rank0(block);

        let rank_in_block_parallel =
            ((offset - cur_rank0) as u64 * crate::util::ONES_STEP_9) as u64;
        let sub_ranks = 64 * crate::util::INV_COUNT_STEP_9 - self.sub_block_ranks(block);
        let sub_block_offset = (crate::util::uleq_step_9(sub_ranks, rank_in_block_parallel)
            .wrapping_mul(crate::util::ONES_STEP_9))
            >> 54
            & 0x7;

        cur_rank0 += sub_ranks >> ((7 - sub_block_offset) * 9) & 0x1FF;

        let word_offset = block_offset + sub_block_offset;
        let ret = word_offset * 64
            + crate::util::select1_in_u64(!self.bit_word(word_offset), offset - cur_rank0);
        ret as usize
    }

    fn block_rank1(&self, block: u64) -> u64 {
        self.block_rank_pairs[block as usize * 2] as u64
    }

    fn block_rank0(&self, block: u64) -> u64 {
        block * BLOCK_SIZE as u64 * 64 - self.block_rank_pairs[block as usize * 2] as u64
    }

    fn sub_block_ranks(&self, sub_block: u64) -> u64 {
        self.block_rank_pairs[sub_block as usize * 2 + 1] as u64
    }

    pub(crate) fn sub_block_rank(&self, sub_block: u64) -> u64 {
        let mut r: u64 = 0;
        let block = sub_block / BLOCK_SIZE as u64;
        r += self.block_rank1(block);
        let left = sub_block % BLOCK_SIZE as u64;
        r += self.sub_block_ranks(block) >> ((7 - left) * 9) & 0x1FF;
        r
    }
}

pub fn build_rank_pairs(bitvec: &BitVec64) -> Vec<u64> {
    let mut block_rank_pairs = Vec::new();
    let mut next_rank = 0u64;
    let mut cur_subrank = 0u64;
    let mut subranks = 0u64;
    block_rank_pairs.push(0);
    let words = bitvec.as_raw_slice();
    for (i, &word) in words.iter().enumerate() {
        let word_pop = word.count_ones();
        let shift = i % BLOCK_SIZE;
        if shift != 0 {
            subranks <<= 9;
            subranks |= cur_subrank;
        }
        next_rank += word_pop as u64;
        cur_subrank += word_pop as u64;

        if shift == BLOCK_SIZE - 1 {
            block_rank_pairs.push(subranks);
            block_rank_pairs.push(next_rank);
            subranks = 0;
            cur_subrank = 0;
        }
    }
    let left = BLOCK_SIZE - words.len() % BLOCK_SIZE;
    for _i in 0..left {
        subranks <<= 9;
        subranks |= cur_subrank;
    }
    block_rank_pairs.push(subranks);

    if words.len() % BLOCK_SIZE != 0 {
        block_rank_pairs.push(next_rank);
        block_rank_pairs.push(0);
    }
    block_rank_pairs
}

pub fn build_select0_hints(block_rank_pairs: &[u64]) -> Vec<u64> {
    let num_blocks = block_rank_pairs.len() / 2 - 1;
    let mut select0_hints = Vec::new();
    let mut cur_zeros_threshold = SELECT_ZEROS_PER_HINT as u64;
    for i in 0..num_blocks {
        let block = i as u64 + 1;
        let block_rank0 =
            block * BLOCK_SIZE as u64 * 64 - block_rank_pairs[block as usize * 2] as u64;
        if block_rank0 > cur_zeros_threshold {
            select0_hints.push(i as u64);
            cur_zeros_threshold += SELECT_ZEROS_PER_HINT as u64;
        }
    }
    select0_hints.push(num_blocks as u64);
    select0_hints
}

#[cfg(test)]
mod tests {
    use proptest::prelude::*;

    use crate::bp::BitVec64;

    proptest! {
        #[test]
        fn select0(
            zero_positions in prop::collection::btree_set(any::<u16>(), 1..1000)
        ) {

            let mut bv = BitVec64::new();
            bv.resize(u16::MAX  as usize + 1, true);

            for one_pos in &zero_positions {
                bv.set(*one_pos as usize, false);
            }

            let bp_vec = super::BpBitVec::from_bitvec(bv);

            for (offset,pos) in zero_positions.into_iter().enumerate() {
                let actual = bp_vec.select0(offset);
                assert_eq!(pos as usize,actual);
            }

        }
    }

    proptest! {
        #[test]
        fn rank1(
            one_positions in prop::collection::btree_set(any::<u16>(), 1..1000)
        ) {
            let mut bv = BitVec64::new();
            bv.resize(u16::MAX  as usize + 1, false);

            for one_pos in &one_positions {
                bv.set(*one_pos as usize, true);
            }

            let bp_vec = super::BpBitVec::from_bitvec(bv);

            for (offset,pos) in one_positions.into_iter().enumerate() {
                let actual = bp_vec.rank1(pos as u64);
                assert_eq!(actual as usize,offset);
            }
        }
    }
}
