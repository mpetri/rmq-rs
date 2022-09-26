use bitvec::prelude::*;
mod excess_tables;
mod rank_select;

use self::excess_tables::{FWD_EXC, FWD_MIN, FWD_MIN_IDX};

pub type BitVec64 = BitVec<u64, Lsb0>;

const BP_BLOCK_SIZE: usize = 4;
const BP_SUPERBLOCK_SIZE: usize = 32;

#[derive(Debug)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct BpBitVec {
    bv: BitVec64,
    select0_hints: Vec<u64>,
    block_rank_pairs: Vec<u64>,
    internal_nodes: u64,
    block_excess_min: Vec<i16>,
    superblock_excess_min: Vec<isize>,
}

impl BpBitVec {
    pub fn from_bitvec(bv: BitVec64) -> Self {
        let block_rank_pairs = rank_select::build_rank_pairs(&bv);
        let select0_hints = rank_select::build_select0_hints(&block_rank_pairs);
        let (internal_nodes, block_excess_min, superblock_excess_min) =
            build_min_tree(&bv, &block_rank_pairs);
        Self {
            bv,
            select0_hints,
            block_rank_pairs,
            internal_nodes,
            block_excess_min,
            superblock_excess_min,
        }
    }

    pub fn excess(&self, offset: usize) -> isize {
        2 * self.rank1(offset as u64) as isize - offset as isize
    }

    pub fn excess_rmq(&self, range: std::ops::RangeInclusive<usize>) -> (usize, isize) {
        let mut cur_excess = self.excess(*range.start());
        let mut min_excess = cur_excess;
        let mut min_excess_idx = *range.start();
        if range.start() == range.end() {
            return (min_excess_idx, min_excess);
        }

        let range_len = range.end() - range.start();
        let word_start_idx = range.start() / 64;
        let word_b_idx = (range.end() - 1) / 64;

        // search in word_a
        let shift_start = *range.start() % 64;
        let shifted_word_start = self.bit_word(word_start_idx as u64) >> shift_start;
        let subword_len_start = range_len.min(64 - shift_start);

        let padded_word_start = if subword_len_start == 64 {
            shifted_word_start
        } else {
            shifted_word_start | (!0 << subword_len_start)
        };

        excess_rmq_in_word(
            padded_word_start,
            *range.start() as u64,
            &mut cur_excess,
            &mut min_excess,
            &mut min_excess_idx,
        );

        if word_start_idx == word_b_idx {
            return (min_excess_idx, min_excess);
        }

        let block_start = word_start_idx / BP_BLOCK_SIZE;
        let block_b = word_b_idx / BP_BLOCK_SIZE;

        cur_excess -= 64 - subword_len_start as isize; // remove padding

        if block_start == block_b {
            // same block
            self.excess_rmq_in_block(
                word_start_idx + 1,
                word_b_idx,
                &mut cur_excess,
                &mut min_excess,
                &mut min_excess_idx,
            );
        } else {
            // search in partial block of word_a
            self.excess_rmq_in_block(
                word_start_idx + 1,
                (block_start + 1) * BP_BLOCK_SIZE,
                &mut cur_excess,
                &mut min_excess,
                &mut min_excess_idx,
            );

            // search in blocks
            let mut block_min_excess = min_excess;
            let mut block_min_idx = usize::MAX;

            let superblock_start = (block_start + 1) / BP_SUPERBLOCK_SIZE;
            let superblock_b = block_b / BP_SUPERBLOCK_SIZE;

            if superblock_start == superblock_b {
                // same superblock
                self.excess_rmq_in_superblock(
                    block_start as u64 + 1,
                    block_b as u64,
                    &mut block_min_excess,
                    &mut block_min_idx,
                );
            } else {
                // partial superblock of a
                self.excess_rmq_in_superblock(
                    block_start as u64 + 1,
                    (superblock_start as u64 + 1) * BP_SUPERBLOCK_SIZE as u64,
                    &mut block_min_excess,
                    &mut block_min_idx,
                );

                // search min superblock in the min tree
                let mut superblock_min_excess = min_excess;
                let mut superblock_min_idx = usize::MAX;
                self.find_min_superblock(
                    superblock_start as u64 + 1,
                    superblock_b as u64,
                    &mut superblock_min_excess,
                    &mut superblock_min_idx,
                );

                if superblock_min_excess < min_excess {
                    self.excess_rmq_in_superblock(
                        superblock_min_idx as u64 * BP_SUPERBLOCK_SIZE as u64,
                        (superblock_min_idx as u64 + 1) * BP_SUPERBLOCK_SIZE as u64,
                        &mut block_min_excess,
                        &mut block_min_idx,
                    );
                }

                // partial superblock of b
                self.excess_rmq_in_superblock(
                    superblock_b as u64 * BP_SUPERBLOCK_SIZE as u64,
                    block_b as u64,
                    &mut block_min_excess,
                    &mut block_min_idx,
                );
            }

            if block_min_excess < min_excess {
                cur_excess = self.get_block_excess(block_min_idx as u64);
                self.excess_rmq_in_block(
                    block_min_idx * BP_BLOCK_SIZE,
                    (block_min_idx + 1) * BP_BLOCK_SIZE,
                    &mut cur_excess,
                    &mut min_excess,
                    &mut min_excess_idx,
                );
            }

            // search in partial block of word_b
            cur_excess = self.get_block_excess(block_b as u64);
            self.excess_rmq_in_block(
                block_b * BP_BLOCK_SIZE,
                word_b_idx,
                &mut cur_excess,
                &mut min_excess,
                &mut min_excess_idx,
            );
        }

        // search in word_b
        let word_b = self.bit_word(word_b_idx as u64);
        let offset_b = range.end() % 64;
        let padded_word_b = if offset_b == 0 {
            word_b
        } else {
            word_b | (!0 << offset_b)
        };

        excess_rmq_in_word(
            padded_word_b,
            word_b_idx as u64 * 64,
            &mut cur_excess,
            &mut min_excess,
            &mut min_excess_idx,
        );

        (min_excess_idx, min_excess)
    }

    pub fn len(&self) -> usize {
        self.bv.len()
    }
}

impl BpBitVec {
    fn bit_word(&self, offset: u64) -> u64 {
        self.bv.as_raw_slice()[offset as usize]
    }

    fn excess_rmq_in_block(
        &self,
        start: usize,
        end: usize,
        excess: &mut isize,
        min_excess: &mut isize,
        min_excess_idx: &mut usize,
    ) {
        for idx in start..end {
            excess_rmq_in_word(
                self.bit_word(idx as u64),
                idx as u64 * 64,
                excess,
                min_excess,
                min_excess_idx,
            );
        }
    }

    fn excess_rmq_in_superblock(
        &self,
        block_start: u64,
        block_end: u64,
        block_min_excess: &mut isize,
        block_min_idx: &mut usize,
    ) {
        if block_start != block_end {
            let superblock = block_start / BP_SUPERBLOCK_SIZE as u64;

            let superblock_excess = self.get_block_excess(superblock * BP_SUPERBLOCK_SIZE as u64);

            for block in block_start..block_end {
                if superblock_excess + (self.block_excess_min[block as usize] as isize)
                    < *block_min_excess
                {
                    *block_min_excess =
                        superblock_excess + self.block_excess_min[block as usize] as isize;
                    *block_min_idx = block as usize;
                }
            }
        }
    }

    fn find_min_superblock(
        &self,
        superblock_start: u64,
        superblock_end: u64,
        superblock_min_excess: &mut isize,
        superblock_min_idx: &mut usize,
    ) {
        if superblock_start == superblock_end {
            return;
        }
        let mut cur_node = self.internal_nodes + superblock_start;
        let mut rightmost_span = superblock_start;

        let mut node_min_excess = self.superblock_excess_min[cur_node as usize];
        let mut node_min_idx = cur_node;

        if superblock_end - superblock_start == 1 {
            *superblock_min_excess = node_min_excess;
            *superblock_min_idx = superblock_start as usize;
            return;
        }

        let mut h = 0;
        loop {
            if (cur_node & 1) == 0 {
                // is a left child
                let right_sibling = cur_node + 1;
                rightmost_span += 1u64 << h;

                if rightmost_span < superblock_end
                    && self.superblock_excess_min[right_sibling as usize] < node_min_excess
                {
                    node_min_excess = self.superblock_excess_min[right_sibling as usize];
                    node_min_idx = right_sibling;
                }

                if rightmost_span >= superblock_end - 1 {
                    cur_node += 1;
                    break;
                }
            }

            cur_node /= 2; // parent
            h += 1;
        }

        while rightmost_span > superblock_end - 1 {
            h -= 1;

            let left_child = cur_node * 2;
            let right_child_span = 1 << h;
            if (rightmost_span - right_child_span) >= (superblock_end - 1) {
                // go to left child
                rightmost_span -= right_child_span;
                cur_node = left_child;
            } else {
                // go to right child and add left subtree to candidate
                // subblocks
                if self.superblock_excess_min[left_child as usize] < node_min_excess {
                    node_min_excess = self.superblock_excess_min[left_child as usize];
                    node_min_idx = left_child;
                }
                cur_node = left_child + 1;
            }
        }

        // check last left-turn
        if rightmost_span < superblock_end
            && self.superblock_excess_min[cur_node as usize] < node_min_excess
        {
            node_min_excess = self.superblock_excess_min[cur_node as usize];
            node_min_idx = cur_node;
        }

        if node_min_excess < *superblock_min_excess {
            cur_node = node_min_idx;
            while cur_node < self.internal_nodes {
                cur_node *= 2;
                // remember that past-the-end nodes are filled with size()
                if self.superblock_excess_min[cur_node as usize + 1]
                    < self.superblock_excess_min[cur_node as usize]
                {
                    cur_node += 1;
                }
            }

            *superblock_min_excess = node_min_excess;
            *superblock_min_idx = cur_node as usize - self.internal_nodes as usize;
        }
    }

    fn get_block_excess(&self, block: u64) -> isize {
        let sub_block_idx = block * BP_BLOCK_SIZE as u64;
        let block_pos = sub_block_idx * 64;
        (2 * self.sub_block_rank(sub_block_idx) - block_pos) as isize
    }
}

fn build_min_tree(bitvec: &BitVec64, block_rank_pairs: &[u64]) -> (u64, Vec<i16>, Vec<isize>) {
    const fn get_block_excess(block: u64, block_rank_pairs: &[u64]) -> isize {
        const fn block_rank1(block: u64, block_rank_pairs: &[u64]) -> u64 {
            block_rank_pairs[block as usize * 2] as u64
        }

        const fn sub_block_ranks(sub_block: u64, block_rank_pairs: &[u64]) -> u64 {
            block_rank_pairs[sub_block as usize * 2 + 1] as u64
        }

        const fn sub_block_rank(sub_block: u64, block_rank_pairs: &[u64]) -> u64 {
            let mut r: u64 = 0;
            let block = sub_block / rank_select::BLOCK_SIZE as u64;
            r += block_rank1(block, block_rank_pairs);
            let left = sub_block % rank_select::BLOCK_SIZE as u64;
            r += sub_block_ranks(block, block_rank_pairs) >> ((7 - left) * 9) & 0x1FF;
            r
        }
        let sub_block_idx = block * BP_BLOCK_SIZE as u64;
        let block_pos = sub_block_idx * 64;

        2 * sub_block_rank(sub_block_idx, block_rank_pairs) as isize - block_pos as isize
    }

    let mut block_excess_min = Vec::new();
    let mut cur_block_min = 0;
    let mut cur_superblock_excess = 0;
    let words = bitvec.as_raw_slice();
    for (sub_block, &word) in words.iter().enumerate() {
        if sub_block % BP_BLOCK_SIZE == 0 {
            if sub_block % (BP_BLOCK_SIZE * BP_SUPERBLOCK_SIZE) == 0 {
                cur_superblock_excess = 0;
            }

            if sub_block != 0 {
                block_excess_min.push(cur_block_min as i16);
                cur_block_min = cur_superblock_excess;
            }
        }
        let mut mask = 1;
        // for last block stop at bit boundary
        let n_bits = match sub_block == words.len() - 1 {
            true if bitvec.len() % 64 != 0 => bitvec.len() % 64,
            _ => 64,
        };
        for _i in 0..n_bits {
            cur_superblock_excess += if (word & mask) != 0 { 1 } else { -1 };
            cur_block_min = cur_block_min.min(cur_superblock_excess);
            mask <<= 1;
        }
    }

    // Flush last block mins
    block_excess_min.push(cur_block_min as i16);

    let n_blocks = words.len() / BP_BLOCK_SIZE + usize::from(words.len() % BP_BLOCK_SIZE != 0);
    let n_superblocks = (n_blocks + BP_SUPERBLOCK_SIZE - 1) / BP_SUPERBLOCK_SIZE;

    assert_eq!(n_blocks, block_excess_min.len());

    let mut n_complete_leaves = 1;
    while n_complete_leaves < n_superblocks {
        n_complete_leaves <<= 1;
    }
    // n_complete_leaves is the smallest power of 2 >= n_superblocks
    let internal_nodes = n_complete_leaves;
    let treesize = internal_nodes + n_superblocks;

    let mut superblock_excess_min = vec![0isize; treesize];

    // Fill in the leaves of the tree
    for superblock in 0..n_superblocks {
        let mut cur_super_min = bitvec.len() as isize;
        let superblock_excess = get_block_excess(
            superblock as u64 * BP_SUPERBLOCK_SIZE as u64,
            block_rank_pairs,
        );

        let start = superblock * BP_SUPERBLOCK_SIZE;
        let stop = n_blocks.min((superblock + 1) * BP_SUPERBLOCK_SIZE);

        for &cur_block_excess in block_excess_min.iter().take(stop).skip(start) {
            cur_super_min = cur_super_min.min(superblock_excess + cur_block_excess as isize);
        }

        superblock_excess_min[internal_nodes + superblock] = cur_super_min;
    }

    // fill in the internal nodes with past-the-boundary values
    // (they will also serve as sentinels in debug)
    for cur_superblock_excess_min in superblock_excess_min.iter_mut().take(internal_nodes) {
        *cur_superblock_excess_min = bitvec.len() as isize;
    }

    // Fill bottom-up the other layers: each node updates the parent
    for node in (1..treesize).rev() {
        let parent = node / 2;
        superblock_excess_min[parent] =
            superblock_excess_min[parent].min(superblock_excess_min[node]);
    }

    (
        internal_nodes as u64,
        block_excess_min,
        superblock_excess_min,
    )
}

fn excess_rmq_in_word(
    word: u64,
    word_start: u64,
    excess: &mut isize,
    min_excess: &mut isize,
    min_excess_idx: &mut usize,
) {
    let mut min_byte_idx = 0;
    let mut min_byte_excess = *min_excess;
    for i in 0..8 {
        let shift = i * 8;
        let byte = (word >> shift) & 0xFF;

        let cur_min = *excess - FWD_MIN[byte as usize] as isize;

        if cur_min < min_byte_excess {
            min_byte_idx = i;
            min_byte_excess = cur_min;
        }

        *excess += FWD_EXC[byte as usize] as isize;
    }

    if min_byte_excess < *min_excess {
        *min_excess = min_byte_excess;
        let shift = min_byte_idx * 8;
        *min_excess_idx = word_start as usize
            + shift as usize
            + FWD_MIN_IDX[((word >> shift) & 0xFF) as usize] as usize;
    }
}
