#![warn(missing_debug_implementations, missing_docs, rust_2018_idioms)]
//!
//! # Range Minimum Query (RMQ)
//!
//! This crate implements a succinct data structure to solve the
//! Range Minimum Query (RMQ) problem in constant time and linear space.
//!
//! This code was derived (almost verbatim) from a succinct version of RMQ, implemented by
//! Giuseppe Ottaviano's succinct c++ library: `https://github.com/ot/succinct`.
//!
//! # What is RMQ (taken from Wikipedia)
//!
//! In computer science, a range minimum query (RMQ) solves the problem of
//! finding the minimal value in a sub-array of an array of comparable objects.
//! Range minimum queries have several use cases in computer science, such as
//! the lowest common ancestor problem and the longest common prefix problem (LCP).
//!
//! For example, when `A = [0,5,2,5,4,3,1,6,3]`, then the answer to the range
//! minimum query for the sub-array `A[2...7] = [2,5,4,3,1,6]` is `6`, as `A[6] = 1`.
//!
//! # Usage
//!
//! ```rust
//! use range_minimum_query::Rmq;
//!
//! let a = [0,5,2,5,4,3,1,6,3];
//! let rmq = Rmq::from_iter(a);
//! let res = rmq.range_minimum(2..=7);
//! assert_eq!(res.unwrap(),6);
//! ```

mod bp;
mod cartesian_tree;
mod util;

use cartesian_tree::CartesianTree;

/// The main RMQ data structure
#[derive(Debug)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct Rmq {
    cartesian_tree: CartesianTree,
}

impl Rmq {
    /// returns the position of the minimum element in `range`
    pub fn range_minimum(&self, range: impl std::ops::RangeBounds<usize>) -> Option<usize> {
        self.cartesian_tree.range_minimum(range)
    }
}

impl<T: Ord> FromIterator<T> for Rmq {
    fn from_iter<I: IntoIterator<Item = T>>(iter: I) -> Self {
        Self {
            cartesian_tree: CartesianTree::from_iter(iter),
        }
    }
}

#[cfg(test)]
mod tests {

    use proptest::prelude::*;

    proptest! {
        #[test]
        fn rmq_works_full_range(
            elems in prop::collection::vec(any::<u16>(), 1..1000)
        ) {
            let other = elems.clone();
            let rmq = super::Rmq::from_iter(other.iter());

            let min_pos = rmq.range_minimum(0..elems.len()).unwrap();

            let min = elems.iter().copied().min().unwrap();
            assert_eq!(elems[min_pos],min);
        }
    }

    proptest! {
        #[test]
        fn rmq_works_partial_range(
            elems in prop::collection::vec(any::<u16>(), 1..1000),
            start in 0usize..1000,
            len in 1usize..100
        ) {
            prop_assume!(start + len < elems.len());
            let other = elems.clone();
            let rmq = super::Rmq::from_iter(other);

            let min_pos = rmq.range_minimum(start..(start+len)).unwrap();

            let min = elems.iter().skip(start).take(len).copied().min().unwrap();
            assert_eq!(elems[min_pos],min);
        }
    }

    proptest! {
        #[test]
        fn rmq_with_max(
            elems in prop::collection::vec(any::<u16>(), 1..1000),
            start in 0usize..1000,
            len in 1usize..100
        ) {
            prop_assume!(start + len < elems.len());
            let other = elems.clone();
            let rmq = super::Rmq::from_iter(other.into_iter().map(std::cmp::Reverse));

            let maxpos = rmq.range_minimum(start..(start+len)).unwrap();

            let max = elems.iter().skip(start).take(len).copied().max().unwrap();
            assert_eq!(elems[maxpos],max);
        }
    }

    proptest! {
        #[test]
        fn rmq_works_prefix_range(
            elems in prop::collection::btree_set(any::<u32>(), 1..1000)
        ) {
            let elem_vec : Vec<u32> = elems.iter().copied().collect();
            let other = elems.clone();
            let rmq = super::Rmq::from_iter(other);

            for end in 0..elems.len() {
                let min_pos = rmq.range_minimum(0..=end).unwrap();
                let min = elems.iter().take(end+1).copied().min().unwrap();
                assert_eq!(elem_vec[min_pos],min);
            }
        }
    }
}
