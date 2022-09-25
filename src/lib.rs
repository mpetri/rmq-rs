mod bp;
mod cartesian_tree;
mod util;

use cartesian_tree::CartesianTree;

pub struct Rmq {
    cartesian_tree: CartesianTree,
}

impl Rmq {
    pub fn range_minimum(&self, range: impl std::ops::RangeBounds<usize>) -> Option<usize> {
        self.cartesian_tree.range_minimum(range)
    }

    pub fn range_minumum_iter(_range: std::ops::Range<usize>) {}
}

impl<T: Ord> FromIterator<T> for Rmq {
    fn from_iter<I: IntoIterator<Item = T>>(iter: I) -> Rmq {
        Rmq {
            cartesian_tree: CartesianTree::from_iter(iter),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use proptest::prelude::*;

    proptest! {
        #[test]
        fn rmq_works_full_range(
            elems in prop::collection::vec(any::<u16>(), 1..1000)
        ) {
            let other = elems.clone();
            let rmq = Rmq::from_iter(other);

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
            let rmq = Rmq::from_iter(other);

            let min_pos = rmq.range_minimum(start..(start+len)).unwrap();

            let min = elems.iter().skip(start).take(len).copied().min().unwrap();
            assert_eq!(elems[min_pos],min);
        }
    }
}