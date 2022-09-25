use crate::bp::BitVec64;
use crate::bp::BpBitVec;

pub(crate) struct CartesianTree {
    bp: BpBitVec,
}

impl CartesianTree {
    fn builder<T: Ord>() -> CartesianTreeBuilder<T> {
        CartesianTreeBuilder::new()
    }

    fn len(&self) -> usize {
        self.bp.len() / 2 - 1
    }

    pub fn range_minimum(&self, range: impl std::ops::RangeBounds<usize>) -> Option<usize> {
        let range_start = match range.start_bound() {
            std::ops::Bound::Included(t) => *t,
            std::ops::Bound::Excluded(t) => *t + 1,
            std::ops::Bound::Unbounded => 0,
        };
        let range_end = match range.end_bound() {
            std::ops::Bound::Included(t) => *t,
            std::ops::Bound::Excluded(t) => *t - 1,
            std::ops::Bound::Unbounded => self.len() - 1,
        };
        let range = range_start..=range_end;
        dbg!(&range);

        if range.is_empty() {
            return None;
        }
        if range.start() == range.end() {
            return Some(*range.start());
        }

        let n = self.len();

        let t = self.bp.select0(n - range.end() - 1);
        let exc_t = (t - 2 * (n - range.end() - 1)) as isize;

        assert!(exc_t - 1 == self.bp.excess(t + 1));

        let x = self.bp.select0(n - range.end());
        let y = self.bp.select0(n - range.start());

        let (w, exc_w) = self.bp.excess_rmq(x..=y);

        dbg!(&w);
        dbg!(&exc_w);

        let rank0_w = (w - exc_w as usize) / 2;

        dbg!(&rank0_w);

        if exc_w >= exc_t - 1 {
            Some(*range.end())
        } else {
            Some(n - rank0_w)
        }
    }
}

impl<T: Ord> FromIterator<T> for CartesianTree {
    fn from_iter<I: IntoIterator<Item = T>>(iter: I) -> CartesianTree {
        let mut builder = CartesianTree::builder();
        for item in iter {
            builder.push(item);
        }
        builder.build()
    }
}

struct CartesianTreeBuilder<T> {
    bp: BitVec64,
    stack: Vec<T>,
}

impl<T: Ord> CartesianTreeBuilder<T> {
    pub fn new() -> CartesianTreeBuilder<T> {
        CartesianTreeBuilder {
            bp: BitVec64::new(),
            stack: Vec::new(),
        }
    }

    pub fn push(&mut self, elem: T) {
        self.bp.push(false);
        while matches!(self.stack.last(), Some(x) if elem.cmp(x) == std::cmp::Ordering::Less) {
            self.stack.pop();
            self.bp.push(true);
        }
        self.stack.push(elem);
    }

    pub fn build(mut self) -> CartesianTree {
        // super-root
        self.bp.push(false);
        while self.stack.pop().is_some() {
            self.bp.push(true);
        }
        self.bp.push(true);

        self.bp.reverse();
        CartesianTree {
            bp: BpBitVec::from_bitvec(self.bp),
        }
    }
}
