use crate::refdata::AlleleId;

/// Dense allele-id bitset used by live call hypotheses.
///
/// The bitset has a fixed universe size so equality is unambiguous:
/// two bitsets over different allele pools are not equal even if their
/// word storage currently contains the same bits.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct AlleleBitSet {
    universe_len: usize,
    words: Vec<u64>,
}

impl AlleleBitSet {
    pub fn empty(universe_len: usize) -> Self {
        Self {
            universe_len,
            words: vec![0; universe_len.div_ceil(64)],
        }
    }

    pub fn full(universe_len: usize) -> Self {
        let mut set = Self {
            universe_len,
            words: vec![u64::MAX; universe_len.div_ceil(64)],
        };
        set.clear_unused_tail_bits();
        set
    }

    pub fn from_ids(universe_len: usize, ids: impl IntoIterator<Item = AlleleId>) -> Self {
        let mut set = Self::empty(universe_len);
        for id in ids {
            set.insert(id);
        }
        set
    }

    pub fn universe_len(&self) -> usize {
        self.universe_len
    }

    pub fn is_empty(&self) -> bool {
        self.words.iter().all(|word| *word == 0)
    }

    pub fn len(&self) -> usize {
        self.words
            .iter()
            .map(|word| word.count_ones() as usize)
            .sum()
    }

    pub fn contains(&self, id: AlleleId) -> bool {
        let index = self.checked_index(id);
        (self.words[index / 64] & (1u64 << (index % 64))) != 0
    }

    pub fn insert(&mut self, id: AlleleId) -> bool {
        let index = self.checked_index(id);
        let word = &mut self.words[index / 64];
        let mask = 1u64 << (index % 64);
        let was_present = (*word & mask) != 0;
        *word |= mask;
        !was_present
    }

    pub fn remove(&mut self, id: AlleleId) -> bool {
        let index = self.checked_index(id);
        let word = &mut self.words[index / 64];
        let mask = 1u64 << (index % 64);
        let was_present = (*word & mask) != 0;
        *word &= !mask;
        was_present
    }

    pub fn union_with(&mut self, other: &Self) {
        self.assert_same_universe(other);
        for (lhs, rhs) in self.words.iter_mut().zip(&other.words) {
            *lhs |= *rhs;
        }
    }

    pub fn intersect_with(&mut self, other: &Self) {
        self.assert_same_universe(other);
        for (lhs, rhs) in self.words.iter_mut().zip(&other.words) {
            *lhs &= *rhs;
        }
    }

    pub fn unioned(&self, other: &Self) -> Self {
        let mut out = self.clone();
        out.union_with(other);
        out
    }

    pub fn intersected(&self, other: &Self) -> Self {
        let mut out = self.clone();
        out.intersect_with(other);
        out
    }

    pub fn iter_ids(&self) -> impl Iterator<Item = AlleleId> + '_ {
        // Walk only the set bits via `trailing_zeros` + clear-lowest-bit
        // (`w & w-1`), skipping all-zero words entirely. The bounds
        // check that was here previously is unneeded because tail bits
        // (positions >= universe_len within the final word) are
        // maintained at zero by every constructor/operator on this
        // type; see `clear_unused_tail_bits`, and `checked_index` in
        // `insert`/`remove`/`contains` which rejects out-of-universe
        // ids before any word write. `union_with` / `intersect_with`
        // preserve the invariant because OR/AND of clean tails stays
        // clean.
        self.words
            .iter()
            .copied()
            .enumerate()
            .flat_map(|(word_index, mut word)| {
                std::iter::from_fn(move || {
                    if word == 0 {
                        return None;
                    }
                    let bit = word.trailing_zeros() as usize;
                    word &= word - 1;
                    Some(AlleleId::new((word_index * 64 + bit) as u32))
                })
            })
    }

    /// Visit each set allele id, in ascending order, without iterator
    /// adapter overhead.
    ///
    /// Slightly faster than [`iter_ids`] for hot inner loops because
    /// the closure inlines cleanly through a single `for` over the
    /// underlying words; iterator-based iteration goes through a
    /// `FlatMap<...>` state machine that LLVM has a harder time fusing.
    pub fn for_each_id<F: FnMut(AlleleId)>(&self, mut f: F) {
        for (word_index, &word) in self.words.iter().enumerate() {
            let mut w = word;
            while w != 0 {
                let bit = w.trailing_zeros() as usize;
                w &= w - 1;
                f(AlleleId::new((word_index * 64 + bit) as u32));
            }
        }
    }

    pub fn to_ids(&self) -> Vec<AlleleId> {
        self.iter_ids().collect()
    }

    fn checked_index(&self, id: AlleleId) -> usize {
        let index = id.as_usize();
        assert!(
            index < self.universe_len,
            "AlleleBitSet: allele id {} outside universe length {}",
            id.index(),
            self.universe_len
        );
        index
    }

    fn assert_same_universe(&self, other: &Self) {
        assert_eq!(
            self.universe_len, other.universe_len,
            "AlleleBitSet operation requires matching universe lengths"
        );
    }

    fn clear_unused_tail_bits(&mut self) {
        let unused = self.words.len() * 64 - self.universe_len;
        if unused == 0 || self.words.is_empty() {
            return;
        }
        let used = 64 - unused;
        let mask = if used == 64 {
            u64::MAX
        } else {
            (1u64 << used) - 1
        };
        if let Some(last) = self.words.last_mut() {
            *last &= mask;
        }
    }
}
