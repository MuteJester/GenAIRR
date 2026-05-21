/// Index into a `NucleotidePool`. Stable for the lifetime of the
/// pool revision it was issued for.
#[derive(Copy, Clone, Eq, PartialEq, Hash, Debug)]
pub struct NucHandle(u32);

impl NucHandle {
    pub const fn new(idx: u32) -> Self {
        Self(idx)
    }
    pub const fn index(self) -> u32 {
        self.0
    }
    pub const fn as_usize(self) -> usize {
        self.0 as usize
    }
}

/// Index into a `Sequence`'s regions list. Stable for the lifetime of
/// the sequence revision it was issued for.
#[derive(Copy, Clone, Eq, PartialEq, Hash, Debug)]
pub struct RegionHandle(u32);

impl RegionHandle {
    pub const fn new(idx: u32) -> Self {
        Self(idx)
    }
    pub const fn index(self) -> u32 {
        self.0
    }
    pub const fn as_usize(self) -> usize {
        self.0 as usize
    }
}
