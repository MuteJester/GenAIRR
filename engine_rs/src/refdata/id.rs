/// Index into an `AllelePool`. Stable for the lifetime of the
/// `RefDataConfig` it was issued for.
///
/// Distinct type from `NucHandle` and `RegionHandle` so the compiler
/// catches handle-confusion at call sites — this is the same
/// discipline as the IR handles in §3 of the design doc.
#[derive(Copy, Clone, Eq, PartialEq, Hash, Debug)]
pub struct AlleleId(u32);

impl AlleleId {
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
