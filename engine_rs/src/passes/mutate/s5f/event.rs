#[derive(Clone, Copy, Debug)]
pub(super) struct S5FMutationEvent {
    pub(super) site: u32,
    pub(super) base: u8,
}

#[derive(Clone, Copy, Debug)]
pub(super) struct WeightedS5FMutationEvent {
    pub(super) event: S5FMutationEvent,
    pub(super) weight: f64,
}
