use super::Pass;

/// An ordered sequence of passes to run as a single simulation.
pub struct PassPlan {
    passes: Vec<Box<dyn Pass>>,
}

impl PassPlan {
    /// Empty plan — runs no passes, leaves the initial IR unchanged.
    pub fn new() -> Self {
        Self { passes: Vec::new() }
    }

    /// Append a pass to the end of the plan.
    pub fn push(&mut self, pass: Box<dyn Pass>) -> &mut Self {
        self.passes.push(pass);
        self
    }

    /// Number of passes in the plan.
    pub fn len(&self) -> usize {
        self.passes.len()
    }

    /// Whether the plan contains zero passes.
    pub fn is_empty(&self) -> bool {
        self.passes.is_empty()
    }

    /// Read-only view of the underlying pass vector.
    pub fn passes(&self) -> &[Box<dyn Pass>] {
        &self.passes
    }
}

impl Default for PassPlan {
    fn default() -> Self {
        Self::new()
    }
}
