use crate::event::EventRecord;
use crate::ir::Simulation;
use crate::trace::Trace;

/// The result of executing a `PassPlan` on an initial `Simulation`.
#[derive(Clone, Debug)]
pub struct Outcome {
    pub revisions: Vec<Simulation>,
    pub pass_names: Vec<String>,
    pub trace: Trace,
    pub events: Vec<EventRecord>,
}

impl Outcome {
    /// The final IR revision after every pass has run.
    pub fn final_simulation(&self) -> &Simulation {
        self.revisions
            .last()
            .expect("Outcome must always contain at least the initial revision")
    }

    /// First revision produced by the pass with the given name.
    pub fn revision_after(&self, name: &str) -> Option<&Simulation> {
        self.pass_names
            .iter()
            .position(|n| n == name)
            .map(|i| &self.revisions[i + 1])
    }

    /// The committed event ledger for this run.
    pub fn events(&self) -> &[EventRecord] {
        &self.events
    }
}
