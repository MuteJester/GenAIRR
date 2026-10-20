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
    /// Attach a run_from continuation to its actual ancestor history. The
    /// shared boundary revision is retained once; addressed choices and
    /// committed event spans are rebased into the combined run.
    pub fn with_continuation(&self, continuation: &Outcome) -> Result<Self, &'static str> {
        let Some(initial) = continuation.revisions.first() else {
            return Err("continuation has no initial simulation");
        };
        let parent = self.final_simulation();
        if parent.pool.as_slice() != initial.pool.as_slice()
            || parent.sequence.regions != initial.sequence.regions
            || parent.assignments != initial.assignments
            || parent.live_call_support != initial.live_call_support
            || parent.coding_end_provenance != initial.coding_end_provenance
            || parent.duplex_provenance != initial.duplex_provenance
            || parent.mutation_count != initial.mutation_count
            || parent.mutation_burden != initial.mutation_burden
        {
            return Err("continuation initial state does not match the parent final state");
        }
        let pass_offset = self.pass_names.len();
        let trace_offset = self.trace.len();
        let mut revisions = self.revisions.clone();
        revisions.extend(continuation.revisions.iter().skip(1).cloned());
        let mut pass_names = self.pass_names.clone();
        pass_names.extend(continuation.pass_names.iter().cloned());
        let mut trace = self.trace.clone();
        trace.append_delta(continuation.trace.clone());
        let mut events = self.events.clone();
        events.extend(continuation.events.iter().cloned().map(|mut event| {
            event.pass_index += pass_offset;
            event.trace_span.start += trace_offset;
            event.trace_span.end += trace_offset;
            event
        }));
        Ok(Self {
            revisions,
            pass_names,
            trace,
            events,
        })
    }

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
