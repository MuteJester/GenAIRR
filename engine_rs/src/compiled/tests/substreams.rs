//! Does a pass's position in the plan decide what it draws?
//!
//! The unit tests under `rng::substream` and `pass::schedule` cover
//! the derivation and the labelling. These cover the run loop: that it
//! hands each pass the stream its label names, and so that a pass's
//! draws survive another pass being inserted in front of it.
//!
//! The pass under test is a stub that draws words and changes nothing,
//! which is the only way to separate the two things a real pass does
//! at once. A `pcr_amplify` in front of a `mutate` moves the molecule
//! as well as the stream position, so a difference downstream would
//! not say which cause produced it.

use super::*;
use std::sync::{Arc, Mutex};

use crate::pass::{PassContext, PassEffect, PassRequirement};
use crate::rng::substream::{StreamRoot, SubstreamLabel};

/// A pass that draws `words` words, records them where the test can
/// read them, and leaves the simulation exactly as it found it.
struct DrawsAndChangesNothing {
    name: String,
    words: u32,
    drawn: Arc<Mutex<Vec<u64>>>,
}

impl DrawsAndChangesNothing {
    fn new(name: &str, words: u32) -> (Self, Arc<Mutex<Vec<u64>>>) {
        let drawn = Arc::new(Mutex::new(Vec::new()));
        (
            Self {
                name: name.to_string(),
                words,
                drawn: Arc::clone(&drawn),
            },
            drawn,
        )
    }
}

impl Pass for DrawsAndChangesNothing {
    fn name(&self) -> &str {
        &self.name
    }

    fn execute(&self, sim: &Simulation, ctx: &mut PassContext) -> Simulation {
        let mut log = self.drawn.lock().expect("test log is not poisoned");
        for _ in 0..self.words {
            log.push(ctx.rng.next_u64());
        }
        sim.clone()
    }

    fn requirements(&self) -> Vec<PassRequirement> {
        Vec::new()
    }

    fn effects(&self) -> Vec<PassEffect> {
        Vec::new()
    }
}

/// A VJ plan whose NP region has both a length and a base to draw, so
/// the molecule is sensitive to whatever stream the NP pass is given.
///
/// `stub` is inserted **before** the NP pass, which is the position
/// that matters: a pass appended at the end has nothing after it to
/// disturb, so a plan built by `push`ing the stub last would pass this
/// file's first test without the engine doing anything at all.
fn drawing_vj_plan(cfg: &RefDataConfig, stub: Option<Box<dyn Pass>>) -> PassPlan {
    let mut plan = PassPlan::new();
    plan.push(Box::new(SampleAllelePass::new(
        Segment::V,
        Box::new(AllelePoolDist::uniform(&cfg.v_pool)),
    )));
    plan.push(Box::new(SampleAllelePass::new(
        Segment::J,
        Box::new(AllelePoolDist::uniform(&cfg.j_pool)),
    )));
    plan.push(Box::new(AssembleSegmentPass::new(Segment::V)));
    if let Some(stub) = stub {
        plan.push(stub);
    }
    plan.push(Box::new(GenerateNPPass::new(
        Segment::Np1,
        Box::new(EmpiricalLengthDist::from_pairs(vec![(2, 0.5), (3, 0.5)])),
        Box::new(UniformBase),
    )));
    plan.push(Box::new(AssembleSegmentPass::new(Segment::J)));
    plan
}

fn molecule(plan: &PassPlan, cfg: &RefDataConfig, seed: u64) -> Vec<Nucleotide> {
    let simulator = compile_test_fixture(plan, Some(cfg), None, ExecutionPolicy::Permissive)
        .expect("fixture compiles");
    simulator
        .run_one(seed)
        .expect("fixture runs")
        .final_simulation()
        .pool
        .as_slice()
        .to_vec()
}

#[test]
fn a_pass_that_consumes_words_does_not_move_the_draws_after_it() {
    // The claim the phase exists for. The stub changes nothing and
    // draws eight words; with one stream shared by the plan those
    // eight words came out from under the NP pass, moving its length
    // draw, its base draws and therefore the molecule.
    let cfg = vj_refdata();
    let without = drawing_vj_plan(&cfg, None);
    let (stub, _log) = DrawsAndChangesNothing::new("test.draws_only", 8);
    let with = drawing_vj_plan(&cfg, Some(Box::new(stub)));

    for seed in 0..32 {
        assert_eq!(
            molecule(&without, &cfg, seed),
            molecule(&with, &cfg, seed),
            "a pass that only consumes words changed the molecule at seed {seed}"
        );
    }
}

#[test]
fn each_pass_draws_the_stream_its_own_label_names() {
    // Ties the run loop to the derivation, for both occurrences of one
    // name: the words a pass sees are the ones `StreamRoot` hands out
    // for its label under the run seed, and nothing else about the
    // plan enters into it. Labelling on the name alone would hand both
    // stubs the same words, which is two `pcr_amplify` steps applying
    // the same errors at the same sites twice.
    let cfg = vj_refdata();
    let mut plan = drawing_vj_plan(&cfg, None);
    let (first, first_log) = DrawsAndChangesNothing::new("test.draws_only", 4);
    let (second, second_log) = DrawsAndChangesNothing::new("test.draws_only", 4);
    plan.push(Box::new(first));
    plan.push(Box::new(second));

    let seed = 4242;
    let _ = molecule(&plan, &cfg, seed);

    let root = StreamRoot::new(seed);
    let expect = |occurrence: u32| -> Vec<u64> {
        let mut rng = root.substream(SubstreamLabel::for_pass("test.draws_only", occurrence));
        (0..4).map(|_| rng.next_u64()).collect()
    };
    let drawn = |log: &Arc<Mutex<Vec<u64>>>| log.lock().expect("not poisoned").clone();

    assert_eq!(drawn(&first_log), expect(0));
    assert_eq!(drawn(&second_log), expect(1));
    assert_ne!(
        drawn(&first_log),
        drawn(&second_log),
        "two passes of the same name drew the same words, so they share a stream"
    );
}
