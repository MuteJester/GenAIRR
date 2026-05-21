use super::common::{vj_plan, vj_refdata, SEED_RANGE};
use genairr_engine::contract::productive;
use genairr_engine::ir::Simulation;
use genairr_engine::junction::compute_junction;
use genairr_engine::pass::PassRuntime;

#[test]
fn property_productive_implies_in_frame_vj() {
    let refdata = vj_refdata();
    let plan = vj_plan(&refdata);
    let contracts = productive();

    for seed in 0..SEED_RANGE {
        let outcome = PassRuntime::execute_with_context(
            &plan,
            Simulation::default(),
            seed,
            Some(&refdata),
            Some(&contracts),
        );
        let junction = compute_junction(outcome.final_simulation(), &refdata)
            .expect("VJ junction should be defined");
        assert!(
            junction.is_in_frame(),
            "seed {} produced out-of-frame junction (length {})",
            seed,
            junction.length
        );
    }
}
