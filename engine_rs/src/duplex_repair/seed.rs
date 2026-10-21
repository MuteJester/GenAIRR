//! Atomic installation of an already validated, single-layout duplex founder.
use super::{
    bindings::DuplexReferenceBindings,
    projection::{analyze_coupled, seed_capability_reasons},
    AnnotatedCoupledProduct, Substrate,
};
use crate::assignment::{AlleleInstance, SegmentOrientation};
use crate::ir::{
    flag, DuplexFounder, DuplexProvenance, NucFlags, NucHandle, Nucleotide, Region, Segment,
    Simulation, SimulationBuilder,
};
use crate::refdata::RefDataConfig;
use std::sync::Arc;

/// A seed is returned only after bindings, physical replay, placement and the
/// entire pool/region layout have been checked. Overlaps stay supported by the
/// standalone projector, but cannot enter this single-owner execution IR.
pub fn seed_from_coupled(
    annotated: AnnotatedCoupledProduct,
    bindings: DuplexReferenceBindings,
    refdata: &RefDataConfig,
) -> Result<Simulation, String> {
    let analysis = analyze_coupled(&annotated, &bindings, refdata)?;
    let seed_reasons = seed_capability_reasons(
        &analysis.paths,
        &analysis.layout_reasons,
        &bindings,
        annotated.product.sequence.len(),
    );
    if !seed_reasons.is_empty() {
        return Err(format!(
            "unsupported duplex seed layout or native IR range: {}",
            seed_reasons.into_iter().collect::<Vec<_>>().join(", ")
        ));
    }
    let substrates = [Substrate::V, Substrate::D, Substrate::J];
    let segments = [Segment::V, Segment::D, Segment::J];
    let query = annotated.product.sequence.as_bytes().to_vec();
    let mut nucleotides: Vec<Option<Nucleotide>> = vec![None; query.len()];
    let mut regions = Vec::new();
    let mut assignments = Vec::new();
    for (i, path) in analysis.paths.iter().enumerate() {
        let source = bindings.for_substrate(substrates[i]);
        let (start, end) = path
            .query_start
            .zip(path.query_end)
            .ok_or("duplex seed body span is unavailable")?;
        let (five, retained_end) = path
            .reference_start
            .zip(path.reference_end)
            .ok_or("duplex seed reference span is unavailable")?;
        let three = source.sequence.len() - retained_end;
        let mut inst = AlleleInstance::new(bindings.allele_id(substrates[i]));
        inst.trim_5 =
            u16::try_from(five).map_err(|_| "duplex seed trim exceeds native IR range")?;
        inst.trim_3 =
            u16::try_from(three).map_err(|_| "duplex seed trim exceeds native IR range")?;
        if i == 1 && annotated.product.d_inverted {
            inst.orientation = SegmentOrientation::ReverseComplement;
        }
        // P scalar caches are not evidence in this model; the typed sidecar
        // retains their explicit unavailability and projection never reads0.
        assignments.push((segments[i], inst));
        regions.push(Region::new(
            segments[i],
            NucHandle::new(start as u32),
            NucHandle::new(end as u32),
        ));
        for q in start..end {
            nucleotides[q] = Some(Nucleotide::synthetic(
                query[q],
                segments[i],
                NucFlags::empty(),
            ));
        }
        for point in &path.body {
            let g = u16::try_from(point.reference_position)
                .map_err(|_| "duplex seed reference position exceeds native IR range")?;
            if g == u16::MAX {
                return Err(
                    "duplex seed reference position uses reserved native IR sentinel".into(),
                );
            }
            let mut n = Nucleotide::germline(query[point.query_position], g, segments[i]);
            if inst.orientation.is_reverse() {
                n.flags = n.flags.with(flag::INVERTED);
            }
            nucleotides[point.query_position] = Some(n);
        }
    }
    for (i, np) in [Segment::Np1, Segment::Np2].into_iter().enumerate() {
        let start = analysis.paths[i].query_end.unwrap();
        let end = analysis.paths[i + 1].query_start.unwrap();
        if start < end {
            regions.push(Region::new(
                np,
                NucHandle::new(start as u32),
                NucHandle::new(end as u32),
            ));
            for q in start..end {
                nucleotides[q] = Some(Nucleotide::synthetic(query[q], np, NucFlags::empty()));
            }
        }
    }
    let founder = Arc::new(DuplexFounder::new(Arc::new(annotated), Arc::new(bindings))?);
    let provenance = Arc::new(DuplexProvenance::new(founder));
    let mut builder = SimulationBuilder::from_simulation(Simulation::with_capacity(query.len()));
    for (segment, assignment) in assignments {
        builder.assign_allele(segment, assignment);
    }
    for nucleotide in nucleotides {
        builder.push_nucleotide(nucleotide.ok_or("duplex seed has an unclassified query column")?);
    }
    regions.sort_by_key(|r| r.start.index());
    for region in regions {
        builder.add_region(region);
    }
    builder.attach_duplex_provenance(provenance)?;
    Ok(builder.seal())
}
