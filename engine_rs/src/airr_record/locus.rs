use crate::refdata::RefDataConfig;

const AIRR_LOCI: [&str; 7] = ["IGH", "IGK", "IGL", "TRA", "TRB", "TRG", "TRD"];

/// refdata-driven locus fallback. Walks each pool's
/// first allele in turn and returns the locus prefix when the name
/// starts with one of the AIRR loci. Used when live-call evidence
/// has wiped every `*_call` (heavy SHM under corruption stack) but
/// the chain is still well-defined by the source data.
pub(super) fn locus_from_refdata(refdata: &RefDataConfig) -> String {
    for entry in [
        refdata.v_pool.iter().next(),
        refdata.j_pool.iter().next(),
        refdata.d_pool.iter().next(),
    ]
    .into_iter()
    .flatten()
    {
        let candidate = locus_prefix(&entry.1.name);
        if !candidate.is_empty() {
            return candidate;
        }
    }
    String::new()
}

fn locus_prefix(name: &str) -> String {
    if name.len() < 3 {
        return String::new();
    }
    let mut prefix = String::with_capacity(3);
    for c in name.chars().take(3) {
        prefix.push(c.to_ascii_uppercase());
    }
    if AIRR_LOCI.contains(&prefix.as_str()) {
        prefix
    } else {
        String::new()
    }
}

pub(super) fn derive_locus(v_call: &str, j_call: &str, d_call: &str) -> String {
    for name in [v_call, j_call, d_call] {
        if name.is_empty() {
            continue;
        }
        let candidate = locus_prefix(name);
        if !candidate.is_empty() {
            return candidate;
        }
    }
    String::new()
}
