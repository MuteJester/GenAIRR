use crate::ir::Segment;
use crate::refdata::{Allele, AllelePool};

mod allele_pool;
mod empirical;
mod filtered;
mod trait_objects;
mod uniform;

/// Build a pool of `n` named alleles for testing. The base
/// sequences are all single-byte `b'A'` since we don't care
/// about content here, only sampling behavior.
fn make_pool(n: usize) -> AllelePool {
    let mut p = AllelePool::new();
    for i in 0..n {
        let _ = p.push(Allele {
            name: format!("a{}*01", i),
            gene: format!("a{}", i),
            seq: b"A".to_vec(),
            segment: Segment::V,
            anchor: None,
            functional_status: None,
            subregions: Vec::new(),
        });
    }
    p
}
