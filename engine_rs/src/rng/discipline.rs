//! Guard: every generator in production code is born in a place that
//! says why it is not a substream.
//!
//! [`Rng::new`](crate::rng::Rng::new) takes a seed word and returns a
//! stream. That is the right primitive and the wrong default. A
//! stream taken straight from a seed is coupled to every other stream
//! taken from the same seed — which is what
//! [`crate::rng::substream`] exists to prevent — and the old way to
//! avoid the coupling was to salt the seed by hand:
//! `Rng::new(seed ^ SOME_SALT)`. That does not work. XOR is
//! invertible, so a salted stream *is* the plain stream of another
//! seed, and for any two salts there is a pair of seeds at which one
//! run's stream is bit-identical to another run's. The measurement is
//! in [`crate::rng::substream`]'s tests.
//!
//! `Rng::new` cannot be made private: integration tests and the
//! derivation in `substream` both need it. So the rule is enforced
//! here instead. Each production file that constructs a generator
//! from a seed word appears in [`SANCTIONED`] with the number of
//! times it does so and the reason that is correct; anything else
//! fails. Adding a site means saying why, in the same place the
//! reader will look to find out whether this engine still has one way
//! to be born a stream.
//!
//! Test code is exempt. A test that pins a draw sequence needs a
//! stream with no history, and `Rng::new(42)` is the clearest way to
//! write one.

use crate::source_guard::mentions_in_production;

/// What the guard looks for: `Rng::new(`, however it is qualified.
const CONSTRUCTOR: &str = "Rng::new(";

/// This file names the constructor in its own needle and in its own
/// synthetic fixtures, neither of which is a call site.
const GUARD_ITSELF: &str = "src/rng/discipline.rs";

/// A production file allowed to construct a generator from a seed
/// word, the number of times it does so, and why a substream is not
/// what that site wants.
struct Sanctioned {
    /// Crate-relative, forward-slashed.
    file: &'static str,
    /// How many production call sites this file is allowed. An exact
    /// budget, not a ceiling: a new site has to be justified here,
    /// and a removed one has to be taken off the list.
    sites: usize,
    why: &'static str,
}

const SANCTIONED: &[Sanctioned] = &[
    Sanctioned {
        file: "src/rng/substream.rs",
        sites: 1,
        why: "`StreamRoot::substream` is the derivation itself -- it \
              folds a label into the run seed and hands the result to \
              `Rng::new`. This is the one place a seed word becomes a \
              stream on purpose.",
    },
    Sanctioned {
        file: "src/lineage/branching.rs",
        sites: 1,
        why: "a child's seed is a fresh word its parent drew, not a \
              position in a stream shared with other roles: one \
              division runs one mutator, so there is nothing for a \
              substream to separate it from. The family's three \
              seed-derived streams are in `lineage::streams`.",
    },
    Sanctioned {
        file: "src/python/lineage.rs",
        sites: 1,
        why: "`sample_clone_sizes` is a one-shot utility whose seed is \
              the caller's own argument and whose draws are the only \
              draws in the call. It is a root, not a role.",
    },
];

#[test]
fn every_production_generator_is_born_in_a_sanctioned_place() {
    let mut counts: std::collections::BTreeMap<String, Vec<String>> = Default::default();
    for mention in mentions_in_production(CONSTRUCTOR) {
        if mention.file == GUARD_ITSELF {
            continue;
        }
        counts
            .entry(mention.file.clone())
            .or_default()
            .push(mention.to_string());
    }

    let mut complaints = Vec::new();
    for (file, sites) in &counts {
        match SANCTIONED.iter().find(|s| s.file == *file) {
            None => complaints.push(format!(
                "{} constructs a generator from a seed word in production code, \
                 and is not in `SANCTIONED`. Draw from a substream \
                 (`StreamRoot::new(seed).substream(label)`), or add the file to \
                 `SANCTIONED` with the reason a substream is wrong here. Sites:\n  {}",
                file,
                sites.join("\n  ")
            )),
            Some(entry) if entry.sites != sites.len() => complaints.push(format!(
                "{} is sanctioned for {} production site(s) and has {}. If the new \
                 one is right, raise the count and extend the reason; if a site went \
                 away, lower it.\nThe reason on record is: {}\nSites:\n  {}",
                file,
                entry.sites,
                sites.len(),
                entry.why,
                sites.join("\n  ")
            )),
            Some(_) => {}
        }
    }

    for entry in SANCTIONED {
        if !counts.contains_key(entry.file) {
            complaints.push(format!(
                "{} is sanctioned for {} production site(s) and has none. Remove it \
                 from `SANCTIONED` -- a stale exemption is an exemption nobody \
                 reviewed.",
                entry.file, entry.sites
            ));
        }
    }

    assert!(
        complaints.is_empty(),
        "The engine has one way to be born a stream and {} exception(s) to it.\n\n{}",
        SANCTIONED.len(),
        complaints.join("\n\n")
    );
}

#[test]
fn the_guard_reads_the_constructor_however_it_is_qualified() {
    use crate::source_guard::mentions_in_body;

    let found = mentions_in_body(
        "src/synthetic.rs",
        "\
fn a() { let mut r = Rng::new(seed); }
fn b() { let mut r = crate::rng::Rng::new(seed); }
fn c() { let mut r = StreamRoot::new(seed).substream(label); }
",
        &[CONSTRUCTOR],
    );
    assert_eq!(
        found.iter().map(|m| m.line).collect::<Vec<_>>(),
        vec![1, 2],
        "both spellings of the constructor count; deriving a substream does not"
    );
}

#[test]
fn no_sanctioned_file_is_listed_twice() {
    let mut files: Vec<&str> = SANCTIONED.iter().map(|s| s.file).collect();
    files.sort_unstable();
    let count = files.len();
    files.dedup();
    assert_eq!(
        files.len(),
        count,
        "a file listed twice would let the second entry's count go unchecked"
    );
}
