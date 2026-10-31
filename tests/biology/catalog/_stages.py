"""The mechanisms, in the order a receptor is actually formed.

Stages follow development rather than the source tree: the inherited locus,
then recombination, then the junction, then whether the result can be
expressed, then what maturation does to it, then clonal history. ``assay`` is
last and is different in kind -- not what the cell did, but what the protocol
did to the molecule afterwards. GenAIRR's oracle claim covers both.

A mechanism may have no engine pass: productivity is an assessment the engine
makes rather than a pass it schedules, and clonal history is driven by the
lineage machinery. The completeness test runs the other way round -- every
pass must be classified -- so a mechanism without a pass is not a gap.
"""
from __future__ import annotations

from ._model import Mechanism

#: Developmental order. Reports and directory listings follow it.
STAGES = ("germline", "recombination", "junction", "expression",
          "maturation", "clonality", "assay", "projection")

#: Passes that are not biology. The engine documents these as reference
#: implementations; classifying them as a mechanism would make the catalog
#: assert something false.
INFRASTRUCTURE = {
    "EchoPass": "the engine's deterministic transform reference, not a biological step",
    "SampleBasePass": "the engine's sampling reference, not a biological step",
}

#: Fields no mechanism is answerable for **as truth**, each with the reason. A
#: field nobody is answerable for is a field nobody has to get right, so the
#: alternative to a claim is an explicit entry here -- not silence.
#:
#: Nine of these are also listed in the ``truth_fields`` of the mechanism they
#: belong to, which is deliberate and is not a contradiction: that is where a
#: reader looks for them. An entry here means the field does not count toward
#: its mechanism's coverage, because nothing a test could do would establish a
#: configuration value echoed back, a label the harness assigns, or a ratio
#: derived from two other fields. See ``Mechanism.verifiable_fields``.
UNCLAIMED_FIELDS = {
    "sequence_id": "a label the harness assigns, not simulated biology; "
                   "clonality.founder claims it for clone identity instead",
    "locus": "claimed by germline.genotype as the only exported field it has; "
             "see that mechanism's truth_gap",
    "productivity_policy": "configuration echoed back, not an outcome",
    "productivity_scope": "configuration echoed back, not an outcome",
    "read_layout": "the configured layout name echoed back -- 'paired_end' "
                   "whenever the pass runs, with no per-record outcome in it. "
                   "insert_size is NOT in this category: it is drawn from a "
                   "distribution per record, and assay.paired_end is "
                   "answerable for it",
    "mutation_rate": "a derived ratio of n_mutations to length, owned via "
                     "maturation.hypermutation rather than separately",
    "c_call": "always empty: no isotype or class-switch recombination is "
              "modelled. See maturation.class_switching's truth_gap",
    # The fourth kind of unclaimed field, and the only one that is unclaimable
    # by the nature of the tool rather than by a modelling decision: an
    # inference output of a search that never happens.
    "v_score": "an alignment score: how well the chosen allele explains the "
               "read, which only a search produces. GenAIRR assembles from a "
               "known allele, so there is no hypothesis to score. Exported as "
               "null because the AIRR column exists",
    "d_score": "see v_score",
    "j_score": "see v_score",
    "v_support": "an alignment support -- an e-value or posterior over the "
                 "call. There is no inference here to be confident about: the "
                 "call names the allele the simulation used. Exported as null "
                 "because the AIRR column exists",
    "d_support": "see v_support",
    "j_support": "see v_support",
}

#: Truth fields more than one mechanism is answerable for, each with the reason.
#:
#: Shared ownership is a real relationship -- several mechanisms genuinely write
#: to the same field -- but an *undeclared* overlap means neither owner is
#: accountable, and a field can read as covered because some other mechanism's
#: test happened to touch it. Declaring them makes a new overlap a failure.
SHARED_FIELDS = {
    "sequence": "the assembled bases are what recombination.assembly produces, "
                "what assay.ambiguous_calls reads to show an N is not a "
                "mutation, and the only exported trace of "
                "assay.duplex_formation; all three are answerable for it",
    "original_v_call": "germline.allele_choice records the V first drawn, and "
                       "expression.receptor_revision is what makes it 'original' "
                       "by replacing it; the field exists because of both",
}

MECHANISMS = (
    # ── germline: what the locus inherited ────────────────────────────────
    Mechanism(
        id="germline.genotype",
        stage="germline",
        summary="Which alleles this individual carries, across both haplotypes, "
                "including whether a haplotype is functional at all.",
        truth_fields=("locus",),
        engine_passes=("SampleGenotypePass",),
        truth_gap="The genotype itself is not exported. There is no haplotype or "
                  "chromosome field, no carried-allele set, and no per-read "
                  "novel-allele or deletion truth, so `locus` is the only field "
                  "available to claim. Until that truth is exported, "
                  "genotype and haplotype inference tools (TIgGER, piglet, "
                  "RAbHIT) cannot be benchmarked against GenAIRR.",
    ),
    Mechanism(
        id="germline.allele_choice",
        stage="germline",
        summary="Which V, D and J allele this rearrangement drew from the "
                "genotype it was given.",
        truth_fields=("v_call", "d_call", "j_call", "original_v_call",
                      "recalled_v_call", "recalled_d_call", "recalled_j_call"),
        engine_passes=("SampleAllelePass",),
    ),

    # ── recombination: RAG cuts and rejoins ───────────────────────────────
    Mechanism(
        id="recombination.trimming",
        stage="recombination",
        summary="Exonuclease loss at each coding end before joining: how many "
                "bases of V, D and J were removed.",
        truth_fields=("v_trim_3", "d_trim_5", "d_trim_3", "j_trim_5",
                      "v_trim_5", "j_trim_3"),
        engine_passes=("TrimPass",),
        notes="v_trim_5 and j_trim_3 are exported for AIRR completeness and are "
              "structurally zero: RAG cleaves only the RSS-adjacent coding end, "
              "so there is no 5' V or 3' J exonuclease loss to record. Measured "
              "zero across 200 records. A non-zero value would mean the model "
              "had acquired biology it does not claim.",
    ),
    Mechanism(
        id="recombination.d_inversion",
        stage="recombination",
        summary="A D segment joined in inverted orientation, which the record "
                "must declare because its germline coordinates run backwards.",
        truth_fields=("d_inverted",),
        engine_passes=("InvertDPass",),
    ),
    Mechanism(
        id="recombination.assembly",
        stage="recombination",
        summary="Copying the surviving germline stretches into the rearranged "
                "sequence, which fixes every germline coordinate the record reports.",
        truth_fields=("sequence", "sequence_length", "v_germline_start", "v_germline_end",
                      "d_germline_start", "d_germline_end", "j_germline_start", "j_germline_end"),
        engine_passes=("AssembleSegmentPass",),
    ),

    # ── junction: hairpin opening and untemplated addition ────────────────
    Mechanism(
        id="junction.p_nucleotide",
        stage="junction",
        summary="Palindromic bases left when a hairpinned coding end is opened "
                "asymmetrically -- templated by the segment's own end, but "
                "reversed and complemented.",
        truth_fields=("p_v_3_length", "p_d_5_length", "p_d_3_length", "p_j_5_length"),
        engine_passes=("PAdditionPass", "EmitCodingEndPPass"),
    ),
    Mechanism(
        id="junction.coding_end_processing",
        stage="junction",
        summary="How a coding end is resolved before joining, which decides "
                "what survives of the conserved anchor codons.",
        truth_fields=("junction_start", "junction_end"),
        engine_passes=("ProcessCodingEndsPass",),
    ),
    Mechanism(
        id="junction.n_addition",
        stage="junction",
        summary="TdT-mediated untemplated nucleotide addition at the V-D and "
                "D-J joins. The np fields report the INTERVAL between the "
                "segment bounds, which is not the same thing: see truth_gap.",
        truth_fields=("np1", "np2", "np1_length", "np2_length", "np1_aa", "np2_aa"),
        truth_gap="The TdT-added amount is not exported separately. np1 and np2 "
                  "are the stretch of the molecule between the V and D (and D "
                  "and J) alignment bounds -- that identity holds exactly, "
                  "measured 100 of 100 records with and without an indel pass "
                  "-- and that stretch contains the P nucleotides as well, "
                  "which the documentation says. It also contains any "
                  "polymerase insertion that lands between the segments: "
                  "measured by running one seed with and without "
                  "polymerase_indels, 24 of 200 records had an np length move "
                  "by the insertion, np1 27 -> 28 and np2 2 -> 3 among them. So "
                  "a consumer reading np1_length as the TdT contribution is "
                  "wrong by the number of insertions that fell in the junction, "
                  "and nothing in the record says by how much. The fields are "
                  "right for what they promise; this mechanism's own claim is "
                  "the part that cannot be established from them.",
        engine_passes=("GenerateNPPass",),
    ),
    Mechanism(
        id="junction.assembly",
        stage="junction",
        summary="The junction as reported: its span, its sequence and its "
                "translation, spanning the conserved anchors at both ends.",
        truth_fields=("junction", "junction_aa", "junction_length"),
    ),

    # ── expression: can the rearrangement be used ─────────────────────────
    Mechanism(
        id="expression.productivity",
        stage="expression",
        summary="Whether the rearrangement can produce a receptor: both "
                "conserved anchors accepted, reading frame across the "
                "junction, absence of a stop codon in the coding region, and "
                "the J tail joining the constant region in frame. Five "
                "predicates under three-valued logic, each of them exported.",
        truth_fields=("productive", "stop_codon", "vj_in_frame",
                      # The five predicates and the span they are assessed
                      # over. Exported so a consumer can reproduce the three
                      # flags above rather than taking them on trust, which
                      # the pipeline-space sweep found they could not: the
                      # coding span was missing, so was the frame's phase, and
                      # so were three of the five predicates.
                      "coding_frame_start", "coding_frame_end",
                      "j_tail_in_frame", "v_anchor_accepted",
                      "j_anchor_accepted",
                      "productivity_policy",
                      "productivity_reasons", "productivity_scope"),
        notes="An assessment the engine makes over a finished sequence rather "
              "than a scheduled pass, so it owns no pass.",
    ),
    Mechanism(
        id="expression.receptor_revision",
        stage="expression",
        summary="Secondary rearrangement replacing an already-joined V, which "
                "leaves the original call recoverable alongside the new one.",
        truth_fields=("receptor_revision_applied", "original_v_call"),
        engine_passes=("ReceptorRevisionPass",),
    ),

    # ── maturation: what happens after expression ─────────────────────────
    Mechanism(
        id="maturation.hypermutation",
        stage="maturation",
        summary="Somatic hypermutation: AID-initiated point substitution, "
                "targeted by sequence context under S5F or spread uniformly.",
        truth_fields=("n_mutations", "n_observable_mutations", "n_v_mutations",
                      "n_d_mutations", "n_j_mutations",
                      "n_np_mutations", "mutation_rate", "n_v_unannotated_mutations",
                      "n_cdr1_mutations", "n_cdr2_mutations", "n_fwr1_mutations",
                      "n_fwr2_mutations", "n_fwr3_mutations",
                      "n_fwr1_observable_mutations",
                      "n_cdr1_observable_mutations",
                      "n_fwr2_observable_mutations",
                      "n_cdr2_observable_mutations",
                      "n_fwr3_observable_mutations",
                      "n_v_unannotated_observable_mutations"),
        engine_passes=("S5FMutationPass", "UniformMutationPass"),
    ),

    Mechanism(
        id="maturation.class_switching",
        stage="maturation",
        summary="Class switch recombination: replacing the constant region so a "
                "matured B cell expresses a different isotype.",
        truth_fields=("c_call",),
        truth_gap="Not modelled. c_call is hardwired empty in the record builder, "
                  "there is no IGHC reference, no switch region and no isotype. "
                  "Isotype-aware clonal inference and switched-versus-naive "
                  "partitioning therefore cannot be benchmarked against GenAIRR, "
                  "and the catalog records that rather than omitting the "
                  "mechanism as though receptor maturation ended at SHM.",
    ),

    # ── clonality: one cell's descendants ─────────────────────────────────
    # Split three ways because the exported truth is three different claims,
    # and one mechanism owning `clone_id` left the five lineage fields -- the
    # ground truth for clonal-inference benchmarking -- owned by nobody.
    Mechanism(
        id="clonality.founder",
        stage="clonality",
        summary="Which rearrangement a clone descends from: every member of a "
                "clone shares one founder, and members of different clones "
                "never do.",
        truth_fields=("clone_id", "sequence_id"),
        truth_gap=(
                  "clone_id is the only always-correct clone key. Membership is "
                  "exactly recoverable from the per-base germline provenance on a "
                  "substitution-only pipeline -- identical within a clone, distinct "
                  "across clones, with 1,500 clones giving 1,500 distinct keys "
                  "where the V/D/J call triple gives only 1,464 -- but a length- "
                  "changing pass breaks it: polymerase_indels splits 5 of 5 lineage "
                  "clones and 14 of 20 repertoire clones, and a variable 3' end "
                  "loss splits 5 of 5. So the provenance key verifies clone_id "
                  "where it holds and is not a general substitute for it."
                  ),
        notes="Driven by the lineage and family-replay machinery rather than "
              "by a scheduled pass.",
    ),
    Mechanism(
        id="clonality.tree_topology",
        stage="clonality",
        summary="The genealogy within a clone: who descends from whom, and at "
                "which generation. This is the truth a clonal-inference tool "
                "is scored against.",
        truth_fields=("lineage_node_id", "lineage_parent_id", "lineage_generation"),
        truth_gap=(
                  "The exported fields do not describe the tree. lineage_parent_id "
                  "names no exported record (0 of about 2,400 measured): sampling "
                  "draws only from the final living generation, so a parent is "
                  "always an unobserved node one generation up. lineage_node_id is "
                  "unique within a clone but not globally -- 872 duplicate ids "
                  "among 1,311 records at 20 clones, and at 1,500 clones it takes "
                  "only the values 1 to 5 -- so a join or dedup on it alone merges "
                  "clones. lineage_generation equals max_generations on every "
                  "record (about 2,500 across 22 configurations), correct by "
                  "construction since only final-generation cells are sampled and "
                  "no internal node is ever emitted, and therefore carrying no per- "
                  "record information about tree position. The topology is exported "
                  "only through result.lineage_trees. Key records on (clone_id, "
                  "lineage_node_id) or on sequence_id, both globally unique."
                  ),
    ),
    Mechanism(
        id="clonality.selection",
        stage="clonality",
        summary="Affinity-driven survival across generations, which decides "
                "which lineages expand rather than die out.",
        truth_fields=("lineage_affinity",),
        truth_gap=(
                  "Inert at the defaults and unverifiable when it is not. "
                  "selection_strength defaults to 0, which makes the field the "
                  "sentinel 0.0 on every record (1,883 of 1,883). With strength 1 "
                  "at the default beta of 1.0 it is exp(-beta*d) as documented, but "
                  "underflows what a consumer would expect: 1.8e-230 to 3.1e-169 on "
                  "exported records at the default max_generations=10, and the 0.0 "
                  "sentinel again at max_generations=30. Selection also does not "
                  "shape the trees at the default beta -- strength 0 and strength 1 "
                  "give byte-identical record sets and tree sizes -- while "
                  "beta=0.01 drives every clone extinct. The selection target is on "
                  "no exported field, so affinity cannot be recomputed from the "
                  "record: given a known target it reproduces 69 of 69 from the "
                  "per-region codon rails and 0 of 69 from sequence_aa."
                  ),
    ),
    Mechanism(
        id="clonality.sampling",
        stage="clonality",
        summary="How many times a given sequence was observed, as distinct "
                "from how many cells carried it.",
        truth_fields=("duplicate_count", "lineage_abundance"),
        truth_gap=(
                  "duplicate_count and lineage_abundance are one number stamped "
                  "twice -- the node's sampling draw count -- and never differ in "
                  "any configuration measured; clonal_repertoire emits no "
                  "lineage_abundance at all. Neither is a read count or a tree-cell "
                  "count, and a single record's value is not a clone size whenever "
                  "corruption or hypermutation is configured. What IS recoverable "
                  "is the configured size, as the per-clone sum: n_sample in "
                  "clonal_lineage, and the drawn size in clonal_repertoire on 150 "
                  "of 150 clones. The drawn size itself is on no record or result "
                  "attribute."
                  ),
    ),

    # ── assay: what the protocol did to the molecule ──────────────────────
    Mechanism(
        id="assay.polymerase_indels",
        stage="assay",
        summary="Insertions and deletions introduced while copying, which shift "
                "every downstream coordinate and must stay reconcilable with "
                "the generating allele.",
        truth_fields=("n_indels", "n_v_indels", "n_d_indels", "n_j_indels"),
        truth_gap=("The three per-segment counters do not partition n_indels. "
                   "They cover V, D and J, and an indel landing in an NP "
                   "region belongs to none of them: there is no n_np_indels "
                   "to match n_np_mutations. Measured on 18 of 200 records "
                   "with a deficit of exactly one, and 0 of 200 once the NP "
                   "lengths are forced to zero, which is what identifies the "
                   "cause. A consumer summing the three and expecting "
                   "n_indels is wrong on about 9 per cent of records at the "
                   "default NP lengths."),
        engine_passes=("IndelPass",),
    ),
    Mechanism(
        id="assay.pcr_error",
        stage="assay",
        summary="Misincorporation during amplification: a substitution that is "
                "not somatic hypermutation and must not be counted as one.",
        truth_fields=("n_pcr_errors",),
        engine_passes=("PCRErrorPass",),
    ),
    Mechanism(
        id="assay.quality_error",
        stage="assay",
        summary="Base-calling error at the sequencer, recorded separately from "
                "amplification error.",
        truth_fields=("n_quality_errors",),
        engine_passes=("QualityErrorPass",),
    ),
    Mechanism(
        id="assay.ambiguous_calls",
        stage="assay",
        summary="Positions the instrument could not call, which appear as N and "
                "must not be read as mutations.",
        truth_fields=("sequence",),
        engine_passes=("NCorruptionPass",),
        truth_gap="No N-count field is exported, and because the sites are drawn "
                  "with replacement the requested count is not recoverable from "
                  "the sequence either (20 requested produced 18 Ns). The claim "
                  "that matters -- ambiguous positions must not be read as "
                  "mutations -- has no field of its own to carry it.",
    ),
    Mechanism(
        id="assay.contamination",
        stage="assay",
        summary="A read that did not come from the sampled repertoire at all.",
        truth_fields=("is_contaminant",),
        engine_passes=("ContaminantPass",),
    ),
    Mechanism(
        id="assay.orientation",
        stage="assay",
        summary="Which strand the read was sequenced from, which decides whether "
                "the recorded sequence needs reverse-complementing to align.",
        truth_fields=("rev_comp",),
        engine_passes=("RevCompPass",),
    ),
    Mechanism(
        id="assay.end_loss",
        stage="assay",
        summary="Sequence missing from either end because the amplicon or the "
                "read did not reach it, as distinct from germline trimming.",
        truth_fields=("end_loss_5_length", "end_loss_3_length"),
        engine_passes=("EndLossPass",),
    ),
    Mechanism(
        id="assay.paired_end",
        stage="assay",
        summary="Sampling two reads from one molecule, with the layout and "
                "insert size needed to reconstruct their relationship.",
        truth_fields=("r1_sequence", "r1_start", "r1_end", "r2_sequence", "r2_start",
                      "r2_end", "read_layout", "insert_size"),
        engine_passes=("PairedEndSamplingPass",),
    ),
    Mechanism(
        id="assay.duplex_formation",
        stage="assay",
        summary="Two molecules repaired into one chimeric duplex, whose halves "
                "come from different rearrangements.",
        truth_fields=("sequence", "duplex_repair"),
        engine_passes=("DuplexRecombinePass",),
        notes="The duplex structure IS exported, as the duplex_repair field: the "
              "product, the observation, the overlaps, the evidence, the current "
              "paths, the bindings and the coordinate system. An earlier version "
              "of this entry declared a truth gap saying the structure was not in "
              "the record, which was false -- an existing test reads "
              "airr_record(refs)['duplex_repair'] directly. The completeness "
              "check could not contradict it because no sampled record mode "
              "produced a duplex record; that mode is now sampled.",
    ),

    # ── projection: what a benchmark actually compares ────────────────────
    # Not a developmental stage. These fields are the projection of whatever
    # every mechanism did -- a trim, an N addition, a mutation and an indel all
    # land in the same CIGAR -- which is why no single stage can own them. They
    # are also precisely what an aligner benchmark compares, so leaving them
    # unowned let the catalog read green while every consumed coordinate could
    # be wrong.
    Mechanism(
        id="projection.segment_coordinates",
        stage="projection",
        summary="Where each segment begins and ends, in the read and in the "
                "alignment. The coordinates a downstream aligner is scored on.",
        truth_fields=("v_sequence_start", "v_sequence_end", "d_sequence_start",
                      "d_sequence_end", "j_sequence_start", "j_sequence_end",
                      "v_alignment_start", "v_alignment_end", "d_alignment_start",
                      "d_alignment_end", "j_alignment_start", "j_alignment_end"),
    ),
    Mechanism(
        id="projection.cigar",
        stage="projection",
        summary="Per-segment CIGAR strings reconciling the read against the "
                "germline, including every edit any mechanism introduced.",
        truth_fields=("v_cigar", "d_cigar", "j_cigar"),
        notes="Under an inverted D the germline coordinates are forward-reference "
              "while the sequence carries the reverse complement, and the CIGAR "
              "carries no strand marker -- so reconstructing germline from "
              "(d_call, d_germline_start/end, d_cigar) without reading d_inverted "
              "yields the wrong strand.",
    ),
    Mechanism(
        id="projection.germline_alignment",
        stage="projection",
        summary="The aligned sequence and its germline counterpart, gap-aware, "
                "with the D region optionally masked.",
        truth_fields=("sequence_alignment", "germline_alignment",
                      "germline_alignment_d_mask"),
    ),
    Mechanism(
        id="projection.translation",
        stage="projection",
        summary="The amino-acid projection of the full rearrangement, which "
                "depends on frame and therefore on everything upstream of it.",
        truth_fields=("sequence_aa",),
    ),
    Mechanism(
        id="projection.alignment_scores",
        stage="projection",
        summary="Per-segment identity for the reported call. Score and support "
                "are exported as the AIRR schema requires the columns to be, "
                "and are always null: see UNCLAIMED_FIELDS.",
        truth_fields=("v_identity", "d_identity", "j_identity",
                      "v_score", "d_score", "j_score",
                      "v_support", "d_support", "j_support"),
        notes="score and support are an aligner's inference outputs -- how "
              "well the chosen allele explains the read, and how confident "
              "the choice is. GenAIRR performs no alignment search: the "
              "allele is the one it used, so there is no hypothesis to score "
              "and nothing to be confident about. The columns stay in the "
              "output as null, which the AIRR schema permits, so a consumer "
              "reading by column name finds them. A benchmark that needs a "
              "score has to compute one from the alignment itself.",
    ),
)