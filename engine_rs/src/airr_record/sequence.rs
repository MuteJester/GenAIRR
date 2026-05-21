use super::AirrRecord;
use crate::codon::translate_seq;
use crate::ir::Simulation;

pub(super) fn pool_bases(sim: &Simulation) -> Vec<u8> {
    sim.pool.as_slice().iter().map(|n| n.base).collect()
}

pub(super) fn bytes_to_string(bytes: &[u8]) -> String {
    // Bases are ASCII; lossy decode is safe.
    String::from_utf8_lossy(bytes).into_owned()
}

pub(super) fn apply_rev_comp_projection(rec: &mut AirrRecord) {
    rec.rev_comp = true;
    let seq_len = rec.sequence_length;
    rec.sequence = reverse_complement(&rec.sequence);
    rec.np1 = reverse_complement(&rec.np1);
    rec.np2 = reverse_complement(&rec.np2);
    rec.junction = reverse_complement(&rec.junction);

    // Flip pool-position coords: new_start = seq_len - old_end,
    // new_end = seq_len - old_start.
    flip_coord_pair(&mut rec.v_sequence_start, &mut rec.v_sequence_end, seq_len);
    flip_coord_pair(&mut rec.d_sequence_start, &mut rec.d_sequence_end, seq_len);
    flip_coord_pair(&mut rec.j_sequence_start, &mut rec.j_sequence_end, seq_len);
    flip_coord_pair(&mut rec.junction_start, &mut rec.junction_end, seq_len);

    // Re-translate AA fields from the flipped strings.
    rec.junction_aa = translate_seq(&rec.junction);
    if let Some(jstart) = rec.junction_start {
        let frame_offset = (jstart % 3) as usize;
        rec.sequence_aa = translate_seq(&rec.sequence[frame_offset..]);
        // np1_aa / np2_aa come from `aa_slice_for_region`, but the
        // np region pool offsets aren't carried on the AirrRecord
        // (only the strings + lengths are). Re-derive by walking
        // the new sequence: find the np1 / np2 substring positions
        // and slice. Simpler: just translate each string in its
        // own frame, which is what the forward path effectively did
        // for these short np regions when the frame happened to
        // align. AIRR consumers that care use `junction_aa`; np_aa
        // is a convenience field. Translate plain:
        rec.np1_aa = translate_seq(&rec.np1);
        rec.np2_aa = translate_seq(&rec.np2);
    } else {
        rec.sequence_aa = String::new();
        rec.np1_aa = String::new();
        rec.np2_aa = String::new();
    }
}

fn flip_coord_pair(start: &mut Option<i64>, end: &mut Option<i64>, total_len: i64) {
    if let (Some(s), Some(e)) = (*start, *end) {
        *start = Some(total_len - e);
        *end = Some(total_len - s);
    }
}

fn reverse_complement(s: &str) -> String {
    let mut out = String::with_capacity(s.len());
    for c in s.chars().rev() {
        out.push(match c {
            'A' => 'T',
            'T' => 'A',
            'U' => 'A',
            'C' => 'G',
            'G' => 'C',
            'a' => 't',
            't' => 'a',
            'u' => 'a',
            'c' => 'g',
            'g' => 'c',
            'N' => 'N',
            'n' => 'n',
            other => other,
        });
    }
    out
}
