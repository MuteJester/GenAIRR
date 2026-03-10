/* ───────────────────────────────────────────────────────────
 *  DissectionBay.tsx — Full sequence dissection viewer.
 *  Assembled color bar, junction bracket, exploded segment
 *  panels, mutation sparklines, P|N|P bars, nucleotide detail.
 *  Ported from v0 template → Docusaurus CSS Modules.
 * ─────────────────────────────────────────────────────────── */

import React, { useState, useMemo, useCallback, useRef, useEffect } from 'react';
import { motion, AnimatePresence } from 'framer-motion';
import { AlertTriangle, X, Dna, Scissors, Zap, FlipVertical, ChevronRight, Copy, Check } from 'lucide-react';
import { SEGMENT_COLORS, type SegmentId, type SegmentMeta, type SequenceRecord } from './types';
import styles from './styles.module.css';

// ── Constants ────────────────────────────────────────────

const SPRING = { type: 'spring' as const, stiffness: 260, damping: 22 };
const STAGGER = 0.06;

function pct(part: number, total: number) {
  return Math.max((part / total) * 100, 1.5);
}

function formatBp(n: number) {
  return `${n} bp`;
}

// ── MutationSparkline ────────────────────────────────────

function MutationSparkline({
  start, end, mutations, height = 24,
}: {
  start: number; end: number;
  mutations: Record<number, string>;
  height?: number;
}) {
  const len = end - start;
  const bins = Math.min(len, 60);
  const binSize = Math.max(1, Math.floor(len / bins));
  const counts = Array.from({ length: bins }, (_, i) => {
    const binStart = start + i * binSize;
    const binEnd = binStart + binSize;
    return Object.keys(mutations).filter(
      (k) => Number(k) >= binStart && Number(k) < binEnd,
    ).length;
  });
  const max = Math.max(...counts, 1);

  return (
    <svg
      width="100%" height={height}
      viewBox={`0 0 ${bins} ${max}`}
      preserveAspectRatio="none"
      className={styles.sparklineWrap}
    >
      {counts.map((c, i) => (
        <rect
          key={i} x={i} y={max - c}
          width={0.85} height={c || 0.15}
          fill={c > 0 ? '#DC2626' : 'var(--pg-border)'}
          opacity={c > 0 ? 0.85 : 0.25}
        />
      ))}
    </svg>
  );
}

// ── PNPBar ───────────────────────────────────────────────

function PNPBar({
  prefix, nRegion, suffix, color,
}: {
  prefix: string; nRegion: string; suffix: string; color: string;
}) {
  const total = prefix.length + nRegion.length + suffix.length;
  if (total === 0) return null;
  const parts = [
    { label: 'P', len: prefix.length, opacity: 0.6 },
    { label: 'N', len: nRegion.length, opacity: 0.25 },
    { label: 'P', len: suffix.length, opacity: 0.6 },
  ];
  return (
    <div className={styles.pnpBar}>
      {parts.map((p, i) =>
        p.len > 0 ? (
          <div
            key={i}
            className={styles.pnpBarPart}
            style={{ width: `${pct(p.len, total)}%`, backgroundColor: color, opacity: p.opacity }}
          >
            <span className={styles.pnpBarLabel}>{p.label}</span>
          </div>
        ) : null,
      )}
    </div>
  );
}

// ── NucleotideStrip ──────────────────────────────────────

function NucleotideStrip({
  sequence, startPos, color, mutations,
}: {
  sequence: string; startPos: number; color: string;
  mutations: Record<number, string>;
}) {
  return (
    <div
      className={styles.ntStrip}
      style={{ backgroundColor: `color-mix(in srgb, ${color} 3%, var(--pg-background))` }}
    >
      {Array.from({ length: Math.ceil(sequence.length / 80) }, (_, row) => {
        const s = row * 80;
        const chunk = sequence.slice(s, s + 80);
        return (
          <div key={row} className={styles.ntStripRow}>
            <span className={styles.ntStripLineNum}>{startPos + s + 1}</span>
            <div className={styles.ntStripChars}>
              {chunk.split('').map((ch, ci) => {
                const absPos = startPos + s + ci;
                const isMut = absPos in mutations;
                return (
                  <span
                    key={ci}
                    className={`${styles.ntStripChar} ${isMut ? styles.ntStripCharMut : ''}`}
                    style={!isMut ? { color: 'var(--pg-foreground)' } : undefined}
                  >
                    {isMut && <span className={styles.ntStripMutDot} />}
                    {ch}
                  </span>
                );
              })}
            </div>
            <span className={styles.ntStripEndNum}>
              {Math.min(startPos + s + 80, startPos + sequence.length)}
            </span>
          </div>
        );
      })}
    </div>
  );
}

// ── TrimDiagram ──────────────────────────────────────────

function TrimDiagram({ label, bases }: { label: string; bases: number }) {
  if (!bases) return null;
  return (
    <div className={styles.trimDiagram}>
      <Scissors size={12} style={{ color: 'var(--pg-muted-fg)' }} />
      <span className={styles.trimText}>
        {label}: <span className={styles.trimValue}>{bases} bp</span> trimmed
      </span>
    </div>
  );
}

// ── MetricCell ───────────────────────────────────────────

function MetricCell({ label, value, color }: { label: string; value: string; color?: string }) {
  return (
    <div className={styles.metricCell}>
      <span className={styles.metricCellLabel}>{label}</span>
      <span className={styles.metricCellValue} style={color ? { color } : undefined}>
        {value}
      </span>
    </div>
  );
}

// ── SegmentPanel ─────────────────────────────────────────

function SegmentPanel({
  segment, result, isActive, index,
}: {
  segment: SegmentMeta; result: SequenceRecord;
  isActive: boolean; index: number;
}) {
  const subSeq = result.sequence.slice(segment.start, segment.end);
  const segLen = segment.end - segment.start;
  const segMutations = useMemo(() => {
    const m: Record<number, string> = {};
    Object.entries(result.mutations).forEach(([k, v]) => {
      const pos = Number(k);
      if (pos >= segment.start && pos < segment.end) m[pos] = v;
    });
    return m;
  }, [result.mutations, segment]);

  const mutCount = Object.keys(segMutations).length;

  return (
    <motion.div
      layout
      initial={{ opacity: 0, y: 24, scale: 0.96 }}
      animate={{ opacity: 1, y: 0, scale: 1 }}
      exit={{ opacity: 0, y: 12, scale: 0.97 }}
      transition={{ ...SPRING, delay: index * STAGGER }}
      className={`${styles.segPanel} ${isActive ? styles.segPanelActive : ''}`}
      style={{ borderTopWidth: 3, borderTopColor: segment.color }}
    >
      {/* Header */}
      <div className={styles.segPanelHeader}>
        <div className={styles.segPanelChip} style={{ backgroundColor: segment.color }}>
          <span className={styles.segPanelChipText}>{segment.label}</span>
        </div>
        <span className={styles.segPanelGene}>{segment.fullLabel}</span>
        <span className={styles.segPanelPos}>
          {segment.start + 1}..{segment.end} ({segLen} nt)
        </span>
      </div>

      {/* Body */}
      <div className={styles.segPanelBody}>
        {/* V segment */}
        {segment.id === 'V' && (
          <>
            <div className={`${styles.metricGrid} ${styles.metricGrid3}`}>
              <MetricCell label="LENGTH" value={formatBp(segLen)} />
              <MetricCell label="MUTATIONS" value={String(mutCount)} color={mutCount > 0 ? '#DC2626' : undefined} />
              <MetricCell label="3' TRIM" value={formatBp(result.v_trim_3)} />
            </div>
            <div>
              <span className={styles.microLabel}>Mutation Density</span>
              <MutationSparkline start={segment.start} end={segment.end} mutations={result.mutations} />
            </div>
            <TrimDiagram label="V-gene 3'" bases={result.v_trim_3} />
          </>
        )}

        {/* D segment */}
        {segment.id === 'D' && (
          <>
            <div className={`${styles.metricGrid} ${styles.metricGrid4}`}>
              <MetricCell label="LENGTH" value={formatBp(segLen)} />
              <MetricCell label="5' TRIM" value={formatBp(result.d_trim_5)} />
              <MetricCell label="3' TRIM" value={formatBp(result.d_trim_3)} />
              <MetricCell label="INVERTED" value={result.d_inverted ? 'YES' : 'NO'} color={result.d_inverted ? '#DC2626' : undefined} />
            </div>
            {result.d_inverted && (
              <div className={styles.invBadge}>
                <FlipVertical size={12} style={{ color: 'var(--pg-destructive)' }} />
                <span className={styles.invBadgeText}>D-gene inverted (reverse complement used)</span>
              </div>
            )}
            <div style={{ display: 'flex', gap: '0.5rem' }}>
              <TrimDiagram label="5'" bases={result.d_trim_5} />
              <TrimDiagram label="3'" bases={result.d_trim_3} />
            </div>
          </>
        )}

        {/* J segment */}
        {segment.id === 'J' && (
          <>
            <div className={`${styles.metricGrid} ${styles.metricGrid3}`}>
              <MetricCell label="LENGTH" value={formatBp(segLen)} />
              <MetricCell label="MUTATIONS" value={String(mutCount)} color={mutCount > 0 ? '#DC2626' : undefined} />
              <MetricCell label="5' TRIM" value={formatBp(result.j_trim_5)} />
            </div>
            <div>
              <span className={styles.microLabel}>Mutation Density</span>
              <MutationSparkline start={segment.start} end={segment.end} mutations={result.mutations} />
            </div>
            <TrimDiagram label="J-gene 5'" bases={result.j_trim_5} />
          </>
        )}

        {/* NP1 segment */}
        {segment.id === 'NP1' && (
          <>
            <div className={`${styles.metricGrid} ${styles.metricGrid3}`}>
              <MetricCell label="P-PREFIX" value={formatBp(result.np1_p_prefix.length)} />
              <MetricCell label="N-ADDITION" value={formatBp(result.np1_n_region.length)} />
              <MetricCell label="P-SUFFIX" value={formatBp(result.np1_p_suffix.length)} />
            </div>
            <div>
              <span className={styles.microLabel}>P | N | P Composition</span>
              <PNPBar prefix={result.np1_p_prefix} nRegion={result.np1_n_region} suffix={result.np1_p_suffix} color={segment.color} />
            </div>
            <div className={styles.npSeqDisplay}>
              <span style={{ color: segment.color, opacity: 0.7 }}>{result.np1_p_prefix}</span>
              <span style={{ color: 'var(--pg-foreground)' }}>{result.np1_n_region}</span>
              <span style={{ color: segment.color, opacity: 0.7 }}>{result.np1_p_suffix}</span>
            </div>
          </>
        )}

        {/* NP2 segment */}
        {segment.id === 'NP2' && (
          <>
            <div className={`${styles.metricGrid} ${styles.metricGrid3}`}>
              <MetricCell label="P-PREFIX" value={formatBp(result.np2_p_prefix.length)} />
              <MetricCell label="N-ADDITION" value={formatBp(result.np2_n_region.length)} />
              <MetricCell label="P-SUFFIX" value={formatBp(result.np2_p_suffix.length)} />
            </div>
            <div>
              <span className={styles.microLabel}>P | N | P Composition</span>
              <PNPBar prefix={result.np2_p_prefix} nRegion={result.np2_n_region} suffix={result.np2_p_suffix} color={segment.color} />
            </div>
            <div className={styles.npSeqDisplay}>
              <span style={{ color: segment.color, opacity: 0.7 }}>{result.np2_p_prefix}</span>
              <span style={{ color: 'var(--pg-foreground)' }}>{result.np2_n_region}</span>
              <span style={{ color: segment.color, opacity: 0.7 }}>{result.np2_p_suffix}</span>
            </div>
          </>
        )}

        {/* Nucleotide detail on active */}
        <AnimatePresence>
          {isActive && (
            <motion.div
              initial={{ height: 0, opacity: 0 }}
              animate={{ height: 'auto', opacity: 1 }}
              exit={{ height: 0, opacity: 0 }}
              transition={{ duration: 0.25 }}
              style={{ overflow: 'hidden' }}
            >
              <span className={styles.microLabel}>
                Nucleotide Sequence
                {mutCount > 0 && (
                  <span style={{ marginLeft: '0.5rem', color: 'var(--pg-destructive)' }}>
                    ({mutCount} mutations highlighted)
                  </span>
                )}
              </span>
              <NucleotideStrip
                sequence={subSeq}
                startPos={segment.start}
                color={segment.color}
                mutations={segMutations}
              />
            </motion.div>
          )}
        </AnimatePresence>
      </div>
    </motion.div>
  );
}

// ── Main DissectionBay ───────────────────────────────────

interface DissectionBayProps {
  result: SequenceRecord;
  onClose: () => void;
}

export default function DissectionBay({ result, onClose }: DissectionBayProps) {
  const [activeSegment, setActiveSegment] = useState<SegmentId | null>(null);
  const [copied, setCopied] = useState(false);
  const bayRef = useRef<HTMLElement>(null);

  const seqLen = result.sequence.length;

  useEffect(() => {
    bayRef.current?.scrollIntoView({ behavior: 'smooth', block: 'start' });
  }, [result.id]);

  const segments: SegmentMeta[] = useMemo(() => [
    { id: 'V', label: 'V', fullLabel: result.v_call, color: SEGMENT_COLORS.V, start: result.v_sequence_start, end: result.v_sequence_end },
    { id: 'NP1', label: 'NP1', fullLabel: 'N/P Region 1', color: SEGMENT_COLORS.NP1, start: result.v_sequence_end, end: result.d_sequence_start },
    { id: 'D', label: 'D', fullLabel: result.d_call, color: SEGMENT_COLORS.D, start: result.d_sequence_start, end: result.d_sequence_end },
    { id: 'NP2', label: 'NP2', fullLabel: 'N/P Region 2', color: SEGMENT_COLORS.NP2, start: result.d_sequence_end, end: result.j_sequence_start },
    { id: 'J', label: 'J', fullLabel: result.j_call, color: SEGMENT_COLORS.J, start: result.j_sequence_start, end: result.j_sequence_end },
  ], [result]);

  const handleSegmentClick = useCallback((id: SegmentId) => {
    setActiveSegment((prev) => (prev === id ? null : id));
  }, []);

  const handleCopy = useCallback(() => {
    navigator.clipboard.writeText(result.sequence);
    setCopied(true);
    setTimeout(() => setCopied(false), 1500);
  }, [result.sequence]);

  const junctionLeftPct = (result.junction_start / seqLen) * 100;
  const junctionWidthPct = ((result.junction_end - result.junction_start) / seqLen) * 100;
  const totalMutations = Object.keys(result.mutations).length;

  return (
    <section ref={bayRef} className={styles.dissectionBay} aria-label={`Sequence dissection: ${result.sequence_id}`}>
      {/* ── Top bar ───────────────────────────────── */}
      <div className={styles.dissectionTopBar}>
        <Dna size={16} style={{ color: 'var(--pg-accent)' }} />
        <div className={styles.dissectionTitleGroup}>
          <h2 className={styles.dissectionTitle}>Sequence Dissection</h2>
          <span className={styles.dissectionSeqId}>{result.sequence_id}</span>
        </div>

        <span className={result.productive ? styles.dissectionBadge : `${styles.dissectionBadge} ${styles.dissectionBadgeNeg}`}>
          {result.productive ? 'Productive' : 'Non-productive'}
        </span>

        <div className={styles.dissectionActions}>
          <button type="button" onClick={handleCopy} className={styles.dissectionBtn}>
            {copied ? <Check size={12} /> : <Copy size={12} />}
            {copied ? 'Copied' : 'Copy Seq'}
          </button>
          <button type="button" onClick={onClose} className={styles.dissectionBtn} aria-label="Close dissection">
            <X size={12} /> Close
          </button>
        </div>
      </div>

      {/* ── Summary metrics ───────────────────────── */}
      <div className={styles.summaryGrid}>
        {[
          { label: 'V-GENE', val: result.v_call, c: SEGMENT_COLORS.V, mono: false },
          { label: 'D-GENE', val: result.d_call, c: SEGMENT_COLORS.D, mono: false },
          { label: 'J-GENE', val: result.j_call, c: SEGMENT_COLORS.J, mono: false },
          { label: 'JUNCTION AA', val: result.junction_aa, c: undefined, mono: true },
          { label: 'TOTAL LENGTH', val: `${seqLen} nt`, c: undefined, mono: false },
          { label: 'MUTATIONS', val: `${totalMutations} (${(result.mutation_rate * 100).toFixed(1)}%)`, c: totalMutations > 0 ? '#DC2626' : undefined, mono: false },
          { label: 'CDR3 LENGTH', val: `${result.cdr3_length} aa`, c: undefined, mono: false },
          { label: 'PRODUCTIVE', val: result.productive ? 'Yes' : 'No', c: result.productive ? undefined : '#DC2626', mono: false },
        ].map((item) => (
          <div key={item.label} className={styles.summaryCell}>
            <span className={styles.summaryCellLabel}>{item.label}</span>
            <span
              className={`${styles.summaryCellValue} ${item.mono ? styles.summaryCellMono : ''}`}
              style={item.c ? { color: item.c } : undefined}
            >
              {item.val}
            </span>
          </div>
        ))}
      </div>

      {/* ── Corruption warnings ───────────────────── */}
      {(!!result.corruption_5prime || !!result.corruption_3prime) && (
        <div className={styles.corruptionRow}>
          {!!result.corruption_5prime && (
            <motion.div initial={{ opacity: 0, x: -12 }} animate={{ opacity: 1, x: 0 }} className={styles.corruptionBadge}>
              <AlertTriangle size={14} style={{ color: 'var(--pg-destructive)' }} />
              <span className={styles.corruptionText}>5&prime; Corruption: +{result.corruption_5prime} nt added</span>
            </motion.div>
          )}
          {!!result.corruption_3prime && (
            <motion.div initial={{ opacity: 0, x: 12 }} animate={{ opacity: 1, x: 0 }} className={styles.corruptionBadge}>
              <AlertTriangle size={14} style={{ color: 'var(--pg-destructive)' }} />
              <span className={styles.corruptionText}>3&prime; Corruption: -{result.corruption_3prime} nt removed</span>
            </motion.div>
          )}
        </div>
      )}

      {/* ── Main visualization ────────────────────── */}
      <div className={styles.assembledArea}>
        <div className={styles.assembledLabel}>
          <span className={styles.microLabel}>Assembled Sequence</span>
          <span className={styles.assembledHint}>Click a segment below to inspect nucleotide detail</span>
        </div>

        {/* Assembled bar */}
        <div
          className={styles.assembledBar}
          style={{
            borderLeft: result.corruption_5prime ? '3px dashed var(--pg-destructive)' : undefined,
            borderRight: result.corruption_3prime ? '3px dashed var(--pg-destructive)' : undefined,
          }}
        >
          {segments.map((seg) => {
            const widthPct = ((seg.end - seg.start) / seqLen) * 100;
            const isActive = activeSegment === seg.id;
            return (
              <motion.div
                key={seg.id}
                className={styles.assembledSegment}
                style={{ width: `${widthPct}%`, backgroundColor: seg.color }}
                initial={false}
                animate={{
                  opacity: activeSegment === null || isActive ? 1 : 0.45,
                  filter: activeSegment === null || isActive ? 'none' : 'grayscale(0.4)',
                }}
                whileHover={{ opacity: 1, filter: 'none' }}
                transition={{ duration: 0.2 }}
                onClick={() => handleSegmentClick(seg.id)}
              >
                {widthPct > 4 && (
                  <span className={styles.assembledSegmentLabel}>{seg.label}</span>
                )}
                <span className={`${styles.assembledSegmentPos} ${styles.assembledSegmentPosLeft}`}>
                  {seg.start + 1}
                </span>
                {widthPct > 8 && (
                  <span className={`${styles.assembledSegmentPos} ${styles.assembledSegmentPosRight}`}>
                    {seg.end}
                  </span>
                )}
              </motion.div>
            );
          })}

          {/* Mutation dots */}
          {Object.keys(result.mutations).map((posStr) => {
            const pos = Number(posStr);
            return (
              <div key={posStr} className={styles.mutDot} style={{ left: `${(pos / seqLen) * 100}%` }} />
            );
          })}
        </div>

        {/* Junction bracket */}
        <div className={styles.junctionBracket} aria-label={`CDR3 junction: ${result.junction_aa}`}>
          <div className={styles.junctionBracketInner} style={{ left: `${junctionLeftPct}%`, width: `${junctionWidthPct}%` }}>
            <svg style={{ width: '100%', height: 12 }} preserveAspectRatio="none">
              <line x1="0" y1="0" x2="0" y2="10" stroke="#9B6FC4" strokeWidth="1.5" />
              <line x1="0" y1="10" x2="100%" y2="10" stroke="#9B6FC4" strokeWidth="1.5" />
              <line x1="100%" y1="0" x2="100%" y2="10" stroke="#9B6FC4" strokeWidth="1.5" />
            </svg>
            <span className={styles.junctionBracketAA}>CDR3: {result.junction_aa}</span>
          </div>
        </div>

        {/* ── Exploded segments ────────────────────── */}
        <div style={{ marginTop: '1rem' }}>
          <div style={{ display: 'flex', alignItems: 'center', gap: '0.5rem', marginBottom: '0.75rem' }}>
            <Zap size={12} style={{ color: 'var(--pg-accent)' }} />
            <span className={styles.microLabel}>Exploded Segments</span>
          </div>

          {/* Leader lines */}
          <div className={styles.leaderRow}>
            {segments.map((seg, i) => (
              <motion.div
                key={seg.id}
                initial={{ opacity: 0 }}
                animate={{ opacity: 1 }}
                transition={{ delay: i * STAGGER + 0.1 }}
                className={styles.leaderLine}
              >
                <div
                  className={styles.leaderLineStem}
                  style={{ backgroundColor: activeSegment === seg.id ? seg.color : undefined }}
                />
                <ChevronRight
                  size={10}
                  style={{
                    transform: 'rotate(90deg)',
                    color: activeSegment === seg.id ? seg.color : 'var(--pg-border)',
                  }}
                />
              </motion.div>
            ))}
          </div>

          {/* Segment panels grid */}
          <AnimatePresence mode="wait">
            <motion.div
              key={result.id}
              className={styles.segmentGrid}
              initial={{ opacity: 0 }}
              animate={{ opacity: 1 }}
              exit={{ opacity: 0 }}
              transition={{ duration: 0.2 }}
            >
              {segments.map((seg, i) => (
                <div key={seg.id} onClick={() => handleSegmentClick(seg.id)} style={{ cursor: 'pointer' }}>
                  <SegmentPanel
                    segment={seg}
                    result={result}
                    isActive={activeSegment === seg.id}
                    index={i}
                  />
                </div>
              ))}
            </motion.div>
          </AnimatePresence>
        </div>

        {/* ── Full-width nt detail for active segment ── */}
        <AnimatePresence>
          {activeSegment && (
            <motion.div
              initial={{ height: 0, opacity: 0 }}
              animate={{ height: 'auto', opacity: 1 }}
              exit={{ height: 0, opacity: 0 }}
              transition={{ duration: 0.25 }}
              style={{ overflow: 'hidden' }}
            >
              {(() => {
                const seg = segments.find((s) => s.id === activeSegment)!;
                const subSeq = result.sequence.slice(seg.start, seg.end);
                const segMuts: Record<number, string> = {};
                Object.entries(result.mutations).forEach(([k, v]) => {
                  const pos = Number(k);
                  if (pos >= seg.start && pos < seg.end) segMuts[pos] = v;
                });
                return (
                  <div className={styles.ntDetailPanel}>
                    <div className={styles.ntDetailHeader}>
                      <div className={styles.ntDetailHeaderLeft}>
                        <div className={styles.ntDetailColorBox} style={{ backgroundColor: seg.color }} />
                        <span className={styles.ntDetailLabel}>
                          {seg.label} Region — Full Nucleotide Sequence
                        </span>
                        <span className={styles.ntDetailPos}>
                          ({seg.start + 1}..{seg.end}, {seg.end - seg.start} nt)
                        </span>
                      </div>
                      {Object.keys(segMuts).length > 0 && (
                        <span className={styles.ntDetailMutCount}>
                          {Object.keys(segMuts).length} mutations
                        </span>
                      )}
                    </div>
                    <NucleotideStrip
                      sequence={subSeq}
                      startPos={seg.start}
                      color={seg.color}
                      mutations={segMuts}
                    />
                  </div>
                );
              })()}
            </motion.div>
          )}
        </AnimatePresence>
      </div>
    </section>
  );
}
