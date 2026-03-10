/* ───────────────────────────────────────────────────────────
 *  PipelineBuilder.tsx — Swiss Lab pipeline with OpCards
 *  and AddOp popover.
 *  Uses our ops.ts (20 real ops) + framer-motion animations.
 * ─────────────────────────────────────────────────────────── */

import React, { useState, useCallback, useRef, useEffect } from 'react';
import ReactDOM from 'react-dom';
import { motion, AnimatePresence } from 'framer-motion';
import {
  Plus, Lock, X, ChevronDown, ChevronRight,
  AlertCircle, Dna, Zap, FlaskConical, Bug,
} from 'lucide-react';
import {
  OP_DEFS, CATEGORY_META, getOpDef, canAddOp, getOpsByCategory,
  type OpDef, type Category,
} from './ops';
import type { PipelineOp } from './pythonRunner';
import styles from './styles.module.css';

// ── Constants ────────────────────────────────────────────

const SPRING = { type: 'spring' as const, stiffness: 300, damping: 26 };
const STAGGER = 0.04;

const CATEGORY_ICONS: Record<Category, typeof Dna> = {
  rearrangement: Dna,
  mutation: Zap,
  biology: FlaskConical,
  artifacts: Bug,
};

// ── Types ────────────────────────────────────────────────

interface PipelineBuilderProps {
  pipeline: PipelineOp[];
  pipelineOpIds: string[];
  selectedInstanceId: string | null;
  onSelectOp: (id: string | null) => void;
  onRemoveOp: (id: string) => void;
  onAddOp: (opId: string) => void;
}

// ── OpCard ────────────────────────────────────────────────

function OpCard({
  op,
  pipelineOp,
  index,
  isSelected,
  onSelect,
  onRemove,
}: {
  op: OpDef;
  pipelineOp: PipelineOp;
  index: number;
  isSelected: boolean;
  onSelect: () => void;
  onRemove: () => void;
}) {
  const meta = CATEGORY_META[op.category];
  const Icon = CATEGORY_ICONS[op.category];
  const isLocked = index === 0;
  const params = pipelineOp.params;

  const modifiedCount = op.params.filter(
    (p) => params[p.name] !== undefined && params[p.name] !== p.default,
  ).length;

  // Check conditions for visible params
  const visibleParams = op.params.filter(
    (p) => !p.condition || p.condition(params),
  );

  return (
    <motion.div
      layout
      initial={{ opacity: 0, scale: 0.92, y: 8 }}
      animate={{ opacity: 1, scale: 1, y: 0 }}
      exit={{ opacity: 0, scale: 0.92, y: -8 }}
      transition={{ ...SPRING, delay: index * STAGGER }}
      className={styles.opCard}
    >
      <button
        type="button"
        onClick={onSelect}
        aria-expanded={isSelected}
        className={`${styles.opCardBtn} ${isSelected ? styles.opCardSelected : ''}`}
      >
        <div
          className={styles.opCardStripe}
          style={{ backgroundColor: meta.color }}
          title={isLocked ? 'Rearrange is required and cannot be removed.' : undefined}
        >
          <Icon size={14} color="white" />
          {isLocked && <Lock size={10} color="rgba(255,255,255,0.7)" />}
        </div>

        <div className={styles.opCardContent}>
          <div className={styles.opCardNameRow}>
            <span className={styles.opCardName}>{op.label}</span>
            {modifiedCount > 0 && (
              <span className={styles.opCardModBadge}>{modifiedCount}</span>
            )}
            {isSelected
              ? <ChevronDown size={12} className={styles.opCardChevron} />
              : <ChevronRight size={12} className={styles.opCardChevron} />}
          </div>
          <span className={styles.opCardDesc}>{op.description}</span>

          {!isSelected && visibleParams.length > 0 && (
            <div className={styles.opCardParamSummary}>
              {visibleParams.slice(0, 3).map((p) => {
                const val = params[p.name] ?? p.default;
                const isModified = val !== p.default;
                return (
                  <span
                    key={p.name}
                    className={isModified ? styles.opCardParamItemModified : styles.opCardParamItem}
                  >
                    {p.label}: {typeof val === 'boolean' ? (val ? 'ON' : 'OFF') : String(val)}
                  </span>
                );
              })}
              {visibleParams.length > 3 && (
                <span className={styles.opCardParamItem}>+{visibleParams.length - 3} more</span>
              )}
            </div>
          )}
        </div>

        {!isLocked && (
          <span
            role="button"
            tabIndex={0}
            aria-label={`Remove ${op.label}`}
            onClick={(e) => { e.stopPropagation(); onRemove(); }}
            onKeyDown={(e) => {
              if (e.key === 'Enter' || e.key === ' ') { e.preventDefault(); e.stopPropagation(); onRemove(); }
            }}
            className={styles.opCardRemove}
          >
            <X size={12} />
          </span>
        )}
      </button>

    </motion.div>
  );
}

// ── FlowArrow ────────────────────────────────────────────

function FlowArrow({ index }: { index: number }) {
  return (
    <motion.div
      initial={{ opacity: 0, scaleX: 0 }}
      animate={{ opacity: 1, scaleX: 1 }}
      exit={{ opacity: 0, scaleX: 0 }}
      transition={{ ...SPRING, delay: index * STAGGER + 0.02 }}
      className={styles.flowArrow}
    >
      <div className={styles.flowArrowLine} />
      <div className={styles.flowArrowHead} />
    </motion.div>
  );
}

// ── AddOpPanel ───────────────────────────────────────────

function AddOpPanel({
  pipelineOpIds,
  onAdd,
}: {
  pipelineOpIds: string[];
  onAdd: (opId: string) => void;
}) {
  const [open, setOpen] = useState(false);
  const btnRef = useRef<HTMLButtonElement>(null);
  const popoverRef = useRef<HTMLDivElement>(null);
  const [pos, setPos] = useState<{ top: number; left: number } | null>(null);

  // Calculate fixed position from button rect
  useEffect(() => {
    if (!open || !btnRef.current) return;
    const rect = btnRef.current.getBoundingClientRect();
    setPos({ top: rect.bottom + 8, left: rect.left });
  }, [open]);

  // Close on click-outside (check both button and portal popover)
  useEffect(() => {
    if (!open) return;
    function handleClick(e: MouseEvent) {
      const target = e.target as Node;
      if (
        btnRef.current && !btnRef.current.contains(target) &&
        popoverRef.current && !popoverRef.current.contains(target)
      ) {
        setOpen(false);
      }
    }
    document.addEventListener('mousedown', handleClick);
    return () => document.removeEventListener('mousedown', handleClick);
  }, [open]);

  // Close on scroll/resize (position would be stale)
  useEffect(() => {
    if (!open) return;
    const close = () => setOpen(false);
    const handleScroll = (e: Event) => {
      // Ignore scroll events originating inside the popover (user scrolling the menu)
      if (popoverRef.current && popoverRef.current.contains(e.target as Node)) return;
      setOpen(false);
    };
    window.addEventListener('resize', close);
    window.addEventListener('scroll', handleScroll, true);
    return () => {
      window.removeEventListener('resize', close);
      window.removeEventListener('scroll', handleScroll, true);
    };
  }, [open]);

  const categories: Category[] = ['rearrangement', 'mutation', 'biology', 'artifacts'];

  const popoverContent = open && pos && ReactDOM.createPortal(
    <div
      ref={popoverRef}
      className={styles.addOpPopover}
      style={{ top: pos.top, left: pos.left }}
    >
      <div className={styles.addOpPopoverHeader}>
        <Plus size={14} style={{ color: 'var(--pg-muted-fg)' }} />
        <span className={styles.addOpPopoverHeaderText}>Add Operation</span>
      </div>
      <div className={styles.addOpPopoverBody}>
        {categories.map((cat) => {
          const ops = getOpsByCategory(cat);
          const catMeta = CATEGORY_META[cat];
          const CatIcon = CATEGORY_ICONS[cat];
          return (
            <div key={cat} className={styles.addOpCatGroup}>
              <div className={styles.addOpCatLabel} style={{ color: catMeta.color }}>
                <CatIcon size={12} />
                <span>{catMeta.label}</span>
              </div>
              {ops.map((op) => {
                const reason = canAddOp(op.id, pipelineOpIds);
                const alreadyIn = pipelineOpIds.includes(op.id);
                const disabled = !!reason || alreadyIn;
                const disabledReason = alreadyIn ? 'Already in pipeline' : reason;
                // Rearrangement ops that aren't already in but would replace the current one
                const isSwap = !alreadyIn && !reason && op.category === 'rearrangement' &&
                  pipelineOpIds.some((id) => {
                    const d = getOpDef(id);
                    return d?.category === 'rearrangement' && d.id !== op.id;
                  });
                return (
                  <button
                    key={op.id}
                    type="button"
                    disabled={disabled}
                    onClick={() => { onAdd(op.id); setOpen(false); }}
                    className={styles.addOpItem}
                  >
                    <div className={styles.addOpItemContent}>
                      <span className={styles.addOpItemName}>{op.label}</span>
                      <p className={styles.addOpItemDesc}>{op.description}</p>
                      {op.params.length > 0 && (
                        <p className={styles.addOpItemParams}>
                          {op.params.length} parameter{op.params.length !== 1 ? 's' : ''}
                        </p>
                      )}
                    </div>
                    {disabled && disabledReason && (
                      <span className={styles.addOpItemReason}>
                        <AlertCircle size={12} />
                        {disabledReason}
                      </span>
                    )}
                    {isSwap && (
                      <span className={styles.addOpItemSwap}>
                        Replaces current rearrangement
                      </span>
                    )}
                  </button>
                );
              })}
            </div>
          );
        })}
      </div>
    </div>,
    document.body,
  );

  return (
    <div className={styles.addOpWrap}>
      <motion.button
        ref={btnRef}
        type="button"
        initial={{ opacity: 0, scale: 0.9 }}
        animate={{ opacity: 1, scale: 1 }}
        className={styles.addOpBtn}
        onClick={() => setOpen(!open)}
        aria-label="Add operation to pipeline"
      >
        <Plus size={16} />
        <span>Add Step</span>
      </motion.button>
      {popoverContent}
    </div>
  );
}

// ── Main Component ───────────────────────────────────────

export default function PipelineBuilder({
  pipeline,
  pipelineOpIds,
  selectedInstanceId,
  onSelectOp,
  onRemoveOp,
  onAddOp,
}: PipelineBuilderProps) {
  const totalModified = pipeline.reduce((acc, pOp) => {
    const def = getOpDef(pOp.opId);
    if (!def) return acc;
    return acc + def.params.filter(
      (p) => pOp.params[p.name] !== undefined && pOp.params[p.name] !== p.default,
    ).length;
  }, 0);

  return (
    <div className={styles.pipelineSection}>
      {/* Section header */}
      <div className={styles.pipelineHeader}>
        <div className={styles.pipelineHeaderDot} />
        <h2 className={styles.pipelineHeaderTitle}>Protocol Pipeline</h2>
        <span className={styles.pipelineHeaderMeta}>
          {pipeline.length} step{pipeline.length !== 1 ? 's' : ''}
        </span>
        {totalModified > 0 && (
          <span className={`${styles.pipelineHeaderMeta} ${styles.pipelineHeaderModified}`}>
            {totalModified} modified
          </span>
        )}
      </div>

      {/* Pipeline flow */}
      <div className={styles.pipelineFlow}>
        <AnimatePresence mode="popLayout">
          {pipeline.map((pOp, i) => {
            const def = getOpDef(pOp.opId);
            if (!def) return null;
            return (
              <React.Fragment key={pOp.instanceId}>
                {i > 0 && <FlowArrow index={i} />}
                <OpCard
                  op={def}
                  pipelineOp={pOp}
                  index={i}
                  isSelected={selectedInstanceId === pOp.instanceId}
                  onSelect={() =>
                    onSelectOp(selectedInstanceId === pOp.instanceId ? null : pOp.instanceId)
                  }
                  onRemove={() => onRemoveOp(pOp.instanceId)}
                />
              </React.Fragment>
            );
          })}
        </AnimatePresence>

        {pipeline.length > 0 && <FlowArrow index={pipeline.length} />}
        <AddOpPanel pipelineOpIds={pipelineOpIds} onAdd={onAddOp} />
      </div>

    </div>
  );
}
