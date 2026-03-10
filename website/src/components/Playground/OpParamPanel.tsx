/* ───────────────────────────────────────────────────────────
 *  OpParamPanel.tsx — Right-side parameter editor panel.
 *  Anchored to results for spatial continuity.
 * ─────────────────────────────────────────────────────────── */

import React from 'react';
import { SlidersHorizontal } from 'lucide-react';
import { CATEGORY_META, getOpDef, type OpDef } from './ops';
import type { PipelineOp } from './pythonRunner';
import styles from './styles.module.css';

interface Props {
  pipeline: PipelineOp[];
  selectedInstanceId: string | null;
  onUpdateParams: (id: string, params: Record<string, any>) => void;
}

export default function OpParamPanel({ pipeline, selectedInstanceId, onUpdateParams }: Props) {
  const selected = pipeline.find((p) => p.instanceId === selectedInstanceId) ?? null;
  const def = selected ? getOpDef(selected.opId) : null;

  if (!selected || !def) {
    return (
      <aside className={styles.paramPanel} aria-live="polite">
        <div className={styles.paramPanelHeader}>
          <div className={styles.paramPanelTitle}>
            <SlidersHorizontal size={12} /> Parameters
          </div>
        </div>
        <div className={styles.paramPanelBody}>
          <div className={styles.paramPanelEmpty}>
            Select an operation in the pipeline to edit parameters.
          </div>
        </div>
      </aside>
    );
  }

  const params = selected.params;
  const visibleParams = def.params.filter((p) => !p.condition || p.condition(params));
  const modifiedCount = def.params.filter(
    (p) => params[p.name] !== undefined && params[p.name] !== p.default,
  ).length;

  const meta = CATEGORY_META[def.category];

  return (
    <aside className={styles.paramPanel} aria-live="polite">
      <div className={styles.paramPanelHeader}>
        <div className={styles.paramPanelTitle}>
          <SlidersHorizontal size={12} /> Parameters
        </div>
        <span className={styles.paramPanelMeta}>{modifiedCount} modified</span>
        <span className={styles.paramPanelCat} style={{ background: meta.color }}>
          {meta.label}
        </span>
      </div>
      <div className={styles.paramPanelBody}>
        <div className={styles.paramPanelOpName}>{def.label}</div>
        <p className={styles.paramPanelOpDesc}>{def.description}</p>

        {visibleParams.length === 0 ? (
          <div className={styles.paramPanelEmpty}>
            This operation has no editable parameters.
          </div>
        ) : (
          <div className={styles.paramPanelFields}>
            {visibleParams.map((p) => (
              <ParamField
                key={p.name}
                param={p}
                value={params[p.name] ?? p.default}
                onChange={(v) => onUpdateParams(selected.instanceId, { ...params, [p.name]: v })}
              />
            ))}
          </div>
        )}
      </div>
    </aside>
  );
}

function ParamField({
  param,
  value,
  onChange,
}: {
  param: OpDef['params'][0];
  value: any;
  onChange: (v: any) => void;
}) {
  switch (param.type) {
    case 'number':
      return (
        <div className={styles.paramField}>
          <div className={styles.paramFieldHeader}>
            <span className={styles.paramFieldLabel}>{param.label}</span>
            <span className={styles.paramFieldDesc}>{param.description}</span>
          </div>
          <div style={{ display: 'flex', alignItems: 'center', gap: '0.5rem' }}>
            <input
              type="number"
              value={value as number}
              min={param.min}
              max={param.max}
              step={param.step ?? 1}
              onChange={(e) => onChange(Number(e.target.value))}
              className={styles.paramInput}
            />
            {param.min !== undefined && param.max !== undefined && (
              <span className={styles.paramInputRange}>{param.min}..{param.max}</span>
            )}
          </div>
        </div>
      );
    case 'select':
      return (
        <div className={styles.paramField}>
          <div className={styles.paramFieldHeader}>
            <span className={styles.paramFieldLabel}>{param.label}</span>
            <span className={styles.paramFieldDesc}>{param.description}</span>
          </div>
          <select
            value={value as string}
            onChange={(e) => onChange(e.target.value)}
            className={styles.paramSelect}
          >
            {param.options?.map((opt) => (
              <option key={opt.value} value={opt.value}>{opt.label}</option>
            ))}
          </select>
        </div>
      );
    case 'range': {
      const numVal = value as number;
      const pctFill =
        param.min !== undefined && param.max !== undefined
          ? ((numVal - param.min) / (param.max - param.min)) * 100
          : 50;
      const decimals = param.step && param.step < 0.01 ? 3 : param.step && param.step < 0.1 ? 2 : 1;
      return (
        <div className={styles.paramField}>
          <div className={styles.paramFieldHeader}>
            <span className={styles.paramFieldLabel}>{param.label}</span>
            <span className={styles.paramFieldDesc}>{param.description}</span>
          </div>
          <div className={styles.rangeWrap}>
            <div className={styles.rangeTrack}>
              <div className={styles.rangeTrackBg} />
              <div className={styles.rangeTrackFill} style={{ width: `${pctFill}%` }} />
              <input
                type="range"
                min={param.min}
                max={param.max}
                step={param.step}
                value={numVal}
                onChange={(e) => onChange(Number(e.target.value))}
                className={styles.rangeInput}
              />
            </div>
            <span className={styles.rangeValue}>{numVal.toFixed(decimals)}</span>
          </div>
        </div>
      );
    }
    case 'boolean':
      return (
        <div className={styles.paramField}>
          <div className={styles.paramFieldHeader}>
            <span className={styles.paramFieldLabel}>{param.label}</span>
            <span className={styles.paramFieldDesc}>{param.description}</span>
          </div>
          <button
            type="button"
            role="switch"
            aria-checked={value as boolean}
            onClick={() => onChange(!(value as boolean))}
            className={`${styles.toggleSwitch} ${value ? styles.toggleSwitchOn : ''}`}
          >
            <div
              style={{
                position: 'absolute',
                width: 12,
                height: 12,
                background: 'var(--pg-background)',
                left: value ? 19 : 3,
                transition: 'left 0.2s ease',
              }}
            />
          </button>
        </div>
      );
    default:
      return null;
  }
}
