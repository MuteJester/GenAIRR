/* ───────────────────────────────────────────────────────────
 *  index.tsx — GenAIRR Playground page component.
 *  Swiss Lab layout: PipelineBuilder → RunBar → CodePreview
 *  → Results + Param Panel, with 4-step loading overlay.
 * ─────────────────────────────────────────────────────────── */

import React, { useState, useCallback, useMemo, useEffect } from 'react';
import { AnimatePresence } from 'framer-motion';
import { AlertTriangle } from 'lucide-react';
import { usePyodide } from './usePyodide';
import { defaultParams, getOpDef } from './ops';
import { buildDisplayScript, type PipelineOp } from './pythonRunner';
import { parseAIRRResults, type SequenceRecord } from './types';
import PipelineBuilder from './PipelineBuilder';
import RunBar from './RunBar';
import DissectionBay from './DissectionBay';
import ResultsTable from './ResultsTable';
import EmptyState from './EmptyState';
import CodePreview from './CodePreview';
import OpParamPanel from './OpParamPanel';
import styles from './styles.module.css';

// ── Loading step mapping ──────────────────────────────────

const LOADING_STEPS = [
  'Loading Python runtime',
  'Initializing Python',
  'Installing GenAIRR',
  'Warming up',
];

function loadingStep(message: string): number {
  if (message.includes('Warming')) return 3;
  if (message.includes('Installing')) return 2;
  if (message.includes('Initializing')) return 1;
  return 0;
}

function loadingProgress(step: number): number {
  return Math.min(((step + 1) / LOADING_STEPS.length) * 100, 100);
}

// ── UID generator ─────────────────────────────────────────

let _nextId = 1;
function uid(): string {
  return `op_${_nextId++}`;
}

// ── Main component ────────────────────────────────────────

export default function Playground() {
  // ── Pyodide ────────────────────────────────────────────
  const {
    pyodideState, loadingMessage, simState,
    results: rawResults, error,
    runSimulation, loadCustomConfig, customConfigName,
  } = usePyodide();

  // ── Pipeline state ─────────────────────────────────────
  const [pipeline, setPipeline] = useState<PipelineOp[]>([
    { instanceId: uid(), opId: 'Rearrange', params: defaultParams('Rearrange') },
  ]);
  const [selectedInstanceId, setSelectedInstanceId] = useState<string | null>(null);
  const [showCode, setShowCode] = useState(false);

  // ── Run config ─────────────────────────────────────────
  const [n, setN] = useState(10);
  const [seed, setSeed] = useState<number | null>(42);
  const [productive, setProductive] = useState(false);
  const [useCustom, setUseCustom] = useState(false);

  // ── Results + selection ────────────────────────────────
  const [selectedId, setSelectedId] = useState<number | null>(null);

  const results: SequenceRecord[] | null = useMemo(() => {
    if (!rawResults || !rawResults.length) return null;
    return parseAIRRResults(rawResults);
  }, [rawResults]);

  const selectedResult = useMemo(
    () => results?.find((r) => r.id === selectedId) ?? null,
    [results, selectedId],
  );

  const pipelineOpIds = useMemo(() => pipeline.map((p) => p.opId), [pipeline]);

  // ── Pipeline mutations ─────────────────────────────────

  const addOp = useCallback((opId: string) => {
    const def = getOpDef(opId);
    if (!def) return;

    if (def.category === 'rearrangement') {
      setPipeline((prev) => [
        { instanceId: uid(), opId, params: defaultParams(opId) },
        ...prev.slice(1),
      ]);
      return;
    }

    setPipeline((prev) => [
      ...prev,
      { instanceId: uid(), opId, params: defaultParams(opId) },
    ]);
  }, []);

  const removeOp = useCallback((instanceId: string) => {
    setPipeline((prev) => prev.filter((op) => op.instanceId !== instanceId));
    setSelectedInstanceId((cur) => (cur === instanceId ? null : cur));
  }, []);

  const updateParams = useCallback((instanceId: string, params: Record<string, any>) => {
    setPipeline((prev) =>
      prev.map((op) => (op.instanceId === instanceId ? { ...op, params } : op)),
    );
  }, []);

  // ── Upload ─────────────────────────────────────────────

  const handleUpload = useCallback(async (bytes: ArrayBuffer, name: string) => {
    await loadCustomConfig(bytes, name);
    setUseCustom(true);
  }, [loadCustomConfig]);

  // ── Generated code for display ─────────────────────────

  const displayCode = useMemo(
    () => buildDisplayScript({ pipeline, useCustomConfig: useCustom, n, seed, productive }),
    [pipeline, useCustom, n, seed, productive],
  );

  // ── Run ────────────────────────────────────────────────

  const handleRun = useCallback(() => {
    setSelectedId(null);
    runSimulation({ pipeline, useCustomConfig: useCustom, n, seed, productive });
  }, [runSimulation, pipeline, useCustom, n, seed, productive]);

  // ── Keyboard shortcut: Cmd/Ctrl + Enter ────────────────

  useEffect(() => {
    const handler = (e: KeyboardEvent) => {
      if (!(e.metaKey || e.ctrlKey) || e.key !== 'Enter') return;
      const target = e.target as HTMLElement | null;
      if (target) {
        const tag = target.tagName.toLowerCase();
        if (['input', 'textarea', 'select'].includes(tag) || target.isContentEditable) {
          return;
        }
      }
      if (pyodideState !== 'ready' || simState === 'running') return;
      e.preventDefault();
      handleRun();
    };
    window.addEventListener('keydown', handler);
    return () => window.removeEventListener('keydown', handler);
  }, [handleRun, pyodideState, simState]);

  // ── Result selection ───────────────────────────────────

  const handleSelect = useCallback((id: number) => {
    setSelectedId((prev) => (prev === id ? null : id));
  }, []);

  const handleCloseDissection = useCallback(() => {
    setSelectedId(null);
  }, []);

  // ── Loading overlay ────────────────────────────────────

  const currentStep = loadingStep(loadingMessage);
  const progress = loadingProgress(currentStep);

  // ── Render ─────────────────────────────────────────────

  return (
    <div className={styles.playground}>
      {/* 4-step loading overlay */}
      {(pyodideState === 'idle' || pyodideState === 'loading') && (
        <div className={styles.loadingOverlay} role="alert" aria-live="polite">
          <div className={styles.loadingInner}>
            {/* Title */}
            <div style={{ textAlign: 'center' }}>
              <h1 className={styles.loadingTitle}>GenAIRR Playground</h1>
              <p className={styles.loadingSubtitle}>Synthetic Antibody Repertoire Generator</p>
            </div>

            {/* Progress bar */}
            <div style={{ width: '100%' }}>
              <div className={styles.loadingProgressTrack}>
                <div className={styles.loadingProgressFill} style={{ width: `${progress}%` }} />
              </div>
              <div className={styles.loadingProgressInfo}>
                <span className={styles.loadingStepLabel}>
                  {LOADING_STEPS[currentStep]}...
                </span>
                <span className={styles.loadingPct}>{Math.round(progress)}%</span>
              </div>
            </div>

            {/* Step indicators */}
            <div className={styles.loadingSteps}>
              {LOADING_STEPS.map((s, i) => (
                <div key={s} className={styles.loadingStep}>
                  <div className={`${styles.loadingStepBar} ${i <= currentStep ? styles.loadingStepBarActive : ''}`} />
                  <span className={`${styles.loadingStepText} ${i <= currentStep ? styles.loadingStepTextActive : ''}`}>
                    {s}
                  </span>
                </div>
              ))}
            </div>

            {/* Hint */}
            <p className={styles.loadingHint}>
              First visit loads the Python runtime — subsequent visits are much faster.
            </p>
          </div>
        </div>
      )}

      {/* Init error */}
      {pyodideState === 'error' && (
        <div className={styles.loadingOverlay}>
          <div className={styles.loadingInner}>
            <AlertTriangle size={36} style={{ color: 'var(--pg-destructive)' }} />
            <h1 className={styles.loadingTitle}>Failed to initialize Python runtime</h1>
            <p className={styles.errorText}>{error}</p>
            <button className={styles.retryBtn} onClick={() => window.location.reload()}>
              Retry
            </button>
          </div>
        </div>
      )}

      {/* Main content */}
      {pyodideState === 'ready' && (
        <div className={styles.playgroundMain}>
          {/* 1. Pipeline builder */}
          <PipelineBuilder
            pipeline={pipeline}
            pipelineOpIds={pipelineOpIds}
            selectedInstanceId={selectedInstanceId}
            onSelectOp={setSelectedInstanceId}
            onRemoveOp={removeOp}
            onAddOp={addOp}
          />

          {/* 2. Run bar */}
          <RunBar
            pyodideState={pyodideState}
            simState={simState}
            n={n}
            seed={seed}
            productive={productive}
            useCustomConfig={useCustom}
            customConfigName={customConfigName}
            onChangeN={setN}
            onChangeSeed={setSeed}
            onChangeProductive={setProductive}
            onRun={handleRun}
            onUploadConfig={handleUpload}
            showCode={showCode}
            onToggleCode={() => setShowCode((v) => !v)}
          />

          <AnimatePresence>
            {showCode && <CodePreview code={displayCode} />}
          </AnimatePresence>

          <div className={styles.playgroundBody}>
            <div className={styles.resultsColumn}>
              {/* Simulation error */}
              {simState === 'error' && error && (
                <div className={styles.simError}>
                  <AlertTriangle size={15} />
                  <pre>{error}</pre>
                </div>
              )}

              {/* Dissection bay (when a sequence is selected) */}
              {selectedResult && (
                <DissectionBay
                  result={selectedResult}
                  onClose={handleCloseDissection}
                />
              )}

              {/* Results / empty / running */}
              {simState === 'running' ? (
                <div className={styles.runningState}>
                  <div className={styles.runningStateInner}>
                    <div className={styles.runningBar}>
                      <div className={styles.runningBarFill} />
                    </div>
                    <span className={styles.runningText}>Generating sequences...</span>
                  </div>
                </div>
              ) : results && results.length > 0 ? (
                <ResultsTable
                  results={results}
                  selectedId={selectedId}
                  onSelect={handleSelect}
                />
              ) : simState === 'idle' || simState === 'done' ? (
                <EmptyState />
              ) : null}
            </div>

            <div className={styles.paramColumn}>
              <OpParamPanel
                pipeline={pipeline}
                selectedInstanceId={selectedInstanceId}
                onUpdateParams={updateParams}
              />
            </div>
          </div>
        </div>
      )}
    </div>
  );
}
