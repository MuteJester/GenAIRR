/* ───────────────────────────────────────────────────────────
 *  RunBar.tsx — Compact run controls bar. Swiss Lab aesthetic.
 *  Config, upload, n, seed, productive toggle, run button.
 * ─────────────────────────────────────────────────────────── */

import React, { useRef } from 'react';
import { Play, Upload, Code } from 'lucide-react';
import type { SimState, PyodideState } from './usePyodide';
import styles from './styles.module.css';

interface Props {
  pyodideState: PyodideState;
  simState: SimState;
  n: number;
  seed: number | null;
  productive: boolean;
  useCustomConfig: boolean;
  customConfigName: string | null;
  onChangeN: (n: number) => void;
  onChangeSeed: (seed: number | null) => void;
  onChangeProductive: (v: boolean) => void;
  onRun: () => void;
  onUploadConfig: (bytes: ArrayBuffer, name: string) => void;
  showCode: boolean;
  onToggleCode: () => void;
}

export default function RunBar({
  pyodideState, simState,
  n, seed, productive,
  useCustomConfig, customConfigName,
  onChangeN, onChangeSeed, onChangeProductive,
  onRun, onUploadConfig, showCode, onToggleCode,
}: Props) {
  const fileRef = useRef<HTMLInputElement>(null);
  const ready = pyodideState === 'ready';
  const running = simState === 'running';

  function handleFile(e: React.ChangeEvent<HTMLInputElement>) {
    const file = e.target.files?.[0];
    if (!file) return;
    const reader = new FileReader();
    reader.onload = () => {
      if (reader.result instanceof ArrayBuffer) {
        onUploadConfig(reader.result, file.name);
      }
    };
    reader.readAsArrayBuffer(file);
    e.target.value = '';
  }

  return (
    <div className={styles.runBar}>
      {/* Config group */}
      <div className={styles.runBarGroup}>
        <span className={styles.runBarLabel}>Config</span>
        <span className={styles.runBarConfigName}>
          {useCustomConfig && customConfigName
            ? customConfigName
            : 'Human IGH (built-in)'}
        </span>
        <button
          type="button"
          className={styles.runBarUpload}
          onClick={() => fileRef.current?.click()}
          disabled={!ready}
          title="Upload your own DataConfig .pkl"
        >
          <Upload size={13} /> Upload .pkl
        </button>
        <input
          ref={fileRef}
          type="file"
          accept=".pkl,.pickle"
          style={{ display: 'none' }}
          onChange={handleFile}
        />
      </div>

      {/* n */}
      <div className={styles.runBarGroup}>
        <span className={styles.runBarLabel}>Sequences</span>
        <input
          type="number"
          className={styles.runBarInput}
          value={n}
          min={1}
          max={500}
          onChange={(e) => onChangeN(Math.max(1, Math.min(500, parseInt(e.target.value) || 1)))}
        />
      </div>

      {/* Seed */}
      <div className={styles.runBarGroup}>
        <span className={styles.runBarLabel}>Seed</span>
        <input
          type="number"
          className={styles.runBarInput}
          value={seed ?? ''}
          placeholder="None"
          onChange={(e) => onChangeSeed(e.target.value === '' ? null : parseInt(e.target.value))}
        />
      </div>

      {/* Productive */}
      <label className={styles.runBarCheck}>
        <input
          type="checkbox"
          checked={productive}
          onChange={(e) => onChangeProductive(e.target.checked)}
        />
        <span>Productive only</span>
      </label>

      <div className={styles.runBarActions}>
        <button
          type="button"
          className={`${styles.runBarCodeToggle} ${showCode ? styles.runBarCodeToggleActive : ''}`}
          onClick={onToggleCode}
          aria-pressed={showCode}
          title="Toggle generated Python"
        >
          <Code size={13} /> Python
        </button>
        <button
          type="button"
          className={styles.runButton}
          onClick={onRun}
          disabled={!ready || running}
          title="Run (Ctrl/Cmd+Enter)"
        >
          {running ? (
            <>
              <span className={styles.runButtonSpinner} />
              Running...
            </>
          ) : (
            <>
              <Play size={15} /> Run Simulation
            </>
          )}
        </button>
      </div>
    </div>
  );
}
