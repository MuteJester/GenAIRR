/* ───────────────────────────────────────────────────────────
 *  usePyodide.ts — React hook that manages the Pyodide Web
 *  Worker. All Python execution happens off the main thread
 *  so the UI stays responsive during simulation.
 * ─────────────────────────────────────────────────────────── */

import { useState, useRef, useEffect, useCallback } from 'react';
import useDocusaurusContext from '@docusaurus/useDocusaurusContext';
import { buildPythonScript, type RunConfig } from './pythonRunner';

const PYODIDE_CDN = 'https://cdn.jsdelivr.net/pyodide/v0.27.3/full/';
const WHEEL_FILENAME = 'genairr-1.0.0-py3-none-any.whl';

export type PyodideState = 'idle' | 'loading' | 'ready' | 'error';
export type SimState = 'idle' | 'running' | 'done' | 'error';

export interface UsePyodideReturn {
  pyodideState: PyodideState;
  loadingMessage: string;
  simState: SimState;
  results: Record<string, any>[] | null;
  error: string | null;
  runSimulation: (cfg: RunConfig) => Promise<void>;
  loadCustomConfig: (fileBytes: ArrayBuffer, fileName: string) => Promise<void>;
  customConfigName: string | null;
  generatedCode: string | null;
}

export function usePyodide(): UsePyodideReturn {
  const { siteConfig } = useDocusaurusContext();
  const baseUrl = siteConfig.baseUrl;
  const workerRef = useRef<Worker | null>(null);
  const [pyodideState, setPyodideState] = useState<PyodideState>('idle');
  const [loadingMessage, setLoadingMessage] = useState('');
  const [simState, setSimState] = useState<SimState>('idle');
  const [results, setResults] = useState<Record<string, any>[] | null>(null);
  const [error, setError] = useState<string | null>(null);
  const [customConfigName, setCustomConfigName] = useState<string | null>(null);
  const [generatedCode, setGeneratedCode] = useState<string | null>(null);

  // Pending promise resolvers for run / loadConfig (one at a time)
  const pendingRun = useRef<{ resolve: () => void; reject: (e: Error) => void } | null>(null);
  const pendingConfig = useRef<{ resolve: () => void; reject: (e: Error) => void } | null>(null);

  // ── Spawn worker & bootstrap Pyodide ─────────────────────

  useEffect(() => {
    let cancelled = false;

    setPyodideState('loading');
    setLoadingMessage('Loading Python runtime...');

    const workerUrl = `${baseUrl}js/pyodideWorker.js`;
    const worker = new Worker(workerUrl);
    workerRef.current = worker;

    worker.onmessage = (e) => {
      if (cancelled) return;
      const msg = e.data;

      switch (msg.type) {
        // ── Init lifecycle ─────────────────────────────────
        case 'status':
          setLoadingMessage(msg.message);
          break;
        case 'ready':
          setPyodideState('ready');
          setLoadingMessage('');
          break;
        case 'initError':
          setPyodideState('error');
          setError(msg.error);
          setLoadingMessage('');
          break;

        // ── Simulation results ─────────────────────────────
        case 'result': {
          const parsed: Record<string, any>[] = JSON.parse(msg.data);
          setResults(parsed);
          setSimState('done');
          pendingRun.current?.resolve();
          pendingRun.current = null;
          break;
        }
        case 'runError':
          setError(msg.error);
          setSimState('error');
          pendingRun.current?.reject(new Error(msg.error));
          pendingRun.current = null;
          break;

        // ── Config upload ──────────────────────────────────
        case 'configLoaded':
          pendingConfig.current?.resolve();
          pendingConfig.current = null;
          break;
        case 'configError':
          setError(`Failed to load config: ${msg.error}`);
          pendingConfig.current?.reject(new Error(msg.error));
          pendingConfig.current = null;
          break;
      }
    };

    // Send init message
    const wheelUrl = `${window.location.origin}${baseUrl}whl/${WHEEL_FILENAME}`;
    worker.postMessage({ type: 'init', indexURL: PYODIDE_CDN, wheelUrl });

    return () => {
      cancelled = true;
      worker.terminate();
      workerRef.current = null;
    };
  }, [baseUrl]);

  // ── Run Simulation ────────────────────────────────────────

  const runSimulation = useCallback(async (cfg: RunConfig) => {
    const worker = workerRef.current;
    if (!worker) return;

    setSimState('running');
    setError(null);
    setResults(null);

    const script = buildPythonScript(cfg);
    setGeneratedCode(script);

    return new Promise<void>((resolve, reject) => {
      pendingRun.current = { resolve, reject };
      worker.postMessage({ type: 'run', script });
    });
  }, []);

  // ── Load Custom Config ────────────────────────────────────

  const loadCustomConfig = useCallback(async (fileBytes: ArrayBuffer, fileName: string) => {
    const worker = workerRef.current;
    if (!worker) return;

    return new Promise<void>((resolve, reject) => {
      pendingConfig.current = {
        resolve: () => {
          setCustomConfigName(fileName);
          setError(null);
          resolve();
        },
        reject: (e) => {
          setCustomConfigName(null);
          reject(e);
        },
      };
      // Transfer the ArrayBuffer to avoid copying
      worker.postMessage({ type: 'loadConfig', bytes: fileBytes }, [fileBytes]);
    });
  }, []);

  return {
    pyodideState,
    loadingMessage,
    simState,
    results,
    error,
    runSimulation,
    loadCustomConfig,
    customConfigName,
    generatedCode,
  };
}
