/* ───────────────────────────────────────────────────────────
 *  pyodideWorker.js — Web Worker that runs Pyodide off the
 *  main thread so the UI stays responsive during simulation.
 *
 *  Protocol (main → worker):
 *    { type: 'init',       indexURL, wheelUrl }
 *    { type: 'run',        script }
 *    { type: 'loadConfig', bytes }   (ArrayBuffer, transferred)
 *
 *  Protocol (worker → main):
 *    { type: 'status',      message }
 *    { type: 'ready' }
 *    { type: 'initError',   error }
 *    { type: 'result',      data }   (JSON string)
 *    { type: 'runError',    error }
 *    { type: 'configLoaded' }
 *    { type: 'configError', error }
 * ─────────────────────────────────────────────────────────── */

let pyodide = null;

self.onmessage = async function (e) {
  const msg = e.data;

  if (msg.type === 'init') {
    try {
      // 1. Load Pyodide script
      self.postMessage({ type: 'status', message: 'Loading Python runtime...' });
      importScripts(msg.indexURL + 'pyodide.js');

      // 2. Initialize Pyodide
      self.postMessage({ type: 'status', message: 'Initializing Python...' });
      pyodide = await self.loadPyodide({ indexURL: msg.indexURL });

      // 3. Install GenAIRR wheel
      self.postMessage({ type: 'status', message: 'Installing GenAIRR...' });
      await pyodide.loadPackage('micropip');
      const micropip = pyodide.pyimport('micropip');
      await micropip.install(msg.wheelUrl);

      // 4. Warm up imports
      self.postMessage({ type: 'status', message: 'Warming up...' });
      await pyodide.runPythonAsync(
        'from GenAIRR import simulate; from GenAIRR.protocol import Protocol'
      );

      self.postMessage({ type: 'ready' });
    } catch (err) {
      self.postMessage({ type: 'initError', error: String(err) });
    }
    return;
  }

  if (msg.type === 'run') {
    try {
      const jsonStr = await pyodide.runPythonAsync(msg.script);
      self.postMessage({ type: 'result', data: jsonStr });
    } catch (err) {
      self.postMessage({ type: 'runError', error: String(err) });
    }
    return;
  }

  if (msg.type === 'loadConfig') {
    try {
      const uint8 = new Uint8Array(msg.bytes);
      pyodide.globals.set('_custom_config_bytes', uint8);
      await pyodide.runPythonAsync(
        "import pickle\n_custom_config = pickle.loads(bytes(_custom_config_bytes))"
      );
      self.postMessage({ type: 'configLoaded' });
    } catch (err) {
      self.postMessage({ type: 'configError', error: String(err) });
    }
    return;
  }
};
