/* ───────────────────────────────────────────────────────────
 *  CodePreview.tsx — Generated Python preview for Playground.
 *  Extracted from PipelineBuilder to allow RunBar grouping.
 * ─────────────────────────────────────────────────────────── */

import React, { useState } from 'react';
import { motion } from 'framer-motion';
import { Code, Copy, Check } from 'lucide-react';
import styles from './styles.module.css';

export default function CodePreview({ code }: { code: string }) {
  const [copied, setCopied] = useState(false);

  const handleCopy = () => {
    navigator.clipboard.writeText(code);
    setCopied(true);
    setTimeout(() => setCopied(false), 2000);
  };

  return (
    <motion.div
      initial={{ height: 0, opacity: 0 }}
      animate={{ height: 'auto', opacity: 1 }}
      exit={{ height: 0, opacity: 0 }}
      transition={{ duration: 0.25 }}
      className={styles.codePreview}
    >
      <div className={styles.codePreviewInner}>
        <div className={styles.codePreviewHeader}>
          <div className={styles.codePreviewLabel}>
            <Code size={12} />
            <span>Generated Python</span>
          </div>
          <button type="button" onClick={handleCopy} className={styles.codePreviewCopy}>
            {copied ? <Check size={12} /> : <Copy size={12} />}
            {copied ? 'Copied' : 'Copy'}
          </button>
        </div>
        <pre className={styles.codePreviewPre}>
          <code>
            {code.split('\n').map((line, i) => (
              <div key={i} style={{ display: 'flex' }}>
                <span className={styles.codePreviewLineNum}>{i + 1}</span>
                <span>
                  {highlightLine(line)}
                </span>
              </div>
            ))}
          </code>
        </pre>
      </div>
    </motion.div>
  );
}

function highlightLine(line: string): React.ReactNode {
  const commentIdx = line.indexOf('#');
  if (commentIdx >= 0) {
    const before = line.slice(0, commentIdx);
    const comment = line.slice(commentIdx);
    return (
      <>
        {before && highlightSegment(before)}
        <span className={styles.codeComment}>{comment}</span>
      </>
    );
  }
  return highlightSegment(line);
}

function highlightSegment(text: string): React.ReactNode {
  const parts = text.split(/(from|import|with|as|for|in|def|class|True|False|None|"[^"]*"|'[^']*')/g);
  return parts.map((part, j) => {
    if (['from', 'import', 'with', 'as', 'for', 'in', 'def', 'class'].includes(part)) {
      return <span key={j} className={styles.codeKw}>{part}</span>;
    }
    if (['True', 'False', 'None'].includes(part)) {
      return <span key={j} className={styles.codeBool}>{part}</span>;
    }
    if (part.startsWith('"') || part.startsWith("'")) {
      return <span key={j} className={styles.codeStr}>{part}</span>;
    }
    return <span key={j}>{part}</span>;
  });
}
