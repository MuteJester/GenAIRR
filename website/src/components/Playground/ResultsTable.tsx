/* ───────────────────────────────────────────────────────────
 *  ResultsTable.tsx — Swiss Lab results table with row
 *  selection, column toggle, CSV export, zebra striping.
 *  Clicking a row triggers DissectionBay display.
 * ─────────────────────────────────────────────────────────── */

import React, { useState, useMemo } from 'react';
import { Download, Columns3 } from 'lucide-react';
import type { SequenceRecord } from './types';
import styles from './styles.module.css';

interface ResultsTableProps {
  results: SequenceRecord[];
  selectedId: number | null;
  onSelect: (id: number) => void;
}

const ALL_COLUMNS = [
  { key: 'sequence_id' as const, label: 'SEQUENCE ID' },
  { key: 'v_call' as const, label: 'V CALL' },
  { key: 'd_call' as const, label: 'D CALL' },
  { key: 'j_call' as const, label: 'J CALL' },
  { key: 'junction_aa' as const, label: 'JUNCTION' },
  { key: 'productive' as const, label: 'PRODUCTIVE' },
  { key: 'mutation_count' as const, label: 'MUTATIONS' },
  { key: 'cdr3_length' as const, label: 'CDR3 LEN' },
] as const;

const MONO_KEYS = new Set(['junction_aa', 'sequence_id']);

export default function ResultsTable({ results, selectedId, onSelect }: ResultsTableProps) {
  const [showAll, setShowAll] = useState(false);
  const columns = showAll ? ALL_COLUMNS : ALL_COLUMNS.slice(0, 6);

  const allRawColumns = useMemo(() => {
    if (!results.length) return [] as string[];
    return Object.keys(results[0]);
  }, [results]);

  function downloadCSV() {
    if (!results.length) return;
    const headers = allRawColumns;
    const rows = results.map((row) =>
      headers.map((h) => {
        const v = (row as any)[h] ?? '';
        const s = typeof v === 'object' ? JSON.stringify(v) : String(v);
        return s.includes(',') || s.includes('"') || s.includes('\n')
          ? `"${s.replace(/"/g, '""')}"` : s;
      }).join(','),
    );
    const csv = [headers.join(','), ...rows].join('\n');
    const blob = new Blob([csv], { type: 'text/csv' });
    const url = URL.createObjectURL(blob);
    const a = document.createElement('a');
    a.href = url;
    a.download = 'genairr_simulation.csv';
    a.click();
    URL.revokeObjectURL(url);
  }

  return (
    <div className={styles.resultsSection}>
      {/* Section header */}
      <div className={styles.resultsHeader}>
        <h2 className={styles.resultsTitle}>Results</h2>
        <span className={styles.resultsCount}>
          {results.length} sequence{results.length !== 1 ? 's' : ''}
        </span>

        {selectedId != null && (
          <span className={styles.resultsSelectedHint}>
            SEQ_{String(selectedId).padStart(4, '0')} selected — scroll up to view dissection
          </span>
        )}

        <div className={styles.resultsActions}>
          <button
            type="button"
            className={styles.resultsBtn}
            onClick={() => setShowAll(!showAll)}
            title={showAll ? 'Show key columns' : 'Show all columns'}
          >
            <Columns3 size={13} />
            {showAll ? 'Fewer columns' : 'All columns'}
          </button>
          <button type="button" className={styles.resultsBtn} onClick={downloadCSV}>
            <Download size={13} /> CSV
          </button>
        </div>
      </div>

      {/* Table */}
      <div className={styles.tableWrapper}>
        <table className={styles.table}>
          <thead>
            <tr>
              <th>#</th>
              {columns.map((col) => (
                <th key={col.key}>{col.label}</th>
              ))}
            </tr>
          </thead>
          <tbody>
            {results.map((r, idx) => {
              const isSelected = selectedId === r.id;
              const rowClass = [
                styles.tableRow,
                isSelected
                  ? styles.tableRowSelected
                  : idx % 2 === 0
                    ? styles.tableRowEven
                    : styles.tableRowOdd,
              ].join(' ');

              return (
                <tr
                  key={r.id}
                  className={rowClass}
                  onClick={() => onSelect(r.id)}
                  aria-selected={isSelected}
                  role="row"
                >
                  <td className={styles.tdIndex}>{r.id}</td>
                  {columns.map((col) => {
                    let val: string;
                    if (col.key === 'productive') {
                      val = r.productive ? 'Yes' : 'No';
                    } else {
                      val = String(r[col.key]);
                    }
                    const isMono = MONO_KEYS.has(col.key);
                    const isAccent = isSelected && col.key === 'sequence_id';
                    const classes = [
                      isMono ? styles.tdMono : '',
                      isAccent ? styles.tdSelectedAccent : '',
                    ].filter(Boolean).join(' ') || undefined;

                    return (
                      <td key={col.key} className={classes}>
                        {val}
                      </td>
                    );
                  })}
                </tr>
              );
            })}
          </tbody>
        </table>
      </div>
    </div>
  );
}
