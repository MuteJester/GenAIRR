/**
 * trace.h — Simulation trace log for live narration.
 *
 * Provides a growable text buffer that pipeline steps can write to
 * via the TRACE() macro. When tracing is disabled (the common case),
 * the macro compiles to a single branch that's never taken — zero
 * overhead on the hot path.
 *
 * Thread-local pointer: set before pipeline execution, cleared after.
 * Step functions access it via trace_current() without signature changes.
 */

#ifndef GENAIRR_TRACE_H
#define GENAIRR_TRACE_H

#include <stdbool.h>
#include <stdarg.h>

/* ── TraceLog: growable text buffer ───────────────────────────── */

typedef struct {
    char  *buf;       /* heap-allocated text buffer */
    int    len;       /* current length (excluding NUL) */
    int    cap;       /* allocated capacity */
} TraceLog;

/** Initialize a trace log (starts with 8KB). */
void  trace_log_init(TraceLog *log);

/** Free the trace log buffer. */
void  trace_log_destroy(TraceLog *log);

/** Clear the buffer (keep allocation). */
void  trace_log_clear(TraceLog *log);

/** Append a formatted line to the trace log.
 *  Grows the buffer as needed. */
void  trace_log_appendf(TraceLog *log, const char *fmt, ...)
#if defined(__GNUC__) || defined(__clang__)
    __attribute__((format(printf, 2, 3)))
#endif
    ;

/* ── Thread-local current trace ───────────────────────────────── */

/** Set the current trace log for this thread.
 *  Pass NULL to disable tracing. */
void  trace_set_current(TraceLog *log);

/** Get the current trace log (NULL if tracing disabled). */
TraceLog *trace_current(void);

/* ── TRACE macro ──────────────────────────────────────────────── */

/**
 * TRACE(fmt, ...) — write a line to the current trace log.
 *
 * Compiles to a single pointer check when tracing is off.
 * When on, appends a formatted line to the trace buffer.
 *
 * Usage:
 *   TRACE("[sample_v] picked %s from pool of %d", name, count);
 */
#define TRACE(fmt, ...) \
    do { \
        TraceLog *_tl = trace_current(); \
        if (_tl) trace_log_appendf(_tl, fmt, ##__VA_ARGS__); \
    } while (0)

#endif /* GENAIRR_TRACE_H */
