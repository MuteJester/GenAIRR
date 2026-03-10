/**
 * trace.c — Simulation trace log implementation.
 */

#include "genairr/trace.h"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

#define TRACE_INITIAL_CAP  8192
#define TRACE_LINE_MAX     1024

/* ── Thread-local trace pointer ───────────────────────────────── */

static _Thread_local TraceLog *g_trace = NULL;

void trace_set_current(TraceLog *log) {
    g_trace = log;
}

TraceLog *trace_current(void) {
    return g_trace;
}

/* ── TraceLog lifecycle ───────────────────────────────────────── */

void trace_log_init(TraceLog *log) {
    log->cap = TRACE_INITIAL_CAP;
    log->buf = (char *)malloc(log->cap);
    log->buf[0] = '\0';
    log->len = 0;
}

void trace_log_destroy(TraceLog *log) {
    free(log->buf);
    log->buf = NULL;
    log->len = 0;
    log->cap = 0;
}

void trace_log_clear(TraceLog *log) {
    log->len = 0;
    if (log->buf) log->buf[0] = '\0';
}

/* ── Append formatted text ────────────────────────────────────── */

void trace_log_appendf(TraceLog *log, const char *fmt, ...) {
    if (!log || !log->buf) return;

    char line[TRACE_LINE_MAX];
    va_list ap;
    va_start(ap, fmt);
    int n = vsnprintf(line, sizeof(line), fmt, ap);
    va_end(ap);

    if (n < 0) return;
    if (n >= (int)sizeof(line)) n = (int)sizeof(line) - 1;

    /* +1 for newline, +1 for NUL */
    int needed = log->len + n + 2;
    if (needed > log->cap) {
        int new_cap = log->cap;
        while (new_cap < needed) new_cap *= 2;
        log->buf = (char *)realloc(log->buf, new_cap);
        log->cap = new_cap;
    }

    memcpy(log->buf + log->len, line, n);
    log->len += n;
    log->buf[log->len++] = '\n';
    log->buf[log->len] = '\0';
}
