/**
 * bench_e2e.c — End-to-end benchmark: rearrange + S5F + serialize.
 *
 * Usage: ./bench_e2e [N]
 *   N = number of sequences to generate (default 10000)
 *
 * Measures wall-clock time for:
 *   Phase 1: Config + correction context setup
 *   Phase 2: N × (pipeline_execute + s5f_mutate + airr_serialize)
 *   Phase 3: Write all records to TSV file
 *
 * Output: timing breakdown + summary statistics.
 */

#include "human_igh_imgt.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

static double now_sec(void) {
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return ts.tv_sec + ts.tv_nsec * 1e-9;
}

int main(int argc, char **argv) {
    int N = 10000;
    if (argc > 1) N = atoi(argv[1]);
    if (N <= 0) N = 10000;

    printf("GenAIRR-C End-to-End Benchmark\n");
    printf("==============================\n");
    printf("Generating %d sequences (HUMAN_IGH_IMGT + S5F)\n\n", N);

    /* ── Phase 1: Setup ──────────────────────────────────────── */
    double t0 = now_sec();

    SimConfig cfg;
    human_igh_imgt_load_config(&cfg);

    S5FModel s5f;
    s5f_model_init(&s5f, 0.05, 0.15, false);
    human_igh_imgt_load_s5f(&s5f);

    Pipeline pl = pipeline_build(&cfg);

    AlleleCorrectionSet ctx = allele_correction_set_build(&cfg);

    double t_setup = now_sec() - t0;
    printf("Phase 1 — Setup:            %8.3f ms\n", t_setup * 1000.0);
    printf("  V alleles: %d, D alleles: %d, J alleles: %d\n",
           cfg.v_alleles.count, cfg.d_alleles.count, cfg.j_alleles.count);
    printf("  Short-D threshold: %d\n\n", ctx.d_bitmap.short_d_threshold);

    /* ── Phase 2: Generate sequences ─────────────────────────── */
    AirrRecord *records = calloc(N, sizeof(AirrRecord));
    if (!records) {
        fprintf(stderr, "Failed to allocate %d records\n", N);
        return 1;
    }

    RngState bench_rng;
    rng_seed(&bench_rng, 42, 0);
    cfg.rng = &bench_rng;

    double t1 = now_sec();

    int success = 0;
    int total_len = 0;
    int total_mutations = 0;
    int total_productive = 0;
    double total_mutation_rate = 0.0;

    for (int i = 0; i < N; i++) {
        ASeq seq;
        aseq_init(&seq);
        SimRecord rec;
        sim_record_init(&rec);

        /* Rearrange */
        pipeline_execute(&pl, &cfg, &seq, &rec, 0, NULL, NULL);
        if (seq.length == 0) {
            aseq_reset(&seq);
            continue;
        }

        /* Mutate */
        S5FResult mres;
        s5f_mutate(&s5f, &seq, &rec, &bench_rng, &mres);

        /* Serialize */
        airr_serialize(&seq, &rec, &cfg, &ctx, &records[success]);

        total_len += records[success].sequence_length;
        total_mutations += records[success].n_mutations;
        total_mutation_rate += records[success].mutation_rate;
        if (records[success].productive) total_productive++;
        success++;

        aseq_reset(&seq);
    }

    double t_generate = now_sec() - t1;

    printf("Phase 2 — Generate:         %8.3f ms  (%d sequences)\n",
           t_generate * 1000.0, success);
    printf("  Per sequence:             %8.3f µs\n",
           t_generate / success * 1e6);
    printf("  Throughput:               %8.0f seq/sec\n",
           success / t_generate);
    printf("  Avg length:               %8.1f bp\n",
           (double)total_len / success);
    printf("  Avg mutations:            %8.1f per seq\n",
           (double)total_mutations / success);
    printf("  Avg mutation rate:        %8.4f\n",
           total_mutation_rate / success);
    printf("  Productive:               %8d / %d (%.1f%%)\n\n",
           total_productive, success,
           100.0 * total_productive / success);

    /* ── Phase 3: Write TSV ──────────────────────────────────── */
    double t2 = now_sec();

    FILE *fp = fopen("/tmp/genairr_bench.tsv", "w");
    if (fp) {
        airr_write_tsv_header(fp);
        for (int i = 0; i < success; i++) {
            airr_write_tsv_row(fp, &records[i]);
        }
        fclose(fp);
    }

    double t_write = now_sec() - t2;
    printf("Phase 3 — Write TSV:        %8.3f ms\n", t_write * 1000.0);

    /* ── Summary ─────────────────────────────────────────────── */
    double t_total = t_setup + t_generate + t_write;
    printf("\n");
    printf("Total wall time:            %8.3f ms\n", t_total * 1000.0);
    printf("  Setup:                    %8.1f%%\n", 100.0 * t_setup / t_total);
    printf("  Generate:                 %8.1f%%\n", 100.0 * t_generate / t_total);
    printf("  Write:                    %8.1f%%\n", 100.0 * t_write / t_total);

    if (fp) {
        /* Get file size */
        fp = fopen("/tmp/genairr_bench.tsv", "r");
        if (fp) {
            fseek(fp, 0, SEEK_END);
            long size = ftell(fp);
            fclose(fp);
            printf("\nOutput: /tmp/genairr_bench.tsv (%.1f MB, %d rows)\n",
                   size / 1e6, success);
        }
    }

    free(records);
    allele_correction_set_destroy(&ctx);
    sim_config_destroy(&cfg);
    return 0;
}
