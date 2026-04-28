cmake_minimum_required(VERSION 3.14)
# Keep this VERSION in sync with `GENAIRR_C_VERSION_STRING` in
# include/genairr/genairr.h and the Python package version in
# setup.py. _load_lib() in src/GenAIRR/_native/__init__.py asserts
# the C library version matches Python at import time (T1-12).
project(genairr-c VERSION 1.0.0 LANGUAGES C)

set(CMAKE_C_STANDARD 11)
set(CMAKE_C_STANDARD_REQUIRED ON)

# ── Build options ────────────────────────────────────────────────
option(GENAIRR_BUILD_TESTS "Build C test executables" ON)
option(GENAIRR_BUILD_BENCHMARKS "Build C benchmark executables" ON)
option(GENAIRR_SANITIZE
       "Build with AddressSanitizer + UndefinedBehaviorSanitizer (GCC/Clang)"
       OFF)

# ── Compiler warnings (GCC/Clang only) ──────────────────────────
if(NOT MSVC)
    add_compile_options(-Wall -Wextra -Wpedantic -Wshadow -Wstrict-prototypes)
else()
    add_compile_options(/W3)
    # MSVC needs _USE_MATH_DEFINES for M_PI in <math.h>
    add_compile_definitions(_USE_MATH_DEFINES)
    # Suppress warnings about POSIX function names (fdopen -> _fdopen etc.)
    add_compile_definitions(_CRT_SECURE_NO_WARNINGS)
endif()

# ── Sanitizers (ASan + UBSan) ────────────────────────────────────
# Set -DGENAIRR_SANITIZE=ON at configure time to build with
# AddressSanitizer + UndefinedBehaviorSanitizer. Catches use-after-free,
# heap/stack overflows, leaks, signed-overflow, alignment, null derefs.
# Adds ~5-10x slowdown — intended for CI and local bug-finding, not
# production. Only supported on GCC/Clang (sanitizer surface on MSVC
# differs and isn't covered here).
if(GENAIRR_SANITIZE AND NOT MSVC)
    message(STATUS "GenAIRR: building with -fsanitize=address,undefined")
    add_compile_options(
        -fsanitize=address,undefined
        -fno-omit-frame-pointer
        -fno-sanitize-recover=all
        -g
    )
    add_link_options(-fsanitize=address,undefined)
endif()

# ── Include paths ─────────────────────────────────────────────────
include_directories(${CMAKE_SOURCE_DIR}/include)
include_directories(${CMAKE_SOURCE_DIR}/vendor/cjson)

# ── Library: libgenairr ───────────────────────────────────────────
set(LIB_SOURCES
    src/core/aseq.c
    src/allele/allele_pool.c
    src/pipeline/sim_config.c
    src/pipeline/pipeline.c
    src/rearrangement/sample.c
    src/rearrangement/trim.c
    src/rearrangement/assemble.c
    src/core/airr.c
    src/mutation/s5f.c
    src/mutation/uniform.c
    src/util/rand_util.c
    src/util/trace.c
    src/corruption/corrupt_5_prime.c
    src/corruption/corrupt_3_prime.c
    src/corruption/quality_errors.c
    src/corruption/pcr_amplification.c
    src/corruption/simulate_umi.c
    src/corruption/reverse_complement.c
    src/corruption/spike_contaminants.c
    src/corruption/insert_indels.c
    src/corruption/insert_ns.c
    src/corruption/paired_end.c
    src/corruption/long_read_errors.c
    src/corruption/trim_to_length.c
    src/corruption/primer_mask.c
    src/simulation/d_inversion.c
    src/simulation/receptor_revision.c
    src/simulation/selection_pressure.c
    src/simulation/simulate_csr.c
    src/core/allele_bitmap.c
    src/core/codon.c
    src/core/functionality.c
    src/io/gdc_write.c
    src/io/gdc_read.c
    src/api/genairr_api.c
    # Anchor resolution subsystem (T2-8)
    src/anchor/codon.c
    src/anchor/result.c
    src/anchor/strategy_explicit.c
    src/anchor/strategy_imgt_gapped.c
    src/anchor/strategy_motif_v.c
    src/anchor/strategy_motif_j.c
    src/anchor/resolver.c
    # Reference loaders (T2-8)
    src/ref/locus_infer.c
    src/ref/build_report.c
    src/ref/loader_base.c
    src/ref/loader_imgt.c
    src/ref/loader_plain.c
    src/ref/loader_airrc.c
    src/ref/loader_ogrdb.c
    src/ref/loader_igblast.c
    src/ref/loader_auto.c
    # Vendored cJSON for AIRR-C / OGRDB JSON loaders
    vendor/cjson/cJSON.c
)

add_library(genairr STATIC ${LIB_SOURCES})
target_include_directories(genairr PUBLIC ${CMAKE_SOURCE_DIR}/include)

# ── Shared library: libgenairr.so / .dylib / .dll (Python bindings) ─
# Note: genairr_api.c is already in LIB_SOURCES (so test executables
# can call the API). The shared lib reuses the same source list — no
# duplicate symbols since each library is built independently.

add_library(genairr_shared SHARED ${LIB_SOURCES})
set_target_properties(genairr_shared PROPERTIES
    OUTPUT_NAME genairr
    POSITION_INDEPENDENT_CODE ON
    SOVERSION 0
    VERSION ${PROJECT_VERSION}
)
target_include_directories(genairr_shared PUBLIC ${CMAKE_SOURCE_DIR}/include)

# Link math library (implicit on MSVC)
if(NOT MSVC)
    target_link_libraries(genairr_shared m)
endif()

# ── Tests ─────────────────────────────────────────────────────────
if(GENAIRR_BUILD_TESTS)
    enable_testing()

    # Helper: link math only on non-MSVC
    if(NOT MSVC)
        set(MATH_LIB m)
    else()
        set(MATH_LIB "")
    endif()

    add_executable(test_aseq tests/test_aseq.c)
    target_link_libraries(test_aseq genairr ${MATH_LIB})
    add_test(NAME test_aseq COMMAND test_aseq)

    add_executable(test_pipeline tests/test_pipeline.c)
    target_link_libraries(test_pipeline genairr ${MATH_LIB})
    add_test(NAME test_pipeline COMMAND test_pipeline)

    add_executable(test_airr tests/test_airr.c)
    target_link_libraries(test_airr genairr ${MATH_LIB})
    add_test(NAME test_airr COMMAND test_airr)

    add_executable(test_corrections tests/test_corrections.c)
    target_link_libraries(test_corrections genairr ${MATH_LIB})
    add_test(NAME test_corrections COMMAND test_corrections)

    add_executable(test_s5f tests/test_s5f.c)
    target_link_libraries(test_s5f genairr ${MATH_LIB})
    add_test(NAME test_s5f COMMAND test_s5f)

    add_executable(test_real_data tests/test_real_data.c data/human_igh_imgt.c)
    target_include_directories(test_real_data PRIVATE ${CMAKE_SOURCE_DIR}/data)
    target_link_libraries(test_real_data genairr ${MATH_LIB})
    add_test(NAME test_real_data COMMAND test_real_data)

    add_executable(test_serialize tests/test_serialize.c)
    target_link_libraries(test_serialize genairr ${MATH_LIB})
    add_test(NAME test_serialize COMMAND test_serialize)

    add_executable(test_ops tests/test_ops.c)
    target_link_libraries(test_ops genairr ${MATH_LIB})
    add_test(NAME test_ops COMMAND test_ops)

    add_executable(test_bitmap tests/test_bitmap.c)
    target_link_libraries(test_bitmap genairr ${MATH_LIB})
    add_test(NAME test_bitmap COMMAND test_bitmap)

    add_executable(test_functionality tests/test_functionality.c)
    target_link_libraries(test_functionality genairr ${MATH_LIB})
    add_test(NAME test_functionality COMMAND test_functionality)

    add_executable(test_codon_rail tests/test_codon_rail.c)
    target_link_libraries(test_codon_rail genairr ${MATH_LIB})
    add_test(NAME test_codon_rail COMMAND test_codon_rail)

    add_executable(test_gdc tests/test_gdc.c)
    target_link_libraries(test_gdc genairr ${MATH_LIB})
    add_test(NAME test_gdc COMMAND test_gdc)

    add_executable(test_rand_util tests/test_rand_util.c)
    target_link_libraries(test_rand_util genairr ${MATH_LIB})
    add_test(NAME test_rand_util COMMAND test_rand_util)

    # Anchor resolution subsystem tests (T2-8)
    add_executable(test_anchor tests/test_anchor.c)
    target_link_libraries(test_anchor genairr ${MATH_LIB})
    add_test(NAME test_anchor COMMAND test_anchor)

    add_executable(test_ref_loader tests/test_ref_loader.c)
    target_link_libraries(test_ref_loader genairr ${MATH_LIB})
    add_test(NAME test_ref_loader COMMAND test_ref_loader)

    add_executable(test_ref_loaders_phase2 tests/test_ref_loaders_phase2.c)
    target_link_libraries(test_ref_loaders_phase2 genairr ${MATH_LIB})
    add_test(NAME test_ref_loaders_phase2 COMMAND test_ref_loaders_phase2)
endif()

# ── Benchmarks ─────────────────────────────────────────────────
if(GENAIRR_BUILD_BENCHMARKS)
    if(NOT MSVC)
        set(MATH_LIB m)
    else()
        set(MATH_LIB "")
    endif()

    add_executable(bench_e2e benchmarks/bench_e2e.c data/human_igh_imgt.c)
    target_include_directories(bench_e2e PRIVATE ${CMAKE_SOURCE_DIR}/data)
    target_link_libraries(bench_e2e genairr ${MATH_LIB})
endif()
