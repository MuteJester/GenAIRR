/**
 * export.h — single source of truth for the GENAIRR_EXPORT macro.
 *
 * Declares which symbols leave the shared library. On Windows we need
 * __declspec(dllexport) explicitly; on ELF / Mach-O the default
 * visibility is "default" so the macro is a no-op.
 *
 * Both genairr_api.h and the anchor/ref subsystem headers include
 * this so the same macro works for every public function.
 */

#ifndef GENAIRR_EXPORT_H
#define GENAIRR_EXPORT_H

#ifdef _WIN32
  #define GENAIRR_EXPORT __declspec(dllexport)
#else
  #define GENAIRR_EXPORT
#endif

#endif /* GENAIRR_EXPORT_H */
