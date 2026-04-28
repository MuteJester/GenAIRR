/**
 * strategies.h — internal declarations of the anchor strategy fns.
 *
 * The resolver in resolver.c iterates a strategy stack of these
 * functions. Each lives in its own .c file for testability and clean
 * extension. This header is private — never installed.
 */

#ifndef GENAIRR_ANCHOR_STRATEGIES_H
#define GENAIRR_ANCHOR_STRATEGIES_H

#include "genairr/anchor.h"
#include "genairr/ref_record.h"

/* All strategy functions share a single signature so the resolver
 * can iterate them through a function-pointer array. Each returns
 * an AnchorResult; CONF_REJECTED tells the resolver to fall through
 * to the next strategy. */

typedef AnchorResult (*AnchorStrategyFn)(const AnchorResolverConfig *cfg,
                                         const LoadedAlleleRecord *rec);

/* V-anchor strategies — order matters in resolver.c */
AnchorResult try_v_explicit(const AnchorResolverConfig *cfg,
                            const LoadedAlleleRecord *rec);
AnchorResult try_v_imgt_gapped(const AnchorResolverConfig *cfg,
                               const LoadedAlleleRecord *rec);
AnchorResult try_v_motif(const AnchorResolverConfig *cfg,
                         const LoadedAlleleRecord *rec);

/* J-anchor strategies — order matters in resolver.c */
AnchorResult try_j_explicit(const AnchorResolverConfig *cfg,
                            const LoadedAlleleRecord *rec);
AnchorResult try_j_motif(const AnchorResolverConfig *cfg,
                         const LoadedAlleleRecord *rec);

#endif /* GENAIRR_ANCHOR_STRATEGIES_H */
