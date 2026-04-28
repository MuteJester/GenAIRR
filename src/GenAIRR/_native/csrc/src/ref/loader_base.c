/**
 * loader_base.c — flat C wrappers around the ReferenceLoader vtable.
 *
 * The vtable is great in C (allows function-pointer dispatch through
 * a struct field) but ctypes-from-Python can't follow function
 * pointers stored in arbitrary struct fields without writing a fragile
 * struct-layout mirror for every concrete loader. Instead we expose
 * `reference_loader_next` / `reference_loader_close` as flat C
 * functions that thinly forward to the vtable. Python only needs to
 * know about the opaque `ReferenceLoader *` handle.
 *
 * Both wrappers are NULL-safe and idempotent in spirit:
 *   - next() on a closed loader returns -1 (consistent with the
 *     "loader has been closed" error path).
 *   - close() on NULL is a no-op.
 */

#include <stddef.h>
#include "genairr/ref_loader.h"

int reference_loader_next(ReferenceLoader *self,
                          LoadedAlleleRecord *out,
                          const char **err_msg) {
    if (!self || !self->vt || !self->vt->next) {
        if (err_msg) *err_msg = "loader is NULL or has no next() method";
        return -1;
    }
    return self->vt->next(self, out, err_msg);
}

void reference_loader_close(ReferenceLoader *self) {
    if (!self) return;
    if (self->vt && self->vt->close) {
        self->vt->close(self);
    }
}
