#define _USE_MATH_DEFINES // for C
#include <math.h>

typedef void corrxor_fun_t(
    uint32_t *sig, uint32_t n_sig, uint32_t *ref, uint32_t n_ref,
    uint32_t *out);

corrxor_fun_t corrxor;
corrxor_fun_t corrxor_popcount;
corrxor_fun_t corrxor_nopop;
corrxor_fun_t corrxor_popcount_3quad;
