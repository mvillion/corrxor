#define _USE_MATH_DEFINES // for C
#include <math.h>

typedef void corrxor_fun_t(
    double *data, long data_len, double k, double *out);

corrxor_fun_t corrxor;
