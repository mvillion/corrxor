#define _USE_MATH_DEFINES // for C
#include <math.h>

typedef void goertzel_fun_t(
    double *data, long data_len, double k, double *out);

typedef void goertzelf_fun_t(
    float *data, long data_len, float k, float *out);

// Goertzel algorithm (for single tone detection)
goertzel_fun_t goertzel;
goertzel_fun_t goertzel_cx;

goertzelf_fun_t goertzelf;
goertzelf_fun_t goertzelf_cx;
