#include <immintrin.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/param.h>
#include "dsp.h"

#define concat2(X, Y) X ## Y
#define concat(X, Y) concat2(X, Y)

#define FLT_TYPE double
#define GOERTZEL_PREFIX(a) goertzel ## a
#include "dsp_simple.c"
#undef GOERTZEL_PREFIX
#undef FLT_TYPE

#define FLT_TYPE float
#define GOERTZEL_PREFIX(a) goertzelf ## a
#include "dsp_simple.c"
#undef GOERTZEL_PREFIX
#undef FLT_TYPE
