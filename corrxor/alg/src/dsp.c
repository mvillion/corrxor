#include <immintrin.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/param.h>
#include "dsp.h"

#define concat2(X, Y) X ## Y
#define concat(X, Y) concat2(X, Y)

void corrxor(double *data, long data_len, double k, double *out)
{
    double omega = 2.0*M_PI*k/data_len;
    double sw = sin(omega);
    double cw = cos(omega);
    double coeff = 2.0*cw;

    double q0 = 0.0;
    double q1 = 0.0;
    double q2 = 0.0;

    long int i;
    for (i = 0; i < data_len/3*3; i += 3)
    {
        q0 = coeff*q1 - q2 + data[i+0];
        q2 = coeff*q0 - q1 + data[i+1];
        q1 = coeff*q2 - q0 + data[i+2];
    }
    for (; i < data_len; i++)
    {
        q0 = coeff*q1 - q2 + data[i];
        q2 = q1;
        q1 = q0;
    }

    // note: dm00446805-the-goertzel-algorithm-to-compute-individual-terms-of-the-discrete-fourier-transform-dft-stmicroelectronics-1.pdf
    // suggests for non-integer k:
//     w2 = 2*pi*k;
//     cw2 = cos(w2);
//     sw2 = sin(w2);
//     I = It*cw2 + Q*sw2;
//     Q = -It*sw2 + Q*cw2;

    out[0] = q1*cw-q2; // real
    out[1] = q1*sw; // imag
}
