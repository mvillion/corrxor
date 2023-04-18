#include <immintrin.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/param.h>
#include "dsp.h"

#define concat2(X, Y) X ## Y
#define concat(X, Y) concat2(X, Y)

void corrxor(
    uint32_t *sig, uint32_t n_sig, uint32_t *ref, uint32_t n_ref,
    uint32_t *out)
{
    uint32_t n_out = n_sig-n_ref+1;

    for (uint32_t i_out = 0; i_out < n_out; i_out++)
    {
        uint32_t acc = 0;
        for (uint32_t k = 0; k < n_ref; k++)
        {
            uint32_t ref_k = ref[k/32];
            ref_k >>= k % 32;
            ref_k &= 1;
            uint32_t sig_k = sig[(k+i_out)/32];
            sig_k >>= (k+i_out) % 32;
            sig_k &= 1;
            acc += ref_k^sig_k;
        }
        out[i_out] = n_ref-acc;
    }
}

void corrxor_popcount(
    uint32_t *sig, uint32_t n_sig, uint32_t *ref, uint32_t n_ref,
    uint32_t *out)
{
    uint32_t n_out = n_sig-n_ref+1;

    for (uint32_t i_out = 0; i_out < n_out; i_out++)
    {
        uint32_t acc = 0;
        uint8_t shift_low = i_out % 32;
        uint8_t shift_high = 32-shift_low;
        if (shift_low == 0)
            for (uint32_t k = 0; k < n_ref/32; k++)
            {
                uint32_t ref_k = ref[k];
                uint32_t sig_k = sig[k+i_out/32];
                acc += __builtin_popcount(ref_k^sig_k);
            }
        else
            for (uint32_t k = 0; k < n_ref/32; k++)
            {
                uint32_t ref_k = ref[k];
                uint32_t sig_low = sig[k+i_out/32];
                uint32_t sig_high = sig[k+i_out/32+1];
                uint32_t sig_k;
                sig_k = sig_low >> shift_low;
                sig_k |= sig_high << shift_high;
                acc += __builtin_popcount(ref_k^sig_k);
            }
        for (uint32_t k = n_ref/32*32; k < n_ref; k++)
        {
            uint32_t ref_k = ref[k/32];
            ref_k >>= k % 32;
            ref_k &= 1;
            uint32_t sig_k = sig[(k+i_out)/32];
            sig_k >>= (k+i_out) % 32;
            sig_k &= 1;
            acc += ref_k^sig_k;
        }
        out[i_out] = n_ref-acc;
    }
}

void corrxor_nopop(
    uint32_t *sig, uint32_t n_sig, uint32_t *ref, uint32_t n_ref,
    uint32_t *out)
{
    uint32_t n_out = n_sig-n_ref+1;

    for (uint32_t i_out = 0; i_out < n_out; i_out++)
    {
        uint32_t acc = 0;
        uint8_t shift_low = i_out % 32;
        uint8_t shift_high = 32-shift_low;
        if (shift_low == 0)
            for (uint32_t k = 0; k < n_ref/32; k++)
            {
                uint32_t ref_k = ref[k];
                uint32_t sig_k = sig[k+i_out/32];
                uint32_t x = ref_k^sig_k;
                x -= ((x >> 1) & 0x55555555);
                x = (x & 0x33333333) + ((x >> 2) & 0x33333333);
                x = (x + (x >> 4)) & 0x0f0f0f0f;
                x += (x >> 8);
                x += (x >> 16);
                acc += x & 0x3f;
            }
        else
            for (uint32_t k = 0; k < n_ref/32; k++)
            {
                uint32_t ref_k = ref[k];
                uint32_t sig_low = sig[k+i_out/32];
                uint32_t sig_high = sig[k+i_out/32+1];
                uint32_t sig_k;
                sig_k = sig_low >> shift_low;
                sig_k |= sig_high << shift_high;
                uint32_t x = ref_k^sig_k;
                x -= ((x >> 1) & 0x55555555);
                x = (x & 0x33333333) + ((x >> 2) & 0x33333333);
                x = (x + (x >> 4)) & 0x0f0f0f0f;
                x += (x >> 8);
                x += (x >> 16);
                acc += x & 0x3f;
            }
        for (uint32_t k = n_ref/32*32; k < n_ref; k++)
        {
            uint32_t ref_k = ref[k/32];
            ref_k >>= k % 32;
            ref_k &= 1;
            uint32_t sig_k = sig[(k+i_out)/32];
            sig_k >>= (k+i_out) % 32;
            sig_k &= 1;
            acc += ref_k^sig_k;
        }
        out[i_out] = n_ref-acc;
    }
}
