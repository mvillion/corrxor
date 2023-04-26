//#include <immintrin.h>
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

#if defined(__ARM_FEATURE_DSP)
static inline uint32_t __SUB_LSR1(uint32_t op1, uint32_t op2)
{
    uint32_t result;
    __asm("sub %0, %1, %2, lsr #1" : "=r" (result) : "r" (op1), "r" (op2));
    return result;
}

static inline uint32_t __ADD_LSR(uint32_t op1, uint32_t op2, uint8_t op3)
{
    uint32_t result;
    __asm("add %0, %1, %2, lsr %3" : "=r" (result) : "r" (op1), "r" (op2), "I" (op3));
    return result;
}
#endif

static inline uint32_t popcount_quad(uint32_t x)
{
    uint32_t x1 = x & 0xaaaaaaaa;
#if defined(__ARM_FEATURE_DSP)
    x = __SUB_LSR1(x, x1);
#else
    x -= x1 >> 1;
#endif
    x1 = x & 0xcccccccc;
    x &= 0x33333333;
#if defined(__ARM_FEATURE_DSP)
    x = __ADD_LSR(x, x1, 2);
#else
    x += x1 >> 2;
#endif
    return x;
}

static inline uint32_t popcount_octet(uint32_t x)
{
    x = popcount_quad(x);
    x = (x + (x >> 4)) & 0x0f0f0f0f;
    return x;
}

static inline uint32_t sum_octet(uint32_t x)
{
#if defined(__ARM_FEATURE_DSP) && 0
    x = __USAD8(x, 0);
#else
    #if 0
        x += (x << 8);
        x += (x << 16);
    #else
        x *= 0x01010101;
    #endif
    x >>= 24;
#endif
    return x;
}

// sum may be over 256
static inline uint32_t sum_octet_over(uint32_t x)
{
#if defined(__ARM_FEATURE_DSP) && 0
    x = __USAD8(x, 0);
#else
    uint32_t x1 = x & 0xff00ff00;
    x &= 0x00ff00ff;
    #if defined(__ARM_FEATURE_DSP)
        x = __ADD_LSR(x, x1, 8);
    #else
        x += x1 >> 8;
    #endif
    x += (x << 16);
    x >>= 16;
#endif
    return x;
}

static inline uint32_t popcount(uint32_t x)
{
    x = popcount_octet(x);
    return sum_octet(x);
}

// method:
// bit 0 is __builtin_popcount or my popcount
// bit 1 (2) is popcount accumulated on 3 quad
// bit 2 (4) is popcount accumulated on 3 quad

static void __attribute__((always_inline)) inline corrxor_popcount_template(
    uint32_t *sig, uint32_t n_sig, uint32_t *ref, uint32_t n_ref,
    uint32_t *out, uint8_t method)
{
    uint32_t acc_k;
    uint32_t k;
    uint32_t n_out = n_sig-n_ref+1;
    uint32_t sig_k;

    for (uint32_t i_out = 0; i_out < n_out; i_out++)
    {
        uint32_t acc = 0;
        uint8_t shift_low = i_out % 32;
        uint8_t shift_high = 32-shift_low;
        if (shift_low == 0 && 0)
            for (k = 0; k < n_ref/32; k++)
            {
                uint32_t sig_k = sig[k+i_out/32];
                uint32_t x = ref[k]^sig_k;
                if ((method & 1) == 0)
                    acc_k = __builtin_popcount(x);
                else // (method == 1)
                    acc_k = popcount(x);
                acc += acc_k;
            }
        else
        {
            k = 0;
            if ((method & 0x2) > 0)
                for (k = 0; k < n_ref/32/3*3; k += 3)
                {
                    uint32_t sig_low = sig[k+0+i_out/32];
                    uint32_t sig_high = sig[k+1+i_out/32];
                    sig_k = sig_low >> shift_low;
                    sig_k |= (uint32_t)((uint64_t)sig_high << shift_high);
                    acc_k = popcount_quad(ref[k+0]^sig_k);

                    sig_low = sig[k+1+i_out/32];
                    sig_high = sig[k+2+i_out/32];
                    sig_k = sig_low >> shift_low;
                    sig_k |= (uint32_t)((uint64_t)sig_high << shift_high);
                    acc_k += popcount_quad(ref[k+1]^sig_k);

                    sig_low = sig[k+2+i_out/32];
                    sig_high = sig[k+3+i_out/32];
                    sig_k = sig_low >> shift_low;
                    sig_k |= (uint32_t)((uint64_t)sig_high << shift_high);
                    acc_k += popcount_quad(ref[k+2]^sig_k);

                    uint32_t x1 = acc_k & 0xf0f0f0f0;
                    acc_k &= 0x0f0f0f0f;
#if defined(__ARM_FEATURE_DSP)
                    acc_k = __ADD_LSR(acc_k, x1, 4);
#else
                    acc_k += x1 >> 4;
#endif
                    acc += sum_octet(acc_k);
                }
            else if ((method & 0x4) > 0)
            {
#define UNROLL_POP8 8
                uint32_t n_ref8 =  n_ref/32/UNROLL_POP8*UNROLL_POP8;
                for (k = 0; k < n_ref8; k += UNROLL_POP8)
                {
                    acc_k = 0;
                    #pragma GCC unroll 8
                    for (uint32_t l = 0; l < UNROLL_POP8; l++)
                    {
                        uint32_t sig_low = sig[k+l+i_out/32];
                        uint32_t sig_high = sig[k+l+1+i_out/32];
                        sig_k = sig_low >> shift_low;
                        sig_k |= (uint32_t)((uint64_t)sig_high << shift_high);
                        acc_k += popcount_octet(ref[k+l]^sig_k);
                    }
                    acc += sum_octet_over(acc_k);
                }
            }
            for (; k < n_ref/32; k++)
            {
                uint32_t sig_low = sig[k+i_out/32];
                uint32_t sig_high = sig[k+i_out/32+1];
                sig_k = sig_low >> shift_low;
                sig_k |= (uint32_t)((uint64_t)sig_high << shift_high);
                uint32_t x = ref[k]^sig_k;
                if ((method & 1) == 0)
                    acc_k = __builtin_popcount(x);
                else // (method == 1)
                    acc_k = popcount(x);
                acc += acc_k;
            }
        }
        for (k = n_ref/32*32; k < n_ref; k++)
        {
            uint32_t ref_k = ref[k/32];
            ref_k >>= k % 32;
            ref_k &= 1;
            sig_k = sig[(k+i_out)/32];
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
    corrxor_popcount_template(sig, n_sig, ref, n_ref, out, 0);
}

void corrxor_nopop(
    uint32_t *sig, uint32_t n_sig, uint32_t *ref, uint32_t n_ref,
    uint32_t *out)
{
    corrxor_popcount_template(sig, n_sig, ref, n_ref, out, 1);
}

void corrxor_popcount_3quad(
    uint32_t *sig, uint32_t n_sig, uint32_t *ref, uint32_t n_ref,
    uint32_t *out)
{
    corrxor_popcount_template(sig, n_sig, ref, n_ref, out, 3);
}

void corrxor_popcount_8octet(
    uint32_t *sig, uint32_t n_sig, uint32_t *ref, uint32_t n_ref,
    uint32_t *out)
{
    corrxor_popcount_template(sig, n_sig, ref, n_ref, out, 5);
}
