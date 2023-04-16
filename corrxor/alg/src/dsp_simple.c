void GOERTZEL_PREFIX()(FLT_TYPE *data, long data_len, FLT_TYPE k, FLT_TYPE *out)
{
    FLT_TYPE omega = 2.0*M_PI*k/data_len;
    FLT_TYPE sw = sin(omega);
    FLT_TYPE cw = cos(omega);
    FLT_TYPE coeff = 2.0*cw;

    FLT_TYPE q0 = 0.0;
    FLT_TYPE q1 = 0.0;
    FLT_TYPE q2 = 0.0;

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

void GOERTZEL_PREFIX(_cx)(
    FLT_TYPE *data, long data_len, FLT_TYPE k, FLT_TYPE *out)
{
    FLT_TYPE omega = 2.0*M_PI*k/data_len;
    FLT_TYPE sw = sin(omega);
    FLT_TYPE cw = cos(omega);
    FLT_TYPE coeff = 2.0*cw;

    FLT_TYPE q0a = 0.0; // for real part
    FLT_TYPE q1a = 0.0;
    FLT_TYPE q2a = 0.0;
    FLT_TYPE q0b = 0.0; // for imaginary part
    FLT_TYPE q1b = 0.0;
    FLT_TYPE q2b = 0.0;

    long int i;
    for (i = 0; i < data_len/3*3; i += 3)
    {
        q0a = coeff*q1a - q2a + data[2*i+0];
        q0b = coeff*q1b - q2b + data[2*i+1];
        q2a = coeff*q0a - q1a + data[2*i+2];
        q2b = coeff*q0b - q1b + data[2*i+3];
        q1a = coeff*q2a - q0a + data[2*i+4];
        q1b = coeff*q2b - q0b + data[2*i+5];
    }
    for (; i < data_len; i++)
    {
        q0a = coeff*q1a - q2a + data[2*i+0];
        q0b = coeff*q1b - q2b + data[2*i+1];
        q2a = q1a;
        q2b = q1b;
        q1a = q0a;
        q1b = q0b;
    }

    out[0] = q1a*cw-q2a; // real
    out[1] = q1a*sw; // imag
    out[0] -= q1b*sw; // real
    out[1] += q1b*cw-q2b; // imag
}
