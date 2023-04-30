#include "include.h"
#include <math.h>

static PyObject* dsp_corrxor_template(
    PyObject* self, PyObject* args, corrxor_fun_t *fun)
{
    PyArrayObject *in_sig;
    uint32_t n_sig;
    PyArrayObject *in_ref;
    uint32_t n_ref;
    PyObject *output;

    int ret = PyArg_ParseTuple(
        args, "O!IO!I", &PyArray_Type, &in_sig, &n_sig, &PyArray_Type,
        &in_ref, &n_ref);
    if (!ret || in_sig == NULL || in_ref == NULL)
        return NULL;

    // check indicated length is compatible with the last dimemssion
    uint32_t n_sig32 = PyArray_DIM(in_sig, PyArray_NDIM(in_sig)-1);
    if (n_sig32*32 < n_sig)
        return NULL;
    uint32_t n_ref32 = PyArray_DIM(in_ref, PyArray_NDIM(in_ref)-1);
    if (n_ref32*32 < n_ref)
        return NULL;

    int n_dim = PyArray_NDIM(in_sig);
    if (n_dim != PyArray_NDIM(in_ref))
        return NULL;

    // create output dimensions
    // last axis is removed, replaced by complex data i&q
    uint32_t n_out = n_sig-n_ref+1;
    npy_intp out_dim[NPY_MAXDIMS];
    memcpy(out_dim, PyArray_DIMS(in_sig), n_dim*sizeof(npy_intp));
    out_dim[n_dim-1] = n_out;

    npy_intp n_data = 1;
    for (int k = 0; k < n_dim-1; k++)
    {
        n_data *= out_dim[k];
        // @todo: dimension broadcasting is not implemented
        if (out_dim[k] != PyArray_DIM(in_ref, k))
            return NULL;
    }

    // Ensure the input array is contiguous.
    // PyArray_GETCONTIGUOUS will increase the reference count.
    in_sig = PyArray_GETCONTIGUOUS(in_sig);
    in_ref = PyArray_GETCONTIGUOUS(in_ref);

    uint32_t *sig = (uint32_t *)PyArray_DATA(in_sig);
    uint32_t *ref = (uint32_t *)PyArray_DATA(in_ref);
    output = PyArray_SimpleNew(n_dim, out_dim, NPY_UINT32);
    uint32_t *out_res = (uint32_t *)PyArray_DATA((PyArrayObject *)output);
    for (npy_intp i_data = 0; i_data < n_data; i_data++)
    {
        fun(sig, n_sig, ref, n_ref, out_res);
        sig += n_sig32;
        ref += n_ref32;
        out_res += n_out;
    }

    // Decrease the reference count
    Py_DECREF(in_ref);
    Py_DECREF(in_sig);
    return output;
}

#define DEF_DSP(name) \
static PyObject* dsp_ ## name (PyObject* self, PyObject* args) \
{ \
    return dsp_corrxor_template(self, args, &name); \
}
DEF_DSP(corrxor)
DEF_DSP(corrxor_popcount)
DEF_DSP(corrxor_nopop)
DEF_DSP(corrxor_popcount_3quad)
DEF_DSP(corrxor_popcount_8octet)
#undef DEF_DSP


#define stringify(x) #x
#define DEF_DSP(radix, arch) \
    { \
        stringify(corrxor_rad ## radix ## _ ## arch), \
        dsp_corrxor_rad ## radix ## _ ## arch, METH_VARARGS, \
        "corrxor radix-" stringify(radix) \
        " algorithm using " stringify(arch) " instructions." \
    },

/* Set up the methods table */
static PyMethodDef methods[] = {
    {
        "corrxor", dsp_corrxor, // Python name, C name
        METH_VARARGS, // input parameters
        "corrxor algorithm." // doc string
    },
    {
        "corrxor_popcount", dsp_corrxor_popcount, METH_VARARGS,
        "corrxor algorithm using popcount instruction."
    },
    {
        "corrxor_nopop", dsp_corrxor_nopop, METH_VARARGS,
        "corrxor algorithm w/o popcount instruction."
    },
    {
        "corrxor_popcount_3quad", dsp_corrxor_popcount_3quad, METH_VARARGS,
        "corrxor algorithm using popcount w/ 3 quad accumulation."
    },
    {
        "corrxor_popcount_8octet", dsp_corrxor_popcount_8octet, METH_VARARGS,
        "corrxor algorithm using popcount w/ 8 octet accumulation."
    },
    {
        "corrxor3_nopop", dsp_corrxor_nopop, METH_VARARGS,
        "corrxor algorithm w/o popcount instruction (n_delay=3)."
    },
    {NULL, NULL, 0, NULL} // sentinel
};

/* Initialize module */
#if PY_VERSION_HEX >= 0x03000000
static struct PyModuleDef moduledef = {
    PyModuleDef_HEAD_INIT,
    "corrxor",
    NULL,
    -1,
    methods,
    NULL,
    NULL,
    NULL,
    NULL
};
PyMODINIT_FUNC PyInit_corrxor(void)
{
    import_array(); // Must be called for NumPy.
    PyObject *m;
    m = PyModule_Create(&moduledef);
    if (!m) return NULL;
    return m;
}
#else
PyMODINIT_FUNC initcorrxor(void)
{
    (void)Py_InitModule("corrxor", methods);
    import_array(); // Must be called for NumPy.
}
#endif
