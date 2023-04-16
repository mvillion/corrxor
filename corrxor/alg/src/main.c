#include "include.h"
#include <math.h>

static PyObject* dsp_corrxor_template(
    PyObject* self, PyObject* args, corrxor_fun_t *fun)
{
    PyArrayObject *in_data;
    PyObject *output;
    double k;

    if (!PyArg_ParseTuple(args, "O!d", &PyArray_Type, &in_data, &k))
        return NULL;
    if (in_data == NULL) return NULL;

    // Ensure the input array is contiguous.
    // PyArray_GETCONTIGUOUS will increase the reference count.
    in_data = PyArray_GETCONTIGUOUS(in_data);

    // create output dimensions
    // last axis is removed, replaced by complex data i&q
    npy_intp out_dim[NPY_MAXDIMS];
    int n_dim = PyArray_NDIM(in_data);
    memcpy(out_dim, PyArray_DIMS(in_data), n_dim*sizeof(npy_intp));
    long int data_len = out_dim[n_dim-1];
    npy_intp n_data = 1;
    for (int k = 0; k < n_dim-1; k++)
        n_data *= out_dim[k];

    double *data = (double *)PyArray_DATA(in_data);
    output = PyArray_SimpleNew(n_dim-1, out_dim, NPY_COMPLEX128);
    double *out_res = (double *)PyArray_DATA((PyArrayObject *)output);
    for (npy_intp i_data = 0; i_data < n_data; i_data++)
    {
        fun(data, data_len, k, out_res);
        data += data_len;
        out_res += 2;
    }

    // Decrease the reference count of ap.
    Py_DECREF(in_data);
    return output;
}

#define DEF_DSP(name) \
static PyObject* dsp_ ## name (PyObject* self, PyObject* args) \
{ \
    return dsp_corrxor_template(self, args, &name); \
}
DEF_DSP(corrxor)
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
