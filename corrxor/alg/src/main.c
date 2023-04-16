#include "include.h"
#include <math.h>

static PyObject* dsp_goertzel_template(
    PyObject* self, PyObject* args, goertzel_fun_t *fun, goertzel_fun_t *fun_cx,
    goertzelf_fun_t *funf, goertzelf_fun_t *funf_cx)
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

    int typenum = PyArray_TYPE(in_data);
    if (typenum == NPY_FLOAT64 || typenum == NPY_COMPLEX128)
    {
        double *data = (double *)PyArray_DATA(in_data);
        output = PyArray_SimpleNew(n_dim-1, out_dim, NPY_COMPLEX128);
        double *out_res = (double *)PyArray_DATA((PyArrayObject *)output);
        if (typenum == NPY_FLOAT64 && fun != NULL)
            for (npy_intp i_data = 0; i_data < n_data; i_data++)
            {
                fun(data, data_len, k, out_res);
                data += data_len;
                out_res += 2;
            }
        else if (typenum == NPY_COMPLEX128 && fun_cx != NULL)
            for (npy_intp i_data = 0; i_data < n_data; i_data++)
            {
                fun_cx(data, data_len, k, out_res);
                data += data_len*2;
                out_res += 2;
            }
    }
    else //if (typenum == NPY_FLOAT32 || typenum == NPY_COMPLEX64)
    {
        float *data = (float *)PyArray_DATA(in_data);
        output = PyArray_SimpleNew(n_dim-1, out_dim, NPY_COMPLEX64);
        float *out_res = (float *)PyArray_DATA((PyArrayObject *)output);
        if (typenum == NPY_FLOAT32 && fun != NULL)
            for (npy_intp i_data = 0; i_data < n_data; i_data++)
            {
                funf(data, data_len, k, out_res);
                data += data_len;
                out_res += 2;
            }
        else if (typenum == NPY_COMPLEX64 && fun_cx != NULL)
            for (npy_intp i_data = 0; i_data < n_data; i_data++)
            {
                funf_cx(data, data_len, k, out_res);
                data += data_len*2;
                out_res += 2;
            }
    }

    // Decrease the reference count of ap.
    Py_DECREF(in_data);
    return output;
}

#define DEF_DSP(name, fun_cx, funf, funf_cx) \
static PyObject* dsp_ ## name (PyObject* self, PyObject* args) \
{ \
    return dsp_goertzel_template(self, args, &name, fun_cx, funf, funf_cx); \
}
DEF_DSP(goertzel, &goertzel_cx, &goertzelf, &goertzelf_cx)
#undef DEF_DSP


#define stringify(x) #x
#define DEF_DSP(radix, arch) \
    { \
        stringify(goertzel_rad ## radix ## _ ## arch), \
        dsp_goertzel_rad ## radix ## _ ## arch, METH_VARARGS, \
        "Goertzel radix-" stringify(radix) \
        " algorithm using " stringify(arch) " instructions." \
    },

/* Set up the methods table */
static PyMethodDef methods[] = {
    {
        "goertzel", dsp_goertzel, // Python name, C name
        METH_VARARGS, // input parameters
        "Goertzel algorithm." // doc string
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
