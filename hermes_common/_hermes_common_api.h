#ifndef __PYX_HAVE_API___hermes_common
#define __PYX_HAVE_API___hermes_common
#include "Python.h"

static PyObject *(*c2numpy_int)(int *, int);
static PyObject *(*c2numpy_double)(double *, int);
static void (*numpy2c_double_inplace)(PyObject *, double **, int *);
static void (*numpy2c_int_inplace)(PyObject *, int **, int *);
static PyObject *(*c2py_AVector)(struct AVector *);
static PyObject *(*c2py_DenseMatrix)(struct DenseMatrix *);
static PyObject *(*c2py_CooMatrix)(struct CooMatrix *);
static PyObject *(*c2py_CSRMatrix)(struct CSRMatrix *);
static PyObject *(*c2py_CSCMatrix)(struct CSCMatrix *);
static PyObject *(*namespace_create)(void);
static void (*namespace_push)(PyObject *, const char*, PyObject *);
static void (*namespace_print)(PyObject *);
static PyObject *(*namespace_pull)(PyObject *, const char*);
static void (*cmd)(const char*);
static void (*set_verbose_cmd)(int);
static void (*insert_object)(const char*, PyObject *);
static PyObject *(*get_object)(const char*);
static PyObject *(*c2py_int)(int);
static int (*py2c_int)(PyObject *);
static char *(*py2c_str)(PyObject *);
static double (*py2c_double)(PyObject *);
static PyObject *(*c2numpy_int_inplace)(int *, int);
static PyObject *(*c2numpy_double_inplace)(double *, int);
static PyObject *(*c2numpy_double_complex_inplace)(__pyx_t_14_hermes_common_cplx *, int);
static void (*numpy2c_double_complex_inplace)(PyObject *, __pyx_t_14_hermes_common_cplx **, int *);
static void (*run_cmd)(const char*, PyObject *);

#ifndef __PYX_HAVE_API_FUNC_import_module
#define __PYX_HAVE_API_FUNC_import_module

#ifndef __PYX_HAVE_RT_ImportModule
#define __PYX_HAVE_RT_ImportModule
static PyObject *__Pyx_ImportModule(const char *name) {
    PyObject *py_name = 0;
    PyObject *py_module = 0;

    #if PY_MAJOR_VERSION < 3
    py_name = PyString_FromString(name);
    #else
    py_name = PyUnicode_FromString(name);
    #endif
    if (!py_name)
        goto bad;
    py_module = PyImport_Import(py_name);
    Py_DECREF(py_name);
    return py_module;
bad:
    Py_XDECREF(py_name);
    return 0;
}
#endif

#endif


#ifndef __PYX_HAVE_RT_ImportFunction
#define __PYX_HAVE_RT_ImportFunction
static int __Pyx_ImportFunction(PyObject *module, const char *funcname, void (**f)(void), const char *sig) {
#if PY_VERSION_HEX < 0x02050000
    char *api = (char *)"__pyx_capi__";
#else
    const char *api = "__pyx_capi__";
#endif
    PyObject *d = 0;
    PyObject *cobj = 0;
    const char *desc;
    const char *s1, *s2;
    union {
        void (*fp)(void);
        void *p;
    } tmp;

    d = PyObject_GetAttrString(module, api);
    if (!d)
        goto bad;
    cobj = PyDict_GetItemString(d, funcname);
    if (!cobj) {
        PyErr_Format(PyExc_ImportError,
            "%s does not export expected C function %s",
                PyModule_GetName(module), funcname);
        goto bad;
    }
    desc = (const char *)PyCObject_GetDesc(cobj);
    if (!desc)
        goto bad;
    s1 = desc; s2 = sig;
    while (*s1 != '\0' && *s1 == *s2) { s1++; s2++; }
    if (*s1 != *s2) {
        PyErr_Format(PyExc_TypeError,
            "C function %s.%s has wrong signature (expected %s, got %s)",
             PyModule_GetName(module), funcname, sig, desc);
        goto bad;
    }
    tmp.p = PyCObject_AsVoidPtr(cobj);
    *f = tmp.fp;
    Py_DECREF(d);
    return 0;
bad:
    Py_XDECREF(d);
    return -1;
}
#endif

static int import__hermes_common(void) {
  PyObject *module = 0;
  module = __Pyx_ImportModule("_hermes_common");
  if (!module) goto bad;
  if (__Pyx_ImportFunction(module, "c2numpy_int", (void (**)(void))&c2numpy_int, "PyObject *(int *, int)") < 0) goto bad;
  if (__Pyx_ImportFunction(module, "c2numpy_double", (void (**)(void))&c2numpy_double, "PyObject *(double *, int)") < 0) goto bad;
  if (__Pyx_ImportFunction(module, "numpy2c_double_inplace", (void (**)(void))&numpy2c_double_inplace, "void (PyObject *, double **, int *)") < 0) goto bad;
  if (__Pyx_ImportFunction(module, "numpy2c_int_inplace", (void (**)(void))&numpy2c_int_inplace, "void (PyObject *, int **, int *)") < 0) goto bad;
  if (__Pyx_ImportFunction(module, "c2py_AVector", (void (**)(void))&c2py_AVector, "PyObject *(struct AVector *)") < 0) goto bad;
  if (__Pyx_ImportFunction(module, "c2py_DenseMatrix", (void (**)(void))&c2py_DenseMatrix, "PyObject *(struct DenseMatrix *)") < 0) goto bad;
  if (__Pyx_ImportFunction(module, "c2py_CooMatrix", (void (**)(void))&c2py_CooMatrix, "PyObject *(struct CooMatrix *)") < 0) goto bad;
  if (__Pyx_ImportFunction(module, "c2py_CSRMatrix", (void (**)(void))&c2py_CSRMatrix, "PyObject *(struct CSRMatrix *)") < 0) goto bad;
  if (__Pyx_ImportFunction(module, "c2py_CSCMatrix", (void (**)(void))&c2py_CSCMatrix, "PyObject *(struct CSCMatrix *)") < 0) goto bad;
  if (__Pyx_ImportFunction(module, "namespace_create", (void (**)(void))&namespace_create, "PyObject *(void)") < 0) goto bad;
  if (__Pyx_ImportFunction(module, "namespace_push", (void (**)(void))&namespace_push, "void (PyObject *, const char*, PyObject *)") < 0) goto bad;
  if (__Pyx_ImportFunction(module, "namespace_print", (void (**)(void))&namespace_print, "void (PyObject *)") < 0) goto bad;
  if (__Pyx_ImportFunction(module, "namespace_pull", (void (**)(void))&namespace_pull, "PyObject *(PyObject *, const char*)") < 0) goto bad;
  if (__Pyx_ImportFunction(module, "cmd", (void (**)(void))&cmd, "void (const char*)") < 0) goto bad;
  if (__Pyx_ImportFunction(module, "set_verbose_cmd", (void (**)(void))&set_verbose_cmd, "void (int)") < 0) goto bad;
  if (__Pyx_ImportFunction(module, "insert_object", (void (**)(void))&insert_object, "void (const char*, PyObject *)") < 0) goto bad;
  if (__Pyx_ImportFunction(module, "get_object", (void (**)(void))&get_object, "PyObject *(const char*)") < 0) goto bad;
  if (__Pyx_ImportFunction(module, "c2py_int", (void (**)(void))&c2py_int, "PyObject *(int)") < 0) goto bad;
  if (__Pyx_ImportFunction(module, "py2c_int", (void (**)(void))&py2c_int, "int (PyObject *)") < 0) goto bad;
  if (__Pyx_ImportFunction(module, "py2c_str", (void (**)(void))&py2c_str, "char *(PyObject *)") < 0) goto bad;
  if (__Pyx_ImportFunction(module, "py2c_double", (void (**)(void))&py2c_double, "double (PyObject *)") < 0) goto bad;
  if (__Pyx_ImportFunction(module, "c2numpy_int_inplace", (void (**)(void))&c2numpy_int_inplace, "PyObject *(int *, int)") < 0) goto bad;
  if (__Pyx_ImportFunction(module, "c2numpy_double_inplace", (void (**)(void))&c2numpy_double_inplace, "PyObject *(double *, int)") < 0) goto bad;
  if (__Pyx_ImportFunction(module, "c2numpy_double_complex_inplace", (void (**)(void))&c2numpy_double_complex_inplace, "PyObject *(__pyx_t_14_hermes_common_cplx *, int)") < 0) goto bad;
  if (__Pyx_ImportFunction(module, "numpy2c_double_complex_inplace", (void (**)(void))&numpy2c_double_complex_inplace, "void (PyObject *, __pyx_t_14_hermes_common_cplx **, int *)") < 0) goto bad;
  if (__Pyx_ImportFunction(module, "run_cmd", (void (**)(void))&run_cmd, "void (const char*, PyObject *)") < 0) goto bad;
  Py_DECREF(module); module = 0;
  return 0;
  bad:
  Py_XDECREF(module);
  return -1;
}

#endif
