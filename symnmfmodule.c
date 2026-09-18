#define PY_SSIZE_T_CLEAN
#include <Python.h>
#include "symnmf.h"

/* Converts a 2D Python list into a dynamically allocated C matrix */
static double **py_to_matrix(PyObject *obj, int *rows, int *cols) {
    double **matrix;
    PyObject *row;
    Py_ssize_t i, j;

    *rows = (int)PyList_Size(obj);
    row = PyList_GetItem(obj, 0);
    *cols = (int)PyList_Size(row);

    matrix = allocate_matrix(*rows, *cols);
    if (!matrix) {
        PyErr_NoMemory();
        return NULL;
    }

    for (i = 0; i < *rows; i++) {
        row = PyList_GetItem(obj, i);
        for (j = 0; j < *cols; j++) {
            matrix[i][j] = PyFloat_AsDouble(PyList_GetItem(row, j));
            if (PyErr_Occurred()) {
                free_matrix(matrix, *rows);
                return NULL;
            }
        }
    }
    return matrix;
}

/* Converts a dynamically allocated C matrix into a 2D Python list */
static PyObject *matrix_to_py(double **matrix, int rows, int cols) {
    PyObject *outer, *row, *value;
    int i, j;

    outer = PyList_New(rows);
    if (!outer) return NULL;

    for (i = 0; i < rows; i++) {
        row = PyList_New(cols);
        if (!row) {
            Py_DECREF(outer);
            return NULL;
        }

        for (j = 0; j < cols; j++) {
            value = PyFloat_FromDouble(matrix[i][j]);
            if (!value) {
                Py_DECREF(row);
                Py_DECREF(outer);
                return NULL;
            }
            PyList_SetItem(row, j, value);
        }
        PyList_SetItem(outer, i, row);
    }
    return outer;
}

/* Python API wrapper for the sym (Similarity Matrix) C function */
static PyObject *py_sym(PyObject *self, PyObject *args) {
    PyObject *input, *result_py;
    double **points, **result;
    int n, d;

    (void)self;
    if (!PyArg_ParseTuple(args, "O", &input)) return NULL;

    points = py_to_matrix(input, &n, &d);
    if (!points) return NULL;

    result = sym(points, n, d);
    free_matrix(points, n);

    if (!result) return PyErr_NoMemory();
    result_py = matrix_to_py(result, n, n);
    free_matrix(result, n);

    return result_py;
}

/* Python API wrapper for the ddg (Degree Matrix) C function */
static PyObject *py_ddg(PyObject *self, PyObject *args) {
    PyObject *input, *result_py;
    double **points, **result;
    int n, d;

    (void)self;
    if (!PyArg_ParseTuple(args, "O", &input)) return NULL;

    points = py_to_matrix(input, &n, &d);
    if (!points) return NULL;

    result = ddg(points, n, d);
    free_matrix(points, n);

    if (!result) return PyErr_NoMemory();
    result_py = matrix_to_py(result, n, n);
    free_matrix(result, n);

    return result_py;
}

/* Python API wrapper for the norm (Normalized Matrix) C function */
static PyObject *py_norm(PyObject *self, PyObject *args) {
    PyObject *input, *result_py;
    double **points, **result;
    int n, d;

    (void)self;
    if (!PyArg_ParseTuple(args, "O", &input)) return NULL;

    points = py_to_matrix(input, &n, &d);
    if (!points) return NULL;

    result = norm(points, n, d);
    free_matrix(points, n);

    if (!result) return PyErr_NoMemory();
    result_py = matrix_to_py(result, n, n);
    free_matrix(result, n);

    return result_py;
}

/* Python API wrapper for the full symnmf optimization algorithm */
static PyObject *py_symnmf(PyObject *self, PyObject *args) {
    PyObject *h_obj, *w_obj, *result_py;
    double **H, **W, **result;
    int n, k, w_rows, w_cols;

    (void)self;
    if (!PyArg_ParseTuple(args, "OO", &h_obj, &w_obj)) return NULL;

    H = py_to_matrix(h_obj, &n, &k);
    if (!H) return NULL;

    W = py_to_matrix(w_obj, &w_rows, &w_cols);
    if (!W) {
        free_matrix(H, n);
        return NULL;
    }
    
    if (w_rows != n || w_cols != n) {
        free_matrix(H, n);
        free_matrix(W, w_rows);
        PyErr_SetString(PyExc_RuntimeError, "An Error Has Occurred");
        return NULL;
    }

    result = symnmf(H, W, n, k);
    free_matrix(H, n);
    free_matrix(W, w_rows);

    if (!result) return PyErr_NoMemory();
    result_py = matrix_to_py(result, n, k);
    free_matrix(result, n);

    return result_py;
}

/* Map Python method names to our C wrapper functions */
static PyMethodDef SymNMFMethods[] = {
    {"sym", py_sym, METH_VARARGS, "Calculate similarity matrix."},
    {"ddg", py_ddg, METH_VARARGS, "Calculate degree matrix."},
    {"norm", py_norm, METH_VARARGS, "Calculate normalized matrix."},
    {"symnmf", py_symnmf, METH_VARARGS, "Run SymNMF."},
    {NULL, NULL, 0, NULL}
};

/* Define the C extension module */
static struct PyModuleDef symnmfmodule = {
    PyModuleDef_HEAD_INIT,
    "symnmfmodule",
    NULL,
    -1,
    SymNMFMethods
};

/* Module initialization function called by Python */
PyMODINIT_FUNC PyInit_symnmfmodule(void) {
    return PyModule_Create(&symnmfmodule);
}