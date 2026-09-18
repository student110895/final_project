#define PY_SSIZE_T_CLEAN
#include <Python.h>
#include "symnmf.h"


static double **py_to_matrix(PyObject *obj, int *rows, int *cols)
{
    double **matrix;
    PyObject *row;
    Py_ssize_t i;
    Py_ssize_t j;

    *rows = (int)PyList_Size(obj);
    row = PyList_GetItem(obj, 0);
    *cols = (int)PyList_Size(row);

    matrix = allocate_matrix(*rows, *cols);
    if (matrix == NULL) {
        PyErr_NoMemory();
        return NULL;
    }

    for (i = 0; i < *rows; i++) {
        row = PyList_GetItem(obj, i);

        for (j = 0; j < *cols; j++) {
            matrix[i][j] =
                PyFloat_AsDouble(PyList_GetItem(row, j));

            if (PyErr_Occurred()) {
                free_matrix(matrix, *rows);
                return NULL;
            }
        }
    }

    return matrix;
}


static PyObject *matrix_to_py(double **matrix, int rows, int cols)
{
    PyObject *outer;
    PyObject *row;
    PyObject *value;
    int i;
    int j;

    outer = PyList_New(rows);

    if (outer == NULL)
        return NULL;

    for (i = 0; i < rows; i++) {
        row = PyList_New(cols);

        if (row == NULL) {
            Py_DECREF(outer);
            return NULL;
        }

        for (j = 0; j < cols; j++) {
            value = PyFloat_FromDouble(matrix[i][j]);

            if (value == NULL) {
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


static PyObject *py_sym(PyObject *self, PyObject *args)
{
    PyObject *input;
    PyObject *result_py;
    double **points;
    double **result;
    int n;
    int d;

    (void)self;

    if (!PyArg_ParseTuple(args, "O", &input))
        return NULL;

    points = py_to_matrix(input, &n, &d);

    if (points == NULL)
        return NULL;

    result = sym(points, n, d);
    free_matrix(points, n);

    if (result == NULL)
        return PyErr_NoMemory();

    result_py = matrix_to_py(result, n, n);
    free_matrix(result, n);

    return result_py;
}


static PyObject *py_ddg(PyObject *self, PyObject *args)
{
    PyObject *input;
    PyObject *result_py;
    double **points;
    double **result;
    int n;
    int d;

    (void)self;

    if (!PyArg_ParseTuple(args, "O", &input))
        return NULL;

    points = py_to_matrix(input, &n, &d);

    if (points == NULL)
        return NULL;

    result = ddg(points, n, d);
    free_matrix(points, n);

    if (result == NULL)
        return PyErr_NoMemory();

    result_py = matrix_to_py(result, n, n);
    free_matrix(result, n);

    return result_py;
}


static PyObject *py_norm(PyObject *self, PyObject *args)
{
    PyObject *input;
    PyObject *result_py;
    double **points;
    double **result;
    int n;
    int d;

    (void)self;

    if (!PyArg_ParseTuple(args, "O", &input))
        return NULL;

    points = py_to_matrix(input, &n, &d);

    if (points == NULL)
        return NULL;

    result = norm(points, n, d);
    free_matrix(points, n);

    if (result == NULL)
        return PyErr_NoMemory();

    result_py = matrix_to_py(result, n, n);
    free_matrix(result, n);

    return result_py;
}


static PyObject *py_symnmf(PyObject *self, PyObject *args)
{
    PyObject *h_obj;
    PyObject *w_obj;
    PyObject *result_py;
    double **H;
    double **W;
    double **result;
    int n;
    int k;
    int w_rows;
    int w_cols;

    (void)self;

    if (!PyArg_ParseTuple(args, "OO", &h_obj, &w_obj))
        return NULL;

    H = py_to_matrix(h_obj, &n, &k);

    if (H == NULL)
        return NULL;

    W = py_to_matrix(w_obj, &w_rows, &w_cols);

    if (W == NULL) {
        free_matrix(H, n);
        return NULL;
    }

    result = symnmf(H, W, n, k);

    free_matrix(H, n);
    free_matrix(W, w_rows);

    if (result == NULL)
        return PyErr_NoMemory();

    result_py = matrix_to_py(result, n, k);
    free_matrix(result, n);

    return result_py;
}


static PyMethodDef SymNMFMethods[] = {
    {"sym", py_sym, METH_VARARGS, "Calculate similarity matrix."},
    {"ddg", py_ddg, METH_VARARGS, "Calculate degree matrix."},
    {"norm", py_norm, METH_VARARGS, "Calculate normalized matrix."},
    {"symnmf", py_symnmf, METH_VARARGS, "Run SymNMF."},
    {NULL, NULL, 0, NULL}
};


static struct PyModuleDef symnmfmodule = {
    PyModuleDef_HEAD_INIT,
    "symnmfmodule",
    NULL,
    -1,
    SymNMFMethods
};


PyMODINIT_FUNC PyInit_symnmfmodule(void)
{
    return PyModule_Create(&symnmfmodule);
}