#define PY_SSIZE_T_CLEAN
#include <Python.h>
#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include "symnmf.h"
#include "matrixutils.h"
#define ERROR "An Error Has Occurred\n"

double** matrix_malloc_py(PyObject *obj, int rows, int cols) {
    double **points;
    int i, j;
    if (!PyList_Check(obj)) {
        printf(ERROR);
        return NULL;
    }
    points = calloc(rows, sizeof(double*));
    if (points == NULL) {
        printf(ERROR);
        return NULL;
    }
    for (i = 0; i < rows; i++) {
        PyObject *row = PyList_GetItem(obj, i);
        if (!PyList_Check(row)) {
            printf(ERROR);
            matrix_memory_free(points, i);
            return NULL;
        }
        points[i] = calloc(cols, sizeof(double));
        if (points[i] == NULL) {
            PyErr_SetString(PyExc_MemoryError, ERROR);
            matrix_memory_free(points, i);
            return NULL;
        }
        for (j = 0; j < cols; j++) {
            PyObject *item = PyList_GetItem(row, j);
            if (!PyFloat_Check(item) && !PyLong_Check(item)) {
                printf(ERROR);
                matrix_memory_free(points, i + 1);
                return NULL;
            }
            points[i][j] = PyFloat_AsDouble(item);
        }
    }
    return points;
}

PyObject* matrix_to_pylist(double **matrix, int rows, int cols) {
    PyObject *PythonResult = PyList_New(rows);
    if (PythonResult == NULL) {
        PyErr_SetString(PyExc_MemoryError, ERROR);
        return NULL;
    }
    for (int i = 0; i < rows; i++) {
        PyObject *PythonRow = PyList_New(cols);
        if (PythonRow == NULL) {
            PyErr_SetString(PyExc_MemoryError, ERROR);
            Py_DECREF(PythonResult);
            return NULL;
        }
        for (int j = 0; j < cols; j++) {
            PyList_SET_ITEM(PythonRow, j, PyFloat_FromDouble(matrix[i][j]));
        }
        PyList_SetItem(PythonResult, i, PythonRow);
    }
    return PythonResult;
}

static PyObject* py_sym(PyObject *self, PyObject *args) {
    PyObject *PythonPoints, *first_row;
    double **points, **result;
    int N, d;
    if (!PyArg_ParseTuple(args, "O", &PythonPoints)) {
        return NULL;
    }
    N = PyList_Size(PythonPoints);
    if(N == 0) {
        printf(ERROR);
        return NULL;
    }
    first_row = PyList_GetItem(PythonPoints, 0);
    d = PyList_Size(first_row);
    points = matrix_malloc_py(PythonPoints, N, d);
    if (points == NULL) {
        return NULL;
    }
    result = sym(points, N, d);
    if (result == NULL) {
        matrix_memory_free(points, N);
        printf(ERROR);
        return NULL;
    }
    PyObject *PythonResult = matrix_to_pylist(result, N, N);
    matrix_memory_free(points, N);
    matrix_memory_free(result, N);
    return PythonResult;
}

static PyObject* py_ddg(PyObject *self, PyObject *args) {
    PyObject *PythonPoints, *first_row;
    double **points, **result;
    int N, d;
    if (!PyArg_ParseTuple(args, "O", &PythonPoints)) {
        return NULL;
    }
    N = PyList_Size(PythonPoints);
    if(N == 0) {
        printf(ERROR);
        return NULL;
    }
    first_row = PyList_GetItem(PythonPoints, 0);
    d = PyList_Size(first_row);
    points = matrix_malloc_py(PythonPoints, N, d);
    if (points == NULL) {
        return NULL;
    }
    result = ddg(points, N, d);
    if (result == NULL) {
        matrix_memory_free(points, N);
        PyErr_SetString(PyExc_RuntimeError, ERROR);
        return NULL;
    }
    PyObject *PythonResult = matrix_to_pylist(result, N, N);
    matrix_memory_free(points, N);
    matrix_memory_free(result, N);
    return PythonResult;
}

static PyObject* py_norm(PyObject *self, PyObject *args) {
    PyObject *PythonPoints, *first_row;
    double **points, **result;
    int N, d;
    if (!PyArg_ParseTuple(args, "O", &PythonPoints)) {
        return NULL;
    }
    N = PyList_Size(PythonPoints);
    if(N == 0) {
        printf(ERROR);
        return NULL;
    }
    first_row = PyList_GetItem(PythonPoints, 0);
    d = PyList_Size(first_row);
    points = matrix_malloc_py(PythonPoints, N, d);
    if (points == NULL) {
        return NULL;
    }
    result = norm(points, N, d);
    if (result == NULL) {
        matrix_memory_free(points, N);
        PyErr_SetString(PyExc_RuntimeError, ERROR);
        return NULL;
    }
    PyObject *PythonResult = matrix_to_pylist(result, N, N);
    matrix_memory_free(points, N);
    matrix_memory_free(result, N);
    return PythonResult;
}

static PyObject* symnmf(PyObject *self, PyObject *args) {
    PyObject *PythonNormMatrix;
    PyObject *PythonInitial;
    int N, d, k, max_iter = 300;
    double **norm_matrix, **initial, **result, eps = 1e-4;
    if (!PyArg_ParseTuple(args, "OiiiO", &PythonNormMatrix, &N, &d, &k, &PythonInitial)) {
        return NULL;
    }
    initial = matrix_malloc_py(PythonInitial, N, k);
    if (initial == NULL) {
        return NULL;
    }
    norm_matrix = matrix_malloc_py(PythonNormMatrix, N, N);
    if (norm_matrix == NULL) {
        matrix_memory_free(initial, N);
        return NULL;
    }
    result = optimize_h(initial, norm_matrix, N, k, max_iter, eps);
    if (result == NULL) {
        matrix_memory_free(initial, N);
        matrix_memory_free(norm_matrix, N);
        PyErr_SetString(PyExc_RuntimeError, ERROR);
        return NULL;
    }
    PyObject *PythonResult = matrix_to_pylist(result, N, k);
    matrix_memory_free(initial, N);
    matrix_memory_free(norm_matrix, N);
    matrix_memory_free(result, N);
    return PythonResult;
}

static PyMethodDef mysymnmf_methods[] = {
    {"symnmf", (PyCFunction)symnmf, METH_VARARGS, "The full NMF algorithm."},
    {"sym", (PyCFunction)py_sym, METH_VARARGS, "Returns the similarity matrix of the given points."},
    {"ddg", (PyCFunction)py_ddg, METH_VARARGS, "Returns the diagonal degree matrix of the given points."},
    {"norm", (PyCFunction)py_norm, METH_VARARGS, "Returns the normalized similarity matrix of the given points."},
    {NULL, NULL, 0, NULL}
};

static struct PyModuleDef mysymnmf_module = {
    PyModuleDef_HEAD_INIT,
    "mysymnmf",
    "The full NMF algorithm.",
    -1,
    mysymnmf_methods
};

PyMODINIT_FUNC PyInit_mysymnmf(void) {
    PyObject *m;
    m = PyModule_Create(&mysymnmf_module);
    if (!m) {
        return NULL;
    }
    return m;
}
