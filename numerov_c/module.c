// This is the module file for the Python module
// NOTE: The interface is supposed to mirror the Python
//       implementation.


#include <Python.h>
#include "eigenvalues.h"

static char module_docstring[] = "This module allows to compute energy eigenvalues for a given potential";
static char energies_docstring[] = "Compute eigenvalues for a given potential and particle";

// Declare the functions we want to expose
static PyObject *energy_numerov_c(PyObject *self, PyObject *args);

// Register the functions with the module
static PyMethodDef module_methods[] = {
  // Test method
  {"energy", energy_numerov_c, METH_VARARGS, energies_docstring},
  // Guard object
  {NULL, NULL, 0, NULL},
};

// Register the module
static struct PyModuleDef numerov_module = {
  PyModuleDef_HEAD_INIT,
  "numerov_c",
  module_docstring,
  -1,
  module_methods,
};

// Initialize the module
PyMODINIT_FUNC PyInit_numerov_c(void) {
  return PyModule_Create(&numerov_module);
}


// Energy wrapper
static PyObject *energy_numerov_c(PyObject *self, PyObject *args) {
  // Get the function inputs
  PyObject *f = NULL; // No need to decref, since this is a function
  double m, x_start, x_end, y_start, y_end;
  int step_count, max_iterations;
  PyArg_ParseTuple(args, "Odddddii",
    &f, &m, &x_start, &x_end, &y_start, &y_end, &step_count, &max_iterations
  );
  // Compute energies and return result (errors are reported in the eigenvalues_energy function)
  return eigenvalues_energy(f, m, x_start, x_end, y_start, y_end, step_count, max_iterations);
}