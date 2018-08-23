#pragma once

#include <Python.h>
#include <math.h>
#include "units.h"
#include "numerov.h"

// Find energy eigenvalues (returns python list)
PyObject *eigenvalues_energy(
  PyObject *f,       // Potential function in Joules (takes meters)
  double m,          // Mass of particle
  double x_start,    // Start of integration range
  double x_end,
  double y_start,    // Boundary condition at x_start
  double y_end,
  int step_count,    // Number of integration steps
  int max_iterations // Max number of energies to return
);