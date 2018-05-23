#import "eigenvalues.h"

#define MASS_L 12l*UNITS_U // Mass for length scaling
#define INT_SHORTHAND(E) integrate(E, f, &m, &x_start, &x_end, &y_start, &y_end, &step_count, &step_width, &norm)

#define ROOT_STEP 1e-3l
#define ROOT_ACCURACY 1e-10l
#define ROOT_MAX_ITERATIONS 1<<10


double callPythonPotential(PyObject *f, double x) {
  // Get value from python defined potential
  PyObject *potential_py = PyObject_CallFunction(f, "d", x);
  if (potential_py == NULL) {
    // Raise an error
    PyErr_BadArgument();
    return 0l;
  }
  // Return the float value (but don't do error checking)
  return PyFloat_AsDouble(potential_py);
}

// Integration function needed by the numerov
// Why this horrid function that gets passed the local variables? Because I need the procedure twice.
double integrate(double E, PyObject *f, double *m, double *x_start, double *x_end, double *y_start, double *y_end, int *step_count, double *step_width, double *norm) {
  // k / kk def as in step fn
  double psi_k_f, psi_kk_f, psi_k_b, psi_kk_b, psi_temp, psi_max, x;
  // Initialize variables
  x = *step_width;
  psi_k_f = *norm;
  psi_k_b = *norm;
  psi_kk_f = *y_start;
  psi_kk_b = *y_end;
  psi_max = 0l;
  // Step forward / backward through integration range
  for (int k = 0; k < *step_count / 2; k ++, x += *step_width) {
    // Forward step
    psi_temp = psi_k_f;
    psi_k_f = numerov_step(
      psi_k_f,
      psi_kk_f,
      units_scaleE(callPythonPotential(f, *x_start + units_unscaleL(x - *step_width, MASS_L)), *m) - E,
      units_scaleE(callPythonPotential(f, *x_start + units_unscaleL(x, MASS_L)), *m) - E,
      units_scaleE(callPythonPotential(f, *x_start + units_unscaleL(x + *step_width, MASS_L)), *m) - E,
      *step_width
    );
    psi_kk_f = psi_temp;
    // Backward step
    psi_temp = psi_k_b;
    psi_k_b = numerov_step(
      psi_k_b,
      psi_kk_b,
      units_scaleE(callPythonPotential(f, *x_end - units_unscaleL(x - *step_width, MASS_L)), *m) - E,
      units_scaleE(callPythonPotential(f, *x_end - units_unscaleL(x, MASS_L)), *m) - E,
      units_scaleE(callPythonPotential(f, *x_end - units_unscaleL(x + *step_width, MASS_L)), *m) - E,
      *step_width
    );
    psi_kk_b = psi_temp;
    // Keep max for normalization (just from first one for convenience, since backwards is scaled after)
    if (fabs(psi_k_f) > psi_max) {
      psi_max = fabs(psi_k_f);
    }
  }
  // Adjust normalization
  *norm = 10l * *norm / psi_max;
  // Return difference in derivatives (and match backwards to be continuous)
  return (psi_k_f - psi_kk_f + (psi_k_b - psi_kk_b) * psi_k_f / psi_k_b) / *step_width;
}

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
) {
  // Test if f is a function
  if (!PyCallable_Check(f)) {
    PyErr_BadArgument();
    return NULL;
  }
  if (max_iterations <= 0) {
    max_iterations = 1<<2;
  }
  // Determine constant values
  double half_length = units_scaleL(x_end-x_start, MASS_L);
  double step_width = 2l * half_length / (double)step_count;

  // Variables for the Numerov function
  double norm = 1l;

  // Variables for the bounds and energies, reference need not be decreased, since
  // we pass on the object
  PyObject* energies_list = PyList_New(0);
  if (energies_list == NULL) {
    return PyErr_NoMemory();
  }
  double lower = 0l, upper = 0l, energy_res;
  double v_upper, v_lower, midpoint;
  int root_success;

  // Variables for the integration, bracketing and root finding
  // NOTE: This won't abort if no bounds are found, but continue to
  //       the next bracketing step.
  for (int i_energy = 0; i_energy < max_iterations; i_energy ++) {
    // Set new bounds
    lower = upper;
    upper += 1e-10l;
    // Expand energy brackets
    v_lower = INT_SHORTHAND(lower);
    root_success = 0;
    for (int i = 0; i < ROOT_MAX_ITERATIONS; i ++) {
      v_upper = INT_SHORTHAND(upper);
      if (v_upper * v_lower <= 0) {
        root_success = 1;
        break;
      } else {
        upper += ROOT_STEP;
      }
    }
    if (!root_success) {
      continue;
    }
    // Bracket energies to find roots
    root_success = 0;
    for (int i = 0; i < ROOT_MAX_ITERATIONS; i ++) {
      v_upper = INT_SHORTHAND(upper);
      midpoint = lower + (upper-lower) / 2l;
      // Adjust bounds
      if (INT_SHORTHAND(midpoint) * v_upper <= 0l) {
        // The zero is between MID <-> UPPER
        lower = midpoint;
      } else {
        // The zero is between LOWER <-> MID
        upper = midpoint;
      }
      // Test if we reached the desired accuracy
      if (upper - lower < ROOT_ACCURACY) {
        root_success = 1;
        energy_res = lower + (upper-lower) / 2l;
      }
    }
    if (!root_success) {
      continue;
    }
    // Store energy
    // Reference need not be decreased, since we pass on the values
    PyObject* energy_py = Py_BuildValue("d", energy_res);
    if (energy_py != NULL) {
      PyList_Append(energies_list, energy_py);
    } else {
      // Raise error with python
      PyErr_NoMemory();
    }
  }
  // We're done
  return energies_list;
}