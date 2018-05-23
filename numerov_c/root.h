#pragma once

#include <math.h>

double root_bracket(
  double (*f)(double), // Function of form f(x) = 0
  double lower,        // Initial lower bound
  double upper,        // Initial upper bound
  double accuracy,     // Desired accuracy
  int max_iteration    // Default is 2^10, if passed a negative value
);

double root_expandBrackets(
  double (*f)(double), // Function of form f(x) = 0
  double lower,        // Initial lower bound, remains fixed
  double upper,        // Initial upper bound, will be adjusted
  double step,         // Used in geometric series to expand bounds
  double step_mod,     // Series used is: upper(n) = upper + step**n * step_mod
  int max_iteration    // Default is 2^10, if passed a negative value
);

// NOTE: This is not used in the final eigenvalues function.