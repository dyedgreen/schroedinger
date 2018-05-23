#include "root.h"


// Returns the resulting zero, NOTE: This expects valid bounds!
// Returns a value of lower - 1, if limit is reached
double root_bracket(
  double (*f)(double), // Function of form f(x) = 0
  double lower,        // Initial lower bound
  double upper,        // Initial upper bound
  double accuracy,     // Desired accuracy
  int max_iteration    // Default is 2^10, if passed a negative value
) {
  double v_upper;
  double midpoint;
  if (max_iteration <= 0) {
    max_iteration = 1<<10;
  }
  for (int i = 0; i < max_iteration; i ++) {
    v_upper = f(upper);
    midpoint = lower + (upper-lower) / 2l;
    // Adjust bounds
    if (f(midpoint) * v_upper <= 0l) {
      // The zero is between MID <-> UPPER
      lower = midpoint;
    } else {
      // The zero is between LOWER <-> MID
      upper = midpoint;
    }
    // Test if we reached the desired accuracy
    if (upper - lower < accuracy) {
      return lower + (upper-lower) / 2l;
    }
  }
  // Not round, return lower bound - 1
  return lower - 1l;
}

// Returns the resulting upper bound (or a double < lower on error)
double root_expandBrackets(
  double (*f)(double), // Function of form f(x) = 0
  double lower,        // Initial lower bound, remains fixed
  double upper,        // Initial upper bound, will be adjusted
  double step,         // Used in geometric series to expand bounds
  double step_mod,     // Series used is: upper(n) = upper + step**n * step_mod
  int max_iteration    // Default is 2^10, if passed a negative value
) {
  double v_upper;
  double v_lower = f(lower);
  double i_upper = upper;
  if (max_iteration <= 0) {
    max_iteration = 1<<10;
  }
  for (int i = 0; i < max_iteration; i ++) {
    v_upper = f(upper);
    if (v_upper * v_lower <= 0) {
      return upper;
    } else {
      upper = i_upper + pow(step, (double)i) * step_mod;
    }
  }
  // Not found, return lower bound - 1
  return lower - 1l;
}