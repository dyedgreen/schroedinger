# Numerov implementation in Python.

import numpy as np

def step(Y_k, Y_kk, f, x_k, step_width):
  """Advance Numerov integration by one step

  @param Y_k        float    (value at k)
  @param Y_kk       float    (value at k-1)
  @param f          callable (function of form Y'' = f(x)*Y)
  @param x_k        float    (x value at k)
  @param step_width float

  @return float (value at k+1)
  """
  pre = 1 - step_width**2 * f(x_k+step_width) / 12
  res = Y_k * (step_width**2 * f(x_k) * 5 / 6 + 2)
  res += Y_kk * (step_width**2 * f(x_k-step_width) / 12 - 1)
  return res / pre

def integrate(Y_0, Y_1, f, x_start, x_end, step_count=1000):
  """Integrate a function using the numerov step
  
  @param Y_0        float    (value at x_start)
  @param Y_1        float    (value at x_1 ie directly after start)
  @param f          callable (see step)
  @param x_start    float    (start of integration range)
  @param x_end      float    (end of integration range)
  @param step_count int

  @return numpy.array
  """
  y_values = np.empty(step_count, dtype=float)
  y_values[0] = Y_0
  y_values[1] = Y_1
  step_width = (x_end - x_start) / step_count
  k = 0
  while k < step_count - 2:
    k += 1
    y_values[k+1] = step(y_values[k], y_values[k-1], f, x_start + k*step_width, step_width)
  return y_values