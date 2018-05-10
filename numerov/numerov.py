# Numerov implementation in Python.

# Imports (if needed later) go here...

# Numerov step
#
# Computes Y(k+1) (Or Y(k-1),
# since it's symmetric)
#
# Takes Yk, Yk-1, fn, xk, step
def step(Y_k, Y_kk, f, x_k, step_width):
  pre = 1 - step_width**2 * f(x_k+step_width) / 12
  res = Y_k * (step_width**2 * f(x_k) * 5 / 6 + 2)
  res += Y_kk * (step_width**2 * f(x_k-step_width) / 12 - 1)
  return res / pre