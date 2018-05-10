# Use the numerov package we are building
# to find energy eigenvalues for arbitrary
# potentials.

import numpy as np
import matplotlib.pyplot as plt
from numerov import numerov


# Test the numerov a bit...
def f(x): 
  return -8.1

for i in range(10, 100):
  if i % 20 == 0:
    print(i)
    x_vals = np.linspace(0, 10, i, dtype=float)
    y_vals = np.empty(len(x_vals), dtype=float)

    # Initial conditions
    y_vals[0] = 0
    y_vals[1] = 1

    for k in range(len(x_vals) - 2):
      k += 1 # This let's us use k as the same index in the formula
      y_vals[k+1] = numerov.step(y_vals[k], y_vals[k-1], f, x_vals[k], x_vals[1]-x_vals[0])

    plt.plot(x_vals, y_vals, ['r:x', 'b:x', 'g:x', 'k:x'][i//20%4], label=str(i)+" steps")
plt.legend()
plt.savefig("results/test.png")