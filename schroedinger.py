# Use the numerov package we are building
# to find energy eigenvalues for arbitrary
# potentials.

import math
import numpy as np
import matplotlib.pyplot as plt
from numerov import numerov


# Test the numerov a bit...
s_c = 100
y_val = numerov.integrate(0, 1, lambda x: -4, 0., 10., step_count=s_c)
x_val = np.linspace(0, 10, s_c)

def sol(w, b, f_0, f_1, n):
  x = np.linspace(0, b, n)
  b = b / n
  f_1 = (f_1 - f_0 * math.cos(w*b)) / math.sin(w*b)
  return f_0 * np.cos(w * x) + f_1 * np.sin(w * x)

plt.plot(x_val, y_val, 'r:x')
plt.plot(x_val, sol(2, 10, 0, 1, s_c), 'b:x')
plt.savefig('results/test.png')

print(np.sum(y_val - sol(2, 10, 0, 1, s_c)))
print(np.sum(y_val - sol(2, 10, 0, 1, s_c)) / s_c)