# Use the numerov package we are building
# to find energy eigenvalues for arbitrary
# potentials.

import math
import numpy as np
import matplotlib.pyplot as plt
from numerov import eigenvalues, units


# Find the energies for the empty well
l = 5e-10
mass = 2*units.u
s_c = 100

def pot(x):
  if x > -2e-10 and x < 2e-10:
    return 0.0
  return 14.0 * 1.6e-19


energies = eigenvalues.energy(pot, mass, -l, l, 0, 0, s_c, 5)
print(energies)

i = 0
col = ['r:','g:','b:','c:','k:']
for en in energies:
  i += 1
  x_val = np.linspace(-l, l, s_c)
  y_val = eigenvalues.psi(pot, mass, en, -l, l, 0, 0, s_c)
  plt.plot(x_val, y_val, col[(i - 1) % len(col)], label=str(i))

  print(units.unscaleE(energies[i-1]), i**2 * units.h**2 / 8 / mass / l**2)

# x_val = np.linspace(-l, l, s_c)
# y_val = eigenvalues.psi(pot, mass, 0.73294830903, -l, l, 0, 0, s_c)
# plt.plot(x_val, y_val, 'r:', label="Test")

plt.grid(True)
plt.legend()
plt.savefig('results/test_energies.png')