# Use the numerov package we are building
# to find energy eigenvalues for arbitrary
# potentials.

import math
import numpy as np
import matplotlib.pyplot as plt
from numerov import eigenvalues, units


# Find the energies for the empty well
l = 2e-10
mass = 2*units.u
s_c = 100

energies = eigenvalues.energy(lambda x: 0, mass, 0, l, 0, 0, s_c, 5)
print(energies)

i = 0
col = ['r:','g:','b:','c:','k:']
for en in energies:
  i += 1
  x_val = np.linspace(0, l, s_c)
  y_val = eigenvalues.psi(lambda x: 0, mass, en, 0, l, 0, s_c)
  plt.plot(x_val, y_val, col[(i - 1) % len(col)], label=str(i))

  print(units.unscaleE(energies[i-1]), i**2 * units.h**2 / 8 / mass / l**2)

plt.grid(True)
plt.legend()
plt.savefig('results/test_energies.png')