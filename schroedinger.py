# Use the numerov package we are building
# to find energy eigenvalues for arbitrary
# potentials.

import math
import numpy as np
import matplotlib.pyplot as plt
from numerov import eigenvalues, units


# Find the energies for the empty well
l = 1
s_c = 100

energies = eigenvalues.energy(lambda x: 0, units.u, 0, l, 0, 0, s_c, 10)
print(energies)
print(units.unscale(energies[0]))

x_val = np.linspace(0, l, s_c)
y_val = eigenvalues.psi(lambda x: 0, units.u, energies[0], 0, l, 0, s_c)
plt.plot(x_val, y_val, 'r:', label='1st')

y_val = eigenvalues.psi(lambda x: 0, units.u, energies[1], 0, l, 0, s_c)
plt.plot(x_val, y_val, 'c:', label='2nd')

y_val = eigenvalues.psi(lambda x: 0, units.u, energies[2], 0, l, 0, s_c)
plt.plot(x_val, y_val, 'g:', label='2nd')

plt.grid(True)
plt.legend()
plt.savefig('results/test_energies.png')