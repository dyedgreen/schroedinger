# Use the numerov package we are building
# to find energy eigenvalues for arbitrary
# potentials.

import math
import numpy as np
import matplotlib.pyplot as plt
from numerov import eigenvalues, units


# Find the energies for the empty well
l = 5e-10
l_print = 5e-10#2.04e-10
mass = 2*units.u
s_c = 300

initial = 0.1

def pot(x):
  #return np.sin(x*5e10) * 8.0e-20
  if x > -2e-10 and x < 2e-10:
    return 0.0
  elif x >= 2e-10:
    return 10.0 * 1.6e-19
  return 14.0 * 1.6e-19


energies = eigenvalues.energy(pot, mass, -l, l, initial, initial, s_c, 5)
print(energies)

correction = (units.h**2 / 8 / mass / (l*2)**2) / units.unscaleE(energies[0], mass)
print(correction)

# x_pot = np.linspace(-l, l_print, s_c)
# y_pot = pot(x_pot)
# y_pot *= 5 / np.max(y_pot)
# plt.plot(x_pot, y_pot, '-k')

i = -1
col = ['r:','g:','b:','c:','k:']
for en in energies:
  i += 1
  # Make graph pretty
  if i == 2:
    continue

  x_val = np.linspace(-l_print, l_print, s_c)
  y_val = eigenvalues.psi(pot, mass, en, -l_print, l_print, initial, initial, s_c)
  plt.plot(x_val, y_val, col[(i) % len(col)], label=str(i))

  print(
    "[" +str(i)+ "]",
    units.unscaleE(energies[i], mass) * correction,
    (i+1)**2 * units.h**2 / 8 / mass / (l*2)**2,
    units.unscaleE(energies[i], mass) * correction - (i+1)**2 * units.h**2 / 8 / mass / (l*2)**2
  )
  #print(y_val)

# x_val = np.linspace(-l, l, s_c)
# y_val = eigenvalues.psi(pot, mass, 0.73294830903, -l, l, 0, 0, s_c)
# plt.plot(x_val, y_val, 'r:', label="Test")

#plt.ylim((-4,4))

plt.grid(True)
plt.legend()
plt.savefig('results/test_energies.png')