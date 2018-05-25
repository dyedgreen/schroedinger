# Use the numerov package we are building
# to find energy eigenvalues for arbitrary
# potentials.

import math
import numpy as np
import matplotlib.pyplot as plt
from numerov import eigenvalues, units
import numerov_c as num_c


# Find the energies for the empty well
l = 10e-10
l_print = l
l_inf = 4e-10
mass = 9.10938e-31 #2*units.u
s_c = 500

initial = 0.1

def pot(x):
  #return 0.;
  #return np.sin(x*5e10) * 8.0e-20
  #if x > -2e-10 and x < -1e-10:
  #  return np.sin(x*5e10) * 14.0e-3 * 1.6e-19
  if x > -2e-10 and x < 2e-10:
    return 0.0
  #elif x >= 2e-10:
  #  return 20.0e-3 * 1.6e-19
  return 14.0 * 1.6e-19




#energies = eigenvalues.energy(pot, mass, -l, l, initial, initial, s_c, 5)
#print(energies)
energies = num_c.energy(pot, mass, -l, l, initial, initial, s_c, 10)
print(energies)

#correction = (units.h**2 / 8 / mass / (l_inf)**2) / units.unscaleE(energies[0], mass)
#print(correction)

x_pot = np.linspace(-l, l_print, s_c)
y_pot = np.empty(len(x_pot))
for i in range(len(x_pot)):
  y_pot[i] = pot(x_pot[i])
y_pot *= .5 / np.max(y_pot)
plt.plot(x_pot, y_pot, 'k-')

i = -1
col = ['r:','g:','b:','c:', 'm:', 'y:']
for en in energies:
  i += 1
  # Make graph pretty
  #if i != 1:
  #  continue

  x_val = np.linspace(-l_print, l_print, s_c)
  y_val = eigenvalues.psi(pot, mass, en, -l_print, l_print, initial, initial, s_c)
  #y_val_c = eigenvalues.psi(pot, mass, energies_c[i], -l_print, l_print, initial, initial, s_c)
  plt.plot(x_val, y_val, col[(i) % len(col)], label=str(i))
  #plt.plot(x_val, y_val_c, col_c[(i) % len(col_c)], label="c " + str(i))

  print(
    "[" +str(i)+ "]",
    units.unscaleE(energies[i], mass) / 1.6e-19,
    (i+1)**2 * units.h**2 / 8 / mass / (l_inf)**2,
    units.unscaleE(energies[i], mass) - (i+1)**2 * units.h**2 / 8 / mass / (l_inf)**2
  )
  #print(y_val)

# # x_val = np.linspace(-l, l, s_c)
# # y_val = eigenvalues.psi(pot, mass, 0.73294830903, -l, l, 0, 0, s_c)
# # plt.plot(x_val, y_val, 'r:', label="Test")

# #plt.ylim((-4,4))

plt.grid(True)
plt.legend()
plt.savefig('results/test_energies.png')