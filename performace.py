# Test the performance of the numerov / numerv_c implementations

from numerov import eigenvalues, units
import numerov_c as num_c
from time import clock

# Test potential (10Amst wide square well, 14eV deep)
def pot(x):
    if x > -1e-9 and x < 1e-9:
        return 0.
    return 14. * units.e

# Settings for both runs
l = 20e-10 # Width of integration range
m = 9.10938e-31 # Electron mass
s = 1000 # Number of integration steps
n = 100 # Number of energy eigenstates to attempt

# Python implementation
print("Starting python implementation")

time_p = clock()
energies_p = eigenvalues.energy(pot, m, -l, l, 0.1, 0.1, s, n)
time_p = clock() - time_p

print("Starting c implementation")

time_c = clock()
energies_c = num_c.energy(pot, m, -l, l, 0.1, 0.1, s, n)
time_c = clock() - time_c

print("\n=======\nResults\n=======\n")
print("Python : " + str(time_p))
print("C      : " + str(time_c))

print("C saved an absolute time of " + str(time_p - time_c) + " s")
print("C is thus " + str(time_p / time_c) + " times as fast!")