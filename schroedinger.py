# Use the numerov package we are building
# to find energy eigenvalues for arbitrary
# potentials.

import math
import numpy as np
import matplotlib.pyplot as plt
from numerov import eigenvalues, units
import numerov_c as num

# General settings
l = 2e-9

# Potential functions
pot = []
mass = []
n = []
names = []

# 14eV deep, 4nm wide square well
def potSquareWell(x):
    if x > -5e-10 and x < -1e-10:
        return 0.
    return 14. * units.e
pot.append(potSquareWell)
mass.append(units.me)
n.append(3)
names.append("Electron in 14eV Square Well")

# Slanted potential well
def potQHM(x):
    return 14. * units.e * units.e * (x/2e-10)**2
pot.append(potQHM)
mass.append(units.me)
n.append(3)
names.append("Electron in Quantum Harmonic Oscillator")

# Slanted potential well
def potSlantedWell(x):
    if x > -5e-10 and x < -1e-10:
        return 14. * units.e * (x - 5e-10) / 4e-10
    return 14. * units.e
pot.append(potSlantedWell)
mass.append(units.me)
n.append(1)
names.append("Electron in Slanted 14eV Well")

########################################
# Solve, plot and print all potentials #
########################################

def underline(length = 0):
    line = ""
    for i in range(length):
        line += "="
    return line

# Print all the graphs
x = np.linspace(-l, l, 1000)
colors = ['r:','g:','b:','c:', 'm:', 'y:']
for i in range(len(pot)):
    # Find energies an print values
    energies = num.energy(pot[i], mass[i], -l, l, 0.1, 0.1, 1000, n[i])
    print("\nResults for " + names[i])
    print(underline(12+len(names[i])))
    print()
    # Graph potential
    pot_y = np.empty(len(x))
    for j in range(len(x)):
        pot_y[j] = pot[i](x[j])
    pot_y *= .5 / np.max(pot_y)
    plt.plot(x, pot_y, 'k-')
    # Graph wave functions
    for j in range(len(energies)):
        print("["+str(j)+"] " + str(units.unscaleE(energies[j], mass[i]) / units.e) + " eV")
        psi = eigenvalues.psi(pot[i], mass[i], energies[j], -l, l, 0.1, 0.1, 1000)
        plt.plot(x, psi, colors[(j) % len(colors)], label=str(units.unscaleE(energies[j], mass[i]) / units.e)[:5] + " eV")
    # Save plot
    plt.grid(True)
    plt.legend()
    plt.xlabel("meters")
    plt.savefig("results/"+"_".join(names[i].lower().split()))
    plt.clf()