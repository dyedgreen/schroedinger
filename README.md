# Numerical Schrödinger Equations

Repository for the 2018 summer project at Imperial. This implements (and uses) a Python module for
computing numerical values of the energy eigenvalues that a particle can have in in a given potential.

## Usage

A simple example of using the module would look like this:

```
from numerov import units
import numerov_c as num

# Test potential (10Amst wide square well, 14eV deep)
def pot(x):
    if x > -1e-9 and x < 1e-9:
        return 0.
    return 14. * units.e

# Settings
l = 20e-10 # Half width of integration range
m = 9.10938e-31 # Electron mass
s = 1000 # Number of integration steps
n = 5 # Number of energy eigenstates to attempt

energies = num.energy(pot, m, -l, l, 0.1, 0.1, s, n)

# Print energies in eV
for en in energies:
    print(units.unscaleE(en) / units.e)
```

For detail documentations of each function, please see the Python source code in the `numerov` folder.
Each function is documented there. The c module only provides the `energy()` function as seen in the
above example. All c module functions have a signature identical to their Python counterparts.

## Building the C extension

Building the C extension should be super simple. Just run `make install` in the `numerov_c` folder. If
you have Python 3 installed as either `python` or `python3`, you should be good to go.

## To Do
- [x] Create repository
- [x] Add Vedin as collaborator
- [x] Write summary on QM course for supervisor meeting
- [x] Meet supervisor
- [x] Formulate goals / road-map
- [x] Get started
- [x] Verify that the eigenvalues function gives correct values for energy
- [x] Get half the energy -> investigate... (discontinuity in inf. case)
- [x] Test on less trivial examples
- [x] Think about optimization to get more energy eigenvalues (current scaling gives lots of values)
- [x] Start thinking about C implementation
- [x] ~Add more fine-grained control over root / bracket expanding functions (would not play nice with c at all...)~
- [ ] Tidy up the `schroedinger.py` file and have it output interesting results

(If you need help using Git, [this](http://try.github.io) is supposedly a nice tutorial.)

## Helpful resources

- [Wikipedia article on numerical method we will use](https://en.wikipedia.org/wiki/Numerov%27s_method)
- Book _Numerical Recipes_

## Research
- [x] Integration of ODEs (numerical recipes)
- [x] Numerov
- [x] Finite potential well
