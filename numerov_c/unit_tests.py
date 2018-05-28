# Unit tests for the numerov_c package
# NOTE: The package must be built before testing

import unittest
import math
import numerov_c

# Constants and conversions (not provided by c package)
h = 6.626070040e-34
pi = 3.14159265358979
u = 1.660539040e-27
e = 1.6021766208e-19
me = 9.10938356e-31
def gammaSquared(m):
  return math.sqrt(8 * pi**2 * m / h**2)
def unscaleE(E, m):
  return E / gammaSquared(m)

class TestEigenvalues(unittest.TestCase):

  def testEnergy(self):
    # 4nm finite 14eV well
    results = [1.47, 5.74, 12.] # in eV
    def potTest(x):
      # shifted to avoid only getting odd states (ie avoid symmetries)
      if x > -5e-10 and x < -1e-10:
        return 0.0
      return 14.0 * e
    energies = numerov_c.energy(potTest, me, -2e-9, 2e-9, .1, .1, 1000, 3)
    error = 0.
    for i in range(len(energies)):
      error += (unscaleE(energies[i], me) / e - results[i])**2
    self.assertTrue(error < 0.01)


if __name__ == '__main__':
    unittest.main()