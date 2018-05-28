# Unit tests for the numerov package

import unittest
import math
import numpy
import numerov
import root
import units
import eigenvalues

class TestNumerov(unittest.TestCase):

  # Test the integrity of the numerov step formula
  def testStep(self):
    results = [
      +0.0,
      -0.9899916597164304,
      -1.9799833194328609,
      -2.9699749791492907,
      -3.9599666388657218,
      -4.9499582985821515,
      -5.939949958298581,
      -6.92994161801501,
      -7.9199332777314435,
      -8.909924937447872
    ]
    for i in range(10):
      self.assertEqual(numerov.step(i, i*3, lambda x: 1.0, 0, 0.1), results[i])

  # Test the integration
  def testIntegrate(self):
    def sol(w, b, f_0, f_1, n):
      x = numpy.linspace(0, b, n)
      b = b / n
      f_1 = (f_1 - f_0 * math.cos(w*b)) / math.sin(w*b)
      return f_0 * numpy.cos(w * x) + f_1 * numpy.sin(w * x)
    s_c = 1000
    res = numerov.integrate(0, 1, lambda x: -4, 0, 10, step_count=s_c)
    res -= sol(2, 10, 0, 1, s_c)
    self.assertTrue(math.fabs(numpy.sum(res)) / s_c < 0.05)


class TestRoot(unittest.TestCase):

  # Test the root finding function
  def testBracket(self):
    def solve(a, b):
      # Note, we always use the + case
      return - a / 2 + math.sqrt(a**2 / 4 - b)
    def bound(a, b):
      return (-a/2, 4*a)
    for i in range(100):
      a = i/0.7
      b = -3*i
      acc = 1e-10
      sol = root.bracket(lambda x: x*x + a*x + b, bound(a,b)[0], bound(a,b)[1], acc)
      self.assertTrue(sol + acc >= solve(a,b) or sol - acc <= solve(a,b))

  # Test the bracket finding function
  def testExpandBrackets(self):
    for i in range(100):
      brack = root.expandBrackets(lambda x: x, -10**(-i), -0.9*10**(-i))
      self.assertTrue(brack[0] <= 0 and brack[1] >= 0)
    self.assertFalse(root.expandBrackets(lambda x: 1.0, 0., 1.))


class TestUnits(unittest.TestCase):

  # Test scaling and un-scaling of energies
  def testScale(self):
    # With custom mass
    for i in range(0, 10000):
      # Energy
      diff = units.unscaleE(units.scaleE(i * 1e-2, 1.4365), 1.4365) - i * 1e-2
      self.assertTrue(math.fabs(diff) < 1e-10)
      # Length
      diff = units.unscaleL(units.scaleL(i * 1e-2, 1.4365), 1.4365) - i * 1e-2
      self.assertTrue(math.fabs(diff) < 1e-10)


class TestEigenvalues(unittest.TestCase):

  def testEnergy(self):
    # 4nm finite 14eV well
    results = [1.47, 5.74, 12.] # in eV
    def potTest(x):
      # shifted to avoid only getting odd states (ie avoid symmetries)
      if x > -5e-10 and x < -1e-10:
        return 0.0
      return 14.0 * units.e
    energies = eigenvalues.energy(potTest, units.me, -2e-9, 2e-9, .1, .1, 1000, 3)
    error = 0.
    for i in range(len(energies)):
      error += (units.unscaleE(energies[i], units.me) / units.e - results[i])**2
    self.assertTrue(error < 0.01)


if __name__ == '__main__':
    unittest.main()