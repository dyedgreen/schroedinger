# Unit tests for the numerov package

import unittest
import numerov

class TestNumerov(unittest.TestCase):

  # Test the integrity of the numerov step formula
  def testStep(self):
    results = [
      0.0,
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

if __name__ == '__main__':
    unittest.main()