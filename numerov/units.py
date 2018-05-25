# Units helps with scaling energies into
# ranges that are more convenient for the
# Numerov algorithm.

import math

h = 6.626070040e-34 # Planck's Constant according to NIST
pi = 3.14159265358979
u = 1.660539040e-27 # Atomic mass unit
e = 1.6021766208e-19 # Fundamental charge 

# We have gamma and gamma squared, as sqrt looses precision
def gammaSquared(m):
  return math.sqrt(8 * pi**2 * m / h**2)
def gamma(m):
  return math.sqrt(gammaSquared(m))

def scaleE(E, m):
  """Scale an energy value into a nice range

  @param E float (energy value in Joules)
  @param m float (mass of particle involved)

  @return float
  """
  return E * gammaSquared(m)

def unscaleE(E, m):
  """Retrieve the energy value in Joules

  @param E float (energy value in scaled units)
  @param m float (see scale)

  @return float
  """
  return E / gammaSquared(m)

def scaleL(L, m):
  """Scale a length value into a nice range

  @param L float (length value in meters)
  @param m float (mass of particle involved)

  @return float
  """
  return L * gamma(m)

def unscaleL(L, m):
  """Retrieve the length value in meters

  @param L float (length value in scaled units)
  @param m float (see scale)

  @return float
  """
  return L / gamma(m)