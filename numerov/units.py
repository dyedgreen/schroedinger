# Units helps with scaling energies into
# ranges that are more convenient for the
# Numerov algorithm.

h = 6.626070040e-34 # Planck's Constant according to NIST
pi = 3.14159265358979
u = 1.660539040e-27 # Atomic mass unit

def scale(E, m=12*u):
  """Scale an energy value into a nice range

  @param E float (energy value in Joules)
  @param m float (mass of particle involved)

  @return float
  """
  return E * 4 * pi * m / h

def unscale(E, m=12*u):
  """Retrieve the energy value in Joules

  @param E float (energy value in scaled units)
  @param m float (see scale)

  @return float
  """
  return E * h / 4 / pi / m