# Algorithm to find numeric values for the energy eigenvalues
# of a given potential function and particle.

if __name__ == 'numerov.eigenvalues':
  from . import numerov
  from . import root
  from . import units
else:
  import numerov
  import root
  import units
import numpy as np
import math

def energy(f, m, x_start, x_end, y_start, y_end, step_count=1000, max_iterations=1000):
  """Find numeric solutions to the energy eigenvalues
     of a given system.
  
  @param f              callable (potential function in Joules)
  @param m              float    (mass of particle in kg)
  @param x_start        float    (start of integration range)
  @param x_end          float
  @param y_start        float    (boundary condition on x_start)
  @param y_end          float
  @param step_count     int      (steps in numerov integration)
  @param max_iterations int      (max number of energies to return)

  @return [float] (energy eigenvalues in scaled units)
  """
  # Find constant values
  half_length = units.scaleL(x_end-x_start, m) / 2
  energies = []
  norm = 1.0
  def f_n(E, back=False):
    if back:
      return lambda x: units.scaleE(f(x_end - units.unscaleL(x, m)), m) - E
    return lambda x: units.scaleE(f(x_start + units.unscaleL(x, m)), m) - E
  def f_r(E):
    nonlocal norm
    # Perform Numerov integration
    y_vals_f = numerov.integrate(y_start, norm, f_n(E, False), 0, half_length, step_count=step_count//2)
    y_vals_b = numerov.integrate(y_start, norm, f_n(E, True), 0, half_length, step_count=step_count//2)
    # Make continuous
    cont_factor = y_vals_f[len(y_vals_f)-1] / y_vals_b[len(y_vals_b)-1]
    # Compute derivatives
    d_f = (y_vals_f[len(y_vals_f)-1] - y_vals_f[len(y_vals_f)-2]) * step_count / 2 / half_length
    d_b = (y_vals_b[len(y_vals_b)-2] - y_vals_b[len(y_vals_b)-1]) * cont_factor * step_count / 2 / half_length
    # Improve normalization
    norm = 10 * norm / np.max(np.abs(y_vals_f))
    return d_f - d_b
  # Find energies using the root / shooting method
  bounds = (0,0)
  i = 0
  while bounds != False and i < max_iterations:
    print("Iteration " + str(i))
    bounds = (bounds[1], bounds[1]+1e-10)
    # Note: Think about better / more intelligent step sizes vs how many steps needed
    bounds = root.expandBrackets(f_r, bounds[0], bounds[1], 1.+1e-3, 1e-1, 1<<15)
    if bounds == False:
      break
    energy = root.bracket(f_r, bounds[0], bounds[1])
    if energy != False:
      energies.append(energy)
    elif i == 0:
      # If this is the first, add one to match c implementation
      max_iterations += 1
    i += 1
  return energies

def psi(f, m, E, x_start, x_end, y_start, y_end, step_count=1000):
  """Find the plot for Psi for a given energy
  
  @param f              callable (potential function in Joules)
  @param m              float    (mass of particle in kg)
  @param E              float    (scaled energy eigenvalue to plot)
  @param x_start        float    (start of integration range)
  @param x_norm         float    (normalization parameter, see energy)
  @param y_start        float    (boundary condition on x_start)
  @param y_end          float
  @param step_count     int      (steps in numerov integration)

  @return list of float (values of psi, not normalized)
  """
  # Find constant values
  half_length = units.scaleL(x_end-x_start, m) / 2
  energies = []
  def f_n(E, back=False):
    if back:
      return lambda x: units.scaleE(f(x_end - units.unscaleL(x, m)), m) - E
    return lambda x: units.scaleE(f(x_start + units.unscaleL(x, m)), m) - E
  # Perform Numerov integration
  y_vals_f = numerov.integrate(y_start, 1.0, f_n(E, False), 0, half_length, step_count=step_count//2)
  y_vals_b = numerov.integrate(y_end, 1.0, f_n(E, True), 0, half_length, step_count=step_count//2)
  # Make continuous
  y_vals_b *= y_vals_f[len(y_vals_f)-1] / y_vals_b[len(y_vals_b)-1]
  # Return normalized result
  return np.concatenate((y_vals_f, np.flip(y_vals_b, 0))) / np.max(np.abs(y_vals_f))