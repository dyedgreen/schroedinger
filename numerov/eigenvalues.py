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

  @return list of float (energy eigenvalues in scaled units)
  """
  # Scale units
  x_start = units.scaleL(x_start, m)
  x_end = units.scaleL(x_end, m)
  # These will be in scaled units
  energies = []
  def f_n(E):
    # Callable for numerov at given energy
    return lambda x: units.scaleE(f(units.unscaleL(x)), m) - E
  def f_r(E):
    # Callable for root
    y_vals = numerov.integrate(y_start, 1.0, f_n(E), x_start, x_end, step_count)
    return y_vals[len(y_vals)-1] - y_end
  # Find energies using the root / shooting method
  bounds = (0,0)
  i = 0
  while bounds != False and i < max_iterations:
    i += 1
    bounds = (bounds[1], bounds[1]+1e-10)
    # Note: Think about better / more intelligent step sizes vs how many steps needed
    bounds = root.expandBrackets(f_r, bounds[0], bounds[1], 1+1e-4, 1e-2, 1<<12)
    if bounds == False:
      break
    energy = root.bracket(f_r, bounds[0], bounds[1])
    if energy != False:
      energies.append(energy)
  return energies

def psi(f, m, E, x_start, x_end, y_start, step_count=1000):
  """Find the plot for Psi for a given energy
  
  @param f              callable (potential function in Joules)
  @param m              float    (mass of particle in kg)
  @param E              float    (energy eigenvalue to plot)
  @param x_start        float    (start of integration range)
  @param x_end          float
  @param y_start        float    (boundary condition on x_start)
  @param step_count     int      (steps in numerov integration)

  @return list of float (values of psi, not normalized)
  """
  # Scale units
  x_start = units.scaleL(x_start, m)
  x_end = units.scaleL(x_end, m)
  # Perform integration
  return numerov.integrate(y_start, 1.0, lambda x: units.scaleE(f(units.unscaleL(x)), m) - E, x_start, x_end, step_count)