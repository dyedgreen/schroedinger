# Root finding implementation in Python

def bracket(f, lower, upper, accuracy=1e-10, max_iteration=1<<10):
  """Find root of function within bracket limits.

  @param f             callable (function of form f(x) = 0)
  @param lower         float    (initial lower bound)
  @parma upper         float    (initial upper bound)
  @param accuracy      float    (desired accuracy, default is 2^-5)
  @param max_iteration int      (default is 2^10)

  @return float or False 
  """
  v_upper = 0.
  v_lower = 0.
  midpoint = 0.
  count = 0
  while count < max_iteration:
    count += 1
    v_upper = f(upper)
    v_lower = f(lower)
    midpoint = lower + (upper-lower) / 2
    v_midpoint = f(midpoint)
    if v_midpoint * v_upper < 0:
      # MID <-> UPPER
      lower = midpoint
    elif v_midpoint * v_lower < 0:
      # LOWER <-> MID
      upper = midpoint
    elif v_midpoint * v_upper == 0:
      return midpoint if v_midpoint == 0 else upper
    elif v_lower == 0:
      return lower
    else:
      # The range was bad
      return False
    if upper - lower <= accuracy:
      return lower + (upper - lower) / 2
  # Reached iteration limit
  return False

def expandBrackets(f, lower, upper, step=1.1, step_mod=1.0, max_iteration=1<<10):
  """Expand brackets upwards using a geometric series

  @param f             callable (function of form f(x) = 0)
  @param lower         float    (initial guess, not varied)
  @param upper         float    (initial guess, is expanded)
  @param step          float    (used in geometric series)
  @param step_mod      float    (used to modify the geometric series: deltaX_n = step_mod * step ** n)
  @param max_iteration int      (default 2^10)

  @return float or False
  """
  v_upper = 0.
  v_lower = f(lower)
  count = 0
  while count < max_iteration:
    count += 1
    v_upper = f(upper)
    if v_upper * v_lower < 0:
      return (lower, upper)
    else:
      upper += step_mod * step ** count
  return False