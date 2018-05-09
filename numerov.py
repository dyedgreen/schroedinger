# An initial attempt at implementing Numerov's
# algorithm in Python. (Once we got it working
# properly, maybe look at using c for speed, or
# getting it to work in 3-D / N-D)

# NOTES:
# This is very crude! Specifically, it only deals
# with the case of linear 'f' and it probably not
# incredibly fast. (Also, Numerov is limited to 1D)

import numpy as np
import math
import matplotlib.pyplot as plt

# Equations need form: y'' = f(x, y)


# Schroedingers EQ (inf well, the 'f' from above -> this is a simple harmonic oscillator)
def schroedinger (x, y, c=1.0):
  return -c*y

# Crude Numerov (approx next step, assuming linear 'f')
def numerov(x, X, y, Y, f, c=1.0):
  h = X - x
  return (h*h / 12 * (10*f(X,Y,c) + f(x,y,c)) + 2*Y - y) / (1 - h*h / 12 * (-c))

# Find result for a range in between two numbers
def find(r, y, Y, f, c=1.0):
  res = np.empty(len(r), float)
  res[0] = y
  res[1] = Y
  for i in range(len(r)-2):
    res[i+2] = numerov(r[i], r[i+1], res[i], res[i+1], f, c)
  return res

# Try it :)
r = np.linspace(0, 10, 100)
res = find(r, 0, 0.1, schroedinger, 8.1)

plt.plot(r, res, 'b:x', label="Solution")
plt.legend()
plt.savefig("test.png")
