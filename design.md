# Algorithm design choices

-> Goal, given a potential function, we want to find the energy eigenvalues

## Numerov
- Use an Numerov step method to carry out the integration
- Use 'shooting' method to match boundary conditions (NR )
- Use 'Van Wijngaarden-Dekker-Brent' method to find root (i.e. boundary) while shooting (NR 9.3)

### Finding the energy (i.e finding roots)
Basically, the function will be

```
  F(E) - B = 0
```
where B is the value at the boundary and F is the integral computed using Numerov. (Note: The second
Psi/Y value is arbitrary as I don't require a normalized wave-function)

-> Start initially using a simpler method that converges linearly (bracketing):
-> Note: the bracketing method only works for points where the function crosses the x-axis (ie wont
work for x^2)


## Notes on Wave function
We have function of the form:

```
  Y'' = (1/c1) * (c2 - f(x)) * Y

  (or as Y'' = F(x) * Y )
```

where c1 = - h_bar^2 / (2m), c2 = Energy, f(x) = potential function

## Notes on energy scale
Energies are scaled using the following scale factor

J' = J * (4π * m / h),

where m is the mass of the particle in the potential
Note: While this makes the energies bigger, its not by too much.
-> Later experiment with different values for the scaling
-> Although: This certainly makes the equations more convenient