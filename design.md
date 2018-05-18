# Algorithm design choices

-> Goal, given a potential function, we want to find the energy eigenvalues

## Numerov
- Use an Numerov step method to carry out the integration
- Use 'shooting' method to match boundary conditions (NR)
- Use 'Van Wijngaarden-Dekker-Brent' method to find root (i.e. boundary) while shooting (NR 9.3)

### Finding the energy (i.e finding roots)
Basically, the function will be

```
  F(E) - B = 0
```
where B is the value at the boundary and F is the integral computed using Numerov. (Note: The second
Psi/Y value is arbitrary as I don't require a normalized wave-function)

-> Start initially using a simpler method that converges linearly (bracketing):
-> Note: the bracketing method only works for points where the function crosses the x-axis (ie won't
work for x^2)


## Notes on Wave function
We have function of the form:

```
  Y'' = (1/c1) * (c2 - f(x)) * Y

  (or as Y'' = F(x) * Y )
```

where c1 = - h_bar^2 / (2m), c2 = Energy, f(x) = potential function

## Notes on energy scale
Energies and lengths are scaled using the following equation

1 = c1,

where we scale the length unit (meters) to achieve this nice form.
Note: This allows us to rewrite the ODE as `Y'' = (f(x)-E) * Y`
Note: For systems substantially bigger than atoms, this scaling makes
      energies to small to be detected

## Energy values for inf square well to small
Multiplying by `2.4` gives pretty exact energies for *ALL* energies computed,
which probably means, that we're missing a constant factor of about that somewhere.
-> Now: 1 order of magnitude to small, but 2.4 factor gives correct digits

## Coming from both sides
-> Check sign of integration
-> Note, use identity of derivatives

## Note on normalization
-> Using adaptive normalization to a max of 10.0, the
   integration is (at least number-range wise) very resilient
   -> Use this to build integrator for symmetric problems using
      gradient condition!