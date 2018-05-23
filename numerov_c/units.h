#pragma once
#import <math.h>

#define UNITS_H 6.626070040e-34l   // Planck's Constant according to NIST
#define UNITS_PI 3.14159265358979l // Circle constant
#define UNITS_U 1.660539040e-27l   // Atomic mass unit
#define UNITS_E 1.6021766208e-19l  // Fundamental charge

// Define macros for scaling / unscaling units or length (energy)
#define units_gammaSquared(m) sqrt(8l * UNITS_PI*UNITS_PI * ((double)m) / UNITS_H / UNITS_H)
#define units_gamma(m) sqrt(units_gammaSquared(m))

#define units_scaleE(E, m) (((double)E) * units_gammaSquared(m))
#define units_unscaleE(E, m) (((double)E) / units_gammaSquared(m))

#define units_scaleL(L, m) (((double)L) * units_gamma(m))
#define units_unscaleL(L, m) (((double)L) / units_gamma(m))