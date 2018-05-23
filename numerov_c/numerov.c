#include "numerov.h"


// Numerov c implementation
double numerov_step(
  double y_k,          // value at k
  double y_kk,         // value at k - 1
  double f_k,          // F(x_k) (Y'' = f(x)*Y)
  double f_kk,         // F(x_k-1)
  double f,            // F(x_k+1)
  double step_width
) {
  double pre = 1l - step_width*step_width * f / 12l;
  double res = y_k * (step_width*step_width * f_k * 5l / 6l + 2l);
  res += y_kk * step_width*step_width * f_kk / 12l - y_kk;
  return res / pre;
}

// TODO: Offer integration (?)