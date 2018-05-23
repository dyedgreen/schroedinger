#pragma once

// Numerov step function
double numerov_step(
  double y_k,          // value at k
  double y_kk,         // value at k - 1
  double f_k,          // F(x_k) (Y'' = f(x)*Y)
  double f_kk,         // F(x_k-1)
  double f,            // F(x_k+1)
  double step_width
);