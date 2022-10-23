#include "fmt.h"

#include <cmath>
#include <iostream>
#include <vector>

void FFT::Convolution(std::vector<Complex>& x, std::vector<Complex>& y) const {
  Rft(x.data());
  Rft(y.data());
  for (int i = 0; i < n_; ++i) {
    x[i] *= y[i];
  }
  IRft(x.data());
}
