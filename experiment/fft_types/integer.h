#pragma once

#include <cstdint>
#include <memory>
#include <vector>

#include "fmt.h"

using int64 = std::int64_t;

class Integer {
 public:
  using Digit = int64;
  static void Convolution(const FMT& fmt,
                          const int n,
                          std::vector<Digit>& x,
                          const std::vector<Digit>& y);
  static double GetFlops(int logn);
  static void Split(std::vector<Complex>& fx, const std::vector<Digit>& x);
  static void Merge(std::vector<Digit>& x, const std::vector<Complex>& fx);

  // Returns MFlops for FFTs, and pseudo MFloops for NTTs.
  static int MeasurePerformance(const FMTFactory& factory, int logn);
};
