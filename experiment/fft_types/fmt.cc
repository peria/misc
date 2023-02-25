#include "fmt.h"

#include <chrono>
#include <cmath>
#include <iostream>
#include <vector>

#define CHECK_DOUBLE_EQ(x, y)                               \
  {                                                         \
    double vx = (x);                                        \
    double vy = (y);                                        \
    if (std::abs((vx) - (vy)) > 1e-6) {                     \
      std::cerr << "Different\n"                            \
                << "Actual: " << #x << " = " << vx << "\n"  \
                << "Expect: " << #y << " = " << vy << "\n"; \
      return false;                                         \
    }                                                       \
  }

void FMT::Convolute(Complex* x, Complex* y) const {
  Rft(x);
  Rft(y);
  for (int i = 0; i < n_; ++i) {
    x[i] *= y[i];
  }
  IRft(x);
}

bool FMT::Test() const {
  std::vector<Complex> x(n_), y(n_);
  for (int i = 0; i < n_; ++i) {
    x[i].real = 1;
    y[i].real = 1;
  }
  Convolute(x.data(), y.data());
  for (int i = 0; i < n_; ++i) {
    CHECK_DOUBLE_EQ(x[i].real, i + 1);
  }
  return true;
}

double FMTFactoryBase::MeasureConvPerf(int logn) const {
  using Clock = std::chrono::steady_clock;
  using MS = std::chrono::milliseconds;

  static constexpr int kMaxLoop = 100000;
  static constexpr int kTimeLimitSec = 1;

  auto&& fmt = Create(logn);
  const int n = 1 << logn;
  std::vector<Complex> x(n), y(n);

  auto start = Clock::now();
  auto due = start + std::chrono::seconds(kTimeLimitSec);
  int count = 0;
  for (count = 0; count < kMaxLoop && Clock::now() < due; ++count) {
    fmt->Convolute(x.data(), y.data());
    x = y;
  }
  auto finish = Clock::now();

  auto duration = finish - start;
  auto duration_sec = std::chrono::duration_cast<MS>(duration).count() * 1e-3;
  return fmt->GetFlops() * 3 / duration_sec * count * 1e-6;
}
