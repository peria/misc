#include <chrono>
#include <iomanip>
#include <iostream>
#include <memory>
#include <vector>

#include "complex.h"
#include "cooley.h"
#include "dpmp.h"
#include "fft.h"
#include "pmp.h"
#include "pmp2.h"
#include "stockham_dif.h"
#include "stockham_dit.h"

// FFTs that do not work with w^(n/4) RFT.
#ifndef REAL_FFT
#include "pmp5.h"
#include "stockham6.h"
#endif

using Clock = std::chrono::system_clock;
using MS = std::chrono::milliseconds;

double MeasurePerformance(const FFT& fft, std::vector<Complex>& data);
bool Verify(const FFT& fft);

int main() {
  const int64 kMaxLogN = 23;

  std::vector<FFTFactoryBase*> factories {
    new FFTFactory<PMP>,
    new FFTFactory<DPMP>,
    new FFTFactory<PMP2>,
    new FFTFactory<Cooley>,
    new FFTFactory<StockhamDIT>,
    new FFTFactory<StockhamDIF>,
#ifndef REAL_FFT
    new FFTFactory<PMP5>,
    new FFTFactory<Stockham6>,
#endif
  };

  static constexpr int64 kColumnWidth = 10;
  std::cout << "| #    |";
  for (auto* factory : factories)
    std::cout << std::setw(kColumnWidth) << factory->name() << " |";
  std::cout << "\n"
            << "|------|";

  for (auto* factory : factories)
    std::cout << "----------:|";
  std::cout << "\n";

  for (int64 logn = 2; logn <= kMaxLogN; ++logn) {
    std::cout << "| 2^" << std::setw(2) << logn << " |";

    const int64 n = 1LL << logn;
    std::vector<Complex> data(n);
    for (int64 i = 0; i < n; ++i)
      data[i] = Complex(2 * i, 2 * i + 1);
    for (auto* factory : factories) {
      std::unique_ptr<FFT> fft(factory->Create(logn));
      if (!Verify(*fft))
        return 0;
      double mflops = MeasurePerformance(*fft, data);
      std::cout << std::setw(kColumnWidth) << std::fixed
                << std::setprecision(3) << mflops << " |";
    }
    std::cout << "\n";
  }

  return 0;
}

double MeasurePerformance(const FFT& fft, std::vector<Complex>& data) {
  static constexpr int64 kTimeLimitSec = 2;
  static constexpr int64 kMaxLoopCount = 10000000;

  auto start = Clock::now();
  auto due_time = start + std::chrono::seconds(kTimeLimitSec);
  int64 loop_count = 0;
  for (loop_count = 0; loop_count < kMaxLoopCount && Clock::now() < due_time; ++loop_count) {
#ifdef REAL_FFT
    fft.rft(data.data(), false);
    fft.rft(data.data(), true);
#else
    fft.dft(data.data(), false);
    fft.dft(data.data(), true);
#endif
  }
  auto end = Clock::now();
  double duration = std::chrono::duration_cast<MS>(end - start).count() * 1e-3;
  double average = duration / loop_count;
  double flop = fft.getFlop() * 2;
  return flop / average * 1e-6;
}

bool Verify(const FFT& fft) {
  static constexpr double kEPS = 1e-4;

  const int64 n = fft.size();
  // For large N, we do not check it.
  if (n >= 256)
    return true;

  std::vector<Complex> x(n);
  for (int64 i = 0; i < n; ++i) {
    x[i] = Complex(2 * i + 1, 2 * i + 2);
  }
#ifdef REAL_FFT
  fft.rft(x.data(), false);
  fft.rft(x.data(), true);
#else
  fft.dft(x.data(), false);
  fft.dft(x.data(), true);
#endif
  bool pass = true;
  for (int64 i = 0; i < n; ++i) {
    auto& d = x[i];
    const double er = 2 * i + 1;
    const double ei = 2 * i + 2;
    if (std::abs(d.real - er) > kEPS ||
        std::abs(d.imag - ei) > kEPS) {
      std::cerr << "Fail at " << i << "/" << n << "\n"
                << "Actual: " << d.real << " + " << d.imag << "i\n"
                << "Expect: " << er << " + " << ei << "i\n";
      pass = false;
    }
  }
  return pass;
}

