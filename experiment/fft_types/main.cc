#include <chrono>
#include <iomanip>
#include <iostream>
#include <memory>
#include <vector>

#include "complex.h"
#include "cooley.h"
#include "pmp.h"
#include "stockham_dit.h"
#include "stockham_dif.h"

using Clock = std::chrono::system_clock;
using MS = std::chrono::milliseconds;

double MeasurePerformance(const FFT& fft, std::vector<Complex>& data);
bool Verify(const std::vector<Complex>& data);

int main() {
  const int64 kMaxLogN = 23;

  std::vector<FFTFactoryBase*> factories {
    new FFTFactory<PMP>,
    new FFTFactory<Cooley>,
    new FFTFactory<StockhamDIT>,
    new FFTFactory<StockhamDIF>,
    // 6StepCT,
    // 6StepStockham
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
    for (int64 i = 0; i < n; ++i) {
      data[i] = Complex(2 * i + 1, 2 * i + 2);
    }

    for (auto* factory : factories) {
      std::unique_ptr<FFT> fft(factory->Create(logn));
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

  auto start = Clock::now();
  auto due_time = start + std::chrono::seconds(kTimeLimitSec);
  int64 loop_count = 0;
  for (loop_count = 0; loop_count < 1e+7 && Clock::now() < due_time; ++loop_count) {
    fft.dft(data.data(), false);
    fft.dft(data.data(), true);
  }
  auto end = Clock::now();
  double duration = std::chrono::duration_cast<MS>(end - start).count() * 1e-3;
  double average = duration / loop_count;
  double flop = fft.getFlop() * 2;
  return flop / average * 1e-6;
}

bool Verify(const std::vector<Complex>& data) {
  const int64 n = data.size();
  // Do not check
  if (n >= 256)
    return true;

  bool pass = true;
  static constexpr double kEPS = 1e-4;
  for (int64 i = 0; i < n; ++i) {
    auto& d = data[i];
    if (std::abs(d.real - (2 * i + 1)) > kEPS ||
        std::abs(d.imag - (2 * i + 2)) > kEPS) {
      std::cerr << i << "/" << n << "\n"
                << "Actual: " << d.real << " + " << d.imag << "i\n"
                << "Expect: " << (2 * i + 1) << " + " << (2 * i + 2)
                << "i\n";
      pass = false;
    }
  }
  return pass;
}

