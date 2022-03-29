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

using Clock = std::chrono::system_clock;
using MS = std::chrono::milliseconds;

static constexpr int64 kMaxLogN = 23;
static constexpr int64 kMaxCheckLogN = 10;

enum class VerifyResult {
  kCorrect,
  kBugInDft,
  kBugInRft,
};

void MeasureFFTs(const std::vector<FFTFactoryBase*>& factories);
double MeasureFFT(const FFT& fft);
bool VerifyFFTs(const std::vector<FFTFactoryBase*>& factories, bool is_check);
VerifyResult VerifyFFT(const FFT& fft);

int main(int argc, const char**) {
  std::vector<FFTFactoryBase*> factories{
      new FFTFactory<PMP>,         new FFTFactory<DPMP>,
      new FFTFactory<PMP2>,        new FFTFactory<Cooley>,
      new FFTFactory<StockhamDIT>, new FFTFactory<StockhamDIF>,
  };

  if (!VerifyFFTs(factories, argc > 1)) {
    return 0;
  }

  MeasureFFTs(factories);
  return 0;
}

void MeasureFFTs(const std::vector<FFTFactoryBase*>& factories) {
  static constexpr int64 kColumnWidth = 10;
  std::cout << "| #    |";
  for (auto&& factory : factories)
    std::cout << std::setw(kColumnWidth) << factory->name() << " |";
  std::cout << "\n";

  std::cout << "|------|";
  for (auto&& factory : factories)
    std::cout << "----------:|";
  std::cout << "\n";

  for (int64 logn = 2; logn <= kMaxLogN; ++logn) {
    std::cout << "| 2^" << std::setw(2) << logn << " |";
    for (auto* factory : factories) {
      std::unique_ptr<FFT> fft(factory->Create(logn));
      double mflops = MeasureFFT(*fft);
      std::cout << std::setw(kColumnWidth) << std::fixed << std::setprecision(3)
                << mflops << " |";
    }
    std::cout << "\n";
  }
}

double MeasureFFT(const FFT& fft) {
  static constexpr int64 kTimeLimitSec = 2;
  static constexpr int64 kMaxLoopCount = 10000000;
  const int64 n = fft.size();

  std::vector<Complex> data(n);
  for (int64 i = 0; i < n; ++i)
    data[i] = Complex(2 * i, 2 * i + 1);

  auto start = Clock::now();
  auto due_time = start + std::chrono::seconds(kTimeLimitSec);
  int64 loop_count = 0;
  while (loop_count < kMaxLoopCount && Clock::now() < due_time) {
    fft.rft(data.data(), false);
    fft.rft(data.data(), true);
    ++loop_count;
  }
  auto end = Clock::now();

  double duration = std::chrono::duration_cast<MS>(end - start).count() * 1e-3;
  double average = duration / loop_count;
  double flop = fft.getFlop() * 2;
  return flop / average * 1e-6;
}

bool VerifyFFTs(const std::vector<FFTFactoryBase*>& factories, bool is_check) {
  for (auto* factory : factories) {
    if (is_check) {
      std::cerr << "Checking " << factory->name() << " ... ";
    }
    for (int64 logn = 2; logn <= kMaxCheckLogN; ++logn) {
      const int64 n = 1LL << logn;
      std::vector<Complex> data(n);
      for (int64 i = 0; i < n; ++i)
        data[i] = Complex(2 * i, 2 * i + 1);
      std::unique_ptr<FFT> fft(factory->Create(logn));
      switch (VerifyFFT(*fft)) {
        case VerifyResult::kCorrect:
          break;
        case VerifyResult::kBugInDft:
          std::cerr << factory->name() << " has a bug in DFT.\n";
          return false;
        case VerifyResult::kBugInRft:
          std::cerr << factory->name() << " has a bug in RFT.\n";
          return false;
      }
    }
    if (is_check) {
      std::cerr << "OK\n";
    }
  }
  return !is_check;
}

VerifyResult VerifyFFT(const FFT& fft) {
  static constexpr double kEPS = 1e-4;

  const int64 n = fft.size();
  std::vector<Complex> x(n);
  for (int64 i = 0; i < n; ++i) {
    x[i] = Complex(2 * i + 1, 2 * i + 2);
  }
  fft.dft(x.data(), false);
  fft.dft(x.data(), true);

  VerifyResult result = VerifyResult::kCorrect;
  for (int64 i = 0; i < n; ++i) {
    auto& d = x[i];
    const double er = 2 * i + 1;
    const double ei = 2 * i + 2;
    if (std::abs(d.real - er) > kEPS || std::abs(d.imag - ei) > kEPS) {
      std::cerr << "Fail at " << i << "/" << n << "\n"
                << "Actual: " << d.real << " + " << d.imag << "i\n"
                << "Expect: " << er << " + " << ei << "i\n";
      result = VerifyResult::kBugInDft;
    }
  }
  if (result != VerifyResult::kCorrect) {
    return result;
  }

  fft.rft(x.data(), false);
  fft.rft(x.data(), true);

  for (int64 i = 0; i < n; ++i) {
    auto& d = x[i];
    const double er = 2 * i + 1;
    const double ei = 2 * i + 2;
    if (std::abs(d.real - er) > kEPS || std::abs(d.imag - ei) > kEPS) {
      std::cerr << "Fail at " << i << "/" << n << "\n"
                << "Actual: " << d.real << " + " << d.imag << "i\n"
                << "Expect: " << er << " + " << ei << "i\n";
      result = VerifyResult::kBugInRft;
    }
  }

  return result;
}
