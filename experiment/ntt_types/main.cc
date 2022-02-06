#include <chrono>
#include <iomanip>
#include <iostream>
#include <memory>
#include <vector>

#include "base.h"
#include "fft.h"  // PMP
#include "mod.h"
#include "ntt.h"
#include "pnt.h"
#include "ref.h"
#include "ref_no_table.h"
#include "refr2.h"

using Clock = std::chrono::system_clock;
using MS = std::chrono::milliseconds;

double GetFFTPerf(int64 logn_complex);
double GetPerf(const NTT& ntt, std::vector<NTT::ElementType>& data);
bool Verify(const NTT& ntt);
int VerifyRoutines(const std::vector<NTTFactoryBase*>& factories);
void MeasurePerformance(const std::vector<NTTFactoryBase*>& factories);

int main(int argc, const char* argv[]) {
  std::vector<NTTFactoryBase*> factories{
      new NTTFactory<RefNTT>, new NTTFactory<PNT>,
      // new NTTFactory<RefNoTable>,
      // new NTTFactory<RefRadix2>,
  };

  if (VerifyRoutines(factories) || argc > 1) {
    return 0;
  }
  MeasurePerformance(factories);

  return 0;
}

void MeasurePerformance(const std::vector<NTTFactoryBase*>& factories) {
  static constexpr int64 kColumnWidth = 10;
  static constexpr int64 kMaxLogN = 22;

  std::cout << "Time compared with FFT [NTT/FFT]. Smaller is better.\n";
  std::cout << "| #    |";
  for (auto* factory : factories)
    std::cout << std::setw(kColumnWidth) << factory->name() << " |";
  std::cout << "\n"
            << "|------|";
  for (auto* factory : factories)
    std::cout << "----------:|";
  std::cout << "\n";

  // Measure performance
  for (int64 logn = 2; logn <= kMaxLogN; ++logn) {
    std::cout << "| 2^" << std::setw(2) << logn << " |";
    const double fft_perf = GetFFTPerf(logn - 1);

    const int64 n = 1LL << logn;
    std::vector<NTT::ElementType> data(n);
    for (int64 i = 0; i < n; ++i)
      data[i] = uint64(i);
    for (auto* factory : factories) {
      std::unique_ptr<NTT> ntt(factory->Create(logn));
      double ntt_perf = GetPerf(*ntt, data);
      std::cout << std::setw(kColumnWidth) << std::fixed << std::setprecision(3)
                << (ntt_perf / fft_perf) << " |";
    }
    std::cout << "\n";
  }
}

double GetFFTPerf(int64 logn) {
  static constexpr int64 kTimeLimitSec = 1;
  static constexpr int64 kMaxLoopCount = 10000000;
  const int64 n = 1LL << logn;
  std::vector<Complex> data(n);
  for (int64 i = 0; i < n; ++i) {
    data[i] = Complex(2 * i, 2 * i + 1);
  }
  FFT fft(logn % 2, logn / 2);

  auto start = Clock::now();
  auto due_time = start + std::chrono::seconds(kTimeLimitSec);
  int64 loop_count = 0;
  for (loop_count = 0; loop_count < kMaxLoopCount && Clock::now() < due_time;
       ++loop_count) {
    fft.rft(data.data(), false);
    fft.rft(data.data(), true);
  }
  auto end = Clock::now();

  double duration_ms = std::chrono::duration_cast<MS>(end - start).count();
  double average = duration_ms / loop_count;
  double perf = average / (n * std::log2(n)) * 1e+6;
  return perf;
}

double GetPerf(const NTT& ntt, std::vector<NTT::ElementType>& data) {
  static constexpr int64 kTimeLimitSec = 1;
  static constexpr int64 kMaxLoopCount = 10000000;

  auto start = Clock::now();
  auto due_time = start + std::chrono::seconds(kTimeLimitSec);
  int64 loop_count = 0;
  for (loop_count = 0; loop_count < kMaxLoopCount && Clock::now() < due_time;
       ++loop_count) {
    ntt.ntt(data.data(), false);
    ntt.ntt(data.data(), true);
  }
  auto end = Clock::now();

  double duration_ms = std::chrono::duration_cast<MS>(end - start).count();
  int64 n = data.size();
  double average = duration_ms / loop_count;
  double perf = average / (n * std::log2(n)) * 1e+6;
  return perf;
}

int VerifyRoutines(const std::vector<NTTFactoryBase*>& factories) {
  static constexpr int64 kMaxCheckLogN = 10;
  for (int64 logn = 2; logn <= kMaxCheckLogN; ++logn) {
    const int64 n = 1LL << logn;
    for (auto* factory : factories) {
      std::unique_ptr<NTT> ntt(factory->Create(logn));
      if (!Verify(*ntt)) {
        std::cerr << "NTT " << factory->name() << " with n = " << n
                  << " has an issue.\n";
        return 1;
      }
    }
  }
  std::cerr << "Routines are verified.\n";
  return 0;
}

bool Verify(const NTT& ntt) {
  const int64 n = ntt.size();
  // For large N, we do not check it.
  if (n >= 256)
    return true;

  std::vector<NTT::ElementType> x(n);
  for (int64 i = 0; i < n; ++i) {
    x[i] = i;
  }
  ntt.ntt(x.data(), false);
  ntt.ntt(x.data(), true);
  bool pass = true;
  for (int64 i = 0; i < n; ++i) {
    auto& xi = x[i];
    if (xi != i) {
      std::cerr << "Fail at " << i << "/" << n << "\n"
                << "Actual: " << xi << " (0x" << std::hex << xi << std::dec
                << ")\n"
                << "Expect: " << i << "\n";
      pass = false;
    }
  }
  return pass;
}
