#include <chrono>
#include <iomanip>
#include <iostream>
#include <memory>
#include <vector>

#include "base.h"
#include "mod.h"
#include "ntt.h"
#include "ref.h"
#include "refr2.h"

using Clock = std::chrono::system_clock;
using MS = std::chrono::milliseconds;

double GetMiops(const NTT& ntt, std::vector<NTT::ElementType>& data);
bool Verify(const NTT& ntt);
int VerifyRoutines(const std::vector<NTTFactoryBase*>& factories);

int main(int argc, const char* argv[]) {
  const int64 kMaxLogN = 20;

  std::vector<NTTFactoryBase*> factories{
      new NTTFactory<RefNTT>,
      new NTTFactory<RefRadix2>,
  };

  if (VerifyRoutines(factories) || argc > 1) {
    return 0;
  }

  static constexpr int64 kColumnWidth = 10;
  std::cout << "| #    |";
  for (auto* factory : factories)
    std::cout << std::setw(kColumnWidth) << factory->name() << " |";
  std::cout << "\n"
            << "|------|";
  for (auto* factory : factories)
    std::cout << "----------:|";
  std::cout << "\n";

  if (argc > 1) {
    return 0;
  }

  // Measure performance
  for (int64 logn = 2; logn <= kMaxLogN; ++logn) {
    std::cout << "| 2^" << std::setw(2) << logn << " |";

    const int64 n = 1LL << logn;
    std::vector<NTT::ElementType> data(n);
    for (int64 i = 0; i < n; ++i)
      data[i] = uint64(i);
    for (auto* factory : factories) {
      std::unique_ptr<NTT> ntt(factory->Create(logn));
      double miops = GetMiops(*ntt, data);
      std::cout << std::setw(kColumnWidth) << std::fixed << std::setprecision(3)
                << miops << " |";
    }
    std::cout << "\n";
  }

  return 0;
}

double GetMiops(const NTT& ntt, std::vector<NTT::ElementType>& data) {
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
  double duration = std::chrono::duration_cast<MS>(end - start).count() * 1e-3;
  double average = duration / loop_count;
  double iop = ntt.getIop() * 2;
  return iop / average * 1e-6;
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
                << "Actual: " << xi << "\n"
                << "Expect: " << i << "\n";
      pass = false;
    }
  }
  return pass;
}
