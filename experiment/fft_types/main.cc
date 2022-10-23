#include <chrono>
#include <iomanip>
#include <iostream>
#include <memory>
#include <vector>

#include "fmt.h"
#include "integer.h"

#include "dpmp.h"

using FactoryVec = std::vector<std::shared_ptr<FMTFactory>>;

bool Test(const FactoryVec& factories);
void MeasurePerformance(const FactoryVec& factories);

int main() {
  FactoryVec factories{
      std::make_unique<DPMPFactory>(),
  };
  if (!Test(factories)) {
    return 0;
  }
  std::cerr << "Tests passed\n";
  MeasurePerformance(factories);

  return 0;
}

bool Test(const FMTFactory& factory) {
  return true;
}

bool Test(const FactoryVec& factories) {
  for (auto&& factory : factories) {
    if (!Test(*factory)) {
      std::cerr << factory->name() << " has some problems.\n";
      return false;
    }
  }
  return true;
}

void MeasurePerformance(const FactoryVec& factories) {
  static constexpr int kMaxLogN = 21;

  // Header (Name)
  std::cerr << "    ";
  for (auto&& factory : factories) {
    std::cerr << " | " << factory->name();
  }
  std::cerr << std::endl;

  std::cerr << "-----";
  for (auto&& factory : factories) {
    std::cerr << "+------";
  }
  std::cerr << std::endl;

  for (int logn = 3; logn <= kMaxLogN; ++logn) {
    std::cerr << "2^" << std::setw(2) << logn;
    for (auto&& factory : factories) {
      int mflops = 3 * Integer::MeasurePerformance(*factory, logn);
      std::cerr << " | " << std::setw(4) << mflops;
    }
    std::cerr << std::endl;
  }
}

int Integer::MeasurePerformance(const FMTFactory& factory, int logn) {
  static constexpr int kMaxLoop = 100000;
  static constexpr int kTimeLimitSec = 1;
  using Clock = std::chrono::steady_clock;
  using MS = std::chrono::milliseconds;

  auto&& fmt = factory.Create(logn);
  int n = 1 << logn;
  std::vector<Digit> x(n, 123), y(n, 123), z(n);

  auto start = Clock::now();
  auto due = start + std::chrono::seconds(kTimeLimitSec);
  int count = 0;
  for (count = 0; count < kMaxLoop && Clock::now() < due; ++count) {
    Convolution(*fmt, n, x, y);
    x = y;
  }
  auto finish = Clock::now();

  auto duration = finish - start;
  auto duration_sec = std::chrono::duration_cast<MS>(duration).count() * 1e-3;
  return Integer::GetFlops(logn) / duration_sec * count * 1e-6;
}
