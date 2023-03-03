#include <iomanip>
#include <iostream>
#include <memory>
#include <vector>

#include "dif.h"
#include "dit.h"
#include "fmt.h"

using FactoryVec = std::vector<std::shared_ptr<FMTFactoryBase>>;

bool TestFMTs(const FactoryVec& factories);
void MeasurePerformance(const FactoryVec& factories);

namespace {
constexpr int kMaxLogN = 23;
}

int main() {
  FactoryVec factories{
      FMTFactory<DIT>::GetFactory(),
      FMTFactory<DIF>::GetFactory(),
  };

  if (!TestFMTs(factories)) {
    return 0;
  }

  MeasurePerformance(factories);

  return 0;
}

bool TestFMTs(const FactoryVec& factories) {
  for (auto&& factory : factories) {
    for (int logn = 2; logn <= kMaxLogN; ++logn) {
      auto&& fmt = factory->Create(logn);
      if (!fmt->Test()) {
        std::cerr << factory->name()
                  << " has some problems. (n = " << (1 << logn) << ")\n";
        return false;
      }
    }
  }
  std::cerr << "Tests passed\n";
  return true;
}

void MeasurePerformance(const FactoryVec& factories) {
  // Header (Name)
  std::cerr << "    ";
  for (auto&& factory : factories) {
    std::cerr << " | " << std::setw(4) << factory->name();
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
      int mflops = factory->MeasureConvPerf(logn);
      std::cerr << " | " << std::setw(4) << mflops;
    }
    std::cerr << std::endl;
  }
}
