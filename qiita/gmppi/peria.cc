#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <memory>
#include <string>
#include <vector>

#include <gmpxx.h>

constexpr int kDefaultDigits = 10000;

class Timer {
  using Clock = std::chrono::system_clock;
  using Ms = std::chrono::milliseconds;

 public:
  Timer(const std::string& name) : name_(name), start_(Clock::now()) {}
  ~Timer() {
    auto duration = Clock::now() - start_;
    std::printf("%7s: %.3f\n", name_.c_str(),
                std::chrono::duration_cast<Ms>(duration).count() * 1e-3);
  }

 private:
  std::string name_;
  Clock::time_point start_;
};

class Formula {
 public:
  enum Type {
    kChudnovsky,
  };
};

void setupFormula(Formula::Type) {}
void compute(const int digits, mpf_class& pi) {}
void output(mpf_class& pi) {}

int main(int argc, char* argv[]) {
  int digits = (argc > 2) ? std::atoi(argv[1]) : kDefaultDigits;
  if (digits <= 0)
    digits = kDefaultDigits;

  setupFormula(Formula::kChudnovsky);
  Timer all("All");
  mpf_class pi;
  {
    Timer timer("Compute");
    compute(digits, pi);
  }
  output(pi);

  return 0;
}
