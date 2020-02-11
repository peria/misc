#include <chrono>
#include <cstdio>
#include <memory>
#include <string>

#include "computer.h"

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

int main(int argc, char* argv[]) {
  int64_t digits = (argc > 1) ? std::atoi(argv[1]) : kDefaultDigits;
  if (digits <= 0)
    digits = kDefaultDigits;
  char* answer_file = (argc > 2) ? argv[2] : nullptr;

  std::unique_ptr<Computer> computer(new Chudnovsky(digits));
  // std::unique_ptr<Computer> computer(new Ramanujan(digits));
  {
    Timer all("All");
    {
      Timer timer("Compute");
      computer->compute();
    }
    computer->output();
  }
  computer->check(answer_file);

  return 0;
}
