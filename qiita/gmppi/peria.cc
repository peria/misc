#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <memory>
#include <optional>
#include <string>

#include "computer.h"

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

Computer::Configuration parseArgs(int argc, char** argv) {
  Computer::Configuration config;
  for (int i = 1, c = 0; i < argc; ++i) {
    char* arg = argv[i];
    if (arg[0] != '-') {
      switch (c++) {
        case 0:
          config.digits = std::strtoll(arg, nullptr, 10);
          break;
        case 1:
          config.compare_file = arg;
          break;
      }
      continue;
    }

    char* key = arg;
    while (*key == '-')
      ++key;
    char* value = (i + 1 < argc) ? argv[++i] : nullptr;
    if (std::strcmp(key, "ramanujan") == 0) {
      config.formula = Computer::Formula::kRamanujan;
    } else if (std::strcmp(key, "no-output") == 0 ||
               std::strcmp(key, "no_output") == 0) {
      config.outputs = false;
    } else if (std::strcmp(key, "threads") == 0) {
      config.number_of_threads = std::atoi(value);
    }
  }

  return config;
}

int main(int argc, char* argv[]) {
  Computer::Configuration config = parseArgs(argc, argv);

  std::unique_ptr<Computer> computer = Computer::Create(config);
  {
    Timer all("All");
    {
      Timer timer("Compute");
      computer->compute();
    }
    computer->output();
  }
  computer->check();

  return 0;
}
