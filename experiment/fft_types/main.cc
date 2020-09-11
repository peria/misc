#include <chrono>
#include <iomanip>
#include <iostream>
#include <memory>
#include <vector>

//#define USE_RADIX8
#include "pmp.h"

using Clock = std::chrono::system_clock;
using MS = std::chrono::milliseconds;

std::vector<std::vector<Complex>> FFT::s_ws;

int main() {
  const int64 kMaxLogN = 26;
  FFT::resizeWs(kMaxLogN + 1);

  std::cout << "| #    | MFLOPS    |\n"
            << "|------|----------:|\n";
  for (int64 logn = 2; logn <= kMaxLogN; ++logn) {
    std::cout << "| 2^" << std::setw(2) << logn << " | ";

    const int64 n = 1LL << logn;
    std::vector<Complex> data(n);
    for (int64 i = 0; i < n; ++i) {
      data[i] = Complex(2 * i + 1, 2 * i + 2);
    }

    FFT fft(logn);
    auto start = Clock::now();
    auto due_time = start + std::chrono::seconds(3);
    int64 loop_count = 0;
    for (loop_count = 0; Clock::now() < due_time; ++loop_count) {
      fft.rft(data.data(), false);
      fft.rft(data.data(), true);
    }
    auto end = Clock::now();

    if (n < 256) {
      bool pass = true;
      for (int64 i = 0; i < n; ++i) {
        auto& d = data[i];
        if (std::abs(d.real - (2 * i + 1)) > 1e-4 ||
            std::abs(d.imag - (2 * i + 2)) > 1e-4) {
          std::cerr << i << "/" << n << "\n"
                    << "Actual: " << d.real << " + " << d.imag << "i\n"
                    << "Expect: " << (2 * i + 1) << " + " << (2 * i + 2)
                    << "i\n";
          pass = false;
        }
      }
      if (!pass)
        return 0;
    }
    auto duration = std::chrono::duration_cast<MS>(end - start).count() * 1e-3;
    std::cout << std::setw(9) << std::fixed << std::setprecision(3)
              << fft.getMFlops(duration / loop_count) << " |\n";
  }

  return 0;
}
