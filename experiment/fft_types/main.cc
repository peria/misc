#include <chrono>
#include <iomanip>
#include <iostream>
#include <memory>
#include <vector>

#include "stockham_dit.h"
#include "stockham_dif.h"
#include "cooley_ditf.h"
#include "cooley_dift.h"

using Clock = std::chrono::system_clock;
using MS = std::chrono::milliseconds;

int main() {
  std::vector<std::unique_ptr<FFT>> ffts;
  ffts.emplace_back(new StockhamDIT);
  ffts.emplace_back(new StockhamDIF);
  ffts.emplace_back(new CooleyDITF);
  ffts.emplace_back(new CooleyDIFT);

  std::cout << "     ";
  for (auto& fft : ffts) {
    std::cout << std::setw(10) << fft->name() << " ";
  }
  std::cout << "\n";

  for (int logn = 2; logn <= 23; ++logn) {
    std::cout << "2^" << std::setw(2) << logn << " ";

    const int n = 1 << logn;
    std::vector<Complex> data(n);
    for (auto& fft : ffts) {
      fft->setUp(n);
      for (int i = 0; i < n; ++i) {
        data[i].real = 2 * i;
        data[i].imag = 2 * i + 1;
      }
      auto start = Clock::now();
      for (int i = 0; i < 10; ++i) {
        fft->dft(data.data());
        fft->idft(data.data());
      }
      auto end = Clock::now();
      fft->tearDown();
      if (n < 256) {
        for (int i = 0; i < n; ++i) {
          if (std::abs(data[i].real - 2 * i) > 1e-4 ||
              std::abs(data[i].imag - (2 * i + 1)) > 1e-4) {
            std::cerr << fft->name() << " : " << i << "/" << n << "\n"
                      << "Actual: " << data[i].real << " + " << data[i].imag << "i\n"
                      << "Expect: " << 2 * i << " + " << 2 * i + 1 << "i\n";
            return 0;
          }
        }
      }
      auto dur = std::chrono::duration_cast<MS>(end - start).count() * 1e-3;
      std::cout << std::setw(10) << std::fixed << std::setprecision(3) << dur << " ";
    }
    std::cout << "\n";
  }
  return 0;
}
