#include "fmt.h"

#include <cmath>
#include <iostream>
#include <vector>
#include <chrono>

#include "fft.h"
#include "ntt.h"
using namespace std;

namespace {

using Clock = chrono::system_clock;

// Settings
const int kMinPower = 10;
const int kMaxPower = 20;
const int kTryTimes = 10;
const int kMaxN = 1 << kMaxPower;

unique_ptr<uint64[]> g_data;
unique_ptr<uint64[]> g_work;
unique_ptr<uint64[]> g_table;
unique_ptr<uint64[]> g_inv_table;

double ToMs(chrono::time_point<Clock>& end,
            chrono::time_point<Clock>& start) {
  auto diff = end - start;
  return chrono::duration_cast<chrono::microseconds>(diff).count() * 1e-3;
}

}  // namespace

int main() {
  g_data.reset(new uint64[kMaxN]);
  g_work.reset(new uint64[kMaxN]);
  g_table.reset(new uint64[kMaxN]);
  g_inv_table.reset(new uint64[kMaxN]);

  double* double_data = reinterpret_cast<double*>(g_data.get());
  uint64* uint64_data = g_data.get();

  RFT::SetGlobals(g_work.get(), g_table.get(), g_inv_table.get());
  NTT::SetGlobals(g_work.get(), g_table.get(), g_inv_table.get());

  if (!RFT::Validate(double_data) || !NTT::Validate(uint64_data))
    return 0;

  // with table initialization
  cout << "Initialize the table everytime in FFT/NTT's call.\n";
  for (int log2n = kMinPower; log2n <= kMaxPower; ++log2n) {
    const int n = 1 << log2n;
    auto fft_start = Clock::now();
    for (int count = 0; count < kTryTimes; ++count) {
      RFT::InitTable(log2n);
      RFT::Forward(log2n, n, double_data);
      RFT::Forward(log2n, n, double_data);
      RFT::Square(double_data, n);
      RFT::Backward(log2n, n, double_data);
    }
    auto fft_end = Clock::now();

    auto ntt_start = Clock::now();
    for (int count = 0; count < kTryTimes; ++count) {
      NTT::InitTable(log2n);
      NTT::Forward(log2n, n, uint64_data);
      NTT::Forward(log2n, n, uint64_data);
      NTT::Square(uint64_data, n);
      NTT::Backward(log2n, n, uint64_data);
    }
    auto ntt_end = Clock::now();
    cout << "2^" << log2n << " "
         << ToMs(fft_end, fft_start) << " "
         << ToMs(ntt_end, ntt_start) << "\n";
  }

  // without table initialization
  cout << "Initialize the table once.\n";
  for (int log2n = kMinPower; log2n <= kMaxPower; ++log2n) {
    const int n = 1 << log2n;
    RFT::InitTable(log2n);
    auto fft_start = Clock::now();
    for (int count = 0; count < kTryTimes; ++count) {
      RFT::Forward(log2n, n, double_data);
      RFT::Forward(log2n, n, double_data);
      RFT::Square(double_data, n);
      RFT::Backward(log2n, n, double_data);
    }
    auto fft_end = Clock::now();

    NTT::InitTable(log2n);
    auto ntt_start = Clock::now();
    for (int count = 0; count < kTryTimes; ++count) {
      NTT::Forward(log2n, n, uint64_data);
      NTT::Forward(log2n, n, uint64_data);
      NTT::Square(uint64_data, n);
      NTT::Backward(log2n, n, uint64_data);
    }
    auto ntt_end = Clock::now();
    cout << "2^" << log2n << " "
         << ToMs(fft_end, fft_start) << " "
         << ToMs(ntt_end, ntt_start) << "\n";
  }

  return 0;
}
