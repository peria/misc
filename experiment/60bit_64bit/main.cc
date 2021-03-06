#include <chrono>
#include <cmath>
#include <cstring>
#include <iostream>
#include <vector>

#include "base.h"
#include "fft.h"

using Clock = std::chrono::system_clock;
using Ms = std::chrono::milliseconds;

enum class NormType {
  Positive,  // [0, B)
  Negative,  // [-B/2, B/2)
};

void DataGenerate(const int64 bits,
                  const int64 n,
                  const NormType norm,
                  std::vector<double>& data) {
  const uint64 value = (norm == NormType::Positive) ?
    ((1ULL << bits) - 1) : ((1ULL << (bits - 1)) - 1);
  for (int64 i = 0; i < n / 2; ++i) {
    data[i] = value;
  }
  if (norm == NormType::Positive)
    return;

  const double kBase = 1ULL << bits;
  for (int64 i = 0; i < n / 2; ++i) {
    if (data[i] >= kBase / 2) {
      data[i] -= kBase;
      data[i + 1] -= 1;
    }
  }
}

double Normalize(const int64 bits, const int64 n, std::vector<double>& c) {
  const double kBase = 1ULL << bits;

  double max_error = 0;
  double carry = 0;
  for (int64 i = 0; i < n; ++i) {
    double ic = std::floor(c[i] + 0.5);
    double e = std::abs(ic - c[i]);
    max_error = std::max(max_error, e);
    carry = std::floor((ic + carry) / kBase);
    c[i] = ic - carry * kBase;
  }
  return max_error;
}

void Convolution(const Rft& rft,
                 std::vector<double>& a,
                 std::vector<double>& b,
                 std::vector<double>& c) {
  rft.run(Fft::Direction::Forward, a);
  rft.run(Fft::Direction::Forward, b);
  c[0] = a[0] * b[0];
  c[1] = a[1] * b[1];
  for (int64 i = 2; i < static_cast<int64>(a.size()); i += 2) {
    c[i] = a[i] * b[i] - a[i + 1] * b[i + 1];
    c[i + 1] = a[i] * b[i + 1] + a[i + 1] * b[i];
  }
  rft.run(Fft::Direction::Backward, c);
}

bool Experiment(const int64 bits,
                const int64 n,
                const int64 log2n,
                const Fft::Radix use,
                const NormType norm) {
  std::vector<double> a(n, 0), b(n, 0), c(n, 0);  // data
  std::vector<double> wk(n), table(n);            // for work
  DataGenerate(bits, n, norm, a);
  DataGenerate(bits, n, norm, b);
  auto start = Clock::now();
  Rft rft(n, log2n, use, wk, table);
  Convolution(rft, a, b, c);
  double max_error = Normalize(bits, n, c);
  auto end = Clock::now();

  if (max_error >= .1)
    return false;

  // total bytes, time[ms], maximum error, Radix, normalize type
  std::cout << (bits * n / 8) << " "
            << std::chrono::duration_cast<Ms>(end - start).count() << " "
            << bits << " " << use << " "
            << (norm == NormType::Positive ? "+" : "-") << " " << max_error
            << std::endl;
  return true;
}

void Test() {
  struct {
    int64 base;
    Fft::Radix radix;
  } data_set[] = {
      {1, Fft::Radix::Two}, {3, Fft::Radix::Three}, {5, Fft::Radix::Five},
  };

  for (auto data : data_set) {
    for (int64 log2n = 1; log2n < 18; ++log2n) {
      const int64 n = data.base << log2n;
      std::vector<double> a(n, 0);
      std::vector<double> wk(n), table(n);  // for work
      for (int64 i = 0; i < n; ++i)
        a[i] = i + 1;
      Rft rft(n, log2n, data.radix, wk, table);
      rft.run(Fft::Direction::Forward, a);
      rft.run(Fft::Direction::Backward, a);
      bool pass = true;
      for (int64 i = 0; i < n; ++i) {
        double diff = std::abs(a[i] - (i + 1));
        if (diff > 1e-3)
          pass = false;
      }
      if (!pass) {
        std::cerr << "Error in " << n << " = "
                  << (data.radix == Fft::Radix::Two
                          ? ""
                          : data.radix == Fft::Radix::Three ? "3 * " : "5 * ")
                  << "2^" << log2n << " in real, "
                  << (data.radix == Fft::Radix::Two
                          ? ""
                          : data.radix == Fft::Radix::Three ? "3 * " : "5 * ")
                  << "2^" << log2n - 1 << " in complex\n";
      }
    }
  }
}

int main(int argc, char* argv[]) {
  if (argc > 1 && std::strcmp(argv[1], "test") == 0) {
    Test();
    return 0;
  }

  struct {
    const int bits;
    const int64 base;
    const Fft::Radix radix;
    const NormType norm;
  } kParameters[] = {
      {16, 1, Fft::Radix::Two, NormType::Positive},
      {15, 1, Fft::Radix::Two, NormType::Positive},
      {20, 3, Fft::Radix::Three, NormType::Positive},
      {12, 5, Fft::Radix::Five, NormType::Positive},
      {16, 1, Fft::Radix::Two, NormType::Negative},
      {15, 1, Fft::Radix::Two, NormType::Negative},
      {20, 3, Fft::Radix::Three, NormType::Negative},
      {12, 5, Fft::Radix::Five, NormType::Negative},
  };

  static constexpr int64 kMaxTotalMemoryMB = 1024;
  for (auto& param : kParameters) {
    for (int64 log2n = 2;; ++log2n) {
      const int64 n = param.base << log2n;
      if ((n >> 20) * 4 > kMaxTotalMemoryMB)
        break;
      if (!Experiment(param.bits, n, log2n, param.radix, param.norm))
        break;
    }
    std::cout << std::endl << std::endl;
  }

  return 0;
}
