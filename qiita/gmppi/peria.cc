#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <memory>
#include <string>

#include <gmpxx.h>

using Integer = mpz_class;
using Float = mpf_class;

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

void drm(const int64_t n0,
         const int64_t n1,
         Integer& x0,
         Integer& y0,
         Integer& z0) {
  if (n0 + 1 == n1) {
    static constexpr int64_t A = 13591409;
    static constexpr int64_t B = 545140134;
    static constexpr int64_t C = 640320;

    if (n0 == 0) {
      x0 = 1;
    } else {
      x0 = n0 * C;
      x0 *= n0 * C;
      x0 *= n0 * C / 24;
    }
    y0 = A + B * n0;
    z0 = -(6 * n0 + 5);
    z0 *= 6 * n0 + 1;
    z0 *= 2 * n0 + 1;
    return;
  }

  Integer x1, y1, z1;
  int64_t m = (n0 + n1) / 2;
  drm(n0, m, x0, y0, z0);
  drm(m, n1, x1, y1, z1);

  // y0 = x1 * y0 + y1 * z0;
  y0 *= x1;
  y1 *= z0;
  y0 += y1;

  x0 *= x1;
  z0 *= z1;
}

void compute(const int64_t digits, Float& pi) {
  const int64_t precision = digits * std::log2(10);
  mpf_set_default_prec(precision);
  pi.set_prec(precision);

  const int64_t n = digits / 14;
  Integer x, y, z;
  drm(0, n, x, y, z);

  x *= 426880;
  mpf_sqrt_ui(pi.get_mpf_t(), 10005);
  pi *= x;
  pi /= y;
}

void output(Float& pi) {
  std::ofstream ofs("pi_peria.out");
  if (!ofs.is_open()) {
    std::cerr << "Could not open the file\n";
    return;
  }
  int64_t precision = pi.get_prec();
  int64_t digits = precision / std::log2(10) + 1;
  ofs << std::setprecision(digits) << pi;
}

int main(int argc, char* argv[]) {
  int64_t digits = (argc > 1) ? std::atoi(argv[1]) : kDefaultDigits;
  if (digits <= 0)
    digits = kDefaultDigits;

  Timer all("All");
  Float pi;
  {
    Timer timer("Compute");
    compute(digits, pi);
  }
  output(pi);

  return 0;
}
