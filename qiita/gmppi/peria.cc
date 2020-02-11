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

class Computer {
 public:
  Computer(int64_t digits) : digits_(digits) {}
  virtual ~Computer() = default;

  void compute();
  void output();

 protected:
  mpf_class pi_;

 private:
  void drm(const int64_t n0,
           const int64_t n1,
           mpz_class& x0,
           mpz_class& y0,
           mpz_class& z0);
  virtual void setXYZ(int64_t k, mpz_class& x, mpz_class& y, mpz_class& z) = 0;
  virtual void postProcess(mpz_class& x, mpz_class& y) = 0;
  virtual int64_t terms(int64_t digits) const = 0;

  const int64_t digits_;
};

class Chudnovsky : public Computer {
 public:
  Chudnovsky(int64_t digits) : Computer(digits) {}
  ~Chudnovsky() override = default;

 private:
  void setXYZ(int64_t k, mpz_class& x, mpz_class& y, mpz_class& z) override;
  void postProcess(mpz_class& x, mpz_class& y) override;
  int64_t terms(int64_t digits) const override { return digits / 14; }
};

void Chudnovsky::setXYZ(int64_t k, mpz_class& x, mpz_class& y, mpz_class& z) {
  static constexpr int64_t A = 13591409;
  static constexpr int64_t B = 545140134;
  static constexpr int64_t C = 640320;

  if (k == 0) {
    x = 1;
  } else {
    x = k * C;
    x *= k * C;
    x *= k * C / 24;
  }
  y = A + B * k;
  z = -(6 * k + 5);
  z *= 6 * k + 1;
  z *= 2 * k + 1;
}

void Chudnovsky::postProcess(mpz_class& x, mpz_class& y) {
  x *= 426880;
  mpf_sqrt_ui(pi_.get_mpf_t(), 10005);
  pi_ *= x;
  pi_ /= y;
}

class Ramanujan : public Computer {
 public:
  Ramanujan(int64_t digits) : Computer(digits) {}
  ~Ramanujan() override = default;

 private:
  void setXYZ(int64_t k, mpz_class& x, mpz_class& y, mpz_class& z) override;
  void postProcess(mpz_class& x, mpz_class& y) override;
  int64_t terms(int64_t digits) const override { return digits / 7.98; }
};

void Ramanujan::setXYZ(int64_t k, mpz_class& x, mpz_class& y, mpz_class& z) {
  static constexpr int64_t A = 1103;
  static constexpr int64_t B = 26390;
  static constexpr int64_t C = 396;

  if (k == 0) {
    x = 1;
  } else {
    x = k * C / 2;
    x *= k * C / 2;
    x *= k * C / 2;
    x *= k * C;
  }
  y = A + B * k;
  z = 4 * k + 1;
  z *= 2 * k + 1;
  z *= 4 * k + 3;
  z *= k + 1;
}

void Ramanujan::postProcess(mpz_class& x, mpz_class& y) {
  x *= 9801;
  y *= 4;
  mpf_sqrt_ui(pi_.get_mpf_t(), 2);
  pi_ *= x;
  pi_ /= y;
}

void Computer::drm(const int64_t n0,
                   const int64_t n1,
                   mpz_class& x0,
                   mpz_class& y0,
                   mpz_class& z0) {
  if (n0 + 1 == n1) {
    setXYZ(n0, x0, y0, z0);
    return;
  }

  mpz_class x1, y1, z1;
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

void Computer::compute() {
  const int64_t precision = digits_ * std::log2(10) + 10;
  pi_.set_prec(precision);

  const int64_t n = terms(digits_);
  mpz_class x, y, z;
  drm(0, n, x, y, z);
  postProcess(x, y);
}

void Computer::output() {
  std::ofstream ofs("pi_peria.out");
  if (!ofs.is_open()) {
    std::cerr << "Could not open the file\n";
    return;
  }
  ofs << std::setprecision(digits_ + 2) << pi_;
}

int main(int argc, char* argv[]) {
  int64_t digits = (argc > 1) ? std::atoi(argv[1]) : kDefaultDigits;
  if (digits <= 0)
    digits = kDefaultDigits;

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

  return 0;
}
