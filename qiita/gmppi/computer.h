#pragma once

#include <cstdint>
#include <memory>
#include <optional>
#include <string>

#include <gmpxx.h>

class Computer {
 public:
  static constexpr int kDefaultDigits = 10000;
  enum class Formula { kChudnovsky, kRamanujan };
  struct Configuration {
    int64_t digits = kDefaultDigits;
    Formula formula = Formula::kChudnovsky;
    bool outputs = true;
    std::optional<std::string> compare_file;
    size_t number_of_threads = 1;
  };

  static std::unique_ptr<Computer> Create(const Configuration& config);
  virtual ~Computer() = default;

  void compute();
  void output();
  void check();

 protected:
  Computer(const Configuration& config) : config_(config) {}

  mpf_class pi_;

 private:
  void drm(const int64_t n0,
           const int64_t n1,
           mpz_class& x0,
           mpz_class& y0,
           mpz_class& z0,
           bool need_z);
  void parallel_drm(const int64_t n0,
                    const int64_t n1,
                    mpz_class& x0,
                    mpz_class& y0,
                    mpz_class& z0,
                    bool need_z,
                    const int number_of_threads);
  virtual void setXYZ(int64_t k, mpz_class& x, mpz_class& y, mpz_class& z) = 0;
  virtual void postProcess(mpz_class& x, mpz_class& y) = 0;
  virtual int64_t terms(int64_t digits) const = 0;
  virtual const char* name() const = 0;

  const Configuration config_;
};

class Chudnovsky : public Computer {
 public:
  Chudnovsky(const Configuration& config) : Computer(config) {}
  ~Chudnovsky() override = default;

 private:
  void setXYZ(int64_t k, mpz_class& x, mpz_class& y, mpz_class& z) override;
  void postProcess(mpz_class& x, mpz_class& y) override;
  int64_t terms(int64_t digits) const override { return digits / 14; }
  const char* name() const override { return "chudnovsky"; }
};

class Ramanujan : public Computer {
 public:
  Ramanujan(const Configuration& config) : Computer(config) {}
  ~Ramanujan() override = default;

 private:
  void setXYZ(int64_t k, mpz_class& x, mpz_class& y, mpz_class& z) override;
  void postProcess(mpz_class& x, mpz_class& y) override;
  int64_t terms(int64_t digits) const override { return digits / 7.98; }
  const char* name() const override { return "ramanujan"; }
};
