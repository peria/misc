#pragma once

#include <cstdint>
#include <memory>
#include <optional>
#include <string>
#include <vector>

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
  struct Parameter {
    mpz_class x;
    mpz_class y;
    mpz_class z;
  };
  struct PrimeFactor {
    int64_t prime = 0;
    int64_t exp = 0;
  };
  class Sieve {
   public:
    explicit Sieve(const int64_t n);

   private:
    struct Element {
      PrimeFactor prime_factor;
      int64_t next = 0;
    };
    std::vector<Element> elements_;
  };

  explicit Computer(const Configuration& config, int64_t sieve_size)
      : config_(config), sieve_(sieve_size) {}

  mpf_class pi_;
  const Configuration config_;

 private:
  void drm(const int64_t n0, const int64_t n1, Parameter& param, bool need_z);
  void parallel_drm(const int64_t n0,
                    const int64_t n1,
                    Parameter& param,
                    bool need_z,
                    const int number_of_threads);
  virtual void setXYZ(int64_t k, Parameter& param) = 0;
  virtual void postProcess(Parameter& param) = 0;
  virtual int64_t terms() const = 0;
  virtual const char* name() const = 0;

  Sieve sieve_;
};

class Chudnovsky : public Computer {
 public:
  explicit Chudnovsky(const Configuration& config);
  ~Chudnovsky() override = default;

 private:
  void setXYZ(int64_t k, Parameter& param) override;
  void postProcess(Parameter& param) override;
  int64_t terms() const override { return terms(config_.digits); }
  const char* name() const override { return "chudnovsky"; }

  static int64_t terms(int64_t digits) { return digits / 14; }
};

class Ramanujan : public Computer {
 public:
  explicit Ramanujan(const Configuration& config);
  ~Ramanujan() override = default;

 private:
  void setXYZ(int64_t k, Parameter& param) override;
  void postProcess(Parameter& param) override;
  int64_t terms() const override { return terms(config_.digits); }
  const char* name() const override { return "ramanujan"; }

  static int64_t terms(int64_t digits) { return digits / 7.98; }
};
