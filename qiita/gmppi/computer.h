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
  struct PrimeFactor {
    int64_t prime = 0;
    int64_t exp = 0;
  };
  struct Factorized {
    void add(const PrimeFactor& f) { factors.push_back(f); }
    Factorized& operator*=(const Factorized& other);
    mpz_class toMpz() const;
    void shrink();
    std::vector<PrimeFactor> factors;

   private:
    mpz_class toMpz(int64_t a, int64_t b) const;
  };
  struct Parameter {
    mpz_class x;
    mpz_class y;
    mpz_class z;
    Factorized factorized_x;
    Factorized factorized_z;
  };

  explicit Computer(const Configuration& config, int64_t max_n)
      : config_(config), factorizer_(max_n) {}
  Factorized factorize(int64_t n, int64_t e = 1) {
    return factorizer_.factorize(n, e);
  }

  mpf_class pi_;
  const Configuration config_;

 private:
  class Factorizer {
   public:
    explicit Factorizer(const int64_t n);
    Factorized factorize(int64_t n, int64_t exp);

   private:
    struct Element {
      PrimeFactor prime_factor;
      int64_t next = 0;
    };
    std::vector<Element> elements_;
  };

  void drm(const int64_t n0, const int64_t n1, Parameter& param, bool need_z);
  void parallelDrm(const int64_t n0,
                   const int64_t n1,
                   Parameter& param,
                   bool need_z,
                   const int number_of_threads);
  void reduceGcd(Parameter& p0, Parameter& p1);
  Factorized splitGcd(Factorized& a, Factorized& b);
  virtual void setXYZ(int64_t k, Parameter& param) = 0;
  virtual void postProcess(Parameter& param) = 0;
  virtual int64_t terms() const = 0;
  virtual const char* name() const = 0;

  Factorizer factorizer_;
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
