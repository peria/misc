#pragma once

#include <cstdint>

#include <gmpxx.h>

class Computer {
 public:
  Computer(int64_t digits) : digits_(digits) {}
  virtual ~Computer() = default;

  void compute();
  void output();
  void check(const char*);

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
  virtual const char* name() const = 0;

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
  const char* name() const override { return "chudnovsky"; }
};

class Ramanujan : public Computer {
 public:
  Ramanujan(int64_t digits) : Computer(digits) {}
  ~Ramanujan() override = default;

 private:
  void setXYZ(int64_t k, mpz_class& x, mpz_class& y, mpz_class& z) override;
  void postProcess(mpz_class& x, mpz_class& y) override;
  int64_t terms(int64_t digits) const override { return digits / 7.98; }
  const char* name() const override { return "ramanujan"; }
};
