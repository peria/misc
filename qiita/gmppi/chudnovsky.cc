#include "computer.h"

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
