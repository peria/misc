#include "computer.h"

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
