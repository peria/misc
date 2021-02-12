#pragma once

#include <cstdint>

using int16 = int16_t;
using int64 = int64_t;
using uint64 = uint64_t;
using uint128 = __uint128_t;

struct Complex {
  Complex() = default;
  Complex(const Complex&) = default;
  Complex(Complex&&) = default;
  Complex(double re, double im) : real(re), imag(im) {}
  Complex& operator=(const Complex&) = default;
  Complex& operator=(Complex&&) = default;
  Complex& operator*=(const Complex& a) {
    double r = real * a.real - imag * a.imag;
    double i = real * a.imag + imag * a.real;
    real = r;
    imag = i;
    return *this;
  }
  Complex& operator*=(const double a) {
    real *= a;
    imag *= a;
    return *this;
  }
  Complex i() const { return Complex{-imag, real}; }
  Complex conj() const { return Complex{real, -imag}; }
  Complex tran() const { return Complex{imag, real}; }

  double real = 0;
  double imag = 0;
};

inline Complex operator+(const Complex& a, const Complex& b) {
  return Complex{a.real + b.real, a.imag + b.imag};
}

inline Complex operator-(const Complex& a, const Complex& b) {
  return Complex{a.real - b.real, a.imag - b.imag};
}

inline Complex operator*(const Complex& a, const Complex& b) {
  return Complex{a.real * b.real - a.imag * b.imag,
                 a.real * b.imag + a.imag * b.real};
}

inline Complex operator*(const Complex& a, const double b) {
  return Complex{a.real * b, a.imag * b};
}
