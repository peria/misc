#pragma once

struct Complex {
  Complex() = default;
  Complex(const Complex&) = default;
  Complex(Complex&&) = default;
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
