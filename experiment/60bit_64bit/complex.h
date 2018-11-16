#pragma once

struct Complex {
  double real;
  double imag;

  inline Complex i() const { return {imag, -real}; }
};

inline Complex operator+(const Complex& lhs, const Complex& rhs) {
  return {lhs.real + rhs.real, lhs.imag + rhs.imag};
}

inline Complex operator-(const Complex& lhs, const Complex& rhs) {
  return {lhs.real - rhs.real, lhs.imag - rhs.imag};
}

inline Complex operator*(const Complex& lhs, const Complex& rhs) {
  return {lhs.real * rhs.real - lhs.imag * rhs.imag,
          lhs.real * rhs.imag + lhs.imag * rhs.real};
}

inline Complex operator*(const double& lhs, const Complex& rhs) {
  return {lhs * rhs.real, lhs * rhs.imag};
}
