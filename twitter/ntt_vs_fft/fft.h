#pragma once

#include "fmt.h"

// RFT class computes DFT in real field.
class RFT {
 public:
  static void SetGlobals(void* work, void* table, void*);
  static bool Validate(double* data);
  static void Square(double* data, int n);

  static void InitTable(const int log2n);
  static void Forward(const int log2n, const int real_n, double* x);
  static void Backward(const int log2n, const int real_n, double* x);
};

// FFT class computes DFT in complex field.
class FFT {
 public:
  static void InitTable(const int log2n);
  static void Forward(const int log2n, const int n, double* data);
  static void Backward(const int log2n, const int n, double* data);

 private:
  static void Core(const int log2n, const int n, double* table, double* work, double* data);
  static void Radix2(const int width, const int height,
                     double* ptr, double* x, double* y);
};
