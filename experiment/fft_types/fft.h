#pragma once

#include "complex.h"

class FFT {
 public:
  virtual ~FFT() {}
  virtual const char* name() const = 0;
  virtual void setUp(int) = 0;
  virtual void tearDown() {}
  virtual void dft(Complex*) = 0;
  virtual void idft(Complex*) = 0;
};
