#pragma once

#include <cmath>
#include <vector>

#include "complex.h"
#include "fft.h"

class Stockham6 final : public FFT {
  static constexpr int64 kThreshold = 10;
  enum class Step {
    kSingle,
    kSix,
    kNine,
  };
  static Step getStep(int64 log2k) {
    if (log2k >= kThreshold)
      return (log2k % 2) ? Step::kNine : Step::kSix;
    return Step::kSingle;
  }

 public:
  Stockham6(int64 log2k)
    : FFT(log2k % 2, log2k / 2),
      step(getStep(log2k)) {
    init();
  }
  ~Stockham6() override = default;

  void dft(Complex* x, bool backward) const override {
    switch (step) {
    case Step::kSingle: {
      if (backward)
        dft<true>(x, log2n_, log4n_);
      else
        dft<false>(x, log2n_, log4n_);
      break;
    }
    case Step::kSix: {
      int64 logm = logn_ / 2;
      int64 log2m = logm % 2;
      int64 log4m = logm / 2;
      if (backward)
        idft6Step(x, log2m, log4m);
      else
        dft6Step(x, log2m, log4m);
      break;
    }
    case Step::kNine: {
      int64 logm = logn_ / 2;
      int64 log2m = logm % 2;
      int64 log4m = logm / 2;
      if (backward)
        idft9Step(x, log2m, log4m);
      else
        dft9Step(x, log2m, log4m);
      break;
    }
    }

    if (backward) {
      double inv = 1.0 / n_;
      for (int i = 0; i < n_; ++i) {
        x[i] *= inv;
      }
    }
  }

  static const char* name() { return "6StepStock"; }

 private:
  void init();

  template <bool backward>
  void dft(Complex* a, const int64 log2m, const int64 log4m) const {
    const int64 logm = log2m + log4m * 2;
    const int64 m = 1LL << logm;

    Complex* x = a;
    Complex* y = const_cast<Complex*>(work.data());
    const Complex* pw = ws.data();
    int64 width = 1;
    int64 height = m;
    for (int64 i = 0; i < log4m - 1; ++i) {
      height /= 4;
      dft4<backward>(x, y, width, height, pw);
      pw += height * 3;
      width *= 4;
      std::swap(x, y);
    }
    {
      height /= 4;
      dft4<backward>(x, log2n_ ? y : a, width, height, pw);
      width *= 4;
    }
    if (log2m) {
      height /= 2;
      dft2(y, a, width);
      width *= 2;
    }
  }

  void dft6Step(Complex* x, const int64 log2m, const int64 log4m) const {
    const int64 logm = log2m + log4m * 2;
    const int64 m = 1LL << logm;
    // Transpose
    for (int64 i = 0; i < m; ++i) {
      for (int64 j = i + 1; j < m; ++j) {
        std::swap(x[i * m + j], x[j * m + i]);
      }
      // Sub FFT
      dft<false>(&x[i * m], log2m, log4m);
    }
    // Transpose, twiddle factor
    for (int64 i = 0; i < m; ++i) {
      x[i * m + i] *= twiddle[i * m + i];
      for (int64 j = i + 1; j < m; ++j) {
        const Complex& tw = twiddle[i * m + j];
        std::swap(x[i * m + j], x[j * m + i]);
        x[i * m + j] *= tw;
        x[j * m + i] *= tw;
      }
    }
    // Sub FFT
    for (int64 i = 0; i < m; ++i) {
      dft<false>(&x[i * m], log2m, log4m);
    }
    // no transpose
  }

  void idft6Step(Complex* x, const int64 log2m, const int64 log4m) const {
    const int64 logm = log2m + log4m * 2;
    const int64 m = 1LL << logm;
    // no transpose
    // Sub FFT
    for (int64 i = 0; i < m; ++i) {
      dft<true>(&x[i * m], log2m, log4m);
    }
    // Transpose, twiddle factor
    for (int64 i = 0; i < m; ++i) {
      x[i * m + i] *= twiddle[i * m + i].conj();
      for (int64 j = i + 1; j < m; ++j) {
        const Complex& tw = twiddle[i * m + j].conj();
        std::swap(x[i * m + j], x[j * m + i]);
        x[i * m + j] *= tw;
        x[j * m + i] *= tw;
      }
    }
    // Sub FFT
    for (int64 i = 0; i < m; ++i) {
      dft<true>(&x[i * m], log2m, log4m);
      // Transpose
      for (int64 j = i + 1; j < m; ++j) {
        std::swap(x[i * m + j], x[j * m + i]);
      }
    }
  }

  void dft9Step(Complex* a, const int64 log2m, const int64 log4m) const {
    const int64 logm = log2m + log4m * 2;
    const int64 m2 = 1LL << (logm * 2);
    // n == 2 * m2

    // twiddle[0:m2] is used for 6 step FFT.
    for (int64 k = 0; k < m2; ++k) {
      Complex w = twiddle[k + m2];
      int64 i0 = k;
      int64 i1 = k + m2;
      Complex a0 = a[i0];
      Complex a1 = a[i1];
      a[i0] = a0 + a1;
      a[i1] = (a0 - a1) * w;
    }
    dft6Step(&a[0], log2m, log4m);
    dft6Step(&a[m2], log2m, log4m);
  }

  void idft9Step(Complex* a, const int64 log2m, const int64 log4m) const {
    const int64 logm = log2m + log4m * 2;
    const int64 m2 = 1LL << (logm * 2);
    // n == 2 * m2

    idft6Step(&a[0], log2m, log4m);
    idft6Step(&a[m2], log2m, log4m);
    // twiddle[0:n/2] is used for 6 step FFT.
    for (int64 k = 0; k < m2; ++k) {
      Complex w = twiddle[k + m2].conj();
      int64 i0 = k;
      int64 i1 = k + m2;
      Complex a0 = a[i0];
      Complex a1 = a[i1] * w;
      a[i0] = a0 + a1;
      a[i1] = a0 - a1;
    }
  }

  void dft2(Complex* x, Complex* y, const int64 width) const {
    for (int64 j = 0; j < width; ++j) {
      int64 i0 = j;
      int64 i1 = j + width;
      Complex x0 = x[i0];
      Complex x1 = x[i1];
      y[i0] = x0 + x1;
      y[i1] = x0 - x1;
    }
  }

  template <bool backward>
  void dft4(Complex* x, Complex* y, const int64 width, const int64 height, const Complex* pw) const {
    {
      for (int64 j = 0; j < width; ++j) {
        int64 ix0 = j;
        int64 ix1 = height * width + j;
        int64 ix2 = height * width * 2 + j;
        int64 ix3 = height * width * 3 + j;
        int64 iy0 = j;
        int64 iy1 = width + j;
        int64 iy2 = 2 * width + j;
        int64 iy3 = 3 * width + j;
        Complex x0 = x[ix0];
        Complex x1 = x[ix1];
        Complex x2 = x[ix2];
        Complex x3 = x[ix3];
        Complex b0 = x0 + x2;
        Complex b1 = x1 + x3;
        Complex b2 = x0 - x2;
        Complex b3 = (x3 - x1).i();
        y[iy0] = b0 + b1;
        y[iy1] = backward ? (b2 - b3) : (b2 + b3);
        y[iy2] = b0 - b1;
        y[iy3] = backward ? (b2 + b3) : (b2 - b3);
      }
    }
    for (int64 k = 1; k < height; ++k) {
      Complex w1 = backward ? pw[k * 3].conj() : pw[k * 3];
      Complex w2 = backward ? pw[k * 3 + 1].conj() : pw[k * 3 + 1];
      Complex w3 = backward ? pw[k * 3 + 2].conj() : pw[k * 3 + 2];
      for (int64 j = 0; j < width; ++j) {
        int64 ix0 = k * width + j;
        int64 ix1 = (k + height) * width + j;
        int64 ix2 = (k + height * 2) * width + j;
        int64 ix3 = (k + height * 3) * width + j;
        int64 iy0 = 4 * k * width + j;
        int64 iy1 = (4 * k + 1) * width + j;
        int64 iy2 = (4 * k + 2) * width + j;
        int64 iy3 = (4 * k + 3) * width + j;
        Complex x0 = x[ix0];
        Complex x1 = x[ix1];
        Complex x2 = x[ix2];
        Complex x3 = x[ix3];
        Complex b0 = x0 + x2;
        Complex b1 = x1 + x3;
        Complex b2 = x0 - x2;
        Complex b3 = (x3 - x1).i();
        y[iy0] = b0 + b1;
        y[iy1] = (backward ? (b2 - b3) : (b2 + b3)) * w1;
        y[iy2] = (b0 - b1) * w2;
        y[iy3] = (backward ? (b2 + b3) : (b2 - b3)) * w3;
      }
    }
  }

  const Step step;
  std::vector<Complex> twiddle;
  std::vector<Complex> ws;
  std::vector<Complex> work;
};

void Stockham6::init() {
  work.resize(n_);

  int64 log4m = log4n_;
  int64 m = n_;
  if (step != Step::kSingle) {
    log4m = logn_ / 4;
    m = 1LL << (logn_ / 2);
  }

  {
    double theta = -2 * M_PI / m;
    int64 width = 1;
    int64 height = m;
    for (int64 i = 0; i < log4m; ++i) {
      height /= 4;
      for (int64 k = 0; k < height; ++k) {
        const double t = theta * width * k;
        ws.push_back(Complex {std::cos(t), std::sin(t)});
        ws.push_back(Complex {std::cos(2*t), std::sin(2*t)});
        ws.push_back(Complex {std::cos(3*t), std::sin(3*t)});
      }
      width *= 4;
    }
  }

  if (step == Step::kSingle)
    return;

  {
    double theta = -2 * M_PI / (m * m);
    for (int64 k = 0; k < m; ++k) {
      for (int64 j = 0; j < m; ++j) {
        double t = theta * k * j;
        twiddle.push_back(Complex {std::cos(t), std::sin(t)});
      }
    }
  }

  if (step == Step::kNine) {
    if (twiddle.size() != n_ / 2) {
      std::cerr << "twiddle.size() != n / 2\n"
                << twiddle.size() << " - " << n_ / 2 << "\n";
    }
    double theta = -2 * M_PI / n_;
    for (int64 k = 0; k < n_ / 2; ++k) {
      double t = theta * k;
      twiddle.push_back(Complex {std::cos(t), std::sin(t)});
    }
  }
}
