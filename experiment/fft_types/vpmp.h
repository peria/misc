#pragma once

#include <cmath>
#include <vector>

#include "complex.h"
#include "fft.h"

// Store the vector of complex numbers as struct-of-vector.
class VPMP : public FFT {
 public:
  VPMP(int64 log2k)
    : FFT(log2k % 2, log2k / 2) {
    init();
  }
  ~VPMP() override = default;

  static const char* name() { return "VPMP"; }

  void dft(Complex* x, bool backward) const override {
    if (backward)
      idft(x);
    else
      dft(x);
  }

 private:
  void init();

  void dft(Complex* a) const {
    double* vreal = reinterpret_cast<double*>(a);
    double* vimag = vreal + n_;

    auto* pw = ws.data();
    int64 l = 1;
    int64 m = n_;
    for (int64 i = 0; i < log4n_; ++i) {
      m /= 4;
      dft4(vreal, vimag, l, m, pw);
      pw += m * 3;
      l *= 4;
    }
    for (int64 i = 0; i < log2n_; ++i) {
      m /= 2;
      dft2(vreal, vimag, l);
      l *= 2;
    }
  }

  void idft(Complex* a) const {
    double* vreal = reinterpret_cast<double*>(a);
    double* vimag = vreal + n_;

    auto* pw = ws.data() + ws.size();
    int64 l = n_;
    int64 m = 1;
    if (log2n_) {
      l /= 2;
      dft2(vreal, vimag, l);
      m *= 2;
    }
    for (int64 i = 0; i < log4n_; ++i) {
      l /= 4;
      pw -= m * 3;
      idft4(vreal, vimag, l, m, pw);
      m *= 4;
    }

    double inv = 1.0 / n_;
    for (int64 i = 0; i < n_; ++i) {
      a[i] *= inv;
    }
  }

  void dft2(double* vr, double* vi, const int64 l) const {
    for (int64 j = 0; j < l; ++j) {
      int64 i0 = 2 * j;
      int64 i1 = 2 * j + 1;
      double a0r = vr[i0];
      double a0i = vi[i0];
      double a1r = vr[i1];
      double a1i = vi[i1];
      vr[i0] = a0r + a1r;
      vr[i1] = a0r - a1r;
      vi[i0] = a0i + a1i;
      vi[i1] = a0i - a1i;
    }
  }

  void dft4(double* vr, double* vi, const int64 l, const int64 m, const Complex* ws) const {
    for (int64 j = 0; j < l; ++j) {
      {
        int64 i0 = 4 * j * m;
        int64 i1 = 4 * j * m + m;
        int64 i2 = 4 * j * m + 2 * m;
        int64 i3 = 4 * j * m + 3 * m;
        double a0r = vr[i0];
        double a0i = vi[i0];
        double a1r = vr[i1];
        double a1i = vi[i1];
        double a2r = vr[i2];
        double a2i = vi[i2];
        double a3r = vr[i3];
        double a3i = vi[i3];
        double b0r = a0r + a2r;
        double b0i = a0i + a2i;
        double b1r = a1r + a3r;
        double b1i = a1i + a3i;
        double b2r = a0r - a2r;
        double b2i = a0i - a2i;
        double b3r = -(a3i - a1i);
        double b3i = a3r - a1r;
        vr[i0] = b0r + b1r;
        vi[i0] = b0i + b1i;
        vr[i1] = b0r - b1r;
        vi[i1] = b0i - b1i;
        vr[i2] = b2r + b3r;
        vi[i2] = b2i + b3i;
        vr[i3] = b2r - b3r;
        vi[i3] = b2i - b3i;
      }
      for (int64 k = 1; k < m; ++k) {
        double w1r = ws[k * 3].real;
        double w1i = ws[k * 3].imag;
        double w2r = ws[k * 3 + 1].real;
        double w2i = ws[k * 3 + 1].imag;
        double w3r = ws[k * 3 + 2].real;
        double w3i = ws[k * 3 + 2].imag;
        int64 i0 = 4 * j * m + k;
        int64 i1 = 4 * j * m + m + k;
        int64 i2 = 4 * j * m + 2 * m + k;
        int64 i3 = 4 * j * m + 3 * m + k;
        double a0r = vr[i0];
        double a0i = vi[i0];
        double a1r = vr[i1];
        double a1i = vi[i1];
        double a2r = vr[i2];
        double a2i = vi[i2];
        double a3r = vr[i3];
        double a3i = vi[i3];
        double b0r = a0r + a2r;
        double b0i = a0i + a2i;
        double b1r = a1r + a3r;
        double b1i = a1i + a3i;
        double b2r = a0r - a2r;
        double b2i = a0i - a2i;
        double b3r = -(a3i - a1i);
        double b3i = (a3r - a1r);
        vr[i0] = b0r + b1r;
        vi[i0] = b0i + b1i;
        double c1r = b0r - b1r;
        double c1i = b0i - b1i;
        double c2r = b2r + b3r;
        double c2i = b2i + b3i;
        double c3r = b2r - b3r;
        double c3i = b2i - b3i;
        vr[i1] = c1r * w2r - c1i * w2i;
        vi[i1] = c1r * w2i + c1i * w2r;
        vr[i2] = c2r * w1r - c2i * w1i;
        vi[i2] = c2r * w1i + c2i * w1r;
        vr[i3] = c3r * w3r - c3i * w3i;
        vi[i3] = c3r * w3i + c3i * w3r;
      }
    }
  }

  void idft4(double* vr, double* vi, const int64 l, const int64 m, const Complex* ws) const {
    for (int64 j = 0; j < l; ++j) {
      {
        int64 i0 = 4 * j * m;
        int64 i1 = 4 * j * m + m;
        int64 i2 = 4 * j * m + 2 * m;
        int64 i3 = 4 * j * m + 3 * m;
        double a0r = vr[i0];
        double a0i = vi[i0];
        double a1r = vr[i1];
        double a1i = vi[i1];
        double a2r = vr[i2];
        double a2i = vi[i2];
        double a3r = vr[i3];
        double a3i = vi[i3];
        double b0r = a0r + a1r;
        double b0i = a0i + a1i;
        double b1r = a0r - a1r;
        double b1i = a0i - a1i;
        double b2r = a2r + a3r;
        double b2i = a2i + a3i;
        double b3r = -(a2i - a3i);
        double b3i = a2r - a3r;
        vr[i0] = b0r + b2r;
        vi[i0] = b0i + b2i;
        vr[i1] = b1r + b3r;
        vi[i1] = b1i + b3i;
        vr[i2] = b0r - b2r;
        vi[i2] = b0i - b2i;
        vr[i3] = b1r - b3r;
        vi[i3] = b1i - b3i;
      }
      for (int64 k = 1; k < m; ++k) {
        double w1r = ws[3 * k].real;
        double w1i = -ws[3 * k].imag;
        double w2r = ws[3 * k + 1].real;
        double w2i = -ws[3 * k + 1].imag;
        double w3r = ws[3 * k + 2].real;
        double w3i = -ws[3 * k + 2].imag;
        int64 i0 = 4 * j * m + k;
        int64 i1 = 4 * j * m + m + k;
        int64 i2 = 4 * j * m + 2 * m + k;
        int64 i3 = 4 * j * m + 3 * m + k;
        double a0r = vr[i0];
        double a0i = vi[i0];
        double a1r = vr[i1];
        double a1i = vi[i1];
        double a2r = vr[i2];
        double a2i = vi[i2];
        double a3r = vr[i3];
        double a3i = vi[i3];
        double b1r = a1r * w2r - a1i * w2i;
        double b1i = a1r * w2i + a1i * w2r;
        double b2r = a2r * w1r - a2i * w1i;
        double b2i = a2r * w1i + a2i * w1r;
        double b3r = a3r * w3r - a3i * w3i;
        double b3i = a3r * w3i + a3i * w3r;
        double c0r = a0r + b1r;
        double c0i = a0i + b1i;
        double c1r = a0r - b1r;
        double c1i = a0i - b1i;
        double c2r = b2r + b3r;
        double c2i = b2i + b3i;
        double c3r = -(b2i - b3i);
        double c3i = b2r - b3r;
        vr[i0] = c0r + c2r;
        vi[i0] = c0i + c2i;
        vr[i1] = c1r + c3r;
        vi[i1] = c1i + c3i;
        vr[i2] = c0r - c2r;
        vi[i2] = c0i - c2i;
        vr[i3] = c1r - c3r;
        vi[i3] = c1i - c3i;
      }
    }
  }

  std::vector<Complex> ws;
};

void VPMP::init() {
  int64 l = 1;
  int64 m = n_;
  const double theta = -2 * M_PI / n_;
  for (int64 i = 0; i < log4n_; ++i) {
    m /= 4;
    for (int64 k = 0; k < m; ++k) {
      double t = theta * l * k;
      ws.push_back(Complex {std::cos(t), std::sin(t)});
      ws.push_back(Complex {std::cos(2*t), std::sin(2*t)});
      ws.push_back(Complex {std::cos(3*t), std::sin(3*t)});
    }
    l *= 4;
  }
}
