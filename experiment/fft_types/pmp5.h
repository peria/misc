#pragma once

#include <cmath>
#include <vector>

#include "complex.h"
#include "fft.h"

// Extend PMP with 6 step fft, but do not transpose at last.
class PMP5 : public FFT {
  static constexpr int64 kThreshold = 13;
  static constexpr int64 kPadding = 4;

 public:
  PMP5(int64 logn)
    : FFT(logn % 2, logn / 2),
      logn1(logn > kThreshold ? logn / 2 : 0),
      log2n1(logn1 % 2),
      log4n1(logn1 / 2),
      n1(1LL << logn1),
      logn2(logn - logn1),
      log2n2(logn2 % 2),
      log4n2(logn2 / 2),
      n2(1LL << logn2) {
    init(log2n2, log4n2);
  }
  ~PMP5() override = default;

  static const char* name() { return "PMP5"; }

  void dft(Complex* x, bool backward) const override {
    if (backward) {
      idft(x);
      double inv = 1.0 / n_;
      for (int64 i = 0; i < n_; ++i)
        x[i] *= inv;
    } else {
      dft(x);
    }
  }

 private:
  void init(int64 log2n, int64 log4n);

  void dft(Complex* a) const {
    const Complex* pw = ws.data();
    if (logn1) {
      const Complex* pw1 = pw + (n2 - log2n2) - (n1 - log2n1);
      Complex* b0 = const_cast<Complex*>(buffer.data());
      Complex* b1 = b0 + n2 + kPadding;
      Complex* b2 = b1 + n2 + kPadding;
      Complex* b3 = b2 + n2 + kPadding;
      for (int64 c = 0; c < n2; c += 2) {
        for (int64 r = 0; r < n1; ++r) {
          b0[r] = a[n2 * r + c];
          b1[r] = a[n2 * r + c + 1];
        }
        dft1d(b0, log2n1, log4n1, pw1);
        dft1d(b1, log2n1, log4n1, pw1);
        for (int64 r = 0; r < n1; ++r) {
          a[n2 * r + c] = b0[r];
          a[n2 * r + c + 1] = b1[r];
        }
      }
      dft1d(a, log2n2, log4n2, pw);
      for (int64 r = 1; r < n1; ++r) {
        for (int64 c = 1; c < n2; ++c) {
          a[n2 * r + c] *= twiddle[n2 * r + c];
        }
        dft1d(a + n2 * r, log2n2, log4n2, pw);
      }
    } else {
      dft1d(a, log2n2, log4n2, pw);
    }
  }

  void dft1d(Complex* a, int64 log2n, int64 log4n, const Complex* pw) const {
    const int64 logn = log2n + log4n * 2;
    const int64 n = 1LL << logn;
    int64 l = 1;
    int64 m = n;
    for (int64 i = 0; i < log4n; ++i) {
      m /= 4;
      dft4(a, l, m, pw);
      pw += m * 3;
      l *= 4;
    }
    for (int64 i = 0; i < log2n; ++i) {
      m /= 2;
      dft2(a, l);
      l *= 2;
    }
  }

  void idft(Complex* a) const {
    const Complex* pw = ws.data() + ws.size();
    if (logn1) {
      idft1d(a, log2n2, log4n2, pw);
      for (int64 r = 1; r < n1; ++r) {
        idft1d(a + n2 * r, log2n2, log4n2, pw);
        for (int64 c = 1; c < n2; ++c) {
          a[n2 * r + c] *= twiddle[n2 * r + c].conj();
        }
      }
      Complex* b0 = const_cast<Complex*>(buffer.data());
      Complex* b1 = b0 + n2 + kPadding;
      for (int64 c = 0; c < n2; c += 2) {
        for (int64 r = 0; r < n1; ++r) {
          b0[r] = a[n2 * r + c];
          b1[r] = a[n2 * r + c + 1];
        }
        idft1d(b0, log2n1, log4n1, pw);
        idft1d(b1, log2n1, log4n1, pw);
        for (int64 r = 0; r < n1; ++r) {
          a[n2 * r + c] = b0[r];
          a[n2 * r + c + 1] = b1[r];
        }
      }
    } else {
      idft1d(a, log2n2, log4n2, pw);
    }
  }

  void idft1d(Complex* a, int64 log2n, int64 log4n, const Complex* pw) const {
    const int64 logn = log2n + log4n * 2;
    const int64 n = 1LL << logn;
    int64 l = n;
    int64 m = 1;
    if (log2n) {
      l /= 2;
      dft2(a, l);
      m *= 2;
    }
    for (int64 i = 0; i < log4n; ++i) {
      l /= 4;
      pw -= m * 3;
      idft4(a, l, m, pw);
      m *= 4;
    }
  }

  void dft2(Complex* a, const int64 l) const {
    for (int64 j = 0; j < l; ++j) {
      int64 i0 = 2 * j;
      int64 i1 = 2 * j + 1;
      Complex a0 = a[i0];
      Complex a1 = a[i1];
      a[i0] = a0 + a1;
      a[i1] = a0 - a1;
    }
  }

  void dft4(Complex* a, const int64 l, const int64 m, const Complex* ws) const {
    for (int64 j = 0; j < l; ++j) {
      {
        int64 i0 = 4 * j * m;
        int64 i1 = 4 * j * m + m;
        int64 i2 = 4 * j * m + 2 * m;
        int64 i3 = 4 * j * m + 3 * m;
        Complex a0 = a[i0];
        Complex a1 = a[i1];
        Complex a2 = a[i2];
        Complex a3 = a[i3];
        Complex b0 = a0 + a2;
        Complex b1 = a1 + a3;
        Complex b2 = a0 - a2;
        Complex b3 = (a3 - a1).i();
        a[i0] = b0 + b1;
        a[i1] = b0 - b1;
        a[i2] = b2 + b3;
        a[i3] = b2 - b3;
      }
      for (int64 k = 1; k < m; ++k) {
        Complex w1 = ws[k * 3];
        Complex w2 = ws[k * 3 + 1];
        Complex w3 = ws[k * 3 + 2];
        int64 i0 = 4 * j * m + k;
        int64 i1 = 4 * j * m + m + k;
        int64 i2 = 4 * j * m + 2 * m + k;
        int64 i3 = 4 * j * m + 3 * m + k;
        Complex a0 = a[i0];
        Complex a1 = a[i1];
        Complex a2 = a[i2];
        Complex a3 = a[i3];
        Complex b0 = a0 + a2;
        Complex b1 = a1 + a3;
        Complex b2 = a0 - a2;
        Complex b3 = (a3 - a1).i();
        a[i0] = b0 + b1;
        a[i1] = (b0 - b1) * w2;
        a[i2] = (b2 + b3) * w1;
        a[i3] = (b2 - b3) * w3;
      }
    }
  }

  void idft4(Complex* a, const int64 l, const int64 m, const Complex* ws) const {
    for (int64 j = 0; j < l; ++j) {
      {
        int64 i0 = 4 * j * m;
        int64 i1 = 4 * j * m + m;
        int64 i2 = 4 * j * m + 2 * m;
        int64 i3 = 4 * j * m + 3 * m;
        Complex a0 = a[i0];
        Complex a1 = a[i1];
        Complex a2 = a[i2];
        Complex a3 = a[i3];
        Complex b0 = a0 + a1;
        Complex b1 = a0 - a1;
        Complex b2 = a2 + a3;
        Complex b3 = (a2 - a3).i();
        a[i0] = b0 + b2;
        a[i1] = b1 + b3;
        a[i2] = b0 - b2;
        a[i3] = b1 - b3;
      }
      for (int64 k = 1; k < m; ++k) {
        Complex w1 = ws[3 * k].conj();
        Complex w2 = ws[3 * k + 1].conj();
        Complex w3 = ws[3 * k + 2].conj();
        int64 i0 = 4 * j * m + k;
        int64 i1 = 4 * j * m + m + k;
        int64 i2 = 4 * j * m + 2 * m + k;
        int64 i3 = 4 * j * m + 3 * m + k;
        Complex a0 = a[i0];
        Complex a1 = a[i1] * w2;
        Complex a2 = a[i2] * w1;
        Complex a3 = a[i3] * w3;
        Complex b0 = a0 + a1;
        Complex b1 = a0 - a1;
        Complex b2 = a2 + a3;
        Complex b3 = (a2 - a3).i();
        a[i0] = b0 + b2;
        a[i1] = b1 + b3;
        a[i2] = b0 - b2;
        a[i3] = b1 - b3;
      }
    }
  }

  std::vector<Complex> ws;
  std::vector<Complex> twiddle;
  std::vector<Complex> buffer;
  const int64 logn1 = 0;
  const int64 log2n1 = 0;
  const int64 log4n1 = 0;
  const int64 n1 = 1;
  const int64 logn2 = 0;
  const int64 log2n2 = 0;
  const int64 log4n2 = 0;
  const int64 n2 = 1;
};

void PMP5::init(int64 log2n, int64 log4n) {
  {
    const int64 logn = log2n + log4n * 2;
    const int64 n = 1LL << logn;
    const double theta = -2 * M_PI / n;
    int64 l = 1;
    int64 m = n;
    for (int64 i = 0; i < log4n; ++i) {
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

  if (logn1) {
    std::vector<int64> rs(n1);
    int64 rr = 0;
    for (int64 i = 1; i < n1; ++i) {
      for(int64 j = (1LL << (logn1 - 1)); j > (rr ^= j); j >>= 1);
      rs[i] = rr;
    }

    // twiddle[r*n2 + c] = w^(bitrev(r)*c)  0<=r<n1, 0<=c<n2
    const double theta = -2 * M_PI / n_;
    for (int64 r = 0; r < n1; ++r) {
      for (int64 c = 0; c < n2; ++c) {
        double t = theta * rs[r] * c;
        twiddle.push_back(Complex {std::cos(t), std::sin(t)});
      }
    }

    buffer.resize((n2 + kPadding) * 4);
  }
}
