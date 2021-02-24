#pragma once

#include <algorithm>
#include <cmath>
#include <vector>

#include "complex.h"
#include "fft.h"

// Extend PMP to repack data with padding.
class PMP2 : public FFT {
  static constexpr int64 kLogBlockSize = 16;
  static constexpr int64 kBlockSize = 1LL << kLogBlockSize;
  static constexpr int64 kPaddingSize = 1;
  static constexpr int64 kStrideSize = kBlockSize + kPaddingSize;
  static_assert(kLogBlockSize % 2 == 0, "kLogBlockSize cannot be odd.");

 public:
  PMP2(int64 log2k)
    : FFT(log2k % 2, log2k / 2),
      log2_block(std::max<int64>(log2k - kLogBlockSize, 0)),
      log2_point(std::min<int64>(log2k, kLogBlockSize)) {
    init();
  }
  ~PMP2() override = default;

  static const char* name() { return "PMP2"; }

  void dft(Complex* x, bool backward) const override {
    if (log2_block) {
      repackWithPadding(x);
      dftWithoutRepack(const_cast<Complex*>(work.data()), backward);
      extendWithPadding(x);
    } else {
      dftWithoutRepack(x, backward);
    }
  }

  void dftWithoutRepack(Complex* x, bool backward) const {
    if (backward) {
      idft(x);
      const double inv = 1.0 / n_;
      for (int64 ib = 0; ib < (1LL << log2_block); ++ib) {
        for (int64 ip = 0; ip < (1LL << log2_point); ++ip) {
          x[ib * kStrideSize + ip] *= inv;
        }
      }
    } else {
      dft(x);
    }
  }

 private:
  void init();
  void repackWithPadding(Complex*) const;
  void extendWithPadding(Complex*) const;

  int64 paddedIndex(int64 i) const {
    return (i >> kLogBlockSize) * kStrideSize + (i & (kBlockSize - 1));
  }

  void dft(Complex* a) const {
    auto* pw = ws.data();
    int64 l = 1;
    int64 m = n_;
    for (int64 i = 0; i < log4n_; ++i) {
      m /= 4;
      dft4(a, l, m, pw);
      pw += m * 3;
      l *= 4;
    }
    for (int64 i = 0; i < log2n_; ++i) {
      m /= 2;
      dft2(a, l);
      l *= 2;
    }
  }

  void idft(Complex* a) const {
    auto* pw = ws.data() + ws.size();
    int64 l = n_;
    int64 m = 1;
    for (int64 i = 0; i < log2n_; ++i) {
      l /= 2;
      dft2(a, l);
      m *= 2;
    }
    for (int64 i = 0; i < log4n_; ++i) {
      l /= 4;
      pw -= m * 3;
      idft4(a, l, m, pw);
      m *= 4;
    }
  }

  void dft2(Complex* a, const int64 l) const {
    for (int64 j = 0; j < l; ++j) {
      int64 i0 = paddedIndex(2 * j);
      int64 i1 = paddedIndex(2 * j + 1);
      Complex a0 = a[i0];
      Complex a1 = a[i1];
      a[i0] = a0 + a1;
      a[i1] = a0 - a1;
    }
  }

  void dft4(Complex* a, const int64 l, const int64 m, const Complex* ws) const {
    for (int64 j = 0; j < l; ++j) {
      {
        int64 k0 = paddedIndex(4 * j * m);
        int64 k1 = paddedIndex(4 * j * m + m);
        int64 k2 = paddedIndex(4 * j * m + 2 * m);
        int64 k3 = paddedIndex(4 * j * m + 3 * m);
        Complex a0 = a[k0];
        Complex a1 = a[k1];
        Complex a2 = a[k2];
        Complex a3 = a[k3];
        Complex b0 = a0 + a2;
        Complex b1 = a1 + a3;
        Complex b2 = a0 - a2;
        Complex b3 = (a3 - a1).i();
        a[k0] = b0 + b1;
        a[k1] = b0 - b1;
        a[k2] = b2 + b3;
        a[k3] = b2 - b3;
      }
      for (int64 k = 1; k < m; ++k) {
        Complex w1 = ws[k * 3];
        Complex w2 = ws[k * 3 + 1];
        Complex w3 = ws[k * 3 + 2];
        int64 k0 = paddedIndex(4 * j * m + k);
        int64 k1 = paddedIndex(4 * j * m + m + k);
        int64 k2 = paddedIndex(4 * j * m + 2 * m + k);
        int64 k3 = paddedIndex(4 * j * m + 3 * m + k);
        Complex a0 = a[k0];
        Complex a1 = a[k1];
        Complex a2 = a[k2];
        Complex a3 = a[k3];
        Complex b0 = a0 + a2;
        Complex b1 = a1 + a3;
        Complex b2 = a0 - a2;
        Complex b3 = (a3 - a1).i();
        a[k0] = b0 + b1;
        a[k1] = (b0 - b1) * w2;
        a[k2] = (b2 + b3) * w1;
        a[k3] = (b2 - b3) * w3;
      }
    }
  }

  void idft4(Complex* a, const int64 l, const int64 m, const Complex* ws) const {
    for (int64 j = 0; j < l; ++j) {
      {
        int64 k0 = paddedIndex(4 * j * m);
        int64 k1 = paddedIndex(4 * j * m + m);
        int64 k2 = paddedIndex(4 * j * m + 2 * m);
        int64 k3 = paddedIndex(4 * j * m + 3 * m);
        Complex a0 = a[k0];
        Complex a1 = a[k1];
        Complex a2 = a[k2];
        Complex a3 = a[k3];
        Complex b0 = a0 + a1;
        Complex b1 = a0 - a1;
        Complex b2 = a2 + a3;
        Complex b3 = (a2 - a3).i();
        a[k0] = b0 + b2;
        a[k1] = b1 + b3;
        a[k2] = b0 - b2;
        a[k3] = b1 - b3;
      }
      for (int64 k = 1; k < m; ++k) {
        Complex w1 = ws[3 * k].conj();
        Complex w2 = ws[3 * k + 1].conj();
        Complex w3 = ws[3 * k + 2].conj();
        int64 k0 = paddedIndex(4 * j * m + k);
        int64 k1 = paddedIndex(4 * j * m + m + k);
        int64 k2 = paddedIndex(4 * j * m + 2 * m + k);
        int64 k3 = paddedIndex(4 * j * m + 3 * m + k);
        Complex a0 = a[k0];
        Complex a1 = a[k1] * w2;
        Complex a2 = a[k2] * w1;
        Complex a3 = a[k3] * w3;
        Complex b0 = a0 + a1;
        Complex b1 = a0 - a1;
        Complex b2 = a2 + a3;
        Complex b3 = (a2 - a3).i();
        a[k0] = b0 + b2;
        a[k1] = b1 + b3;
        a[k2] = b0 - b2;
        a[k3] = b1 - b3;
      }
    }
  }

  std::vector<Complex> ws;
  std::vector<Complex> work;
  int64 log2_block;
  int64 log2_point;
};

void PMP2::init() {
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

  if (log2_block) {
    work.resize((kStrideSize << log2_block) - kPaddingSize);
  }
}

void PMP2::repackWithPadding(Complex* a) const {
  Complex* x = const_cast<Complex*>(work.data());
  for (int64 ib = 0; ib < (1LL << log2_block); ++ib) {
    for (int64 ip = 0; ip < kBlockSize; ++ip) {
      x[ib * kStrideSize + ip] = a[ib * kBlockSize + ip];
    }
  }
}

void PMP2::extendWithPadding(Complex* a) const {
  for (int64 ib = 0; ib < (1LL << log2_block); ++ib) {
    for (int64 ip = 0; ip < (1LL << log2_point); ++ip) {
      a[ib * kBlockSize + ip] = work[ib * kStrideSize + ip];
    }
  }
}
