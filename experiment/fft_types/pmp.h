#pragma once
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <vector>

using int16 = int16_t;
using int64 = int64_t;
using uint64 = uint64_t;
using uint128 = __uint128_t;

struct Complex {
  Complex() = default;
  Complex(const Complex&) = default;
  Complex(Complex&&) = default;
  Complex(double r, double i) : real(r), imag(i) {}
  Complex& operator=(const Complex&) = default;
  Complex& operator=(Complex&&) = default;
  Complex& operator*=(const Complex& a) {
    double r = this->real * a.real - this->imag * a.imag;
    double i = this->real * a.imag + this->imag * a.real;
    this->real = r;
    this->imag = i;
    return *this;
  }
  Complex& operator*=(double a) {
    real *= a;
    imag *= a;
    return *this;
  }
  Complex i() const { return Complex{-imag, real}; }
  Complex conj() const { return Complex{real, -imag}; }
  static constexpr double kSqrt05 = 0.7071067811865476;
  Complex w8() const {
    return Complex{(real + imag) * kSqrt05, (imag - real) * kSqrt05};
  }

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

class FFT {
 public:
#if defined(USE_RADIX8)
  FFT(int64 logn_)
      : n(1 << logn_),
        logn(logn_),
        log2n(0),
        log4n(logn_ * 2 % 3),
        log8n((logn_ - log4n * 2) / 3) {}
#else
  FFT(int64 logn_)
      : n(1 << logn_),
        logn(logn_),
        log2n(logn_ % 2),
        log4n(logn_ / 2),
        log8n(0) {}
#endif

  void dft(Complex* x, bool backward) {
    if (backward)
      idft(x);
    else
      dft(x);
  }

#if defined(USE_RADIX8)
  void rft(Complex* x, bool backward) {
    if (backward) {
      idft(x);
      double theta = 2 * M_PI / n / 4;
      Complex qw1{std::cos(theta), std::sin(theta)};
      Complex qw2 = qw1 * qw1;
      Complex qw3 = qw2 * qw1;
      auto& ws = GetWs();
      for (int64 i = 0; i < n / 8; ++i) {
        Complex w = ws[i].conj();
        x[4 * i] *= w;
        x[4 * i + 1] *= w * qw1;
        x[4 * i + 2] *= w * qw2;
        x[4 * i + 3] *= w * qw3;
      }
      for (int64 i = n / 8; i < n / 4; ++i) {
        Complex w = ws[i - n / 8].w8().conj();
        x[4 * i] *= w;
        x[4 * i + 1] *= w * qw1;
        x[4 * i + 2] *= w * qw2;
        x[4 * i + 3] *= w * qw3;
      }
    } else {
      double theta = -2 * M_PI / n / 4;
      Complex qw1{std::cos(theta), std::sin(theta)};
      Complex qw2 = qw1 * qw1;
      Complex qw3 = qw2 * qw1;
      auto& ws = GetWs();
      for (int64 i = 0; i < n / 8; ++i) {
        Complex w = ws[i];
        x[4 * i] = x[4 * i] * w;
        x[4 * i + 1] *= w * qw1;
        x[4 * i + 2] *= w * qw2;
        x[4 * i + 3] *= w * qw3;
      }
      for (int64 i = n / 8; i < n / 4; ++i) {
        Complex w = ws[i - n / 8].w8();
        x[4 * i] *= w;
        x[4 * i + 1] *= w * qw1;
        x[4 * i + 2] *= w * qw2;
        x[4 * i + 3] *= w * qw3;
      }
      dft(x);
    }
  }
#else
  void rft(Complex* x, bool backward) {
    if (backward) {
      idft(x);
      double theta = 2 * M_PI / n / 4;
      Complex qw1{std::cos(theta), std::sin(theta)};
      Complex qw2 = qw1 * qw1;
      Complex qw3 = qw2 * qw1;
      auto& ws = GetWs();
      for (int64 i = 0; i < n / 4; ++i) {
        Complex w = ws[i].conj();
        x[4 * i] *= w;
        x[4 * i + 1] *= w * qw1;
        x[4 * i + 2] *= w * qw2;
        x[4 * i + 3] *= w * qw3;
      }
    } else {
      double theta = -2 * M_PI / n / 4;
      Complex qw1{std::cos(theta), std::sin(theta)};
      Complex qw2 = qw1 * qw1;
      Complex qw3 = qw2 * qw1;
      auto& ws = GetWs();
      for (int64 i = 0; i < n / 4; ++i) {
        Complex w = ws[i];
        x[4 * i] *= w;
        x[4 * i + 1] *= w * qw1;
        x[4 * i + 2] *= w * qw2;
        x[4 * i + 3] *= w * qw3;
      }
      dft(x);
    }
  }
#endif

  double getMFlops(double seconds) {
    double flop = 98 * log8n * n / 8 + 34 * log4n * n / 4 + 4 * log2n * n / 2;
    double flops = flop / seconds;
    return flops * 1e-6;
  }

  static void resizeWs(int64 n) { s_ws.resize(n); }

 private:
  const std::vector<Complex>& GetWs();

  void dft(Complex* a) {
    auto& ws = GetWs();
    const Complex* pw = ws.data();
    int64 l = 1;
    int64 m = n;
    for (int64 i = 0; i < log8n; ++i) {
      m /= 8;
      dft8(a, l, m, pw);
      pw += m;
      l *= 8;
    }
    for (int64 i = 0; i < log4n; ++i) {
      m /= 4;
      dft4(a, l, m, pw);
      pw += m;
      l *= 4;
    }
    for (int64 i = 0; i < log2n; ++i) {
      m /= 2;
      dft2(a, l, m, pw);
      pw += m;
      l *= 2;
    }
  }

  void idft(Complex* a) {
    auto pw = GetWs().data() + GetWs().size();
    int64 l = n;
    int64 m = 1;
    for (int64 i = 0; i < log2n; ++i) {
      l /= 2;
      pw -= m;
      idft2(a, l, m, pw);
      m *= 2;
    }
    for (int64 i = 0; i < log4n; ++i) {
      l /= 4;
      pw -= m;
      idft4(a, l, m, pw);
      m *= 4;
    }
    for (int64 i = 0; i < log8n; ++i) {
      l /= 8;
      pw -= m;
      idft8(a, l, m, pw);
      m *= 8;
    }
    double inv = 1.0 / n;
    for (int64 i = 0; i < n; ++i)
      a[i] *= inv;
  }

  void dft2(Complex* a, const int64 l, const int64 m, const Complex* ws) {
    for (int64 j = 0; j < l; ++j) {
      for (int64 k = 0; k < m; ++k) {
        int64 k0 = 2 * j * m + k;
        int64 k1 = k0 + m;
        Complex a0 = a[k0];
        Complex a1 = a[k1];
        a[k0] = a0 + a1;
        a[k1] = a0 - a1;
      }
    }
  }

  void dft4(Complex* a, const int64 l, const int64 m, const Complex* ws) {
    for (int64 j = 0; j < l; ++j) {
      for (int64 k = 0; k < m; ++k) {
        Complex w1 = ws[k];
        Complex w2 = w1 * w1;
        Complex w3 = w2 * w1;
        int64 k0 = 4 * j * m + k;
        int64 k1 = 4 * j * m + m + k;
        int64 k2 = 4 * j * m + 2 * m + k;
        int64 k3 = 4 * j * m + 3 * m + k;
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

  void dft8(Complex* a, const int64 l, const int64 m, const Complex* ws) {
    for (int64 j = 0; j < l; ++j) {
      int64 k0 = 8 * j * m;
      int64 k1 = 8 * j * m + m;
      int64 k2 = 8 * j * m + 2 * m;
      int64 k3 = 8 * j * m + 3 * m;
      int64 k4 = 8 * j * m + 4 * m;
      int64 k5 = 8 * j * m + 5 * m;
      int64 k6 = 8 * j * m + 6 * m;
      int64 k7 = 8 * j * m + 7 * m;
      {
        Complex a0 = a[k0];
        Complex a1 = a[k1];
        Complex a2 = a[k2];
        Complex a3 = a[k3];
        Complex a4 = a[k4];
        Complex a5 = a[k5];
        Complex a6 = a[k6];
        Complex a7 = a[k7];
        Complex b0 = a0 + a4;
        Complex b1 = a1 + a5;
        Complex b2 = a2 + a6;
        Complex b3 = a3 + a7;
        Complex b4 = a0 - a4;
        Complex b5 = (a1 - a5).w8();
        Complex b6 = (a6 - a2).i();
        Complex b7 = (a7 - a3).w8().i();
        Complex c0 = b0 + b2;
        Complex c1 = b1 + b3;
        Complex c2 = b0 - b2;
        Complex c3 = (b3 - b1).i();
        Complex c4 = b4 + b6;
        Complex c5 = b5 + b7;
        Complex c6 = b4 - b6;
        Complex c7 = (b7 - b5).i();
        a[k0] = c0 + c1;
        a[k1] = c0 - c1;
        a[k2] = c2 + c3;
        a[k3] = c2 - c3;
        a[k4] = c4 + c5;
        a[k5] = c4 - c5;
        a[k6] = c6 + c7;
        a[k7] = c6 - c7;
      }
      for (int64 k = 1; k < m; ++k) {
        Complex w1 = ws[k];
        Complex w2 = w1 * w1;
        Complex w3 = w2 * w1;
        Complex w4 = w2 * w2;
        Complex w5 = w3 * w2;
        Complex w6 = w3 * w3;
        Complex w7 = w4 * w3;
        ++k0;
        ++k1;
        ++k2;
        ++k3;
        ++k4;
        ++k5;
        ++k6;
        ++k7;
        Complex a0 = a[k0];
        Complex a1 = a[k1];
        Complex a2 = a[k2];
        Complex a3 = a[k3];
        Complex a4 = a[k4];
        Complex a5 = a[k5];
        Complex a6 = a[k6];
        Complex a7 = a[k7];
        Complex b0 = a0 + a4;
        Complex b1 = a1 + a5;
        Complex b2 = a2 + a6;
        Complex b3 = a3 + a7;
        Complex b4 = a0 - a4;
        Complex b5 = (a1 - a5).w8();
        Complex b6 = (a6 - a2).i();
        Complex b7 = (a7 - a3).w8().i();
        Complex c0 = b0 + b2;
        Complex c1 = b1 + b3;
        Complex c2 = b0 - b2;
        Complex c3 = (b3 - b1).i();
        Complex c4 = b4 + b6;
        Complex c5 = b5 + b7;
        Complex c6 = b4 - b6;
        Complex c7 = (b7 - b5).i();
        a[k0] = c0 + c1;
        a[k1] = (c0 - c1) * w4;
        a[k2] = (c2 + c3) * w2;
        a[k3] = (c2 - c3) * w6;
        a[k4] = (c4 + c5) * w1;
        a[k5] = (c4 - c5) * w5;
        a[k6] = (c6 + c7) * w3;
        a[k7] = (c6 - c7) * w7;
      }
    }
  }

  void idft2(Complex* a, const int64 l, const int64 m, const Complex* ws) {
    for (int64 j = 0; j < l; ++j) {
      for (int64 k = 0; k < m; ++k) {
        int64 k0 = 2 * j * m + k;
        int64 k1 = k0 + m;
        Complex a0 = a[k0];
        Complex a1 = a[k1];
        a[k0] = a0 + a1;
        a[k1] = a0 - a1;
      }
    }
  }

  void idft4(Complex* a, const int64 l, const int64 m, const Complex* ws) {
    for (int64 j = 0; j < l; ++j) {
      for (int64 k = 0; k < m; ++k) {
        Complex w1 = ws[k].conj();
        Complex w2 = w1 * w1;
        Complex w3 = w2 * w1;
        int64 k0 = 4 * j * m + k;
        int64 k1 = 4 * j * m + m + k;
        int64 k2 = 4 * j * m + 2 * m + k;
        int64 k3 = 4 * j * m + 3 * m + k;
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

  void idft8(Complex* a, const int64 l, const int64 m, const Complex* ws) {
    for (int64 j = 0; j < l; ++j) {
      int64 k0 = 8 * j * m;
      int64 k1 = 8 * j * m + m;
      int64 k2 = 8 * j * m + 2 * m;
      int64 k3 = 8 * j * m + 3 * m;
      int64 k4 = 8 * j * m + 4 * m;
      int64 k5 = 8 * j * m + 5 * m;
      int64 k6 = 8 * j * m + 6 * m;
      int64 k7 = 8 * j * m + 7 * m;
      {
        Complex a0 = a[k0];
        Complex a1 = a[k1];
        Complex a2 = a[k2];
        Complex a3 = a[k3];
        Complex a4 = a[k4];
        Complex a5 = a[k5];
        Complex a6 = a[k6];
        Complex a7 = a[k7];
        Complex b0 = a0 + a1;
        Complex b1 = a0 - a1;
        Complex b2 = a2 + a3;
        Complex b3 = (a2 - a3).i();
        Complex b4 = a4 + a5;
        Complex b5 = a4 - a5;
        Complex b6 = a6 + a7;
        Complex b7 = (a6 - a7).i();
        Complex c0 = b0 + b2;
        Complex c1 = b1 + b3;
        Complex c2 = b0 - b2;
        Complex c3 = b1 - b3;
        Complex c4 = b4 + b6;
        Complex c5 = (b5 + b7).w8().i();
        Complex c6 = (b4 - b6).i();
        Complex c7 = (b7 - b5).w8();
        a[k0] = c0 + c4;
        a[k1] = c1 + c5;
        a[k2] = c2 + c6;
        a[k3] = c3 + c7;
        a[k4] = c0 - c4;
        a[k5] = c1 - c5;
        a[k6] = c2 - c6;
        a[k7] = c3 - c7;
      }
      for (int64 k = 1; k < m; ++k) {
        Complex w1 = ws[k].conj();
        Complex w2 = w1 * w1;
        Complex w3 = w2 * w1;
        Complex w4 = w2 * w2;
        Complex w5 = w3 * w2;
        Complex w6 = w3 * w3;
        Complex w7 = w4 * w3;
        ++k0;
        ++k1;
        ++k2;
        ++k3;
        ++k4;
        ++k5;
        ++k6;
        ++k7;
        Complex a0 = a[k0];
        Complex a1 = a[k1] * w4;
        Complex a2 = a[k2] * w2;
        Complex a3 = a[k3] * w6;
        Complex a4 = a[k4] * w1;
        Complex a5 = a[k5] * w5;
        Complex a6 = a[k6] * w3;
        Complex a7 = a[k7] * w7;
        Complex b0 = a0 + a1;
        Complex b1 = a0 - a1;
        Complex b2 = a2 + a3;
        Complex b3 = (a2 - a3).i();
        Complex b4 = a4 + a5;
        Complex b5 = a4 - a5;
        Complex b6 = a6 + a7;
        Complex b7 = (a6 - a7).i();
        Complex c0 = b0 + b2;
        Complex c1 = b1 + b3;
        Complex c2 = b0 - b2;
        Complex c3 = b1 - b3;
        Complex c4 = b4 + b6;
        Complex c5 = (b5 + b7).w8().i();
        Complex c6 = (b4 - b6).i();
        Complex c7 = (b7 - b5).w8();
        a[k0] = c0 + c4;
        a[k1] = c1 + c5;
        a[k2] = c2 + c6;
        a[k3] = c3 + c7;
        a[k4] = c0 - c4;
        a[k5] = c1 - c5;
        a[k6] = c2 - c6;
        a[k7] = c3 - c7;
      }
    }
  }

  static std::vector<std::vector<Complex>> s_ws;

  const int64 n;
  const int64 logn;
  const int64 log2n;
  const int64 log4n;
  const int64 log8n;
  const int64 logn_all = 0;
};

const std::vector<Complex>& FFT::GetWs() {
  std::vector<Complex>& ws = s_ws[logn];
  if (ws.size())
    return ws;

  int64 l = 1;
  int64 m = 1 << logn;
  ws.reserve(m);
  const double theta = -2 * M_PI / n;
  for (int64 i = 0; i < log8n; ++i) {
    m /= 8;
    for (int64 k = 0; k < m; ++k) {
      double t = theta * l * k;
      Complex w{std::cos(t), std::sin(t)};
      ws.push_back(w);
    }
    l *= 8;
  }
  for (int64 i = 0; i < log4n; ++i) {
    m /= 4;
    for (int64 k = 0; k < m; ++k) {
      double t = theta * l * k;
      Complex w{std::cos(t), std::sin(t)};
      ws.push_back(w);
    }
    l *= 4;
  }
  if (log2n == 1) {
    m /= 2;
    for (int64 k = 0; k < m; ++k) {
      const double t = theta * l * k;
      Complex w{std::cos(t), std::sin(t)};
      ws.push_back(w);
    }
    l *= 2;
  }
  return ws;
}
