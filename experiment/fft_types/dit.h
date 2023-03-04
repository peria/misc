#include "fmt.h"

#include <cmath>

class DIT : public FMT {
 public:
  DIT(int logn) : FMT(logn) { Init(); }

  static const char* name() { return "DIT"; }

 private:
  void Rft(Complex*) const override;
  void IRft(Complex*) const override;

  void Init();
  void Dft(Complex*) const;
  void IDft(Complex*) const;

  void Dft2(Complex*, const int, const int) const;
  void IDft2(Complex*, const int, const int) const;
  void Dft4(Complex*, const int, const int) const;
  void IDft4(Complex*, const int, const int) const;

  std::vector<Complex> rws_;
  std::vector<Complex> ws_;
  Complex qw_;
};

void DIT::Init() {
  const double theta = 2 * M_PI / n_;
  for (int i = 0; i < n_ / 4; ++i) {
    double t = theta * i;
    Complex w(std::cos(t), std::sin(t));
    rws_.push_back(w);
  }

  for (int j = 0, r = 0; j < n_ / 2; ++j) {
    double t = theta * r;
    Complex w1(std::cos(t), std::sin(t));
    ws_.push_back(w1);
    for (int b = n_ / 4; b; b >>= 1) {
      r ^= b;
      if (r & b) {
        break;
      }
    }
  }

  qw_ = Complex(std::cos(theta / 4), std::sin(theta / 4));
}

void DIT::Rft(Complex* x) const {
  const double theta = 2 * M_PI / (4 * n_);
  const Complex& qw1 = qw_;
  const Complex qw2 = qw1 * qw1;
  const Complex qw3 = qw2 * qw1;
  for (int i = 0; i < n_ / 4; ++i) {
    const Complex& w0 = rws_[i];
    x[4 * i + 0] *= w0;
    x[4 * i + 1] *= w0 * qw1;
    x[4 * i + 2] *= w0 * qw2;
    x[4 * i + 3] *= w0 * qw3;
  }
  Dft(x);
}

void DIT::IRft(Complex* x) const {
  const double theta = 2 * M_PI / (4 * n_);
  const Complex qw1 = qw_.conj();
  const Complex qw2 = qw1 * qw1;
  const Complex qw3 = qw2 * qw1;
  IDft(x);
  for (int i = 0; i < n_ / 4; ++i) {
    const Complex w0 = rws_[i].conj();
    x[4 * i + 0] *= w0;
    x[4 * i + 1] *= w0 * qw1;
    x[4 * i + 2] *= w0 * qw2;
    x[4 * i + 3] *= w0 * qw3;
  }
}

void DIT::Dft(Complex* x) const {
  int m = n_;
  int l = 1;
  for (int i = 0; i < log2n_; ++i) {
    m /= 2;
    Dft2(x, m, l);
    l *= 2;
  }
  for (int i = 0; i < log4n_; ++i) {
    m /= 4;
    Dft4(x, m, l);
    l *= 4;
  }
}

void DIT::IDft(Complex* x) const {
  int m = 1;
  int l = n_;
  for (int i = 0; i < log4n_; ++i) {
    l /= 4;
    IDft4(x, m, l);
    m *= 4;
  }
  for (int i = 0; i < log2n_; ++i) {
    l /= 2;
    IDft2(x, m, l);
    m *= 2;
  }

  double inverse = 1.0 / n_;
  for (int i = 0; i < n_; ++i) {
    x[i] *= inverse;
  }
}

void DIT::Dft2(Complex* x, const int m, const int l) const {
  for (int j = 0; j < l; ++j) {
    const Complex& w1 = ws_[j];
    for (int k = 0; k < m; ++k) {
      int k0 = 2 * j * m + k;
      int k1 = 2 * j * m + m + k;
      Complex x0 = x[k0];
      Complex x1 = x[k1] * w1;
      x[k0] = x0 + x1;
      x[k1] = x0 - x1;
    }
  }
}

void DIT::IDft2(Complex* x, const int m, const int l) const {
  for (int j = 0; j < l; ++j) {
    const Complex w1 = ws_[j].conj();
    for (int k = 0; k < m; ++k) {
      int k0 = 2 * j * m + k;
      int k1 = 2 * j * m + m + k;
      Complex x0 = x[k0];
      Complex x1 = x[k1];
      x[k0] = x0 + x1;
      x[k1] = (x0 - x1) * w1;
    }
  }
}

void DIT::Dft4(Complex* x, const int m, const int l) const {
  for (int j = 0; j < l; ++j) {
    const Complex& w1 = ws_[2 * j];
    const Complex& w2 = ws_[j];
    const Complex w3 = w2 * w1;
    for (int k = 0; k < m; ++k) {
      int k0 = 4 * j * m + k;
      int k1 = 4 * j * m + m + k;
      int k2 = 4 * j * m + 2 * m + k;
      int k3 = 4 * j * m + 3 * m + k;
      const Complex& x0 = x[k0];
      const Complex x1 = x[k1] * w1;
      const Complex x2 = x[k2] * w2;
      const Complex x3 = x[k3] * w3;
      const Complex y0 = x0 + x2;
      const Complex y1 = x1 + x3;
      const Complex y2 = x0 - x2;
      const Complex y3 = (x1 - x3).i();
      x[k0] = y0 + y1;
      x[k1] = y0 - y1;
      x[k2] = y2 + y3;
      x[k3] = y2 - y3;
    }
  }
}

void DIT::IDft4(Complex* x, const int m, const int l) const {
  for (int j = 0; j < l; ++j) {
    const Complex w1 = ws_[2 * j].conj();
    const Complex w2 = ws_[j].conj();
    const Complex w3 = w1 * w2;
    for (int k = 0; k < m; ++k) {
      int k0 = 4 * j * m + k;
      int k1 = 4 * j * m + m + k;
      int k2 = 4 * j * m + 2 * m + k;
      int k3 = 4 * j * m + 3 * m + k;
      const Complex& x0 = x[k0];
      const Complex& x1 = x[k1];
      const Complex& x2 = x[k2];
      const Complex& x3 = x[k3];
      const Complex y0 = x0 + x1;
      const Complex y1 = x0 - x1;
      const Complex y2 = x2 + x3;
      const Complex y3 = (x3 - x2).i();
      x[k0] = y0 + y2;
      x[k1] = (y1 + y3) * w1;
      x[k2] = (y0 - y2) * w2;
      x[k3] = (y1 - y3) * w3;
    }
  }
}
