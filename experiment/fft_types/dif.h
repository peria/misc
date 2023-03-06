#include "fmt.h"

#include <cmath>

class DIF : public FMT {
 public:
  DIF(int logn) : FMT(logn) { Init(); }
  double GetMemory() const override;

  static const char* name() { return "DIF"; }

 private:
  void Rft(Complex*) const override;
  void IRft(Complex*) const override;

  void Init();
  void Dft(Complex*) const;
  void IDft(Complex*) const;

  void Dft2(Complex*, const int, const int, const Complex*) const;
  void IDft2(Complex*, const int, const int, const Complex*) const;
  void Dft4(Complex*, const int, const int, const Complex*) const;
  void IDft4(Complex*, const int, const int, const Complex*) const;

  std::vector<Complex> ws_;
  Complex qw_;
};

double DIF::GetMemory() const {
  return (ws_.size() + 1) * sizeof(Complex) / 1024.0 / 1024.0;
}

void DIF::Init() {
  const double theta = 2 * M_PI / n_;
  int l = 1;
  int m = n_;
  for (int i = 0; i < logn_; ++i) {
    m /= 2;
    for (int k = 0; k < m; ++k) {
      double t = theta * k * l;
      Complex w1(std::cos(t), std::sin(t));
      ws_.push_back(w1);
    }
    l *= 2;
  }

  qw_ = Complex(std::cos(theta / 4), std::sin(theta / 4));
}

void DIF::Rft(Complex* x) const {
  const double theta = 2 * M_PI / (4 * n_);
  const Complex& qw1 = qw_;
  const Complex qw2 = qw1 * qw1;
  const Complex qw3 = qw2 * qw1;
  for (int i = 0; i < n_ / 4; ++i) {
    const Complex& w0 = ws_[i];
    x[4 * i + 0] *= w0;
    x[4 * i + 1] *= w0 * qw1;
    x[4 * i + 2] *= w0 * qw2;
    x[4 * i + 3] *= w0 * qw3;
  }
  Dft(x);
}

void DIF::IRft(Complex* x) const {
  const double theta = 2 * M_PI / (4 * n_);
  const Complex qw1 = qw_.conj();
  const Complex qw2 = qw1 * qw1;
  const Complex qw3 = qw2 * qw1;
  IDft(x);
  for (int i = 0; i < n_ / 4; ++i) {
    const Complex w0 = ws_[i].conj();
    x[4 * i + 0] *= w0;
    x[4 * i + 1] *= w0 * qw1;
    x[4 * i + 2] *= w0 * qw2;
    x[4 * i + 3] *= w0 * qw3;
  }
}

void DIF::Dft(Complex* x) const {
  const Complex* wp = ws_.data();
  int m = n_;
  int l = 1;
  for (int i = 0; i < log4n_; ++i) {
    m /= 4;
    Dft4(x, m, l, wp);
    wp += 3 * m;
    l *= 4;
  }
  for (int i = 0; i < log2n_; ++i) {
    m /= 2;
    Dft2(x, m, l, wp);
    wp += m;
    l *= 2;
  }
}

void DIF::IDft(Complex* x) const {
  const Complex* wp = ws_.data() + ws_.size();
  int m = 1;
  int l = n_;
  for (int i = 0; i < log2n_; ++i) {
    l /= 2;
    wp -= m;
    Dft2(x, m, l, wp);
    m *= 2;
  }
  for (int i = 0; i < log4n_; ++i) {
    l /= 4;
    wp -= 3 * m;
    IDft4(x, m, l, wp);
    m *= 4;
  }

  double inverse = 1.0 / n_;
  for (int i = 0; i < n_; ++i) {
    x[i] *= inverse;
  }
}

void DIF::Dft2(Complex* x, const int m, const int l, const Complex* wp) const {
  for (int j = 0; j < l; ++j) {
    const Complex& w1 = wp[0];
    int k0 = 2 * j;
    int k1 = 2 * j + 1;
    Complex x0 = x[k0];
    Complex x1 = x[k1];
    x[k0] = x0 + x1;
    x[k1] = x0 - x1;
  }
}

void DIF::Dft4(Complex* x, const int m, const int l, const Complex* wp) const {
  const Complex* wp2 = wp + 2 * m;
  for (int j = 0; j < l; ++j) {
    for (int k = 0; k < m; ++k) {
      const Complex& w1 = wp[k];
      const Complex& w2 = wp2[k];
      const Complex& w3 = w1 * w2;
      int k0 = 4 * j * m + k;
      int k1 = 4 * j * m + m + k;
      int k2 = 4 * j * m + 2 * m + k;
      int k3 = 4 * j * m + 3 * m + k;
      const Complex& x0 = x[k0];
      const Complex& x1 = x[k1];
      const Complex& x2 = x[k2];
      const Complex& x3 = x[k3];
      const Complex y0 = x0 + x2;
      const Complex y1 = x1 + x3;
      const Complex y2 = x0 - x2;
      const Complex y3 = (x1 - x3).i();
      x[k0] = y0 + y1;
      x[k1] = (y0 - y1) * w2;
      x[k2] = (y2 + y3) * w1;
      x[k3] = (y2 - y3) * w3;
    }
  }
}

void DIF::IDft4(Complex* x, const int m, const int l, const Complex* wp) const {
  const Complex* wp2 = wp + 2 * m;
  for (int j = 0; j < l; ++j) {
    for (int k = 0; k < m; ++k) {
      const Complex w1 = wp[k].conj();
      const Complex w2 = wp2[k].conj();
      const Complex w3 = w1 * w2;
      int k0 = 4 * j * m + k;
      int k1 = 4 * j * m + m + k;
      int k2 = 4 * j * m + 2 * m + k;
      int k3 = 4 * j * m + 3 * m + k;
      Complex x0 = x[k0];
      Complex x1 = x[k1] * w2;
      Complex x2 = x[k2] * w1;
      Complex x3 = x[k3] * w3;
      Complex y0 = x0 + x1;
      Complex y1 = x0 - x1;
      Complex y2 = x2 + x3;
      Complex y3 = (x3 - x2).i();
      x[k0] = y0 + y2;
      x[k1] = y1 + y3;
      x[k2] = y0 - y2;
      x[k3] = y1 - y3;
    }
  }
}
