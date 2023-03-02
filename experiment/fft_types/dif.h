#include "fmt.h"

#include <cmath>

class DIF : public FMT {
 public:
  DIF(int logn) : FMT(logn) { Init(); }

  static const char* name() { return "DIF"; }

 private:
  void Rft(Complex*) const override;
  void IRft(Complex*) const override;

  void Init();
  void Dft(Complex*) const;
  void IDft(Complex*) const;

  std::vector<Complex> ws_;
  Complex qw_;
};

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
  for (int i = 0; i < logn_; ++i) {
    m /= 2;
    for (int j = 0; j < l; ++j) {
      for (int k = 0; k < m; ++k) {
        const Complex& w1 = wp[k];
        int k0 = 2 * j * m + k;
        int k1 = 2 * j * m + m + k;
        Complex x0 = x[k0];
        Complex x1 = x[k1];
        x[k0] = x0 + x1;
        x[k1] = (x0 - x1) * w1;
      }
    }
    wp += m;
    l *= 2;
  }
}

void DIF::IDft(Complex* x) const {
  const double theta = 2 * M_PI / n_;

  const Complex* wp = ws_.data() + ws_.size();
  int m = 1;
  int l = n_;
  for (int i = 0; i < logn_; ++i) {
    l /= 2;
    wp -= m;
    for (int j = 0; j < l; ++j) {
      for (int k = 0; k < m; ++k) {
        const Complex& w1 = wp[k];
        int k0 = 2 * j * m + k;
        int k1 = 2 * j * m + m + k;
        Complex x0 = x[k0];
        Complex x1 = x[k1] * w1.conj();
        x[k0] = x0 + x1;
        x[k1] = x0 - x1;
      }
    }
    m *= 2;
  }

  double inverse = 1.0 / n_;
  for (int i = 0; i < n_; ++i) {
    x[i] *= inverse;
  }
}