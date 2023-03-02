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
  for (int i = 0; i < logn_; ++i) {
    m /= 2;
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
    l *= 2;
  }
}

void DIT::IDft(Complex* x) const {
  int m = 1;
  int l = n_;
  for (int i = 0; i < logn_; ++i) {
    l /= 2;
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
    m *= 2;
  }

  double inverse = 1.0 / n_;
  for (int i = 0; i < n_; ++i) {
    x[i] *= inverse;
  }
}