#include "fmt.h"

#include <cmath>

class DDIT : public FMT {
 public:
  DDIT(int logn) : FMT(logn) { Init(); }
  double GetMemory() const override;

  static const char* name() { return "DDIT"; }

 private:
  void Rft(Complex*) const override;
  void IRft(Complex*) const override;

  void Init();
  void Dft(Complex*) const;
  void IDft(Complex*) const;

  void Dft2(Complex*, const int, const int, const int, const int) const;
  void IDft2(Complex*, const int, const int, const int, const int) const;
  void Dft4(Complex*, const int, const int, const int, const int) const;
  void IDft4(Complex*, const int, const int, const int, const int) const;

  std::vector<Complex> rws_;
  std::vector<Complex> ws_;
  Complex qw_;
};

double DDIT::GetMemory() const {
  return (rws_.size() + ws_.size() + 1) * sizeof(Complex) / 1024.0 / 1024.0;
}

void DDIT::Init() {
  {
    double theta = 2 * M_PI / 8;
    ws_.push_back(Complex(1, 0));
    ws_.push_back(Complex(std::cos(theta), std::sin(theta)));
    for (int i = 1; i < logn_ - 2; ++i) {
      int m = 1 << i;
      theta *= 0.5;
      Complex w(std::cos(theta), std::sin(theta));
      for (int j = 0; j < m; ++j) {
        ws_.push_back(w * ws_[j]);
      }
    }
  }

  {
    const double theta = 2 * M_PI / n_;
    for (int i = 0; i < n_ / 4; ++i) {
      double t = theta * i;
      Complex w(std::cos(t), std::sin(t));
      rws_.push_back(w);
    }
    const double qtheta = theta / 4;
    qw_ = Complex(std::cos(qtheta), std::sin(qtheta));
  }
}

void DDIT::Rft(Complex* x) const {
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

void DDIT::IRft(Complex* x) const {
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

void DDIT::Dft(Complex* x) const {
  static constexpr int kCacheSize = 1 << 18;
  int m = n_;
  int l = 1;
  for (int i = 0; i < log2n_; ++i) {
    m /= 2;
    Dft2(x, m, l, 0, l);
    l *= 2;
  }
  if (n_ <= kCacheSize) {
    for (int i = 0; i < log4n_; ++i) {
      m /= 4;
      Dft4(x, m, l, 0, l);
      l *= 4;
    }
  } else {
    int i = 0;
    for (i = 0; i < log4n_ && m > kCacheSize; ++i) {
      m /= 4;
      Dft4(x, m, l, 0, l);
      l *= 4;
    }
    int n_block = n_ / kCacheSize;
    for (int j = 0; j < n_block; ++j) {
      int m_block = m;
      int l_block = l;
      for (int i_block = i, j_block = l_block / n_block; i_block < log4n_;
           ++i_block) {
        int j0 = j * j_block;
        int j1 = (j + 1) * j_block;
        m_block /= 4;
        Dft4(x, m_block, l_block, j0, j1);
        l_block *= 4;
        j_block *= 4;
      }
    }
  }
}

void DDIT::IDft(Complex* x) const {
  int m = 1;
  int l = n_;
  for (int i = 0; i < log4n_; ++i) {
    l /= 4;
    IDft4(x, m, l, 0, l);
    m *= 4;
  }
  for (int i = 0; i < log2n_; ++i) {
    l /= 2;
    Dft2(x, m, l, 0, l);
    m *= 2;
  }

  double inverse = 1.0 / n_;
  for (int i = 0; i < n_; ++i) {
    x[i] *= inverse;
  }
}

void DDIT::Dft2(Complex* x,
                const int m,
                const int l,
                const int,
                const int) const {
  for (int k = 0; k < m; ++k) {
    int k0 = k;
    int k1 = m + k;
    Complex x1 = x[k0] - x[k1];
    x[k0] += x[k1];
    x[k1] = x1;
  }
}

void DDIT::Dft4(Complex* x,
                const int m,
                const int l,
                const int j0,
                const int jn) const {
  for (int j = j0; j < jn; ++j) {
    const Complex& w1 = ws_[j];
    const Complex w2 = w1 * w1;
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

void DDIT::IDft4(Complex* x,
                 const int m,
                 const int l,
                 const int j0,
                 const int jn) const {
  for (int j = j0; j < jn; ++j) {
    const Complex w1 = ws_[j].conj();
    const Complex w2 = w1 * w1;
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
