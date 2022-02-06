#pragma

#include "ntt.h"

class PNT final : public NTT {
  using Type = ElementType;

 public:
  static const char* name() { return "PNT"; }

  PNT(int64 log2n) : NTT(log2n % 2, log2n / 2) { init(); }

  void ntt(ElementType* x, bool backward) const {
    if (backward) {
      intt(x);
      Mod<P> inv = pow(Mod<P>(n_), P - 2);
      for (int64 i = 0; i < n_; ++i) {
        x[i] *= inv;
      }
    } else {
      ntt(x);
    }
  }

 private:
  void ntt(ElementType* x) const {
    const Type* wkp = ws_.data();
    int64 m = n_;
    for (int64 i = 0; i < log4n_; ++i) {
      m /= 4;
      ntt4(x, m, wkp);
      wkp += 3 * (m - 1);
    }
    if (log2n_) {
      m /= 2;
      ntt2(x, m);
    }
  }

  void intt(ElementType* x) const {
    const Type* wkp = iws_.data();
    int64 m = 1;
    if (log2n_) {
      intt2(x, m);
      m *= 2;
    }
    for (int64 i = 0; i < log4n_; ++i) {
      intt4(x, m, wkp);
      wkp += 3 * (m - 1);
      m *= 4;
    }
  }

  void ntt2(ElementType* x, int64 m) const {
    for (int64 k = 0; k < m; ++k) {
      for (int64 j = 0; j < n_; j += 2 * m) {
        Type& x0 = x[k + j];
        Type& x1 = x[k + j + m];
        Type x1w = x0 - x1;
        x0 = x0 + x1;
        x1 = x1w;
      }
    }
  }

  void ntt4(Type* x, int64 m, const Type* wkp) const {
    static const Type wq = pow(Type(W), (P - 1) / 4);
    for (int64 j = 0; j < n_; j += 4 * m) {
      Type& x0 = x[j];
      Type& x1 = x[j + m];
      Type& x2 = x[j + 2 * m];
      Type& x3 = x[j + 3 * m];
      Type y0 = x0 + x2;
      Type y1 = x1 + x3;
      Type y2 = x0 - x2;
      Type y3 = (x1 - x3) * wq;
      x0 = y0 + y1;
      x2 = y2 + y3;
      x1 = y0 - y1;
      x3 = y2 - y3;
    }
    for (int64 k = 1; k < m; ++k) {
      const Type& wk = *wkp++;
      const Type& wk2 = *wkp++;
      const Type& wk3 = *wkp++;
      for (int64 j = 0; j < n_; j += 4 * m) {
        Type& x0 = x[k + j];
        Type& x1 = x[k + j + m];
        Type& x2 = x[k + j + 2 * m];
        Type& x3 = x[k + j + 3 * m];
        Type y0 = x0 + x2;
        Type y1 = x1 + x3;
        Type y2 = x0 - x2;
        Type y3 = (x1 - x3) * wq;
        x0 = y0 + y1;
        x2 = (y2 + y3) * wk;
        x1 = (y0 - y1) * wk2;
        x3 = (y2 - y3) * wk3;
      }
    }
  }

  void intt2(ElementType* x, int64 m) const {
    for (int64 k = 0; k < m; ++k) {
      for (int64 j = 0; j < n_; j += 2 * m) {
        Type& x0 = x[k + j];
        Type& x1 = x[k + j + m];
        Type x1w = x1;
        x1 = x0 - x1w;
        x0 = x0 + x1w;
      }
    }
  }

  void intt4(Type* x, int64 m, const Type* wkp) const {
    static const Type wq = pow(Type(W), (P - 1) / 4 * 3);
    for (int64 j = 0; j < n_; j += 4 * m) {
      Type& x0 = x[j];
      Type& x1 = x[j + m];
      Type& x2 = x[j + 2 * m];
      Type& x3 = x[j + 3 * m];
      Type y0 = x0 + x1;
      Type y1 = x0 - x1;
      Type y2 = x2 + x3;
      Type y3 = (x2 - x3) * wq;
      x0 = y0 + y2;
      x1 = y1 + y3;
      x2 = y0 - y2;
      x3 = y1 - y3;
    }
    for (int64 k = 1; k < m; ++k) {
      const Type& wk = *wkp++;
      const Type& wk2 = *wkp++;
      const Type& wk3 = *wkp++;
      for (int64 j = 0; j < n_; j += 4 * m) {
        Type& x0 = x[k + j];
        Type& x1 = x[k + j + m];
        Type& x2 = x[k + j + 2 * m];
        Type& x3 = x[k + j + 3 * m];
        Type x1w = x1 * wk2;
        Type x2w = x2 * wk;
        Type x3w = x3 * wk3;
        Type y0 = x0 + x1w;
        Type y1 = x0 - x1w;
        Type y2 = x2w + x3w;
        Type y3 = (x2w - x3w) * wq;
        x0 = y0 + y2;
        x1 = y1 + y3;
        x2 = y0 - y2;
        x3 = y1 - y3;
      }
    }
  }

  void init() {
    int64 m = n_;
    for (int64 i = 0; i < log4n_; ++i) {
      m /= 4;
      Type w = pow(Type(W), (P - 1) / (4 * m));
      Type w2 = w * w;
      Type w3 = w2 * w;
      Type wk(w);
      Type wk2(w2);
      Type wk3(w3);
      for (int64 k = 1; k < m; ++k) {
        ws_.push_back(wk);
        ws_.push_back(wk2);
        ws_.push_back(wk3);
        wk *= w;
        wk2 *= w2;
        wk3 *= w3;
      }
    }

    m = 1;
    if (log2n_) {
      m *= 2;
    }
    for (int64 i = 0; i < log4n_; ++i) {
      Type w = pow(Type(W), P - 1 - (P - 1) / (4 * m));
      Type w2 = w * w;
      Type w3 = w2 * w;
      Type wk(w);
      Type wk2(w2);
      Type wk3(w3);
      for (int64 k = 1; k < m; ++k) {
        iws_.push_back(wk);
        iws_.push_back(wk2);
        iws_.push_back(wk3);
        wk *= w;
        wk2 *= w2;
        wk3 *= w3;
      }
      m *= 4;
    }
  }

  std::vector<Mod<P>> ws_;
  std::vector<Mod<P>> iws_;
};
