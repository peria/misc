#pragma

#include "ntt.h"

class RefNTT final : public NTT {
  using Type = ElementType;

 public:
  RefNTT(int64 log2n) : NTT(log2n % 2, log2n / 2) {}

  void ntt(ElementType* x, bool backward) const {
    if (backward) {
      dit(x);
      Mod<P> inv = pow(Mod<P>(n_), P - 2);
      for (int64 i = 0; i < n_; ++i) {
        x[i] *= inv;
      }
    } else {
      dif(x);
    }
  }

  void dif(ElementType* x) const {
    int64 m = n_;
    for (int64 i = 0; i < log4n_; ++i) {
      m /= 4;
      dif4(x, m);
    }
    if (log2n_) {
      m /= 2;
      dif2(x, m);
    }
  }

  void dit(ElementType* x) const {
    for (int64 m = 1; m < n_; m *= 2) {
      dit2(x, m);
    }
  }

  static const char* name() { return "Refer"; }

 private:
  void dif2(ElementType* x, int64 m) const {
    Type w = pow(Type(W), (P - 1) / (2 * m));
    for (int64 k = 0; k < m; ++k) {
      Type wk = pow(w, k);
      for (int64 j = 0; j < n_; j += 2 * m) {
        Type& x0 = x[k + j];
        Type& x1 = x[k + j + m];
        Type x1w = (x0 - x1) * wk;
        x0 = x0 + x1;
        x1 = x1w;
      }
    }
  }

  void dif4(ElementType* x, int64 m) const {
    static const Type wq = pow(Type(W), (P - 1) / 4);
    Type w = pow(Type(W), (P - 1) / (4 * m));
    Type w2 = w * w;
    Type w3 = w2 * w;
    Type wk(1);
    Type wk2(1);
    Type wk3(1);
    for (int64 k = 0; k < m; ++k) {
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
      wk *= w;
      wk2 *= w2;
      wk3 *= w3;
    }
  }

  void dit2(ElementType* x, int64 m) const {
    Type w = pow(Type(W), P - 1 - (P - 1) / (2 * m));
    for (int64 k = 0; k < m; ++k) {
      Type wk = pow(w, k);
      for (int64 j = 0; j < n_; j += 2 * m) {
        Type& x0 = x[k + j];
        Type& x1 = x[k + j + m];
        Type x1w = x1 * wk;
        x1 = x0 - x1w;
        x0 = x0 + x1w;
      }
    }
  }
};
