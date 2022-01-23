#pragma

#include "ntt.h"

class RefRadix2 final : public NTT {
  using Type = ElementType;

 public:
  RefRadix2(int64 log2n, int64 log4n = 0) : NTT(log2n, log4n) {}

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
    Type w = pow(Type(W), (P - 1) / n_);
    int64 m = n_;
    for (int64 i = 0; i < logn_; ++i) {
      m /= 2;
      Type wk(1);
      for (int64 k = 0; k < m; ++k) {
        for (int64 j = 0; j < n_; j += 2 * m) {
          Type& x0 = x[k + j];
          Type& x1 = x[k + j + m];
          Type x1w = (x0 - x1) * wk;
          x0 = x0 + x1;
          x1 = x1w;
        }
        wk *= w;
      }
      w *= w;
    }
  }

  void dit(ElementType* x) const {
    int64 m = 1;
    for (int64 i = 0; i < logn_; ++i) {
      Type w = pow(Type(W), P - 1 - (P - 1) / (2 * m));
      Type wk(1);
      for (int64 k = 0; k < m; ++k) {
        for (int64 j = 0; j < n_; j += 2 * m) {
          Type& x0 = x[k + j];
          Type& x1 = x[k + j + m];
          Type x1w = x1 * wk;
          x1 = x0 - x1w;
          x0 = x0 + x1w;
        }
        wk *= w;
      }
      m *= 2;
    }
  }

  static const char* name() { return "ReferR2"; }
};
