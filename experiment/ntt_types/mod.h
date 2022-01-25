#pragma once

#include "mod_base.h"
// #include "mod_specific.h"
// #include "montgomery.h"

template <uint64 P>
Mod<P> pow(Mod<P> a, uint64 e) {
  Mod<P> r(1);
  if (e & 1) {
    r *= a;
  }
  e >>= 1;
  while (e) {
    a *= a;
    if (e & 1) {
      r *= a;
    }
    e >>= 1;
  }
  return r;
}
