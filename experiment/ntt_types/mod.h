#pragma once

template <uint64 P>
struct ModTrait {
  static constexpr uint64 kIP = 1;                // P^(-1) mod 2^64
  static constexpr uint64 kR2 = uint128(-P) % P;  // uint128(-P) % P;
};

template <>
struct ModTrait<0xffffffff00000001ULL> {
  static constexpr uint64 kIP = 0x100000001ULL;
  static constexpr uint64 kR2 = 0xfffffffe00000001ULL;
  static constexpr uint64 kW = 7;
};

#if 0
#include "mod_base.h"
#include "mod_specific.h"
#else
#include "montgomery.h"
#endif

template <uint64 P>
Mod<P> pow(Mod<P> a, uint64 e) {
  Mod<P> r(1);
  if (e & 1) {
    r *= a;
  }
  while (e >>= 1) {
    a *= a;
    if (e & 1) {
      r *= a;
    }
  }
  return r;
}
