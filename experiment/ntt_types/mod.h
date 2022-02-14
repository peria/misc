#pragma once

template <uint64 P>
struct ModTrait {
  static constexpr uint64 kW = 1;   // Primitive root
  static constexpr uint64 kIW = 1;  // W * IW = 1 % P
  static constexpr uint64 kR2 = 1;  // 2^128 % P
  static constexpr uint64 kIP = 1;  // P^(-1) mod 2^64
};

template <>
struct ModTrait<0xffffffff00000001ULL> {
  static constexpr uint64 kW = 7;
  static constexpr uint64 kIW = 0x249249246db6db6eULL;
  static constexpr uint64 kR2 = 0xfffffffe00000001ULL;
  static constexpr uint64 kIP = 0x100000001ULL;
};

template <>
struct ModTrait<0x7fffffffc8000001ULL> {
  static constexpr uint64 kW = 5;
  static constexpr uint64 kIW = 0x4cccccccab333334ULL;
  static constexpr uint64 kR2 = 0x30fffffe40000004ULL;
  static constexpr uint64 kIP = 0x8c40000038000001ULL;
};

template <>
struct ModTrait<0x7fffffff60000001ULL> {
  static constexpr uint64 kW = 0x3ULL;
  static constexpr uint64 kIW = 0x2aaaaaaa75555556ULL;
  static constexpr uint64 kR2 = 0xffffffce0000001ULL;
  static constexpr uint64 kIP = 0xe4000000a0000001ULL;
};

template <>
struct ModTrait<0x7ffffffef0000001ULL> {
  static constexpr uint64 kW = 0x5;
  static constexpr uint64 kIW = 0x33333332c6666667ULL;
  static constexpr uint64 kR2 = 0x40000010ffffffbULL;
  static constexpr uint64 kIP = 0xa100000110000001ULL;
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
