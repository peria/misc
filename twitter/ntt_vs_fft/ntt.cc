#include "ntt.h"

#include <cstring>
#include <iostream>

using namespace std;

namespace {

uint64* g_work = nullptr;
uint64* g_table = nullptr;
uint64* g_inv_table = nullptr;

}

namespace Mod {

const uint64 kLowerMask = (1ULL << 32) - 1;
const uint64 P = 0xffffffff00000001ULL; // 2^64-2^32+1
const uint64 kRoot = 0x185629DCDA58878CULL; // 7^((P-1)/2^32) % P

inline uint64 sub(uint64 a, uint64 b) {
  uint64 c = a - b;
  if (a < b)
    c += P;
  return c;
}

inline uint64 add(uint64 a, uint64 b) {
  return sub(a, P - b);
}

inline uint64 reduce(uint64 a1, uint64 a0) {
  uint64 a10 = a1 & kLowerMask, a11 = a1 >> 32;
  if (a0 >= P)
    a0 -= P;
  a0 = sub(a0, a10);
  a0 = sub(a0, a11);
  return add(a0, a10 << 32);
}

inline uint64 mul(uint64 a, uint64 b) {
  // 64bit * 64bit -> 128bit
#if defined(UINT128)
  uint128 d = static_cast<uint128>(a) * static_cast<uint128>(b);
  uint64 d1 = d >> 64;
  uint64 d0 = d;
#else
  uint64 a0 = a & kLowerMask, a1 = a >> 32;
  uint64 b0 = b & kLowerMask, b1 = b >> 32;
  uint64 c00 = a0 * b0, c11 = a1 * b1;
  uint64 c01 = a0 * b1, c10 = a1 * b0;
  c01 += c10;
  if (c01 < c10)
    c11 += 1ULL << 32;
  uint64 d0 = c00 + (c01 << 32);
  if (d0 < c00)
    ++c11;
  uint64 d1 = c11 + (c01 >> 32);
#endif

  // [d1, d0] mod P -> e
  return reduce(d1, d0);
}

inline uint64 pow(uint64 a, uint64 e) {
  uint64 res = (e & 1) ? a : 1;
  for (; e >>= 1;) {
    a = mul(a, a);
    if (e & 1)
      res = mul(res, a);
  }
  return res;
}

inline uint64 inv(uint64 a) {
  return pow(a, P - 2);
}

}  // namespace Mod

void NTT::SetGlobals(void* work, void* table, void* inv_table) {
  g_work = static_cast<uint64*>(work);
  g_table = static_cast<uint64*>(table);
  g_inv_table = static_cast<uint64*>(inv_table);
}

bool NTT::Validate(uint64* data) {
  for (int log2n = 1; log2n < 10; ++log2n) {
    const int n = 1 << log2n;
    InitTable(log2n);
    for (int i = 0; i < n; ++i)
      data[i] = i;
    Forward(log2n, n, data);
    Backward(log2n, n, data);
    for (int i = 0; i < n; ++i) {
      if (data[i] != i) {
        cout << "2^" << log2n << " [" << i << "] " << data[i] << "\n";
        return false;
      }
    }
  }
    return true;
}

void NTT::Square(uint64* data, int n) {
  for (int i = 0; i < n; ++i)
    data[i] = Mod::mul(data[i], data[i]);
}

void NTT::InitTable(int log2n) {
  uint64* ptr0 = g_table;
  uint64* ptr1 = g_inv_table;

  int width = 1;
  for (int k = 0; k < log2n; ++k) {
    uint64 w = Mod::pow(Mod::kRoot, (1ULL << 31) / width);
    uint64 m = Mod::inv(w);
    uint64 wi = 1, mi = 1;
    for (int i = 0; i < width; ++i) {
      ptr0[i] = wi;
      ptr1[i] = mi;
      wi = Mod::mul(wi, w);
      mi = Mod::mul(mi, m);
    }
    ptr0 += width;
    ptr1 += width;
    width *= 2;
  }
}

void NTT::Forward(int log2n, int n, uint64* data) {
  Core(log2n, n, g_table, g_work, data);
}

void NTT::Backward(int log2n, int n, uint64* data) {
  Core(log2n, n, g_inv_table, g_work, data);

  uint64 inv = Mod::inv(n);
  for (int i = 0; i < n; ++i)
    data[i] = Mod::mul(data[i], inv);
}

void NTT::Core(int log2n, int n, uint64* table, uint64* y, uint64* x) {
  int width = 1, height = n;
  for (int i = 0; i < log2n; ++i) {
    height /= 2;
    if (i % 2) {
      Radix2(width, height, table, y, x);
    } else {
      Radix2(width, height, table, x, (i == log2n - 1) ? x : y);
    }
    table += width;
    width *= 2;
  }
}

void NTT::Radix2(const int width, const int height,
                 uint64* ptr, uint64* x, uint64* y) {
  for (int j = 0; j < height; ++j) {
    for (int i = 0; i < width; ++i) {
      int ix0 = j * width + i, ix1 = ix0 + height * width;
      int iy0 = j * 2 * width + i, iy1 = iy0 + width;
      uint64 t = Mod::mul(x[ix1], ptr[i]);
      y[iy1] = Mod::sub(x[ix0], t);
      y[iy0] = Mod::add(x[ix0], t);
    }
  }
}
