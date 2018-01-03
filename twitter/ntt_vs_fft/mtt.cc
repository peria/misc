#include "mtt.h"

#include <cstring>
#include <iostream>

using namespace std;

#if !defined(UINT128)
#error "Run on the environment that supports 128-bit integer"
#endif

namespace {
uint64* g_work = nullptr;
uint64* g_table = nullptr;
uint64* g_inv_table = nullptr;

using Montgomery = uint64;
}

namespace montgomery {

// TODO: Use smaller P for modular to reduce computing cost of reduce().
constexpr uint64 P = 0xffffffff00000001ULL;  // 2^64-2^32+1

inline uint64 reduce(const uint128& a) {
  static constexpr uint64 kInverse = 0xfffffffeffffffffULL;
  uint64 a64 = a;
  // r128 can be > 2^128 after the addition.
  uint128 r128 = (static_cast<uint128>(a64 * kInverse) * P + a);
  uint64 r = r128 >> 64;
  return (r >= P || r128 < a) ? (r - P) : r;
}

inline Montgomery convert(uint64 a) {
  static constexpr uint64 kR2 = 0xfffffffe00000001ULL;  // (2^64)^2 % P
  return reduce(static_cast<uint128>(a) * kR2);
}

inline Montgomery sub(Montgomery a, Montgomery b) {
  Montgomery c = a - b;
  return (a < b) ? (c + P) : c;
}

inline Montgomery add(Montgomery a, Montgomery b) {
  return sub(a, P - b);
}

}  // namespace montgomery

void MTT::SetGlobals(void* work, void* table, void* inv_table) {
  g_work = static_cast<uint64*>(work);
  g_table = static_cast<uint64*>(table);
  g_inv_table = static_cast<uint64*>(inv_table);
}

bool MTT::Validate(uint64* data) {
  for (int log2n = 1; log2n < 10; ++log2n) {
    const int n = 1 << log2n;
    InitTable(log2n);
    for (int i = 0; i < n; ++i)
      data[i] = i;
    Forward(log2n, n, data);
    Backward(log2n, n, data);
    for (int i = 0; i < n; ++i) {
      if (data[i] != i) {
        cout << "2^" << log2n << " [" << i << "] " << std::hex << data[i] << "\n";
        return false;
      }
    }
  }
  return true;
}

void MTT::Square(uint64* data, int n) {
  for (int i = 0; i < n; ++i)
    data[i] = montgomery::reduce(static_cast<uint128>(data[i]) * data[i]);
}

void MTT::InitTable(int log2n) {
  static constexpr Montgomery kRoots[] = {
    0xfffffffe00000002, 0xfffffffeffff0001, 0xfeffffff01000001, 0x0000000010000000,
    0xffffbfff00000001, 0xfffffffeffffff81, 0x07fffffffffff800, 0xe60ca9645a7a425e,
    0x5c411f4d8ab91088, 0x8bfed970d671fbb7, 0x1da1c8cedc0a82b1, 0x959dfcb4779eb1b1,
    0x35d17996b4e99746, 0x10bba1e10e56548b, 0x2306baaae6467556, 0xbf79450ceba724c2,
    0xaa3d8a0ca9f1cf0a, 0x05f9beab78de26d9, 0x8caa33007781b093, 0x5e93e76c70b1e9c6,
    0x32322652d8cb2ab7, 0xe67246b3ce63a09e, 0x36fbc989de66dc62, 0xc307e16fb62a525e,
    0x6ecfefd745751a91, 0x78d6e28499e74d1f, 0x915a171c5dce5b0b, 0x004a4484a6b1267b,
    0xa46d26647bea105f, 0xb86a0843c8fa27b2, 0x5588e6586a6c9a32, 0xda58878b0d514e98,
  };
  static constexpr Montgomery kInvRoots[] = {
    0xfffffffe00000002, 0x0000000000010000, 0xfffffeff00000001, 0xfffffffefffffff1,
    0xfffbffff00040001, 0x0000000002000000, 0x00000010000ffff0, 0xa4a4eeaf3df82e41,
    0xcf496a72f2d4f925, 0xd8d87587eb8ca0a6, 0x66a8e50c22b30ceb, 0x1fbd6e9f4552731d,
    0x9eb031492a3d5599, 0x895adfb5360defa9, 0xc8241328224d0c6f, 0x9a04ed17f8c6f7cb,
    0xb8343b3d7b2d6a6a, 0x3f8b95d59d926e3a, 0x234c7df8961145bf, 0x9fe6ed797090bbef,
    0x564a23670495d3ec, 0x36458c0ff54c19a3, 0x72cf655d4176ebbf, 0x63814db68b932df7,
    0x71bb0b99ce256069, 0xbb2075985a8c94d3, 0x38e675d8352604e8, 0xd8857182e8d2a3c4,
    0x4454645bd6da5192, 0xeec4f0092ef7c6ec, 0x7ab506ffef8a0c65, 0xb6fc8717d24cc2b2,
  };

  uint64* ptr0 = g_table;
  uint64* ptr1 = g_inv_table;
  int width = 1;
  for (int k = 0; k < log2n; ++k) {
    const uint128 w = kRoots[k];
    const uint128 m = kInvRoots[k];
    Montgomery wi = montgomery::convert(1), mi = montgomery::convert(1);
    for (int i = 0; i < width; ++i) {
      ptr0[i] = wi;
      ptr1[i] = mi;
      wi = montgomery::reduce(wi * w);
      mi = montgomery::reduce(mi * m);
    }
    ptr0 += width;
    ptr1 += width;
    width *= 2;
  }
}

void MTT::Forward(int log2n, int n, uint64* data) {
  for (int i = 0; i < n; ++i)
    data[i] = montgomery::convert(data[i]);
  Core(log2n, n, g_table, g_work, data);
}

void MTT::Backward(int log2n, int n, uint64* data) {
  Core(log2n, n, g_inv_table, g_work, data);

  static constexpr uint64 kInverses[] = {
    0,
    0x7FFFFFFF80000001, 0xBFFFFFFF40000001, 0xDFFFFFFF20000001, 0xEFFFFFFF10000001,
    0xF7FFFFFF08000001, 0xFBFFFFFF04000001, 0xFDFFFFFF02000001, 0xFEFFFFFF01000001,
    0xFF7FFFFF00800001, 0xFFBFFFFF00400001, 0xFFDFFFFF00200001, 0xFFEFFFFF00100001,
    0xFFF7FFFF00080001, 0xFFFBFFFF00040001, 0xFFFDFFFF00020001, 0xFFFEFFFF00010001,
    0xFFFF7FFF00008001, 0xFFFFBFFF00004001, 0xFFFFDFFF00002001, 0xFFFFEFFF00001001,
    0xFFFFF7FF00000801, 0xFFFFFBFF00000401, 0xFFFFFDFF00000201, 0xFFFFFEFF00000101,
    0xFFFFFF7F00000081, 0xFFFFFFBF00000041, 0xFFFFFFDF00000021, 0xFFFFFFEF00000011,
    0xFFFFFFF700000009, 0xFFFFFFFB00000005, 0xFFFFFFFD00000003, 0xFFFFFFFE00000002,
    0x7FFFFFFF00000001, 0xBFFFFFFF00000001, 0xDFFFFFFF00000001, 0xEFFFFFFF00000001,
    0xF7FFFFFF00000001, 0xFBFFFFFF00000001, 0xFDFFFFFF00000001, 0xFEFFFFFF00000001,
    0xFF7FFFFF00000001, 0xFFBFFFFF00000001, 0xFFDFFFFF00000001, 0xFFEFFFFF00000001,
    0xFFF7FFFF00000001, 0xFFFBFFFF00000001, 0xFFFDFFFF00000001, 0xFFFEFFFF00000001,
    0xFFFF7FFF00000001, 0xFFFFBFFF00000001, 0xFFFFDFFF00000001, 0xFFFFEFFF00000001,
    0xFFFFF7FF00000001, 0xFFFFFBFF00000001, 0xFFFFFDFF00000001, 0xFFFFFEFF00000001,
    0xFFFFFF7F00000001, 0xFFFFFFBF00000001, 0xFFFFFFDF00000001, 0xFFFFFFEF00000001,
    0xFFFFFFF700000001, 0xFFFFFFFB00000001, 0xFFFFFFFD00000001, 0xFFFFFFFE00000001
  };

  const uint128 inv = kInverses[log2n];
  for (int i = 0; i < n; ++i) {
    data[i] = montgomery::reduce(data[i] * inv);
  }
}

void MTT::Core(int log2n, int n, uint64* table, uint64* y, uint64* x) {
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

void MTT::Radix2(const int width, const int height,
                 uint64* ptr, uint64* x, uint64* y) {
  for (int j = 0; j < height; ++j) {
    for (int i = 0; i < width; ++i) {
      int ix0 = j * width + i, ix1 = ix0 + height * width;
      int iy0 = j * 2 * width + i, iy1 = iy0 + width;

      Montgomery t = montgomery::reduce(static_cast<uint128>(x[ix1]) * ptr[i]);
      y[iy1] = montgomery::sub(x[ix0], t);
      y[iy0] = montgomery::add(x[ix0], t);
    }
  }
}
