#include <chrono>
#include <cstdint>
#include <cstdio>
#include <iostream>
#include <memory>

using uint128_t = __uint128_t;
using Clock = std::chrono::steady_clock;

constexpr int64_t kComputeSize = 4;
constexpr int64_t kNumThreads = 4;
static_assert(kNumThreads % 2 == 0, "Number of threads must be even");

// A = s * \sum_k 2^(-a*k+b)/(c*k+d)         (if f == kNoFlip)
//     s * \sum_k (-1)^k 2^(-a*k+b)/(c*k+d)  (if f == kFlip)
struct Term {
  enum class Sign : uint8_t { kPositive, kNegative };
  enum class Flip : uint8_t { kNoFlip, kFlip };
  int64_t a;
  int64_t b;
  int64_t c;
  int64_t d;
  Sign s;
  Flip f;
};

inline uint64_t umul64hi(uint64_t x, uint64_t y) {
#if 1
  // Asm
  uint64_t z;
  asm("mul %2" : "=d"(z) : "a"(x), "r"(y));
  return z;
#elif 0
  // Uint128
  return (uint128_t(x) * y) >> 64;
#else
  // Pure C++
  static constexpr uint64_t kMask = 0xffffffffULL;
  uint64_t x0 = x & kMask;
  uint64_t x1 = x >> 32;
  uint64_t y0 = y & kMask;
  uint64_t y1 = y >> 32;
  uint64_t t = ((x0 * y0) >> 32) + ((x1 * y0) & kMask) + ((x0 * y1) & kMask);
  uint64_t r = x1 * y1 + ((x1 * y0) >> 32) + ((x0 * y1) >> 32) + (t >> 32);
  return r;
#endif
}

// Make 2^63 <= |d| < 2^64
int64_t NormalizeDivider(uint64_t& d) {
  int64_t shift = 0;
  if ((d & (-1ULL << 32)) == 0) {
    shift += 32;
    d <<= 32;
  }
  if ((d & (-1ULL << 48)) == 0) {
    shift += 16;
    d <<= 16;
  }
  if ((d & (-1ULL << 56)) == 0) {
    shift += 8;
    d <<= 8;
  }
  if ((d & (-1ULL << 60)) == 0) {
    shift += 4;
    d <<= 4;
  }
  if ((d & (-1ULL << 62)) == 0) {
    shift += 2;
    d <<= 2;
  }
  if ((d & (-1ULL << 63)) == 0) {
    shift += 1;
    d <<= 1;
  }
  return shift;
}

uint64_t Div2ByNormalized1(uint64_t n0,
                           uint64_t n1,
                           const int64_t shift,
                           const uint64_t d,
                           uint64_t& rem) {
  // Nomalize numerator ==> remain
  uint64_t r0 = n0 << shift;
  uint64_t r1 = (n1 << shift) | (n0 >> (64 - shift));

  uint64_t d1 = d >> 32;
  uint64_t qt = r1 / d1;

  uint64_t p0 = d * qt;
  uint64_t p1 = umul64hi(d, qt);
  {
    const uint64_t u0 = (r1 << 32) | (r0 >> 32);
    const uint64_t u1 = r1 >> 32;
    while (p1 > u1 || (p1 == u1 && p0 > u0)) {
      if (p0 < d)
        --p1;
      p0 -= d;
      --qt;
    }
  }

  uint64_t q = qt << 32;
  uint64_t w0 = p0 << 32;
  uint64_t w1 = (p1 << 32) | (p0 >> 32);
  if (w0 > r0)
    --w1;
  r0 -= w0;
  r1 -= w1;

  qt = ((r1 << 32) | (r0 >> 32)) / d1;
  p0 = d * qt;
  p1 = umul64hi(d, qt);
  while (p1 > r1 || (p1 == r1 && p0 > r0)) {
    if (p0 < d)
      --p1;
    p0 -= d;
    --qt;
  }
  r0 -= p0;
  q |= qt;

  rem = r0 >> shift;
  return q;
}

// Computes a^e mod m using Montgomery Reduction.
// Reference; https://min-25.hatenablog.com/entry/2017/08/20/171214
uint64_t PowMod(uint64_t a, uint64_t e, const uint64_t m) {
  uint64_t inv = m;
  for (int i = 0; i < 5; ++i)
    inv *= 2 - inv * m;
  uint64_t mm = m;
  int64_t shift = NormalizeDivider(mm);
  uint64_t r2;
  Div2ByNormalized1(-m, (~0ULL) % m, shift, mm, r2);

  uint64_t r = -umul64hi(r2 * inv, m);
  r = (r & (1ULL << 63)) ? (r + m) : r;
  a *= r;
  for (; e; e >>= 1) {
    if (e & 1) {
      r = umul64hi(r, a) - umul64hi(r * a * inv, m);
      r = (r & (1ULL << 63)) ? (r + m) : r;
    }
    a = umul64hi(a, a) - umul64hi(a * a * inv, m);
    a = (a & (1ULL << 63)) ? (a + m) : a;
  }
  r = -umul64hi(r * inv, m);
  return (r & (1ULL << 63)) ? (r + m) : r;
}

void DivMod1(const uint64_t a,
             const uint64_t b,
             const int64_t size,
             uint64_t* dst) {
#if 0
  // Uint128
  uint128_t t = a;
  for (int64_t i = size - 1; i >= 0; --i) {
    t <<= 64;
    dst[i] = t / b;
    t %= b;
  }
#else
  // Pure C++
  uint64_t nb = b;
  int64_t shift = NormalizeDivider(nb);
  uint64_t b1 = nb >> 32;
  uint64_t r0 = a << shift;
  for (int64_t i = size - 1; i >= 0; --i) {
    uint64_t r1 = r0;
    r0 = 0;
    uint64_t qt = r1 / b1;
    uint64_t p0 = b * qt;
    uint64_t p1 = umul64hi(nb, qt);
    {
      const uint64_t u0 = (r1 << 32) | (r0 >> 32);
      const uint64_t u1 = r1 >> 32;
      while (p1 > u1 || (p1 == u1 && p0 > u0)) {
        if (p0 < nb)
          --p1;
        p0 -= nb;
        --qt;
      }
    }
    uint64_t q = qt << 32;
    uint64_t w0 = p0 << 32;
    uint64_t w1 = (p1 << 32) | (p0 >> 32);
    if (w0 > r0)
      --w1;
    r0 -= w0;
    r1 -= w1;

    qt = ((r1 << 32) | (r0 >> 32)) / b1;
    p0 = nb * qt;
    p1 = umul64hi(nb, qt);
    while (p1 > r1 || (p1 == r1 && p0 > r0)) {
      if (p0 < nb)
        --p1;
      p0 -= nb;
      --qt;
    }
    r0 -= p0;
    dst[i] = q | qt;
  }
#endif
}

void Div(const uint64_t* a,
         const uint64_t b,
         const int64_t size,
         uint64_t* dst) {
#if 0
  uint128_t t;
  for (int64_t i = size - 1; i >= 0; --i) {
    t = (t << 64) + a[i];
    dst[i] = t / b;
    t %= b;
  }
#else
  // Pure C++
  uint64_t nb = b;
  int64_t shift = NormalizeDivider(nb);
  uint64_t b1 = nb >> 32;
  uint64_t r0 = 0;
  for (int64_t i = size - 1; i >= 0; --i) {
    uint64_t r1 = r0;
    r0 = a[i];
    r1 |= r0 >> (64 - shift);
    r0 <<= shift;
    uint64_t qt = r1 / b1;
    uint64_t p0 = b * qt;
    uint64_t p1 = umul64hi(nb, qt);
    {
      const uint64_t u0 = (r1 << 32) | (r0 >> 32);
      const uint64_t u1 = r1 >> 32;
      while (p1 > u1 || (p1 == u1 && p0 > u0)) {
        if (p0 < nb)
          --p1;
        p0 -= nb;
        --qt;
      }
    }
    uint64_t q = qt << 32;
    uint64_t w0 = p0 << 32;
    uint64_t w1 = (p1 << 32) | (p0 >> 32);
    if (w0 > r0)
      --w1;
    r0 -= w0;
    r1 -= w1;

    qt = ((r1 << 32) | (r0 >> 32)) / b1;
    p0 = nb * qt;
    p1 = umul64hi(nb, qt);
    while (p1 > r1 || (p1 == r1 && p0 > r0)) {
      if (p0 < nb)
        --p1;
      p0 -= nb;
      --qt;
    }
    r0 -= p0;
    dst[i] = q | qt;
  }
#endif
}

void Add(const uint64_t* a,
         const uint64_t* b,
         const int64_t size,
         uint64_t* dst) {
  uint64_t carry = 0;
  for (int64_t i = 0; i < size; ++i) {
    uint64_t x = a[i];
    uint64_t y = b[i];
    uint64_t z = x + carry;
    carry = (z < x) ? 1 : 0;
    uint64_t w = z + y;
    if (w < z)
      carry = 1;
    dst[i] = w;
  }
}

void Subtract(const uint64_t* a,
              const uint64_t* b,
              const int64_t size,
              uint64_t* dst) {
  uint64_t carry = 0;
  for (int64_t i = 0; i < size; ++i) {
    uint64_t x = a[i];
    uint64_t y = b[i];
    uint64_t z = x - carry;
    carry = (z > x) ? 1 : 0;
    uint64_t w = z - y;
    if (w > z)
      carry = 1;
    dst[i] = w;
  }
}

void ComputeIntegralKernel(const int64_t thread_id,
                           const int64_t n,
                           const int64_t bit_shift,
                           const int64_t a,
                           const int64_t c,
                           const int64_t d,
                           uint64_t* thread_sum) {
  for (int64_t i = 0; i < kComputeSize; ++i)
    thread_sum[i] = 0;

  for (int64_t i = thread_id; i < n; i += kNumThreads) {
    const int64_t mod = c * i + d;
    const int64_t pow2 = bit_shift - a * i;
    uint64_t operand[kComputeSize];
    const uint64_t rem = PowMod(2, pow2, mod);
    DivMod1(rem, mod, kComputeSize, operand);
    Add(thread_sum, operand, kComputeSize, thread_sum);
  }
}

void ComputeIntegralPart(const int64_t n,
                         const int64_t bit_shift,
                         const int64_t a,
                         const int64_t c,
                         const int64_t d,
                         uint64_t term_sum[2][kComputeSize]) {
  const int64_t n_in_parallel = n / kNumThreads * kNumThreads;

  std::unique_ptr<uint64_t[]> thread_sums(
      new uint64_t[kComputeSize * kNumThreads]);
  for (int64_t thread_id = 0; thread_id < kNumThreads; ++thread_id) {
    uint64_t* thread_sum = thread_sums.get() + thread_id * kComputeSize;
    ComputeIntegralKernel(thread_id, n_in_parallel, bit_shift, a, c, d,
                          thread_sum);
  }

  // single thread in CPU
  for (int64_t i = n_in_parallel; i < n; ++i) {
    const int64_t mod = c * i + d;
    const int64_t pow2 = bit_shift - a * i;
    uint64_t operand[kComputeSize];
    const uint64_t rem = PowMod(2, pow2, mod);
    DivMod1(rem, mod, kComputeSize, operand);
    Add(term_sum[i & 1], operand, kComputeSize, term_sum[i & 1]);
  }

  // sync here

  for (int64_t i = 0; i < kNumThreads; ++i) {
    uint64_t* thread_sum = thread_sums.get() + i * kComputeSize;
    Add(term_sum[i & 1], thread_sum, kComputeSize, term_sum[i & 1]);
  }
}

void Dump(const uint64_t* a, const int64_t size) {
  for (int64_t i = size - 1; i >= 0; --i) {
    std::printf("%016lx", a[i]);
  }
  std::printf("\n");
}

int main() {
  // XxxIndex are 1-origin indecies.
  static constexpr int64_t kHexIndex = 10000000;
  static constexpr int64_t kBitIndex = kHexIndex * 4 - 3;

  const Term kTerms[] = {
      {10, -1, 4, 1, Term::Sign::kNegative, Term::Flip::kFlip},
      {10, -6, 4, 3, Term::Sign::kNegative, Term::Flip::kFlip},
      {10, 2, 10, 1, Term::Sign::kPositive, Term::Flip::kFlip},
      {10, 0, 10, 3, Term::Sign::kNegative, Term::Flip::kFlip},
      {10, -4, 10, 5, Term::Sign::kNegative, Term::Flip::kFlip},
      {10, -4, 10, 7, Term::Sign::kNegative, Term::Flip::kFlip},
      {10, -6, 10, 9, Term::Sign::kPositive, Term::Flip::kFlip},
  };

  uint64_t part_pi[kComputeSize]{};
  auto start_time = Clock::now();
  for (auto&& term : kTerms) {
    auto term_start_time = Clock::now();

    uint64_t term_sum[2][kComputeSize]{};  // [0]: k is even, [1]: k is odd.
    const int64_t bit_shift = kBitIndex + term.b - 1;
    const int64_t integer_n = (bit_shift >= 0) ? (bit_shift / term.a) : 0LL;
    const int64_t zero_n = (bit_shift + 64 * kComputeSize) / term.a;

    ComputeIntegralPart(integer_n, bit_shift, term.a, term.c, term.d, term_sum);

    for (int64_t i = integer_n; i < zero_n; ++i) {
      int64_t shift = bit_shift + 64 * kComputeSize - i * term.a;
      const uint64_t mod = term.c * i + term.d;
      uint64_t operand[kComputeSize + 1]{};
      operand[shift / 64] = 1ULL << (shift % 64);
      Div(operand, mod, kComputeSize + 1, operand);
      Add(term_sum[i & 1], operand, kComputeSize, term_sum[i & 1]);
    }

    if (term.f == Term::Flip::kFlip) {
      Subtract(term_sum[0], term_sum[1], kComputeSize, term_sum[0]);
    } else {
      Add(term_sum[0], term_sum[1], kComputeSize, term_sum[0]);
    }

    if (term.s == Term::Sign::kPositive) {
      Add(part_pi, term_sum[0], kComputeSize, part_pi);
    } else {
      Subtract(part_pi, term_sum[0], kComputeSize, part_pi);
    }
    auto term_end_time = Clock::now();

    double term_time = std::chrono::duration_cast<std::chrono::milliseconds>(
                           term_end_time - term_start_time)
                           .count() *
                       1e-3;
    std::printf("%.3fsec Term[1/(%2ldk%+ld)] : ", term_time, term.c, term.d);
    Dump(part_pi, kComputeSize);
  }
  auto end_time = Clock::now();

  double total_time = std::chrono::duration_cast<std::chrono::milliseconds>(
                          end_time - start_time)
                          .count() *
                      1e-3;
  std::cout << (kComputeSize * 16) << " hex digits of pi from " << kHexIndex
            << " th hex digit are\n";
  for (int64_t i = kComputeSize - 1; i >= 0; --i) {
    std::printf("%016lX", part_pi[i]);
  }
  std::printf("\nCompute Time: %.3f sec.\n", total_time);

  return 0;
}
