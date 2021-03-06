#include <chrono>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <iostream>
#include <memory>

#ifdef __NVCC__
#include <cuda.h>
#include <cuda_runtime.h>
#endif

#ifdef __NVCC__
#define HOST   __host__
#define DEVICE __device__
#define GLOBAL __global__
#else
#define HOST
#define DEVICE
#define GLOBAL
#endif

using uint128_t = __uint128_t;
using Clock = std::chrono::steady_clock;

constexpr int64_t kComputeSize = 4;
constexpr int64_t kNumBlocks = 128;
constexpr int64_t kNumGrids = 128;
constexpr int64_t kNumThreads = kNumBlocks * kNumGrids;
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

struct Output {
  Output(const int64_t size, const int64_t num_threads) {
    part_pi = new uint64_t[size];
    for (int64_t i = 0; i < size; ++i)
      part_pi[i] = 0;
    term_sums = new uint64_t[2 * size];
#ifdef __NVCC__
    cudaMalloc(&thread_sums, sizeof(uint64_t) * size * num_threads);
#else
    thread_sums = new uint64_t[size * num_threads];
#endif
  }
  ~Output() {
    delete[](part_pi);
    delete[](term_sums);
#ifdef __NVCC__
    cudaFree(&thread_sums);
#else
    delete[](thread_sums);
#endif
  }

  uint64_t* part_pi = nullptr;
  uint64_t* term_sums = nullptr;
  uint64_t* thread_sums = nullptr;
};

DEVICE HOST
inline uint64_t umul64hi(const uint64_t x, const uint64_t y) {
#ifdef __CUDA_ARCH__
  return __umul64hi(x, y);
#elif 1
  // Asm
  uint64_t z;
  asm("movq %1,%%rax\n\t"
      "mul %2\n\t"
      "movq %%rdx,%0\n\t"
      : "=r"(z)
      : "r"(x), "r"(y)
      : "rax", "rdx");
  return z;
#elif 1
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
DEVICE HOST
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

// https://gmplib.org/~tege/division-paper.pdf
DEVICE HOST
uint64_t GetInverse(uint64_t d) {
#if 0
  // Uint128
  static constexpr uint128_t kFull128 = -uint128_t(1);
  return uint64_t(kFull128 / d);
#else
  // Pure C++
  uint64_t d0 = d & 1;
  uint64_t d9 = d >> 55;
  uint64_t d40 = (d >> 24) + 1;
  uint64_t d63 = (d >> 1) + (d & 1);
  uint64_t v0 = ((1ULL << 19) - (3ULL << 8)) / d9;  // Table
  uint64_t v1 = (v0 << 11) - ((v0 * v0 * d40) >> 40) - 1;
  uint64_t v2 = (v1 << 13) + ((v1 * ((1ULL << 60) - v1 * d40)) >> 47);
  uint64_t e = -v2 * d63 + (v2 / 2) * d0;
  uint64_t v3 = (v2 << 31) + (umul64hi(v2, e) >> 1);
  uint64_t v4 = v3 - (umul64hi(v3 + 1, d) + d);
  return v4;
#endif
}

DEVICE HOST
uint64_t Div2ByNormalized1(uint64_t r0,
                           uint64_t r1,
                           const uint64_t d,
                           uint64_t& rem) {
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
    --r1;
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

  rem = r0;
  return q;
}

DEVICE HOST
uint64_t Div2ByNormalized1(const uint64_t u0,
                           const uint64_t u1,
                           const uint64_t d,
                           const uint64_t v,
                           uint64_t& rem) {
#if 0
  // Uint128
  uint128_t q = uint128_t(v) * u1;
  q += (uint128_t(u1) << 64) + u0;
  uint64_t q1 = q >> 64;
  ++q1;
  uint64_t r = u0 - q1 * d;
  if (r > uint64_t(q)) {
    --q1;
    r += d;
  }
  if (r >= d) {
    ++q1;
    r -= d;
  }
  rem = r;
  return q1;
#else
  // Pure C++
  uint64_t q0 = v * u1;
  uint64_t q1 = umul64hi(v, u1);
  q1 += u1;
  q0 += u0;
  if (q0 < u0)
    ++q1;
  ++q1;
  uint64_t r = u0 - q1 * d;
  if (r > q0) {
    --q1;
    r += d;
  }
  if (r >= d) {
    ++q1;
    r -= d;
  }
  rem = r;
  return q1;
#endif
}

// Computes a^e mod m using Montgomery Reduction.
// Reference; https://min-25.hatenablog.com/entry/2017/08/20/171214
DEVICE HOST
uint64_t PowMod(uint64_t a, uint64_t e, const uint64_t m) {
  uint64_t inv = m;
  for (int i = 0; i < 5; ++i)
    inv *= 2 - inv * m;
  uint64_t mm = m;
  int64_t shift = NormalizeDivider(mm);
  uint64_t r2;
#if 0
  // Uint128
  r2 = -uint128_t(m) % m;
#else
  // Pure C++
  {
    uint64_t m0 = -m;
    uint64_t m1 = (~0ULL) % m;
    m1 = (m1 << shift) + (m0 >> (64 - shift));
    m0 <<= shift;
    uint64_t inv = GetInverse(mm);
    Div2ByNormalized1(m0, m1, mm, inv, r2);
  }
  r2 >>= shift;
#endif

  uint64_t r = -umul64hi(r2 * inv, m);
  r += (r >> 63) * m;
  a *= r;

  if (e & 1) {
    r = umul64hi(r, a) - umul64hi(r * a * inv, m);
    r += (r >> 63) * m;
  }
  while (e >>= 1) {
    a = umul64hi(a, a) - umul64hi(a * a * inv, m);
    a += (a >> 63) * m;
    if (e & 1) {
      r = umul64hi(r, a) - umul64hi(r * a * inv, m);
      r += (r >> 63) * m;
    }
  }
  r = -umul64hi(r * inv, m);
  r += (r >> 63) * m;
  return r;
}

DEVICE HOST
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
  uint64_t inv = GetInverse(nb);
  uint64_t r1 = a << shift;
  for (int64_t i = size - 1; i >= 0; --i) {
    uint64_t r0 = 0;
    dst[i] = Div2ByNormalized1(r0, r1, nb, inv, r1);
  }
#endif
}

DEVICE HOST
void Div(const uint64_t* a,
         const uint64_t b,
         const int64_t size,
         uint64_t* dst) {
#if 0
  // Uint128
  uint128_t r;
  for (int64_t i = size - 1; i >= 0; --i) {
    r = (r << 64) + a[i];
    dst[i] = r / b;
    r %= b;
  }
#else
  // Pure C++
  uint64_t nb = b;
  int64_t shift = NormalizeDivider(nb);
  uint64_t inv = GetInverse(nb);
  uint64_t r1 = 0;
  for (int64_t i = size - 1; i >= 0; --i) {
    r1 |= a[i] >> (64 - shift);
    uint64_t r0 = a[i] << shift;
    dst[i] = Div2ByNormalized1(r0, r1, nb, inv, r1);
  }
#endif
}

DEVICE HOST
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

DEVICE HOST
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

GLOBAL
void ComputeIntegralKernel(
#ifndef __NVCC__
    const int64_t thread_id,
#endif
    const int64_t n,
    const int64_t bit_shift,
    const int64_t a,
    const int64_t c,
    const int64_t d,
    uint64_t* thread_sums) {
#ifdef __NVCC__
  const int64_t thread_id = blockIdx.x * blockDim.x + threadIdx.x;
#endif
  uint64_t* thread_sum = thread_sums + thread_id * kComputeSize;
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
                         const Term& term,
                         Output& output) {
  const int64_t n_in_parallel = n / kNumThreads * kNumThreads;

#ifdef __NVCC__
  ComputeIntegralKernel<<<kNumGrids, kNumBlocks>>>(
      n_in_parallel, bit_shift, term.a, term.c, term.d, output.thread_sums);
#else
  for (int64_t thread_id = 0; thread_id < kNumThreads; ++thread_id) {
    ComputeIntegralKernel(thread_id, n_in_parallel, bit_shift, term.a, term.c,
                          term.d, output.thread_sums);
  }
#endif

  uint64_t* term_sums[] = {output.term_sums, output.term_sums + kComputeSize};
  // single thread in CPU
  for (int64_t i = n_in_parallel; i < n; ++i) {
    const int64_t mod = term.c * i + term.d;
    const int64_t pow2 = bit_shift - term.a * i;
    uint64_t operand[kComputeSize];
    const uint64_t rem = PowMod(2, pow2, mod);
    DivMod1(rem, mod, kComputeSize, operand);
    Add(term_sums[i & 1], operand, kComputeSize, term_sums[i & 1]);
  }

  std::unique_ptr<uint64_t[]> thread_sums(
      new uint64_t[kComputeSize * kNumThreads]);
  // Copy GPU memory to main memory.
#ifdef __NVCC__
  cudaMemcpy(thread_sums.get(), output.thread_sums,
             sizeof(uint64_t) * kComputeSize * kNumThreads,
             cudaMemcpyDeviceToHost);
#else
  std::memcpy(thread_sums.get(), output.thread_sums,
              sizeof(uint64_t) * kComputeSize * kNumThreads);
#endif
  for (int64_t i = 0; i < kNumThreads; ++i) {
    uint64_t* thread_sum = thread_sums.get() + i * kComputeSize;
    Add(term_sums[i & 1], thread_sum, kComputeSize, term_sums[i & 1]);
  }
}

void ComputeTerm(const Term& term, const int64_t hex_index, Output& output) {
  const int64_t bit_index = hex_index * 4 - 3;
  const int64_t bit_shift = bit_index + term.b - 1;
  const int64_t integer_n = (bit_shift >= 0) ? (bit_shift / term.a) : 0LL;
  const int64_t zero_n = (bit_shift + 64 * kComputeSize) / term.a;

  for (int64_t i = 0; i < 2 * kComputeSize; ++i)
    output.term_sums[i] = 0;

  ComputeIntegralPart(integer_n, bit_shift, term, output);

  uint64_t* term_sums[] = {output.term_sums, output.term_sums + kComputeSize};
  for (int64_t i = integer_n; i < zero_n; ++i) {
    int64_t shift = bit_shift + 64 * kComputeSize - i * term.a;
    const uint64_t mod = term.c * i + term.d;
    uint64_t operand[kComputeSize + 1]{};
    operand[shift / 64] = 1ULL << (shift % 64);
    Div(operand, mod, kComputeSize + 1, operand);
    Add(term_sums[i & 1], operand, kComputeSize, term_sums[i & 1]);
  }

  if (term.f == Term::Flip::kFlip) {
    Subtract(term_sums[0], term_sums[1], kComputeSize, term_sums[0]);
  } else {
    Add(term_sums[0], term_sums[1], kComputeSize, term_sums[0]);
  }

  if (term.s == Term::Sign::kPositive) {
    Add(output.part_pi, term_sums[0], kComputeSize, output.part_pi);
  } else {
    Subtract(output.part_pi, term_sums[0], kComputeSize, output.part_pi);
  }
}

void Dump(const uint64_t* a, const int64_t size) {
  for (int64_t i = size - 1; i >= 0; --i) {
    std::printf("%016lx", a[i]);
  }
  std::printf("\n");
}

double DurationInSec(const Clock::time_point& a, const Clock::time_point& b) {
  if (a > b)
    return DurationInSec(b, a);

  using MSec = std::chrono::milliseconds;
  auto diff = b - a;
  return std::chrono::duration_cast<MSec>(diff).count() * 1e-3;
}

int main() {
  // HexIndex is 1-origin index.
  static constexpr int64_t kHexIndex = 10000000;

  const Term kTerms[] = {
      {10, -1, 4, 1, Term::Sign::kNegative, Term::Flip::kFlip},
      {10, -6, 4, 3, Term::Sign::kNegative, Term::Flip::kFlip},
      {10, 2, 10, 1, Term::Sign::kPositive, Term::Flip::kFlip},
      {10, 0, 10, 3, Term::Sign::kNegative, Term::Flip::kFlip},
      {10, -4, 10, 5, Term::Sign::kNegative, Term::Flip::kFlip},
      {10, -4, 10, 7, Term::Sign::kNegative, Term::Flip::kFlip},
      {10, -6, 10, 9, Term::Sign::kPositive, Term::Flip::kFlip},
  };

  Output output(kComputeSize, kNumThreads);
  auto start_time = Clock::now();
  for (auto&& term : kTerms) {
    auto term_start_time = Clock::now();
    ComputeTerm(term, kHexIndex, output);
    auto term_end_time = Clock::now();
    double term_time = DurationInSec(term_end_time, term_start_time);
    std::printf("%.3fsec Term[1/(%2ldk%+ld)] : ", term_time, term.c, term.d);
    Dump(output.part_pi, kComputeSize);
  }
  auto end_time = Clock::now();

  double total_time = DurationInSec(end_time, start_time);
  std::cout << (kComputeSize * 16) << " hex digits of pi from " << kHexIndex
            << " th hex digit are\n";
  for (int64_t i = kComputeSize - 1; i >= 0; --i) {
    std::printf("%016lX", output.part_pi[i]);
  }
  std::printf("\nCompute Time: %.3f sec.\n", total_time);

  return 0;
}
