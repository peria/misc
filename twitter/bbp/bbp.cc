#include <cstdint>
#include <cstdio>
#include <iostream>

using uint128_t = __uint128_t;

uint64_t PowMod(uint64_t a, uint64_t e, const uint64_t m) {
  uint128_t a2 = a;
  uint64_t r = 1;
  for (; e; e >>= 1) {
    if (e & 1) {
      r = a * a2 % m;
    }
    a2 = a2 * a2 % m;
  }
  return r;
}

void DivMod1(const uint64_t a, const uint64_t b, const int64_t size, uint64_t* dst) {
  uint128_t t = a;
  for (int64_t i = 0; i < size; ++i) {
    t <<= 64;
    dst[i] = t / b;
    t %= b;
  }
}

void Div(const uint64_t* a, const uint64_t b, const int64_t size, uint64_t* dst) {
  uint128_t t;
  for (int64_t i = 0; i < size; ++i) {
    t = (t << 64) + a[i];
    dst[i] = t / b;
    t %= b;
  }
}

void Add(const uint64_t* a, const uint64_t* b, const int64_t size, uint64_t* dst) {
  uint64_t carry = 0;
  for (int64_t i = size - 1; i >= 0; --i) {
    uint64_t x = a[i];
    uint64_t y = b[i];
    uint64_t z = x + carry;
    if (z < x)
      carry = 1;
    uint64_t w = y + z;
    if (w < z)
      carry = 1;
    dst[i] = w;
  }
}

void Subtract(const uint64_t* a, const uint64_t* b, const int64_t size, uint64_t* dst) {
  uint64_t carry = 0;
  for (int64_t i = size - 1; i >= 0; --i) {
    uint64_t x = a[i];
    uint64_t y = b[i];
    uint64_t z = x - carry;
    if (z > x)
      carry = 1;
    uint64_t w = z - y;
    if (w > z)
      carry = 1;
    dst[i] = w;
  }
}

int main() {
  // XxxIndex are 1-origin indecies.
  static constexpr int64_t kHexIndex = 1;
  static constexpr int64_t kBitIndex = kHexIndex * 4 - 3;
  static constexpr int64_t kComputeSize = 4;

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
  const Term kTerms[] = {
    {10, -1, 4, 1, Term::Sign::kNegative, Term::Flip::kFlip},
    {10, -6, 4, 3, Term::Sign::kNegative, Term::Flip::kFlip},
    {10, 2, 10, 1, Term::Sign::kPositive, Term::Flip::kFlip},
    {10, 0, 10, 3, Term::Sign::kNegative, Term::Flip::kFlip},
    {10, -4, 10, 5, Term::Sign::kNegative, Term::Flip::kFlip},
    {10, -4, 10, 7, Term::Sign::kNegative, Term::Flip::kFlip},
    {10, -6, 10, 9, Term::Sign::kPositive, Term::Flip::kFlip},
  };

  uint64_t part_pi[kComputeSize] {};
  for (auto&& term : kTerms) {
    uint64_t term_sum[2][kComputeSize] {}; // [0]: k is even, [1]: k is odd.
    const int64_t bit_shift = kBitIndex + term.b - 1;
    const int64_t integer_n = (bit_shift >= 0) ? (bit_shift / term.a) : -1LL;
    const int64_t zero_n = (bit_shift + 64 * kComputeSize) / term.a;

    for (int64_t i = 0; i < integer_n; ++i) {
      const int64_t mod = term.c * i + term.d;
      const int64_t pow2 = bit_shift - term.a * i;
      uint64_t operand[kComputeSize];
      const uint64_t rem = PowMod(2, pow2, mod);
      DivMod1(rem, mod, kComputeSize, operand);
      Add(term_sum[i & 1], operand, kComputeSize, term_sum[i & 1]);
    }

    for (int64_t i = integer_n; i < zero_n; ++i) {
      int64_t shift = bit_shift + 64 * kComputeSize - i * term.a;
      const uint64_t mod = term.c * i + term.d;
      uint64_t operand[kComputeSize] {};
      operand[shift / 64] = 1ULL << (shift % 64);
      Div(operand, mod, kComputeSize, operand);
      Add(term_sum[i & 1], operand, kComputeSize, term_sum[i & 1]);
    }
    if (term.f == Term::Flip::kFlip) {
      Subtract(term_sum[0], term_sum[1], kComputeSize, term_sum[0]);
    } else {
      Add(term_sum[0], term_sum[1], kComputeSize, term_sum[0]);
    }
    Add(part_pi, term_sum[0], kComputeSize, part_pi);
  }

  for (int64_t i = 0; i < kComputeSize; ++i) {
    std::printf("%016lX", part_pi[i]);
  }
  std::puts("");

  return 0;
}
