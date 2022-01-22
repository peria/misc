#pragma once

#include <cassert>
#include <cstdint>

using int64 = std::int64_t;
using uint64 = std::uint64_t;
using uint128 = __uint128_t;

#ifndef LOCAL
#define ASSERT(x) assert(x)
#else
#define ASSERT(x) \
  {}
#endif