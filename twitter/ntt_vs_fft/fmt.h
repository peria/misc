#pragma once

#include <cstdint>
#include <memory>

using uint64 = std::uint64_t;

#if __GNUC__
#define UINT128
using uint128 = __uint128_t;
#endif
