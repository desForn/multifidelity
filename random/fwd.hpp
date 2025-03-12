#pragma once

#include "core/core.hpp"

// **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //

namespace Random
{
    using Arithmetic::index_t;
    using Arithmetic::integer_t;
    using Arithmetic::real_t;
    using Arithmetic::complex_t;
    using Arithmetic::negative_1;

#ifndef RANDOM_SEED
    constexpr index_t default_seed = negative_1; // -1 to use time()
#else
    constexpr index_t default_seed = RANDOM_SEED;
#endif
}
