#pragma once

#include "fwd.hpp"

// **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //

namespace Arithmetic
{
    using index_t = std::size_t;
    using integer_t = std::make_signed_t<index_t>;

#ifdef ARITHMETIC_REAL_T
#if ARITHMETIC_REAL_T > 2
#error "Unknown value for ARITHMETIC_REAL_T"
#endif
#endif
#if ARITHMETIC_REAL_T == 1
    using real_t = double;
#elif ARITHMETIC_REAL_T == 2
    using real_t = long double;
#else
    using real_t = float;
#endif
    using complex_t = std::complex<real_t>;

    extern constexpr index_t negative_1 = static_cast<index_t>(-1);
}

// **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //
