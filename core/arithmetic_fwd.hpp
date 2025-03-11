#pragma once

#include "fwd.hpp"

// **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //

namespace Arithmetic
{
    using index_t = std::size_t;
    using integer_t = std::make_signed_t<index_t>;

    using real_t = double;
    using complex_t = std::complex<real_t>;

    extern constexpr index_t negative_1 = static_cast<index_t>(-1);
}

// **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //
