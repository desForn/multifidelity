#pragma once

#include "core/core.hpp"

// **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //

namespace Matrix
{
    using Arithmetic::index_t;
    using Arithmetic::integer_t;
    using Arithmetic::real_t;
    using Arithmetic::complex_t;
    using Arithmetic::negative_1;

    // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //

    // Exception type
    struct invalid_index { std::array<index_t, 2> index; };
    
    // Maximum size of arrays. Larger values are stored as vectors
    constexpr index_t maximum_array_storage = 128;

    template<class field_type = real_t>
    class matrix_base;

    // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //

    template<class matrix_type>
    concept matrix_c =
        std::derived_from<matrix_type, matrix_base<typename matrix_type::field_type>>;

    template<class matrix_type, class field_type>
    concept matrix_field_type_c = matrix_c<matrix_type> and
        std::same_as<typename matrix_type::field_type, field_type>;
}

// **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //

