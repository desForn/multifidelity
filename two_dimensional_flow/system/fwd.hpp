#pragma once

#include "sampling/chebyshev_lobatto.hpp"
#include "sampling/equispaced.hpp"
#include "sampling/exponential_sampling.hpp"

namespace Two_dimensional_flow
{
    using Arithmetic::index_t;
    using Arithmetic::integer_t;
    using Arithmetic::real_t;
    using Arithmetic::complex_t;
    using Arithmetic::negative_1;

    // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //
    // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //

    using fidelity_type = std::array<index_t, 1>;
    using coordinate_type = std::array<real_t, 1>;
    using coordinate_index_type = std::array<index_t, 1>;

    using fidelity_coordinate_type = std::tuple<fidelity_type, coordinate_type>;
    using fidelity_coordinate_index_type = std::tuple<fidelity_type, coordinate_index_type>;

    using output_type = std::array<real_t, 2>;
    using pair_type = std::pair<fidelity_coordinate_type, output_type>;
}

