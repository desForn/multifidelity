#pragma once

#include "core/core.hpp"

#include "sampling/exponential_sampling.hpp"
#include "sampling/chebyshev_lobatto.hpp"
#include "sampling/exponential_refinement.hpp"
#include "sampling/sampling_base.hpp"

#include "smolyak/smolyak_approximation.hpp"
#include "smolyak/smolyak_traits.hpp"

#include "matrix/dense_matrix.hpp"

// **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //

namespace Platoon
{
    using Arithmetic::index_t;
    using Arithmetic::integer_t;
    using Arithmetic::real_t;
    using Arithmetic::complex_t;
    using Arithmetic::negative_1;

    // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //

    namespace Domain
    {
        constexpr index_t n_parameters = 4;
        constexpr index_t n_fidelity_levels = 4;

        constexpr index_t min_geometry = 0;
        constexpr index_t max_geometry = 3;

        constexpr real_t min_spacing = 5.4375;
        constexpr real_t max_spacing = 10.875;

        constexpr real_t min_yaw = 0;
        constexpr real_t max_yaw = static_cast<real_t>(10 * M_PI) / 180;

        constexpr real_t min_reynolds = 3.84E5;
        constexpr real_t max_reynolds = 7.68E5;
    }

    // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //

    namespace Settings
    {
        constexpr std::array<index_t, Domain::n_fidelity_levels> max_jobs{10, 10, 10, 10};

        const std::array<index_t, Domain::n_fidelity_levels> partition{1, 16, 28, 64};
        const std::array<index_t, Domain::n_fidelity_levels> n_processors{1, 16, 56, 320};
        const std::array<real_t, Domain::n_fidelity_levels> start_average{1500, 1500, 1500, 1500};
    }

    // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //

    using sampling_type = Sampling::exponential<Sampling::chebyshev_lobatto>;
    using richardson_extrapolation_sampling_type = Sampling::exponential_refinement<0.5, 1>;

    using fidelity_traits = Smolyak::Smolyak_traits::multifidelity<>;
    using richardson_extrapolation_traits = Smolyak::Smolyak_traits::
            multifidelity_richardson_extrapolation<0.5>;
    using interpolation_traits = Smolyak::Smolyak_traits::polynomial_interpolation<sampling_type>;
    using matrix_traits = Smolyak::Smolyak_traits::matrix<true>;
    using matrix_traits_initialiser = Smolyak::Smolyak_traits::matrix_initialiser;

    // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //

    //                            component, fidelity, geometry, spacing, yaw,    reynolds
    using input_type = std::tuple<index_t,   index_t,  index_t,  real_t,  real_t, real_t>;
    using output_type = std::array<real_t, 2>; // {force, cost}
    using pair_type = std::pair<input_type, output_type>;
}

