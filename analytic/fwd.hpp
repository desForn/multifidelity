#pragma once

#include "core/core.hpp"

#include "invocable/invocable_function_hash_table.hpp"

#include "sampling/exponential_sampling.hpp"
#include "sampling/chebyshev_lobatto.hpp"
#include "sampling/equispaced.hpp"
#include "sampling/exponential_refinement.hpp"

#include "polynomial/polynomial_point_value.hpp"
#include "polynomial/chebyshev_polynomial.hpp"
#include "polynomial/multivariate_polynomial_point_value.hpp"

#include "smolyak/smolyak_approximation.hpp"
#include "smolyak/smolyak_traits.hpp"

namespace Analytic
{
    using Arithmetic::index_t;
    using Arithmetic::integer_t;
    using Arithmetic::real_t;
    using Arithmetic::complex_t;
    using Arithmetic::negative_1;

    // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //

    constexpr real_t rho(real_t epsilon)
    {
        ASSERT_ASSUME(0 < epsilon and epsilon < 0.5);

        epsilon = 1 / (1 - epsilon);
        return epsilon + std::sqrt(epsilon * epsilon - 1);
    }

    auto f = []<index_t n_variates>(const std::array<real_t, n_variates> &arg) -> real_t
    {
        real_t ret = 1;

        for (real_t i : arg)
            ret *= (1 - i * i);

        ret = std::sqrt(ret);

        return ret;
    };

    // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //
    // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //

    using sampling_type = Sampling::exponential<Sampling::chebyshev_lobatto>;
    using polynomial_type = Polynomial::polynomial_point_value<sampling_type>;
    using interpolation_traits = Smolyak::Smolyak_traits::polynomial_interpolation<sampling_type>;
    using fidelity_traits = Smolyak::Smolyak_traits::multifidelity<>;
    using richardson_traits =
            Smolyak::Smolyak_traits::multifidelity_richardson_extrapolation<0.5>;
}

// **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //

