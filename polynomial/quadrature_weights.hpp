#pragma once

#include "fourier_transform.hpp"

// **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //
// **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //

namespace Polynomial
{
    template<>
    std::vector<real_t> compute_quadrature_weights
        <Sampling::exponential<Sampling::chebyshev_lobatto>>(index_t sampling_level)
    {
        if (sampling_level == static_cast<index_t>(-1))
            return {};

        if (sampling_level == 0)
            return {2};

        integer_t n = Sampling::n_points_exponential_sampling(sampling_level);

        /* p[i] == int_{-1}^{1} chebyshev_polynomial_i_kind::basis_polynomial(i);
         * w is the result of applying the adjoint of the discrete cosine transform to p */

        chebyshev_polynomial_i_kind p;
        p.set_extents(n);

        p[0] = 2;

        // Note the extra factor of 2 needed for normalisation of the adjoint operator
        for (integer_t i = 2; i < n - 1; i += 2)
            p[i] = static_cast<real_t>(4) / static_cast<real_t>(1 - i * i);

        p[n - 1] = static_cast<real_t>(2) / static_cast<real_t>(1 - (n - 1) * (n - 1));

        chebyshev_polynomial_point_value w = discrete_cosine_transform(p);

        std::vector<real_t> ret = std::move(w).f();

        index_t unfolded_n = 2 * (n - 1);

        ret.front() /= unfolded_n;

        std::for_each(std::begin(ret) + 1, std::end(ret) - 1,
                      [unfolded_n](real_t & i) { i *= static_cast<real_t>(2) / unfolded_n; });

        ret.back() /= unfolded_n;

        return ret;
    }
}
