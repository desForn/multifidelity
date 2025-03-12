#pragma once

#include "fwd.hpp"

namespace Polynomial
{
    void discrete_fourier_transform(std::span<complex_t>);
    void inverse_discrete_fourier_transform(std::span<complex_t>);

    void multidimensional_discrete_fourier_transform
            (std::span<complex_t>, std::span<const index_t>);
    void inverse_multidimensional_discrete_fourier_transform
            (std::span<complex_t>, std::span<const index_t>);

    chebyshev_polynomial_i_kind discrete_cosine_transform(const chebyshev_polynomial_point_value &);
    chebyshev_polynomial_point_value discrete_cosine_transform(const chebyshev_polynomial_i_kind &);

    template<index_t n_variates>
    multivariate_polynomial<chebyshev_polynomial_i_kind, n_variates> discrete_cosine_transform
            (const multivariate_polynomial<chebyshev_polynomial_point_value, n_variates> &);

    template<index_t n_variates>
    multivariate_polynomial<chebyshev_polynomial_point_value, n_variates> discrete_cosine_transform
            (const multivariate_polynomial<chebyshev_polynomial_i_kind, n_variates> &);
}
