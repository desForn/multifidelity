#pragma once

#include "multivariate_polynomial.hpp"

// **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //

/* Implementation of polynomials in the canonical representation:
 *      p = sum_{j=0}^n c_j X^j */

namespace Polynomial
{
    namespace Apparatus
    {
        template<class coefficients_type_>
        struct polynomial_traits
        {
            using coefficients_type = coefficients_type_;
        };

        // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //
        // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //

        template<class coefficients_type, bool is_proxy, class variate_type>
        auto evaluate
                (const polynomial_base
                        <polynomial_traits<coefficients_type>, is_proxy> &p,
                 const variate_type &x)
                 -> decltype(std::declval<coefficients_type>() * x)
        {
            // Horner's method
            std::span<const coefficients_type> coefficients = p.coefficients();

            if (p.empty())
                return static_cast<variate_type>(0);

            variate_type ret = static_cast<variate_type>(coefficients.back());

            for (auto it = std::crbegin(coefficients) + 1; it != std::crend(coefficients); ++it)
                ret = *it + ret * x;

            return ret;
        }

        template<class coefficients_type, bool is_proxy_lhs, bool is_proxy_rhs>
        polynomial<coefficients_type> operator*
                (const polynomial_base<polynomial_traits<coefficients_type>,
                        is_proxy_lhs> &lhs,
                 const polynomial_base<polynomial_traits<coefficients_type>,
                         is_proxy_rhs> &rhs)
        {
            if (lhs.empty() or rhs.empty())
                return {};

            index_t lhs_extents = lhs.extents();
            index_t rhs_extents = rhs.extents();

            index_t ret_extents = lhs_extents + rhs_extents - 1;

            polynomial<coefficients_type> ret;
            ret.set_extents(ret_extents);

            for (index_t i = 0; i != ret_extents; ++i)
                for (index_t j = (i < rhs_extents ? 0 : i - rhs_extents + 1);
                     j != std::min(i + 1, lhs_extents); ++j)
                    ret[i] += lhs[j] * rhs[i - j];

            return ret;
        }

        template<class coefficients_type, class variate_type>
        auto evaluate_basis_functions
        (const variate_type &x, index_t n, tag<polynomial<coefficients_type>>)
        -> std::vector<decltype(std::declval<coefficients_type>() * x)>
        {
            using ret_type = decltype(std::declval<coefficients_type>() * x);

            std::vector<ret_type> ret;

            if (n == 0)
                return ret;

            ret.reserve(n);
            ret.emplace_back(static_cast<coefficients_type>(1));

            for(index_t i = 1; i != n; ++i)
                ret.emplace_back(ret.back() * x);

            return ret;
        }
    }

    template<index_t n_variates, class coefficients_type>
    multivariate_polynomial<polynomial<coefficients_type>, n_variates> operator*
            (const multivariate_polynomial<polynomial<coefficients_type>, n_variates> &lhs,
             const multivariate_polynomial<polynomial<coefficients_type>, n_variates> &rhs)
    {
        using polynomial_type = multivariate_polynomial<polynomial<coefficients_type>, n_variates>;

        if (lhs.empty() or rhs.empty())
            return {};

        std::array<index_t, n_variates> lhs_extents = lhs.extents();
        std::array<index_t, n_variates> rhs_extents = rhs.extents();
        std::array<index_t, n_variates> ret_extents =
                Utility::subtract_array(Utility::sum_array(lhs_extents, rhs_extents), 1);

        polynomial_type ret;
        ret.set_extents(ret_extents);

        for (const std::array<index_t, n_variates> &i: Utility::counter{lhs_extents})
            for (const std::array<index_t, n_variates> &j: Utility::counter{rhs_extents})
                ret[Utility::sum_array(i, j)] += lhs[i] * rhs[j];

        return ret;
    }
}
