#pragma once

#include "polynomial_recurrence_relation.hpp"

// **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //

// https://en.wikipedia.org/wiki/Legendre_polynomials

namespace Polynomial
{
    namespace Apparatus
    {
        template<class coefficients_type_ = real_t>
        struct legendre_recurrence_relation
        {
            using coefficients_type = coefficients_type_;

            static void gamma_0(index_t) {}
            static real_t gamma_1(index_t);
            static real_t delta(index_t);

            static void base_case_10() {}
            static real_t base_case_11();
        };

        // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //

        template<class coefficients_type>
        real_t legendre_recurrence_relation<coefficients_type>::gamma_1(index_t i)
        {
            return 2 - 1 / static_cast<real_t>(i + 1);
        }

        template<class coefficients_type>
        real_t legendre_recurrence_relation<coefficients_type>::delta(index_t i)
        {
            return -1 + 1 / static_cast<real_t>(i + 1);
        }

        template<class coefficients_type>
        real_t legendre_recurrence_relation<coefficients_type>::base_case_11()
        {
            return 1;
        }

        // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //
        // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //

        template<class coefficients_type, bool is_proxy_lhs, bool is_proxy_rhs>
        legendre_polynomial<coefficients_type> operator*
                (const polynomial_base<polynomial_recurrence_relation_traits
                        <legendre_recurrence_relation<coefficients_type>>, is_proxy_lhs> &lhs,
                 const polynomial_base<polynomial_recurrence_relation_traits
                         <legendre_recurrence_relation<coefficients_type>>, is_proxy_rhs> &rhs)
        {
            if (lhs.empty() or rhs.empty())
                return {};

            index_t lhs_extents = lhs.extents();
            index_t rhs_extents = rhs.extents();

            index_t ret_extents = lhs_extents + rhs_extents - 1;

            legendre_polynomial<coefficients_type> ret;
            ret.set_extents(ret_extents);

            index_t max_order = std::max(lhs_extents, rhs_extents) - 1;
            index_t min_order = std::min(lhs_extents, rhs_extents) - 1;

            index_t ret_order = ret_extents - 1;

            std::vector<real_t> A(ret_order + 1);

            A[0] = 1;

            for (index_t i = 1; i != ret_order + 1; ++i)
                A[i] = A[i - 1] * static_cast<real_t>(2 * i - 1) / static_cast<real_t>(i);

            for (index_t i = 0; i <= max_order; ++i)
                for (index_t j = 0; j <= std::min(i, min_order); ++j)
                    for (index_t k = 0; k <= j; ++k)
                    {
                        coefficients_type c = static_cast<coefficients_type>(0);

                        bool skip = true;

                        if (i < lhs_extents and j < rhs_extents)
                        {
                            c = lhs[i] * rhs[j];
                            skip = false;
                        }

                        if (i != j and j < lhs_extents and i < rhs_extents)
                        {
                            c += lhs[j] * rhs[i];
                            skip = false;
                        }

                        if (skip)
                            continue;

                        ret[i + j - 2 * k] +=
                                A[i - k] * A[k] * A[j - k] * (2 * i + 2 * j - 4 * k + 1) /
                                (A[i + j - k] * (2 * i + 2 * j - 2 * k + 1)) * c;
                    }

            return ret;
        }
    }

    // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //
    // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //

    template<class coefficients_type, bool is_proxy_lhs, bool is_proxy_rhs>
    requires((std::is_same_v<coefficients_type, real_t> or
              std::is_same_v<coefficients_type, complex_t>))
    coefficients_type inner_product
            (const Apparatus::polynomial_base<Apparatus::polynomial_recurrence_relation_traits
                    <Apparatus::legendre_recurrence_relation<coefficients_type>>,
                    is_proxy_lhs> &lhs,
             const Apparatus::polynomial_base<Apparatus::polynomial_recurrence_relation_traits
                     <Apparatus::legendre_recurrence_relation<coefficients_type>>,
                     is_proxy_rhs> &rhs)
    {
        coefficients_type ret = 0;

        for (index_t i = 0; i != std::min(lhs.extents(), rhs.extents()); ++i)
            ret += lhs[i] * rhs[i] * static_cast<real_t>(2) / static_cast<real_t>(2 * i + 1);

        return ret;
    }

    template<class coefficients_type, bool is_proxy>
    requires(std::is_same_v<coefficients_type, real_t> and
             std::is_same_v<coefficients_type, complex_t>)
    real_t norm2
            (const Apparatus::polynomial_base<Apparatus::polynomial_recurrence_relation_traits
                    <Apparatus::legendre_recurrence_relation<coefficients_type>>, is_proxy> &rhs)
    {
        if constexpr (std::is_same_v<coefficients_type, real_t>)
            return inner_product(rhs, rhs);

        else
            return inner_product(rhs, rhs).real();
    }

    template<class coefficients_type, bool is_proxy>
    requires(std::is_same_v<coefficients_type, real_t> and
             std::is_same_v<coefficients_type, complex_t>)
    real_t norm
            (const Apparatus::polynomial_base<Apparatus::polynomial_recurrence_relation_traits
                    <Apparatus::legendre_recurrence_relation<coefficients_type>>, is_proxy> &rhs)
    {
        return std::sqrt(norm2(rhs));
    }

    template<class coefficients_type, bool is_proxy>
    coefficients_type integrate
            (const Apparatus::polynomial_base<Apparatus::polynomial_recurrence_relation_traits
                    <Apparatus::legendre_recurrence_relation<coefficients_type>>, is_proxy> &rhs)
    {
        if (rhs.empty())
            return static_cast<coefficients_type>(0);

        return rhs[0] * 2;
    }
}
