#pragma once

#include "polynomial_recurrence_relation.hpp"
#include "../miscellany/hypergeometric_function.hpp"

// **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //

//https://en.wikipedia.org/wiki/Jacobi_polynomials

namespace Polynomial
{
    namespace Apparatus
    {
        template<real_t alpha, real_t beta, class coefficients_type_>
        struct jacobi_recurrence_relation
        {
            using coefficients_type = coefficients_type_;

            static std::conditional_t<alpha == beta, void, real_t> gamma_0(index_t);
            static real_t gamma_1(index_t);
            static real_t delta(index_t);

            static std::conditional_t<alpha == beta, void, real_t> base_case_10();
            static real_t base_case_11();
        };

        // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //

        template<real_t alpha, real_t beta, class coefficients_type>
        std::conditional_t<alpha == beta, void, real_t>
        jacobi_recurrence_relation<alpha, beta, coefficients_type>::gamma_0(index_t i)
        {
            if constexpr (alpha == beta)
                return;

            else
            {
                if constexpr (alpha + beta == 0)
                {
                    if (i == 0)
                        return (alpha - beta) / 2;

                    return 0;
                }

                if constexpr (alpha + beta + 1 == 0)
                    return (alpha + beta) * (alpha - beta) / ((2 * i - 1) * (i + 1));

                return (alpha + beta) * (alpha - beta) * (2 * i + alpha + beta + 1)
                       / (2 * (i + 1) * (i + alpha + beta + 1) * (2 * i + alpha + beta));
            }
        }

        template<real_t alpha, real_t beta, class coefficients_type>
        real_t jacobi_recurrence_relation<alpha, beta, coefficients_type>::gamma_1(index_t i)
        {
            if constexpr (alpha + beta == 0 or alpha + beta + 1 == 0)
                return 2 - 1 / static_cast<real_t>(i + 1);

            return (2 * i + alpha + beta + 1) * (2 * i + alpha + beta + 2)
                   / (2 * (i + 1) * (i + alpha + beta + 1));
        }

        template<real_t alpha, real_t beta, class coefficients_type>
        real_t jacobi_recurrence_relation<alpha, beta, coefficients_type>::delta(index_t i)
        {
            ASSERT_ASSUME(i != 0);

            if constexpr (alpha == 0 and beta == 0)
                return -1 + 1 / static_cast<real_t>(i + 1);

            return -(i + alpha) * (i + beta) * (2 * i + alpha + beta + 2)
                   / ((i + 1) * (i + alpha + beta + 1) * (2 * i + alpha + beta));
        }

        template<real_t alpha, real_t beta, class coefficients_type>
        std::conditional_t<alpha == beta, void, real_t>
        jacobi_recurrence_relation<alpha, beta, coefficients_type>::base_case_10()
        {
            if constexpr (alpha == beta)
                return;

            return 0.5 * (alpha - beta);
        }

        template<real_t alpha, real_t beta, class coefficients_type>
        real_t jacobi_recurrence_relation<alpha, beta, coefficients_type>::base_case_11()
        {
            return 1 + 0.5 * (alpha + beta);
        }

        // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //
        // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //

        template<real_t alpha, real_t beta,
                 class coefficients_type, bool is_proxy_lhs, bool is_proxy_rhs>
        jacobi_polynomial<alpha, beta> operator*
                (const polynomial_base<polynomial_recurrence_relation_traits
                        <jacobi_recurrence_relation<alpha, beta, coefficients_type>>,
                        is_proxy_lhs> &lhs,
                 const polynomial_base<polynomial_recurrence_relation_traits
                         <jacobi_recurrence_relation<alpha, beta, coefficients_type>>,
                         is_proxy_rhs> &rhs)
        {
            // doi.org/10.1186/s13662-015-0509-4 -- https://rdcu.be/dr0yH (eqs. 3, 9, 10)

            constexpr real_t ab1 = alpha + beta + 1;

            if (lhs.empty() or rhs.empty())
                return {};

            index_t lhs_extents = lhs.extents();
            index_t rhs_extents = rhs.extents();

            index_t ret_extents = lhs_extents + rhs_extents - 1;

            jacobi_polynomial<alpha, beta, coefficients_type> ret;
            ret.set_extents(ret_extents);

            index_t max_order = std::max(lhs_extents, rhs_extents) - 1;
            index_t min_order = std::min(lhs_extents, rhs_extents) - 1;

            std::vector<real_t>
                    g_a(max_order + min_order + 1);  // g_a[i] == std::tgamma(i + alpha + 1)
            std::vector<real_t>
                    g_b(max_order + min_order + 1);  // g_b[i] == std::tgamma(i + beta + 1)
            std::vector<real_t>
                    g_ab(2 * max_order + 2);         // g_ab[i] == std::tgamma(i + alpha + beta + 1)
            std::vector<real_t>
                    p_a(max_order + min_order + 1);  // p_a[i] = pochhammer(alpha + 1, i)

            g_a[0] = std::tgamma(alpha + 1);
            for (index_t i = 1; i != std::size(g_a); ++i)
                g_a[i] = (alpha + i) * g_a[i - 1];

            g_b[0] = std::tgamma(beta + 1);
            for (index_t i = 1; i != std::size(g_b); ++i)
                g_b[i] = (beta + i) * g_b[i - 1];

            // alpha + beta + 1 may be 0
            // In this case g_ab[0] is undefined, but it is well handed below
            g_ab[1] = std::tgamma(alpha + beta + 2);

            if constexpr (ab1 != 0)
                g_ab[0] = g_ab[1] / (alpha + beta + 1);

            for (index_t i = 2; i != std::size(g_ab); ++i)
                g_ab[i] = (alpha + beta + i) * g_ab[i - 1];

            p_a[0] = 1;
            for (index_t i = 1; i != std::size(p_a); ++i)
                p_a[i] = (alpha + i) * p_a[i - 1];

            for (index_t i = 0, factorial_i = 1; i <= max_order; factorial_i *= ++i)

                // factorial_imj == factorial(i - j)
                for (index_t j = 0, factorial_j = 1, factorial_imj = factorial_i;
                     j <= std::min(i, min_order);
                     factorial_imj /= (i != j) ? (i - j) : 1, factorial_j *= ++j)

                    // factorial_imjpk == factorial(i - j + k)
                    for (index_t k = 0, factorial_k = 1, factorial_imjpk = factorial_imj;
                         k <= 2 * j;
                         factorial_k *= ++k, factorial_imjpk *= (i - j + k))
                    {
                        index_t imjpk = i - j + k;

                        /* If one of {i, j, imjpk} is 0 and the two others disctinct,
                         * the contribution is 0. This is by orthogonality of the Jacobi
                         * polynomials. */

                        if (std::min(j, imjpk) == 0 and i + j + imjpk != 2 * std::max(i, imjpk))
                            continue;

                        real_t c = 0;

                        if (i < lhs_extents and j < rhs_extents)
                            c = lhs[i] * rhs[j];

                        if (i != j and j < lhs_extents and i < rhs_extents)
                            c += lhs[j] * rhs[i];

                        if (c == 0)
                            continue;

                        real_t a = 0;

                        real_t ri = i;
                        real_t rj = j;
                        real_t rk = k;

                        index_t i0 = 2 * i - 2 * j + k;
                        real_t r0 = i0;

                        real_t p0 = 1; // pochhammer(-rk, l)
                        real_t p1 = 1; // pochhammer(rk + 2 * ri - 2 * rj + alpha + beta + 1, l)
                        real_t p2 = 1; // pochhammer(2 * ri + alpha + beta + 2, l)

                        for (index_t l = 0, factorial_l = 1;
                             l <= k;
                             p0 *= (-rk + l),
                             p1 *= (r0 + ab1 + l),
                             p2 *= (2 * ri + ab1 + 1 + l),
                             factorial_l *= ++l)
                        {
                            real_t sum_a = (p0 * p1) / (factorial_l * p2);

                            if (l != 0)
                                sum_a *= Miscellany::hypergeometric_function
                                        (std::array<real_t, 4>
                                                 {-static_cast<real_t>(l), l + 2 * ri - 2 * rj + 1,
                                                  -rj, -rj - beta},
                                         std::array<real_t, 3>
                                                 {ri - rj + 1, ri - rj + alpha + 1,
                                                  -2 * rj - alpha - beta},
                                         real_t{1});

                            a += sum_a;
                        }

                        if (a == 0)
                            continue;

                        a *= factorial_i * g_a[0];

                        // lim (z -> 0) (z * std::tgamma(z)) == 1
                        if (ab1 != 0 or r0 != 0)
                            a *= (r0 + rk + ab1) * g_ab[i0];

                        a /= factorial_imj * factorial_k * g_ab[2 * i + 1];
                        a *= g_a[imjpk] * g_b[i];
                        a /= g_b[imjpk] * g_a[i - j] * g_a[j];

                        if (j != 0)
                            a *= g_ab[2 * j] / g_ab[j];

                        a *= factorial_imjpk / p_a[imjpk];
                        a *= p_a[i] / factorial_i;
                        a *= p_a[j] / factorial_j;

                        ret[imjpk] += a * c;
                    }

            return ret;
        }
    }

    // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //
    // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //

    template<real_t alpha, real_t beta,
             class coefficients_type, bool is_proxy_lhs, bool is_proxy_rhs>
    requires((std::is_same_v<coefficients_type, real_t> or
              std::is_same_v<coefficients_type, complex_t>))
    coefficients_type inner_product
            (const Apparatus::polynomial_base<Apparatus::polynomial_recurrence_relation_traits
                    <Apparatus::jacobi_recurrence_relation<alpha, beta, coefficients_type>>,
                    is_proxy_lhs> &lhs,
             const Apparatus::polynomial_base<Apparatus::polynomial_recurrence_relation_traits
                     <Apparatus::jacobi_recurrence_relation<alpha, beta, coefficients_type>>,
                     is_proxy_rhs> &rhs)
    {
        index_t n = std::min(lhs.extents(), rhs.extents());

        coefficients_type ret = static_cast<coefficients_type>(0);

        if (n == 0)
            return ret;

        constexpr real_t ab1 = alpha + beta + 1;
        constexpr real_t pow = std::pow(2, ab1);

        real_t g_a = std::tgamma(alpha + 1);
        real_t g_b = std::tgamma(beta + 1);
        real_t g_ab;

        index_t i = 0;
        index_t i_factorial = 1;

        if constexpr (ab1 != 0)
            g_ab = std::tgamma(ab1);

        else
        {
            ret = lhs[0] * rhs[0] * pow * g_a * g_b;

            g_a *= alpha + 1;
            g_b *= beta + 1;
            g_ab = std::tgamma(ab1 + 1);
        }

        for (; i != n;
               g_a *= alpha + 1 + i,
               g_b *= beta + 1 + i,
               g_ab *= ab1 + i,
               i_factorial *= ++i)
            ret += lhs[i] * rhs[i] * pow / (2 * i + ab1) * g_a * g_b / (g_ab * i_factorial);

        return ret;
    }

    template<real_t alpha, real_t beta,
             class coefficients_type, bool is_proxy_rhs>
    real_t norm2(const Apparatus::polynomial_base<Apparatus::polynomial_recurrence_relation_traits
            <Apparatus::jacobi_recurrence_relation<alpha, beta, coefficients_type>>, is_proxy_rhs> &rhs)
    requires((std::is_same_v<coefficients_type, real_t> or
              std::is_same_v<coefficients_type, complex_t>))
    {
        if constexpr (std::is_same_v<coefficients_type, real_t>)
            return inner_product(rhs, rhs);

        else
            return inner_product(rhs, rhs).real();
    }

    template<real_t alpha, real_t beta,
             class coefficients_type, bool is_proxy_rhs>
    real_t norm(const Apparatus::polynomial_base<Apparatus::polynomial_recurrence_relation_traits
            <Apparatus::jacobi_recurrence_relation<alpha, beta, coefficients_type>>, is_proxy_rhs> &rhs)
    requires((std::is_same_v<coefficients_type, real_t> or
              std::is_same_v<coefficients_type, complex_t>))
    {
        return std::sqrt(norm2(rhs));
    }
}
