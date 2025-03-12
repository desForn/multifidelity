#pragma once

#include "polynomial_canonical_basis.hpp"

// **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //

/* Here it is implemented the Clenshaw algorithm [https://en.wikipedia.org/wiki/Clenshaw_algorithm]
 * as a way to evaluate polynomials basis that are expressed recurrently.
 *      Note that we have modified the wikipedia notation:
 *      alpha is substituted by gamma and beta by delta.
 *      This is to avoid conflicts with alpha and beta from jacobi polynomials.
 * The formulation is as follows:
 *      phi_0 = 1,
 *      phi_1 = base_case_10() + base_case_11() * X,
 *
 *      for all k > 1:
 *      phi_k = (gamma_0(k - 1) + gamma_1(k - 1) * X) * phi_{k - 1} + delta(k - 1) * phi_{k - 2};
 *
 * All functions in the previos formulae are provided as public methods of the template argument
 * class. If "gamma_0" or "base_case_10" should always return 0, one can set its return type to
 * "void". See "jacobi_polynomial.hpp" for an example. */

namespace Polynomial::Apparatus
{
    template<class recurrence_relation_traits, bool is_proxy, class variate_type>
    auto evaluate(const polynomial_base<polynomial_recurrence_relation_traits
                  <recurrence_relation_traits>, is_proxy> &p,
                  const variate_type &x)
    -> decltype(std::declval<typename recurrence_relation_traits::coefficients_type>() * x)
    {
        using ret_type = decltype
        (std::declval<typename recurrence_relation_traits::coefficients_type>() * x);

        if (p.empty())
            return static_cast<ret_type>(0);

        ret_type b0{0}, b1{0}, b2{0};

        for (index_t i = p.extents() - 1; i != 0; b2 = b1, b1 = b0, --i)
        {
            b0 = static_cast<ret_type>(recurrence_relation_traits::gamma_1(i) * x);

            if constexpr (!std::is_void_v<decltype(recurrence_relation_traits::gamma_0(i))>)
                b0 += recurrence_relation_traits::gamma_0(i);

            b0 *= b1;

            if constexpr (!std::is_void_v<decltype(recurrence_relation_traits::delta(i))>)
                b0 += recurrence_relation_traits::delta(i + 1) * b2;

            b0 += p[i];
        }

        /* After termination the value that should store b1 (in wikipedia notation), is in fact
         * stored at b0 */

        ret_type ret = static_cast<ret_type>(recurrence_relation_traits::base_case_11() * x);

        if constexpr (!std::is_void_v<decltype(recurrence_relation_traits::base_case_10())>)
            ret += recurrence_relation_traits::base_case_10();

        ret *= b0;

        if constexpr (!std::is_void_v<decltype(recurrence_relation_traits::delta(1))>)
            ret += recurrence_relation_traits::delta(1) * b2;

        ret += p[0];

        return ret;
    }

    // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //
    // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //

    template<class recurrence_relation_traits>
    class conversion_coefficients<polynomial_recurrence_relation<recurrence_relation_traits>,
            polynomial<typename recurrence_relation_traits::coefficients_type>>
    {
    public:
        using coefficients_type = recurrence_relation_traits::coefficients_type;
        static constexpr bool implemented = true;
        static coefficients_type evaluate(index_t n)
        {
            if (n == 0)
                return static_cast<coefficients_type>(1);

            if (n == 1)
            {
                if constexpr (!std::is_void_v<decltype
                (recurrence_relation_traits::base_case_10())>)
                    return recurrence_relation_traits::base_case_10();

                else
                    return static_cast<coefficients_type>(0);
            }

            if (n == 2)
                return recurrence_relation_traits::base_case_11();

            /* if i == 0, then j * j + j - 2 * n == 0, so j == (sqrt(8 * n + 1) - 1) / 2
             * if i != 0, as i <= j, then j == trunc((sqrt(8 * n + 1) - 1) / 2)
             * to avoid rounding errors issues, substitute n by n + 0.5 */

            index_t j = static_cast<index_t>
                    ((std::sqrt(8 * static_cast<real_t>(n + 0.5) + 1) - 1) / 2);

            index_t i = n - j * (j + 1) / 2;

            coefficients_type ret = static_cast<coefficients_type>(0);

            const auto &cache_invocable = forward_conversion
                    <polynomial_recurrence_relation<recurrence_relation_traits>,
                            polynomial<coefficients_type>>::invocable_;

            /* The relations above come from the following facts:
             * vector_index(i, j) - vector_index(i - 1, j) == 1
             * vector_index(i, j) - vector_index(i, j - 1) == j */

            if constexpr (!std::is_void_v<decltype(recurrence_relation_traits::delta(i))>)
                if (j - i >= 2) // vector_index(i, j - 2) == n - 2 * j + 1
                    ret += recurrence_relation_traits::delta(j - 1)
                           * cache_invocable(n - 2 * j + 1);

            if constexpr (!std::is_void_v<decltype(recurrence_relation_traits::gamma_0(i))>)
                if (j != i) // vector_index(i, j - 1) == n - j
                    ret += recurrence_relation_traits::gamma_0(j - 1)
                           * cache_invocable(n - j);

            if (i != 0) // vector_index(i - 1, j - 1) == n - j - 1
                ret += recurrence_relation_traits::gamma_1(j - 1)
                       * cache_invocable(n - j - 1);

            return ret;
        }
    };
}
