#pragma once

#include "fwd.hpp"

// **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //

namespace Miscellany
{
    template<Arithmetic::real_complex_c lhs_type, index_t lhs_n,
             Arithmetic::real_complex_c rhs_type, index_t rhs_n,
             Arithmetic::real_complex_c z_type>
    std::common_type_t<lhs_type, rhs_type, z_type> hypergeometric_function
            (const std::array<lhs_type, lhs_n> &lhs, const std::array<rhs_type, rhs_n> &rhs,
             z_type z, real_t tol = 0)
    {
        /* https://en.wikipedia.org/wiki/Hypergeometric_function
         * The behaviour is only well-defined if the return value should be a finite number.
         * Note the introduction of "if (num == 0) break;" to enable simpler handling of
         * singular cases */

        using ret_type = std::common_type_t<lhs_type, rhs_type, z_type>;

        ASSERT_ASSUME(tol >= 0);

        auto less = Arithmetic::less(tol);

        ret_type ret_prev = 0;
        ret_type ret = 1;
        real_t n_factorial = 1;
        lhs_type num = 1;
        rhs_type den = 1;
        z_type z_power = 1;

        for (index_t n = 1; less(std::abs(ret_prev * tol), std::abs(ret_prev - ret)); ++n)
        {
            for (const lhs_type &i: lhs)
                num *= i + (n - 1);

            if (num == 0)
                break;

            for (const rhs_type &i: rhs)
                den *= i + (n - 1);

            n_factorial *= n;

            z_power *= z;

            ret_prev = ret;
            ret += (num * z_power) / (den * n_factorial);
        }

        return ret;
    }
}
