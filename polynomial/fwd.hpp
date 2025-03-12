#pragma once

#include "core/core.hpp"

#include "../sampling/chebyshev_lobatto.hpp"
#include "../sampling/exponential_sampling.hpp"

// **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //

namespace Polynomial
{
    using Arithmetic::index_t;
    using Arithmetic::integer_t;
    using Arithmetic::real_t;
    using Arithmetic::complex_t;
    using Arithmetic::negative_1;

    // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //
    // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //

    namespace Apparatus
    {
        template<class type>
        struct multivariate_polynomial_c_apparatus;
    }

    template<class polynomial_type>
    concept polynomial_c =
    requires()
    {
        requires polynomial_type::polynomial_tag;
        requires std::is_same_v<decltype(polynomial_type::n_variates), const index_t>;
    };

    template<class polynomial_type>
    concept univariate_polynomial_c =
    requires()
    {
        requires polynomial_c<polynomial_type>;
        requires polynomial_type::n_variates == 1;
    };

    template<class polynomial_type>
    concept multivariate_polynomial_c =
            Apparatus::multivariate_polynomial_c_apparatus<polynomial_type>::value;

    // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //
    // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //

    namespace Apparatus
    {
        template<class type>
        class tag {};

        template<class polynomial_traits, bool is_proxy = false>
        class polynomial_base;

        template<class coefficients_type>
        struct polynomial_traits;

        template<class recurrence_relation_traits_>
        struct polynomial_recurrence_relation_traits
        {
            using recurrence_relation_traits = recurrence_relation_traits_;
            using coefficients_type = recurrence_relation_traits::coefficients_type;
        };

        template<class recurrence_relation_traits>
        using polynomial_recurrence_relation = polynomial_base
                <polynomial_recurrence_relation_traits<recurrence_relation_traits>>;

        template<class coefficients_type>
        struct legendre_recurrence_relation;

        template<index_t kind, class coefficients_type = real_t>
        struct chebyshev_recurrence_relation;

        template<real_t alpha, real_t beta, class coefficients_type = real_t>
        struct jacobi_recurrence_relation;
    }

    template<class coefficients_type>
    using polynomial = Apparatus::polynomial_base<Apparatus::polynomial_traits<coefficients_type>>;

    template<class coefficients_type = real_t>
    using legendre_polynomial = Apparatus::polynomial_base
            <Apparatus::polynomial_recurrence_relation_traits<Apparatus::legendre_recurrence_relation
                    <coefficients_type>>>;

    template<class coefficients_type = real_t, index_t kind = 1>
    using chebyshev_polynomial = Apparatus::polynomial_recurrence_relation
            <Apparatus::chebyshev_recurrence_relation<kind, coefficients_type>>;

    using chebyshev_polynomial_i_kind = chebyshev_polynomial<real_t, 1>;
    using chebyshev_polynomial_ii_kind = chebyshev_polynomial<real_t, 2>;
    using chebyshev_polynomial_iii_kind = chebyshev_polynomial<real_t, 3>;
    using chebyshev_polynomial_iv_kind = chebyshev_polynomial<real_t, 4>;

    template<real_t alpha, real_t beta, class coefficients_type = real_t>
    requires(alpha > -1 and beta > -1)
    using jacobi_polynomial = Apparatus::polynomial_base
            <Apparatus::polynomial_recurrence_relation_traits
                    <Apparatus::jacobi_recurrence_relation<alpha, beta, coefficients_type>>>;

    template<Sampling::sampling_c sampling_type,
             class codomain_field = sampling_type::coordinates_type,
             class container_type = std::vector<codomain_field>>
    class polynomial_point_value;

    template<univariate_polynomial_c, index_t>
    class multivariate_polynomial;

    using chebyshev_polynomial_point_value = polynomial_point_value
            <Sampling::exponential<Sampling::chebyshev_lobatto>>;

    namespace Apparatus
    {
        template<class type>
        struct multivariate_polynomial_c_apparatus
        {
            static constexpr bool value = false;
        };

        template<univariate_polynomial_c univariate_polynomial_type, index_t n_variates>
        struct multivariate_polynomial_c_apparatus<multivariate_polynomial
                <univariate_polynomial_type, n_variates>>
        {
            static constexpr bool value = true;
        };
    }
}
