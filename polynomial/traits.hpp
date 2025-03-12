#pragma once

#include "fwd.hpp"

// **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //

namespace Polynomial
{
    namespace Apparatus
    {
        template<polynomial_c>
        struct univariate_polynomial_type_apparatus;

        template<class polynomials_traits>
        struct univariate_polynomial_type_apparatus<polynomial_base<polynomials_traits>>
        {
            using type = polynomial_base<polynomials_traits>;
        };

        template<index_t n_variates, univariate_polynomial_c univariate_polynomial_type>
        struct univariate_polynomial_type_apparatus
                <multivariate_polynomial<univariate_polynomial_type, n_variates>>
        {
            using type = univariate_polynomial_type;
        };

        // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //

        template<polynomial_c>
        struct polynomial_proxy_type_apparatus;

        template<class polynomial_traits>
        struct polynomial_proxy_type_apparatus<polynomial_base<polynomial_traits, true>>
        {
            using type = polynomial_base<polynomial_traits, true>;
        };

        template<class polynomial_traits>
        struct polynomial_proxy_type_apparatus<polynomial_base<polynomial_traits, false>>
        {
            using type = polynomial_base<polynomial_traits, true>;
        };

        template<class univariate_polynomial, index_t n_variates>
        struct polynomial_proxy_type_apparatus
                <multivariate_polynomial<univariate_polynomial, n_variates>>
        {
            using type = multivariate_polynomial
                    <typename polynomial_proxy_type_apparatus<univariate_polynomial>::type, n_variates>;
        };

        template<polynomial_c polynomial_type>
        struct polynomial_remove_proxy_type_apparatus
        {
            using type = polynomial_type;
        };

        template<class polynomial_traits>
        struct polynomial_remove_proxy_type_apparatus<polynomial_base<polynomial_traits, true>>
        {
            using type = polynomial_base<polynomial_traits, false>;
        };

        template<class polynomial_traits>
        struct polynomial_remove_proxy_type_apparatus<polynomial_base<polynomial_traits, false>>
        {
            using type = polynomial_base<polynomial_traits, false>;
        };

        template<class univariate_polynomial, index_t n_variates>
        struct polynomial_remove_proxy_type_apparatus
                <multivariate_polynomial<univariate_polynomial, n_variates>>
        {
            using type = multivariate_polynomial
                    <typename polynomial_remove_proxy_type_apparatus<univariate_polynomial>::type,
                     n_variates>;
        };

        // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //

        template<polynomial_c, class>
        struct convert_polynomial_coefficients_type_apparatus;

        template<index_t n_variates, univariate_polynomial_c univariate_polynomial_type,
                class destination_coefficients_type>
        struct convert_polynomial_coefficients_type_apparatus<multivariate_polynomial
                <univariate_polynomial_type, n_variates>, destination_coefficients_type>
        {
            using type = multivariate_polynomial<typename
                    convert_polynomial_coefficients_type_apparatus
                    <univariate_polynomial_type, destination_coefficients_type>::type, n_variates>;
        };

        template<class origin_coefficients_type, class destination_coefficients_type>
        struct convert_polynomial_coefficients_type_apparatus
                <polynomial<origin_coefficients_type>, destination_coefficients_type>
        {
            using type = polynomial<destination_coefficients_type>;
        };

        template<class origin_coefficients_type, class destination_coefficients_type>
        struct convert_polynomial_coefficients_type_apparatus
                <legendre_polynomial<origin_coefficients_type>, destination_coefficients_type>
        {
            using type = legendre_polynomial<destination_coefficients_type>;
        };

        template<index_t kind, class origin_coefficients_type, class destination_coefficients_type>
        struct convert_polynomial_coefficients_type_apparatus
                <chebyshev_polynomial<origin_coefficients_type, kind>,
                        destination_coefficients_type>
        {
            using type = chebyshev_polynomial<destination_coefficients_type, kind>;
        };

        template<real_t alpha, real_t beta,
                 class origin_coefficients_type, class destination_coefficients_type>
        struct convert_polynomial_coefficients_type_apparatus
                <jacobi_polynomial<alpha, beta, origin_coefficients_type>,
                        destination_coefficients_type>
        {
            using type = jacobi_polynomial<alpha, beta, destination_coefficients_type>;
        };
    }

    // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //

    namespace Traits
    {
        template<polynomial_c polynomial_type>
        using univariate_polynomial_type =
                Apparatus::univariate_polynomial_type_apparatus<polynomial_type>::type;

        template<polynomial_c polynomial_type>
        using polynomial_proxy_type =
                Apparatus::polynomial_proxy_type_apparatus<polynomial_type>::type;

        template<polynomial_c polynomial_type>
        using polynomial_remove_proxy_type =
                Apparatus::polynomial_remove_proxy_type_apparatus<polynomial_type>::type;

        template<polynomial_c polynomial_type, class destination_coefficients_type>
        using convert_polynomial_coefficients_type =
                Apparatus::convert_polynomial_coefficients_type_apparatus
                        <polynomial_type, destination_coefficients_type>::type;
    }
}
