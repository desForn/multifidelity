#pragma once

#include "polynomial_recurrence_relation.hpp"
#include "polynomial_point_value.hpp"
#include "fwd_fourier_transform.hpp"

// **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //

/* https://en.wikipedia.org/wiki/Chebyshev_polynomials
 * The four kinds of chebyshev polynomials are implemented. They correspond to rescaling of
 * the Jacobi polynomials with (alpha, beta) = (-0.5, -0.5), (0.5, 0.5), (-0.5, 0.5), (0.5, -0.5) */


namespace Polynomial
{
    namespace Apparatus
    {
        template<index_t kind, class coefficients_type_>
        struct chebyshev_recurrence_relation
        {
            static_assert(kind > 0 and kind < 5);

            using coefficients_type = coefficients_type_;

            static void gamma_0(index_t) {}
            static real_t gamma_1(index_t);
            static real_t delta(index_t);

            static std::conditional_t<kind <= 2, void, real_t> base_case_10();
            static real_t base_case_11();
        };

        // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //

        template<index_t kind, class coefficients_type>
        real_t chebyshev_recurrence_relation<kind, coefficients_type>::gamma_1(index_t)
        {
            return 2;
        }

        template<index_t kind, class coefficients_type>
        real_t chebyshev_recurrence_relation<kind, coefficients_type>::delta(index_t)
        {
            return -1;
        }

        template<index_t kind, class coefficients_type>
        std::conditional_t<kind <= 2, void, real_t>
        chebyshev_recurrence_relation<kind, coefficients_type>::base_case_10()
        {
            if constexpr (kind <= 2)
                return;

            if constexpr (kind == 3)
                return -1;

            if constexpr (kind == 4)
                return 1;
        }

        template<index_t kind, class coefficients_type>
        real_t chebyshev_recurrence_relation<kind, coefficients_type>::base_case_11()
        {
            if constexpr (kind == 1)
                return 1;

            return 2;
        }

        // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //
        // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //

        template<class coefficients_type, bool is_proxy_lhs, bool is_proxy_rhs>
        chebyshev_polynomial<coefficients_type, 1> operator*
                (const polynomial_base<polynomial_recurrence_relation_traits
                        <chebyshev_recurrence_relation<1, coefficients_type>>, is_proxy_lhs> &lhs,
                 const polynomial_base<polynomial_recurrence_relation_traits
                         <chebyshev_recurrence_relation<1, coefficients_type>>, is_proxy_rhs> &rhs)
        {
            // https://en.wikipedia.org/wiki/Chebyshev_polynomials#Products_of_Chebyshev_polynomials

            if (lhs.empty() or rhs.empty())
                return {};

            index_t lhs_extents = lhs.extents();
            index_t rhs_extents = rhs.extents();

            index_t ret_extents = lhs_extents + rhs_extents - 1;

            chebyshev_polynomial<coefficients_type, 1> ret;
            ret.set_extents(ret_extents);

            index_t max_order = std::max(lhs_extents, rhs_extents) - 1;
            index_t min_order = std::min(lhs_extents, rhs_extents) - 1;

            for (index_t i = 0; i <= max_order; ++i)
                for (index_t j = 0; j <= std::min(i, min_order); ++j)
                {
                    real_t c = 0;

                    if (i < lhs_extents and j < rhs_extents)
                        c = lhs[i] * rhs[j];

                    if (i != j and j < lhs_extents and i < rhs_extents)
                        c += lhs[j] * rhs[i];

                    if (c == 0)
                        continue;

                    ret[i + j] += 0.5 * c;
                    ret[i - j] += 0.5 * c;
                }

            return ret;
        }

        template<class coefficients_type, bool is_proxy_lhs, bool is_proxy_rhs>
        chebyshev_polynomial<coefficients_type, 2> operator*
                (const polynomial_base<polynomial_recurrence_relation_traits
                        <chebyshev_recurrence_relation<2, coefficients_type>>, is_proxy_lhs> &lhs,
                 const polynomial_base<polynomial_recurrence_relation_traits
                         <chebyshev_recurrence_relation<2, coefficients_type>>, is_proxy_rhs> &rhs)
        {
            // https://en.wikipedia.org/wiki/Chebyshev_polynomials#Products_of_Chebyshev_polynomials

            if (lhs.empty() or rhs.empty())
                return {};

            index_t lhs_extents = lhs.extents();
            index_t rhs_extents = rhs.extents();

            index_t ret_extents = lhs_extents + rhs_extents - 1;

            chebyshev_polynomial<coefficients_type, 2> ret;
            ret.set_extents(ret_extents);

            index_t max_order = std::max(lhs_extents, rhs_extents) - 1;
            index_t min_order = std::min(lhs_extents, rhs_extents) - 1;

            for (index_t i = 0; i <= max_order; ++i)
                for (index_t j = 0; j <= std::min(i, min_order); ++j)
                {
                    real_t c = 0;

                    if (i < lhs_extents and j < rhs_extents)
                        c = lhs[i] * rhs[j];

                    if (i != j and j < lhs_extents and i < rhs_extents)
                        c += lhs[j] * rhs[i];

                    if (c == 0)
                        continue;

                    for (index_t k = 0; k <= j; ++k)
                        ret[i - j + 2 * k] += c;
                }

            return ret;
        }

        template<class coefficients_type, bool is_proxy_lhs, bool is_proxy_rhs>
        chebyshev_polynomial<coefficients_type, 3> operator*
                (const polynomial_base<polynomial_recurrence_relation_traits
                        <chebyshev_recurrence_relation<3, coefficients_type>>, is_proxy_lhs> &lhs,
                 const polynomial_base<polynomial_recurrence_relation_traits
                         <chebyshev_recurrence_relation<3, coefficients_type>>, is_proxy_rhs> &rhs)
        {
            // https://www.jstor.org/stable/26413520

            if (lhs.empty() or rhs.empty())
                return {};

            index_t lhs_extents = lhs.extents();
            index_t rhs_extents = rhs.extents();

            index_t ret_extents = lhs_extents + rhs_extents - 1;

            chebyshev_polynomial<coefficients_type, 3> ret;
            ret.set_extents(ret_extents);

            index_t max_order = std::max(lhs_extents, rhs_extents) - 1;
            index_t min_order = std::min(lhs_extents, rhs_extents) - 1;

            for (index_t i = 0; i <= max_order; ++i)
                for (index_t j = 0; j <= std::min(i, min_order); ++j)
                {
                    real_t c = 0;

                    if (i < lhs_extents and j < rhs_extents)
                        c = lhs[i] * rhs[j];

                    if (i != j and j < lhs_extents and i < rhs_extents)
                        c += lhs[j] * rhs[i];

                    if (c == 0)
                        continue;

                    for (index_t k = 0; k < j; ++k)
                    {
                        ret[i - j + 2 * k] += c;
                        ret[i - j + 2 * k + 1] -= c;
                    }

                    ret[i + j] += c;
                }

            return ret;
        }

        template<class coefficients_type, bool is_proxy_lhs, bool is_proxy_rhs>
        chebyshev_polynomial<coefficients_type, 4> operator*
                (const polynomial_base<polynomial_recurrence_relation_traits
                        <chebyshev_recurrence_relation<4, coefficients_type>>, is_proxy_lhs> &lhs,
                 const polynomial_base<polynomial_recurrence_relation_traits
                         <chebyshev_recurrence_relation<4, coefficients_type>>, is_proxy_rhs> &rhs)
        {
            // https://www.jstor.org/stable/26413520

            if (lhs.empty() or rhs.empty())
                return {};

            index_t lhs_extents = lhs.extents();
            index_t rhs_extents = rhs.extents();

            index_t ret_extents = lhs_extents + rhs_extents - 1;

            chebyshev_polynomial<coefficients_type, 4> ret;
            ret.set_extents(ret_extents);

            index_t max_order = std::max(lhs_extents, rhs_extents) - 1;
            index_t min_order = std::min(lhs_extents, rhs_extents) - 1;

            for (index_t i = 0; i <= max_order; ++i)
                for (index_t j = 0; j <= std::min(i, min_order); ++j)
                {
                    real_t c = 0;

                    if (i < lhs_extents and j < rhs_extents)
                        c = lhs[i] * rhs[j];

                    if (i != j and j < lhs_extents and i < rhs_extents)
                        c += lhs[j] * rhs[i];

                    if (c == 0)
                        continue;

                    for (index_t k = 0; k <= 2 * j; ++k)
                        ret[i - j + k] += c;
                }

            return ret;
        }
    }

    // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //

    template<index_t kind, class coefficients_type, bool is_proxy_lhs, bool is_proxy_rhs>
    requires((std::is_same_v<coefficients_type, real_t> or
              std::is_same_v<coefficients_type, complex_t>))
    coefficients_type inner_product
            (const Apparatus::polynomial_base<Apparatus::polynomial_recurrence_relation_traits
                    <Apparatus::chebyshev_recurrence_relation<kind, coefficients_type>>,
                    is_proxy_lhs> &lhs,
             const Apparatus::polynomial_base<Apparatus::polynomial_recurrence_relation_traits
                     <Apparatus::chebyshev_recurrence_relation<kind, coefficients_type>>,
                     is_proxy_rhs> &rhs)
    {
        if (lhs.empty() or rhs.empty())
            return static_cast<coefficients_type>(0);

        constexpr integer_t c = (kind <= 2) ? 1 : 2;

        coefficients_type ret = 2 * lhs[0] * rhs[0];

        for (index_t i = 0; i != std::min(lhs.extents(), rhs.extents()); ++i)
            ret += c * lhs[i] * rhs[i];

        return ret * M_PI_2;
    }

    template<index_t kind, class coefficients_type, bool is_proxy_rhs>
    requires(std::is_same_v<coefficients_type, real_t> and
             std::is_same_v<coefficients_type, complex_t>)
    real_t norm2(const Apparatus::polynomial_base<Apparatus::polynomial_recurrence_relation_traits
            <Apparatus::chebyshev_recurrence_relation<kind, coefficients_type>>, is_proxy_rhs> &rhs)
    {
        if constexpr (std::is_same_v<coefficients_type, real_t>)
            return inner_product(rhs, rhs);

        else
            return inner_product(rhs, rhs).real();
    }

    template<index_t kind, class coefficients_type, bool is_proxy_rhs>
    requires(std::is_same_v<coefficients_type, real_t> and
             std::is_same_v<coefficients_type, complex_t>)
    real_t norm(const Apparatus::polynomial_base<Apparatus::polynomial_recurrence_relation_traits
            <Apparatus::chebyshev_recurrence_relation<kind, coefficients_type>>, is_proxy_rhs> &rhs)
    {
        return std::sqrt(norm2(rhs));
    }

    template<class coefficients_type, bool is_proxy_rhs>
    coefficients_type integrate
            (const Apparatus::polynomial_base<Apparatus::polynomial_recurrence_relation_traits
                    <Apparatus::chebyshev_recurrence_relation<1, coefficients_type>>, is_proxy_rhs> &rhs)
    {
        coefficients_type ret = 0;

        for (index_t i = 0; i < rhs.extents(); i += 2)
            ret += rhs[i] * 2 / (1 - static_cast<real_t>(i * i));

        return ret;
    }

    template<index_t n_variates>
    real_t integrate
            (const multivariate_polynomial<chebyshev_polynomial_i_kind, n_variates> &rhs)
    {
        std::array<index_t, n_variates> begin{Utility::uniform_array<index_t, n_variates>(0)};
        std::array<integer_t, n_variates> strides
            {Utility::uniform_array<integer_t , n_variates>(2)};

        real_t ret = 0;

        for (const auto &i : Utility::counter{begin, rhs.extents(), strides})
        {
            real_t c = 1;

            for (index_t v = 0; v != n_variates; ++v)
                c *= 2 / (1 - static_cast<real_t>(i[v] * i[v]));

            ret += rhs[i] * c;
        }

        return ret;
    }

    // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //
    // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //

    template<>
    class conversion<chebyshev_polynomial_i_kind, chebyshev_polynomial_point_value>
    {
        // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //
    public:
        using origin_polynomial_type = chebyshev_polynomial_i_kind;
        using destination_polynomial_type = chebyshev_polynomial_point_value;

        // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //
    public:
        // It has to be used statically
        conversion() = delete;
        ~conversion() = delete;

        // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //
    public:
        static destination_polynomial_type convert(const origin_polynomial_type &arg)
        {
            return discrete_cosine_transform(arg);
        }
    };

    // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //

    template<>
    class conversion<chebyshev_polynomial_point_value, chebyshev_polynomial_i_kind>
    {
        // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //
    public:
        using origin_polynomial_type = chebyshev_polynomial_point_value;
        using destination_polynomial_type = chebyshev_polynomial_i_kind;

        // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //
    public:
        // It has to be used statically
        conversion() = delete;
        ~conversion() = delete;

        // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //
    public:
        static destination_polynomial_type convert(const origin_polynomial_type &arg)
        {
            return discrete_cosine_transform(arg);
        }
    };

    // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //

    template<index_t n_variates>
    class conversion<multivariate_polynomial<chebyshev_polynomial_point_value, n_variates>,
                     multivariate_polynomial<chebyshev_polynomial_i_kind, n_variates>>
    {
        // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //
    public:
        using origin_polynomial_type =
                multivariate_polynomial<chebyshev_polynomial_point_value, n_variates>;
        using destination_polynomial_type =
                multivariate_polynomial<chebyshev_polynomial_i_kind, n_variates>;

        // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //
    public:
        // It has to be used statically
        conversion() = delete;
        ~conversion() = delete;

        // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //
    public:
        static destination_polynomial_type convert(const origin_polynomial_type &arg)
        {
            return discrete_cosine_transform(arg);
        }
    };

    // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //

    template<index_t n_variates>
    class conversion<multivariate_polynomial<chebyshev_polynomial_i_kind, n_variates>,
                     multivariate_polynomial<chebyshev_polynomial_point_value, n_variates>>
    {
        // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //
    public:
        using origin_polynomial_type =
                multivariate_polynomial<chebyshev_polynomial_i_kind, n_variates>;
        using destination_polynomial_type =
                multivariate_polynomial<chebyshev_polynomial_point_value, n_variates>;

        // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //
    public:
        // It has to be used statically
        conversion() = delete;
        ~conversion() = delete;

        // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //
    public:
        static destination_polynomial_type convert(const origin_polynomial_type &arg)
        {
            return discrete_cosine_transform(arg);
        }
    };
}

#include "fourier_transform.hpp"
