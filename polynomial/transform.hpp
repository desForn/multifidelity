#pragma once

#include "traits.hpp"

#include "invocable/invocable_function.hpp"
#include "invocable/invocable_function_sequential_vector.hpp"
#include "invocable/invocable_function_hash_table.hpp"

// **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //

/*
 * Let p be a polynomial expressed in "origin_polynomial_type" with coefficients c_o.
 * The "conversion" class allows to express p as a "destination_polynomial_type" with
 * coefficients c_d.
 *
 * The linear relation c_d = K c_o holds with an upper triangular matrix K.
 * K is represented as a vector v with numbering as follows
 *
 *     | v[0] v[1] v[3] v[6] ... |
 *     |  0   v[2] v[4] v[7] ... |                             | v[i + j * (j + 1) / 2] if i <= j
 * K = |  0    0   v[5] v[8] ... | ; or, in general, K[i, j] = |
 *     |  0    0    0   v[9] ... |                             | 0                      if i > j
 *     | ...  ...  ...  ...  ... |
 *
 * For conversions between canonical or Jacobi polynomials The vector v can be computed inductively.
 *
 * The Jacobi-canonical transform is computed using the recurrence relation
 * https://en.wikipedia.org/wiki/Jacobi_polynomials#Recurrence_relations.
 *
 * The canonical-Jacobi transform is computed by inverting K (using an inverse_conversion)
 *
 * See https://doi.org/10.1090/mcom/3377 for the Jacobi-Jacobi transform algorithm.
 *
 * Additional algorithms relying on other techniques rather than matrix multiplication can be
 * implemented via specializing the "conversion" class.
 */

namespace Polynomial
{
    template<class origin_polynomial_type, class destination_polynomial_type>
    class conversion
    {
        static_assert(false, "Unimplemented");
    };

    // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //
    // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //

    namespace Apparatus
    {
        template<class origin_polynomial_type, class destination_polynomial_type>
        class conversion_coefficients
        {
        public:
            static constexpr bool implemented = false;

            static void evaluate(index_t)
            {
                static_assert(false, "Specialise this class to provide an implementation");
            }
        };

        template<class origin_polynomial_type, class destination_polynomial_type>
        concept implemented_conversion_coefficients_c =
        requires()
        {
            requires conversion_coefficients<origin_polynomial_type, destination_polynomial_type>::
            implemented;
        };

        template<class origin_polynomial_type, class destination_polynomial_type>
        struct implemented_conversion
        {
            static constexpr bool value = false;
        };

        template<univariate_polynomial_c origin_polynomial_type,
                univariate_polynomial_c destination_polynomial_type>
        requires(implemented_conversion_coefficients_c
                         <Traits::polynomial_remove_proxy_type<origin_polynomial_type>,
                          destination_polynomial_type> or
                 implemented_conversion_coefficients_c
                         <destination_polynomial_type,
                          Traits::polynomial_remove_proxy_type<origin_polynomial_type>>)
        struct implemented_conversion<origin_polynomial_type, destination_polynomial_type>
        {
            static constexpr bool value = true;
        };

        template<>
        struct implemented_conversion<chebyshev_polynomial_i_kind, chebyshev_polynomial_point_value>
        {
            static constexpr bool value = true;
        };

        template<>
        struct implemented_conversion<chebyshev_polynomial_point_value, chebyshev_polynomial_i_kind>
        {
            static constexpr bool value = true;
        };

        template<univariate_polynomial_c origin_univariate_polynomial_type,
                univariate_polynomial_c destination_univariate_polynomial_type,
                index_t n_variates>
        requires(implemented_conversion<origin_univariate_polynomial_type,
                                        destination_univariate_polynomial_type>::value)
        struct implemented_conversion
                <multivariate_polynomial<origin_univariate_polynomial_type, n_variates>,
                 multivariate_polynomial<destination_univariate_polynomial_type, n_variates>>
        {
            static constexpr bool value = true;
        };

        template<class origin_polynomial_type, class destination_polynomial_type>
        concept implemented_conversion_c =
        requires()
        {
            requires implemented_conversion
                    <origin_polynomial_type, destination_polynomial_type>::value;
        };

        // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //
        // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //

        template<class, class>
        struct conversion_type
        {
            using type = void;
        };

        template<univariate_polynomial_c origin_polynomial_type_,
                univariate_polynomial_c destination_polynomial_type_>
        class forward_conversion;

        template<univariate_polynomial_c origin_polynomial_type_,
                univariate_polynomial_c destination_polynomial_type_>
        class backward_conversion;

        template<multivariate_polynomial_c origin_polynomial_type_,
                 multivariate_polynomial_c destination_polynomial_type_>
        class multivariate_forward_conversion;

        // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //
        // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //

        template<univariate_polynomial_c origin_polynomial_type_,
                univariate_polynomial_c destination_polynomial_type_>
        class forward_conversion
        {
            // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //

            template<univariate_polynomial_c origin_polynomial_type_bis,
                    univariate_polynomial_c destination_polynomial_type_bis>
            friend
            class backward_conversion;

            template<class origin_polynomial_type_bis, class destination_polynomial_type_bis>
            friend
            class conversion_coefficients;

            template<multivariate_polynomial_c origin_polynomial_type_bis,
                     multivariate_polynomial_c destination_polynomial_type_bis>
            friend
            class multivariate_forward_conversion;

            // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //
        public:
            using origin_polynomial_type =
                    Traits::polynomial_remove_proxy_type<origin_polynomial_type_>;
            using origin_polynomial_type_proxy =
                    Traits::polynomial_proxy_type<origin_polynomial_type>;
            using destination_polynomial_type = destination_polynomial_type_;

            // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //
        private:
            static inline Invocable::invocable_function_sequential_vector invocable_
                    {conversion_coefficients
                     < origin_polynomial_type, destination_polynomial_type > ::evaluate};

            // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //
        public:
            // It has to be used statically
            forward_conversion() = delete;
            ~forward_conversion() = delete;

            // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //
        public:
            template<class type>
            static destination_polynomial_type convert(const type &)
            requires(std::is_same_v<type, origin_polynomial_type> or
                     std::is_same_v<type, origin_polynomial_type_proxy>);
            static void clear_cache();
        };

        /* The matrix K is computed using a recursive algorithm like in
         * [https://doi.org/10.1090/mcom/3377, Theorem 2.4]. The notation there is followed except
         * the fact that there is a replacement of j + 1 by j.
         *
         * Class forward_conversion_coefficients provides 4 methods to compute epsilon_{1...4}.
         * It also must provide a method named base_case that returns the values of the base case.
         * This class must be specialized for every conversion formula.
         * If any of the methods should return always zero return a void type. */

        template<univariate_polynomial_c origin_polynomial_type_,
                univariate_polynomial_c destination_polynomial_type_>
        struct forward_conversion_coefficients;

        template<univariate_polynomial_c origin_polynomial_type_,
                univariate_polynomial_c destination_polynomial_type_>
        class backward_conversion
        {
            template<multivariate_polynomial_c origin_polynomial_type_bis,
                     multivariate_polynomial_c destination_polynomial_type_bis>
            friend
            class multivariate_forward_conversion;

            // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //
        public:
            using origin_polynomial_type =
                    Traits::polynomial_remove_proxy_type<origin_polynomial_type_>;
            using origin_polynomial_type_proxy =
                    Traits::polynomial_proxy_type<origin_polynomial_type>;
            using destination_polynomial_type = destination_polynomial_type_;

            using forward_conversion_type = forward_conversion
                    <destination_polynomial_type, origin_polynomial_type>;

            using coefficients_type = destination_polynomial_type::coefficients_type;

            // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //
        private:
            static coefficients_type conversion_coefficients(index_t);

            // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //
        private:
            static inline Invocable::invocable_function_sequential_vector
                    invocable_{conversion_coefficients};

            static inline std::vector<coefficients_type> cache{}; // See "conversion_coefficients"

            // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //
        public:
            // It has to be used statically
            backward_conversion() = delete;
            ~backward_conversion() = delete;

            // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //
        public:
            template<class type>
            static destination_polynomial_type convert(const type &)
            requires(std::is_same_v<type, origin_polynomial_type> or
                     std::is_same_v<type, origin_polynomial_type_proxy>);
            static void clear_cache();
            static void clear_forward_transform_cache();
        };

        // The mutivariate case is implemented with a tensor product approach

        template<multivariate_polynomial_c origin_polynomial_type_,
                 multivariate_polynomial_c destination_polynomial_type_>
        class multivariate_forward_conversion
        {
        public:
            using origin_polynomial_type =
                    Traits::polynomial_remove_proxy_type<origin_polynomial_type_>;
            using origin_polynomial_type_proxy =
                    Traits::polynomial_proxy_type<origin_polynomial_type>;

            using destination_polynomial_type = destination_polynomial_type_;

            using origin_univariate_polynomial_type =
                    Traits::univariate_polynomial_type<origin_polynomial_type>;

            using destination_univariate_polynomial_type =
                    Traits::univariate_polynomial_type<destination_polynomial_type>;

            using conversion_type = conversion_type
                    <origin_univariate_polynomial_type,
                     destination_univariate_polynomial_type>::type;

            static_assert(!std::is_same_v<conversion_type, void>);

            static constexpr index_t n_variates = origin_polynomial_type::n_variates;
            using indices_type = std::array<index_t, n_variates>;

            using coefficients_type = destination_polynomial_type::coefficients_type;

            // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //
        public:
            // It has to be used statically
            multivariate_forward_conversion() = delete;
            ~multivariate_forward_conversion() = delete;

            // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //
        public:
            template<class type>
            static destination_polynomial_type convert(const type &)
            requires(std::is_same_v<type, origin_polynomial_type> or
                     std::is_same_v<type, origin_polynomial_type_proxy>);
        };

        // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //
        // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //

        template<univariate_polynomial_c origin_polynomial_type,
                univariate_polynomial_c destination_polynomial_type>
        template<class type>
        destination_polynomial_type
        forward_conversion<origin_polynomial_type, destination_polynomial_type>::
        convert(const type &arg)
        requires(std::is_same_v<type, origin_polynomial_type> or
                 std::is_same_v<type, origin_polynomial_type_proxy>)
        {
            destination_polynomial_type ret;
            ret.set_extents(arg.extents());

            for (index_t i = 0, c = 0; i != arg.extents(); ++i)
                for (index_t j = 0; j != i + 1; ++j, ++c)
                    ret[j] += arg[i] * invocable_(c);

            return ret;
        }

        template<univariate_polynomial_c origin_polynomial_type,
                univariate_polynomial_c destination_polynomial_type>
        void forward_conversion<origin_polynomial_type, destination_polynomial_type>::clear_cache()
        {
            invocable_.clear();
            return;
        }

        // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //
        // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //

        template<univariate_polynomial_c origin_polynomial_type,
                univariate_polynomial_c destination_polynomial_type>
        auto backward_conversion<origin_polynomial_type, destination_polynomial_type>::
        conversion_coefficients(index_t n) -> coefficients_type
        {
            /* This function computes inverse(K)[i, j] where K is the matrix associated to the
             * forward transform.
             *
             * To invert the triangular matrix, the algorithm must procede in a diferent order e.g.
             * to compute inverse(K)[0, 2], the coefficients [1, 2] and [2, 2] must be computed.
             * so when this function is called to compute [0, j], all coefficients up to [j, j]
             * are computed and stored inside the vector "cache".
             *
             * In general, when this function is called to compute [i, j] first it checks if "cache"
             * holds the values of [0, j] ... [j, j] by checking if j == size(cache) + 1.
             * If it is not the case, "cache" is computed */

            /* if i == 0, then j * j + j - 2 * n == 0, so j == (sqrt(8 * n + 1) - 1) / 2
             * if i != 0, as i <= j, then j == trunc((sqrt(8 * n + 1) - 1) / 2)
             * to avoid rounding errors issues, substitute n by n + 0.5 */

            index_t j = static_cast<index_t>
            ((std::sqrt(8 * static_cast<real_t>(n + 0.5) + 1) - 1) / 2);

            index_t i = n - j * (j + 1) / 2;

            if (j == std::size(cache) + 1)
                return cache[i];

            cache.resize(j + 1);

            for (index_t k = j; k != static_cast<index_t>(-1); --k)
            {
                coefficients_type aux = static_cast<coefficients_type> (k == j);

                for (index_t l = k + 1; l != j + 1; ++l)
                    aux -= cache[l] *
                           forward_conversion_type::invocable_(k + l * (l + 1) / 2);

                cache[k] = aux / forward_conversion_type::invocable_(k + k * (k + 1) / 2);
            }

            return cache[i];
        }

        // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //

        template<univariate_polynomial_c origin_polynomial_type,
                univariate_polynomial_c destination_polynomial_type>
        template<class type>
        destination_polynomial_type backward_conversion
                <origin_polynomial_type, destination_polynomial_type>::convert(const type &arg)
        requires(std::is_same_v<type, origin_polynomial_type> or
                 std::is_same_v<type, origin_polynomial_type_proxy>)
        {
            destination_polynomial_type ret;
            ret.set_extents(arg.extents());

            for (index_t i = 0, c = 0; i != arg.extents(); ++i)
                for (index_t j = 0; j != i + 1; ++j, ++c)
                    ret[j] += arg[i] * invocable_(c);

            return ret;
        }

        template<univariate_polynomial_c origin_polynomial_type,
                univariate_polynomial_c destination_polynomial_type>
        void backward_conversion<origin_polynomial_type, destination_polynomial_type>::
        clear_cache()
        {
            invocable_.clear();
            return;
        }

        template<univariate_polynomial_c origin_polynomial_type,
                univariate_polynomial_c destination_polynomial_type>
        void backward_conversion<origin_polynomial_type, destination_polynomial_type>::
        clear_forward_transform_cache()
        {
            forward_conversion_type::clear_cache();

            return;
        }


        // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //
        // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //

        template<multivariate_polynomial_c origin_polynomial_type,
                 multivariate_polynomial_c destination_polynomial_type>
        template<class type>
        destination_polynomial_type multivariate_forward_conversion
                <origin_polynomial_type, destination_polynomial_type>::convert(const type &arg)
        requires(std::is_same_v<type, origin_polynomial_type> or
                 std::is_same_v<type, origin_polynomial_type_proxy>)
        {
            const indices_type &extents = arg.extents();

            destination_polynomial_type ret;
            ret.set_extents(extents);

            for (const indices_type &i: Utility::counter(extents))
                for (const indices_type &j: Utility::counter(Utility::sum_array(i, 1)))
                {
                    coefficients_type c{1};

                    for (index_t v = 0; v != n_variates; ++v)
                        c *= conversion_type::invocable_(j[v] + i[v] * (i[v] + 1) / 2);

                    ret[j] += arg[i] * c;
                }

            return ret;
        }


        // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //
        // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //

        template<univariate_polynomial_c origin_polynomial_type,
                 univariate_polynomial_c destination_polynomial_type>
        requires(implemented_conversion_coefficients_c
                    <origin_polynomial_type, destination_polynomial_type>)
        struct conversion_type<origin_polynomial_type, destination_polynomial_type>
        {
            using type = forward_conversion<origin_polynomial_type, destination_polynomial_type>;
        };

        template<univariate_polynomial_c origin_polynomial_type,
                univariate_polynomial_c destination_polynomial_type>
        requires(!implemented_conversion_coefficients_c
                    <origin_polynomial_type, destination_polynomial_type> and
                 implemented_conversion_coefficients_c
                    <destination_polynomial_type, origin_polynomial_type>)
        struct conversion_type<origin_polynomial_type, destination_polynomial_type>
        {
            using type = backward_conversion<origin_polynomial_type, destination_polynomial_type>;
        };

        template<multivariate_polynomial_c origin_polynomial_type,
                 multivariate_polynomial_c destination_polynomial_type>
        requires(implemented_conversion_c<origin_polynomial_type, destination_polynomial_type>)
        struct conversion_type<origin_polynomial_type, destination_polynomial_type>
        {
            using type = multivariate_forward_conversion
                    <origin_polynomial_type, destination_polynomial_type>;
        };
    }

    // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //
    // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //

    template<class origin_polynomial_type_, class destination_polynomial_type_>
    requires(Apparatus::implemented_conversion_c<origin_polynomial_type_, destination_polynomial_type_>)
    class conversion<origin_polynomial_type_, destination_polynomial_type_>
    {
        // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //
    public:
        using origin_polynomial_type =
                Traits::polynomial_remove_proxy_type<origin_polynomial_type_>;
        using origin_polynomial_proxy_type =
                Traits::polynomial_proxy_type<origin_polynomial_type>;
        using destination_polynomial_type = destination_polynomial_type_;

    private:
        using conversion_type =
                Apparatus::conversion_type<origin_polynomial_type, destination_polynomial_type>::type;

        static_assert(!std::is_same_v<conversion_type, void>);

        // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //
    public:
        // It has to be used statically
        conversion() = delete;
        ~conversion() = delete;

        // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //
    public:
        template<class type>
        static destination_polynomial_type convert(const type &arg)
        requires(std::is_same_v<type, origin_polynomial_type> or
                 std::is_same_v<type, origin_polynomial_proxy_type>)
        {
            return conversion_type::convert(arg);
        }
    };
}
