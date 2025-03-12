#pragma once

#include "transform.hpp"

// **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //

/* This class represent an abstract polynomial p expressed as
 *      sum_{j = 0}^n (c_j * phi_j(X)).
 *
 * c_j represent the coefficients over which the polynomial is taken.
 * Their type is provided by typedef "coefficients_type" that must be public in the
 * class provided as template argument - in fact this is the only requirement of the
 * template argument class.
 * It is assumed that "coefficients_type" represents a commutative ring, with the ring
 * operations been given by the typical operators +, -, * (and / if it is in fact a field).
 * In addition, it is required that the elements 0, 1, and -1 of the ring are constructible
 * as static_cast<coefficeints_type>(0) and similarly for 1, and -1.
 * Finally, they must be comparable to 0 using the syntax c == 0 and c != 0.
 *
 * phi_j represent the polynomial basis. It has two requirements:
 *      phi_0 = 1
 *      phi_n is of exact order n
 *
 * To enable the evaluation of the polynomial by substituting X, one must overload the
 * function evaluate (see below). Similarly, to enable multiplication of two polynomials
 * one must overload the operator* function. Note that the exponentiation (operator^)
 * uses the operator* but it is not required to implement it. Also note that if operator* is
 * provided, composition of polynomials p and q of the same type can be computed as p(q).
 *
 * The polynomial p is represented simply with a "std::vector<coefficients_type>" which
 * holds the values of c_j in the natural ordering. Access to this vector is provided by the
 * coefficients() method. The values of c_j can also be accessed and modified by operator[].
 * The vector is resized automatically by operator[]. It can also be resized manually by the
 * methods "set_order", "set_extents" and "clear". The "extents" method returns the size of
 * the coefficients vector and "order()" returns "extents() - 1". (Note the overflow in case
 * of having an empty vector). Finally, note that "order()" or "extents() - 1" may not
 * return the exact order of the polynomial since there may be trailing zeros in the
 * coefficients vector. The methods "set_order()" or "set_extents()" (which are equivalent)
 * remove these trailing zeros. */

namespace Polynomial
{
    namespace Apparatus
    {
        template<class polynomial_traits, bool is_proxy, class variate_type>
        void evaluate(const polynomial_base<polynomial_traits, is_proxy> &, const variate_type &)
        {
            // This function must return arg0(arg1)
            static_assert(false, "Overload this function to provide an implementation");
        }

        template<class polynomial_traits, bool is_proxy_arg0, bool is_proxy_arg1>
        void operator*
                (const polynomial_base<polynomial_traits, is_proxy_arg0> &,
                 const polynomial_base<polynomial_traits, is_proxy_arg1> &)
        {
            // This function must return arg0 * arg1
            static_assert(false, "Overload this function to provide an implementation");
        }

        // Overload this function to optimise
        template<polynomial_c polynomial_type, class variate_type>
        auto evaluate_basis_functions(const variate_type &x, index_t n, tag<polynomial_type>)
        -> std::vector<decltype(polynomial_type::basis_polynomial(n)(x))>
        {
            using coefficients_type = polynomial_type::coefficients_type;
            using ret_type = decltype(polynomial_type::basis_polynomial(n)(x));
            std::vector<ret_type> ret;

            ret.reserve(n);

            polynomial_type p;
            p.set_extents(n);

            for (index_t it = 0; it != n; ++it)
            {
                p[it] = static_cast<coefficients_type>(1);
                ret.emplace_back(p(x));
                p[it] = static_cast<coefficients_type>(0);
            }

            return ret;
        }

        template<class polynomial_traits_, bool is_proxy_>
        class polynomial_base
        {
        public:
            static constexpr bool polynomial_tag = true;

            using polynomial_traits = polynomial_traits_;
            static constexpr bool is_proxy = is_proxy_;

            using coefficients_type = polynomial_traits::coefficients_type;
            static constexpr index_t n_variates = 1;

        private:
            using coefficients_container_type = std::conditional_t
                    <is_proxy, std::span<const coefficients_type>, std::vector<coefficients_type>>;

            // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //
        private:
            coefficients_container_type coefficients_{};

            // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //
        public:
            polynomial_base()requires(!is_proxy) = default;
            ~polynomial_base() = default;

            polynomial_base(const polynomial_base &) = default;
            polynomial_base &operator=(const polynomial_base &) = default;

            polynomial_base(polynomial_base &&) noexcept = default;
            polynomial_base &operator=(polynomial_base &&) noexcept = default;

            polynomial_base
                    (const polynomial_base<polynomial_traits, true> &)requires(!is_proxy);
            polynomial_base &operator=
                    (const polynomial_base<polynomial_traits, true> &)requires(!is_proxy);

            template<class type>
            explicit polynomial_base(const type &)
            requires(!is_proxy and std::convertible_to<type, coefficients_type> and
                     !Apparatus::implemented_conversion_c<type, polynomial_base>);

            template<class type>
            explicit polynomial_base(const type &arg)
            requires(!is_proxy and Apparatus::implemented_conversion_c<type, polynomial_base>);

            explicit polynomial_base(coefficients_container_type);

            // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //
        public:
            [[nodiscard]] coefficients_type &operator[](index_t) requires(!is_proxy);
            [[nodiscard]] const coefficients_type &operator[](index_t) const;

            [[nodiscard]] std::span<coefficients_type> coefficients() &requires(!is_proxy);
            [[nodiscard]] std::span<const coefficients_type> coefficients() const &;
            [[nodiscard]] std::vector<coefficients_type> &&coefficients() &&requires(!is_proxy);

            [[nodiscard]] bool empty() const;
            polynomial_base &clear()requires(!is_proxy);

            [[nodiscard]] index_t order() const;
            polynomial_base &set_order()requires(!is_proxy);
            polynomial_base &set_order(index_t)requires(!is_proxy);

            [[nodiscard]] index_t extents() const;
            polynomial_base &set_extents()requires(!is_proxy);
            polynomial_base &set_extents(index_t)requires(!is_proxy);

            // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //
        public:
            template<class variate_type>
            auto operator()(const variate_type &x) const -> decltype(evaluate(*this, x));

            template<bool is_proxy_arg0>
            polynomial_base &operator+=
                    (const polynomial_base<polynomial_traits, is_proxy_arg0> &)requires(!is_proxy);

            template<bool is_proxy_arg0>
            polynomial_base &operator-=
                    (const polynomial_base<polynomial_traits, is_proxy_arg0> &)requires(!is_proxy);

            template<bool is_proxy_arg0>
            polynomial_base &operator*=
                    (const polynomial_base<polynomial_traits, is_proxy_arg0> &)requires(!is_proxy);

            template<class type>
            polynomial_base &operator+=(const type &)requires(!is_proxy);

            template<class type>
            polynomial_base &operator-=(const type &)requires(!is_proxy);

            template<class type>
            polynomial_base &operator*=(const type &)requires(!is_proxy);

            template<class type>
            polynomial_base &operator/=(const type &)requires(!is_proxy);

            polynomial_base &operator^=(index_t)requires(!is_proxy);

            template<bool is_proxy_arg0>
            bool operator==(const polynomial_base<polynomial_traits, is_proxy_arg0> &) const;
            template<bool is_proxy_arg0>
            bool operator!=(const polynomial_base<polynomial_traits, is_proxy_arg0> &) const;

            // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //
        public:
            static polynomial_base basis_polynomial(index_t)requires(!is_proxy);
        };

        // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //
        // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //

        template<class polynomial_traits, bool is_proxy>
        polynomial_base<polynomial_traits, is_proxy>::polynomial_base
                (const polynomial_base<polynomial_traits, true> &arg0)requires(!is_proxy) :
                coefficients_{std::cbegin(arg0.coefficients()), std::cend(arg0.coefficients())} {}

        template<class polynomial_traits, bool is_proxy>
        auto polynomial_base<polynomial_traits, is_proxy>::operator=
                (const polynomial_base<polynomial_traits, true> &arg0)
        -> polynomial_base &requires(!is_proxy)
        {
            coefficients_ = coefficients_container_type{arg0.coefficients()};
        }

        template<class polynomial_traits, bool is_proxy>
        template<class type>
        polynomial_base<polynomial_traits, is_proxy>::polynomial_base
                (const type &arg0)
        requires(!is_proxy and
                 std::convertible_to<type, coefficients_type> and
                 !Apparatus::implemented_conversion_c<type, polynomial_base>) :
                coefficients_{static_cast<coefficients_type>(arg0)} {}

        template<class polynomial_traits, bool is_proxy>
        template<class type>
        polynomial_base<polynomial_traits, is_proxy>::polynomial_base
                (const type &arg)
        requires(!is_proxy and Apparatus::implemented_conversion_c<type, polynomial_base>)
        {
            *this = conversion<type, polynomial_base>::convert(arg);

            return;
        }


        template<class polynomial_traits, bool is_proxy>
        polynomial_base<polynomial_traits, is_proxy>::polynomial_base
                (coefficients_container_type arg0) : coefficients_{std::move(arg0)} {}

        // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //

        template<class polynomial_traits, bool is_proxy>
        auto polynomial_base<polynomial_traits, is_proxy>::operator[]
                (index_t i) -> coefficients_type &requires(!is_proxy)
        {
            if (i >= extents())
                coefficients_.resize(i + 1, static_cast<coefficients_type>(0));

            return coefficients_[i];
        }

        template<class polynomial_traits, bool is_proxy>
        auto polynomial_base<polynomial_traits, is_proxy>::operator[](index_t i) const
        -> const coefficients_type &
        {
            ASSERT_ASSUME(i < extents());
            return coefficients_[i];
        }

        template<class polynomial_traits, bool is_proxy>
        auto polynomial_base<polynomial_traits, is_proxy>::coefficients() &
        -> std::span<coefficients_type>requires(!is_proxy)
        {
            return {coefficients_};
        }

        template<class polynomial_traits, bool is_proxy>
        auto polynomial_base<polynomial_traits, is_proxy>::coefficients() const &
        -> std::span<const coefficients_type>
        {
            return {coefficients_};
        }

        template<class polynomial_traits, bool is_proxy>
        auto polynomial_base<polynomial_traits, is_proxy>::coefficients() &&
        -> std::vector<coefficients_type> &&requires(!is_proxy)
        {
            return std::move(coefficients_);
        }

        template<class polynomial_traits, bool is_proxy>
        bool polynomial_base<polynomial_traits, is_proxy>::empty() const
        {
            return coefficients_.empty();
        }

        template<class polynomial_traits, bool is_proxy>
        auto polynomial_base<polynomial_traits, is_proxy>::clear()
        -> polynomial_base &requires(!is_proxy)
        {
            coefficients_.clear();
            return *this;
        }

        template<class polynomial_traits, bool is_proxy>
        index_t polynomial_base<polynomial_traits, is_proxy>::order() const
        {
            return extents() - 1;
        }

        template<class polynomial_traits, bool is_proxy>
        auto
        polynomial_base<polynomial_traits, is_proxy>::set_order()
        -> polynomial_base &requires(!is_proxy)
        {
            set_extents();
            return *this;
        }

        template<class polynomial_traits, bool is_proxy>
        auto polynomial_base<polynomial_traits, is_proxy>::set_order(index_t arg0)
        -> polynomial_base &requires(!is_proxy)
        {
            set_extents(arg0 + 1);
            return *this;
        }

        template<class polynomial_traits, bool is_proxy>
        index_t polynomial_base<polynomial_traits, is_proxy>::extents() const
        {
            return std::size(coefficients_);
        }

        template<class polynomial_traits, bool is_proxy>
        auto polynomial_base<polynomial_traits, is_proxy>::set_extents()
        -> polynomial_base &requires(!is_proxy)
        {
            index_t i = extents() - 1;

            for (; i != negative_1 and coefficients_[i] == static_cast<coefficients_type>(0); --i);

            coefficients_.resize(++i);
            return *this;
        }

        template<class polynomial_traits, bool is_proxy>
        auto polynomial_base<polynomial_traits, is_proxy>::set_extents(index_t arg0)
        -> polynomial_base &requires(!is_proxy)
        {
            coefficients_.resize(arg0, static_cast<coefficients_type>(0));
            return *this;
        }

        // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //

        template<class polynomial_traits, bool is_proxy>
        template<class variate_type>
        auto polynomial_base<polynomial_traits, is_proxy>::operator()(const variate_type &x) const
        -> decltype(evaluate(*this, x))
        {
            return evaluate(*this, x);
        }

        template<class polynomial_traits, bool is_proxy>
        template<bool is_proxy_arg0>
        auto polynomial_base<polynomial_traits, is_proxy>::operator+=
                (const polynomial_base<polynomial_traits, is_proxy_arg0> &arg0)
        -> polynomial_base &requires(!is_proxy)
        {
            if (arg0.empty())
                return *this;

            if (this->empty())
                return *this = arg0;

            if (arg0.extents() > this->extents())
                this->set_extents(arg0.extents());

            auto it_this = std::begin(coefficients_);
            auto it_arg0 = std::cbegin(arg0.coefficients());

            while (it_arg0 != std::cend(arg0.coefficients()))
                *it_this++ += *it_arg0++;

            return *this;
        }

        template<class polynomial_traits, bool is_proxy>
        template<bool is_proxy_arg0>
        auto polynomial_base<polynomial_traits, is_proxy>::operator-=
                (const polynomial_base<polynomial_traits, is_proxy_arg0> &arg0)
        -> polynomial_base &requires(!is_proxy)
        {
            if (arg0.empty())
                return *this;

            if (this->empty())
                return *this = -arg0;

            if (arg0.extents() > this->extents())
                this->set_extents(arg0.extents());

            auto it_this = std::begin(coefficients_);
            auto it_arg0 = std::cbegin(arg0.coefficients());

            while (it_arg0 != std::cend(arg0.coefficients()))
                *it_this++ -= *it_arg0++;

            return *this;
        }

        template<class polynomial_traits, bool is_proxy>
        template<bool is_proxy_arg0>
        auto polynomial_base<polynomial_traits, is_proxy>::operator*=
                (const polynomial_base<polynomial_traits, is_proxy_arg0> &arg0)
        -> polynomial_base &requires(!is_proxy)
        {
            return *this = *this * arg0;
        }

        template<class polynomial_traits, bool is_proxy>
        template<class type>
        auto polynomial_base<polynomial_traits, is_proxy>::operator+=
                (const type &arg0) -> polynomial_base &requires(!is_proxy)
        {
            (*this)[0] += arg0;

            return *this;
        }

        template<class polynomial_traits, bool is_proxy>
        template<class type>
        auto polynomial_base<polynomial_traits, is_proxy>::operator-=
                (const type &arg0) -> polynomial_base &requires(!is_proxy)
        {
            (*this)[0] -= arg0;

            return *this;
        }

        template<class polynomial_traits, bool is_proxy>
        template<class type>
        auto polynomial_base<polynomial_traits, is_proxy>::operator*=
                (const type &arg0) -> polynomial_base &requires(!is_proxy)
        {
            std::ranges::for_each
                    (coefficients_, [&arg0](coefficients_type &arg_bis) { arg_bis *= arg0; });

            return *this;
        }

        template<class polynomial_traits, bool is_proxy>
        template<class type>
        auto polynomial_base<polynomial_traits, is_proxy>::operator/=
                (const type &arg0) -> polynomial_base &requires(!is_proxy)
        {
            std::ranges::for_each
                    (coefficients_, [&arg0](coefficients_type &arg_bis) { arg_bis /= arg0; });

            return *this;
        }

        template<class polynomial_traits, bool is_proxy>
        auto polynomial_base<polynomial_traits, is_proxy>::operator^=
                (index_t arg0) -> polynomial_base &requires(!is_proxy)
        {
            return *this = *this ^ arg0;
        }

        template<class polynomial_traits, bool is_proxy>
        template<bool is_proxy_arg0>
        bool polynomial_base<polynomial_traits, is_proxy>::operator==
                (const polynomial_base<polynomial_traits, is_proxy_arg0> &arg0) const
        {
            if (this->extents() != arg0.extents())
                return false;

            return std::ranges::equal(coefficients_, arg0.coefficients());
        }

        template<class polynomial_traits, bool is_proxy>
        template<bool is_proxy_arg0>
        bool polynomial_base<polynomial_traits, is_proxy>::operator!=
                (const polynomial_base<polynomial_traits, is_proxy_arg0> &arg0) const
        {
            return !(*this == arg0);
        }

        // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //

        template<class polynomial_traits, bool is_proxy>
        auto polynomial_base<polynomial_traits, is_proxy>::basis_polynomial(index_t arg0)
        -> polynomial_base requires(!is_proxy)
        {
            polynomial_base ret;
            ret[arg0] = static_cast<coefficients_type>(1);

            return ret;
        }

        // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //
        // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //

        template<class polynomial_traits>
        polynomial_base<polynomial_traits> operator-(polynomial_base<polynomial_traits> arg0)
        {
            using coefficients_type = typename polynomial_traits::coefficients_type;

            std::ranges::for_each
                    (arg0.coefficients(),
                     [](coefficients_type &arg_bis) { arg_bis = -std::move(arg_bis); });

            return arg0;
        }

        template<class polynomial_traits, bool is_proxy_arg1>
        polynomial_base<polynomial_traits> operator+
                (polynomial_base<polynomial_traits> arg0,
                 const polynomial_base<polynomial_traits, is_proxy_arg1> &arg1)
        {
            return arg0 += arg1;
        }

        template<class polynomial_traits, bool is_proxy_arg1>
        polynomial_base<polynomial_traits> operator-
                (polynomial_base<polynomial_traits> arg0,
                 const polynomial_base<polynomial_traits, is_proxy_arg1> &arg1)
        {
            return arg0 -= arg1;
        }

        template<class polynomial_traits>
        polynomial_base<polynomial_traits> operator+
                (polynomial_base<polynomial_traits> arg0,
                 const typename polynomial_traits::coefficients_type &arg1)
        {
            return arg0 += arg1;
        }

        template<class polynomial_traits>
        polynomial_base<polynomial_traits> operator-
                (polynomial_base<polynomial_traits> arg0,
                 const typename polynomial_traits::coefficients_type &arg1)
        {
            return arg0 -= arg1;
        }

        template<class polynomial_traits>
        polynomial_base<polynomial_traits> operator*
                (polynomial_base<polynomial_traits> arg0,
                 const typename polynomial_traits::coefficients_type &arg1)
        {
            return arg0 *= arg1;
        }

        template<class polynomial_traits>
        polynomial_base<polynomial_traits> operator/
                (polynomial_base<polynomial_traits> arg0,
                 const typename polynomial_traits::coefficients_type &arg1)
        {
            return arg0 /= arg1;
        }

        template<class polynomial_traits>
        polynomial_base<polynomial_traits> operator+
                (const typename polynomial_traits::coefficients_type &arg0,
                 polynomial_base<polynomial_traits> arg1)
        {
            return arg1 += arg0;
        }

        template<class polynomial_traits>
        polynomial_base<polynomial_traits> operator-
                (const typename polynomial_traits::coefficients_type &arg0,
                 polynomial_base<polynomial_traits> arg1)
        {
            std::ranges::for_each
                (arg1.coefficients(),
                     [&arg0](auto &arg_bis) { arg_bis = arg0 - std::move(arg_bis); });

            return arg1;
        }

        template<class polynomial_traits>
        polynomial_base<polynomial_traits> operator*
                (const typename polynomial_traits::coefficients_type &arg0,
                 polynomial_base<polynomial_traits> arg1)
        {
            return arg1 *= arg0;
        }

        template<class polynomial_traits, bool is_proxy_arg0>
        polynomial_base<polynomial_traits> operator^
                (const polynomial_base<polynomial_traits, is_proxy_arg0> &arg0, index_t arg1)
        {
            if (arg1 == 0)
                return static_cast<polynomial_base<polynomial_traits>>(1);

            if (arg1 == 1)
                return arg0;

            polynomial_base<polynomial_traits> ret = arg0 ^ (arg1 / 2);
            ret *= ret;

            if (arg1 % 2)
                ret *= arg0;

            return ret;
        }

        // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //

        template<class polynomial_traits>
        polynomial_base<polynomial_traits> operator-
                (polynomial_base<polynomial_traits, true> arg0)
        {
            return -static_cast<polynomial_base<polynomial_traits>>(arg0);
        }

        template<class polynomial_traits, bool is_proxy_arg1>
        polynomial_base<polynomial_traits> operator+
                (polynomial_base<polynomial_traits, true> arg0,
                 const polynomial_base<polynomial_traits, is_proxy_arg1> &arg1)
        {
            return static_cast<polynomial_base<polynomial_traits>>(arg0) + arg1;
        }

        template<class polynomial_traits, bool is_proxy_arg1>
        polynomial_base<polynomial_traits> operator-
                (polynomial_base<polynomial_traits, true> arg0,
                 const polynomial_base<polynomial_traits, is_proxy_arg1> &arg1)
        {
            return static_cast<polynomial_base<polynomial_traits>>(arg0) - arg1;
        }

        template<class polynomial_traits>
        polynomial_base<polynomial_traits> operator+
                (polynomial_base<polynomial_traits, true> arg0,
                 const typename polynomial_traits::coefficients_type &arg1)
        {
            return static_cast<polynomial_base<polynomial_traits>>(arg0) + arg1;
        }

        template<class polynomial_traits>
        polynomial_base<polynomial_traits> operator-
                (polynomial_base<polynomial_traits, true> arg0,
                 const typename polynomial_traits::coefficients_type &arg1)
        {
            return static_cast<polynomial_base<polynomial_traits>>(arg0) - arg1;
        }

        template<class polynomial_traits>
        polynomial_base<polynomial_traits> operator*
                (polynomial_base<polynomial_traits, true> arg0,
                 const typename polynomial_traits::coefficients_type &arg1)
        {
            return static_cast<polynomial_base<polynomial_traits>>(arg0) * arg1;
        }

        template<class polynomial_traits>
        polynomial_base<polynomial_traits> operator/
                (polynomial_base<polynomial_traits, true> arg0,
                 const typename polynomial_traits::coefficients_type &arg1)
        {
            return static_cast<polynomial_base<polynomial_traits>>(arg0) / arg1;
        }

        template<class polynomial_traits>
        polynomial_base<polynomial_traits> operator+
                (const typename polynomial_traits::coefficients_type &arg0,
                 polynomial_base<polynomial_traits, true> arg1)
        {
            return arg0 + static_cast<polynomial_base<polynomial_traits>>(arg1);
        }

        template<class polynomial_traits>
        polynomial_base<polynomial_traits> operator-
                (const typename polynomial_traits::coefficients_type &arg0,
                 polynomial_base<polynomial_traits, true> arg1)
        {
            return arg0 - static_cast<polynomial_base<polynomial_traits>>(arg1);
        }

        template<class polynomial_traits>
        polynomial_base<polynomial_traits> operator*
                (const typename polynomial_traits::coefficients_type &arg0,
                 polynomial_base<polynomial_traits, true> arg1)
        {
            return arg0 * static_cast<polynomial_base<polynomial_traits>>(arg1);
        }
    }
}
