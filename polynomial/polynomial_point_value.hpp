#pragma once

#include "transform.hpp"

#include "invocable/invocable_function_hash_table.hpp"
#include "invocable/evaluate.hpp"
#include "../sampling/exponential_sampling.hpp"
#include "../sampling/chebyshev_lobatto.hpp"

// **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //

namespace Polynomial
{
    template<Sampling::sampling_c sampling_type_>
    class barycentric_interpolation_weights
    {
        // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //
    public:
        using sampling_type = sampling_type_;
        using domain_field = sampling_type::coordinates_type;

    private:
        class factory_apparatus;

        // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //
    private:
        std::vector<domain_field> w_;

        // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //
    public:
        barycentric_interpolation_weights() = delete;
        ~barycentric_interpolation_weights() = default;

        barycentric_interpolation_weights
                (const barycentric_interpolation_weights &) = delete;
        barycentric_interpolation_weights &operator=
                (const barycentric_interpolation_weights &) = delete;

        barycentric_interpolation_weights
                (barycentric_interpolation_weights &&) noexcept = default;
        barycentric_interpolation_weights &operator=
                (barycentric_interpolation_weights &&) noexcept = default;

    private:
        explicit barycentric_interpolation_weights(index_t);

    public:
        static const barycentric_interpolation_weights &factory(index_t);
        static void clear();

        // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //
    public:
        domain_field operator[](index_t) const;
    };

    // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //

    template<class domain_field>
    std::vector<domain_field> compute_barycentric_interpolation_weights
            (std::span<const domain_field>);

    // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //
    // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //

    template<Sampling::sampling_c sampling_type_>
    class quadrature_weights
    {
        // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //
    public:
        using sampling_type = sampling_type_;
        using domain_field = sampling_type::coordinates_type;

    private:
        class factory_apparatus;

        // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //
    private:
        std::vector<domain_field> w_;

        // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //
    public:
        quadrature_weights() = delete;
        ~quadrature_weights() = default;

        quadrature_weights(const quadrature_weights &) = delete;
        quadrature_weights &operator=(const quadrature_weights &) = delete;

        quadrature_weights(quadrature_weights &&) noexcept = default;
        quadrature_weights &operator=(quadrature_weights &&) noexcept = default;

    private:
        explicit quadrature_weights(index_t);

    public:
        static const quadrature_weights &factory(index_t);
        static void clear();

        // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //
    public:
        domain_field operator[](index_t) const;
    };

    // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //

    template<Sampling::sampling_c sampling_type>
    std::vector<typename sampling_type::coordinates_type> compute_quadrature_weights(index_t)
    {
        static_assert(false, "Unimplemented");
    }

    // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //
    // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //

    template<Sampling::sampling_c sampling_type_, class codomain_field_, class container_type_>
    class polynomial_point_value
    {
    public:
        using sampling_type = sampling_type_;
        using weights_type = barycentric_interpolation_weights<sampling_type>;

        using domain_field = sampling_type::coordinates_type;
        using codomain_field = codomain_field_;
        using container_type = container_type_;
        static constexpr bool is_proxy = 
            not std::same_as<container_type, std::vector<codomain_field>>;

        static constexpr bool polynomial_tag = true;
        static constexpr index_t n_variates = 1;

        // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //
    private:
        index_t sampling_level_{negative_1};
        index_t order_
        {sampling_level_ == negative_1 ?
            negative_1 :
            sampling_type::n_points(sampling_level_) - 1};
        container_type f_{};

        // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //
    public:
        polynomial_point_value()requires(not is_proxy) = default;
        ~polynomial_point_value() = default;

        polynomial_point_value(const polynomial_point_value &) = default;
        polynomial_point_value &operator=(const polynomial_point_value &) = default;

        polynomial_point_value(polynomial_point_value &&) noexcept = default;
        polynomial_point_value &operator=(polynomial_point_value &&) noexcept = default;

        template<class container_type_other>
        polynomial_point_value
                (const polynomial_point_value
                        <sampling_type, codomain_field, container_type_other> &)
        requires(not is_proxy and
                 not std::same_as<container_type_other, std::vector<codomain_field>>);

        template<class container_type_other>
        polynomial_point_value &operator=
                (const polynomial_point_value
                        <sampling_type, codomain_field, container_type_other> &)
        requires(not is_proxy and
                 not std::same_as<container_type_other, std::vector<codomain_field>>);

        explicit polynomial_point_value(const codomain_field &)requires(not is_proxy);

        template<class type>
        explicit polynomial_point_value(const type &)
        requires(not is_proxy and
                 Apparatus::implemented_conversion_c<type, polynomial_point_value>);

        template<class function_type>
        polynomial_point_value(const function_type &, index_t)
        requires(Invocable::invocable_c<function_type, codomain_field(domain_field)> and
                 not is_proxy);

        template<class function_type>
        polynomial_point_value(const function_type &, index_t, index_t)
        requires(Invocable::invocable_c<function_type, codomain_field(domain_field)> and
                 not is_proxy);

        polynomial_point_value(container_type, index_t);
        polynomial_point_value(container_type, index_t, index_t);

        // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //
    public:
        [[nodiscard]] index_t sampling_level() const;
        [[nodiscard]] index_t order() const;
        [[nodiscard]] const sampling_type &sampling() const;
        [[nodiscard]] const weights_type &weights() const;

        polynomial_point_value &shrink()requires(not is_proxy);
        polynomial_point_value &set_order(index_t)requires(not is_proxy);

        [[nodiscard]] index_t extents() const;

        [[nodiscard]] bool empty() const;
        polynomial_point_value &clear()requires(not is_proxy);

        [[nodiscard]] decltype(auto) f(this auto &&self) requires(is_proxy) { return (self.f_); }
        [[nodiscard]] std::span<codomain_field> f(this polynomial_point_value &self)
            requires(not is_proxy)
            { return {std::begin(self.f_), std::end(self.f_)}; }
        [[nodiscard]] std::span<const codomain_field> f(this const polynomial_point_value &self)
            requires(not is_proxy)
            { return {std::cbegin(self.f_), std::cend(self.f_)}; }
        [[nodiscard]] std::vector<codomain_field> &&f(this polynomial_point_value&& self)
            requires(not is_proxy)
            { return std::move(self.f_); }

        [[nodiscard]] codomain_field &f(index_t);
        [[nodiscard]] const codomain_field &f(index_t) const;
        [[nodiscard]] domain_field x(index_t) const;
        [[nodiscard]] domain_field w(index_t) const;

        // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //
    public:
        template<Arithmetic::real_complex_c variate_type>
        auto operator()(const variate_type &) const
        -> std::common_type_t<domain_field, codomain_field, variate_type>;

        template<class container_type_other>
        polynomial_point_value &operator+=
                (const polynomial_point_value
                    <sampling_type, codomain_field, container_type_other> &)requires(not is_proxy);

        template<class container_type_other>
        polynomial_point_value &operator-=
                (const polynomial_point_value
                    <sampling_type, codomain_field, container_type_other> &)requires(not is_proxy);

        template<class container_type_other>
        polynomial_point_value &operator*=
                (const polynomial_point_value
                    <sampling_type, codomain_field, container_type_other> &)requires(not is_proxy);

        polynomial_point_value &operator+=(const codomain_field &)requires(not is_proxy);
        polynomial_point_value &operator-=(const codomain_field &)requires(not is_proxy);
        polynomial_point_value &operator*=(const codomain_field &)requires(not is_proxy);
        polynomial_point_value &operator/=(const codomain_field &)requires(not is_proxy);

        polynomial_point_value &operator^=(index_t)requires(not is_proxy);

    private:
        template<class function_type, class container_type_other>
        void operator_assign
            (const polynomial_point_value<sampling_type, codomain_field, container_type_other> &,
             const function_type &, index_t)requires(not is_proxy);
    };

    template<Sampling::sampling_c sampling_type, class codomain_field, class container_type>
    codomain_field integrate
        (const polynomial_point_value<sampling_type, codomain_field, container_type> &);

    // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //
    // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //

    template<Sampling::sampling_c sampling_type, class codomain_field, class container_type>
    template<class container_type_other>
    polynomial_point_value<sampling_type, codomain_field, container_type>::polynomial_point_value
        (const polynomial_point_value <sampling_type, codomain_field, container_type_other> &arg)
    requires(not is_proxy and not std::same_as<container_type_other, std::vector<codomain_field>>) :
            sampling_level_{arg.sampling_level()},
            order_{arg.order()},
            f_{std::cbegin(arg.f()), std::cend(arg.f())} {}

    template<Sampling::sampling_c sampling_type, class codomain_field, class container_type>
    template<class container_type_other>
    auto polynomial_point_value<sampling_type, codomain_field, container_type>::operator=
        (const polynomial_point_value <sampling_type, codomain_field, container_type_other> &arg)
    -> polynomial_point_value &
    requires(not is_proxy and not std::same_as<container_type_other, std::vector<codomain_field>>)
    {
        sampling_level_ = arg.sampling_level();
        order_ = arg.order();
        f_ = container_type{std::begin(arg.f()), std::cend(arg.f())};

        return *this;
    }

    template<Sampling::sampling_c sampling_type, class codomain_field, class container_type>
    polynomial_point_value<sampling_type, codomain_field, container_type>::polynomial_point_value
            (const codomain_field &arg) requires(not is_proxy):
            sampling_level_{0},
            f_{arg} {}

    template<Sampling::sampling_c sampling_type, class codomain_field, class container_type>
    template<class type>
    polynomial_point_value<sampling_type, codomain_field, container_type>::polynomial_point_value
            (const type &arg)
    requires(not is_proxy and Apparatus::implemented_conversion_c<type, polynomial_point_value>)
    {
        *this = conversion<type, polynomial_point_value>::convert(arg);

        return;
    }

    template<Sampling::sampling_c sampling_type, class codomain_field, class container_type>
    template<class function_type>
    polynomial_point_value<sampling_type, codomain_field, container_type>::polynomial_point_value
            (const function_type &function, index_t sampling_level)
    requires(Invocable::invocable_c<function_type, codomain_field(domain_field)> and not is_proxy) :
            sampling_level_{sampling_level},
            f_{Invocable::evaluate(function, sampling().coordinates())} {}

    template<Sampling::sampling_c sampling_type, class codomain_field, class container_type>
    template<class function_type>
    polynomial_point_value<sampling_type, codomain_field, container_type>::polynomial_point_value
            (const function_type &function, index_t sampling_level, index_t order)
    requires(Invocable::invocable_c<function_type, codomain_field(domain_field)> and not is_proxy) :
            polynomial_point_value{function, sampling_level}
    {
        set_order(order);
        return;
    }

    template<Sampling::sampling_c sampling_type, class codomain_field, class container_type>
    polynomial_point_value<sampling_type, codomain_field, container_type>::polynomial_point_value
            (container_type function, index_t sampling_level):
            sampling_level_{sampling_level},
            f_{std::move(function)}
    {
        ASSERT_ASSUME(std::size(f_) == std::size(sampling().coordinates()));

        return;
    }

    template<Sampling::sampling_c sampling_type, class codomain_field, class container_type>
    polynomial_point_value<sampling_type, codomain_field, container_type>::polynomial_point_value
            (container_type values, index_t sampling_level, index_t order) :
            sampling_level_{sampling_level},
            order_{order},
            f_{std::move(values)}
    {
        ASSERT_ASSUME(std::size(f_) == std::size(sampling().coordinates()));
        ASSERT_ASSUME(std::size(f_) > order_);

        return;
    }

    // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //

    template<Sampling::sampling_c sampling_type, class codomain_field, class container_type>
    index_t polynomial_point_value<sampling_type, codomain_field, container_type>::
        sampling_level() const
    {
        return sampling_level_;
    }

    template<Sampling::sampling_c sampling_type, class codomain_field, class container_type>
    index_t polynomial_point_value<sampling_type, codomain_field, container_type>::order() const
    {
        return order_;
    }

    template<Sampling::sampling_c sampling_type, class codomain_field, class container_type>
    auto polynomial_point_value<sampling_type, codomain_field, container_type>::sampling() const
    -> const sampling_type &
    {
        ASSERT_ASSUME(sampling_level_ != negative_1);
        return sampling_type::factory(sampling_level_);
    }

    template<Sampling::sampling_c sampling_type, class codomain_field, class container_type>
    auto polynomial_point_value<sampling_type, codomain_field, container_type>::weights() const
    -> const weights_type &
    {
        ASSERT_ASSUME(sampling_level_ != negative_1);
        return weights_type::factory(sampling_level_);
    }

    template<Sampling::sampling_c sampling_type, class codomain_field, class container_type>
    auto polynomial_point_value<sampling_type, codomain_field, container_type>::shrink()
    -> polynomial_point_value &requires(not is_proxy)
    {
        if (order_ == negative_1)
        {
            this->clear();
            return *this;
        }

        index_t new_sampling_level = sampling_type::required_level(order_ + 1);

        if (new_sampling_level != sampling_level_)
            *this = polynomial_point_value{*this, sampling_level_};

        return *this;
    }

    template<Sampling::sampling_c sampling_type, class codomain_field, class container_type>
    auto polynomial_point_value<sampling_type, codomain_field, container_type>::set_order
            (index_t new_order) -> polynomial_point_value &requires(not is_proxy)
    {
        // Set manually the order if the exact order (or an upper bound) is known
        ASSERT_ASSUME(new_order == negative_1 or new_order < extents());

        order_ = new_order;

        return *this;
    }

    template<Sampling::sampling_c sampling_type, class codomain_field, class container_type>
    index_t polynomial_point_value<sampling_type, codomain_field, container_type>::extents() const
    {
        return std::size(f_);
    }

    template<Sampling::sampling_c sampling_type, class codomain_field, class container_type>
    bool polynomial_point_value<sampling_type, codomain_field, container_type>::empty() const
    {
        return std::empty(f_);
    }

    template<Sampling::sampling_c sampling_type, class codomain_field, class container_type>
    auto polynomial_point_value<sampling_type, codomain_field, container_type>::clear()
    -> polynomial_point_value &requires(not is_proxy)
    {
        order_ = negative_1;
        sampling_level_ = negative_1;
        f_.clear();

        return *this;
    }

    template<Sampling::sampling_c sampling_type, class codomain_field, class container_type>
    auto polynomial_point_value<sampling_type, codomain_field, container_type>::f(index_t i)
    -> codomain_field &
    {
        return f_[i];
    }

    template<Sampling::sampling_c sampling_type, class codomain_field, class container_type>
    auto polynomial_point_value<sampling_type, codomain_field, container_type>::f(index_t i) const
    -> const codomain_field &
    {
        return f_[i];
    }

    template<Sampling::sampling_c sampling_type, class codomain_field, class container_type>
    auto polynomial_point_value<sampling_type, codomain_field, container_type>::x(index_t i) const
    -> domain_field
    {
        return sampling()[i];
    }

    template<Sampling::sampling_c sampling_type, class codomain_field, class container_type>
    auto polynomial_point_value<sampling_type, codomain_field, container_type>::w(index_t i) const
    -> domain_field
    {
        return weights()[i];
    }

    // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //

    template<Sampling::sampling_c sampling_type, class codomain_field, class container_type>
    template<Arithmetic::real_complex_c variate_type>
    auto polynomial_point_value<sampling_type, codomain_field, container_type>::operator()
            (const variate_type &x) const
    -> std::common_type_t<domain_field, codomain_field, variate_type>
    {
        using ret_type = std::common_type_t<domain_field, codomain_field, variate_type>;

        const sampling_type &sampling_ = sampling();
        const weights_type &weights_ = weights();

        if (order_ == negative_1)
            return static_cast<variate_type>(0);

        ret_type num = static_cast<variate_type>(0);
        ret_type den = static_cast<variate_type>(0);

        auto is_valid = [](const ret_type &arg_bis) -> bool
        {
            if constexpr (std::is_same_v<ret_type, real_t>)
                return (std::isfinite(arg_bis));

            else
                return std::isfinite(arg_bis.real()) and std::isfinite(arg_bis.imag());
        };

        for (index_t i = 0; i != extents(); ++i)
        {
            ret_type a = weights_[i] / (x - sampling_[i]);

            if (!is_valid(a))
                return f_[i];

            num += a * f_[i];
            den += a;
        }

        ret_type ret = num / den;

        return ret;
    }

    template<Sampling::sampling_c sampling_type, class codomain_field, class container_type>
    template<class container_type_other>
    auto polynomial_point_value<sampling_type, codomain_field, container_type>::operator+=
            (const polynomial_point_value
                    <sampling_type, codomain_field, container_type_other> &rhs)
    -> polynomial_point_value &requires(not is_proxy)
    {
        if (rhs.order() == negative_1)
            return *this;

        if (this->order() == negative_1)
            return *this = rhs;

        operator_assign(rhs, std::plus{}, std::max(this->order(), rhs.order()));

        return *this;
    }

    template<Sampling::sampling_c sampling_type, class codomain_field, class container_type>
    template<class container_type_other>
    auto polynomial_point_value<sampling_type, codomain_field, container_type>::operator-=
            (const polynomial_point_value
                    <sampling_type, codomain_field, container_type_other> &rhs)
    -> polynomial_point_value &requires(not is_proxy)
    {
        if (rhs.order() == negative_1)
            return *this;

        if (this->order() == negative_1)
            return *this = -rhs;

        operator_assign(rhs, std::minus{}, std::max(this->order(), rhs.order()));

        return *this;
    }

    template<Sampling::sampling_c sampling_type, class codomain_field, class container_type>
    template<class container_type_other>
    auto polynomial_point_value<sampling_type, codomain_field, container_type>::operator*=
            (const polynomial_point_value
                    <sampling_type, codomain_field, container_type_other> &rhs)
    -> polynomial_point_value &requires(not is_proxy)
    {
        if (this->order() == negative_1)
            return *this;

        if (rhs.order() == negative_1)
            return *this->clear();

        operator_assign(rhs, std::multiplies{}, this->order() + rhs.order());

        return *this;
    }

    template<Sampling::sampling_c sampling_type, class codomain_field, class container_type>
    auto polynomial_point_value<sampling_type, codomain_field, container_type>::operator+=
            (const codomain_field &rhs) -> polynomial_point_value &requires(not is_proxy)
    {
        std::ranges::for_each(f_, [&rhs](codomain_field &lhs) { lhs += rhs; });

        return *this;
    }

    template<Sampling::sampling_c sampling_type, class codomain_field, class container_type>
    auto polynomial_point_value<sampling_type, codomain_field, container_type>::operator-=
            (const codomain_field &rhs) -> polynomial_point_value &requires(not is_proxy)
    {
        std::ranges::for_each(f_, [&rhs](codomain_field &lhs) { lhs -= rhs; });

        return *this;
    }

    template<Sampling::sampling_c sampling_type, class codomain_field, class container_type>
    auto polynomial_point_value<sampling_type, codomain_field, container_type>::operator*=
            (const codomain_field &rhs) -> polynomial_point_value &requires(not is_proxy)
    {
        std::ranges::for_each(f_, [&rhs](codomain_field &lhs) { lhs *= rhs; });

        return *this;
    }

    template<Sampling::sampling_c sampling_type, class codomain_field, class container_type>
    auto polynomial_point_value<sampling_type, codomain_field, container_type>::operator/=
            (const codomain_field &rhs) -> polynomial_point_value &requires(not is_proxy)
    {
        std::ranges::for_each(f_, [&rhs](codomain_field &lhs) { lhs /= rhs; });

        return *this;
    }

    template<Sampling::sampling_c sampling_type, class codomain_field, class container_type>
    auto polynomial_point_value<sampling_type, codomain_field, container_type>::operator^=
            (index_t rhs) -> polynomial_point_value &requires(not is_proxy)
    {
        return *this = *this ^ rhs;
    }

    template<Sampling::sampling_c sampling_type, class codomain_field, class container_type>
    template<class function_type, class container_type_other>
    void polynomial_point_value<sampling_type, codomain_field, container_type>::operator_assign
            (const polynomial_point_value
                    <sampling_type, codomain_field, container_type_other> &rhs,
             const function_type &function, index_t new_order)requires(not is_proxy)
    {
        index_t new_sampling_level = sampling_type::required_level(new_order + 1);

        if (new_sampling_level <= sampling_level_)
        {
            if (this->sampling_level_ == rhs.sampling_level())
                std::transform(std::cbegin(f_), std::cend(f_),
                               std::cbegin(rhs.f()), std::begin(f_), function);

            else
                for (index_t i = 0; i != this->extents(); ++i)
                    f_[i] = function(f_[i], rhs(sampling()[i]));
        }

        else
        {
            auto lambda = [&lhs = *this, &rhs = rhs, &function = function](const domain_field &x)
                    -> codomain_field
            {
                return function(lhs(x), rhs(x));
            };

            *this = polynomial_point_value{lambda, new_sampling_level, new_order};
        }

        order_ = new_order;
        return;
    }

    // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //
    // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //

    template<Sampling::sampling_c sampling_type, class codomain_field>
    polynomial_point_value<sampling_type, codomain_field> operator-
            (polynomial_point_value<sampling_type, codomain_field> arg)
    {
        std::ranges::for_each(arg.f(), [](codomain_field &a) { a = -a; });

        return arg;
    }

    template<Sampling::sampling_c sampling_type, class codomain_field, class container_type_arg1>
    polynomial_point_value<sampling_type, codomain_field> operator+
            (polynomial_point_value<sampling_type, codomain_field> arg0,
             const polynomial_point_value<sampling_type, codomain_field, container_type_arg1> &arg1)
    {
        return arg0 += arg1;
    }

    template<Sampling::sampling_c sampling_type, class codomain_field, class container_type_arg1>
    polynomial_point_value<sampling_type, codomain_field> operator-
            (polynomial_point_value<sampling_type, codomain_field> arg0,
             const polynomial_point_value<sampling_type, codomain_field, container_type_arg1> &arg1)
    {
        return arg0 -= arg1;
    }

    template<Sampling::sampling_c sampling_type, class codomain_field, class container_type_arg1>
    polynomial_point_value<sampling_type, codomain_field> operator*
            (polynomial_point_value<sampling_type, codomain_field> arg0,
             const polynomial_point_value<sampling_type, codomain_field, container_type_arg1> &arg1)
    {
        return arg0 *= arg1;
    }

    template<Sampling::sampling_c sampling_type, class codomain_field>
    polynomial_point_value<sampling_type, codomain_field> operator+
            (polynomial_point_value<sampling_type, codomain_field> arg0,
             const codomain_field &arg1)
    {
        return arg0 += arg1;
    }

    template<Sampling::sampling_c sampling_type, class codomain_field>
    polynomial_point_value<sampling_type, codomain_field> operator-
            (polynomial_point_value<sampling_type, codomain_field> arg0,
             const codomain_field &arg1)
    {
        return arg0 -= arg1;
    }

    template<Sampling::sampling_c sampling_type, class codomain_field>
    polynomial_point_value<sampling_type, codomain_field> operator*
            (polynomial_point_value<sampling_type, codomain_field> arg0,
             const codomain_field &arg1)
    {
        return arg0 *= arg1;
    }

    template<Sampling::sampling_c sampling_type, class codomain_field>
    polynomial_point_value<sampling_type, codomain_field> operator/
            (polynomial_point_value<sampling_type, codomain_field> arg0,
             const codomain_field &arg1)
    {
        return arg0 /= arg1;
    }

    template<Sampling::sampling_c sampling_type, class codomain_field>
    polynomial_point_value<sampling_type, codomain_field> operator+
            (const codomain_field &arg0,
             polynomial_point_value<sampling_type, codomain_field> arg1)
    {
        return arg1 += arg0;
    }

    template<Sampling::sampling_c sampling_type, class codomain_field>
    polynomial_point_value<sampling_type, codomain_field> operator-
            (const codomain_field &arg0,
             polynomial_point_value<sampling_type, codomain_field> arg1)
    {
        std::ranges::for_each
            (arg1.f(), [&arg0](codomain_field &arg_bis) { arg_bis = arg0 - arg_bis; });

        return arg1;
    }

    template<Sampling::sampling_c sampling_type, class codomain_field>
    polynomial_point_value<sampling_type, codomain_field> operator*
            (const codomain_field &arg0,
             polynomial_point_value<sampling_type, codomain_field> arg1)
    {
        return arg1 *= arg0;
    }

    template<Sampling::sampling_c sampling_type, class codomain_field, class container_type_arg0>
    polynomial_point_value<sampling_type, codomain_field> operator^
            (const polynomial_point_value<sampling_type, codomain_field, container_type_arg0> &arg0,
             index_t arg1)
    {
        using ret_type = polynomial_point_value<sampling_type, codomain_field>;

        if (arg1 == 0)
            return static_cast<ret_type>(1);

        if (arg1 == 1)
            return arg0;

        ret_type ret = arg0 ^ (arg1 / 2);
        ret *= ret;

        if (arg1 % 2)
            ret *= arg0;

        return ret;
    }

    // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //

    template<Sampling::sampling_c sampling_type, class codomain_field, class container_type>
    polynomial_point_value<sampling_type, codomain_field> operator-
            (const polynomial_point_value<sampling_type, codomain_field, container_type> &arg)
        requires(not std::same_as<container_type, std::vector<codomain_field>>)
    {
        return -static_cast<polynomial_point_value<sampling_type, codomain_field>>(arg);
    }

    template<Sampling::sampling_c sampling_type, class codomain_field,
             class container_type_arg0, class container_type_arg1>
    polynomial_point_value<sampling_type, codomain_field> operator+
            (const polynomial_point_value<sampling_type, codomain_field, container_type_arg0> &arg0,
             const polynomial_point_value<sampling_type, codomain_field, container_type_arg1> &arg1)
        requires(not std::same_as<container_type_arg0, std::vector<codomain_field>>)
    {
        return static_cast<polynomial_point_value<sampling_type, codomain_field>>(arg0) + arg1;
    }

    template<Sampling::sampling_c sampling_type, class codomain_field,
             class container_type_arg0, class container_type_arg1>
    polynomial_point_value<sampling_type, codomain_field> operator-
            (const polynomial_point_value<sampling_type, codomain_field, container_type_arg0> &arg0,
             const polynomial_point_value<sampling_type, codomain_field, container_type_arg1> &arg1)
        requires(not std::same_as<container_type_arg0, std::vector<codomain_field>>)
    {
        return static_cast<polynomial_point_value<sampling_type, codomain_field>>(arg0) - arg1;
    }

    template<Sampling::sampling_c sampling_type, class codomain_field,
             class container_type_arg0, class container_type_arg1>
    polynomial_point_value<sampling_type, codomain_field> operator*
            (const polynomial_point_value<sampling_type, codomain_field, container_type_arg0> &arg0,
             const polynomial_point_value<sampling_type, codomain_field, container_type_arg1> &arg1)
        requires(not std::same_as<container_type_arg0, std::vector<codomain_field>>)
    {
        return static_cast<polynomial_point_value<sampling_type, codomain_field>>(arg0) * arg1;
    }

    template<Sampling::sampling_c sampling_type, class codomain_field, class container_type>
    polynomial_point_value<sampling_type, codomain_field> operator+
            (polynomial_point_value<sampling_type, codomain_field, container_type> arg0,
             const codomain_field &arg1)
        requires(not std::same_as<container_type, std::vector<codomain_field>>)
    {
        return static_cast<polynomial_point_value<sampling_type, codomain_field>>(arg0) + arg1;
    }

    template<Sampling::sampling_c sampling_type, class codomain_field, class container_type>
    polynomial_point_value<sampling_type, codomain_field> operator-
            (polynomial_point_value<sampling_type, codomain_field, container_type> arg0,
             const codomain_field &arg1)
        requires(not std::same_as<container_type, std::vector<codomain_field>>)
    {
        return static_cast<polynomial_point_value<sampling_type, codomain_field>>(arg0) - arg1;
    }

    template<Sampling::sampling_c sampling_type, class codomain_field, class container_type>
    polynomial_point_value<sampling_type, codomain_field> operator*
            (polynomial_point_value<sampling_type, codomain_field, container_type> arg0,
             const codomain_field &arg1)
        requires(not std::same_as<container_type, std::vector<codomain_field>>)
    {
        return static_cast<polynomial_point_value<sampling_type, codomain_field>>(arg0) * arg1;
    }

    template<Sampling::sampling_c sampling_type, class codomain_field, class container_type>
    polynomial_point_value<sampling_type, codomain_field> operator/
            (polynomial_point_value<sampling_type, codomain_field, container_type> arg0,
             const codomain_field &arg1)
        requires(not std::same_as<container_type, std::vector<codomain_field>>)
    {
        return static_cast<polynomial_point_value<sampling_type, codomain_field>>(arg0) / arg1;
    }

    template<Sampling::sampling_c sampling_type, class codomain_field, class container_type>
    polynomial_point_value<sampling_type, codomain_field> operator+
            (const codomain_field &arg0,
             polynomial_point_value<sampling_type, codomain_field, container_type> arg1)
        requires(not std::same_as<container_type, std::vector<codomain_field>>)
    {
        return arg0 + static_cast<polynomial_point_value<sampling_type, codomain_field>>(arg1);
    }

    template<Sampling::sampling_c sampling_type, class codomain_field, class container_type>
    polynomial_point_value<sampling_type, codomain_field> operator-
            (const codomain_field &arg0,
             polynomial_point_value<sampling_type, codomain_field, container_type> arg1)
        requires(not std::same_as<container_type, std::vector<codomain_field>>)
    {
        return arg0 - static_cast<polynomial_point_value<sampling_type, codomain_field>>(arg1);
    }

    template<Sampling::sampling_c sampling_type, class codomain_field, class container_type>
    polynomial_point_value<sampling_type, codomain_field> operator*
            (const codomain_field &arg0,
             polynomial_point_value<sampling_type, codomain_field, container_type> arg1)
        requires(not std::same_as<container_type, std::vector<codomain_field>>)
    {
        return arg0 * static_cast<polynomial_point_value<sampling_type, codomain_field>>(arg1);
    }

    // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //
    // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //

    template<Sampling::sampling_c sampling_type>
    class barycentric_interpolation_weights<sampling_type>::factory_apparatus
    {
    private:
        friend class barycentric_interpolation_weights<sampling_type>;

        inline static auto lambda_ = [](index_t level)
        {
            return std::make_unique<barycentric_interpolation_weights>
                    (barycentric_interpolation_weights{level - 1});
        };

        inline static Invocable::invocable_function_hash_table invocable_
                {lambda_, Invocable::invocable<const std::unique_ptr
                        <barycentric_interpolation_weights> &(index_t)>{}};
    };

    // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //

    template<Sampling::sampling_c sampling_type>
    barycentric_interpolation_weights<sampling_type>::
            barycentric_interpolation_weights(index_t level) :
            w_{compute_barycentric_interpolation_weights
                (sampling_type::factory(level).coordinates())} {}

    template<Sampling::sampling_c sampling_type>
    auto barycentric_interpolation_weights<sampling_type>::factory(index_t level)
    -> const barycentric_interpolation_weights &
    {
        return *factory_apparatus::invocable_(level + 1);
    }

    template<Sampling::sampling_c sampling_type>
    void barycentric_interpolation_weights<sampling_type>::clear()
    {
        factory_apparatus::invocable_.clear();
        return;
    }

    // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //

    template<Sampling::sampling_c sampling_type>
    auto barycentric_interpolation_weights<sampling_type>::operator[](index_t arg) const
            -> domain_field
    {
        ASSERT_ASSUME(arg < std::size(w_));

        return w_[arg];
    }

    // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //

    template<class domain_field>
    std::vector<domain_field> compute_barycentric_interpolation_weights
        (std::span<const domain_field> x)
    {
        /* This function computes the barycentric_interpolation_weights
         * w[i] = 1 / product_j (x[i] - x[j]).
         * First, a direct approach is tried. If no underflow or overflow is detected, the
         * function returns. The results are discarded and a more sofisticated result is tried.
         * Note that the barycentric_interpolation_weights are centered around 1 in logarithmic
         * scale. */

        if (std::empty(x))
            return {};

        index_t n = std::size(x);
        std::vector<domain_field> w(n, static_cast<domain_field>(1));

        auto i_w = std::begin(w);

        for (auto i_x = std::cbegin(x); i_x != std::cend(x); ++i_w, ++i_x)
            for (auto j_x = std::cbegin(x); j_x != std::cend(x); ++j_x)
                if (i_x != j_x)
                    *i_w *= *i_x - *j_x;

        real_t min, max;

        i_w = std::begin(w);

        min = std::abs(*i_w);
        max = std::abs(*i_w);

        for (++i_w; i_w != std::end(w); ++i_w)
        {
            min = std::min(min, std::abs(*i_w));
            max = std::max(max, std::abs(*i_w));
        }

        std::ranges::for_each(w, [w_ref = std::sqrt(min * max)](domain_field &i)
                                 {
                                     i = w_ref / i;
                                 });

        if (std::ranges::all_of(w,[](const domain_field &i) -> bool
                                  {
                                      if constexpr (std::is_same_v<domain_field, real_t>)
                                          return std::isfinite(i);
                                      else
                                          return std::isfinite(i.real()) and
                                                 std::isfinite(i.imag());
                                  }))
            return w;

        /* If overflow is detected the barycentric_interpolation_weights are computed in a more
         * cautious way */

        std::ranges::fill(w, static_cast<domain_field>(1));

        std::vector<domain_field> delta(n * n);
        std::vector<std::array<typename std::vector<domain_field>::iterator, 2>> it(n);

        {
            auto it_delta = std::begin(delta);

            for (index_t i = 0; i != n; ++i, it_delta += n)
            {
                it[i][0] = it_delta;
                it[i][1] = it_delta + n - 1;

                std::transform(std::cbegin(x), std::cend(x), it_delta,
                               [&x_i = x[i]](const domain_field &j)
                               {
                                   return x_i - j;
                               });

                *(it_delta + i) = static_cast<domain_field>(1);

                std::sort(it_delta, it_delta + n,
                          [](const domain_field &i,const domain_field &j)
                          {
                              return std::abs(i) < std::abs(j);
                          });
            }
        }

        bool done = false;

        while (!done)
        {
            done = true;

            for (index_t i = 0; i != n; ++i)
            {
                if (it[i][0] > it[i][1])
                    continue;

                done = false;
                std::array<domain_field, 2> a = {w[i] * *it[i][0], w[i] * *it[i][1]};

                if (std::abs(a[0]* a[1]) > 1)
                {
                    w[i] = a[0];
                    ++it[i][0];

                    a[0] = w[i] * *it[i][0];
                    a[1] = w[i] * *it[i][1];

                    min = std::min(min, std::abs(w[i]));
                    max = std::max(max, std::abs(w[i]));
                }

                else
                {
                    w[i] = a[1];
                    --it[i][1];

                    a[0] = w[i] * *it[i][0];
                    a[1] = w[i] * *it[i][1];

                    min = std::min(min, std::abs(w[i]));
                    max = std::max(max, std::abs(w[i]));
                }

                while (it[i][0] <= it[i][1] and std::abs(a[0]) > min and std::abs(a[1]) < max)
                {
                    if (std::abs(a[0]* a[1]) > 1)
                    {
                        w[i] = a[0];

                        ++it[i][0];

                        a[0] = w[i] * *it[i][0];
                        a[1] = w[i] * *it[i][1];

                        min = std::min(min, std::abs(w[i]));
                        max = std::max(max, std::abs(w[i]));
                    }

                    else
                    {
                        w[i] = a[1];
                        --it[i][1];

                        a[0] = w[i] * *it[i][0];
                        a[1] = w[i] * *it[i][1];

                        min = std::min(min, std::abs(w[i]));
                        max = std::max(max, std::abs(w[i]));
                    }
                }

                while (it[i][0] <= it[i][1] and std::abs(a[0]) > min and std::abs(a[0]) < max)
                {
                    w[i] = a[0];
                    ++it[i][0];

                    a[0] = w[i] * *it[i][0];
                    a[1] = w[i] * *it[i][1];

                    min = std::min(min, std::abs(w[i]));
                    max = std::max(max, std::abs(w[i]));
                }

                while (it[i][0] <= it[i][1] and std::abs(a[1]) > min and std::abs(a[1]) < max)
                {
                    w[i] = a[1];
                    --it[i][1];

                    a[0] = w[i] * *it[i][0];
                    a[1] = w[i] * *it[i][1];

                    min = std::min(min, std::abs(w[i]));
                    max = std::max(max, std::abs(w[i]));
                }
            }

            i_w = std::begin(w);

            min = std::abs(*i_w);
            max = std::abs(*i_w);

            for (++i_w; i_w != std::end(w); ++i_w)
            {
                min = std::min(min, std::abs(*i_w));
                max = std::max(max, std::abs(*i_w));
            }

            std::ranges::for_each
            (w, [w_ref = 1 / std::abs(std::sqrt(min) * std::sqrt(max))] (domain_field &i)
                {
                    i *= w_ref;
                });
        }

        std::ranges::for_each(w, [](domain_field &i) { i = 1 / i; });

        return w;
    }

    // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //
    // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //

    template<Sampling::sampling_c sampling_type>
    class quadrature_weights<sampling_type>::factory_apparatus
    {
    private:
        friend class quadrature_weights<sampling_type>;

        inline static auto lambda_ = [](index_t level)
        {
            return std::make_unique<quadrature_weights>(quadrature_weights{level - 1});
        };

        inline static Invocable::invocable_function_hash_table invocable_
                {lambda_, Invocable::invocable<const std::unique_ptr
                        <quadrature_weights> &(index_t)>{}};
    };

    // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //

    template<Sampling::sampling_c sampling_type>
    quadrature_weights<sampling_type>::quadrature_weights(index_t level) :
    w_{compute_quadrature_weights<sampling_type>(level)} {}

    template<Sampling::sampling_c sampling_type>
    auto quadrature_weights<sampling_type>::factory(index_t level) -> const quadrature_weights &
    {
        return *factory_apparatus::invocable_(level + 1);
    }

    template<Sampling::sampling_c sampling_type>
    void quadrature_weights<sampling_type>::clear()
    {
        factory_apparatus::invocable_.clear();
        return;
    }

    template<Sampling::sampling_c sampling_type>
    auto quadrature_weights<sampling_type>::operator[](index_t arg) const -> domain_field
    {
        ASSERT_ASSUME(arg < std::size(w_));
        return w_[arg];
    }

    // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //

    template<>
    class barycentric_interpolation_weights<Sampling::exponential<Sampling::chebyshev_lobatto>>
    {
        // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //
    public:
        using sampling_type = Sampling::exponential<Sampling::chebyshev_lobatto>;
        using domain_field = real_t;

    private:
        class factory_apparatus;

        // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //
    private:
        index_t n_;

        // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //
    public:
        barycentric_interpolation_weights() = delete;
        ~barycentric_interpolation_weights() = default;

        barycentric_interpolation_weights
                (const barycentric_interpolation_weights &) = delete;
        barycentric_interpolation_weights &operator=
                (const barycentric_interpolation_weights &) = delete;

        barycentric_interpolation_weights
                (barycentric_interpolation_weights &&) noexcept = default;
        barycentric_interpolation_weights &operator=
                (barycentric_interpolation_weights &&)noexcept = default;

    private:
        explicit barycentric_interpolation_weights(index_t);

    public:
        static const barycentric_interpolation_weights &factory(index_t);
        static void clear();

        // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //
    public:
        real_t operator[](index_t) const;
    };

    // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //

    class barycentric_interpolation_weights
            <Sampling::exponential<Sampling::chebyshev_lobatto>>::factory_apparatus
    {
    private:
        friend class barycentric_interpolation_weights
                <Sampling::exponential<Sampling::chebyshev_lobatto>>;

        inline static auto lambda_ = [](index_t level)
        {
            return std::make_unique<barycentric_interpolation_weights>
                    (barycentric_interpolation_weights{level - 1});
        };

        inline static Invocable::invocable_function_hash_table invocable_
                {lambda_, Invocable::invocable<const std::unique_ptr
                        <barycentric_interpolation_weights> &(index_t)>{}};
    };

    // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //

    barycentric_interpolation_weights<Sampling::exponential<Sampling::chebyshev_lobatto>>::
    barycentric_interpolation_weights(index_t arg) :
            n_{Sampling::n_points_exponential_sampling(arg)} {}

    auto barycentric_interpolation_weights<Sampling::exponential<Sampling::chebyshev_lobatto>>::
        factory(index_t level)
            -> const barycentric_interpolation_weights &
    {
        return *factory_apparatus::invocable_(level + 1);
    }

    void barycentric_interpolation_weights<Sampling::exponential<Sampling::chebyshev_lobatto>>::
        clear()
    {
        factory_apparatus::invocable_.clear();
        return;
    }

    // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //

    real_t barycentric_interpolation_weights<Sampling::exponential
            <Sampling::chebyshev_lobatto>>::operator[](index_t arg) const
    {
        ASSERT_ASSUME(arg < n_);

        if (arg == 0 or arg == n_ - 1)
            return 1;

        if (arg % 2 == 0)
            return 2;

        return -2;
    }

    // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //
    // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //

    template<Sampling::sampling_c sampling_type, class codomain_field, class container_type>
    codomain_field integrate
            (const polynomial_point_value<sampling_type, codomain_field, container_type> &arg)
    {
        using quadrature_weights_type = quadrature_weights<sampling_type>;
        const quadrature_weights_type &weights =
                quadrature_weights_type::factory(arg.sampling_level());

        codomain_field ret = 0;

        for (index_t i = 0; i != std::size(arg.f()); ++i)
            ret += weights[i] * arg.f(i);

        return ret;
    }
}

// **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //

