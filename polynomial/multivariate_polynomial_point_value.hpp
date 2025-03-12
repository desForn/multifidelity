#pragma once

#include "polynomial_point_value.hpp"

#include "invocable/fwd.hpp"
#include "../sampling/exponential_sampling.hpp"

// **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //

namespace Polynomial
{
    template<class sampling_type_, class codomain_field_, index_t n_variates_>
    class multivariate_polynomial
            <polynomial_point_value<sampling_type_, codomain_field_>, n_variates_>
    {
        static_assert(n_variates_ != 0);

        // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //
    public:
        using univariate_polynomial_type = polynomial_point_value<sampling_type_, codomain_field_>;

        using sampling_type = univariate_polynomial_type::sampling_type;
        using weights_type = univariate_polynomial_type::weights_type;

        using domain_field = univariate_polynomial_type::domain_field;
        using codomain_field = univariate_polynomial_type::codomain_field;

        static constexpr bool polynomial_tag = true;
        static constexpr index_t n_variates = n_variates_;

        using indices_type = std::array<index_t, n_variates>;

        // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //
    private:
        indices_type sampling_level_
                {Utility::uniform_array<index_t, n_variates>(negative_1)};

        indices_type order_
                {Utility::transform_array
                         ([](index_t arg_bis) -> index_t
                          {
                              return arg_bis == negative_1 ?
                                        negative_1 :
                                        sampling_type::n_points(arg_bis) - 1;
                          },
                          sampling_level_)};

        std::vector<codomain_field> f_{};

        // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //
    public:
        multivariate_polynomial() = default;
        ~multivariate_polynomial() = default;

        multivariate_polynomial(const multivariate_polynomial &) = default;
        multivariate_polynomial &operator=(const multivariate_polynomial &) = default;

        multivariate_polynomial(multivariate_polynomial &&) noexcept = default;
        multivariate_polynomial &operator=(multivariate_polynomial &&) noexcept = default;

        explicit multivariate_polynomial(const codomain_field &);

        template<class type>
        explicit multivariate_polynomial(const type &)
        requires(Apparatus::implemented_conversion_c<type, multivariate_polynomial>);

        explicit multivariate_polynomial(univariate_polynomial_type)requires(n_variates == 1);
        multivariate_polynomial(univariate_polynomial_type, index_t);

        template<class function_type>
        multivariate_polynomial(const function_type &, const indices_type &)
        requires(Invocable::invocable_c
                <function_type, codomain_field(std::array<domain_field, n_variates>)>);

        template<class function_type>
        multivariate_polynomial(const function_type &, const indices_type &, const indices_type &)
        requires(Invocable::invocable_c
                <function_type, codomain_field(std::array<domain_field, n_variates>)>);

        multivariate_polynomial(std::vector<codomain_field>, const indices_type &);
        multivariate_polynomial
                (std::vector<codomain_field>, const indices_type &, const indices_type &);

        // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //
    public:
        [[nodiscard]] const indices_type &sampling_level() const;
        [[nodiscard]] index_t sampling_level(index_t) const;

        [[nodiscard]] const indices_type &order() const;
        [[nodiscard]] index_t order(index_t) const;

        [[nodiscard]] std::array<const sampling_type *, n_variates> sampling() const;
        [[nodiscard]] const sampling_type &sampling(index_t) const;

        [[nodiscard]] std::array<const weights_type *, n_variates> weights() const;
        [[nodiscard]] const weights_type &weights(index_t) const;

        [[nodiscard]] indices_type extents() const;
        [[nodiscard]] index_t extents(index_t) const;

        [[nodiscard]] indices_type strides() const;
        [[nodiscard]] index_t strides(index_t) const;

        [[nodiscard]] std::span<codomain_field> f() &;
        [[nodiscard]] std::span<const codomain_field> f() const &;
        [[nodiscard]] std::vector<codomain_field> &&f() &&;

        multivariate_polynomial &shrink();
        multivariate_polynomial &set_order(const indices_type &);

        template<std::integral indices_types>
        multivariate_polynomial &set_order(const std::array<indices_types, n_variates> &);

        template<std::integral... indices_types>
        requires(sizeof...(indices_types) == n_variates_)
        multivariate_polynomial &set_order(indices_types...);

        [[nodiscard]] bool empty() const;
        multivariate_polynomial &clear();

        [[nodiscard]] codomain_field &f(const indices_type &);

        template<std::integral indices_types>
        [[nodiscard]] codomain_field &f(const std::array<indices_types, n_variates> &);

        template<std::integral... indices_types>
        requires(sizeof...(indices_types) == n_variates_)
        [[nodiscard]] codomain_field &f(indices_types...);

        [[nodiscard]] const codomain_field &f(const indices_type &) const;

        template<std::integral indices_types>
        [[nodiscard]] const codomain_field &f(const std::array<indices_types, n_variates> &) const;

        template<std::integral... indices_types>
        requires(sizeof...(indices_types) == n_variates_)
        [[nodiscard]] const codomain_field &f(indices_types...) const;

        [[nodiscard]] std::array<domain_field, n_variates> x(const indices_type &) const;

        template<std::integral indices_types>
        [[nodiscard]] std::array<domain_field, n_variates> x
                (const std::array<indices_types, n_variates> &) const;

        template<std::integral... indices_types>
        requires(sizeof...(indices_types) == n_variates_)
        [[nodiscard]] std::array<domain_field, n_variates> x(indices_types...) const;

        [[nodiscard]] std::array<domain_field, n_variates> w(const indices_type &) const;

        template<std::integral indices_types>
        [[nodiscard]] std::array<domain_field, n_variates> w
                (const std::array<indices_types, n_variates> &) const;

        template<std::integral... indices_types>
        requires(sizeof...(indices_types) == n_variates_)
        [[nodiscard]] std::array<domain_field, n_variates> w(indices_types...) const;

        [[nodiscard]] index_t offset(const indices_type &) const;

        // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //
    public:
        template<Arithmetic::real_complex_c variates_type>
        auto operator()(const std::array<variates_type, n_variates> &) const
        -> std::common_type_t<domain_field, codomain_field, variates_type>;

        template<Arithmetic::real_complex_c... variates_types>
        auto operator()(const variates_types &...) const
        -> std::common_type_t<domain_field, codomain_field, variates_types...>;

        multivariate_polynomial &operator+=(const multivariate_polynomial &);
        multivariate_polynomial &operator-=(const multivariate_polynomial &);
        multivariate_polynomial &operator*=(const multivariate_polynomial &);

        multivariate_polynomial &operator+=(const codomain_field &);
        multivariate_polynomial &operator-=(const codomain_field &);
        multivariate_polynomial &operator*=(const codomain_field &);
        multivariate_polynomial &operator/=(const codomain_field &);

        multivariate_polynomial &operator^=(index_t);

    private:
        template<class function_type>
        void operator_assign
                (const multivariate_polynomial &, const function_type &, const indices_type &);
    };

    // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //

    template<class sampling_type, class codomain_field, index_t n_variates>
    codomain_field integrate(const multivariate_polynomial
            <polynomial_point_value<sampling_type, codomain_field>, n_variates> &);

    // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //
    // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //

    template<class sampling_type, class codomain_field, index_t n_variates>
    multivariate_polynomial<polynomial_point_value<sampling_type, codomain_field>, n_variates>::
    multivariate_polynomial(const codomain_field &arg) :
            sampling_level_{Utility::uniform_array<index_t, n_variates>(0)},
            f_{arg} {}

    template<class sampling_type, class codomain_field, index_t n_variates>
    template<class type>
    multivariate_polynomial<polynomial_point_value<sampling_type, codomain_field>, n_variates>::
    multivariate_polynomial(const type &arg)
    requires(Apparatus::implemented_conversion_c<type, multivariate_polynomial>)
    {
        *this = conversion<type, multivariate_polynomial>::convert(arg);

        return;
    }

    template<class sampling_type, class codomain_field, index_t n_variates>
    multivariate_polynomial<polynomial_point_value<sampling_type, codomain_field>, n_variates>::
    multivariate_polynomial(univariate_polynomial_type p)requires(n_variates == 1) :
            sampling_level_{std::array{p.sampling_level()}},
            order_{std::array{p.order()}},
            f_{std::move(p).f()} {}

    template<class sampling_type, class codomain_field, index_t n_variates>
    multivariate_polynomial<polynomial_point_value<sampling_type, codomain_field>, n_variates>::
    multivariate_polynomial(univariate_polynomial_type p, index_t variate) :
            sampling_level_{Utility::uniform_array<index_t, n_variates>(0)}
    {
        if (p.order() == negative_1)
            this->clear();

        else
        {
            order_[variate] = p.order();
            sampling_level_[variate] = p.sampling_level();

            f_ = std::move(p).f();

            ASSERT_ASSUME(std::size(f_) == extents(0) * strides(0));
        }

        return;
    }

    template<class sampling_type, class codomain_field, index_t n_variates>
    template<class function_type>
    multivariate_polynomial<polynomial_point_value<sampling_type, codomain_field>, n_variates>::
    multivariate_polynomial(const function_type &f, const indices_type &sampling_level)
    requires(Invocable::invocable_c
             <function_type, codomain_field(std::array<domain_field, n_variates>)>) :
            sampling_level_{sampling_level}
    {
        if (std::ranges::find(order_, negative_1) != std::cend(order_))
            this->clear();
        
        else
        {
            const auto sampling_ = sampling();

            std::vector<std::array<domain_field, n_variates>> coordinates;
            coordinates.reserve(extents(0) * strides(0));

            std::array<domain_field, n_variates> x;
            for (const auto &i: Utility::counter(extents()))
            {
                std::transform(std::cbegin(sampling_), std::cend(sampling_),
                               std::cbegin(i),
                               std::begin(x),
                               [](const sampling_type *arg0, index_t arg1) -> domain_field
                               {
                                   return (*arg0)[arg1];
                               });

                coordinates.emplace_back(x);
            }

            f_ = Invocable::evaluate(f, coordinates);
        }

        return;
    }

    template<class sampling_type, class codomain_field, index_t n_variates>
    template<class function_type>
    multivariate_polynomial<polynomial_point_value<sampling_type, codomain_field>, n_variates>::
    multivariate_polynomial(const function_type &f, const indices_type &sampling_level,
                            const indices_type &order)
    requires(Invocable::invocable_c
            <function_type, codomain_field(std::array<domain_field, n_variates>)>) :
            sampling_level_{sampling_level},
            order_{order}
    {
        ASSERT_ASSUME(std::ranges::equal(order_, sampling_level_,
                [](index_t o, index_t l) -> bool
                {
                    return o == negative_1 or o < sampling_type::n_points(l);
                }));

        if (std::ranges::find(order_, negative_1) != std::cend(order_))
            this->clear();

        else
        {
            const auto sampling_ = sampling();

            std::vector<std::array<domain_field, n_variates>> coordinates;
            coordinates.reserve(extents(0) * strides(0));

            std::array<domain_field, n_variates> x;
            for (const auto &i: Utility::counter(extents()))
            {
                std::transform(std::cbegin(sampling_), std::cend(sampling_),
                               std::cbegin(i),
                               std::begin(x),
                               [](const sampling_type *arg0, index_t arg1) -> domain_field
                               {
                                   return (*arg0)[arg1];
                               });

                coordinates.emplace_back(x);
            }

            f_ = Invocable::evaluate(f, coordinates);
        }

        return;
    }

    template<class sampling_type, class codomain_field, index_t n_variates>
    multivariate_polynomial<polynomial_point_value<sampling_type, codomain_field>, n_variates>::
    multivariate_polynomial(std::vector<codomain_field> f, const indices_type &sampling_level) :
            sampling_level_{sampling_level},
            f_{std::move(f)}
    {
        ASSERT_ASSUME(std::size(f_) == extents(0) * strides(0));

        return;
    }

    template<class sampling_type, class codomain_field, index_t n_variates>
    multivariate_polynomial<polynomial_point_value<sampling_type, codomain_field>, n_variates>::
    multivariate_polynomial(std::vector<codomain_field> f, const indices_type &sampling_level,
                            const indices_type &order) :
            sampling_level_{sampling_level},
            order_{order},
            f_{std::move(f)}
    {
        ASSERT_ASSUME(std::size(f_) == extents(0) * strides(0));
        ASSERT_ASSUME(std::ranges::equal(order_, sampling_level_,
                                          [](index_t i, index_t j)
                                          {
                                              return i == negative_1 or
                                                     i < sampling_type::n_points(j);
                                          }));

        return;
    }

    // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //

    template<class sampling_type, class codomain_field, index_t n_variates>
    auto multivariate_polynomial
            <polynomial_point_value<sampling_type, codomain_field>, n_variates>::sampling_level()
    const -> const indices_type &
    {
        return sampling_level_;
    }

    template<class sampling_type, class codomain_field, index_t n_variates>
    index_t multivariate_polynomial
            <polynomial_point_value<sampling_type, codomain_field>, n_variates>::sampling_level
            (index_t arg) const
    {
        return sampling_level_[arg];
    }

    template<class sampling_type, class codomain_field, index_t n_variates>
    auto multivariate_polynomial
            <polynomial_point_value<sampling_type, codomain_field>, n_variates>::order() const
    -> const indices_type &
    {
        return order_;
    }

    template<class sampling_type, class codomain_field, index_t n_variates>
    index_t multivariate_polynomial
            <polynomial_point_value<sampling_type, codomain_field>, n_variates>::order
            (index_t arg) const
    {
        return order_[arg];
    }

    template<class sampling_type, class codomain_field, index_t n_variates>
    auto multivariate_polynomial
            <polynomial_point_value<sampling_type, codomain_field>, n_variates>::sampling() const
    -> std::array<const sampling_type *, n_variates>
    {
        ASSERT_ASSUME(sampling_level_.front() != negative_1);

        return Utility::transform_array
                ([](index_t arg_bis) -> const sampling_type *
                 {
                     return &sampling_type::factory(arg_bis);
                 },
                 sampling_level_);
    }

    template<class sampling_type, class codomain_field, index_t n_variates>
    auto multivariate_polynomial
            <polynomial_point_value<sampling_type, codomain_field>, n_variates>::sampling
            (index_t arg) const -> const sampling_type &
    {
        ASSERT_ASSUME(sampling_level_[arg] != negative_1);
        return sampling_type::factory(sampling_level_[arg]);
    }

    template<class sampling_type, class codomain_field, index_t n_variates>
    auto multivariate_polynomial
            <polynomial_point_value<sampling_type, codomain_field>, n_variates>::weights() const
    -> std::array<const weights_type *, n_variates>
    {
        ASSERT_ASSUME(sampling_level_.front() != negative_1);

        return Utility::transform_array
                ([](index_t arg_bis) -> const weights_type *
                 {
                     return &weights_type::factory(arg_bis);
                 },
                 sampling_level_);
    }

    template<class sampling_type, class codomain_field, index_t n_variates>
    auto multivariate_polynomial
            <polynomial_point_value<sampling_type, codomain_field>, n_variates>::weights
            (index_t arg) const -> const weights_type &
    {
        ASSERT_ASSUME(sampling_level_[arg] != negative_1);
        return weights_type::factory(sampling_level_[arg]);
    }

    template<class sampling_type, class codomain_field, index_t n_variates>
    auto multivariate_polynomial
            <polynomial_point_value<sampling_type, codomain_field>, n_variates>::extents() const
    -> indices_type
    {
        if(sampling_level_.front() == negative_1)
            return Utility::uniform_array<index_t, n_variates>(0);

        return Utility::transform_array
                        ([](index_t arg_bis) -> index_t
                         {
                             return sampling_type::n_points(arg_bis);
                         },
                         sampling_level_);
    }

    template<class sampling_type, class codomain_field, index_t n_variates>
    index_t multivariate_polynomial
            <polynomial_point_value<sampling_type, codomain_field>, n_variates>::extents
            (index_t arg) const
    {
        if (sampling_level_[arg] == negative_1)
            return 0;

        return sampling_type::n_points(sampling_level_[arg]);
    }

    template<class sampling_type, class codomain_field, index_t n_variates>
    auto multivariate_polynomial
            <polynomial_point_value<sampling_type, codomain_field>, n_variates>::strides() const
    -> indices_type
    {
        return Utility::strides(extents());
    }

    template<class sampling_type, class codomain_field, index_t n_variates>
    index_t multivariate_polynomial
            <polynomial_point_value<sampling_type, codomain_field>, n_variates>::strides
            (index_t arg) const
    {
        return strides()[arg];
    }

    template<class sampling_type, class codomain_field, index_t n_variates>
    auto multivariate_polynomial
            <polynomial_point_value<sampling_type, codomain_field>, n_variates>::f() &
    -> std::span<codomain_field>
    {
        return {std::begin(f_), std::end(f_)};
    }

    template<class sampling_type, class codomain_field, index_t n_variates>
    auto multivariate_polynomial
            <polynomial_point_value<sampling_type, codomain_field>, n_variates>::f() const &
    -> std::span<const codomain_field>
    {
        return {std::cbegin(f_), std::cend(f_)};
    }

    template<class sampling_type, class codomain_field, index_t n_variates>
    auto multivariate_polynomial
            <polynomial_point_value<sampling_type, codomain_field>, n_variates>::f() &&
    -> std::vector<codomain_field> &&
    {
        return std::move(f_);
    }

    template<class sampling_type, class codomain_field, index_t n_variates>
    auto multivariate_polynomial
            <polynomial_point_value<sampling_type, codomain_field>, n_variates>::shrink()
    -> multivariate_polynomial &
    {
        if (order_.front() == negative_1)
            return this->clear();

        indices_type new_sampling_level =
                Utility::transform_array
                        ([](index_t arg_bis) -> index_t
                         {
                             return sampling_type::required_level(arg_bis + 1);
                         },
                         order_);

        if (new_sampling_level != sampling_level_)
            *this = multivariate_polynomial{*this, new_sampling_level};

        return *this;
    }

    template<class sampling_type, class codomain_field, index_t n_variates>
    auto multivariate_polynomial
            <polynomial_point_value<sampling_type, codomain_field>, n_variates>::set_order
            (const indices_type &arg) -> multivariate_polynomial &
    {
        if (std::ranges::find(arg, negative_1) != std::cend(arg))
            order_ = std::ranges::fill(order_, negative_1);

        else
        {
            ASSERT_ASSUME(std::ranges::equal(arg, extents(),
                                             [](index_t i, index_t j) { return i < j; }));

            order_ = arg;
        }

        return *this;
    }

    template<class sampling_type, class codomain_field, index_t n_variates>
    template<std::integral indices_types>
    auto multivariate_polynomial
            <polynomial_point_value<sampling_type, codomain_field>, n_variates>::set_order
            (const std::array<indices_types, n_variates> &arg) -> multivariate_polynomial &
    {
        return this->set_order(Utility::static_cast_array<index_t>(arg));
    }

    template<class sampling_type, class codomain_field, index_t n_variates>
    template<std::integral... indices_types>
    requires(sizeof...(indices_types) == n_variates)
    auto multivariate_polynomial
            <polynomial_point_value<sampling_type, codomain_field>, n_variates>::set_order
            (indices_types... args) -> multivariate_polynomial &
    {
        return this->set_order(std::array<index_t, n_variates>(args...));
    }

    template<class sampling_type, class codomain_field, index_t n_variates>
    bool multivariate_polynomial
            <polynomial_point_value<sampling_type, codomain_field>, n_variates>::empty() const
    {
        return std::empty(f_);
    }

    template<class sampling_type, class codomain_field, index_t n_variates>
    auto multivariate_polynomial
            <polynomial_point_value<sampling_type, codomain_field>, n_variates>::clear()
    -> multivariate_polynomial &
    {
        sampling_level_ = {Utility::uniform_array<index_t, n_variates>(negative_1)};
        order_ = {Utility::uniform_array<index_t, n_variates>(negative_1)};

        f_.clear();

        return *this;
    }

    template<class sampling_type, class codomain_field, index_t n_variates>
    auto multivariate_polynomial
            <polynomial_point_value<sampling_type, codomain_field>, n_variates>::f
            (const indices_type &arg) -> codomain_field &
    {
        return f_[offset(arg)];
    }

    template<class sampling_type, class codomain_field, index_t n_variates>
    template<std::integral indices_types>
    auto multivariate_polynomial
            <polynomial_point_value<sampling_type, codomain_field>, n_variates>::f
            (const std::array<indices_types, n_variates> &arg) -> codomain_field &
    {
        return (*this)[Utility::static_cast_array<index_t>(arg)];
    }

    template<class sampling_type, class codomain_field, index_t n_variates>
    template<std::integral... indices_types>
    requires(sizeof...(indices_types) == n_variates)
    auto multivariate_polynomial
            <polynomial_point_value<sampling_type, codomain_field>, n_variates>::f
            (indices_types... args) -> codomain_field &
    {
        return (*this)[indices_type{args...}];
    }

    template<class sampling_type, class codomain_field, index_t n_variates>
    auto multivariate_polynomial
            <polynomial_point_value<sampling_type, codomain_field>, n_variates>::f
            (const indices_type &arg) const -> const codomain_field &
    {
        return f_[offset(arg)];
    }

    template<class sampling_type, class codomain_field, index_t n_variates>
    template<std::integral indices_types>
    auto multivariate_polynomial
            <polynomial_point_value<sampling_type, codomain_field>, n_variates>::f
            (const std::array<indices_types, n_variates> &arg) const -> const codomain_field &
    {
        return (*this)[Utility::static_cast_array<index_t>(arg)];
    }

    template<class sampling_type, class codomain_field, index_t n_variates>
    template<std::integral... indices_types>
    requires(sizeof...(indices_types) == n_variates)
    auto multivariate_polynomial
            <polynomial_point_value<sampling_type, codomain_field>, n_variates>::f
            (indices_types... args) const -> const codomain_field &
    {
        return (*this)[indices_type{args...}];
    }

    template<class sampling_type, class codomain_field, index_t n_variates>
    auto multivariate_polynomial
            <polynomial_point_value<sampling_type, codomain_field>, n_variates>::x
            (const indices_type &arg) const -> std::array<domain_field, n_variates>
    {
        ASSERT_ASSUME(std::ranges::equal
                              (arg, extents(),
                               [](index_t arg_bis0, index_t arg_bis1)
                               {
                                   return arg_bis0 < arg_bis1;
                               }));

       const auto sampling_ = sampling();

        std::array<domain_field, n_variates> ret;

        std::transform(std::cbegin(sampling_), std::cend(sampling_),
                       std::cbegin(arg),
                       std::begin(ret),
                       [](const sampling_type *arg0, index_t arg1) -> domain_field
                       {
                           return (*arg0)[arg1];
                       });

        return ret;
    }

    template<class sampling_type, class codomain_field, index_t n_variates>
    template<std::integral indices_types>
    auto multivariate_polynomial
            <polynomial_point_value<sampling_type, codomain_field>, n_variates>::x
            (const std::array<indices_types, n_variates> &arg) const
    -> std::array<domain_field, n_variates>
    {
        return this->x(Utility::static_cast_array<index_t>(arg));
    }

    template<class sampling_type, class codomain_field, index_t n_variates>
    template<std::integral... indices_types>
    requires(sizeof...(indices_types) == n_variates)
    auto multivariate_polynomial
            <polynomial_point_value<sampling_type, codomain_field>, n_variates>::x
            (indices_types... args) const -> std::array<domain_field, n_variates>
    {
        return this->x(indices_type{args...});
    }

    template<class sampling_type, class codomain_field, index_t n_variates>
    auto multivariate_polynomial
            <polynomial_point_value<sampling_type, codomain_field>, n_variates>::w
            (const indices_type &arg) const -> std::array<domain_field, n_variates>
    {
        ASSERT_ASSUME(std::ranges::equal
                              (arg, extents(),
                               [](index_t arg_bis0, index_t arg_bis1)
                               {
                                   return arg_bis0 < arg_bis1;
                               }));

        const auto weights_ = weights();

        domain_field ret = 1;

        for (index_t i = 0; i != n_variates; ++i)
            ret *= (*weights_[i])[arg[i]];

        return ret;
    }

    template<class sampling_type, class codomain_field, index_t n_variates>
    template<std::integral indices_types>
    auto multivariate_polynomial
            <polynomial_point_value<sampling_type, codomain_field>, n_variates>::w
            (const std::array<indices_types, n_variates> &arg) const
    -> std::array<domain_field, n_variates>
    {
        return this->w(Utility::static_cast_array<index_t>(arg));
    }

    template<class sampling_type, class codomain_field, index_t n_variates>
    template<std::integral... indices_types>
    requires(sizeof...(indices_types) == n_variates)
    auto multivariate_polynomial
            <polynomial_point_value<sampling_type, codomain_field>, n_variates>::w
            (indices_types... args) const -> std::array<domain_field, n_variates>
    {
        return this->w(indices_type{args...});
    }

    template<class sampling_type, class codomain_field, index_t n_variates>
    index_t multivariate_polynomial
            <polynomial_point_value<sampling_type, codomain_field>, n_variates>::offset
            (const indices_type &arg) const
    {
        ASSERT_ASSUME(std::ranges::equal
                              (arg, extents(),
                               [](index_t arg_bis0, index_t arg_bis1)
                               {
                                   return arg_bis0 < arg_bis1;
                               }));

        indices_type strides_ = strides();

        return std::inner_product
                (std::cbegin(strides_), std::cend(strides_), std::cbegin(arg),
                 static_cast<index_t>(0));
    }

    // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //

    template<class sampling_type, class codomain_field, index_t n_variates>
    template<Arithmetic::real_complex_c variates_type>
    auto multivariate_polynomial
            <polynomial_point_value<sampling_type, codomain_field>, n_variates>::operator()
            (const std::array<variates_type, n_variates> &x) const
    -> std::common_type_t<domain_field, codomain_field, variates_type>
    {
        using ret_type = std::common_type_t<domain_field, codomain_field, variates_type>;

        if (order_.front() == negative_1)
            return static_cast<ret_type>(0);

        const indices_type extents_ = extents();

        std::vector<ret_type> a;

        a.reserve(std::size(f_) / extents_.back());

        for (auto i = std::cbegin(f_); i != std::cend(f_); i += extents_.back())
        {
            polynomial_point_value<sampling_type, codomain_field, std::span<const codomain_field>> p
                    {std::span<const codomain_field>{i, i + extents_.back()},
                     sampling_level_.back(), order_.back()};

            a.emplace_back(p(x.back()));
        }

        if constexpr (n_variates != 1)
        {
            std::vector<ret_type> b;
            b.reserve(std::size(a) / extents_[n_variates - 2]);

            for (index_t d = n_variates - 2; d != negative_1; --d)
            {
                for (auto i = std::cbegin(a); i != std::cend(a); i += extents_[d])
                {
                    polynomial_point_value<sampling_type, ret_type, std::span<const ret_type>> p
                            {std::span<const ret_type>{i, i + extents_[d]},
                             sampling_level_[d], order_[d]};

                    b.emplace_back(p(x[d]));
                }

                std::swap(a, b);
                b.clear();
            }
        }

        return a.front();
    }

    template<class sampling_type, class codomain_field, index_t n_variates>
    template<Arithmetic::real_complex_c... variates_types>
    auto multivariate_polynomial
            <polynomial_point_value<sampling_type, codomain_field>, n_variates>::operator()
            (const variates_types &...args) const
    -> std::common_type_t<domain_field, codomain_field, variates_types...>
    {
        using type = std::common_type_t<variates_types...>;

        return (*this)(std::array<type, n_variates>{args...});
    }

    template<class sampling_type, class codomain_field, index_t n_variates>
    auto multivariate_polynomial
            <polynomial_point_value<sampling_type, codomain_field>, n_variates>::operator+=
            (const multivariate_polynomial &arg) -> multivariate_polynomial &
    {
        if (std::empty(arg.f_))
            return *this;

        if (std::empty(f_))
            return *this = arg;

        operator_assign(arg, std::plus{}, Utility::max_array(order_, arg.order_));

        return *this;
    }

    template<class sampling_type, class codomain_field, index_t n_variates>
    auto multivariate_polynomial
            <polynomial_point_value<sampling_type, codomain_field>, n_variates>::operator-=
            (const multivariate_polynomial &arg) -> multivariate_polynomial &
    {
        if (std::empty(arg.f_))
            return *this;

        if (std::empty(f_))
            return *this = -arg;

        operator_assign(arg, std::minus{}, Utility::max_array(order_, arg.order_));

        return *this;
    }

    template<class sampling_type, class codomain_field, index_t n_variates>
    auto multivariate_polynomial
            <polynomial_point_value<sampling_type, codomain_field>, n_variates>::operator*=
            (const multivariate_polynomial &arg) -> multivariate_polynomial &
    {
        if (std::empty(f_))
            return *this;

        if (std::empty(arg.f_))
            return this->clear();

        operator_assign(arg, std::multiplies{}, Utility::sum_array(order_, arg.order_));

        return *this;
    }

    template<class sampling_type, class codomain_field, index_t n_variates>
    auto multivariate_polynomial
            <polynomial_point_value<sampling_type, codomain_field>, n_variates>::operator+=
            (const codomain_field &arg) -> multivariate_polynomial &
    {
        std::ranges::for_each(f_, [&arg](codomain_field &arg_bis) { arg_bis += arg; });

        return *this;
    }

    template<class sampling_type, class codomain_field, index_t n_variates>
    auto multivariate_polynomial
            <polynomial_point_value<sampling_type, codomain_field>, n_variates>::operator-=
            (const codomain_field &arg) -> multivariate_polynomial &
    {
        std::ranges::for_each(f_, [&arg](codomain_field &arg_bis) { arg_bis -= arg; });

        return *this;
    }

    template<class sampling_type, class codomain_field, index_t n_variates>
    auto multivariate_polynomial
            <polynomial_point_value<sampling_type, codomain_field>, n_variates>::operator*=
            (const codomain_field &arg) -> multivariate_polynomial &
    {
        std::ranges::for_each(f_, [&arg](codomain_field &arg_bis) { arg_bis *= arg; });

        return *this;
    }

    template<class sampling_type, class codomain_field, index_t n_variates>
    auto multivariate_polynomial
            <polynomial_point_value<sampling_type, codomain_field>, n_variates>::operator/=
            (const codomain_field &arg) -> multivariate_polynomial &
    {
        std::ranges::for_each(f_, [&arg](codomain_field &arg_bis) { arg_bis /= arg; });

        return *this;
    }

    template<class sampling_type, class codomain_field, index_t n_variates>
    auto multivariate_polynomial
            <polynomial_point_value<sampling_type, codomain_field>, n_variates>::operator^=
            (index_t arg) -> multivariate_polynomial &
    {
        return *this = *this ^ arg;
    }

    template<class sampling_type, class codomain_field, index_t n_variates>
    template<class function_type>
    void multivariate_polynomial
            <polynomial_point_value<sampling_type, codomain_field>, n_variates>::operator_assign
            (const multivariate_polynomial &p, const function_type &function,
             const indices_type &new_order)
    {
        indices_type new_sampling_level =
                {Utility::transform_array
                         ([](index_t arg_bis) -> index_t
                          {
                              return sampling_type::required_level(arg_bis + 1);
                          },
                          new_order)};

        bool sampling_level_is_enough =
                std::ranges::equal(new_sampling_level, sampling_level_,
                                   [](index_t arg_bis0, index_t arg_bis1)
                                   {
                                       return arg_bis0 <= arg_bis1;
                                   });

        if (sampling_level_is_enough)
        {
            if (this->sampling_level_ == p.sampling_level_)
                std::transform(std::cbegin(f_), std::cend(f_),
                               std::cbegin(p.f_), std::begin(f_), function);

            else
                for (const auto &i: Utility::counter(extents()))
                {
                    const indices_type strides_ = strides();

                    index_t offset = std::inner_product(std::cbegin(strides_), std::cend(strides_),
                                                        std::cbegin(i), static_cast<index_t>(0));

                    std::array<domain_field, n_variates> x_ = x(i);

                    f_[offset] = function((*this)(x_), p(x_));
                }

            order_ = new_order;
        }

        else
        {
            auto lambda = [&arg_bis0 = *this, &arg_bis1 = p, &function = function]
                    (const std::array<domain_field, n_variates> &x) -> codomain_field
            {
                return function(arg_bis0(x), arg_bis1(x));
            };

            *this = multivariate_polynomial{lambda, new_sampling_level, new_order};
        }

        return;
    }

    // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //
    // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //

    template<class sampling_type, class codomain_field, index_t n_variates>
    multivariate_polynomial<polynomial_point_value<sampling_type, codomain_field>, n_variates>
    operator-(multivariate_polynomial
              <polynomial_point_value<sampling_type, codomain_field>, n_variates> arg)
    {
        std::ranges::for_each
                (arg.f(), [](codomain_field &arg_bis) { return arg_bis = -arg_bis; });

        return arg;
    }

    template<class sampling_type, class codomain_field, index_t n_variates>
    multivariate_polynomial<polynomial_point_value<sampling_type, codomain_field>, n_variates>
    operator+(multivariate_polynomial
              <polynomial_point_value<sampling_type, codomain_field>, n_variates> arg0,
              const multivariate_polynomial
              <polynomial_point_value<sampling_type, codomain_field>, n_variates> &arg1)
    {
        return arg0 += arg1;
    }

    template<class sampling_type, class codomain_field, index_t n_variates>
    multivariate_polynomial<polynomial_point_value<sampling_type, codomain_field>, n_variates>
    operator-(multivariate_polynomial
              <polynomial_point_value<sampling_type, codomain_field>, n_variates> arg0,
              const multivariate_polynomial
              <polynomial_point_value<sampling_type, codomain_field>, n_variates> &arg1)
    {
        return arg0 -= arg1;
    }

    template<class sampling_type, class codomain_field, index_t n_variates>
    multivariate_polynomial<polynomial_point_value<sampling_type, codomain_field>, n_variates>
    operator*(multivariate_polynomial
              <polynomial_point_value<sampling_type, codomain_field>, n_variates> arg0,
              const multivariate_polynomial
              <polynomial_point_value<sampling_type, codomain_field>, n_variates> &arg1)
    {
        return arg0 *= arg1;
    }

    template<class sampling_type, class codomain_field, index_t n_variates>
    multivariate_polynomial<polynomial_point_value<sampling_type, codomain_field>, n_variates>
    operator+(multivariate_polynomial
              <polynomial_point_value<sampling_type, codomain_field>, n_variates> arg0,
              const typename multivariate_polynomial
                      <polynomial_point_value<sampling_type, codomain_field>, n_variates>::
              codomain_field &arg1)
    {
        return arg0 += arg1;
    }

    template<class sampling_type, class codomain_field, index_t n_variates>
    multivariate_polynomial<polynomial_point_value<sampling_type, codomain_field>, n_variates>
    operator-(multivariate_polynomial
              <polynomial_point_value<sampling_type, codomain_field>, n_variates> arg0,
              const typename multivariate_polynomial
                      <polynomial_point_value<sampling_type, codomain_field>, n_variates>::
              codomain_field &arg1)
    {
        return arg0 -= arg1;
    }

    template<class sampling_type, class codomain_field, index_t n_variates>
    multivariate_polynomial<polynomial_point_value<sampling_type, codomain_field>, n_variates>
    operator*(multivariate_polynomial
              <polynomial_point_value<sampling_type, codomain_field>, n_variates> arg0,
              const typename multivariate_polynomial
                      <polynomial_point_value<sampling_type, codomain_field>, n_variates>::
                      codomain_field &arg1)
    {
        return arg0 *= arg1;
    }

    template<class sampling_type, class codomain_field, index_t n_variates>
    multivariate_polynomial<polynomial_point_value<sampling_type, codomain_field>, n_variates>
    operator/(multivariate_polynomial
              <polynomial_point_value<sampling_type, codomain_field>, n_variates> arg0,
              const typename multivariate_polynomial
                      <polynomial_point_value<sampling_type, codomain_field>, n_variates>::
              codomain_field &arg1)
    {
        return arg0 /= arg1;
    }

    template<class sampling_type, class codomain_field, index_t n_variates>
    multivariate_polynomial<polynomial_point_value<sampling_type, codomain_field>, n_variates>
    operator+(const typename multivariate_polynomial
                      <polynomial_point_value<sampling_type, codomain_field>, n_variates>::
              codomain_field &arg0,
              multivariate_polynomial
              <polynomial_point_value<sampling_type, codomain_field>, n_variates> arg1)
    {
        return arg1 += arg0;
    }

    template<class sampling_type, class codomain_field, index_t n_variates>
    multivariate_polynomial<polynomial_point_value<sampling_type, codomain_field>, n_variates>
    operator-(const typename multivariate_polynomial
                      <polynomial_point_value<sampling_type, codomain_field>, n_variates>::
              codomain_field &arg0,
              multivariate_polynomial
              <polynomial_point_value<sampling_type, codomain_field>, n_variates> arg1)
    {
        std::ranges::for_each
                (arg1.f(), [&arg0](codomain_field &arg_bis) { arg_bis = arg0 - arg_bis; });

        return arg1;
    }

    template<class sampling_type, class codomain_field, index_t n_variates>
    multivariate_polynomial<polynomial_point_value<sampling_type, codomain_field>, n_variates>
    operator*(const typename multivariate_polynomial
                      <polynomial_point_value<sampling_type, codomain_field>, n_variates>::
              codomain_field &arg0,
              multivariate_polynomial
              <polynomial_point_value<sampling_type, codomain_field>, n_variates> arg1)
    {
        return arg1 *= arg0;
    }

    template<class sampling_type, class codomain_field, index_t n_variates>
    multivariate_polynomial<polynomial_point_value<sampling_type, codomain_field>, n_variates>
    operator^(multivariate_polynomial
              <polynomial_point_value<sampling_type, codomain_field>, n_variates> arg0,
              index_t arg1)
    {
        using ret_type = multivariate_polynomial
                <polynomial_point_value<sampling_type, codomain_field>, n_variates>;

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
    // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //

    template<class sampling_type, class codomain_field, index_t n_variates>
    codomain_field integrate(const multivariate_polynomial
            <polynomial_point_value<sampling_type, codomain_field>, n_variates> &arg)
    {
        using quadrature_weights = quadrature_weights<sampling_type>;

        if (arg.sampling_level(0) == negative_1)
            return 0;

        std::array<const quadrature_weights *, n_variates> weights{Utility::transform_array
        ([](index_t sampling_level) -> const quadrature_weights *
         {
            return &quadrature_weights::factory(sampling_level);
         },
        arg.sampling_level())};

        codomain_field ret = 0;

        for (const auto& i : Utility::counter(arg.extents()))
        {
            codomain_field c = arg.f(i);

            for (index_t j = 0; j != n_variates; ++j)
                c *= (*weights[j])[i[j]];

            ret += c;
        }

        return ret;
    }
}
