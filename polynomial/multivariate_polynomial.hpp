#pragma once

#include "polynomial_base.hpp"
#include "invocable/invocable_function_hash_table.hpp"

// **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //

namespace Polynomial
{
    template<class univariate_polynomial_traits, index_t n_variates_>
    class multivariate_polynomial<Apparatus::polynomial_base<univariate_polynomial_traits>, n_variates_>
    {
        static_assert(n_variates_ != 0);

        // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //
    public:
        using univariate_polynomial_type = Apparatus::polynomial_base<univariate_polynomial_traits>;

        using coefficients_type = univariate_polynomial_type::coefficients_type;

        static constexpr bool polynomial_tag = true;
        static constexpr index_t n_variates = n_variates_;

        using indices_type = std::array<index_t, n_variates>;

        // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //
    private:
        std::vector<coefficients_type> coefficients_{};
        indices_type extents_{Utility::uniform_array<index_t, n_variates>(0)};
        indices_type strides_{Utility::uniform_array<index_t, n_variates>(0)};

        // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //
    public:
        multivariate_polynomial() = default;
        ~multivariate_polynomial() = default;

        multivariate_polynomial(const multivariate_polynomial &) = default;
        multivariate_polynomial &operator=(const multivariate_polynomial &) = default;

        multivariate_polynomial(multivariate_polynomial &&) noexcept = default;
        multivariate_polynomial &operator=(multivariate_polynomial &&) noexcept = default;

        template<class type>
        explicit multivariate_polynomial(const type &)
        requires(std::convertible_to<type, coefficients_type> and
                 !Apparatus::implemented_conversion_c<type, multivariate_polynomial>);

        template<class type>
        explicit multivariate_polynomial(const type &arg)
        requires(Apparatus::implemented_conversion_c<type, multivariate_polynomial>);

        explicit multivariate_polynomial(univariate_polynomial_type)requires(n_variates == 1);

        multivariate_polynomial(univariate_polynomial_type, index_t);

        explicit multivariate_polynomial(std::vector<coefficients_type>, const indices_type &);

        // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //
    public:
        [[nodiscard]] std::span<coefficients_type> coefficients() &;
        [[nodiscard]] std::span<const coefficients_type> coefficients() const &;
        [[nodiscard]] std::vector<coefficients_type> &&coefficients() &&;

        [[nodiscard]] coefficients_type &operator[](const indices_type &);
        [[nodiscard]] const coefficients_type &operator[](const indices_type &) const;

        template<std::integral indices_types>
        [[nodiscard]] coefficients_type &operator[](const std::array<indices_types, n_variates> &);

        template<std::integral indices_types>
        [[nodiscard]] const coefficients_type &operator[]
            (const std::array<indices_types, n_variates> &) const;

        template<std::integral... indices_types>
        requires(sizeof...(indices_types) == n_variates_)
        [[nodiscard]] coefficients_type &operator[](indices_types...);

        template<std::integral... indices_types>
        requires(sizeof...(indices_types) == n_variates_)
        [[nodiscard]] const coefficients_type &operator[](indices_types...) const;

        bool empty() const;
        multivariate_polynomial &clear();

        [[nodiscard]] indices_type order() const;
        [[nodiscard]] index_t order(index_t) const;
        multivariate_polynomial &set_order();
        multivariate_polynomial &set_order(const indices_type &);

        template<std::integral indices_types>
        multivariate_polynomial &set_order(const std::array<indices_types, n_variates> &);

        template<class... indices_types>
        requires(sizeof...(indices_types) == n_variates_ and (std::integral<indices_types> and ...))
        multivariate_polynomial &set_order(indices_types...);

        multivariate_polynomial &set_extents();
        multivariate_polynomial &set_extents(const indices_type &);

        template<std::integral indices_types>
        multivariate_polynomial &set_extents(const std::array<indices_types, n_variates> &);

        template<class... indices_types>
        requires(sizeof...(indices_types) == n_variates_ and (std::integral<indices_types> and ...))
        multivariate_polynomial &set_extents(indices_types...);

        [[nodiscard]] const indices_type &extents() const;
        [[nodiscard]] index_t extents(index_t) const;
        [[nodiscard]] const indices_type &strides() const;
        [[nodiscard]] index_t strides(index_t) const;

        [[nodiscard]] index_t offset(const indices_type &) const;

    private:
        bool check_order
            (const indices_type &extents, index_t i, indices_type aux = {}, index_t j = 0) const;

        void reset_coefficients
            (std::vector<coefficients_type> &new_coefficients, const indices_type &new_extents,
             const indices_type &new_strides, index_t i_old = 0, index_t i_new = 0,
             index_t i_var = 0);

        // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //
    public:
        template<class... variate_types>
        auto operator()(const variate_types &...x) const -> decltype(Apparatus::evaluate(*this, x...));

        multivariate_polynomial &operator+=(const multivariate_polynomial &);
        multivariate_polynomial &operator-=(const multivariate_polynomial &);
        multivariate_polynomial &operator*=(const multivariate_polynomial &);

        multivariate_polynomial &operator+=(const coefficients_type &);
        multivariate_polynomial &operator-=(const coefficients_type &);
        multivariate_polynomial &operator*=(const coefficients_type &);
        multivariate_polynomial &operator/=(const coefficients_type &);

        multivariate_polynomial &operator^=(index_t);

        bool operator==(const multivariate_polynomial &) const = default;
        bool operator!=(const multivariate_polynomial &) const = default;

    private:
        template<class invocable_type>
        void elementwise_op(const multivariate_polynomial &, const invocable_type &);

        template<class invocable_type>
        void elementwise_op_apparatus
            (const multivariate_polynomial &arg1, const invocable_type &invocable,
             index_t i_var = 0, index_t offset_this = 0, index_t offset_arg1 = 0);

        // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //
    public:
        static multivariate_polynomial basis_polynomial(const indices_type &);

        template<std::integral indices_types>
        static multivariate_polynomial basis_polynomial
            (const std::array<indices_types, n_variates> &);

        template<class... indices_types>
        requires(sizeof...(indices_types) == n_variates_ and
                 (std::is_convertible_v<indices_types, index_t> and ...))
        static multivariate_polynomial basis_polynomial(indices_types...);
    };

    // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //

    namespace Apparatus
    {
        template<index_t n_variates, class univariate_polynomial_traits, class... variates_types>
        requires(sizeof...(variates_types) == n_variates and
                 (std::is_same_v<variates_types,
                                 typename polynomial_base<univariate_polynomial_traits>::
                                 coefficients_type> and ...))
        polynomial_base<univariate_polynomial_traits>::coefficients_type evaluate
            (const multivariate_polynomial
                   <polynomial_base<univariate_polynomial_traits>, n_variates> &p,
             const variates_types &...x)
        {
            using univariate_polynomial_type = polynomial_base<univariate_polynomial_traits>;

            using coefficients_type = typename multivariate_polynomial
                    <polynomial_base<univariate_polynomial_traits>, n_variates>::coefficients_type;

            if (p.empty())
                return static_cast<coefficients_type>(0);

            const std::array<index_t, n_variates> &extents = p.extents();

            std::vector<coefficients_type> a;
            std::span<const coefficients_type> c = p.coefficients();

            a.reserve(std::size(c) / extents.back());

            for (auto it = std::cbegin(c); it != std::cend(c); it += extents.back())
            {
                Traits::polynomial_proxy_type<univariate_polynomial_type> q
                        {std::span<const coefficients_type>{it, it + extents.back()}};

                a.emplace_back(q(Utility::get_value(n_variates - 1, x...)));
            }

            if constexpr (n_variates != 1)
            {
                std::vector<coefficients_type> b;

                b.reserve(std::size(a) / extents[n_variates - 2]);

                for (index_t d = n_variates - 2; d != negative_1; --d)
                {
                    for (auto it = std::cbegin(a); it != std::cend(a); it += extents[d])
                    {
                        Traits::polynomial_proxy_type<univariate_polynomial_type> q
                                {std::span{it, it + extents[d]}};

                        b.emplace_back(q(Utility::get_value(d, x...)));
                    }

                    std::swap(a, b);
                    b.clear();
                }
            }

            return a.front();
        }

        template<index_t n_variates, class univariate_polynomial_traits, class... variates_types>
        requires(sizeof...(variates_types) == n_variates and
                 (std::is_same_v
                    <variates_types, typename polynomial_base<univariate_polynomial_traits>::
                             coefficients_type> and ...))
        polynomial_base<univariate_polynomial_traits>::coefficients_type
        evaluate_with_basis_functions
                (const multivariate_polynomial
                        <polynomial_base<univariate_polynomial_traits>, n_variates> &p,
                const variates_types &...x)
        {
            using univariate_polynomial_type = polynomial_base<univariate_polynomial_traits>;
            using coefficients_type = typename multivariate_polynomial
                    <polynomial_base<univariate_polynomial_traits>, n_variates>::coefficients_type;

            if (p.empty())
                return static_cast<coefficients_type>(0);

            const std::array<index_t, n_variates> &extents = p.extents();

            std::array<std::vector<coefficients_type>, n_variates> basis_polynomials_evaluations;

            for (index_t i = 0; i != n_variates; ++i)
                basis_polynomials_evaluations[i] = evaluate_basis_functions
                        (Utility::get_value(i, x...), extents[i],
                         tag<univariate_polynomial_type>{});

            std::array<coefficients_type, n_variates> aux;

            aux.front() = basis_polynomials_evaluations.front().front();

            for (index_t i = 1; i != n_variates; ++i)
                aux[i] = aux[i - 1] * basis_polynomials_evaluations[i].front();

            coefficients_type ret = static_cast<coefficients_type>(0);

            std::array<index_t, n_variates> i{};
            while (i.front() != extents.front())
            {
                ret += p[i] * aux.back();

                ++i.back();

                for (index_t j = n_variates - 1; j != 0; --j)
                {
                    bool update = i[j] == extents[j];
                    if (update)
                    {
                        i[j] = 0;
                        ++i[j - 1];
                    }

                    if ((!update or j == 1) and i.front() != extents.front())
                    {
                        for (index_t k = j - 1; k != n_variates; ++k)
                        {
                            if (k == 0)
                                aux.front() = basis_polynomials_evaluations.front()[i.front()];

                            else
                                aux[k] = aux[k - 1] * basis_polynomials_evaluations[k][i[k]];
                        }

                        break;
                    }
                }
            }

            return ret;
        }

        template<index_t i_variate, index_t n_variates,
                univariate_polynomial_c polynomial_type, class... variates_types>
        auto evaluate_apparatus
                (std::span<const typename polynomial_type::coefficients_type> coefficients,
                 const std::array<index_t, n_variates> extents,
                 const variates_types &...x)
        {
            using variate_type = std::remove_cvref_t
                    <Utility::get_type<i_variate, variates_types...>>;

            using coefficients_type_next = decltype
            ((std::declval<polynomial_type>())(std::declval<variate_type>()));

            std::vector<coefficients_type_next> ret;
            index_t ext = extents[i_variate];

            ret.reserve(std::size(coefficients) / ext);

            for (auto it = std::cbegin(coefficients); it != std::cend(coefficients); it += ext)
            {
                Traits::polynomial_proxy_type<polynomial_type> p{std::span{it, it + ext}};

                ret.emplace_back(p(Utility::get_value(i_variate, x...)));
            }

            if constexpr (i_variate == 0)
                return ret.front();

            else
            {
                using next_polynomial_type = Traits::convert_polynomial_coefficients_type
                        <polynomial_type, coefficients_type_next>;

                std::span<const coefficients_type_next> span{ret.cbegin(), ret.cend()};

                return evaluate_apparatus<i_variate - 1, n_variates, next_polynomial_type>
                        (span, extents, x...);
            }
        }

        template<index_t n_variates, class univariate_polynomial_traits, class... variates_types>
        requires(sizeof...(variates_types) == n_variates and
                 !(std::is_same_v
                         <variates_types,
                          typename polynomial_base<univariate_polynomial_traits>::
                          coefficients_type> and ...))
        auto evaluate(const multivariate_polynomial
                <polynomial_base<univariate_polynomial_traits>, n_variates> &p,
                      const variates_types &...x)
        -> decltype(std::declval<typename polynomial_base<univariate_polynomial_traits>::
        coefficients_type>() * (x * ...))
        {
            return evaluate_apparatus<n_variates - 1, n_variates,
                                 polynomial_base<univariate_polynomial_traits>>
                    (p.coefficients(), p.extents(), x...);
        }

        template<index_t n_variates, class univariate_polynomial_traits, class variates_type>
        auto evaluate
                (const multivariate_polynomial
                        <polynomial_base<univariate_polynomial_traits>, n_variates> &p,
                 const std::array<variates_type, n_variates> &x)
        {
            auto lambda = [&p]<class... variates_types>(const variates_types &... x)
            {
                return evaluate(p, x...);
            };

            return std::apply(lambda, x);
        }

        template<index_t n_variates, class univariate_polynomial_traits, class variates_type>
        auto evaluate_with_basis_functions
                (const multivariate_polynomial
                        <polynomial_base<univariate_polynomial_traits>, n_variates> &p,
                 const std::array<variates_type, n_variates> &x)
        {
            auto lambda = [&p]<class... variates_types>(const variates_types &... x)
            {
                return evaluate_with_basis_functions(p, x...);
            };

            return std::apply(lambda, x);
        }

        template<index_t n_variates, class univariate_polynomial_traits, class... variates_types>
        requires(sizeof...(variates_types) == n_variates)
        auto evaluate
                (const multivariate_polynomial
                        <polynomial_base<univariate_polynomial_traits>, n_variates> &p,
                 const std::tuple<variates_types...> &x)
        {
            auto lambda = [&p](const variates_types &... x)
            {
                return evaluate(p, x...);
            };

            return std::apply(lambda, x);
        }
    }

    template<class univariate_polynomial_traits, index_t n_variates>
    multivariate_polynomial<Apparatus::polynomial_base<univariate_polynomial_traits>, n_variates>
    operator*(const multivariate_polynomial
                     <Apparatus::polynomial_base<univariate_polynomial_traits>, n_variates> &arg0,
              const multivariate_polynomial
                     <Apparatus::polynomial_base<univariate_polynomial_traits>, n_variates> &arg1)
    {
        using univariate_polynomial_type = Apparatus::polynomial_base<univariate_polynomial_traits>;
        using polynomial_type = multivariate_polynomial
                <Apparatus::polynomial_base<univariate_polynomial_traits>, n_variates>;
        using coefficients_type = polynomial_type::coefficients_type;

        if (arg0.empty() or arg1.empty())
            return {};

        auto lambda = [](index_t i, index_t j) -> std::vector<coefficients_type>
        {
            return (univariate_polynomial_type::basis_polynomial(i) *
                    univariate_polynomial_type::basis_polynomial(j)).coefficients();
        };

        Invocable::invocable_function_hash_table invocable
                {lambda,
                 Invocable::invocable<const std::vector<coefficients_type> &(index_t, index_t)>{}};

        std::array<index_t, n_variates> arg0_extents = arg0.extents();
        std::array<index_t, n_variates> arg1_extents = arg1.extents();
        std::array<index_t, n_variates> ret_extents =
                Utility::subtract_array(Utility::sum_array(arg0_extents, arg1_extents), 1);

        polynomial_type ret;
        ret.set_extents(ret_extents);

        for (const std::array<index_t, n_variates> &i: Utility::counter{arg0_extents})
            for (const std::array<index_t, n_variates> &j: Utility::counter{arg1_extents})
            {
                std::array<std::span<const coefficients_type>, n_variates> a;
                coefficients_type b = arg0[i] * arg1[j];

                for (index_t v = 0; v != n_variates; ++v)
                    a[v] = invocable(std::max(i[v], j[v]), std::min(i[v], j[v]));

                for (const std::array<index_t, n_variates> &k:
                        Utility::counter{Utility::sum_array(Utility::sum_array(i, j), 1)})
                {
                    coefficients_type c = b;

                    for (index_t v = 0; v != n_variates; ++v)
                        c *= a[v][k[v]];

                    ret[k] += c;
                }
            }

        return ret;
    }

    // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //
    // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //

    template<class univariate_polynomial_traits, index_t n_variates>
    template<class type>
    multivariate_polynomial<Apparatus::polynomial_base<univariate_polynomial_traits>, n_variates>::
    multivariate_polynomial(const type &arg0)
    requires(std::convertible_to<type, coefficients_type> and
             !Apparatus::implemented_conversion_c<type, multivariate_polynomial>) :
            coefficients_{static_cast<coefficients_type>(arg0)},
            extents_{Utility::uniform_array<index_t, n_variates>(1)},
            strides_{Utility::uniform_array<index_t, n_variates>(1)} {}

    template<class univariate_polynomial_traits, index_t n_variates>
    template<class type>
    multivariate_polynomial<Apparatus::polynomial_base<univariate_polynomial_traits>, n_variates>::
    multivariate_polynomial
            (const type &arg)requires(Apparatus::implemented_conversion_c<type, multivariate_polynomial>)
    {
        *this = conversion<type, multivariate_polynomial>::convert(arg);

        return;
    }

    template<class univariate_polynomial_traits, index_t n_variates>
    multivariate_polynomial<Apparatus::polynomial_base<univariate_polynomial_traits>, n_variates>::
    multivariate_polynomial(univariate_polynomial_type arg0)requires(n_variates == 1) :
            coefficients_{std::move(arg0).coefficients()},
            extents_{std::size(coefficients_)},
            strides_{1} {}

    template<class univariate_polynomial_traits, index_t n_variates>
    multivariate_polynomial<Apparatus::polynomial_base<univariate_polynomial_traits>, n_variates>::
    multivariate_polynomial
            (univariate_polynomial_type arg0, index_t var):
            coefficients_{std::move(arg0).coefficients()},
            extents_{Utility::uniform_array<index_t, n_variates>(1)},
            strides_{Utility::uniform_array<index_t, n_variates>(1)}
    {
        ASSERT_ASSUME(var < n_variates);
        extents_[var] = std::size(coefficients_);
        return;
    }

    template<class univariate_polynomial_traits, index_t n_variates>
    multivariate_polynomial<Apparatus::polynomial_base<univariate_polynomial_traits>, n_variates>::
    multivariate_polynomial
            (std::vector<coefficients_type> coefficients, const indices_type &extents) :
            coefficients_{std::move(coefficients)},
            extents_{extents}
    {
        strides_.back() = 1;

        for (index_t i = n_variates - 1; i != 0; --i)
            strides_[i - 1] = strides_[i] * extents_[i];

        ASSERT_ASSUME(strides_.front() * extents_.front() == std::size(coefficients_));

        return;
    }

    // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //

    template<class univariate_polynomial_traits, index_t n_variates>
    auto multivariate_polynomial<Apparatus::polynomial_base<univariate_polynomial_traits>, n_variates>::
    coefficients() &
    -> std::span<coefficients_type>
    {
        return {coefficients_};
    }

    template<class univariate_polynomial_traits, index_t n_variates>
    auto multivariate_polynomial<Apparatus::polynomial_base<univariate_polynomial_traits>, n_variates>::
    coefficients() const &
    -> std::span<const coefficients_type>
    {
        return {coefficients_};
    }

    template<class univariate_polynomial_traits, index_t n_variates>
    auto multivariate_polynomial<Apparatus::polynomial_base<univariate_polynomial_traits>, n_variates>::
    coefficients() &&
    -> std::vector<coefficients_type> &&
    {
        return std::move(coefficients_);
    }

    template<class univariate_polynomial_traits, index_t n_variates>
    auto multivariate_polynomial<Apparatus::polynomial_base<univariate_polynomial_traits>, n_variates>::
    operator[](const indices_type &i) ->
    coefficients_type &
    {
        if (!std::ranges::equal(i, extents_, [](index_t i, index_t j) { return i < j; }))
        {
            indices_type new_extents;

            std::transform(std::cbegin(extents_), std::cend(extents_),
                           std::cbegin(i), std::begin(new_extents),
                           [](index_t i, index_t j) { return std::max(i, j + 1); });

            set_extents(new_extents);
        }

        return coefficients_[offset(i)];
    }

    template<class univariate_polynomial_traits, index_t n_variates>
    auto multivariate_polynomial<Apparatus::polynomial_base<univariate_polynomial_traits>, n_variates>::
    operator[]
            (const indices_type &i) const -> const coefficients_type &
    {
        ASSERT_ASSUME(std::ranges::equal(i, extents_, [](index_t i, index_t j) { return i < j; }));

        return coefficients_[offset(i)];
    }

    template<class univariate_polynomial_traits, index_t n_variates>
    template<std::integral indices_types>
    auto multivariate_polynomial<Apparatus::polynomial_base<univariate_polynomial_traits>, n_variates>::
    operator[]
            (const std::array<indices_types, n_variates> &i) -> coefficients_type &
    {
        return (*this)[Utility::static_cast_array<index_t>(i)];
    }

    template<class univariate_polynomial_traits, index_t n_variates>
    template<std::integral indices_types>
    auto multivariate_polynomial<Apparatus::polynomial_base<univariate_polynomial_traits>, n_variates>::
    operator[]
            (const std::array<indices_types, n_variates> &i) const -> const coefficients_type &
    {
        return (*this)[Utility::static_cast_array<index_t>(i)];
    }

    template<class univariate_polynomial_traits, index_t n_variates>
    template<std::integral... indices_types>
    requires(sizeof...(indices_types) == n_variates)
    auto multivariate_polynomial<Apparatus::polynomial_base<univariate_polynomial_traits>, n_variates>::
    operator[](indices_types... i) -> coefficients_type &
    {
        return (*this)[std::array{static_cast<index_t>(i)...}];
    }

    template<class univariate_polynomial_traits, index_t n_variates>
    template<std::integral... indices_types>
    requires(sizeof...(indices_types) == n_variates)
    auto multivariate_polynomial<Apparatus::polynomial_base<univariate_polynomial_traits>, n_variates>::
    operator[](indices_types... i) const -> const coefficients_type &
    {
        return (*this)[std::array{static_cast<index_t>(i)...}];
    }

    template<class univariate_polynomial_traits, index_t n_variates>
    bool multivariate_polynomial<Apparatus::polynomial_base<univariate_polynomial_traits>, n_variates>::
    empty() const
    {
        return coefficients_.empty();
    }

    template<class univariate_polynomial_traits, index_t n_variates>
    auto multivariate_polynomial<Apparatus::polynomial_base<univariate_polynomial_traits>, n_variates>::
    clear() -> multivariate_polynomial &
    {
        coefficients_.clear();

        extents_ = Utility::uniform_array<index_t, n_variates>(0);
        strides_ = Utility::uniform_array<index_t, n_variates>(0);

        return *this;
    }

    template<class univariate_polynomial_traits, index_t n_variates>
    auto multivariate_polynomial<Apparatus::polynomial_base<univariate_polynomial_traits>, n_variates>::
    order() const -> indices_type
    {
        indices_type ret;

        auto it_ret = std::begin(ret);
        auto it_ext = std::cbegin(extents_);

        while (it_ret != std::end(ret))
            *it_ret++ = (*it_ext++) - 1;

        return ret;
    }

    template<class univariate_polynomial_traits, index_t n_variates>
    index_t multivariate_polynomial
            <Apparatus::polynomial_base<univariate_polynomial_traits>, n_variates>::order
            (index_t arg) const
    {
        ASSERT_ASSUME(arg < n_variates);
        return extents(arg) - 1;
    }

    template<class univariate_polynomial_traits, index_t n_variates>
    auto multivariate_polynomial<Apparatus::polynomial_base<univariate_polynomial_traits>, n_variates>::
    set_order() -> multivariate_polynomial &
    {
        return this->set_extents();
    }

    template<class univariate_polynomial_traits, index_t n_variates>
    auto multivariate_polynomial<Apparatus::polynomial_base<univariate_polynomial_traits>, n_variates>::
    set_order(const indices_type &order) -> multivariate_polynomial &
    {
        indices_type new_extents = order;
        std::ranges::for_each(new_extents, [](index_t &i){ ++i; });

        set_extents(new_extents);

        return *this;
    }

    template<class univariate_polynomial_traits, index_t n_variates>
    template<std::integral indices_types>
    auto multivariate_polynomial<Apparatus::polynomial_base<univariate_polynomial_traits>, n_variates>::
    set_order(const std::array<indices_types, n_variates> &order) -> multivariate_polynomial &
    {
        return set_order(Utility::static_cast_array<index_t>(order));
    }

    template<class univariate_polynomial_traits, index_t n_variates>
    template<class... indices_types>
    requires(sizeof...(indices_types) == n_variates and (std::integral<indices_types> and ...))
    auto multivariate_polynomial<Apparatus::polynomial_base<univariate_polynomial_traits>, n_variates>::
    set_order(indices_types... indices) -> multivariate_polynomial &
    {
        return this->set_order(std::array{static_cast<index_t>(indices)...});
    }

    template<class univariate_polynomial_traits, index_t n_variates>
    auto multivariate_polynomial<Apparatus::polynomial_base<univariate_polynomial_traits>, n_variates>::
    set_extents() -> multivariate_polynomial &
    {
        indices_type new_extents = extents_;

        for (index_t i = 0; i != n_variates; ++i)
            for (; new_extents[i] != 0 and !check_order(new_extents, i); --new_extents[i])
                ;

        set_extents(new_extents);

        return *this;
    }

    template<class univariate_polynomial_traits, index_t n_variates>
    auto multivariate_polynomial<Apparatus::polynomial_base<univariate_polynomial_traits>, n_variates>::
    set_extents(const indices_type &new_extents) -> multivariate_polynomial &
    {
        if (new_extents == extents_)
            return *this;

        if (std::ranges::find(new_extents,0) != std::end(new_extents))
        {
            this->clear();

            return *this;
        }

        indices_type new_strides;

        new_strides.back() = 1;

        for (index_t j = n_variates - 1; j != 0; --j)
            new_strides[j - 1] = new_strides[j] * new_extents[j];

        index_t n = new_strides.front() * new_extents.front();

        std::vector<coefficients_type> new_coefficients(n, static_cast<coefficients_type>(0));

        if (!coefficients_.empty())
            reset_coefficients(new_coefficients, new_extents, new_strides);

        coefficients_ = std::move(new_coefficients);
        extents_ = std::move(new_extents);
        strides_ = std::move(new_strides);

        return *this;
    }

    template<class univariate_polynomial_traits, index_t n_variates>
    template<std::integral indices_types>
    auto multivariate_polynomial<Apparatus::polynomial_base<univariate_polynomial_traits>, n_variates>::
    set_extents(const std::array<indices_types, n_variates> &new_extents)
    -> multivariate_polynomial &
    {
        return this->set_extents(Utility::static_cast_array<index_t>(new_extents));
    }

    template<class univariate_polynomial_traits, index_t n_variates>
    template<class... indices_types>
    requires(sizeof...(indices_types) == n_variates and (std::integral<indices_types> and ...))
    auto multivariate_polynomial<Apparatus::polynomial_base<univariate_polynomial_traits>, n_variates>::
    set_extents(indices_types... indices) -> multivariate_polynomial &
    {
        return this->set_extents(std::array{static_cast<index_t>(indices)...});
    }

    template<class univariate_polynomial_traits, index_t n_variates>
    auto multivariate_polynomial<Apparatus::polynomial_base<univariate_polynomial_traits>, n_variates>::
    extents() const -> const indices_type &
    {
        return extents_;
    }

    template<class univariate_polynomial_traits, index_t n_variates>
    index_t multivariate_polynomial
            <Apparatus::polynomial_base<univariate_polynomial_traits>, n_variates>::extents
            (index_t arg) const
    {
        ASSERT_ASSUME(arg < n_variates);
        return extents_[arg];
    }

    template<class univariate_polynomial_traits, index_t n_variates>
    auto multivariate_polynomial<Apparatus::polynomial_base<univariate_polynomial_traits>, n_variates>::
    strides() const -> const indices_type &
    {
        return strides_;
    }

    template<class univariate_polynomial_traits, index_t n_variates>
    index_t multivariate_polynomial
            <Apparatus::polynomial_base<univariate_polynomial_traits>, n_variates>::
            strides(index_t arg) const
    {
        ASSERT_ASSUME(arg < n_variates);
        return strides_[arg];
    }

    template<class univariate_polynomial_traits, index_t n_variates>
    bool multivariate_polynomial<Apparatus::polynomial_base<univariate_polynomial_traits>, n_variates>::
    check_order(const indices_type &extents, index_t i, indices_type aux, index_t j) const
    {
        if (j == n_variates)
            return (*this)[aux] != 0;

        if (j == i)
        {
            aux[j] = extents[i] - 1;
            return check_order(extents, i, aux, j + 1);
        }

        for (index_t k = 0; k != extents_[j]; ++k)
        {
            ASSERT_ASSUME(j < n_variates);
            aux[j] = k;

            if (check_order(extents, i, aux, j + 1))
                return true;
        }

        return false;
    }

    template<class univariate_polynomial_traits, index_t n_variates>
    void multivariate_polynomial<Apparatus::polynomial_base<univariate_polynomial_traits>, n_variates>::
    reset_coefficients
            (std::vector<coefficients_type> &new_coefficients,
             const indices_type &new_extents,
             const indices_type &new_strides, index_t i_old,
             index_t i_new, index_t i_var)
    {
        if (i_var == n_variates)
        {
            new_coefficients[i_new] = std::move(coefficients_[i_old]);

            return;
        }

        for (index_t i = 0; i != std::min(new_extents[i_var], extents_[i_var]); ++i)
        {
            reset_coefficients(new_coefficients, new_extents, new_strides, i_old, i_new, i_var + 1);

            i_old += strides_[i_var];
            i_new += new_strides[i_var];
        }

        return;
    }

    template<class univariate_polynomial_traits, index_t n_variates>
    index_t
    multivariate_polynomial<Apparatus::polynomial_base<univariate_polynomial_traits>, n_variates>::offset
        (const indices_type &arg) const
    {
        return std::inner_product(std::cbegin(strides_), std::cend(strides_),
                                  std::cbegin(arg), static_cast<index_t>(0));
    }

    // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //

    template<class univariate_polynomial_traits, index_t n_variates>
    template<class... variate_types>
    auto multivariate_polynomial<Apparatus::polynomial_base<univariate_polynomial_traits>, n_variates>::
    operator()(const variate_types &...x) const -> decltype(Apparatus::evaluate(*this, x...))
    {
        return Apparatus::evaluate(*this, x...);
    }

    template<class univariate_polynomial_traits, index_t n_variates>
    auto multivariate_polynomial<Apparatus::polynomial_base<univariate_polynomial_traits>, n_variates>::
    operator+=(const multivariate_polynomial &arg0) -> multivariate_polynomial &
    {
        if (empty())
            return *this = arg0;

        if (arg0.empty())
            return *this;

        elementwise_op(arg0, [](coefficients_type &i,const coefficients_type &j) { i += j; });

        return *this;
    }

    template<class univariate_polynomial_traits, index_t n_variates>
    auto multivariate_polynomial<Apparatus::polynomial_base<univariate_polynomial_traits>, n_variates>::
    operator-=(const multivariate_polynomial &arg0) -> multivariate_polynomial &
    {
        if (coefficients_.empty())
        {
            *this = -arg0;
            return *this;
        }

        if (arg0.coefficients().empty())
            return *this;

        elementwise_op(arg0, [](coefficients_type &i, const coefficients_type &j){ i -= j; });

        return *this;
    }

    template<class univariate_polynomial_traits, index_t n_variates>
    auto multivariate_polynomial<Apparatus::polynomial_base<univariate_polynomial_traits>, n_variates>::
    operator*=(const multivariate_polynomial &arg0) -> multivariate_polynomial &
    {
        return *this = *this * arg0;
    }

    template<class univariate_polynomial_traits, index_t n_variates>
    auto multivariate_polynomial<Apparatus::polynomial_base<univariate_polynomial_traits>, n_variates>::
    operator+=(const coefficients_type &arg0) -> multivariate_polynomial &
    {
        if (coefficients_.empty())
            set_extents(Utility::uniform_array<index_t, n_variates>(1));

        coefficients_[0] += arg0;

        return *this;
    }

    template<class univariate_polynomial_traits, index_t n_variates>
    auto multivariate_polynomial<Apparatus::polynomial_base<univariate_polynomial_traits>, n_variates>::
    operator-=(const coefficients_type &arg0) -> multivariate_polynomial &
    {
        if (coefficients_.empty())
            set_extents(Utility::uniform_array<index_t, n_variates>(1));

        coefficients_[0] -= arg0;

        return *this;
    }

    template<class univariate_polynomial_traits, index_t n_variates>
    auto multivariate_polynomial<Apparatus::polynomial_base<univariate_polynomial_traits>, n_variates>::
    operator*=(const coefficients_type &arg0) -> multivariate_polynomial &
    {
        std::ranges::for_each(coefficients_,
                              [&arg0](coefficients_type &arg_bis) { arg_bis *= arg0; });

        return *this;
    }

    template<class univariate_polynomial_traits, index_t n_variates>
    auto multivariate_polynomial<Apparatus::polynomial_base<univariate_polynomial_traits>, n_variates>::
    operator/=(const coefficients_type &arg0) -> multivariate_polynomial &
    {
        std::ranges::for_each(coefficients_,
                              [&arg0](coefficients_type &arg_bis) { arg_bis /= arg0; });

        return *this;
    }

    template<class univariate_polynomial_traits, index_t n_variates>
    auto multivariate_polynomial<Apparatus::polynomial_base<univariate_polynomial_traits>, n_variates>::
    operator^=(index_t exp) -> multivariate_polynomial &
    {
        return *this = *this ^ exp;
    }

    template<class univariate_polynomial_traits, index_t n_variates>
    template<class invocable_type>
    void multivariate_polynomial<Apparatus::polynomial_base<univariate_polynomial_traits>, n_variates>::
    elementwise_op(const multivariate_polynomial &arg0, const invocable_type &invocable)
    {
        indices_type new_extents;

        auto it_ex_this = std::cbegin(this->extents_);
        auto it_ex_arg0 = std::cbegin(arg0.extents_);
        auto it_ex_new = std::begin(new_extents);

        while (it_ex_new != std::end(new_extents))
            *it_ex_new++ = std::max(*it_ex_this++, *it_ex_arg0++);

        set_extents(new_extents);

        elementwise_op_apparatus(arg0, invocable);

        return;
    }

    template<class univariate_polynomial_traits, index_t n_variates>
    template<class invocable_type>
    void multivariate_polynomial<Apparatus::polynomial_base<univariate_polynomial_traits>, n_variates>::
    elementwise_op_apparatus
            (const multivariate_polynomial &arg0,
             const invocable_type &invocable,
             index_t i_var, index_t offset_this, index_t offset_arg0)
    {
        if (i_var == n_variates)
        {
            invocable(coefficients_[offset_this], arg0.coefficients_[offset_arg0]);
            return;
        }

        for (index_t i = 0; i != std::min(extents_[i_var], arg0.extents_[i_var]); ++i)
        {
            elementwise_op_apparatus(arg0, invocable, i_var + 1, offset_this, offset_arg0);

            offset_this += strides_[i_var];
            offset_arg0 += arg0.strides_[i_var];
        }

        return;
    }

    // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //

    template<class univariate_polynomial_traits, index_t n_variates>
    auto multivariate_polynomial<Apparatus::polynomial_base<univariate_polynomial_traits>, n_variates>::
    basis_polynomial(const indices_type &i) -> multivariate_polynomial
    {
        multivariate_polynomial ret;
        ret[i] = static_cast<coefficients_type>(1);

        return ret;
    }

    template<class univariate_polynomial_traits, index_t n_variates>
    template<std::integral indices_types>
    auto multivariate_polynomial<Apparatus::polynomial_base<univariate_polynomial_traits>, n_variates>::
    basis_polynomial(const std::array<indices_types, n_variates> &i) -> multivariate_polynomial
    {
        return basis_polynomial(Utility::static_cast_array<index_t>(i));
    }

    template<class univariate_polynomial_traits, index_t n_variates>
    template<class... indices_types>
    requires(sizeof...(indices_types) == n_variates and
             (std::is_convertible_v<indices_types, index_t> and ...))
    auto multivariate_polynomial<Apparatus::polynomial_base<univariate_polynomial_traits>, n_variates>::
    basis_polynomial(indices_types... i) -> multivariate_polynomial
    {
        return basis_polynomial(std::array<index_t, n_variates>{static_cast<index_t>(i)...});
    }

// **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //

    template<class univariate_polynomial_traits, index_t n_variates>
    multivariate_polynomial<Apparatus::polynomial_base<univariate_polynomial_traits>, n_variates>
    operator-(multivariate_polynomial
              <Apparatus::polynomial_base<univariate_polynomial_traits>, n_variates> arg0)
    {
        std::ranges::for_each(arg0.coefficients(),
                              [](auto &arg_bis)
                              {
                                  arg_bis = -std::move(arg_bis);
                              });

        return arg0;
    }

    template<class univariate_polynomial_traits, index_t n_variates>
    multivariate_polynomial<Apparatus::polynomial_base<univariate_polynomial_traits>, n_variates>
    operator+(multivariate_polynomial
              <Apparatus::polynomial_base<univariate_polynomial_traits>, n_variates> arg0,
              const multivariate_polynomial
                      <Apparatus::polynomial_base<univariate_polynomial_traits>, n_variates> &arg1)
    {
        return arg0 += arg1;
    }

    template<class univariate_polynomial_traits, index_t n_variates>
    multivariate_polynomial<Apparatus::polynomial_base<univariate_polynomial_traits>, n_variates>
    operator-(multivariate_polynomial
              <Apparatus::polynomial_base<univariate_polynomial_traits>, n_variates> arg0,
              const multivariate_polynomial
                      <Apparatus::polynomial_base<univariate_polynomial_traits>, n_variates> &arg1)
    {
        return arg0 -= arg1;
    }

    template<index_t n_variates, class univariate_polynomial_traits>
    multivariate_polynomial<Apparatus::polynomial_base<univariate_polynomial_traits>, n_variates>
    operator+(multivariate_polynomial
              <Apparatus::polynomial_base<univariate_polynomial_traits>, n_variates> arg0,
             const typename univariate_polynomial_traits::coefficients_type &arg1)
    {
        return arg0 += arg1;
    }

    template<index_t n_variates, class univariate_polynomial_traits>
    multivariate_polynomial<Apparatus::polynomial_base<univariate_polynomial_traits>, n_variates>
    operator-(multivariate_polynomial
              <Apparatus::polynomial_base<univariate_polynomial_traits>, n_variates> arg0,
             const typename univariate_polynomial_traits::coefficients_type &arg1)
    {
        return arg0 -= arg1;
    }

    template<index_t n_variates, class univariate_polynomial_traits>
    multivariate_polynomial<Apparatus::polynomial_base<univariate_polynomial_traits>, n_variates>
    operator*(multivariate_polynomial
                      <Apparatus::polynomial_base<univariate_polynomial_traits>, n_variates> arg0,
              const typename univariate_polynomial_traits::coefficients_type &arg1)
    {
        return arg0 *= arg1;
    }

    template<index_t n_variates, class univariate_polynomial_traits>
    multivariate_polynomial<Apparatus::polynomial_base<univariate_polynomial_traits>, n_variates>
    operator/(multivariate_polynomial
              <Apparatus::polynomial_base<univariate_polynomial_traits>, n_variates> arg0,
             const typename univariate_polynomial_traits::coefficients_type &arg1)
    {
        return arg0 /= arg1;
    }

    template<index_t n_variates, class univariate_polynomial_traits>
    multivariate_polynomial<Apparatus::polynomial_base<univariate_polynomial_traits>, n_variates>
    operator+(const typename univariate_polynomial_traits::coefficients_type &arg0,
              multivariate_polynomial
                      <Apparatus::polynomial_base<univariate_polynomial_traits>, n_variates> arg1)
    {
        return arg1 += arg0;
    }

    template<index_t n_variates, class univariate_polynomial_traits>
    multivariate_polynomial<Apparatus::polynomial_base<univariate_polynomial_traits>, n_variates>
    operator-(const typename univariate_polynomial_traits::coefficients_type &arg0,
            multivariate_polynomial
              <Apparatus::polynomial_base<univariate_polynomial_traits>, n_variates> arg1)
    {
        std::ranges::for_each(arg1.coefficients(),
                              [&arg0](auto &arg_bis)
                              {
                                  arg_bis = arg0 - std::move(arg_bis);
                              });

        return arg1;
    }

    template<index_t n_variates, class univariate_polynomial_traits>
    multivariate_polynomial<Apparatus::polynomial_base<univariate_polynomial_traits>, n_variates>
    operator*(const typename univariate_polynomial_traits::coefficients_type &arg0,
            multivariate_polynomial
              <Apparatus::polynomial_base<univariate_polynomial_traits>, n_variates> arg1)
    {
        return arg1 *= arg0;
    }

    template<class univariate_polynomial_traits, index_t n_variates>
    multivariate_polynomial<Apparatus::polynomial_base<univariate_polynomial_traits>, n_variates>
    operator^(const multivariate_polynomial
            <Apparatus::polynomial_base<univariate_polynomial_traits>, n_variates> &arg0, index_t arg1)
    {
        using ret_type = multivariate_polynomial
                <Apparatus::polynomial_base<univariate_polynomial_traits>, n_variates>;

        if (arg1 == 0)
            return ret_type{1};

        if (arg1 == 1)
            return arg0;

        ret_type ret = arg0 ^ (arg1 / 2);
        ret *= ret;

        if (arg1 % 2)
            ret *= arg0;

        return ret;
    }
}
