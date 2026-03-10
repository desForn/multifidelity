#pragma once

#include "fwd.hpp"
#include "polynomial/polynomial_point_value.hpp"
#include "polynomial/quadrature_weights.hpp"

// **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //

namespace Smolyak::Smolyak_traits
{
    template<Sampling::sampling_c sampling_type, bool incremental_ = true>
    class polynomial_interpolation
    {
    public:
        using polynomial_type = Polynomial::polynomial_point_value<sampling_type>;
        using proxy_polynomial_type = Polynomial::polynomial_point_value<sampling_type, real_t,
              std::ranges::stride_view<std::span<const real_t>>>;

        using field_type = polynomial_type::codomain_field;
        using primary_domain = polynomial_type::domain_field;
        using secondary_domain = primary_domain;
        using initialiser_type = void;
        static constexpr bool incremental = incremental_;

        class embedding_type;
        class dual_embedding_type;
        template<class function_type>
        class integrate_homomorphism_type;

    private:
        static constexpr auto one_function = [](real_t) -> real_t { return 1; };

        // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //
    private:
        const index_t sampling_level_;

        // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //
    public:
        polynomial_interpolation() = delete;
        ~polynomial_interpolation() = default;

        polynomial_interpolation(const polynomial_interpolation &) = default;
        polynomial_interpolation &operator=(const polynomial_interpolation &) = default;

        polynomial_interpolation(polynomial_interpolation &&) noexcept = default;
        polynomial_interpolation &operator=(polynomial_interpolation &&) noexcept = default;

        polynomial_interpolation(index_t arg) :
            sampling_level_{arg} {}

        // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //
    public:
        void linear_operator(index_t, std::span<const field_type>, std::span<field_type>) const
            requires(incremental);
        index_t linear_operator_output_size() const requires(incremental)
            { return sampling_type::n_points(sampling_level_); }
        field_type evaluate
            (std::span<const field_type> input_span, const secondary_domain &x) const
        {
            std::ranges::stride_view<std::span<const field_type>> a =
                std::ranges::stride_view{input_span, 1};
            return proxy_polynomial_type{a, sampling_level_}(x);
        }
        std::span<const primary_domain> coordinates() const
            { return sampling_type::factory(sampling_level_).coordinates(); }
        std::span<const primary_domain> new_coordinates() const
            { return sampling_type::factory(sampling_level_).new_coordinates(); }
        embedding_type embedding(index_t destination_level) const
            { return embedding_type{destination_level}; }
        dual_embedding_type dual_embedding(index_t destination_level) const
            { return dual_embedding_type{destination_level};}
        integrate_homomorphism_type<decltype(one_function)> integrate_homomorphism() const
            { return integrate_homomorphism_type<decltype(one_function)>{one_function, 0}; }
        integrate_homomorphism_type<decltype(one_function)> integrate_homomorphism
            (index_t increased_accuracy) const
            { return integrate_homomorphism_type<decltype(one_function)>
                {one_function, increased_accuracy}; }
        template<class invocable_type>
        requires(Invocable::invocable_c<invocable_type, real_t(real_t)>)
        integrate_homomorphism_type<invocable_type> integrate_homomorphism
            (invocable_type &&invocable) const
            { return integrate_homomorphism_type{std::forward<invocable_type>(invocable)}; }
        template<class invocable_type>
        requires(Invocable::invocable_c<invocable_type, real_t(real_t)>)
        integrate_homomorphism_type<invocable_type> integrate_homomorphism
        (invocable_type &&invocable, index_t increased_accuracy) const
        {
            return integrate_homomorphism_type
                {std::forward<invocable_type>(invocable), increased_accuracy};
        }
    };

    // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //
    // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //

    template<Sampling::sampling_c sampling_type, bool incremental>
    class polynomial_interpolation<sampling_type, incremental>::embedding_type
    {
    private:
        const index_t destination_level_;
        const index_t output_size_;
        const sampling_type &sampling_;

        // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //
    public:
        embedding_type() = delete;
        ~embedding_type() = default;

        embedding_type(const embedding_type &) = default;
        embedding_type &operator=(const embedding_type &) = default;

        embedding_type(embedding_type &&) noexcept = default;
        embedding_type &operator=(embedding_type &&) = default;

        explicit embedding_type(index_t);

        // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //
    public:
        void operator()
            (index_t, index_t, std::span<const field_type>, std::span<field_type>) const;
        index_t input_size(index_t level) const { return sampling_.n_points(level); }
        index_t output_size(index_t) const { return output_size_; }
    };

    // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //
    // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //

    template<Sampling::sampling_c sampling_type, bool incremental>
    class polynomial_interpolation<sampling_type, incremental>::dual_embedding_type
    {
    private:
        using quadrature_type = Polynomial::quadrature_weights<sampling_type>;

        // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //
    private:
        const index_t destination_level_;
        const index_t output_size_;
        const sampling_type &sampling_;
        const quadrature_type &quadrature_;

        // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //
    public:
        dual_embedding_type() = delete;
        ~dual_embedding_type() = default;

        dual_embedding_type(const dual_embedding_type &) = default;
        dual_embedding_type &operator=(const dual_embedding_type &) = default;

        dual_embedding_type(dual_embedding_type &&) noexcept = default;
        dual_embedding_type &operator=(dual_embedding_type &&) = default;

        explicit dual_embedding_type(index_t);

        // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //
    public:
        void operator()
            (index_t, index_t, std::span<const field_type>, std::span<field_type>) const;
        index_t input_size(index_t level) const { return sampling_type::n_points(level); }
        index_t output_size(index_t) const { return output_size_; }
    };

    // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //
    // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //

    template<Sampling::sampling_c sampling_type, bool incremental>
    template<class function_t>
    class polynomial_interpolation<sampling_type, incremental>::integrate_homomorphism_type
    {
        static_assert(Invocable::invocable_c<function_t, real_t(real_t)>);
    public:
        using quadrature_type = Polynomial::quadrature_weights<sampling_type>;
        using function_type = function_t;

    private:
        const function_type function_;
        mutable std::unordered_map<index_t, polynomial_type> interpolated_function_{};
        index_t increased_accuracy_;

    public:
        integrate_homomorphism_type() = delete;
        ~integrate_homomorphism_type() = default;

        integrate_homomorphism_type(const integrate_homomorphism_type &) = default;
        integrate_homomorphism_type &operator=(const integrate_homomorphism_type &) = default;

        integrate_homomorphism_type(integrate_homomorphism_type &&) noexcept = default;
        integrate_homomorphism_type &operator=(integrate_homomorphism_type &&) noexcept = default;

        integrate_homomorphism_type(function_type function, index_t increased_accuracy = 0) :
            function_{std::move(function)}, increased_accuracy_{increased_accuracy} {}

        // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //
    public:
        void operator()
            (index_t, index_t, std::span<const field_type>, std::span<field_type>) const;
        index_t input_size(index_t level) const { return sampling_type::n_points(level); }
        index_t output_size(index_t) const { return 1; }
    };

    // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //
    // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //

    template<Sampling::sampling_c sampling_type, bool incremental>
    void polynomial_interpolation<sampling_type, incremental>::linear_operator
        (index_t stride, std::span<const field_type> input_span, std::span<field_type> output_span)
        const requires(incremental)
    {
        ASSERT_ASSUME(std::size(input_span) == stride * (linear_operator_output_size() - 1) + 1);
        ASSERT_ASSUME(std::size(input_span) == std::size(output_span));

        if (sampling_level_ == 0)
        {
            auto it_i = std::cbegin(input_span);
            auto it_o = std::begin(output_span);

            for (;it_i < std::cend(input_span); it_i += stride, it_o += stride)
                *it_o = *it_i;
            
            return;
        }

        if constexpr (Sampling::exponential_sampling_c<sampling_type>)
        {
            auto current_f = std::ranges::stride_view{input_span, static_cast<long int>(stride)};
            auto current_p = proxy_polynomial_type{current_f, sampling_level_};

            if (sampling_level_ == 1)
            {
                field_type previous_value = current_f[1];
                auto it_o = std::cbegin(output_span);
                auto it_c = std::cbegin(coordinates());

                for(; it_o < std::cend(output_span); it_o += stride, ++it_c)
                    *it_o = current_p(*it_c) - previous_value;
            }

            else
            {
                auto previous_f = std::ranges::stride_view{input_span,
                    2 * static_cast<long int>(stride)};
                auto previous_p = proxy_polynomial_type{previous_f, sampling_level_ - 1};

                auto it_o = std::cbegin(output_span);
                auto it_c = std::cbegin(coordinates());

                for(; it_o < std::cend(output_span); it_o += stride, ++it_c)
                    *it_o = current_p(*it_c) - previous_p(*it_c);
            }

            return;
        }

        else
        {
            auto current_f = std::ranges::stride_view{input_span, stride};
            auto current_p = proxy_polynomial_type{current_f, sampling_level_};
            auto previous_p = polynomial_type{current_p, sampling_level_ - 1};
            
            auto it_o = std::cbegin(output_span);
            auto it_c = std::cbegin(coordinates());

            for(; it_o < std::cend(output_span); it_o += stride, ++it_c)
                *it_o = current_p(*it_c) - previous_p(*it_c);

            return;
        }
    }

    // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //
    // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //

    template<Sampling::sampling_c sampling_type, bool incremental>
    polynomial_interpolation<sampling_type, incremental>::embedding_type::embedding_type
            (index_t destination_level) :
            destination_level_
                    {sampling_type::required_level(
                            2 * sampling_type::n_points(destination_level) - 1)},
            output_size_{sampling_type::n_points(destination_level_)},
            sampling_{sampling_type::factory(destination_level_)} {}

    // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //

    template<Sampling::sampling_c sampling_type, bool incremental>
    void polynomial_interpolation<sampling_type, incremental>::embedding_type::operator()
        (index_t level, index_t stride, std::span<const field_type> input_span,
        std::span<field_type> output_span) const
    {
        ASSERT_ASSUME(level <= destination_level_);
        ASSERT_ASSUME(std::size(input_span) == (input_size(level) - 1) * stride + 1);
        ASSERT_ASSUME(std::size(output_span) == (output_size(level) - 1) * stride + 1);

        auto a = std::ranges::stride_view{input_span, static_cast<long int>(stride)};
        proxy_polynomial_type p{a, level};

        auto it = std::begin(output_span);
        for (const auto &i: sampling_)
        {
            *it = p(i);
            it += stride;
        }

        return;
    }

    // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //
    // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //

    template<Sampling::sampling_c sampling_type, bool incremental>
    polynomial_interpolation<sampling_type, incremental>::dual_embedding_type::dual_embedding_type
            (index_t destination_level) :
            destination_level_
                    {sampling_type::required_level(
                            2 * sampling_type::n_points(destination_level) - 1)},
            output_size_{sampling_type::n_points(destination_level_)},
            sampling_{sampling_type::factory(destination_level_)},
            quadrature_{quadrature_type::factory(destination_level_)} {}

    // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //

    template<Sampling::sampling_c sampling_type, bool incremental>
    void polynomial_interpolation<sampling_type, incremental>::dual_embedding_type::operator()
        (index_t level, index_t stride, std::span<const field_type> input_span,
        std::span<field_type> output_span) const
    {
        ASSERT_ASSUME(level <= destination_level_);
        ASSERT_ASSUME(std::size(input_span) == (input_size(level) - 1) * stride + 1);
        ASSERT_ASSUME(std::size(output_span) == (output_size(level) - 1) * stride + 1);

        auto a = std::ranges::stride_view{input_span, static_cast<long int>(stride)};
        proxy_polynomial_type p{a, level};

        auto it = std::begin(output_span);
        index_t c = 0;

        for (const auto &i: sampling_)
        {
            if constexpr (std::same_as<field_type, complex_t>)
                *it = std::conj(p(i)) * quadrature_[c++];
            else
                *it = p(i) * quadrature_[c++];
            it += stride;
        }

        return;
    }

    template<Sampling::sampling_c sampling_type, bool incremental>
    template<class function_t>
    void polynomial_interpolation<sampling_type, incremental>::
        integrate_homomorphism_type<function_t>::operator()
        (index_t level, index_t stride, std::span<const field_type> input_span,
        std::span<field_type> output_span) const
    {
        ASSERT_ASSUME(std::size(input_span) == (input_size(level) - 1) * stride + 1);
        ASSERT_ASSUME(std::size(output_span) == 1);

        auto interpolated_function = interpolated_function_.find(level + increased_accuracy_);

        if (interpolated_function == std::cend(interpolated_function_))
        {
            auto emplaced = interpolated_function_.emplace(
                level + increased_accuracy_,
                polynomial_type{function_, level + increased_accuracy_});
            ASSERT_ASSUME(std::get<1>(emplaced));
            interpolated_function = std::get<0>(emplaced);
        }

        auto a = std::ranges::stride_view{input_span, static_cast<long int>(stride)};
        proxy_polynomial_type q{a, level};

        const quadrature_type &quadrature = quadrature_type::factory(level + increased_accuracy_);

        auto it = std::begin(output_span);

        const auto &f = std::get<1>(*interpolated_function).f();
        if (increased_accuracy_ == 0)
        {
            const auto &qf = q.f();
            ASSERT_ASSUME(std::size(f) == std::size(qf) and std::size(f) == input_size(level));

            *it = 0;
            for (index_t i = 0; i != std::size(f); ++i)
            {
                if constexpr (std::same_as<field_type, complex_t>)
                    *it += f[i] * quadrature[i] * std::conj(qf[i]);
                else
                    *it += f[i] * quadrature[i] * qf[i];
            }
        }

        else
        {
            sampling_type sampling = sampling_type::factory(level + increased_accuracy_);
            ASSERT_ASSUME(sampling.n_points() == std::size(f));
            *it = 0;
            index_t c = 0;

            for (const auto &i : sampling)
            {
                if constexpr (std::same_as<field_type, complex_t>)
                    *it += f[c] * quadrature[c] * std::conj(q(i));
                else
                    *it += f[c] * quadrature[c] * q(i);
                    
                ++c;
            }
        }

        return;
    }

    // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //
    // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //

}

// **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //
// **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //

