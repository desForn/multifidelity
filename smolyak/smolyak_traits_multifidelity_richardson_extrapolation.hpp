#pragma once

#include "fwd.hpp"
#include "polynomial/polynomial_point_value.hpp"
#include "sampling/exponential_refinement.hpp"

// **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //

namespace Smolyak::Smolyak_traits
{
    template<real_t growth_ratio_, bool incremental_ = true>
    class multifidelity_richardson_extrapolation
    {
    public:
        static constexpr real_t growth_ratio = growth_ratio_;
        using sampling_type = Sampling::exponential_refinement<growth_ratio, 1>;
        using polynomial_type = Polynomial::polynomial_point_value<sampling_type>;
        
        using primary_domain = index_t;
        using secondary_domain = void;
        using initialiser_type = void;
        static constexpr bool incremental = incremental_;

        class embedding_type;
        class dual_embedding_type;

        // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //
    private:
        std::vector<index_t> coordinates_;

        // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //
    public:
        multifidelity_richardson_extrapolation() = delete;
        ~multifidelity_richardson_extrapolation() = default;

        multifidelity_richardson_extrapolation
            (const multifidelity_richardson_extrapolation &) = default;
        multifidelity_richardson_extrapolation &operator=
            (const multifidelity_richardson_extrapolation &) = default;

        multifidelity_richardson_extrapolation
            (multifidelity_richardson_extrapolation &&) noexcept = default;
        multifidelity_richardson_extrapolation &operator=
            (multifidelity_richardson_extrapolation &&) noexcept = default;

        multifidelity_richardson_extrapolation(index_t fidelity_level) :
            coordinates_(fidelity_level + 1) { std::ranges::iota(coordinates_, 0); }

        // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //
    public:
        index_t fidelity_level() const { return std::size(coordinates_) - 1; }
        void linear_operator(index_t, std::span<const real_t>, std::span<real_t>) const;
        index_t linear_operator_output_size() const { return 1; }
        real_t evaluate(std::span<const real_t>) const;
        std::span<const primary_domain> coordinates() const { return {coordinates_}; }
        std::span<const primary_domain> new_coordinates() const
            { return {std::end(coordinates_) - 1, 1}; }
        embedding_type embedding(index_t destination_level) const
            { return embedding_type{destination_level}; }
        dual_embedding_type dual_embedding(index_t destination_level) const
            { return dual_embedding_type{destination_level}; }
    };

    // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //
    // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //

    template<real_t growth_ratio, bool increment>
    class multifidelity_richardson_extrapolation<growth_ratio, increment>::embedding_type
    {
    public:
        embedding_type() = delete;
        ~embedding_type() = default;

        embedding_type(const embedding_type &) = default;
        embedding_type &operator=(const embedding_type &) = default;

        embedding_type(embedding_type &&) noexcept = default;
        embedding_type &operator=(embedding_type &&) noexcept = default;

        explicit embedding_type(index_t) {}

        // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //
    public:
        void operator()(index_t, index_t, std::span<const real_t>, std::span<real_t>) const;
        index_t input_size(index_t) const { return 1; }
        index_t output_size(index_t) const { return 1; }

    };

    // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //
    // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //

    template<real_t growth_ratio, bool increment>
    class multifidelity_richardson_extrapolation<growth_ratio, increment>::dual_embedding_type
    {
    public:
        dual_embedding_type() = delete;
        ~dual_embedding_type() = default;

        dual_embedding_type(const dual_embedding_type &) = default;
        dual_embedding_type &operator=(const dual_embedding_type &) = default;

        dual_embedding_type(dual_embedding_type &&) noexcept = default;
        dual_embedding_type &operator=(dual_embedding_type &&) noexcept = default;

        explicit dual_embedding_type(index_t) {}

        // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //
    public:
        void operator()(index_t, index_t, std::span<const real_t>, std::span<real_t>) const;
        index_t input_size(index_t) const { return 1; }
        index_t output_size(index_t) const { return 1; }
    };

    // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //
    // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //

    template<real_t growth_ratio, bool increment>
    void multifidelity_richardson_extrapolation<growth_ratio, increment>::
        linear_operator(index_t stride, std::span<const real_t> input_span,
                        std::span<real_t> output_span) const
    {
        index_t fidelity_level_ = fidelity_level();

        ASSERT_ASSUME(std::size(input_span) == stride * fidelity_level_ + 1);
        ASSERT_ASSUME(std::size(output_span) == 1);

        if (fidelity_level_ == 0)
        {
            output_span.front() = input_span.front();
            return;
        }

        if constexpr (incremental)
        {
            std::vector<real_t> a;
            for (auto it = std::cbegin(input_span); it < std::cend(input_span) - 1; it += stride)
                a.emplace_back(*it);
            real_t previous = polynomial_type{a, fidelity_level_ - 1}
                (static_cast<real_t>(0));

            a.emplace_back(input_span.back());
            real_t current = polynomial_type{a, fidelity_level_}(static_cast<real_t>(0));

            output_span.front() = current - previous;
        }
        else
        {
            std::vector<real_t> a;
            for (auto it = std::cbegin(input_span); it < std::cend(input_span); it += stride)
                a.emplace_back(*it);

            output_span.front() = polynomial_type{a, fidelity_level_}(static_cast<real_t>(0));
        }

        return;
    }

    template<real_t growth_ratio, bool increment>
    real_t multifidelity_richardson_extrapolation<growth_ratio, increment>::
        evaluate(std::span<const real_t> input_span) const
    {
        ASSERT_ASSUME(std::size(input_span) == 1);
        return input_span.front();
    }

    // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //
    // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //

    template<real_t growth_ratio, bool increment>
    void multifidelity_richardson_extrapolation<growth_ratio, increment>::embedding_type::
        operator()(index_t, index_t, std::span<const real_t> input_span,
            std::span<real_t> output_span) const
    {
        ASSERT_ASSUME(std::size(input_span) == 1 and std::size(output_span) == 1);
        output_span.front() = input_span.front();
        return;
    }

    // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //
    // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //

    template<real_t growth_ratio, bool increment>
    void multifidelity_richardson_extrapolation<growth_ratio, increment>::dual_embedding_type::
        operator()(index_t, index_t, std::span<const real_t> input_span,
            std::span<real_t> output_span) const
    {
        ASSERT_ASSUME(std::size(input_span) == 1 and std::size(output_span) == 1);
        output_span.front() = input_span.front();
        return;
    }
}

// **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //

