#pragma once

#include "fwd.hpp"

// **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //

namespace Smolyak::Smolyak_traits
{
    template<bool incremental_ = true>
    class multifidelity
    {
    public:
        using primary_domain = index_t;
        using secondary_domain = void;
        using initialiser_type = void;
        static constexpr bool incremental = incremental_;

        class embedding_type;
        class dual_embedding_type;

        // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //
    private:
        const std::array<index_t, 2> fidelity_level_;

        // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //
    public:
        multifidelity() = delete;
        ~multifidelity() = default;

        multifidelity(const multifidelity &) = default;
        multifidelity &operator=(const multifidelity &) = default;

        multifidelity(multifidelity &&) noexcept = default;
        multifidelity &operator=(multifidelity &&) noexcept = default;
        multifidelity(index_t fidelity_level) :
            fidelity_level_{fidelity_level - 1, fidelity_level} {}

        // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //
    public:
        void linear_operator(index_t, std::span<const real_t>, std::span<real_t>) const
            requires(incremental);
        index_t linear_operator_output_size() const requires(incremental) { return 1; }
        real_t evaluate(std::span<const real_t>) const;
        std::span<const primary_domain> coordinates() const;
        std::span<const primary_domain> new_coordinates() const
            { return {std::cbegin(fidelity_level_) + 1, std::cend(fidelity_level_)}; }
        embedding_type embedding(index_t) const { return {}; }
        dual_embedding_type dual_embedding(index_t) const { return {}; }
    };

    // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //
    // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //

    template<bool incremental>
    class multifidelity<incremental>::embedding_type
    {
    public:
        embedding_type() = default;
        ~embedding_type() = default;

        embedding_type(const embedding_type &) = default;
        embedding_type &operator=(const embedding_type &) = default;

        embedding_type(embedding_type &&) noexcept = default;
        embedding_type &operator=(embedding_type &&) noexcept = default;

        // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //
    public:
        void operator()(index_t, index_t, std::span<const real_t>, std::span<real_t>) const;
        index_t input_size(index_t) const { return 1; }
        index_t output_size(index_t) const { return 1; }

    };

    // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //
    // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //

    template<bool incremental>
    class multifidelity<incremental>::dual_embedding_type
    {
    public:
        dual_embedding_type() = default;
        ~dual_embedding_type() = default;

        dual_embedding_type(const dual_embedding_type &) = default;
        dual_embedding_type &operator=(const dual_embedding_type &) = default;

        dual_embedding_type(dual_embedding_type &&) noexcept = default;
        dual_embedding_type &operator=(dual_embedding_type &&) noexcept = default;

        // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //
    public:
        void operator()(index_t, index_t, std::span<const real_t>, std::span<real_t>) const;
        index_t input_size(index_t) const { return 1; }
        index_t output_size(index_t) const { return 1; }
    };

    // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //
    // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //

    template<bool incremental>
    void multifidelity<incremental>::linear_operator
        (index_t stride, std::span<const real_t> input_span, std::span<real_t> output_span)
        const requires(incremental)
    {
        ASSERT_ASSUME(std::size(input_span) == (fidelity_level_.back() == 0 ? 1 : stride + 1));
        ASSERT_ASSUME(std::size(output_span) == 1);

        if (fidelity_level_.back() == 0)
            output_span.front() = input_span.front();
        else
            output_span.front() = input_span.back() - input_span.front();
        return;
    }

    template<bool incremental>
    real_t multifidelity<incremental>::evaluate(std::span<const real_t> input_span) const
    {
        ASSERT_ASSUME(std::size(input_span) == 1);
        return input_span.front();
    }

    template<bool incremental>
    auto multifidelity<incremental>::coordinates() const -> std::span<const primary_domain>
    {
        if constexpr (incremental)
        {
            if (fidelity_level_.back() == 0)
                return {std::cbegin(fidelity_level_) + 1, std::cend(fidelity_level_)};
            else
                return {std::cbegin(fidelity_level_), std::cend(fidelity_level_)};
        }
        else
            return {std::cbegin(fidelity_level_) + 1, std::cend(fidelity_level_)};
    }

    // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //
    // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //

    template<bool incremental>
    void multifidelity<incremental>::embedding_type::operator()
        (index_t, index_t, std::span<const real_t> input_span, std::span<real_t> output_span) const
    {
        ASSERT_ASSUME(std::size(input_span) == 1 and std::size(output_span) == 1);

        output_span.front() = input_span.front();
        return;
    }

    // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //
    // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //

    template<bool incremental>
    void multifidelity<incremental>::dual_embedding_type::operator()
        (index_t, index_t, std::span<const real_t> input_span, std::span<real_t> output_span) const
    {
        ASSERT_ASSUME(std::size(input_span) == 1 and std::size(output_span) == 1);

        output_span.front() = input_span.front();
        return;
    }
}

// **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //

