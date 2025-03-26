#pragma once

#include "core/core.hpp"
#include <thread>

// **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //

namespace Smolyak
{
    using Arithmetic::index_t;
    using Arithmetic::integer_t;
    using Arithmetic::real_t;
    using Arithmetic::complex_t;
    using Arithmetic::negative_1;

    // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //

    namespace Apparatus
    {
        template<class smolyak_traits_type>
        concept smolyak_traits_void_initialiser =
        requires(index_t level)
        {
            requires std::same_as<typename smolyak_traits_type::initialiser_type, void>;
            { smolyak_traits_type{level} } -> std::same_as<smolyak_traits_type>;
        };

        template<class smolyak_traits_type>
        concept smolyak_traits_non_void_initialiser =
        requires(index_t level)
        {
            requires not std::same_as<typename smolyak_traits_type::initialiser_type, void>;
            { smolyak_traits_type
                {level, std::declval<typename smolyak_traits_type::initialiser_type>()} }
            -> std::same_as<smolyak_traits_type>;
        };

        template<class smolyak_traits_type>
        concept smolyak_traits_void_evaluate =
        requires(smolyak_traits_type smolyak_traits, std::span<real_t> const_span)
        {
            requires std::same_as<typename smolyak_traits_type::secondary_domain, void>;
            { smolyak_traits.evaluate(const_span) } -> std::same_as<real_t>;
        };

        template<class smolyak_traits_type>
        concept smolyak_traits_non_void_evaluate =
        requires(smolyak_traits_type smolyak_traits,
                 std::span<const real_t> const_span,
                smolyak_traits_type::secondary_domain x)
        {
            requires not std::same_as<typename smolyak_traits_type::secondary_domain, void>;
            { smolyak_traits.evaluate(const_span, x) } -> std::same_as<real_t>;
        };
    }

    template<class smolyak_traits_type>
    concept smolyak_traits_c =
    requires(const smolyak_traits_type smolyak_traits)
    {
        typename smolyak_traits_type::primary_domain;
        typename smolyak_traits_type::secondary_domain;
        typename smolyak_traits_type::initialiser_type;

        { smolyak_traits.coordinates() }
        -> Traits::range_c<typename smolyak_traits_type::primary_domain>;
        { smolyak_traits.new_coordinates() }
        -> Traits::range_c<typename smolyak_traits_type::primary_domain>;
    }
    and (Apparatus::smolyak_traits_void_initialiser<smolyak_traits_type> or
        Apparatus::smolyak_traits_non_void_initialiser<smolyak_traits_type>)
    and (Apparatus::smolyak_traits_void_evaluate<smolyak_traits_type> or
        Apparatus::smolyak_traits_non_void_evaluate<smolyak_traits_type>);

    template<class homomorphism_type>
    concept homomorphism_c =
    requires(const homomorphism_type homomorphism,
             index_t level, index_t stride,
             std::span<const real_t> input_span,
             std::span<real_t> output_span)
    {
        { homomorphism(level, stride, input_span, output_span) } -> std::same_as<void>;
        { homomorphism.input_size(level) } -> std::same_as<index_t>;
        { homomorphism.output_size(level) } -> std::same_as<index_t>;
    };

    namespace Apparatus
    {
        template<class smolyak_traits_type>
        concept smolyak_traits_non_trivial_linear_operator_c =
                smolyak_traits_c<smolyak_traits_type> and
        requires(const smolyak_traits_type smolyak_traits, index_t stride,
                 std::span<const real_t> const_span,
                 std::span<real_t> span)
        {
            { smolyak_traits.linear_operator_output_size() } -> std::same_as<index_t>;
            { smolyak_traits.linear_operator(stride, const_span, span) } -> std::same_as<void>;
        };

        template<class smolyak_traits_type>
        concept smolyak_traits_trivial_linear_operator_c =
                smolyak_traits_c<smolyak_traits_type> and
                not smolyak_traits_non_trivial_linear_operator_c<smolyak_traits_type>;

        template<class smolyak_traits_type>
        concept smolyak_traits_inner_product_c =
                smolyak_traits_c<smolyak_traits_type> and
        requires(smolyak_traits_type smolyak_traits, index_t level,
                std::span<const real_t> const_span)
        {
            requires homomorphism_c<typename smolyak_traits_type::embedding_type>;

            requires homomorphism_c<typename smolyak_traits_type::dual_embedding_type>;

            { smolyak_traits.embedding(level) }
            -> std::same_as<typename smolyak_traits_type::embedding_type>;
            { smolyak_traits.dual_embedding(level) }
            -> std::same_as<typename smolyak_traits_type::dual_embedding_type>;
        };
    }
}

// **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //

