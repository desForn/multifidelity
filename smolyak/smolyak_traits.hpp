#pragma once

#include "smolyak_traits_polynomial_interpolation.hpp"
#include "smolyak_traits_multifidelity.hpp"
#include "smolyak_traits_multifidelity_richardson_extrapolation.hpp"
#include "smolyak_traits_matrix.hpp"

// **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //

namespace Smolyak::Smolyak_traits
{
    template<class function_t>
    class identity_homomorphism_type
    {
    public:
        using function_type = function_t;

    private:
        const function_type function_;

    public:
        identity_homomorphism_type() = delete;
        ~identity_homomorphism_type() = default;

        identity_homomorphism_type(const identity_homomorphism_type &) = default;
        identity_homomorphism_type &operator=(const identity_homomorphism_type &) = default;

        identity_homomorphism_type(identity_homomorphism_type &&) noexcept = default;
        identity_homomorphism_type &operator=(identity_homomorphism_type &&) noexcept = default;

        identity_homomorphism_type(function_type &&function) : function_{std::move(function)} {}

    // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //
    public:
        void operator()(index_t level, index_t stride, std::span<const real_t> input_span,
                std::span<real_t> output_span) const
        {
            ASSERT_ASSUME(std::size(input_span) == (input_size(level) - 1) * stride + 1);
            ASSERT_ASSUME(std::size(output_span) == (output_size(level) - 1) * stride + 1);

            auto i = std::cbegin(input_span);
            auto j = std::begin(output_span);

            while (i < std::cend(input_span))
            {
                *j = *i;
                i += stride;
                j += stride;
            }

            return;
        }

        index_t input_size(index_t level) const { return function_(level); }
        index_t output_size(index_t level) const { return function_(level); }
    };
}

// **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //
// **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //

