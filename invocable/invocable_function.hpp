#pragma once

#include "fwd.hpp"

// **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //

namespace Invocable
{
    template<class function_type, class output_type, class... input_types>
    class invocable_function<output_type(input_types...), function_type>
    {
        // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //
    public:
        using signature_type = output_type(input_types...);

        // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //
    protected:
        const function_type function_;

        // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //
    public:
        invocable_function() = delete;
        ~invocable_function() = default;

        invocable_function(const invocable_function &) = default;
        invocable_function &operator=(const invocable_function &) = default;

        invocable_function(invocable_function &&) noexcept = default;
        invocable_function &operator=(invocable_function &&) noexcept = default;

        explicit invocable_function(function_type function) :
                function_{std::move(function)} {}

        invocable_function(function_type function, const invocable<signature_type> &) :
                function_{std::move(function)} {}

        // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //
    public:
        template<class... call_operator_input_types>
        requires(is_callable_c<signature_type, call_operator_input_types...>)
        output_type operator()(call_operator_input_types &&... inputs) const
        {
            return function_(std::forward<call_operator_input_types>(inputs)...);
        }

        // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //
    public:
        const function_type &function() const & { return function_; }
        function_type &&function() && { return std::move(function_); }
    };

    // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //

    template<class function_type>
    invocable_function(function_type)
        -> invocable_function<invocable_signature<function_type>, function_type>;

    template<class function_type, class function_signature>
    invocable_function(function_type, invocable<function_signature>)
        -> invocable_function<function_signature, function_type>;
}
