#pragma once

#include "fwd.hpp"

// **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //

namespace Invocable
{
    template<class function_type, class output_type, class input_type>
    class invocable_function_sequential_vector<output_type(input_type), function_type>
    {
        static_assert(std::is_integral_v<std::remove_cvref_t<input_type>>);

        // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //
    public:
        using signature_type = output_type(input_type);
        using table_type = std::vector<std::remove_cvref_t<output_type>>;

        // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //
    private:
        mutable table_type table_{};
        function_type function_;

        // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //
    public:
        invocable_function_sequential_vector() = delete;
        ~invocable_function_sequential_vector() = default;

        invocable_function_sequential_vector
                (const invocable_function_sequential_vector &) = delete;
        invocable_function_sequential_vector &operator=
                (const invocable_function_sequential_vector &) = delete;

        invocable_function_sequential_vector
                (invocable_function_sequential_vector &&) noexcept = default;
        invocable_function_sequential_vector &operator=
                (invocable_function_sequential_vector &&) noexcept = default;

        explicit invocable_function_sequential_vector(function_type function) :
                function_(std::move(function)) {}

        invocable_function_sequential_vector
                (function_type function, const invocable<signature_type>) :
                function_(std::move(function)) {}

        // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //
    public:
        output_type operator()(input_type input) const
        {
            if constexpr (std::signed_integral<input_type>)
                ASSERT_ASSUME(input >= 0);

            if (static_cast<std::size_t>(input) >= std::size(table_))
            {
                table_.reserve(input + 1);

                for (index_t i = std::size(table_); i <= static_cast<index_t>(input); ++i)
                    table_.emplace_back(function_(i));
            }

            return table_[input];
        }

        // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //
    public:
        const function_type &function() const & { return function_; }
        function_type &&function() && { return std::move(function_); }

        const table_type &table() const & { return table_; }
        table_type &&table() && { return std::move(table_); }

        invocable_function_sequential_vector &clear() { table_.clear(); return *this; }
    };

    // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //

    template<class function_type>
    invocable_function_sequential_vector(function_type) ->
    invocable_function_sequential_vector<invocable_signature<function_type>, function_type>;

    template<class function_type, class signature_type>
    invocable_function_sequential_vector(function_type, invocable<signature_type>) ->
    invocable_function_sequential_vector<signature_type, function_type>;
}
