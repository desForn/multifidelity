#pragma once

#include "fwd.hpp"

// **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //

namespace Invocable
{
    template<class function_type, class output_type, class... input_types>
    class invocable_function_binary_tree<output_type(input_types...), function_type>
    {
        static_assert(sizeof...(input_types) != 0);

        // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //
    public:
        using signature_type = output_type(input_types...);

        using key_type = std::tuple<std::remove_cvref_t<input_types>...>;
        using table_type = std::map
                <key_type, std::remove_cvref_t<output_type>, std::less<key_type>>;

        // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //
    private:
        const function_type function_;
        mutable table_type table_{};

        // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //
    public:
        invocable_function_binary_tree() = delete;
        ~invocable_function_binary_tree() = default;

        invocable_function_binary_tree(const invocable_function_binary_tree &) = delete;
        invocable_function_binary_tree &operator=(const invocable_function_binary_tree &) = delete;

        invocable_function_binary_tree(invocable_function_binary_tree &&) noexcept = default;
        invocable_function_binary_tree &
        operator=(invocable_function_binary_tree &&) noexcept = default;

        explicit invocable_function_binary_tree(function_type function) :
                function_(std::move(function)) {}

        explicit invocable_function_binary_tree(function_type function, invocable<signature_type>) :
                function_(std::move(function)) {}

        // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //
    public:
        template<class... input_types_bis>
        requires(is_callable_c<signature_type, input_types_bis...>)
        output_type operator()(input_types_bis &&... inputs) const
        {
            std::tuple inputs_tuple = std::make_tuple(std::forward<input_types_bis>(inputs)...);
            auto it = table_.find(inputs_tuple);

            if (it != std::end(table_))
                return std::get<1>(*it);

            output_type ret = std::apply(function_, inputs_tuple);
            table_.emplace(std::move(std::make_pair(inputs_tuple, ret)));

            return ret;
        }

        // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //
    public:
        const function_type &function() const & { return function_; }
        function_type &&function() && { return std::move(function_); }

        const table_type &table() const & { return table_; }
        table_type &&table() && { return std::move(table_); }

        invocable_function_binary_tree &clear() { table_.clear(); return *this; }
    };

    // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //
    // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //

    template<class function_type>
    invocable_function_binary_tree(function_type) ->
    invocable_function_binary_tree<invocable_signature<function_type>, function_type>;

    template<class function_type, class signature_type>
    invocable_function_binary_tree(function_type, invocable<signature_type>) ->
    invocable_function_binary_tree<signature_type, function_type>;
}
