#pragma once

#include "fwd.hpp"

// **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //

namespace Invocable
{
    template<class function_type, class output_type, class hash_type, class... input_types>
    class invocable_function_hash_table<output_type(input_types...), function_type, hash_type>
    {
        static_assert(sizeof...(input_types) != 0);

        // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //
    public:
        using signature_type = output_type(input_types...);

    private:
        static constexpr index_t inputs_size = sizeof...(input_types);
        using first_input_type = Utility::get_type<0, input_types...>;

    public:
        using key_type = std::tuple<std::remove_cvref_t<input_types>...>;
        using table_type = std::unordered_map
                <key_type, std::remove_cvref_t<output_type>, hash_type, std::equal_to<key_type>>;

        // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //
    private:
        mutable table_type table_{};
        function_type function_;

        // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //
    public:
        invocable_function_hash_table() = delete;
        ~invocable_function_hash_table() = default;

        invocable_function_hash_table(const invocable_function_hash_table &) = delete;
        invocable_function_hash_table &operator=(const invocable_function_hash_table &) = delete;

        invocable_function_hash_table(invocable_function_hash_table &&) noexcept = default;
        invocable_function_hash_table &
        operator=(invocable_function_hash_table &&) noexcept = default;

        explicit invocable_function_hash_table(function_type function) :
                function_{std::move(function)} {}

        invocable_function_hash_table(function_type function, const invocable<signature_type> &) :
                function_{std::move(function)} {}

        invocable_function_hash_table(function_type function, const hash_type &hash) :
                function_{std::move(function)},
                table_{0, hash} {}

        invocable_function_hash_table(function_type function, const hash_type &hash,
                                      const invocable<signature_type> &) :
                function_{std::move(function)},
                table_{0, hash} {}

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

            return std::get<1>(*std::get<0>(table_.emplace
                (std::make_pair(inputs_tuple, std::apply(function_, inputs_tuple)))));
        }

        // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //
    public:
        const function_type &function() const & { return function_; }
        function_type &&function() && { return std::move(function_); }

        const table_type &table() const & { return table_; }
        table_type &&table() && { return std::move(table_); }

        invocable_function_hash_table &clear() { table_.clear(); return *this; }
    };

    // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //

    template<class function_type>
    invocable_function_hash_table(function_type) ->
    invocable_function_hash_table<invocable_signature<function_type>, function_type, Core::hash>;

    template<class function_type, class signature_type>
    invocable_function_hash_table(function_type, invocable<signature_type>) ->
    invocable_function_hash_table<signature_type, function_type, Core::hash>;

    template<class function_type, class hash_type>
    invocable_function_hash_table(function_type, hash_type) ->
    invocable_function_hash_table<invocable_signature<function_type>, function_type, hash_type>;

    template<class function_type, class hash_type, class signature_type>
    invocable_function_hash_table(function_type, hash_type, invocable<signature_type>) ->
    invocable_function_hash_table<signature_type, function_type, hash_type>;
}
