#pragma once

#include "fwd.hpp"

// **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //

namespace Invocable
{
    constexpr auto n_ary_operator(const auto &n_ary_operator, auto &&... operator_arguments)
    {
        return [&]<class... input_types>(input_types &&... inputs) ->
            decltype(n_ary_operator(operator_arguments(std::forward<input_types>(inputs)...)...))
        {
            return n_ary_operator(operator_arguments(std::forward<input_types>(inputs)...)...);
        };
    }

    template<class n_ary_operator_type, class... operator_argument_types>
    constexpr auto n_ary_operator_copy
        (n_ary_operator_type &&n_ary_operator, operator_argument_types &&... operator_arguments)
    {
        return [n_ary_operator = std::forward<n_ary_operator_type>(n_ary_operator),
                ... operator_arguments = std::forward<operator_argument_types>(operator_arguments)]
                <class... input_types>(input_types &&... inputs) ->
            decltype(n_ary_operator(operator_arguments(std::forward<input_types>(inputs)...)...))
        {
            return n_ary_operator(operator_arguments(std::forward<input_types>(inputs)...)...);
        };
    }

    // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //
    namespace Operators
    {
        constexpr auto plus =
                []<class argument_type>(argument_type &&argument) ->
                        decltype(+(std::forward<argument_type>(argument)))
                {
                    return +(std::forward<argument_type>(argument));
                };

        constexpr auto minus =
                []<class argument_type>(argument_type &&argument) ->
                        decltype(-(std::forward<argument_type>(argument)))
                {
                    return -(std::forward<argument_type>(argument));
                };

        constexpr auto sum =
                []<class... argument_types>(argument_types &&... arguments) ->
                        decltype((std::forward<argument_types>(arguments) + ...))
                {
                    return (std::forward<argument_types>(arguments) + ...);
                };

        constexpr auto subtract =
                []<class lhs_type, class rhs_type>(lhs_type &&lhs, rhs_type &&rhs) ->
                        decltype(std::forward<lhs_type>(lhs) - std::forward<rhs_type>(rhs))
                {
                    return std::forward<lhs_type>(lhs) - std::forward<rhs_type>(rhs);
                };

        constexpr auto multiply =
                []<class... argument_types>(argument_types &&... arguments) ->
                        decltype((std::forward<argument_types>(arguments) * ...))
                {
                    return (std::forward<argument_types>(arguments) * ...);
                };

        constexpr auto divide =
                []<class lhs_type, class rhs_type>(lhs_type &&lhs, rhs_type &&rhs) ->
                        decltype(std::forward<lhs_type>(lhs) / std::forward<rhs_type>(rhs))
                {
                    return std::forward<lhs_type>(lhs) / std::forward<rhs_type>(rhs);
                };

        constexpr auto square =
                []<class argument_type>(argument_type &&argument) -> decltype(argument * argument)
                {
                    return argument * argument;
                };

        constexpr auto sqrt =
                []<class argument_type>(argument_type &&argument) ->
                        decltype(std::sqrt(std::forward<argument_type>(argument)))
                {
                    return std::sqrt(std::forward<argument_type>(argument));
                };

        constexpr auto abs =
                []<class argument_type>(argument_type &&argument) ->
                        decltype(std::abs(std::forward<argument_type>(argument)))
                {
                    return std::abs(std::forward<argument_type>(argument));
                };
    }

    // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //

    namespace Apparatus
    {
        template<class type>
        struct compose_apparatus
        {
            compose_apparatus() = delete;
            ~compose_apparatus() = default;

            compose_apparatus(const compose_apparatus &) = delete;
            compose_apparatus &operator=(const compose_apparatus &) = delete;

            explicit compose_apparatus(const type &value) : value_(value) {}

            // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //

            type &value_;

            template<class argument_type>
            constexpr auto
            operator^(argument_type &&argument) ->
            decltype(value_(std::forward<argument_type>(argument)))
            {
                return value_(std::forward<argument_type>(argument));
            }
        };
    }

    // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //

    template<class... argument_types>
    constexpr auto sum(argument_types &&... arguments) ->
    decltype(n_ary_operator(Operators::sum, std::forward<argument_types>(arguments)...))
    {
        return n_ary_operator(Operators::sum, std::forward<argument_types>(arguments)...);
    }

    template<class lhs_type, class rhs_type>
    constexpr auto subtract(lhs_type &&lhs, rhs_type &&rhs) ->
    decltype(n_ary_operator
        (Operators::subtract, std::forward<lhs_type>(lhs), std::forward<rhs_type>(rhs)))
    {
        return n_ary_operator
            (Operators::subtract, std::forward<lhs_type>(lhs), std::forward<rhs_type>(rhs));
    }

    template<class... argument_types>
    constexpr auto multiply(argument_types &&... arguments) ->
    decltype(n_ary_operator(Operators::multiply,std::forward<argument_types>(arguments)...))
    {
        return n_ary_operator(Operators::multiply, std::forward<argument_types>(arguments)...);
    }

    template<class lhs_type, class rhs_type>
    constexpr auto divide(lhs_type &&lhs, rhs_type &&rhs) ->
    decltype(n_ary_operator(Operators::divide,
                            std::forward<lhs_type>(lhs),
                            std::forward<rhs_type>(rhs)))
    {
        return n_ary_operator
            (Operators::divide, std::forward<lhs_type>(lhs), std::forward<rhs_type>(rhs));
    }

    template<class argument_type>
    constexpr auto square(argument_type &&argument) ->
    decltype(n_ary_operator(Operators::square, std::forward<argument_type>(argument)))
    {
        return n_ary_operator(Operators::square, std::forward<argument_type>(argument));
    }

    template<class argument_type>
    constexpr auto sqrt(argument_type &&argument) ->
    decltype(n_ary_operator(Operators::sqrt, std::forward<argument_type>(argument)))
    {
        return n_ary_operator(Operators::sqrt, std::forward<argument_type>(argument));
    }

    template<class argument_type>
    constexpr auto abs(argument_type &&argument) ->
    decltype(n_ary_operator(Operators::abs, std::forward<argument_type>(argument)))
    {
        return n_ary_operator(Operators::abs, std::forward<argument_type>(argument));
    }

    template<class... argument_types>
    constexpr auto compose(argument_types &&... arguments)
    {
        return [&arguments...]<class input_type>(input_type &&input) ->
                decltype((Apparatus::compose_apparatus<argument_types>{arguments} ^ ... ^
                            std::forward<input_type>(input)))
        {
            return (Apparatus::compose_apparatus<argument_types>{arguments} ^ ... ^
                        std::forward<input_type>(input));
        };
    }

    template<class... argument_types, index_t... sequence_indices>
    constexpr auto tensor_product_apparatus
        (std::integer_sequence<index_t, sequence_indices...>, const argument_types &... arguments)
    {
        return [&arguments...]<class... input_types>(input_types &&...inputs) ->
            decltype((... * Utility::get_value<sequence_indices>(arguments...)
                    (Utility::get_value<sequence_indices>(std::forward<input_types>(inputs)...))))
        requires(sizeof...(input_types) == sizeof...(arguments))
        {
            return (... * Utility::get_value<sequence_indices>(arguments...)
                    (Utility::get_value<sequence_indices>(std::forward<input_types>(inputs)...)));
        };
    }

    template<class... argument_types>
    constexpr auto tensor_product(const argument_types &...arguments)
    {
        return tensor_product_apparatus(std::make_integer_sequence<index_t, sizeof...(arguments)>{},
                                   arguments...);
    }
}
