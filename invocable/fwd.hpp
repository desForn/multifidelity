#pragma once

#include "core/core.hpp"

namespace Invocable
{
    using Arithmetic::index_t;
    using Arithmetic::integer_t;

    using Arithmetic::real_t;
    using Arithmetic::complex_t;

    // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //

    template<class>
    struct invocable {};

    template<class, class>
    class invocable_function;

    template<class, class>
    class invocable_function_sequential_vector;

    template<class, class>
    class invocable_function_disk_sequential_vector;

    template<class, class, class>
    class invocable_function_hash_table;

    template<class, class>
    class invocable_function_binary_tree;

    template<class>
    class invocable_disk_binary_tree;

    template<class, class>
    class invocable_function_disk_binary_tree;

    template<class, class, class...>
    class piecewise_invocable;

    // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //

    namespace Apparatus
    {
        template<index_t, class, class...>
        struct is_callable_apparatus : std::false_type {};

        template<index_t i, class signature_type, class... call_types>
        requires(Utility::count_signature_input_types<signature_type> == sizeof...(call_types))
        struct is_callable_apparatus<i, signature_type, call_types...>
        {
            static constexpr bool value =
                    std::is_convertible_v<
                            Utility::get_type<i, call_types...>,
                            Utility::get_signature_input_type<i, signature_type>> and
                    is_callable_apparatus<i + 1, signature_type, call_types...>::value;
        };

        template<class signature_type, class... call_types>
        requires(Utility::count_signature_input_types<signature_type> == sizeof...(call_types))
        struct is_callable_apparatus<sizeof...(call_types) - 1, signature_type, call_types...>
        {
            static constexpr bool value =
                    std::is_convertible_v<
                            Utility::get_type<sizeof...(call_types) - 1, call_types...>,
                            Utility::get_signature_input_type
                                    <sizeof...(call_types) - 1, signature_type>>;
        };

        template<class signature_type>
        requires(Utility::count_signature_input_types<signature_type> == 0)
        struct is_callable_apparatus<0, signature_type>
        {
            static constexpr bool value = true;
        };
    }

    template<class signature_type, class... call_types>
    concept is_callable_c = Apparatus::is_callable_apparatus<0, signature_type, call_types...>::value;

    // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //

    namespace Apparatus
    {
        template<class type_>
        struct wrapper
        {
            using type = type_;
        };
    
        template<class type>
        concept particular_invocable_c = requires { typename type::signature_type; };

        // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //

        std::false_type convert_to_function(...);

        template<particular_invocable_c invocable_type>
        std::function<typename invocable_type::signature_type>
        convert_to_function(const invocable_type &);

        template<class invocable_type,
                 class converted_function = decltype(std::function(std::declval<invocable_type>()))>
        requires(!particular_invocable_c<invocable_type>)
        converted_function convert_to_function(const invocable_type &);

        // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //

        template<class output_type, class... input_types>
        wrapper<output_type(input_types...)> invocable_signature_apparatus
            (const std::function<output_type(input_types...)> &);

        std::false_type invocable_signature_apparatus(const std::false_type &);
    }

    template<class type>
    concept general_invocable_c =
            !std::is_same_v
                    <decltype(Apparatus::convert_to_function(std::declval<type>())),
                     std::false_type>;

    namespace Apparatus
    {
        template<general_invocable_c invocable_type>
        using non_decayed_invocable_signature =
                typename decltype(Apparatus::invocable_signature_apparatus
                    (Apparatus::convert_to_function(std::declval<invocable_type>())))::type;
    }

    // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //
    namespace Apparatus
    {
        template<class non_signature_type>
        struct invocable_signature_c_apparatus : std::false_type {};

        template<class output_type, class... input_types>
        struct invocable_signature_c_apparatus<output_type(input_types...)> : std::true_type {};
    }

    template<class type>
    concept invocable_signature_c = Apparatus::invocable_signature_c_apparatus<type>::value;

    // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //
    namespace Apparatus
    {
        template<class>
        struct decayed_invocable_signature_apparatus;

        template<class non_invocable_type>
        struct decay_invocable_apparatus
        {
            using type = non_invocable_type;
        };

        template<general_invocable_c invocable_type>
        struct decay_invocable_apparatus<invocable_type>
        {
            using type = invocable
                <typename decayed_invocable_signature_apparatus
                        <non_decayed_invocable_signature<invocable_type>>::type>;
        };

        template<class non_signature_type>
        struct decayed_invocable_signature_apparatus
        {
            using type = typename decay_invocable_apparatus<non_signature_type>::type;
        };

        template<class output_type, class... input_types>
        struct decayed_invocable_signature_apparatus<output_type(input_types...)>
        {
            using type =
                    output_type(typename decayed_invocable_signature_apparatus<input_types>::type...);
        };
    }

    template<class invocable_type>
    using decay_invocable = typename Apparatus::decay_invocable_apparatus<invocable_type>::type;

    namespace Apparatus
    {
        template<invocable_signature_c signature_type>
        using decayed_invocable_signature =
                typename Apparatus::decayed_invocable_signature_apparatus<signature_type>::type;
    }

    template<general_invocable_c invocable_type>
    using invocable_signature =
            Apparatus::decayed_invocable_signature
                <Apparatus::non_decayed_invocable_signature<invocable_type>>;

    // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //
    namespace Apparatus
    {
        template<class, class>
        struct invocable_c_apparatus : std::false_type {};

        template<class invocable_type, class output_type, class... input_types>
        struct invocable_c_apparatus<invocable_type, output_type(input_types...)>
        {
            static constexpr bool value =
                std::is_invocable_r_v<output_type, invocable_type, input_types...>;
        };

        template<class invocable_type, class output_type, class... input_types>
        requires(invocable_c_apparatus<invocable_type, output_type(input_types...)>::value)
        std::true_type invocable_multiple_inputs_r_c_apparatus
            (const invocable_type &, const output_type &,
             const Utility::homotype_pack<input_types...> &);

        std::false_type invocable_multiple_inputs_r_c_apparatus(...);

        template<class invocable_type, class... input_types>
        requires(std::is_invocable_v<invocable_type, input_types...>)
        std::true_type invocable_multiple_inputs_c_apparatus
            (const invocable_type &, const Utility::homotype_pack<input_types...> &);

        std::false_type invocable_multiple_inputs_c_apparatus(...);
    }

    template<class invocable_type, class signature_type>
    concept invocable_c = Apparatus::invocable_c_apparatus<invocable_type, signature_type>::value;

    template<class invocable_type, class input_type, index_t n>
    concept invocable_multiple_inputs_c =
            decltype(Apparatus::invocable_multiple_inputs_c_apparatus
                (std::declval<invocable_type>(),
                 Utility::generate_homotype_pack<input_type, n>{}))::value;

    template<class invocable_type, class output_type, class input_type, index_t n>
    concept invocable_multiple_inputs_r_c =
            decltype(Apparatus::invocable_multiple_inputs_r_c_apparatus
            (std::declval<invocable_type>(),
             std::declval<output_type>(),
             Utility::generate_homotype_pack<input_type, n>{}))::value;
}
