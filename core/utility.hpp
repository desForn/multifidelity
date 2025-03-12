#pragma once

#include "arithmetic_fwd.hpp"

namespace Utility
{
    namespace Apparatus
    {
        template<auto first_argument, auto... other_arguments>
        struct sum_apparatus
        {
            static constexpr auto value = first_argument + sum_apparatus<other_arguments...>::value;
        };

        template<auto first_argument>
        struct sum_apparatus<first_argument>
        {
            static constexpr auto value = first_argument;
        };
    }

    template<auto... arguments>
    constexpr auto sum = Apparatus::sum_apparatus<arguments...>::value;

    // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //

    namespace Apparatus
    {
        template<Arithmetic::index_t i, class first_type, class... other_types>
        struct get_type_apparatus
        {
            static_assert(i < sizeof...(other_types) + 1, "Utility::get_type: out of range");

            using type = typename get_type_apparatus<i - 1, other_types...>::type;
        };

        template<class first_type, class... other_types>
        struct get_type_apparatus<0, first_type, other_types...>
        {
            using type = first_type;
        };
    }

    template<Arithmetic::index_t i, class... types>
    using get_type = typename Apparatus::get_type_apparatus<i, types...>::type;

    // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //

    template<class first_input_type, class... input_types>
    requires((std::is_same_v<first_input_type, input_types> and ...))
    auto &&get_value
    (Arithmetic::index_t i, first_input_type &&first_input, input_types &&...inputs) noexcept
    {
        if (i == 0)
            return std::forward<first_input_type>(first_input);

        if constexpr (sizeof...(input_types) != 0)
            return get_value(i - 1, std::forward<input_types>(inputs)...);

        throw std::out_of_range("Utility::get_value");
    }

    template<Arithmetic::index_t i, class first_type, class... other_types>
    constexpr get_type<i, first_type, other_types...> &&
    get_value(first_type &&first, other_types &&... others) noexcept
    {
        static_assert(i < sizeof...(other_types) + 1, "Utility::get_value: out of range");

        if constexpr (i == 0)
            return std::forward<first_type>(first);

        else
            return get_value<i - 1, other_types...>(std::forward<other_types>(others)...);
    }

    // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //

    namespace Apparatus
    {
        template<class type, class first_type, class... other_types>
        constexpr Arithmetic::index_t count_value_apparatus
                (Arithmetic::index_t counter, const type &value_to_count,
                 const first_type &first_arg, const other_types &... other_args)
        {
            if (value_to_count == first_arg)
                ++counter;

            if constexpr (sizeof...(other_types) == 0)
                return counter;

            else
                return count_value_apparatus(counter, value_to_count, other_args...);
        }
    }

    template<class type, class... types>
    requires ((... and std::equality_comparable_with<type, types>))
    constexpr Arithmetic::index_t count_value(const type &value_to_count, const types &... args)
    {
        return Apparatus::count_value_apparatus(0, value_to_count, args...);
    }

    // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //

    namespace Apparatus
    {
        template<class type_to_count, class first_type_to_search, class... other_types_to_search>
        consteval Arithmetic::index_t count_type_apparatus(Arithmetic::index_t current_counter)
        {
            if constexpr (std::is_same_v<type_to_count, first_type_to_search>)
                ++current_counter;

            if constexpr (sizeof...(other_types_to_search) != 0)
                current_counter = count_type_apparatus<type_to_count, other_types_to_search...>
                        (current_counter);

            return current_counter;
        }

        template<class... Args>
        consteval Arithmetic::index_t count_type_apparatus()
        {
            if constexpr (sizeof...(Args) == 1)
                return 0;

            else
                return count_type_apparatus<Args...>(0);
        }
    }

    template<class type_to_count, class... list_of_types>
    constexpr Arithmetic::index_t count_type = Apparatus::count_type_apparatus
            <type_to_count, list_of_types...>();

    // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //

    namespace Apparatus
    {
        template<class type_to_count, class first_type_to_search, class... other_types_to_search>
        consteval Arithmetic::index_t count_convertible_type_apparatus
            (Arithmetic::index_t current_counter)
        {
            if constexpr (std::is_convertible_v<type_to_count, first_type_to_search>)
                ++current_counter;

            if constexpr (sizeof...(other_types_to_search) != 0)
                current_counter = count_convertible_type_apparatus
                        <type_to_count, other_types_to_search...>(current_counter);

            return current_counter;
        }

        template<class... Args>
        consteval Arithmetic::index_t count_convertible_type_apparatus()
        {
            if constexpr (sizeof...(Args) == 1)
                return 0;

            else
                return count_convertible_type_apparatus<Args...>(0);
        }
    }

    template<class type_to_count, class... list_of_types>
    constexpr Arithmetic::index_t count_convertible_type =
            Apparatus::count_convertible_type_apparatus<type_to_count, list_of_types...>();

    // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //

    namespace Apparatus
    {
        template<Arithmetic::index_t i, class type_to_search, class... list_of_types>
        consteval Arithmetic::index_t find_type_apparatus()
        {
            if constexpr (i == 0)
                return Arithmetic::negative_1;

            if constexpr (std::is_same_v<type_to_search, get_type<i, list_of_types...>>)
                return i;

            else
                return find_type_apparatus<i + 1, type_to_search, list_of_types...>();
        }
    }

    template<class type_to_search, class... list_of_types>
    constexpr Arithmetic::index_t find_type = Apparatus::find_type_apparatus
            <0, type_to_search, list_of_types...>();

    // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //

    namespace Apparatus
    {
        template<Arithmetic::index_t i, class type_to_search, class... list_of_types>
        consteval Arithmetic::index_t find_convertible_type_apparatus()
        {
            if constexpr (i == 0)
                return Arithmetic::negative_1;

            if constexpr (std::is_convertible_v<type_to_search, get_type<i, list_of_types...>>)
                return i;

            else
                return find_convertible_type_apparatus
                        <i + 1, type_to_search, list_of_types...>();
        }
    }

    template<class type_to_search, class... list_of_types>
    constexpr Arithmetic::index_t find_convertible_type =
            Apparatus::find_convertible_type_apparatus<0, type_to_search, list_of_types...>();

    // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //

    namespace Apparatus
    {
        template<Arithmetic::index_t i, Arithmetic::index_t j,
                 class type_to_search, class... list_of_types>
        consteval Arithmetic::index_t find_nth_type_apparatus()
        {
            if constexpr (i == 0)
                return Arithmetic::negative_1;

            if constexpr (std::is_same_v<type_to_search, get_type<i, list_of_types...>>)
                if constexpr (j == 0)
                    return i;
                else
                    return find_nth_type_apparatus<
                            i + 1, j - 1, type_to_search, list_of_types...>();

            else
                return find_nth_type_apparatus<i + 1, j, type_to_search, list_of_types...>();
        }
    }

    template<Arithmetic::index_t nth, class type_to_search, class... list_of_types>
    constexpr Arithmetic::index_t find_nth_type =
            Apparatus::find_nth_type_apparatus<0, nth, type_to_search, list_of_types...>();

    // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //

    namespace Apparatus
    {
        template<Arithmetic::index_t i, Arithmetic::index_t j,
                 class type_to_search, class... list_of_types>
        consteval Arithmetic::index_t find_nth_convertible_type_apparatus()
        {
            if constexpr (i == 0)
                return Arithmetic::negative_1;

            if constexpr (std::is_convertible_v<type_to_search, get_type<i, list_of_types...>>)
                if constexpr (j == 0)
                    return i;
                else
                    return find_nth_convertible_type_apparatus<
                            i + 1, j - 1, type_to_search, list_of_types...>();

            else
                return find_nth_convertible_type_apparatus<
                        i + 1, j, type_to_search, list_of_types...>();
        }
    }

    template<Arithmetic::index_t nth, class type_to_search, class... list_of_types>
    constexpr Arithmetic::index_t find_nth_convertible_type =
            Apparatus::find_nth_convertible_type_apparatus
                    <0, nth, type_to_search, list_of_types...>();

    // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //

    namespace Apparatus
    {
        template<Arithmetic::index_t i, class first_type, class... other_types>
        struct get_non_void_type_apparatus
        {
            using type = void;
        };

        template<Arithmetic::index_t i, class first_type, class... other_types>
        requires(i < sizeof...(other_types) + 1 - count_type<void, other_types...>)
        struct get_non_void_type_apparatus<i, first_type, other_types...>
        {
            using type = typename get_non_void_type_apparatus<i - 1, other_types...>::type;
        };

        template<Arithmetic::index_t i, class... other_types>
        requires(i < sizeof...(other_types) - count_type<void, other_types...>)
        struct get_non_void_type_apparatus<i, void, other_types...>
        {
            using type = typename get_non_void_type_apparatus<i, other_types...>::type;
        };

        template<class first_type, class... other_types>
        struct get_non_void_type_apparatus<0, first_type, other_types...>
        {
            using type = first_type;
        };

        template<class... other_types>
        requires(sizeof...(other_types) != count_type<void, other_types...>)
        struct get_non_void_type_apparatus<0, void, other_types...>
        {
            using type = typename get_non_void_type_apparatus<0, other_types...>::type;
        };
    }

    template<Arithmetic::index_t i, class...types>
    using get_non_void_type = typename Apparatus::get_non_void_type_apparatus<i, types...>::type;

    // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //

    namespace Apparatus
    {
        template<class...>
        struct non_void_tuple_apparatus;

        template<>
        struct non_void_tuple_apparatus<>
        {
            using type = void;
        };

        template<class first_type, class...other_types>
        struct non_void_tuple_apparatus<first_type, other_types...>
        {
            using type = non_void_tuple_apparatus<std::tuple<first_type>, other_types...>::type;
        };

        template<class...other_types>
        struct non_void_tuple_apparatus<void, other_types...>
        {
            using type = non_void_tuple_apparatus<other_types...>::type;
        };

        template<class...types_in_tuple>
        struct non_void_tuple_apparatus<std::tuple<types_in_tuple...>>
        {
            using type = std::tuple<types_in_tuple...>;
        };

        template<class...types_in_tuple, class first_to_append, class... others_to_append>
        struct non_void_tuple_apparatus<std::tuple<types_in_tuple...>,
                first_to_append, others_to_append...>
        {
            using type = non_void_tuple_apparatus
                    <std::tuple<types_in_tuple..., first_to_append>, others_to_append...>::type;
        };

        template<class...types_in_tuple, class... others_to_append>
        struct non_void_tuple_apparatus<std::tuple<types_in_tuple...>, void, others_to_append...>
        {
            using type = non_void_tuple_apparatus
                    <std::tuple<types_in_tuple...>, others_to_append...>::type;
        };
    }

    template<class... types>
    using non_void_tuple = Apparatus::non_void_tuple_apparatus<types...>::type;

    // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //

    namespace Apparatus
    {
        template<Arithmetic::index_t i, class... types>
        constexpr auto head_tuple_apparatus(const std::tuple<types...> arg)
        {
            if constexpr (i == 0)
                return std::tuple{std::get<i>(arg)};
            else
                return std::tuple_cat(head_tuple_apparatus<i - 1, types...>(arg), std::tuple{std::get<i>(arg)});
        }

        template<Arithmetic::index_t i, class... types>
        constexpr auto tail_tuple_apparatus(const std::tuple<types...> arg)
        {
            if constexpr (i == sizeof...(types) - 2)
                return std::tuple{std::get<i + 1>(arg)};
            else
                return std::tuple_cat(std::tuple{std::get<i + 1>(arg)},
                        tail_tuple_apparatus<i + 1, types...>(arg));
        }
    }

    template<Arithmetic::index_t i, class... types>
    constexpr auto head_tuple(const std::tuple<types...> &arg)
    {
        static_assert(i < sizeof...(types) - 1, "Utility::head_tuple: out of range");
        return Apparatus::head_tuple_apparatus<i, types...>(arg);
    }

    template<Arithmetic::index_t i, class... types>
    constexpr auto tail_tuple(const std::tuple<types...> &arg)
    {
        static_assert(i < sizeof...(types) - 1, "Utility::tail_tuple: out of range");
        return Apparatus::tail_tuple_apparatus<i, types...>(arg);
    }

    template<Arithmetic::index_t i, class... types>
    constexpr auto split_tuple(const std::tuple<types...> &arg)
    {
        static_assert(i < sizeof...(types) - 1, "Utility::split: out of range");
        return std::tuple{head_tuple<i, types...>(arg), tail_tuple<i, types...>(arg)};
    }

    // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //

    namespace Apparatus
    {
        template<class signature_type>
        struct get_signature_output_type_apparatus
        {
            static_assert(false, "Utility::get_signature_outputtype: invalid type");
        };

        template<class output_type, class... input_types>
        struct get_signature_output_type_apparatus<output_type(input_types...)>
        {
            using type = output_type;
        };

        template<Arithmetic::index_t i, class signature_type>
        struct get_signature_input_type_apparatus
        {
            static_assert(false, "Utility::get_signature_input_type: invalid type");
        };

        template<Arithmetic::index_t i, class output_type, class... input_types>
        struct get_signature_input_type_apparatus<i, output_type(input_types...)>
        {
            using type = Utility::get_type<i, input_types...>;
        };

        template<class signature_type>
        struct count_signature_input_types_apparatus
        {
            static_assert(false, "Utility::count_signature_input_types: invalid type");
        };

        template<class output_type, class... input_types>
        struct count_signature_input_types_apparatus<output_type(input_types...)>
        {
            static constexpr Arithmetic::index_t value = sizeof...(input_types);
        };
    }

    template<class signature_type>
    using get_signature_output_type = typename Apparatus::get_signature_output_type_apparatus
            <signature_type>::type;

    template<Arithmetic::index_t i, class signature_type>
    using get_signature_input_type = typename Apparatus::get_signature_input_type_apparatus
            <i, signature_type>::type;

    template<class signature_type>
    constexpr Arithmetic::index_t
            count_signature_input_types = Apparatus::count_signature_input_types_apparatus
            <signature_type>::value;

    // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //

    namespace Apparatus
    {
        template<Arithmetic::index_t n, class output_type, class... input_types>
        struct multiple_inputs_signature_apparatus
        {
            using type = typename multiple_inputs_signature_apparatus<n - 1, output_type,
                    Utility::get_type<0, input_types...>,
                    input_types...>::type;
        };

        template<class output_type, class input_type>
        struct multiple_inputs_signature_apparatus<0, output_type, input_type>
        {
            using type = output_type();
        };

        template<class output_type, class... input_types>
        struct multiple_inputs_signature_apparatus<1, output_type, input_types...>
        {
            using type = output_type(input_types...);
        };
    }

    template<class output_type, class input_type, Arithmetic::index_t n_inputs>
    using multiple_inputs_signature =
            typename Apparatus::multiple_inputs_signature_apparatus
                    <n_inputs, output_type, input_type>::type;

    // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //

    template<class...>
    struct homotype_pack {};

    namespace Apparatus
    {
        template<Arithmetic::index_t i, class... types>
        struct generate_homotype_pack_apparatus
        {
            using type = typename generate_homotype_pack_apparatus
                    <i - 1, get_type<0, types...>, types...>::type;
        };

        template<class... types>
        struct generate_homotype_pack_apparatus<1, types...>
        {
            using type = homotype_pack<types...>;
        };

        template<class... types>
        struct generate_homotype_pack_apparatus<0, types...>
        {
            using type = homotype_pack<>;
        };
    }

    template<class type, Arithmetic::index_t i>
    using generate_homotype_pack = typename Apparatus::generate_homotype_pack_apparatus
            <i, type>::type;

    // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //

    namespace Apparatus
    {
        template<class>
        struct is_array_apparatus
        {
            static constexpr bool value = false;
        };

        template<class type, std::size_t size>
        struct is_array_apparatus<std::array<type, size>>
        {
            static constexpr bool value = true;
        };

        template<class>
        struct is_tuple_apparatus
        {
            static constexpr bool value = false;
        };

        template<class ...types>
        struct is_tuple_apparatus<std::tuple<types...>>
        {
            static constexpr bool value = true;
        };
    }

    template<class type>
    concept is_array = Apparatus::is_array_apparatus<type>::value;

    template<class type>
    concept is_tuple = Apparatus::is_tuple_apparatus<type>::value;

    // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //

    template<class type, size_t size>
    constexpr std::array<type, size> uniform_array(const type &value)
    {
        std::array<type, size> ret;
        std::ranges::fill(ret, value);

        return ret;
    }

    namespace Apparatus
    {
        struct uniform_tuple_type_apparatus
        {
            template<class... types>
            static std::tuple<types...> function(homotype_pack<types...>);
        };
    }

    template<class type, size_t size>
    using uniform_tuple_type = decltype(
        Apparatus::uniform_tuple_type_apparatus::function(generate_homotype_pack<type, size>{}));

    namespace Apparatus
    {
        template<class... types, class type, size_t i = 0>
        void uniform_tuple_apparatus(std::tuple<types...> &ret, const type &value)
        {
            std::get<i>(ret) = value;

            if constexpr (i < sizeof...(types) - 1)
                uniform_tuple_apparatus<types..., type, i + 1>(ret, value);

            return;
        }
    }

    template<class type, size_t size>
    constexpr uniform_tuple_type<type, size> uniform_tuple(const type &value)
    {
        uniform_tuple_type<type, size> ret;
        return Apparatus::uniform_tuple_apparatus(ret, value);
    }

    namespace Apparatus
    {
        template<size_t i, class... types>
        void tuple_to_array_apparatus(std::array<get_type<0, types...>, sizeof...(types)> &ret,
                                 const std::tuple<types...> &arg)
        {
            std::get<i>(ret) = std::get<i>(arg);

            if constexpr (i < sizeof...(types) - 1)
                tuple_to_array_apparatus<i + 1, types...>(ret, arg);

            return;
        }

        template<size_t i, class type, size_t size>
        void array_to_tuple_apparatus(uniform_tuple_type<type, size> &ret,
                                 const std::array<type, size> &arg)
        {
            std::get<i>(ret) = std::get<i>(arg);

            if constexpr (i < size - 1)
                array_to_tuple_apparatus<i + 1, type, size>(ret, arg);

            return;
        }
    }

    template<class... types>
    requires((std::same_as<get_type<0, types...>, types> && ...))
    std::array<get_type<0, types...>, sizeof...(types)> tuple_to_array
        (const std::tuple<types...> &arg)
    {
        std::array<get_type<0, types...>, sizeof...(types)> ret;
        Apparatus::tuple_to_array_apparatus<0, types...>(ret, arg);
        return ret;
    }

    template<class type, size_t size>
    uniform_tuple_type<type, size> array_to_tuple(const std::array<type, size> &arg)
    {
        uniform_tuple_type<type, size> ret;
        Apparatus::array_to_tuple_apparatus<0, type, size>(ret, arg);
        return ret;
    }

    // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //
    
    template<class function_type, class type, std::size_t size>
    constexpr auto transform_array
            (const function_type &arg0, const std::array<type, size> &arg1)
    -> std::array<decltype(arg0(arg1.front())), size>
    {
        return std::apply([&arg0](const auto &...args_bis)
                          {
                              return std::array{arg0(args_bis)...};
                          },
                          arg1);
    }

    template<class function_type, class... args_types>
    constexpr auto transform_array
            (const function_type &arg0, const args_types &...args)
    -> std::array
            <std::common_type_t<decltype(arg0(std::declval<get_type<0, args_types>>()))...>,
                    sizeof...(args_types)>
    {
        using ret_type = std::common_type_t
                <decltype(arg0(std::declval<get_type<0, args_types>>()))...>;

        return std::array<ret_type, sizeof...(args_types)>{arg0(args)...};
    }

    namespace Apparatus
    {
        template<Arithmetic::index_t... sequence>
        constexpr auto transform_tuple_apparatus
                (const auto &arg0, const auto &arg1,
                 std::integer_sequence<Arithmetic::index_t, sequence...>)
        {
            return std::tuple(arg0(std::get<sequence>(arg1))...);
        }
    }

    template<class function_type, class... arg_types>
    constexpr auto transform_tuple
            (const function_type &arg0, const std::tuple<arg_types...> &arg1)
    -> std::tuple<decltype(arg0(std::declval<arg_types>()...))>
    {
        return Apparatus::transform_tuple_apparatus
            (arg0, arg1, std::make_integer_sequence<Arithmetic::integer_t, sizeof...(arg_types)>{});
    }

    template<class function_type, class type, std::size_t size>
    constexpr auto transform_tuple
            (const function_type &arg0, const std::array<type, size> &arg1)
    {
        return std::apply(
                [&arg0](const auto &...args_bis)
                {
                    return std::tuple{arg0(args_bis)...};
                },
                arg1);
    }

    namespace Apparatus
    {
        template<Arithmetic::index_t... sequence>
        constexpr auto transform_tuple
                (const auto &arg0, const auto &arg1, const auto &arg2,
                 std::integer_sequence<Arithmetic::index_t, sequence...>)
        {
            return std::tuple{arg0(std::get<sequence>(arg1), std::get<sequence>(arg2))...};
        }
    }

    template<class function_type, class... arg1_types, class... arg2_types>
    requires(sizeof...(arg1_types) == sizeof...(arg2_types))
    constexpr auto transform_tuple
            (const function_type &arg0,
             const std::tuple<arg1_types...> &arg1,
             const std::tuple<arg2_types...> &arg2)
    {
        return Apparatus::transform_tuple
            (arg0, arg1, arg2,
             std::make_integer_sequence<Arithmetic::index_t, sizeof...(arg1_types)>{});
    }

    template<class function_type, class... arg1_types, class arg2_type, std::size_t arg2_size>
    requires(sizeof...(arg1_types) == arg2_size)
    constexpr auto transform_tuple
            (const function_type &arg0,
             const std::tuple<arg1_types...> &arg1,
             const std::array<arg2_type, arg2_size> &arg2)
    {
        return Apparatus::transform_tuple
                (arg0, arg1, arg2,
                 std::make_integer_sequence<Arithmetic::index_t, sizeof...(arg1_types)>{});
    }

    template<class function_type, class arg1_type, std::size_t arg1_size, class... arg2_types>
    requires(sizeof...(arg2_types) == arg1_size)
    constexpr auto transform_tuple
            (const function_type &arg0,
             const std::array<arg1_type, arg1_size> &arg1,
             const std::tuple<arg2_types...> &arg2)
    {
        return Apparatus::transform_tuple
                (arg0, arg1, arg2,
                 std::make_integer_sequence<Arithmetic::index_t, sizeof...(arg2_types)>{});
    }

    template<class function_type, class arg1_type, class arg2_type, std::size_t size>
    constexpr auto transform_tuple
            (const function_type &arg0,
             const std::array<arg1_type, size> &arg1,
             const std::array<arg2_type, size> &arg2)
    {
        return Apparatus::transform_tuple
                (arg0, arg1, arg2,
                 std::make_integer_sequence<Arithmetic::index_t, size>{});
    }

    // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //

    template<class type, size_t size, class type_bis = type>
    requires(requires(type a, const type_bis b) {{a += b};})
    constexpr std::array<type, size> &increment
        (std::array<type, size> &arg0, Arithmetic::index_t arg1, const type_bis &arg2 = 1)
    {
        ASSERT_ASSUME(arg1 < size);
        arg0[arg1] += arg2;

        return arg0;
    }

    // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //

    template<class type, size_t size, class type_bis = type>
    requires(requires(type a, const type_bis b) {{a += b};})
    constexpr std::array<type, size> increment_copy
        (const std::array<type, size> &arg0, Arithmetic::index_t arg1,
         const type_bis &arg2 = 1)
    {
        std::array<type, size> ret{arg0};
        increment(ret, arg1, arg2);

        return ret;
    }

    // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //

    template<class type, size_t size>
    constexpr std::array<type, size> sum_array
            (std::array<type, size> lhs, const std::array<type, size> &rhs)
    {
        std::transform(std::cbegin(lhs), std::cend(lhs), std::cbegin(rhs), std::begin(lhs),
                       std::plus{});

        return lhs;
    }

    template<class type, size_t size>
    constexpr std::array<type, size> subtract_array
            (std::array<type, size> lhs, const std::array<type, size> &rhs)
    {
        std::transform(std::cbegin(lhs), std::cend(lhs), std::cbegin(rhs), std::begin(lhs),
                       std::minus{});

        return lhs;
    }

    template<class type, size_t size>
    constexpr std::array<type, size> multiply_array
            (std::array<type, size> lhs, const std::array<type, size> &rhs)
    {
        std::transform(std::cbegin(lhs), std::cend(lhs), std::cbegin(rhs), std::begin(lhs),
                       std::multiplies{});

        return lhs;
    }

    template<class type, size_t size>
    constexpr std::array<type, size> divide_array
            (std::array<type, size> lhs, const std::array<type, size> &rhs)
    {
        std::transform(std::cbegin(lhs), std::cend(lhs), std::cbegin(rhs), std::begin(lhs),
                       std::divides{});

        return lhs;
    }

    template<class type, size_t size>
    constexpr std::array<type, size> modulus_array
            (std::array<type, size> lhs, const std::array<type, size> &rhs)
    {
        std::transform(std::cbegin(lhs), std::cend(lhs), std::cbegin(rhs), std::begin(lhs),
                       std::modulus{});

        return lhs;
    }

    template<class type, size_t size>
    constexpr std::array<type, size> negate_array(std::array<type, size> lhs)
    {
        std::transform(std::cbegin(lhs), std::cend(lhs), std::begin(lhs), std::negate{});

        return lhs;
    }

    // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //

    template<class type, size_t size>
    constexpr std::array<type, size> min_array
            (std::array<type, size> lhs, const std::array<type, size> &rhs)
    {
        std::transform(std::cbegin(lhs), std::cend(lhs), std::cbegin(rhs), std::begin(lhs),
                       [](const type &i, const type &j) { return std::min(i, j); });

        return lhs;
    }

    template<class type, size_t size>
    constexpr std::array<type, size> max_array
            (std::array<type, size> lhs, const std::array<type, size> &rhs)
    {
        std::transform(std::cbegin(lhs), std::cend(lhs), std::cbegin(rhs), std::begin(lhs),
                       [](const type &i, const type &j) { return std::max(i, j); });

        return lhs;
    }

    // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //

    template<class type, size_t size, class rhs_type>
    constexpr std::array<type, size> sum_array(std::array<type, size> lhs, const rhs_type &rhs)
    {
        std::ranges::for_each(lhs, [&rhs = rhs](type &i) { i += rhs; });

        return lhs;
    }

    template<class type, size_t size, class rhs_type>
    constexpr std::array<type, size> subtract_array(std::array<type, size> lhs, const rhs_type &rhs)
    {
        std::ranges::for_each(lhs, [&rhs = rhs](type &i) { i -= rhs; });

        return lhs;
    }

    template<class type, size_t size, class rhs_type>
    constexpr std::array<type, size> multiply_array(std::array<type, size> lhs, const rhs_type &rhs)
    {
        std::ranges::for_each(lhs, [&rhs = rhs](type &i) { i *= rhs; });

        return lhs;
    }

    template<class type, size_t size, class rhs_type>
    constexpr std::array<type, size> divide_array(std::array<type, size> lhs, const rhs_type &rhs)
    {
        std::ranges::for_each(lhs, [&rhs = rhs](type &i) { i /= rhs; });

        return lhs;
    }

    template<class type, size_t size, class rhs_type>
    constexpr std::array<type, size> modulus_array(std::array<type, size> lhs, const rhs_type &rhs)
    {
        std::ranges::for_each(lhs, [&rhs = rhs](type &i) { i %= rhs; });

        return lhs;
    }

    // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //

    template<class lhs_type, std::convertible_to<lhs_type> rhs_type, size_t size>
    constexpr std::array<lhs_type, size> static_cast_array(const std::array<rhs_type ,size> &rhs)
    {
        std::array<lhs_type, size> ret;

        auto it_lhs = std::cbegin(rhs);
        auto it_ret = std::begin(ret);

        while (it_lhs != std::cend(rhs))
            *it_ret++ = static_cast<lhs_type>(*it_lhs++);

        return ret;

    }

    template<size_t size>
    constexpr std::array<Arithmetic::integer_t, size> signed_array
        (const std::array<Arithmetic::index_t, size> &lhs)
    {
        return static_cast_array<Arithmetic::integer_t>(lhs);
    }

    template<size_t size>
    constexpr std::array<Arithmetic::index_t, size> unsigned_array
        (const std::array<Arithmetic::integer_t, size> &lhs)
    {
        return static_cast_array<Arithmetic::index_t>(lhs);
    }
    
    // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //
    // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //
    
    template<size_t size>
    constexpr std::array<Arithmetic::index_t, size> strides
        (const std::array<Arithmetic::index_t, size> &arg)
    {
        std::array<Arithmetic::index_t, size> ret;
        ret.back() = 1;
        
        std::transform(std::crbegin(ret), std::crend(ret) - 1,
                       std::crbegin(arg), std::rbegin(ret) + 1,
                       std::multiplies<Arithmetic::index_t>{});
            
        return ret;
    }
}

// **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //
