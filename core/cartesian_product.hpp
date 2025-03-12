#pragma once

#include "utility.hpp"
#include "traits.hpp"

// **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //

namespace Utility
{
    namespace Apparatus
    {
        template<class function_type, class... arguments_types>
        using cartesian_product_return_type =
                decltype(std::declval<function_type>()
                        (std::declval<typename arguments_types::value_type>()...));

        template<Arithmetic::index_t dim = 0, class function_type, class... arguments_types>
        void cartesian_product_apparatus
                (auto &ret,
                 const function_type &function,
                 const std::tuple<arguments_types const *...> &arguments,
                 std::tuple<typename arguments_types::value_type const *...> current = {})
        {
            if constexpr (dim == sizeof...(arguments_types))
                ret.emplace_back(std::apply(
                        [&function]<class... input_types>(const input_types &...inputs)
                                requires((std::is_same_v
                                        <input_types,
                                         typename arguments_types::value_type const *> and ...))
                        {
                            return function(*inputs...);
                        },
                        current));

            else
                for (auto it_bis = std::cbegin(*std::get<dim>(arguments));
                     it_bis != std::end(*std::get<dim>(arguments)); ++it_bis)
                {
                    std::get<dim>(current) = &*it_bis;
                    cartesian_product_apparatus<dim + 1>(ret, function, arguments, current);
                }

            return;
        }

        template<class function_type, class... arguments_types, class output_type
        = cartesian_product_return_type<function_type, arguments_types...>>
        std::vector<output_type> cartesian_product_apparatus
                (const function_type &function, arguments_types *const ... arguments)
        {
            Arithmetic::index_t n = (std::size(*arguments) * ...);
            std::vector<output_type> ret;
            ret.reserve(n);

            cartesian_product_apparatus(ret, function, std::tuple{arguments...});

            return ret;
        }
    }

    template<class function_type, Traits::container_c... arguments_types, class output_type
    = Apparatus::cartesian_product_return_type<function_type, arguments_types...>>
    std::vector<output_type> cartesian_product
            (const function_type &function, const arguments_types &... arguments)
    {
        return Apparatus::cartesian_product_apparatus(function, &arguments...);
    }

    template<Traits::container_c... arguments_types>
    std::vector<std::tuple<typename arguments_types::value_type...>>
    cartesian_product(const arguments_types &... arguments)
    {
        static constexpr auto lambda =
                []<class... input_types>(const input_types &...inputs)
                        requires((std::is_same_v
                                <input_types, typename arguments_types::value_type> and ...))
                {
                    return std::tuple{inputs...};
                };

        return Apparatus::cartesian_product_apparatus(lambda, &arguments...);
    }

    template<Traits::container_c... arguments_types>
    requires((std::is_same_v<typename arguments_types::value_type,
             typename Utility::get_type<0, arguments_types...>::value_type> and ...))
    std::vector<std::array<typename Utility::get_type<0, arguments_types...>::value_type,
                           sizeof...(arguments_types)>>
    cartesian_product_array(const arguments_types &... arguments)
    {
        static constexpr auto lambda =
                []<class... input_types>(const input_types &...inputs)
                        requires((std::is_same_v
                                <input_types, typename arguments_types::value_type> and ...))
                {
                    return std::array{inputs...};
                };

        return Apparatus::cartesian_product_apparatus(lambda, &arguments...);
    }

    template<Traits::container_c arguments_type, size_t size>
    std::vector<std::array<typename arguments_type::value_type, size>>
    cartesian_product(const std::array<arguments_type, size> &arg)
    {
        auto lambda = [](const auto &...args)
        {
            return cartesian_product(args...);
        };

        return std::apply(lambda, arg);
    }

    template<Traits::container_c... arguments_types>
    std::vector<std::tuple<typename arguments_types::value_type...>>
    cartesian_product(const std::tuple<arguments_types...> &arg)
    {
        auto lambda = [](const auto &...args)
        {
            return cartesian_product(args...);
        };

        if constexpr
            (!(std::same_as<Utility::get_type<0, arguments_types...>, arguments_types> and ...))
            return std::apply(lambda, arg);

        else // convert std::vector<std::array> to std::vector<std::tuple>
        {
            auto vec = std::apply(lambda, arg);

            std::vector<std::tuple<typename arguments_types::value_type...>> ret;
            ret.reserve(std::size(vec));

            for (const auto &i : vec)
                ret.emplace_back(Utility::transform_tuple([](const auto &j){ return j; }, i));

            return ret;
        }
    }
}

// **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //
