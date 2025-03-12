#pragma once

#include "core/core.hpp"

// **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //

namespace Sampling
{
    using Arithmetic::index_t;
    using Arithmetic::integer_t;
    using Arithmetic::real_t;
    using Arithmetic::complex_t;
    using Arithmetic::negative_1;

    // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //

    template<class sampling_type>
    concept sampling_c =
    requires(sampling_type sampling, index_t i)
    {
        typename sampling_type::coordinates_type;
        typename sampling_type::sampling_traits;

        { sampling_type{i} } -> std::same_as<sampling_type>;
        { sampling_type::factory(i) } -> std::same_as<const sampling_type &>;
        { sampling_type::clear_factory() } -> std::same_as<void>;

        { sampling_type::n_points(i) } -> std::same_as<index_t>;
        { sampling_type::required_level(i) } -> std::same_as<index_t>;

        { sampling.level() } -> std::same_as<index_t>;
        { sampling.n_points() } -> std::same_as<index_t>;
        { sampling.coordinates() }
            -> std::same_as<std::span<const typename sampling_type::coordinates_type>>;
        { std::move(sampling).coordinates() }
        -> std::same_as<std::vector<typename sampling_type::coordinates_type> &&>;
        { sampling[i] } -> std::same_as<const typename sampling_type::coordinates_type &>;

        { std::begin(sampling) } ->
        std::same_as<decltype(std::begin(std::move(sampling).coordinates()))>;
        { std::end(sampling) } ->
        std::same_as<decltype(std::end(std::move(sampling).coordinates()))>;

        { std::cbegin(sampling) } ->
        std::same_as<decltype(std::cbegin(std::move(sampling).coordinates()))>;
        { std::cend(sampling) } ->
        std::same_as<decltype(std::cend(std::move(sampling).coordinates()))>;

        { std::rbegin(sampling) } ->
        std::same_as<decltype(std::rbegin(std::move(sampling).coordinates()))>;
        { std::rend(sampling) } ->
        std::same_as<decltype(std::rend(std::move(sampling).coordinates()))>;

        { std::crbegin(sampling) } ->
        std::same_as<decltype(std::crbegin(std::move(sampling).coordinates()))>;
        { std::crend(sampling) } ->
        std::same_as<decltype(std::crend(std::move(sampling).coordinates()))>;

    };

    template<class sampling_traits>
    concept sampling_traits_c =
    requires(index_t i)
    {
        typename sampling_traits::coordinates_type;

        { sampling_traits::coordinates(i) } ->
        std::same_as<std::vector<typename sampling_traits::coordinates_type>>;

        { sampling_traits::n_points(i) } -> std::same_as<index_t>;
        { sampling_traits::required_level(i) } -> std::same_as<index_t>;
    };

    template<class sampling_traits>
    concept sampling_traits_new_level_c = sampling_traits_c<sampling_traits> and
    requires(index_t i)
    {
        { sampling_traits::n_new_points(i) } -> std::same_as<index_t>;
        { sampling_traits::is_new_point(i, i) } -> std::same_as<bool>;
    };

    template<class sampling_type>
    concept exponential_sampling_c = sampling_c<sampling_type> and
    requires()
    {
        requires(sampling_type::sampling_traits::exponential_sampling);
    };
}

// **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //

