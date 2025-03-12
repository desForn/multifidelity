#pragma once

#include "fwd.hpp"
#include "sampling_base.hpp"

// **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //

namespace Sampling
{
    constexpr index_t n_points_exponential_sampling(index_t);
    constexpr index_t n_new_points_exponential_sampling(index_t);
    constexpr index_t required_level_exponential_sampling(index_t);

    template<sampling_traits_c simple_sampling_traits_>
    class exponential_traits
    {
    public:
        using simple_sampling_traits = simple_sampling_traits_;
        using coordinates_type = simple_sampling_traits::coordinates_type;
        static constexpr bool exponential_sampling = true;

        // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //
    private:
        // It has to be used statically
        exponential_traits() = delete;
        ~exponential_traits() = delete;

        // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //
    public:
        static std::vector<coordinates_type> coordinates(index_t);
        static index_t n_points(index_t);
        static index_t required_level(index_t);

        static index_t n_new_points(index_t);
        static bool is_new_point(index_t, index_t);
    };

    // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //
    // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //

    template<sampling_traits_c simple_sampling_traits>
    auto exponential_traits<simple_sampling_traits>::coordinates(index_t level) ->
    std::vector<coordinates_type>
    {
        return simple_sampling_traits::coordinates(n_points(level));
    }

    template<sampling_traits_c simple_sampling_traits>
    index_t exponential_traits<simple_sampling_traits>::n_points(index_t level)
    {
        return n_points_exponential_sampling(level);
    }

    template<sampling_traits_c simple_sampling_traits>
    index_t exponential_traits<simple_sampling_traits>::required_level(index_t n_points)
    {
        return required_level_exponential_sampling(n_points);
    }

    template<sampling_traits_c simple_sampling_traits>
    index_t exponential_traits<simple_sampling_traits>::n_new_points(index_t level)
    {
        return n_new_points_exponential_sampling(level);
    }

    template<sampling_traits_c simple_sampling_traits>
    bool exponential_traits<simple_sampling_traits>::is_new_point(index_t i, index_t level)
    {
        if (level < 2)
            return i != 1;

        return i % 2;
    }

    // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //
    // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //

    constexpr index_t n_points_exponential_sampling(index_t level)
    {
        return Arithmetic::exp2(level) + (level != 0);
    }

    constexpr index_t n_new_points_exponential_sampling(index_t level)
    {
        return n_points_exponential_sampling(level) / 2 + (level < 2);
    }

    constexpr index_t required_level_exponential_sampling(index_t n_points)
    {
        if (n_points < 3)
            return n_points > 1;

        return Arithmetic::log2_ceil(n_points - 1);
    }

    // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //

    template<sampling_c sampling_type>
    using exponential = sampling<exponential_traits<typename sampling_type::sampling_traits>>;
}

// **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //

