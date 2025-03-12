#pragma once

#include "sampling_base.hpp"

// **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //

namespace Sampling
{
    template<real_t growth_ratio_>
    class exponential_refinement_traits
    {
    public:
        using coordinates_type = real_t;
        static constexpr real_t growth_ratio = growth_ratio_;

        // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //
    private:
        // It has to be used statically
        exponential_refinement_traits() = delete;
        ~exponential_refinement_traits() = delete;

        // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //
    public:
        static std::vector<coordinates_type> coordinates(index_t);
        static index_t n_points(index_t);
        static index_t required_level(index_t);
    };

    // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //
    // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //

    template<real_t growth_ratio>
    auto exponential_refinement_traits<growth_ratio>::coordinates(index_t level)
    -> std::vector<coordinates_type>
    {
        if (level == 0)
            return {};

        std::vector<coordinates_type> ret(level);

        auto it = std::begin(ret);
        *it = 1;

        while (it != std::end(ret) - 1)
            *++it = *it * growth_ratio;

        return ret;
    }

    template<real_t growth_ratio>
    index_t exponential_refinement_traits<growth_ratio>::n_points(index_t level)
    {
        return level;
    }

    template<real_t growth_ratio>
    index_t exponential_refinement_traits<growth_ratio>::required_level(index_t level)
    {
        return level;
    }

    template<real_t growth_ratio, index_t offset = 0>
    using exponential_refinement = sampling<exponential_refinement_traits<growth_ratio>, offset>;
}