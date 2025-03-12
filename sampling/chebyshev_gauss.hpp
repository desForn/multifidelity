#pragma once

#include "sampling_base.hpp"

// **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //

namespace Sampling
{
    class chebyshev_gauss_traits
    {
    public:
        using coordinates_type = real_t;

        // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //
    private:
        // It has to be used statically
        chebyshev_gauss_traits() = delete;
        ~chebyshev_gauss_traits() = delete;

        // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //
    public:
        static std::vector<coordinates_type> coordinates(index_t);
        static index_t n_points(index_t);
        static index_t required_level(index_t);
    };

    // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //
    // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //

    auto chebyshev_gauss_traits::coordinates(index_t level) -> std::vector<coordinates_type>
    {
        if (level == 1)
            return {0};

        std::vector<coordinates_type> ret(level);

        for (index_t i = 0; i != level; ++i)
            ret[i] = std::sin(((2 * i + 1) / static_cast<real_t>(2 * level) - (0.5)) * M_PI);

        return ret;
    }

    index_t chebyshev_gauss_traits::n_points(index_t level)
    {
        return level;
    }

    index_t chebyshev_gauss_traits::required_level(index_t level)
    {
        return level;
    }

    using chebyshev_gauss = sampling<chebyshev_gauss_traits>;
}
