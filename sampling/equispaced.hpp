#pragma once

#include "sampling_base.hpp"
#include "../miscellany/totient_function.hpp"

// **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //

namespace Sampling
{
    class equispaced_traits
    {
    public:
        using coordinates_type = real_t;

        // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //
    private:
        // It has to be used statically
        equispaced_traits() = delete;
        ~equispaced_traits() = delete;

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

    auto equispaced_traits::coordinates(index_t level) -> std::vector<coordinates_type>
    {
        if (level == 1)
            return {0};

        std::vector<coordinates_type> ret(level);

        for (index_t i = 0; i != level; ++i)
            ret[i] = 2 * i / static_cast<real_t>(level - 1) - 1;

        return ret;
    }

    index_t equispaced_traits::n_points(index_t level)
    {
        return level;
    }

    index_t equispaced_traits::required_level(index_t level)
    {
        return level;
    }

    index_t equispaced_traits::n_new_points(index_t level)
    {
        if (level < 3)
            return level;

        if (level == 3)
            return 0;

        return Miscellany::totient_function(level - 1);
    }

    bool equispaced_traits::is_new_point(index_t i, index_t level)
    {
        ASSERT_ASSUME(i < level);

        if (level < 3)
            return true;

        if (level == 3)
            return false;

        return std::gcd(i, level - 1) == 1;
    }

    using equispaced = sampling<equispaced_traits>;
}
