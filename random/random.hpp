#pragma once

#include "fwd.hpp"

// **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //

namespace Random
{
    namespace Apparatus
    {
        unsigned int generate_seed()
        {
            if constexpr (default_seed == negative_1)
                return static_cast<unsigned int>(time(nullptr));

            else
                return static_cast<unsigned int>(default_seed);
        }

        static std::mt19937 random_generator_{generate_seed()};
        static std::normal_distribution<real_t> standard_normal_distribution_{};
    }

    void set_seed(index_t seed)
    {
        Apparatus::random_generator_.seed(static_cast<unsigned int>(seed));
    }

    void set_seed_time()
    {
        set_seed(time(nullptr));
    }

    index_t engine()
    {
        return Apparatus::random_generator_();
    }

    integer_t integer(integer_t min, integer_t max)
    {
        integer_t d = max - min + 1;

        ASSERT_ASSUME(d > 0);

        return engine() % d + min;
    }

    real_t real()
    {
        return static_cast<real_t>(engine() - std::mt19937::min()) /
               static_cast<real_t>(std::mt19937::max() - std::mt19937::min());
    }

    real_t real(real_t min, real_t max)
    {
        real_t d = max - min;

        ASSERT_ASSUME(d > static_cast<real_t>(0));

        return real() * d + min;
    }

    template<index_t n_variates>
    std::array<real_t, n_variates> real_array(real_t min, real_t max)
    {
        std::array<real_t, n_variates> ret;

        std::for_each(std::begin(ret), std::end(ret),
                      [min, max](auto &i) -> void
                      {
                          i = real(min, max);
                          return;
                      });

        return ret;
    }

    template<index_t n_variates>
    std::array<real_t, n_variates> real_array()
    {
        return real_array<n_variates>(0, 1);
    }

    real_t standard_normal_distribution()
    {
        return Apparatus::standard_normal_distribution_(Apparatus::random_generator_);
    }

    real_t normal_distribution(real_t mean, real_t standard_deviation)
    {
        return standard_normal_distribution() * standard_deviation + mean;
    }
}
