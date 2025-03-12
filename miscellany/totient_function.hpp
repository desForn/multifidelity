#pragma once

#include "fwd.hpp"

namespace Miscellany
{
    index_t totient_function(index_t arg)
    {
        index_t ret = 1;

        for (index_t i = 2; i < arg; ++i)
            ret += std::gcd(i, arg) == 1;

        return ret;
    }
}
