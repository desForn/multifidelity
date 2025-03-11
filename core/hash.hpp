#pragma once

#include "fwd.hpp"
#include "boost/functional/hash.hpp"

// **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //

namespace Core
{
    struct hash
    {
        template<class... types>
        static std::size_t operator()(const types &...inputs)
        {
            return boost::hash_value(inputs...);
        }
    };
}

// **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //


