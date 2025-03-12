#pragma once

#include "fwd.hpp"
#include "core/missing_data.hpp"

// **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //

namespace Invocable
{
    template<class function_type, class inputs_type, std::size_t size, class outputs_type =
             decltype(std::declval<function_type>()(std::declval<inputs_type>()))>
    std::array<outputs_type, size> evaluate
    (const function_type &function, const std::array<inputs_type, size> &args)
    {
        std::array<outputs_type, size> ret;
        
        try
        {
            auto it_ret = std::begin(ret);
            for (const inputs_type &i : args)
                *it_ret++ = function(i);
        }
        
        catch (Core::missing_data_base &)
        {
            Core::throw_missing_points(function, args);
        }
        
        return ret;
    }

    template<class function_type, class args_type, class outputs_type =
             decltype(std::declval<function_type>()(*std::cbegin(std::declval<args_type>())))>
    std::vector<outputs_type> evaluate
    (const function_type &function, const args_type &args)
    {
        std::vector<outputs_type> ret;
        ret.reserve(std::size(args));
        
        try
        {
            for (const auto &i : args)
                ret.emplace_back(function(i));
        }
        
        catch (Core::missing_data_base &)
        {
            Core::throw_missing_points(function, args);
        }
        
        return ret;
    }
}
