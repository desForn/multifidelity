#pragma once

#include "fwd.hpp"
#include "core/missing_data.hpp"

// **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //

namespace Invocable
{
    template<class function_type, class args_type,
        class exception_container_type =
            std::unordered_set<typename args_type::value_type, Core::hash>,
        class outputs_type =
            decltype(std::declval<function_type>()(std::declval<typename args_type::value_type>()))>
    std::vector<outputs_type> evaluate
    (const function_type &function, const args_type &args, const index_t n_threads = 1)
    {
        ASSERT_ASSUME(n_threads != 0);

        std::vector<outputs_type> ret(std::size(args));

        auto work = [&function, &args, &ret, n_threads](const index_t arg) ->
            Core::missing_data<exception_container_type>
        {
            Core::missing_data<exception_container_type> missing_data_;

            auto it_args = std::cbegin(args) + arg;
            auto it_ret = std::begin(ret) + arg;

            for (; it_ret < std::cend(ret);
                    std::advance(it_args, n_threads), it_ret += n_threads)
            {
                try { *it_ret = function(*it_args); }
                catch (Core::missing_data_base &)
                {
                    missing_data_.data().emplace(*it_args);
                }
            }

            return missing_data_;
        };

        std::vector<decltype(std::async(std::launch::async, work, 0))> futures;

        for (index_t i = 0; i != n_threads and i != std::size(args); ++i)
            futures.emplace_back(std::async(std::launch::async, work, i));

        Core::missing_data<exception_container_type> missing_data_;

        for (auto &i : futures)
            missing_data_.data().insert_range(i.get().data());

        if (not std::empty(missing_data_.data()))
                throw missing_data_;

        return ret;
    }
}
