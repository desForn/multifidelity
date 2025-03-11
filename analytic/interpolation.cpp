#include "fwd.hpp"
#include "random/random.hpp"
#include "invocable/invocable_function_hash_table.hpp"
#include "invocable/invocable_function_sequential_vector.hpp"

using namespace Analytic;
using namespace Smolyak;

namespace Apparatus
{
    struct smolyak_approximation_struct
    {
        template<class... types>
        static smolyak_approximation_factory<types...> function(Utility::homotype_pack<types...>);
    };
}

#ifndef N_VARIATES
static_assert(false, "Must define N_VARIATES");
#endif

constexpr index_t n_variates = N_VARIATES;

using array_type = std::array<real_t, n_variates>;
using tuple_type = Utility::uniform_tuple_type<real_t, n_variates>;
using smolyak_approximation_factory_type = decltype(
        ::Apparatus::smolyak_approximation_struct::function(
                Utility::generate_homotype_pack<interpolation_traits, n_variates>{}));

[[noreturn]] void exit()
{
    Print::println("Invalid arguments.");
    Print::println("\t[executable name] minimum_error n_test_points");

    std::abort();
}

int main(int argc, char *argv[])
{
    if (argc != 3)
        exit();

    real_t min_e = std::atof(argv[1]);
    index_t m = std::atoll(argv[2]);

    constexpr real_t epsilon = 0.03;

    auto function = [epsilon](const tuple_type &arg) -> real_t
    {
        array_type arg_f = Utility::transform_array
                ([epsilon](real_t arg_secundum){ return arg_secundum * (1 - epsilon); },
                 Utility::tuple_to_array(arg));

        return f(arg_f);
    };

    auto cost_function = [](const tuple_type &) -> real_t
    {
        return 1;
    };

    auto e_w_component_lambda = [](index_t arg) -> real_t
    {
        static constexpr real_t log_rho = std::log(rho(epsilon));
        static constexpr real_t log_2 = std::log(2);

        real_t a = log_2 * (static_cast<integer_t>(arg) * 2);
        return log_rho * std::exp(a) + a;
    };

    Invocable::invocable_function_sequential_vector e_w_component{e_w_component_lambda};

    auto e_w = [&e_w_component](const std::array<index_t, n_variates> &arg)
    {
        real_t ret = 0;
        std::ranges::for_each(arg, [&ret, &e_w_component](index_t a) { ret += e_w_component(a); });
        return ret;
    };

    auto smolyak_approximation = smolyak_approximation_factory_type::
            factory(function, cost_function);

    std::vector<real_t> error, cost;

    std::vector<tuple_type> test_points;

    Random::set_seed(0);
    for (index_t i = 0; i != m; ++i)
        test_points.emplace_back(Utility::array_to_tuple(Random::real_array<n_variates>(-1, 1)));

    real_t max_f = 0;
    for (const tuple_type &j: test_points)
        max_f = std::max(max_f, std::abs(function(j)));

    std::unordered_set<std::array<index_t, n_variates>, Core::hash> levels;
    std::multimap<real_t, std::array<index_t, n_variates>> ordered_levels;

    real_t L = e_w(std::array<index_t, n_variates>{});
    ordered_levels.insert(std::pair{L, std::array<index_t, n_variates>{}});
    auto ordered_levels_it = std::cbegin(ordered_levels);

    real_t allocated_L = L;

    std::ofstream data{std::string{"data_interpolation_"} + std::to_string(n_variates) + ".dat"};
    ASSERT_ASSUME(data.is_open());
    Disk::write(data, "cost error\n");

    for (index_t i = 0; true; ++i)
    {
        Print::println("Iteration: ", i);

        if (ordered_levels_it == std::cend(ordered_levels))
        {
            std::array<index_t, n_variates> level{};
            std::array<real_t, n_variates> e_w_level;
            std::ranges::transform(level, std::begin(e_w_level),
                    [&e_w_component](index_t arg) { return e_w_component(arg); });
            real_t L_level = std::accumulate(std::cbegin(e_w_level), std::cend(e_w_level),
                    static_cast<real_t>(0));

            auto next_level = [&level, &e_w_level, &L_level, L, &allocated_L, &e_w_component]()
            {
                auto set = [&level, &e_w_level, &L_level, &e_w_component]
                        (index_t i, index_t value) -> void
                {
                    level[i] = value;
                    L_level -= e_w_level[i];
                    e_w_level[i] = e_w_component(level[i]);
                    L_level += e_w_level[i];

                    return;
                };

                if (level == std::array<index_t, n_variates>{})
                {
                    while (L_level <= L)
                        set(n_variates - 1, level.back() + 1);

                    allocated_L = L_level * 2;

                    return true;
                }

                set(n_variates - 1, level.back() + 1);

                if (L_level <= allocated_L)
                    return true;

                for (index_t i = n_variates - 1; i != 0; --i)
                {
                    if (level[i - 1] >= level[i])
                        continue;

                    set(i - 1, level[i - 1] + 1);
                    for (index_t j = i; j != n_variates; ++j)
                        set(j, level[i - 1]);

                    while (L_level <= L)
                        set(n_variates - 1, level.back() + 1);
                    if (L_level <= allocated_L)
                        return true;
                }

                return false;
            };

            while (next_level())
                ordered_levels.insert(std::pair{L_level, level});

            ordered_levels_it = ordered_levels.find(L);
            ASSERT_ASSUME(ordered_levels_it != std::cend(ordered_levels));

            do
            {
                ++ordered_levels_it;
                ASSERT_ASSUME(ordered_levels_it != std::cend(ordered_levels));
            } while (levels.find(ordered_levels_it->second) != std::cend(levels));
        }

        while(ordered_levels_it != std::end(ordered_levels)) {
        L = ordered_levels_it->first;
        std::array<index_t, n_variates> level{ordered_levels_it->second};
        do
        {
            levels.insert(level);
        } while (std::next_permutation(std::begin(level), std::end(level)));

        ++ordered_levels_it;}

        smolyak_approximation.activate_nodes(levels);

        real_t e = 0;
        for (const tuple_type &j: test_points)
            e = std::max(e, std::abs(smolyak_approximation(j) - function(j)));
        e /= max_f;

        real_t w = smolyak_approximation.cost();

        error.emplace_back(e);
        cost.emplace_back(w);

        Disk::write(data, std::to_string(cost.back()), " ", std::to_string(error.back()), "\n");
        data << std::flush;

        Print::println("\tCost = ", cost.back(), "; Error = ", error.back());

        if (error.back() < min_e)
            break;
    }

    data.close();

    return 0;
}
