#include "snapshot_invocable.hpp"

#include "smolyak/smolyak_approximation.hpp"
#include "smolyak/smolyak_traits.hpp"

using namespace Polynomial;
using namespace Smolyak;
using namespace Two_dimensional_flow;

using polynomial_point_value_type =
        multivariate_polynomial<chebyshev_polynomial_point_value, 1>;
using polynomial_coefficients_type = multivariate_polynomial<chebyshev_polynomial_i_kind, 1>;

using sampling_type = Sampling::exponential<Sampling::chebyshev_lobatto>;

// **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //

int main(int argc, char* argv[])
{
    if (argc != 2 and argc != 3)
    {
        Print::println("Must be called as:");
        Print::println("\t[executable name] [path to data file]");
        Print::println("\t[executable name] [path to data file] n_levels");
        std::abort();
    }

    std::string path{argv[1]};
    index_t n_levels = negative_1;

    if (argc == 3)
        n_levels = std::atol(argv[2]);

    snapshot_invocable snapshot_invocable_{path};

    auto qoi_function = [&snapshot_invocable_]
        (const std::tuple<index_t, real_t> &arg)
        {
            return snapshot_invocable_.qoi_function()
                (std::array<index_t, 1>{std::get<0>(arg)},
                 std::array<real_t, 1>{std::get<1>(arg)});
        };
    auto cost_function = [&snapshot_invocable_]
        (const std::tuple<index_t, real_t> &arg)
        {
            return snapshot_invocable_.cost_function()
                (std::array<index_t, 1>{std::get<0>(arg)},
                 std::array<real_t, 1>{std::get<1>(arg)});
        };

    auto smolyak = smolyak_approximation_factory
        <Smolyak::Smolyak_traits::multifidelity<>,
        Smolyak::Smolyak_traits::polynomial_interpolation<sampling_type>>::factory
        (qoi_function, cost_function);

    smolyak.maximum_level() = std::array<index_t, 2>{5, negative_1};

    index_t i = 0;
    try
    {
        for (; i != n_levels; ++i)
            smolyak.activate_node_adaptively(verbosity::all);
    }

    catch(Core::missing_data<std::vector<std::tuple<index_t, real_t>>> &e)
    {
        Print::println("Missing data: ");
        Print::println(e.data());

        std::ofstream missing_data{"missing_data.dat", std::ios::out};
        Disk::write(missing_data, "Fidelity_level angle_index\n");

        for (const std::tuple<index_t, real_t> &i : e.data())
            Disk::write(missing_data, std::to_string(std::get<0>(i)), " ",
                    std::to_string(sampling_type::index(std::get<1>(i))), "\n");
    }

    Print::println("Fitted ", i, " levels.");

    for (const auto &i : smolyak.nodes())
        if (i.second.state() != Smolyak::node_state::unfitted)
            smolyak.activate_node(i.first);

    std::ofstream file{"surrogate.dat"};
    Disk::write(file, "theta F_x\n");
    
    index_t n_data = 2001;
    
    for (index_t i = 0; i != n_data; ++i)
    {
        real_t x = static_cast<real_t>(i) / (n_data - 1) * 2 - 1;
        real_t y = smolyak(x);
        real_t z = x * std::numbers::pi_v<real_t> / 4;

        Disk::write(file, std::to_string(z), " ",
                          std::to_string(y), "\n");
    }

    return 0;
}

