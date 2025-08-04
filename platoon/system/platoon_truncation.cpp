#include "snapshot_invocable.hpp"
#include "generate_jobs.hpp"

using namespace Platoon;
using namespace Smolyak;

int main(int argc, char* argv[])
{
    ASSERT_ASSUME(argc < 3);
    index_t n_iter = negative_1;
    if (argc == 2)
        n_iter = std::atoll(argv[1]);

    snapshot_invocable invocable_{};

    auto function = invocable_.qoi_function();
    auto cost_function = invocable_.cost_function();

    Matrix::dense_matrix<real_t> m0(3, 1), m1(3, 3), m2(3, 3);

    for (index_t i = 0; i != 3; ++i)
    {
        m0[i, 0] = 1;
        m1[i, i] = 1;
    }

    m2[1, 0] = -1;
    m2[2, 0] = -1;
    m2[1, 1] = 1;
    m2[2, 2] = 1;

    std::reference_wrapper<const Matrix::matrix_base<real_t>> r0{m0}, r1{m1}, r2{m2};

    Smolyak_traits::matrix_initialiser geometry_initialiser{std::vector{r0, r2}};
    Smolyak_traits::matrix_initialiser component_initialiser{std::vector{r1}};

    try
    {
        for (index_t i = 0; i != 3; ++i)
        {
            input_type a{i, 0, 0, 0, 0, 0};
            (void) function(a);
            (void) cost_function(a);
        }
    }
    catch (Core::missing_data_base &)
    {
        input_type a{};
        if (not exists_data(a))
            generate_jobs(std::vector{input_type{0, 0, 0, 0, 0, 0}});
        else
            generate_jobs(std::vector{input_type{0, 1, 0, 0, 0, 0}});
        return 0;
    }

    auto smolyak_approximation_ = smolyak_approximation_factory
            <matrix_traits, fidelity_traits, matrix_traits,
            interpolation_traits, interpolation_traits, interpolation_traits>::factory
        (function, cost_function, std::tuple{component_initialiser, geometry_initialiser});

    smolyak_approximation_.maximum_level() = std::array<index_t, 6>
            {0, 2, 0, negative_1, negative_1, negative_1};

    std::vector<real_t> norms, error_estimation, cost;
    try
    {
        for(index_t i = 0; i != n_iter; ++i)
        {
            auto adaptivity_information = 
                smolyak_approximation_.activate_node_adaptively(verbosity::all);

            cost.emplace_back(smolyak_approximation_.cost());
            error_estimation.emplace_back(adaptivity_information.error_estimation);
        }
    }

    catch (Core::missing_data<std::vector<input_type>> &e)
    {
        std::vector<input_type> v;

        for (input_type i : e.data())
        {
            index_t &fidelity = std::get<1>(i);

            if (fidelity != 0 or exists_data(i))
                ++fidelity;
            v.emplace_back(i);
        }
        generate_jobs(v);
    }

    real_t norm = smolyak_approximation_.norm();

    {
        std::ofstream file{"convergence.dat"};
        Disk::write(file, "error cost\n");

        for (index_t i = 0; i != std::size(cost); ++i)
            Disk::write(file, std::to_string(error_estimation[i] / norm), " ",
                              std::to_string(cost[i] / (3 * 3600)), "\n");
    }

    // Surfaces
    /*
    std::array<char, 3> components_c{'d', 's', 'l'};
    constexpr index_t n_elements = 30;

    auto normalisation_factors = snapshot_invocable::normalisation_factors();

    for (index_t i = 0; i != 3; ++i)
        for (index_t j = 0; j != 3; ++j)
        {
            std::ofstream f{std::string{"platoon_"} + components_c[i] + std::to_string(j) + ".dat"};
            ASSERT_ASSUME(f.is_open());
            
            for (index_t m = 0; m != n_elements + 1; ++m)
            {
                for (index_t n = 0; n != n_elements + 1; ++n)
                {
                    real_t x = static_cast<real_t>(m) / n_elements * 2 - 1;
                    real_t y = static_cast<real_t>(n) / n_elements * 2 - 1;

                    std::tuple<index_t, index_t, real_t, real_t, real_t> a{i, j, x, y, 0};
                    real_t z = smolyak_approximation_(a) * normalisation_factors[i][0] +
                        normalisation_factors[i][1];


                    Disk::write(f, std::to_string(x) + ' ' + std::to_string(y) + ' ' +
                            std::to_string(z) + '\n');
                }
                Disk::write(f, '\n');
            }
            f.close();
        }

    */
    // Histograms
    /*const auto &nodes{smolyak_approximation_.nodes()};

    std::array<std::array<std::vector<index_t>, 3>, 5> levels, levels_n, levels_cost;

    for (const auto &i_nodes : nodes)
    {
        if (i_nodes.second.state() == node_state::unfitted)
            continue;

        const std::array<index_t, 6> &node = i_nodes.first;
        const index_t &fidelity_level{node[1]};
        const index_t &geometry_level{node[2]};
        index_t n = geometry_level == 0 ? 1 : 2;

        for (index_t i = 0; i != 3; ++i)
            n *= sampling_type::n_new_points(node[3 + i]);

        real_t cost = 0;

        std::array<index_t, 5> begin, end;

        begin[0] = fidelity_level;
        end[0] = fidelity_level + 1;

        begin[1] = geometry_level;
        end[1] = 1 + 2 * geometry_level;

        for (index_t i = 0; i != 3; ++i)
        {
            begin[2 + i] = node[3 + i] == 0 ? 0 : sampling_type::n_points(node[3 + i] - 1);
            end[2 + i] = sampling_type::n_points(node[3 + i]);
        }

        for (const std::array<index_t, 5> &i : Utility::counter{begin, end})
        {
            std::tuple<index_t, index_t, index_t, real_t, real_t, real_t> tuple;
            std::get<0>(tuple) = 0;
            std::get<1>(tuple) = i[0];
            std::get<2>(tuple) = i[1];
            std::get<3>(tuple) = sampling_type::coordinate_from_index(i[2]);
            std::get<4>(tuple) = sampling_type::coordinate_from_index(i[3]);
            std::get<5>(tuple) = sampling_type::coordinate_from_index(i[4]);

            cost += cost_function(tuple);
        }

        cost /= 3600;

        for (index_t i = 0; i != 5; ++i)
        {
            if (node[i + 1] >= std::size(levels[i][fidelity_level]))
                levels[i][fidelity_level].resize(node[i + 1] + 1, 0);
            ++levels[i][fidelity_level][node[i + 1]];

            if (node[i + 1] >= std::size(levels_n[i][fidelity_level]))
                levels_n[i][fidelity_level].resize(node[i + 1] + 1, 0);
            levels_n[i][fidelity_level][node[i + 1]] += n;

            if (node[i + 1] >= std::size(levels_cost[i][fidelity_level]))
                levels_cost[i][fidelity_level].resize(node[i + 1] + 1, 0);
            levels_cost[i][fidelity_level][node[i + 1]] += cost;
        }
    }

    for (index_t i = 0; i != std::size(levels); ++i)
    {
        for (index_t j = 0; j != 3; ++j)
        {
            std::ofstream file{std::string{"platoon_histogram_levels_"} + std::to_string(i) +
            "_" + std::to_string(j) + ".dat"};
            Disk::write(file, "x y\n");

            for (index_t k = 0; k != std::size(levels[i][j]); ++k)
                Disk::write(file, std::to_string(k), " ", std::to_string(levels[i][j][k]), "\n");
        }
    }

    for (index_t i = 0; i != std::size(levels_cost); ++i)
    {
        for (index_t j = 0; j != 3; ++j)
        {
            std::ofstream file{std::string{"platoon_histogram_levels_cost_"} + std::to_string(i) +
                                "_" + std::to_string(j) + ".dat"};
            Disk::write(file, "x y\n");

            for (index_t k = 0; k != std::size(levels_cost[i][j]); ++k)
                Disk::write(file, std::to_string(k), " ",
                            std::to_string(levels_cost[i][j][k]), "\n");
        }
    }

    for (index_t i = 0; i != std::size(levels_n); ++i)
    {
        for (index_t j = 0; j != 3; ++j)
        {
            std::ofstream file{std::string{"platoon_histogram_levels_n_"} + std::to_string(i) +
                               "_" + std::to_string(j) + ".dat"};
            Disk::write(file, "x y\n");

            for (index_t k = 0; k != std::size(levels_n[i][j]); ++k)
                Disk::write(file, std::to_string(k), " ",
                            std::to_string(levels_n[i][j][k]), "\n");
        }
    }*/

    return 0;
}

