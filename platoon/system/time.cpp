#include "snapshot_invocable.hpp"
#include "generate_jobs.hpp"
#include <chrono>

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
            {0, 2, 1, negative_1, negative_1, negative_1};

    std::chrono::steady_clock clock;

    auto t0 = clock.now();

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
    }

    auto t1 = clock.now();

    real_t norm = smolyak_approximation_.norm();

    auto t2 = clock.now();

    Print::println(norm);
    Print::println(std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0) / 1000.0);
    Print::println(std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1) / 1000.0);
    Print::println(std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t0) / 1000.0);

    return 0;
}

