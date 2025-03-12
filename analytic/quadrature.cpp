#include "fwd.hpp"
#include "random/random.hpp"
#include "invocable/invocable_function_hash_table.hpp"

#ifndef N_VARIATES
    #error "Must define N_VARIATES"
#endif

using namespace Analytic;
using namespace Smolyak;

namespace Apparatus
{
    struct smolyak_approximation_truncate_struct
    {
        template<class... types>
        static smolyak_approximation_factory<fidelity_traits, types...>
                function(Utility::homotype_pack<types...>);
    };

    struct smolyak_approximation_richardson_struct
    {
        template<class... types>
        static smolyak_approximation_factory<richardson_traits, types...>
        function(Utility::homotype_pack<types...>);
    };
}

constexpr index_t n_variates = N_VARIATES;

using array_type = std::array<real_t, n_variates>;
using tuple_type = Utility::uniform_tuple_type<real_t, n_variates>;
using tuple_extnended_type = decltype(std::tuple_cat(std::tuple<index_t>{}, tuple_type{}));

using smolyak_approximation_factory_truncate_type = decltype(
    ::Apparatus::smolyak_approximation_truncate_struct::function(
        Utility::generate_homotype_pack<interpolation_traits, n_variates>{}));
using smolyak_approximation_factory_richardson_type = decltype(
    ::Apparatus::smolyak_approximation_richardson_struct::function(
        Utility::generate_homotype_pack<interpolation_traits, n_variates>{}));

[[noreturn]] void exit()
{
    Print::println("Invalid arguments.");
    Print::println("\t[executable name] n_iter n_test_points quadrature_level fidelity");
    Print::println("fidelity:");
    Print::println("\tclen");
    Print::println("\trich");
    Print::println("\ttrap");

    std::abort();
}

int main(int argc, char *argv[])
{
    if (argc != 5)
        exit();

    index_t n = std::atoll(argv[1]);
    index_t m = std::atoll(argv[2]);
    index_t l = std::atoll(argv[3]);
    std::string fidelity_string{argv[4]};

    auto l_array = Utility::uniform_array<index_t, n_variates>(l);

    enum fidelity_t {clenshaw_curtis, richardson, trapezoidal};
    fidelity_t fidelity;

    if (fidelity_string == "clen")
        fidelity = fidelity_t::clenshaw_curtis;

    else if (fidelity_string == "rich")
        fidelity = fidelity_t::richardson;

    else if (fidelity_string == "trap")
        fidelity = fidelity_t::trapezoidal;

    else
        exit();

    constexpr real_t epsilon = 0.03;

    // integral(sqrt(1-x*x) dx, x in (-(1 - epsilon), 1 - epsilon)) / (1 - epsilon) =
    // sqrt(1-(1-epsilon)^2) * (1 - epsilon) + asin(1 - epsilon) / (1 - epsilon)
    const real_t univariate_quadrature = (std::sqrt(1 - std::pow(1 - epsilon, 2)) * (1 - epsilon) +
            std::asin(1 - epsilon)) / (1 - epsilon);

    auto function = [epsilon, fidelity](const tuple_extnended_type &arg) -> real_t
    {
        array_type arg_f = Utility::transform_array
                ([epsilon](real_t arg_bis){ return arg_bis * (1 - epsilon); },
                 Utility::tuple_to_array(Utility::tail_tuple<0>(arg)));

        real_t x = f(arg_f);

        auto integrand = [x, epsilon](real_t arg) -> real_t
        {
            arg *= (1 - epsilon);
            return x * f(std::array<real_t, 1>{arg});
        };

        if (fidelity == fidelity_t::clenshaw_curtis)
            return Polynomial::integrate
                    (Polynomial::polynomial_point_value<sampling_type>{integrand, std::get<0>(arg)});

        Sampling::exponential<Sampling::equispaced> sampling(std::get<0>(arg));
        std::vector<real_t> y = Invocable::evaluate(integrand, sampling.coordinates());

        if (std::size(y) == 1)
            return y.front() * 2;

        y.front() /= 2;
        y.back() /= 2;
        real_t ret = std::accumulate(std::cbegin(y), std::cend(y), static_cast<real_t>(0));
        ret *= 2 / static_cast<real_t>(std::size(y) - 1);
        return ret;
    };

    auto exact_quadrature = [epsilon, univariate_quadrature](const tuple_type &arg) -> real_t
    {
        array_type arg_f = Utility::transform_array
                ([epsilon](real_t arg_bis){ return arg_bis * (1 - epsilon); },
                 Utility::tuple_to_array(arg));

        return f(arg_f) * univariate_quadrature;
    };

    auto exact_quadrature_2 = [&exact_quadrature](const array_type &arg) -> real_t
    {
        real_t ret = exact_quadrature(Utility::array_to_tuple(arg));
        return ret * ret;
    };

    auto cost_function = [](const tuple_extnended_type &arg) -> real_t
    {
        return Sampling::n_points_exponential_sampling(std::get<0>(arg));
    };

    auto smolyak_approximation_truncate = smolyak_approximation_factory_truncate_type::
        factory(function, cost_function);

    auto smolyak_approximation_richardson = smolyak_approximation_factory_richardson_type::
        factory(function, cost_function);

    std::vector<tuple_type> test_points;

    Random::set_seed(0);
    for (index_t i = 0; i != m; ++i)
        test_points.emplace_back(Utility::array_to_tuple(Random::real_array<n_variates>(-1, 1)));

    real_t norm_u_f = 0;
    for (const tuple_type &j: test_points)
        norm_u_f = std::max(norm_u_f, std::abs(exact_quadrature(j)));

    Polynomial::multivariate_polynomial
        <Polynomial::chebyshev_polynomial_point_value, n_variates>
            p_ref{exact_quadrature_2, l_array};

    real_t norml2_f = std::sqrt(Polynomial::integrate(p_ref));
    Print::println(norml2_f);
    Print::println();

    std::vector<real_t> cost, error_u, error_l2, error_estimation;

    for (index_t i = 0; i != n; ++i)
    {
        real_t w, estimation;

        if (fidelity == fidelity_t::richardson)
        {
            auto info = smolyak_approximation_richardson.activate_node_adaptively
                                                    (Smolyak::verbosity::all);
            w = smolyak_approximation_richardson.cost();
            estimation = info.error_estimation;
        }

        else
        {
            auto info = smolyak_approximation_truncate.activate_node_adaptively
                                                  (Smolyak::verbosity::all);
            w = smolyak_approximation_truncate.cost();
            estimation = info.error_estimation;
        }

        cost.emplace_back(w);
        error_estimation.emplace_back(estimation);

        real_t e = 0;
        for (const tuple_type &j: test_points)
        {
            if (fidelity == fidelity_t::richardson)
                e = std::max(e,
                             std::abs(smolyak_approximation_richardson(j) - exact_quadrature(j)));
            else
                e = std::max(e,
                             std::abs(smolyak_approximation_truncate(j) - exact_quadrature(j)));
        }
        e /= norm_u_f;

        error_u.emplace_back(e);

        auto error_invocable =
            [fidelity, &exact_quadrature, &smolyak_approximation_richardson,
            &smolyak_approximation_truncate](const array_type &arg) -> real_t
                {
                    auto tuple = Utility::array_to_tuple(arg);
                    real_t ret;
                    
                    if (fidelity == fidelity_t::richardson)
                        ret = smolyak_approximation_richardson(tuple);
                    else
                        ret = smolyak_approximation_truncate(tuple);

                    ret -= exact_quadrature(tuple);

                    return ret * ret;
                };

        Polynomial::multivariate_polynomial
            <Polynomial::chebyshev_polynomial_point_value, n_variates>
                p{error_invocable, l_array};

        error_l2.emplace_back(std::sqrt(Polynomial::integrate(p)));

        std::ofstream data{std::string{"data_quadrature_"} + std::to_string(n_variates + 1) +
                           "_" + fidelity_string};
        ASSERT_ASSUME(data.is_open());
        Disk::write(data, "cost error_u error_l2 error_estimation\n");
        for (index_t i = 0; i != std::size(cost); ++i)
            Disk::write(data, std::to_string(cost[i]) + ' ' + std::to_string(error_u[i]) + ' ' +
                              std::to_string(error_l2[i]) + ' ' +
                              std::to_string(error_estimation[i]) + '\n');
            
        data.close();
    }

    return 0;
}
