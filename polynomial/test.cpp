#include "polynomial.hpp"

#include "../sampling/sampling.hpp"

#include "../random/random.hpp"

using Arithmetic::index_t;
using Arithmetic::real_t;
using Arithmetic::complex_t;

auto equal_to = Arithmetic::equal_to(1E-6);

using namespace Polynomial;

// **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //

void print_separator()
{
    std::cout
            << "// **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //\n";
}

// **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //

int main()
{
#if 1 // Check basic operations in univariate polynomials
    {
        {
            using polynomial_type = jacobi_polynomial<-0.2, 0.7>;

            index_t order = 4;
            polynomial_type p, q, r;

            for (index_t i = 0; i != order + 1; ++i)
            {
                p[i] = Random::real();
                q[i] = Random::real();
            }

            real_t x = Random::real();

            r = p(q);
            assert(equal_to(r(x), p(q(x))));
            Print::println("Correct cmp");

            r = p + q;
            assert(equal_to(r(x), p(x) + q(x)));
            Print::println("Correct sum");

            r = p - q;
            assert(equal_to(r(x), p(x) - q(x)));
            Print::println("Correct sub");

            r = p * q;
            assert(equal_to(r(x), p(x) * q(x)));
            Print::println("Correct mul");

            r = p ^ 3;
            assert(equal_to(r(x), std::pow(p(x), 3)));
            Print::println("Correct pow");
        }
    }
#endif

#if 1 // Check basic operations in multivariate polynomials
    {
        {
            constexpr index_t n_variates = 3;
            std::array<index_t, n_variates> extents =
                    Utility::uniform_array<index_t, n_variates>(4);


            std::array<index_t, n_variates> extents_arr =
                    Utility::uniform_array<index_t, n_variates>(4);

            using p_t = multivariate_polynomial<polynomial<real_t>, n_variates>;

            p_t p, q, r;
            std::array<p_t, n_variates> arr;
            {
                std::array<index_t, n_variates> i =
                        Utility::uniform_array<index_t, n_variates>(0);

                while (i[0] != extents[0])
                {
                    p[i] = Random::real();
                    q[i] = Random::real();

                    ++i.back();

                    for (index_t j = n_variates - 1; j != 0; --j)
                    {
                        if (i[j] == extents[j])
                        {
                            i[j] = 0;
                            ++i[j - 1];
                        }

                        else
                            break;
                    }
                }

                 i = Utility::uniform_array<index_t, n_variates>(0);

                while (i[0] != extents_arr[0])
                {
                    for (index_t j = 0; j != n_variates; ++j)
                        arr[j][i] = Random::real();

                    ++i.back();

                    for (index_t j = n_variates - 1; j != 0; --j)
                    {
                        if (i[j] == extents_arr[j])
                        {
                            i[j] = 0;
                            ++i[j - 1];
                        }

                        else
                            break;
                    }
                }
            }

            std::array<real_t, n_variates> x;
            std::array<real_t, n_variates> arr_x;

            for (index_t i = 0; i != n_variates; ++i)
                x[i] = Random::real();

            for (index_t i = 0; i != n_variates; ++i)
                arr_x[i] = arr[i](x);

            r = p(arr);

            assert(equal_to(r(x), p(arr_x)));
            Print::println("Correct cmp");

            r = p + q;
            assert(equal_to(r(x), p(x) + q(x)));
            Print::println("Correct sum");

            r = p - q;
            assert(equal_to(r(x), p(x) - q(x)));
            Print::println("Correct sub");

            r = p * q;
            assert(equal_to(r(x), p(x) * q(x)));
            Print::println("Correct mul");

            r = p ^ 3;
            assert(equal_to(r(x), std::pow(p(x), 3)));
            Print::println("Correct pow");
        }
    }
#endif

#if 1 // Check univariate legendre polynomial
    {
        Polynomial::legendre_polynomial<real_t> p, q, pq;

        static_assert(Apparatus::implemented_conversion_coefficients_c
                              <legendre_polynomial<real_t>, polynomial<real_t>> or
                              Apparatus::implemented_conversion_coefficients_c
                                      <polynomial<real_t>, legendre_polynomial<real_t>> );

        for (index_t i = 0; i <= 5; ++i)
        {
            p[i] = Random::real();
            q[i] = Random::real();
        }

        pq = p * q;

        real_t x = Random::real();

        assert(equal_to(pq(x), p(x) * q(x)));

        pq = q * p;
        assert(equal_to(pq(x), p(x) * q(x)));


        pq = p(q);
        assert(equal_to(pq(x), p(q(x))));

        static_assert(Apparatus::implemented_conversion_c
                <legendre_polynomial<real_t>, polynomial<real_t>> and
                              Apparatus::implemented_conversion_c
                                      <polynomial<real_t>, legendre_polynomial<real_t>>);

        for (index_t i = 0; i != 13; ++i)
        {
            legendre_polynomial<real_t> l = legendre_polynomial<real_t>::basis_polynomial(i);
            polynomial<real_t> p{l};
            Polynomial::Traits::polynomial_proxy_type<polynomial<real_t>> pr{p.coefficients()};
            legendre_polynomial<real_t> l_bis{pr};

            assert(std::ranges::equal(l.coefficients(), l_bis.coefficients(), equal_to));
        }

        Print::println("Successful test univariate Legendre");
    }

#endif

#if 1 // Check univariate chebyshev polynomial
    {
        Polynomial::chebyshev_polynomial_i_kind p, q, pq;

        for (index_t i = 0; i <= 5; ++i)
        {
            p[i] = Random::real();
            q[i] = Random::real();
        }

        pq = p * q;

        real_t x = Random::real();

        assert(equal_to(pq(x), p(x) * q(x)));

        pq = q * p;

        assert(equal_to(pq(x), p(x) * q(x)));

        pq = p(q);
        assert(equal_to(pq(x), p(q(x))));

        for (index_t i = 0; i != 13; ++i)
        {
            chebyshev_polynomial_i_kind l = chebyshev_polynomial_i_kind::basis_polynomial(i);
            polynomial<real_t> p{l};
            Polynomial::Traits::polynomial_proxy_type<polynomial<real_t>> pr{p.coefficients()};
            chebyshev_polynomial_i_kind l_bis{pr};

            assert(std::ranges::equal(l.coefficients(), l_bis.coefficients(), equal_to));
        }

        Print::println("Successful test univariate Chebyshev");
    }

#endif

#if 1 // Check multivariate legendre polynomial
    {
        using namespace Polynomial;

        constexpr index_t n_variates = 3;
        index_t extents = 5;

        using polynomial_type = polynomial<real_t>;
        using multivariate_polynomial_type = multivariate_polynomial<polynomial_type, n_variates>;

        using legendre_type = legendre_polynomial<real_t>;
        using multivariate_legendre_type = multivariate_polynomial<legendre_type, n_variates>;

        std::array<polynomial_type, n_variates> a;
        std::array<multivariate_polynomial_type, n_variates> a_bis;
        multivariate_polynomial_type p, q;

        for (index_t i = 0; i != n_variates; ++i)
        {
            for (index_t j = 0; j != extents; ++j)
                a[i][j] = Random::real();

            a_bis[i] = multivariate_polynomial_type{a[i], i};
        }

        for (const auto &i : Utility::counter(Utility::uniform_array<index_t, n_variates>(extents)))
        {
            p[i] = Random::real();
            q[i] = Random::real();
        }

        multivariate_polynomial_type pa{p(a_bis)};

        std::array<real_t, n_variates> x = Random::real_array<n_variates>();
        std::array<real_t, n_variates> y;

        for (index_t i = 0; i != n_variates; ++i)
            y[i] = a[i](x[i]);

        assert(equal_to(pa(x), p(y)));

        std::array<legendre_type , n_variates> a_l;
        std::array<multivariate_legendre_type, n_variates> a_l_bis;
        multivariate_legendre_type p_l{p}, q_l{q};

        for (index_t i = 0; i != n_variates; ++i)
        {
            a_l[i] = legendre_type{a[i]};
            a_l_bis[i] = multivariate_legendre_type{a_l[i], i};
        }

        multivariate_legendre_type pa_l = p_l(a_l_bis);

        assert(equal_to(p(x), p_l(x)));
        assert(equal_to(q(x), q_l(x)));

        p = static_cast<multivariate_polynomial_type>(p_l);

        assert(equal_to(p(x), p_l(x)));

        assert(equal_to(p(x) * q(x), (p * q)(x)));
        assert(equal_to(p(x) * q(x), (p_l * q_l)(x)));
        assert(equal_to(pa_l(x), pa(x)));

        Print::println("Successful test multivariate Legendre");
    }
#endif

#if 1 // Check multivariate chebyshev polynomial
    {
        using namespace Polynomial;

        constexpr index_t n_variates = 3;
        index_t extents = 5;

        using polynomial_type = polynomial<real_t>;
        using multivariate_polynomial_type = multivariate_polynomial<polynomial_type, n_variates>;

        using chebyshev_type = chebyshev_polynomial_i_kind;
        using multivariate_chebyshev_type = multivariate_polynomial<chebyshev_type , n_variates>;

        std::array<polynomial_type, n_variates> a;
        std::array<multivariate_polynomial_type, n_variates> a_bis;
        multivariate_polynomial_type p, q;

        for (index_t i = 0; i != n_variates; ++i)
        {
            for (index_t j = 0; j != extents; ++j)
                a[i][j] = Random::real();

            a_bis[i] = multivariate_polynomial_type{a[i], i};
        }

        for (const auto &i : Utility::counter(Utility::uniform_array<index_t, n_variates>(extents)))
        {
            p[i] = Random::real();
            q[i] = Random::real();
        }

        multivariate_polynomial_type pa{p(a_bis)};

        std::array<real_t, n_variates> x = Random::real_array<n_variates>();
        std::array<real_t, n_variates> y;

        for (index_t i = 0; i != n_variates; ++i)
            y[i] = a[i](x[i]);

        assert(equal_to(pa(x), p(y)));

        std::array<chebyshev_type , n_variates> a_l;
        std::array<multivariate_chebyshev_type , n_variates> a_l_bis;
        multivariate_chebyshev_type p_l{p}, q_l{q};

        for (index_t i = 0; i != n_variates; ++i)
        {
            a_l[i] = chebyshev_type{a[i]};
            a_l_bis[i] = multivariate_chebyshev_type{a_l[i], i};
        }

        multivariate_chebyshev_type pa_l = p_l(a_l_bis);

        assert(equal_to(p(x), p_l(x)));
        assert(equal_to(q(x), q_l(x)));

        p = static_cast<multivariate_polynomial_type>(p_l);

        assert(equal_to(p(x), p_l(x)));

        assert(equal_to(p(x) * q(x), (p * q)(x)));
        assert(equal_to(p(x) * q(x), (p_l * q_l)(x)));
        assert(equal_to(pa_l(x), pa(x)));

        Print::println("Successful test multivariate Chebyshev");
    }
#endif

#if 1 // Check multivariate jacobi polynomial
    {
        using namespace Polynomial;

        constexpr index_t n_variates = 3;
        index_t extents = 5;

        using polynomial_type = polynomial<real_t>;
        using multivariate_polynomial_type = multivariate_polynomial<polynomial_type, n_variates>;

        using jacobi_type = jacobi_polynomial<0.3, -0.5>;
        using multivariate_jacobi_type = multivariate_polynomial<jacobi_type , n_variates>;

        std::array<polynomial_type, n_variates> a;
        std::array<multivariate_polynomial_type, n_variates> a_bis;
        multivariate_polynomial_type p, q;

        for (index_t i = 0; i != n_variates; ++i)
        {
            for (index_t j = 0; j != extents; ++j)
                a[i][j] = Random::real();

            a_bis[i] = multivariate_polynomial_type{a[i], i};
        }

        for (const auto &i : Utility::counter(Utility::uniform_array<index_t, n_variates>(extents)))
        {
            p[i] = Random::real();
            q[i] = Random::real();
        }

        multivariate_polynomial_type pa{p(a_bis)};

        std::array<real_t, n_variates> x = Random::real_array<n_variates>();
        std::array<real_t, n_variates> y;

        for (index_t i = 0; i != n_variates; ++i)
            y[i] = a[i](x[i]);

        assert(equal_to(pa(x), p(y)));

        std::array<jacobi_type , n_variates> a_l;
        std::array<multivariate_jacobi_type , n_variates> a_l_bis;
        multivariate_jacobi_type p_l{p}, q_l{q};

        for (index_t i = 0; i != n_variates; ++i)
        {
            a_l[i] = jacobi_type{a[i]};
            a_l_bis[i] = multivariate_jacobi_type{a_l[i], i};
        }

        multivariate_jacobi_type pa_l = p_l(a_l_bis);

        assert(equal_to(p(x), p_l(x)));
        assert(equal_to(q(x), q_l(x)));

        p = static_cast<multivariate_polynomial_type>(p_l);

        assert(equal_to(p(x), p_l(x)));

        assert(equal_to(p(x) * q(x), (p * q)(x)));
        assert(equal_to(p(x) * q(x), (p_l * q_l)(x)));
        assert(equal_to(pa_l(x), pa(x)));

        Print::println("Successful test multivariate Jacobi");
    }
#endif

#if 1 // Check multivariate chebyshev_point_value polynomial
    {
        constexpr index_t n_variates = 4;
        index_t sampling_level = 4;

        using polynomial_type = Polynomial::polynomial<real_t>;
        using polynomial_point_value_type = chebyshev_polynomial_point_value;
        using sampling_type = polynomial_point_value_type::sampling_type;

        index_t extents = sampling_type::n_points(sampling_level);

        using multivariate_polynomial_type =
                Polynomial::multivariate_polynomial<polynomial_type, n_variates>;

        using multivariate_polynomial_point_value_type =
                Polynomial::multivariate_polynomial<polynomial_point_value_type, n_variates>;

        polynomial_type p;

        for (index_t i = 0; i != extents; ++i)
            p[i] = Random::real();

        polynomial_point_value_type pv{p, sampling_level, p.order()};
        polynomial_point_value_type pb{pv};

        pb -= pb;
        pb += pv;

        real_t x = Random::real();
        complex_t xc{x};

        assert(equal_to(p(x), pv(x)));
        assert(equal_to(p(x), pb(x)));

        assert(equal_to(p(x), p(xc)));
        assert(equal_to(p(x), pv(xc)));
        assert(equal_to(p(x), pb(xc)));

        assert(equal_to(p.order(), pv.order()));
        assert(equal_to(p.order(), pb.order()));

        assert(equal_to(sampling_level, pv.sampling_level()));
        assert(equal_to(sampling_level, pb.sampling_level()));

        multivariate_polynomial_type mp;

        for (const auto &i : Utility::counter(Utility::uniform_array<index_t, n_variates>(extents)))
            mp[i] = Random::real();

        multivariate_polynomial_point_value_type mpv{mp, Utility::uniform_array<index_t, n_variates>(sampling_level)};

        auto mx = Random::real_array<n_variates>();
        auto mxc = Utility::static_cast_array<complex_t>(mx);

        assert(equal_to(mp(mx), mp(mxc)));
        assert(equal_to(mp(mx), mpv(mx)));
        assert(equal_to(mp(mx), mpv(mxc)));
        assert(equal_to(mp(mx), mpv(mx[0], mxc[1], mx[2], mxc[3])));

        Print::println("Successful test multivariate point-value Chebyshev");
    }
#endif

#if 1 // Check quadrature

    {
        using polynomial_type = Polynomial::chebyshev_polynomial_i_kind;
        using polynomial_point_value_type = Polynomial::chebyshev_polynomial_point_value;

        auto func = [](real_t x) -> real_t
        {
            return std::log(x + 4);
        };

        const index_t n = 3;

        polynomial_point_value_type p(func, n);
        polynomial_type q = static_cast<polynomial_type>(p);

        real_t x = 0.637;

        assert(equal_to(p(x), q(x)));

        p = static_cast<polynomial_point_value_type>(q);
        assert(equal_to(p(x), q(x)));
        assert(equal_to(Polynomial::integrate(p), Polynomial::integrate(q)));

        Print::println("Successful test univariate quadrature");
    }

#endif

#if 1 // Check multivariate quadrature

    {
        constexpr index_t n_variates = 3;

        using polynomial_type = Polynomial::multivariate_polynomial
                <Polynomial::chebyshev_polynomial_i_kind, n_variates>;
        using polynomial_point_value_type = Polynomial::multivariate_polynomial
                <Polynomial::chebyshev_polynomial_point_value, n_variates>;

        auto func = [](const std::array<real_t, n_variates> &x) -> real_t
        {
            real_t ret = n_variates + 1;

            for (index_t i = 0; i != n_variates; ++i)
                ret += x[i];

            ret = std::log(ret);

            return std::sqrt(ret);
        };

        std::array<real_t, n_variates> x;

        x.front() = 0.637;

        for (index_t i = 1; i != n_variates; ++i)
            x[i] = x[i - 1] / 1;

        for (index_t n = 0; n != 4; ++n)
        {
            std::array<index_t, n_variates>
                    sampling_level{Utility::uniform_array<index_t, n_variates>(n)};

            polynomial_point_value_type p(func, sampling_level);
            polynomial_type q = static_cast<polynomial_type>(p);

            p *= p;
            q *= q;

            assert(equal_to(p(x), q(x)));

            p = static_cast<polynomial_point_value_type>(q);
            assert(equal_to(p(x), q(x)));
            assert(equal_to(Polynomial::integrate(p), Polynomial::integrate(q)));

            Print::println("Successful test multivariate quadrature");
        }
    }

#endif

#if 1 // Check fourier transform

    {
        const index_t n = 8;

        for (index_t i = 0; i != n ; ++i)
        {
            std::vector<complex_t> f(n);
            real_t alpha = 2 * M_PI / n;

            for (index_t j = 0; j != n; ++j)
                f[j] = complex_t{std::cos(i * j * alpha), std::sin(i * j * alpha)};

            Polynomial::discrete_fourier_transform(std::span{f});

            for (index_t j = 0; j != n; ++j)
                assert(equal_to(f[j], static_cast<real_t>(i == j)));
        }

        Print::println("Successful test univariate FFT");
    }

#endif

#if 1 // Check multivariate fourier transform

    {
        constexpr index_t n_variates = 3;
        constexpr std::array<index_t, n_variates> n{8, 4, 8};
        constexpr index_t s = std::accumulate(std::begin(n), std::end(n), 1, std::multiplies{});

        index_t c = 0;
        for (const auto &i : Utility::counter(n))
        {
            std::vector<complex_t> f;

            for (const auto &j : Utility::counter(n))
            {
                real_t angle = 0;

                for (index_t k = 0; k != n_variates; ++k)
                    angle += i[k] * j[k] / static_cast<real_t>(n[k]);

                angle *= 2 * M_PI;

                f.emplace_back(complex_t{std::cos(angle), std::sin(angle)});
            }

            Polynomial::multidimensional_discrete_fourier_transform(f, n);

            for (index_t j = 0; j != s; ++j)
                assert(equal_to(f[j], static_cast<real_t>(c == j)));

            ++c;
        }

        Print::println("Successful test multivariate FFT");
    }

#endif
}

