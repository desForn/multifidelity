#include "fwd.hpp"

using namespace Platoon;

int main(int argc, const char* argv[])
{
    if(argc != 3)
    {
        Print::println("in qoi: invalid arguments.");
        std::abort();
    }

    index_t fidelity = std::atoll(argv[1]);
    real_t cost = std::atof(argv[2]);

    std::ifstream file{"postProcessing/force/0/force.dat"};

    if(!file.is_open())
    {
        Print::println("in qoi: cannot open force.dat");
        std::ofstream abort{"abort"};
        std::abort();
    }

    for (std::string line, to_search{"# Time"}; std::getline(file, line); )
    {
        if (std::size(line) < 6)
            continue;

        std::string begin{std::cbegin(line), std::cbegin(line) + 6};

        if (begin == to_search)
            break;
    }

    std::vector<std::array<real_t, 3>> f;
    std::vector<real_t> t;

    for (std::string line; std::getline(file, line); )
    {
        auto split = std::ranges::views::split(line, ' ');
        auto it = std::begin(split);

        std::string word{std::string_view(*it)};
        while (std::empty(word))
            word = static_cast<std::string>(std::string_view(*++it));

        t.emplace_back(std::stod(word));
        ++it;

        std::array<real_t, 3> a;

        for (auto it_a = std::begin(a); it != std::end(split) and it_a != std::end(a); ++it)
        {
            word = static_cast<std::string>(std::string_view(*it));

            if (std::empty(word))
                continue;

            *it_a++ = std::stod(word);
        }

        f.emplace_back(a);
    }
    file.close();

    auto it_t = std::begin(t);
    while(*it_t < Settings::start_average[fidelity])
        ++it_t;
    auto it_f = std::begin(f) + std::distance(std::begin(t), it_t);

    std::array<std::array<real_t, 2>, 3> min_max;
    std::array<real_t, 3> threshold;
    real_t e = 0.2;

    for (index_t i = 0; i != 3; ++i)
    {
        min_max[i][0] = (*std::min_element(it_f, std::end(f),
            [i](const std::array<real_t, 3> &arg0, const std::array<real_t, 3> &arg1)
            {
                return arg0[i] < arg1[i];
            }))[i];

        min_max[i][1] = (*std::min_element(it_f, std::end(f),
            [i](const std::array<real_t, 3> &arg0, const std::array<real_t, 3> &arg1)
            {
                return arg0[i] > arg1[i];
            }))[i];

        threshold[i] = min_max[i][1] * e + min_max[i][0] * (1 - e);
    }

#ifdef PLOT_HISTORY
    std::array<std::vector<real_t>, 3> plot_f;
    std::array<std::vector<real_t>, 3> plot_t;
    std::array<std::vector<real_t>, 3> plot_a;
    std::array<std::vector<real_t>, 3> plot_a_t;
    real_t start;
#endif

    std::array<real_t, 3> average_f;
    std::array<real_t, 3> average_t;

    for (index_t i = 0; i != 3; ++i)
    {
        auto start_f = it_f;
        auto start_t = it_t;
        auto end_f = std::end(f);
        auto end_t = std::end(t);

        ASSERT_ASSUME(end_f > start_f and end_t > start_t);
#ifdef PLOT_HISTORY
        for (auto it = start_f; it != end_f; ++it)
            (plot_f[i]).emplace_back((*it)[i]);

        for (auto it = start_t; it != end_t; ++it)
            (plot_t[i]).emplace_back(*it);

        start = *start_t;
#endif

        average_t[i] = *(end_t - 1) - *start_t;
        average_f[i] = (*start_f)[i] * (*(start_t + 1) - *start_t) / 2;

#ifdef PLOT_HISTORY
        plot_a_t[i].emplace_back((*(start_t + 1) + *start_t) / 2);
        plot_a[i].emplace_back(average_f[i] / (plot_a_t[i].back() - start));
#endif

        for (; start_f < end_f - 1; ++start_f, ++start_t)
        {
            average_f[i] += (*start_f)[i] * (*(start_t + 1) - *(start_t - 1)) / 2;
#ifdef PLOT_HISTORY
            plot_a_t[i].emplace_back((*(start_t + 1) + *(start_t - 1)) / 2);
            plot_a[i].emplace_back(average_f[i] / (plot_a_t[i].back() - start));
#endif
        }

        average_f[i] += (*start_f)[i] * (*start_t - *(start_t - 1)) / 2;

#ifdef PLOT_HISTORY
        plot_a_t[i].emplace_back(*start_t);
        plot_a[i].emplace_back(average_f[i] / (plot_a_t[i].back() - start));
#endif

        average_f[i] /= average_t[i];
    }

    std::ofstream data{"data", std::ios::out | std::ios::binary};
    ASSERT_ASSUME(data.is_open());

    Disk::write(data, average_f);
    Disk::write(data, cost);

    Print::println(average_f, ' ', cost);

#ifdef PLOT_HISTORY
    Plot::plot_2d pl;

    for (index_t i = 0; i != 3; ++i)
    {
        pl.filename() = "temp" + std::to_string(i);
        pl.title() = "Force_" + std::to_string(i);
        pl(plot_t[i], plot_f[i]);

        std::ranges::for_each(plot_a[i], [average = average_f[i]](real_t &arg)
        {
            arg = std::log10(std::abs(arg - average));
        });

        pl.filename() = "temp" + std::to_string(i + 3);
        pl.title() = "Average_" + std::to_string(i);
        pl(plot_a_t[i], plot_a[i]);
    }
#endif
    return 0;
}
