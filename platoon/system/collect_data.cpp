#include "snapshot_invocable.hpp"

using namespace Platoon;

std::array<index_t, 5> increase_fidelity_level(std::array<index_t, 5> arg)
{
    ++std::get<0>(arg);
    return arg;
}

int main()
{
    std::vector<std::pair<std::tuple<index_t, index_t, real_t, real_t, real_t>,
                          std::array<real_t, 4>>> v;

    std::array<index_t, 5> l{}, h{4, 4, negative_1, negative_1, negative_1};
    std::array<index_t, 5> i{l};

    if (not exists_data(increase_fidelity_level(i)))
        goto end_loop;

    begin_loop:
    {
        std::ifstream data{file_name_data(increase_fidelity_level(i))};
        ASSERT_ASSUME(data.is_open());

        std::array<real_t, 4> x;
        Disk::read(data, x);
        data.close();

        Print::println(i, x);

        std::tuple<index_t, index_t, real_t, real_t, real_t> input
                {i[0],
                 i[1],
                 sampling_type::coordinate_from_index(i[2]),
                 sampling_type::coordinate_from_index(i[3]),
                 sampling_type::coordinate_from_index(i[4])};

        v.emplace_back(input, x);

        ++i.back();

        for (index_t j = 4; j != 0; --j)
        {
            if (std::ranges::equal(i, h, std::less{}) and exists_data(increase_fidelity_level(i)))
                goto begin_loop;

            i[j] = 0;
            ++i[j - 1];
        }

        if (std::ranges::equal(i, h, std::less{}) and exists_data(increase_fidelity_level(i)))
            goto begin_loop;
    }
    end_loop:

    std::ofstream file{"data", std::ios::out | std::ios::binary};

    ASSERT_ASSUME(file.is_open());

    Print::println(std::size(v));
    Disk::write(file, std::size(v));

    for (const auto &i : v)
        Disk::write(file, i);

    file.close();

    return 0;
}
