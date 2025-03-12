#include "snapshot_invocable.hpp"

using namespace Two_dimensional_flow;

std::string path;

bool exists(const std::array<index_t, 2> &fi_an)
{
    return std::filesystem::exists(path + std::string{"data/"} + std::to_string(fi_an[0]) + "_" +
                                   std::to_string(fi_an[1]) + "/data");
}

int main(int argc, char* argv[])
{
    ASSERT_ASSUME(argc == 2);

    path = argv[1];

    using sampling_type = Sampling::exponential<Sampling::chebyshev_lobatto>;

    std::vector<snapshot_data<sampling_type>> v;

    for (std::array<index_t, 2> fi_an{0, 0}; exists(fi_an);
         fi_an = exists(Utility::increment(fi_an, 1)) ? fi_an :
                 std::array<index_t, 2>{fi_an[0] + 1, 0})
    {
        std::ifstream data{path + std::string{"data/"} + std::to_string(fi_an[0]) + "_" +
                          std::to_string(fi_an[1]) + "/data", std::ios::in | std::ios::binary};

        ASSERT_ASSUME(data.is_open());

        real_t qoi;
        real_t cost;

        Disk::read(data, qoi);
        Disk::read(data, cost);
        data.close();

        v.emplace_back(fi_an[0], sampling_type::coordinate_from_index(fi_an[1]), qoi, cost);
    }

    std::ofstream file{path + "system/data", std::ios::out | std::ios::binary};

    ASSERT_ASSUME(file.is_open());

    Disk::write(file, std::size(v));

    for(const auto &i : v)
        Disk::write(file, i);

    file.close();

    return 0;
}

