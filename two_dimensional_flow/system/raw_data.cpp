#include "snapshot_invocable.hpp"

using namespace Two_dimensional_flow;

using sampling_type = Sampling::exponential<Sampling::chebyshev_lobatto>;

// **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //

void error()
{
    Print::println("Call this program as:");
    Print::println("\t [executable name] [path to data file] fidelity_level_to_show");

    std::abort();
}

int main(int argc, char* argv[])
{
    if (argc != 3)
        error();

    std::string path{argv[1]};
    index_t fidelity_level = std::atoll(argv[2]);

    snapshot_invocable snapshot_invocable_{path};

    auto qoi_function = [&snapshot_invocable_, fidelity_level]
        (real_t arg)
        {
            return snapshot_invocable_.qoi_function()
                (std::array<index_t, 1>{fidelity_level},
                 std::array<real_t, 1>{arg});
        };

    std::vector<real_t> x, y;
    std::ofstream file{"raw_data.dat"};
    Disk::write(file, "angle force\n");

    try
    {
        for (index_t i = 0; true; ++i)
        {
            real_t a = sampling_type::coordinate_from_index(i);
            real_t b = qoi_function(a); // Will throw
            x.emplace_back(a);
            y.emplace_back(b);
            Disk::write(file, std::to_string(x.back()), " ");
            Disk::write(file, std::to_string(y.back()), "\n");
        }
    }

    catch (Core::missing_data_base &)
    {
        file.close();
    }

    return 0;
}
