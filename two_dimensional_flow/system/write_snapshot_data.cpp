#include "core/core.hpp"

using namespace Arithmetic;

int main(int argc, char *argv[])
{
    ASSERT_ASSUME(argc == 5);

    index_t fidelity = std::atoll(argv[1]);
    index_t angle_id = std::atoll(argv[2]);

    real_t qoi = std::atof(argv[3]);
    real_t cost = std::atof(argv[4]);

    std::string path{"../data/"};
    path += std::to_string(fidelity) + '_' + std::to_string(angle_id) + "/data";

    std::ofstream file{path};

    ASSERT_ASSUME(file.is_open());

    Disk::write(file, qoi);
    Disk::write(file, cost);

    return 0;
}

