#include "core/core.hpp"

#include "sampling/chebyshev_lobatto.hpp"
#include "sampling/exponential_sampling.hpp"

using namespace Arithmetic;

int main(int argc, const char* argv[])
{
    using sampling_type = Sampling::exponential<Sampling::chebyshev_lobatto>;

    ASSERT_ASSUME(argc == 2);

    index_t sampling_index = std::stoull(argv[1]);

    real_t c = sampling_type::coordinate_from_index(sampling_index);

    Print::println(c * M_PI_4);

    return 0;
}
