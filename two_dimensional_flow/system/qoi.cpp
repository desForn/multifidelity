#include "core/core.hpp"

using namespace Arithmetic;

int main(int argc, const char* argv[])
{
    if(argc != 3)
    {
        Print::println("in qoi: invalid arguments. Must be called as");
        Print::println("\t[executable name] fidelity_level sampling_index.");
        std::abort();
    }

    std::string fidelity_level{argv[1]};
    std::string sampling_index{argv[2]};

    std::string path{"../data/"};
    path += fidelity_level + "_" + sampling_index + "/postProcessing/force/0/force.dat";

    std::ifstream file{path};

    if(!file.is_open())
    {
        Print::println("in qoi: cannot open ", path);
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

    real_t value = 0;

    for (std::string line; std::getline(file, line); )
    {
        auto split = std::ranges::views::split(line, ' ');
        auto it = std::begin(split);

        ++it;

        for (; it != std::end(split); ++it)
        {
            std::string word{std::string_view(*it)};

            if (std::empty(word))
                continue;

            value = std::stod(word);
            break;
        }
    }

    Print::println(value);

    return 0;
}

