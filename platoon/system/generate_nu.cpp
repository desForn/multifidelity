#include "fwd.hpp"

using namespace Arithmetic;

int main(int argc, char *argv[])
{
    ASSERT_ASSUME(argc == 2);
    real_t reynolds = std::atof(argv[1]);
    real_t nu = 1 / reynolds;

    std::string str = "FoamFile\n{\n\t\tversion\t\t2.0;\n\t\tformat\t\tascii;\n";
    str += "\t\tclass\t\tdictionary;\n\t\tobject\t\ttransportProperties;\n}\n\n";
    str += "transportModel\t\tNewtonian;\n";
    str += "nu\t\t\t\t\t[0 2 -1 0 0 0 0]\t\t";

    const char *c_str = str.c_str();
    std::ofstream file{"constant/transportProperties"};
    ASSERT_ASSUME(file.is_open());
    file << std::scientific;
    file.write(c_str, std::strlen(c_str));

    file << nu << ";\n";

    file.close();

    return 0;
}

