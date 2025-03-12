#include "fwd.hpp"

using namespace Arithmetic;

int main(int argc, char *argv[])
{
    ASSERT_ASSUME(argc == 2);
    real_t yaw = std::atof(argv[1]);

    std::string str = "FoamFile\n{\n\t\tversion\t\t2.0;\n\t\tformat\t\tascii;\n";
    str += "\t\tclass\t\tvolVectorField;\n\t\tobject\t\tU;\n}\n\n";
    str += "dimensions\t\t\t[0 1 -1 0 0 0 0];\n";
    str += "internalField\t\tuniform (" + std::to_string(std::cos(yaw)) + " " +
            std::to_string(std::sin(yaw)) + " 0);\n\n";
    str += "boundaryField\n{\n";
    str += "\t\"(inlet|left)\"\n\t{\n\t\ttype\tfixedValue;\n\t\tvalue\t$internalField;\n\t}\n\n";
    str += "\t\"(outlet|right)\"\n\t{\n\t\ttype\tzeroGradient;\n\t}\n\n";
    str += "\t\"(top|ground)\"\n\t{\n\t\ttype\tsymmetry;\n\t}\n\n";
    str += "\t\"(body|bodyRef|stilts|stiltsRef)\"\n\t{\n\t\ttype\tfixedValue;\n";
    str += "\t\tvalue\tuniform (0 0 0);\n}\t\n\n";
    str += "\t\"proc.*\"\n\t{\n\t\ttype\tprocessor;\n\t}\n}\n\n";

    const char *c_str = str.c_str();
    std::ofstream file{"0/U"};
    ASSERT_ASSUME(file.is_open());
    file.write(c_str, std::strlen(c_str));
    file.close();

    return 0;
}

