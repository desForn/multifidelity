#include "argList.H"
#include "IOmanip.H"
#include <fstream>
#include <iomanip>
#include "dynamicFvMesh.H"

using namespace Foam;

class myDynamicFvMesh : public dynamicFvMesh
{
	public:
		myDynamicFvMesh(const IOobject& io) : dynamicFvMesh(io)
		{
		};


		bool update() {return true;}
};

point rotate(point p, double angle)
{
    double x = p[0] * std::cos(angle) + p[1] * std::sin(angle);
    p[1] = -p[0] * std::sin(angle) + p[1] * std::cos(angle);
    p[0] = x;
    return p;
}

const double s0 = 9.5;
const double x0 = 3.875;
const double x1 = 7.625;
const double x2 = 9;
const double x3 = 9.125;
double s, dx, dy, dz;

point move(point p)
{
    if (s > s0)
    {
        if (p[0] >= x2)
            p[0] += s - s0;

        else if (p[0] > x1)
            p[0] += (s - s0) * (p[0] - x1) / (x2 - x1);
    }

    else
    {
        if (p[0] >= x3)
            p[0] += s - s0;

        else if (p[0] > x0)
            p[0] += (s - s0) * (p[0] - x0) / (x3 - x0);
    }

    return p;
}

int main(int argc, char *argv[])
{
    argList::addOption("s", "s", "spacing");
    argList::addOption("g", "g", "geometry");

    #include "setRootCase.H"
    #include "createTime.H"

    if (not (args.found("s") and args.found("g")))
    {
        std::cout << "Invalids arguments" << std::endl;
        return 0;
    }

    s = std::stod(args["s"]);
    int g = std::stoi(args["g"]);

    if (g == 1)
        s -= 0.0097;

    else if (g == 2)
        s -= 0.0556;

    else if (g == 3)
        s += 0.01389;

    else if (g != 0)
    {
        std::cout << "Invalid geometry" << std::endl;
		return 0;
    }

	autoPtr<myDynamicFvMesh> meshPtr
        (
                        new myDynamicFvMesh
                        (
                                IOobject
                                (
                                        myDynamicFvMesh::defaultRegion,
                                        runTime.timeName(),
                                        runTime,
                                        IOobject::MUST_READ,
					IOobject::AUTO_WRITE
                                )
                        )
	);

	myDynamicFvMesh& mesh = meshPtr();

	pointField points(mesh.nPoints(), point(0, 0, 0));

	forAll(mesh.points(), iPoint)
        points[iPoint] = move(mesh.points()[iPoint]);
	
	mesh.movePoints(points);

	mesh.write();

	std::cout << "End" << std::endl;

	return 0;
}
