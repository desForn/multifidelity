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

double angle(double r, double angle_max)
{
    if (r < 2)
        return angle_max;
        
    r = (16 - r) / 14;
    
    return (3 * std::pow(r, 2) - 2 * std::pow(r, 3)) * angle_max;
}

int main(int argc, char *argv[])
{
	argList::addOption("a", "a", "angle of rotation");

	#include "setRootCase.H"
	#include "createTime.H"

	if (!(args.found("a")))
	{
		std::cout << "The option a must be given" << std::endl;
		return 0;
	}

	double a = std::stod(args["a"]);
	
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
	{
	    double r = std::sqrt(std::pow(mesh.points()[iPoint][0], 2) + std::pow(mesh.points()[iPoint][1], 2));
		points[iPoint] = rotate(mesh.points()[iPoint], angle(r, a));
	}
	
	mesh.movePoints(points);

	mesh.write();

	std::cout << "End" << std::endl;

	return 0;
}

