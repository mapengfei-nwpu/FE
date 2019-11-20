#include <iostream>
#include <dolfin.h>
#include "Poisson.h"
#include "BoxAdjacents.h"
#include "IBInterpolation.h"
using namespace dolfin;

class xplusy : public Expression
{
public:
	xplusy() : Expression(2) { }
	void eval(Array<double>& values, const Array<double>& x) const
	{
		values[0] = x[0];
		values[1] = x[1];
	}
};

int main()
{
	Point p0(0, 0, 0);
	Point p1(1, 1, 0);
	BoxAdjacents ba({ p0, p1 }, { 8, 8 }, CellType::Type::quadrilateral);
	auto V = std::make_shared<Poisson::FunctionSpace>(ba.mesh());
	auto g = std::make_shared<xplusy>();
	Function u(V);
	u.interpolate(*g);
	std::vector<double> body(8);
	std::vector<std::array<double, 2>> body_coordinates(4);
	body_coordinates[0][0] = 0.1; body_coordinates[0][1] = 0.5;
	body_coordinates[1][0] = 0.2; body_coordinates[1][1] = 0.6;
	body_coordinates[2][0] = 0.3; body_coordinates[2][1] = 0.7;
	body_coordinates[3][0] = 0.4; body_coordinates[3][1] = 0.8;
	DeltaInterplation::fun1(u, ba, body, body_coordinates);
	for (size_t i = 0; i < 8; i++)
	{
		std::cout << body[i] << std::endl;
	}


	File file("xplusy.pvd");
	file << u;
}
