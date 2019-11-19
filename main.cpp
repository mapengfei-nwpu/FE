#include <iostream>
#include <dolfin.h>
#include "Poisson.h"
#include "BoxAdjacents.h"
#include "IBInterpolation.h"
using namespace dolfin;

class xplusy : public Expression
{
	void eval(Array<double>& values, const Array<double>& x) const
	{
		values[0] = 2.0 + x[0] + x[1];
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
	std::vector<double> body(4);
	std::vector<std::array<double, 2>> body_coordinates(4);
	body_coordinates[0][0] = 0.25; body_coordinates[0][1] = 0.25;
	body_coordinates[1][0] = 0.25; body_coordinates[1][1] = 0.75;
	body_coordinates[2][0] = 0.75; body_coordinates[2][1] = 0.25;
	body_coordinates[3][0] = 0.75; body_coordinates[3][1] = 0.75;
	DeltaInterplation::fun1(u, ba, body, body_coordinates);


	File file("xplusy.pvd");
	file << u;
}
