#include <iostream>
#include <dolfin.h>
#include "Poisson.h"
#include "BoxAdjacents.h"
#include "IBInterpolation.h"
using namespace dolfin;

class xplusy : public Expression
{
public:
	xplusy() : Expression(2) {}
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
	BoxAdjacents ba({ p0, p1 }, { 16, 16 }, CellType::Type::quadrilateral);
	auto V = std::make_shared<Poisson::FunctionSpace>(ba.mesh());
	auto g = std::make_shared<xplusy>();
	Function v(V);
	v.interpolate(*g);

	Point p2(0.25, 0.25, 0);
	Point p3(0.75, 0.75, 0);
	BoxAdjacents bb({ p2, p3 }, { 16, 16 }, CellType::Type::quadrilateral);
	auto U = std::make_shared<Poisson::FunctionSpace>(bb.mesh());
	Function u(U);

	DeltaInterplation di(ba);
	di.fluid_to_solid(v, u);

	File file("xplusy.pvd");
	file << u;
}
