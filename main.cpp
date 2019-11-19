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
		values[0] = x[0] + x[1];
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

	File file("xplusy.pvd");
	file << u;

	// BoxAdjacents ba({p0,p1}, {8,8,8}, CellType::Type::hexahedron);
	// ba.check();
}
