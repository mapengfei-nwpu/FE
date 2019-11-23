#include <iostream>
#include <dolfin.h>
#include "Chanel.h"
#include "Circle.h"
#include "BoxAdjacents.h"
#include "IBInterpolation.h"
using namespace dolfin;

class xplusy : public Expression
{
public:
	xplusy() : Expression(2) {}
	void eval(Array<double> &values, const Array<double> &x) const
	{
		values[0] = x[0];
		values[1] = x[1];
	}
};

int main()
{
	/// Construct the mesh of the box.
	Point p0(0, 0, 0);
	Point p1(1, 1, 0);
	BoxAdjacents box({p0, p1}, {256, 256}, CellType::Type::quadrilateral);
	auto V = std::make_shared<Chanel::FunctionSpace>(box.mesh());
	auto g = std::make_shared<xplusy>();
	Function v(V);

	/// Load mesh from file.
	auto circle = std::make_shared<Mesh>("./circle.xml.gz");
	auto U = std::make_shared<Circle::FunctionSpace>(circle);
	Function u(U);
	u.interpolate(*g);

	/// print mesh size
	/// solid mesh should be finer than chanel mesh
	std::cout << "chanel size : " << box.mesh()->hmax() << std::endl
			  << "circle size : " << circle->hmax() << std::endl;

	/// Interpolate the mesh.
	DeltaInterplation interpolation(box);
	interpolation.solid_to_fluid(v, u);

	File file("fluid.pvd");
	file << v;
}
