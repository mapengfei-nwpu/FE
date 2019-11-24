#include <iostream>
#include <dolfin.h>
#include <dolfin/geometry/SimplexQuadrature.h>

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
		values[0] = 1.0;
		values[1] = 1.0;
	}
};

int main()
{
	/// Construct the mesh of the box.
	Point p0(0, 0, 0);
	Point p1(0.4, 0.4, 0);
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
	// u.interpolate(*g);

	File file("fluid.pvd");
	file << v;

	SimplexQuadrature gq(2, 9);
	Point p2(0.4, 0.0, 0);
	std::vector<Point> tri = {p0, p1, p2};
	auto gsqr = gq.compute_quadrature_rule(tri, 2);
	for (size_t i = 0; i < gsqr.second.size(); i++)
	{
		std::cout << gsqr.second[i] << std::endl;
	}
}
