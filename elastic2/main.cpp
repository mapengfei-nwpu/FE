#include <iostream>
#include <dolfin.h>

#include "ElasticStructure.h"
using namespace dolfin;

class Source : public Expression
{
public:
	Source() : Expression(2) {}
	void eval(Array<double> &values, const Array<double> &x) const
	{
		values[0] = 0.0;
		values[1] = -0.2;
	}
};

class Traction : public Expression
{
public:
	Traction() : Expression(2) {}
	void eval(Array<double> &values, const Array<double> &x) const
	{
		values[0] = 0.0;
		values[1] = 10.0;
	}
};

// Define noslip domain
class FixedBoundary : public SubDomain
{
	bool inside(const Array<double> &x, bool on_boundary) const
	{
		return near(x[0], 0);
	}
};

int main()
{
	/// Load mesh from file.
	auto LL = 1.0;
	auto WW = 0.2;
	auto _mesh = RectangleMesh::create({Point(0.0, 0), Point(LL, WW)}, {100, 30}, CellType::Type::triangle);
	auto mesh = std::make_shared<Mesh>(/*"../circle.xml.gz"); */_mesh);
	auto V = std::make_shared<ElasticStructure::FunctionSpace>(mesh);

	// Create velocity term.
	auto f = std::make_shared<Constant>(0.0, 0.0);
	auto p_f = std::make_shared<Constant>(0.0);
	auto u_f = std::make_shared<Traction>();

	/* no boundary condition*/
	auto fixed_boundary = std::make_shared<FixedBoundary>();
	auto u_fixed = std::make_shared<Constant>(0.0, 0.0);
	DirichletBC fixed(V, u_fixed, fixed_boundary);
	DirichletBC bc(V, u_fixed, fixed_boundary);
	/**/

	// Define variational problem
	ElasticStructure::BilinearForm a(V, V);
	ElasticStructure::LinearForm L(V);
	L.f = f;
	L.p_f = p_f;
	L.u_f = u_f;

	// Solve PDE

	Function u(V);
	solve(a == L, u, bc);

	ALE::move(*mesh, u);
	// Save solution in VTK format
	File file_u("force.pvd");
	file_u << u;

	return 0;
}
