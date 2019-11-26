#include <iostream>
#include <dolfin.h>

#include "ElasticStructure.h"
using namespace dolfin;

class velocity : public Expression
{
public:
	velocity() : Expression(2) {}
	void eval(Array<double> &values, const Array<double> &x) const
	{
		values[0] = 0.0;
		values[1] = 0.0;
	}
};

int main()
{
	/// Load mesh from file.
	auto circle = std::make_shared<Mesh>("../circle.xml.gz");
	auto V = std::make_shared<ElasticStructure::FunctionSpace>(circle);

	// Create velocity term.
	auto f = std::make_shared<velocity>();

	// Define variational problem
	ElasticStructure::BilinearForm a(V, V);
	ElasticStructure::LinearForm L(V);
	L.u = f;

	// Solve PDE
	Function u(V);
	solve(a == L, u);

	// Save solution in VTK format
	File file_u("force.pvd");
	std::cout<<"yes"<<std::endl;
	file_u << u;

	return 0;
}
