#include <iostream>
#include "MeshElementsCollection.h"
#include "GaussQuadrature.h"
int main() {
	/********** test 1 ************
	MeshElementsCollection a;
	a.print();
	*******************************/
	///* test 2
	MeshElementsCollection a;
	auto b = a.get_mesh_element(3);
	a.print();
	b.print();
	GaussQuadrature c;
	auto d = c.quadrature_weights_nodes(b);
		
//≤‚ ‘ x+y;

	double sum = 0.0;
	std::cout << "gauss_nodes" << std::endl;
	for (std::size_t i = 0; i < 4; ++i) {
		auto gp = d[i].first;
		sum += (1.0) * d[i].second;
		std::cout << gp[0] << "," << gp[1] << std::endl;
	}
	auto jacobi_det = [&b]() {
		auto p = b.coordinates;
		return (p[1][0] - p[0][0]) * (p[2][1] - p[0][1]) - (p[2][0] - p[0][0]) * (p[1][1] - p[0][1]);
	};
	std::cout << jacobi_det()*sum << std::endl;

	std::cout << "hello world!" << std::endl;
}