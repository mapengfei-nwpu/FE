#include <iostream>
#include <stdio.h>
#include "MeshElementsCollection.h"
#include "GaussQuadrature.h"
#include "BasisFunction.h"
int main() {
	/********** test 1 ************
	MeshElementsCollection a;
	a.print();
	*******************************/
	///* test 2
	MeshElementsCollection a;
	auto b = a.get_mesh_element(4);
	a.print();
	b.print();
	GaussQuadrature c;
	auto d = c.quadrature_weights_nodes(b);

	//quadrature test

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

	// basis function test
	BasisFunction f(b);
	auto basis_table = f.basis;
	std::cout << f.get_jacobian_determian() << std::endl;
	std::cout << (f.*basis_table[2]["dy"])(b.coordinates[0][0], b.coordinates[0][1]) << std::endl;
	std::cout << "hello world!" << std::endl;
}