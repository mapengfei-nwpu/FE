#include <iostream>
#include <stdio.h>
#include "MeshElementsCollection.h"
#include "FiniteElementsCollection.h"
#include "GaussQuadrature.h"
#include "BasisFunction.h"
#include "Assembler.h"



bool near(double x, double y) {
	if (std::abs(x - y) < 0.000001)
		return true;
	else return false;
}
void dirichlet(MeshElementsCollection& mesh, Eigen::MatrixXd& A) {
	std::size_t I = A.rows();
	for (size_t i = 0; i < I; i++) {
		auto p = mesh.get_node(i);
		if (near(p[0], 0) || near(p[1], 1))
			A(i,i) = 10000000000000000.0;
	}
}


int main() {
	/********** test 1 ************
	MeshElementsCollection a;
	a.print();
	*******************************/
	///* test 2
	/*
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
	*/
	/// test the matrix assembler.
	MeshElementsCollection mesh;
	FiniteElementsCollection functionspace(mesh);
	Assembler a(mesh, functionspace);
	auto A = a.matrixGlobalAssembler();
	std::size_t I = A.rows();
	std::size_t J = A.rows();
	dirichlet(mesh, A);
	for (size_t i = 0; i < I; i++) {
		for (size_t j = 0; j < J; j++)
		{
			//std::cout << A(i, j) << "  ";
		}
		printf("\n");
	}

	//MeshElementsCollection mesh;
	//FiniteElementsCollection functionspace(mesh);
	//Assembler a(mesh, functionspace);
	auto F = a.vectorGlobalAssembler();

	//std::size_t I = F.rows();
	for (size_t i = 0; i < I; i++) {

		std::cout <<F(i)<< "  ";
		printf("\n");
	}
	//MeshElementsCollection mesh;
	//FiniteElementsCollection functionspace(mesh);
	//Assembler a(mesh, functionspace);
	//auto A = a.matrixGlobalAssembler();
	//auto F = a.vectorGlobalAssembler();
	//dirichlet(mesh, F);
	Eigen::VectorXd X = A.lu().solve(F);
	
	//std::size_t I = X.rows();
	//std::size_t J = X.cols();
	for (size_t i = 0; i < I; i++) {
		auto p = mesh.get_node(i); 
		std::cout << p[0] << "        " << p[1] << "        " << X(i) << "     " << std::endl;
	}

	std::cout << F.sum() << std::endl;

}
