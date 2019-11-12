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
void dirichlet(FiniteElementsCollection& functionspace, Eigen::SparseMatrix<double>& A) {
	std::size_t I = A.rows();
	for (size_t i = 0; i < I; i++) {
		auto p = functionspace.get_node(i);
		if (near(p[0], 0) || near(p[1], 1))
			A.coeffRef(i, i) = 1000000000;
	}
}


int main() {
	/// assemble the matrix.
	MeshElementsCollection mesh;
	FiniteElementsCollection functionspace(mesh);
	auto temp1 = functionspace.get_finite_element(0);
	temp1.print();
	
	Assembler a(mesh, functionspace);
	auto A = a.matrixGlobalAssembler();
	/// std::cout << A << std::endl;
	std::size_t I = A.rows();
	std::size_t J = A.rows();
	dirichlet(functionspace, A);
	for (size_t i = 0; i < I; i++) {
		for (size_t j = 0; j < J; j++)
		{
			//std::cout << A.coeffRef(i, j) << "  ";
		}
		/// printf("\n");
	}
	/// assemble the right side vector.
	auto F = a.vectorGlobalAssembler();
	auto temp = 0.0;
	for (size_t i = 0; i < I; i++) {
		temp += F(i);
	}
	std::cout << temp << std::endl;

	/// solve the linear system.
	Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
	solver.compute(A);
	auto X = solver.solve(F);
	for (size_t i = 0; i < I; i++) {
		auto p = functionspace.get_node(i); 
		std::cout << p[0] << "        " << p[1] << "        " << X(i) << "     " << std::endl;
	}
}
