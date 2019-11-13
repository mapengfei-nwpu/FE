#pragma once
#include<vector>
#include "FiniteElementsCollection.h"
#include <Eigen/Sparse>
#include <Eigen/Dense>

typedef std::vector<std::array<size_t, 4>> Boundary;
class BoundaryCondition {
public:
	BoundaryCondition(FiniteElementsCollection& functionspace, Boundary boundary)
		:functionspace(functionspace),
		boundary(boundary) {
		;
	}
	void dirichlet(Eigen::SparseMatrix<double>& A, Eigen::VectorXd& b,/* a callable object */) {
		auto elements = functionspace._T;
		for (size_t i = 0; i < boundary.size(); i++)
		{
			if (boundary[i][0] == 0) /// if it is dirichlet boundary condition
			{
				auto elem = elements[boundary[i][1]];
				for (size_t j = 0; j < elem.size(); j++)
				{
					A.coeffRef(boundary[i][2], elem[j]) = 0.0;
				}
				A.coeffRef(boundary[i][2], boundary[i][2]) = 1.0;
				b.coeffRef(boundary[i][2]) = 0.0;///zero boundary condition.
			}
		}
	}
private:

	FiniteElementsCollection functionspace;
	Boundary boundary;
};
