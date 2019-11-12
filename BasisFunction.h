#pragma once
#include <vector>
#include <iostream>
#include <Eigen/Dense>
#include "MeshElement.h"
#include "FiniteElement.h"
#define M 6


class BasisFunction {
public:
	BasisFunction(const FiniteElement& e) {
		if (e.coordinates.size() != M) std::cout << "wrong" << std::endl;
		std::array<std::array<double, 2>, M> coord;
		for (size_t i = 0; i < M; i++)
			for (size_t j = 0; j < 2; j++)
			{
				coord[i][j] = e.coordinates[i][j];
			}
		kappa_calculation(coord);
		/// _j = (coord[1][0] - coord[0][0]) * (coord[2][1] - coord[0][1])
		///	- (coord[2][0] - coord[0][0]) * (coord[1][1] - coord[0][1]);
		_j = (coord[2][0] - coord[0][0]) * (coord[4][1] - coord[0][1])
			- (coord[4][0] - coord[0][0]) * (coord[2][1] - coord[0][1]);

	}


	double get_jacobian_determian() {
		return _j;
	}

	double psi(double x, double y, size_t i) {
		if (i < 0 || i >= M) {
			std::cout << "wrong" << std::endl; return 1.0;
		}
		// return kappa[i][0] * x + kappa[i][1] * y + kappa[i][2];
		return kappa[i][0] * x * x + kappa[i][1] * y * y + kappa[i][2] * x * y
			+ kappa[i][3] * x + kappa[i][4] * y + kappa[i][5];
	}
	double psi_x(double x, double y, size_t i) {
		// return kappa[i][0];
		return kappa[i][0] * x * 2.0 + kappa[i][2] * y + kappa[i][3];
	}
	double psi_y(double x, double y, size_t i) {
		return kappa[i][1] * y * 2.0 + kappa[i][2] * x + kappa[i][4];
	}
	void kappa_calculation(std::array<std::array<double, 2>, M> coordinates) {

		/// set matrix entries.
		Eigen::MatrixXd A(M, M);
		for (size_t j = 0; j < M; j++)
		{
			auto c = coordinates[j];
			A(j, 0) = c[0] * c[0];
			A(j, 1) = c[1] * c[1];
			A(j, 2) = c[0] * c[1];
			A(j, 3) = c[0];
			A(j, 4) = c[1];
			A(j, 5) = 1.0;
		}

		/// set vector entries and solve parameters.
		Eigen::VectorXd e(M);
		for (size_t i = 0; i < M; ++i)
		{
			e.setZero();
			e(i) = 1.0;
			Eigen::VectorXd x = A.lu().solve(e);
			/// std::cout << A << std::endl;
			/// std::cout << e << std::endl;
			/// std::cout << x << std::endl;
			/// copy result to parameters matrix.
			for (size_t j = 0; j < M; j++)
			{
				auto temp = x(j);
				kappa[i][j] = x(j);
			}
		}

		/// output parameter matrix.
		/*for (size_t i = 0; i < M; i++)
		{
			for (size_t j = 0; j < M; j++)
			{
				std::cout << kappa[i][j] << "   ";
			}
			std::cout << std::endl;

		}*/

	}

private:
	/// parameters matrix of basis functions.
	double kappa[M][M];
	/// the determinant of jacobian matrix.
	double _j;
};