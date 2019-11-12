#pragma once
#include <functional>
#include <Eigen/Dense>
#include <map>
#include <string>
#include <vector>
#include "MeshElement.h"
#define M 6
class BasisFunction;
typedef double(BasisFunction::* basis_type)(double, double);


class BasisFunction {
public:
	std::vector<std::map<std::string, basis_type>> basis;
	BasisFunction(const MeshElement& e) {
		auto c = e.coordinates;
		x1 = c[0][0]; y1 = c[0][1];
		x2 = c[1][0]; y2 = c[1][1];
		x3 = c[2][0]; y3 = c[2][1];
		_j = (x2 - x1) * (y3 - y1) - (x3 - x1) * (y2 - y1);

		double A[M][M];
		double e[M];
		/// copy to kappa. 
	}


	double get_jacobian_determian() {
		return _j;
	}

	double psi(double x, double y, size_t i) {
		if (i < 0 || i >= M) {
			std::cout << "wrong" << std::endl; return 1.0;
		}
		return kappa[i][0] * x * x + kappa[i][1] * y * y + kappa[i][2] * x * y
			+ kappa[i][3] * x + kappa[i][4] * y + kappa[i][5];
	}
	double psi_x(double x, double y, size_t i) {
		return kappa[i][0] * x * 2.0 + kappa[i][2] * y + kappa[i][3];
	}
	double psi_y(double x, double y, size_t i) {
		return kappa[i][1] * y * 2.0 + kappa[i][2] * x + kappa[i][4];
	}
	void kappa_calculation(std::array<std::array<double,2>,3> coordinates,size_t i) {
		coordinates[0][0];
		Eigen::MatrixXd A(3, 3);
		Eigen::VectorXd e = Eigen::VectorXd::Zero(3);
		e(i) = 1.0;
		for (size_t j = 0; j < 3; ++j)
		{
			auto c = coordinates[j];
			A(j, 0) = c[0] * c[0];
			A(j, 1) = c[1] * c[1];
			A(j, 1) = c[1] * c[1];
		}
	}

private:


	double kappa[M][M];
	// the determinant of jacobian matrix.
	double _j;
	// coordinates
	double x1, x2, x3, y1, y2, y3;

	// reference coordinates.
	double x_reference(double x, double y) {
		return ((y3 - y1) * (x - x1) - (x3 - x1) * (y - y1)) / _j;
	};
	double y_reference(double x, double y) {
		return (-(y2 - y1) * (x - x1) + (x2 - x1) * (y - y1)) / _j;
	};
};