#pragma once
#include <functional>
#include <map>
#include "MeshElement.h"
class BasisFunction {
public:
	std::map<std::size_t, std::function<double(double, double)>> binops;
	BasisFunction(const MeshElement& e) {
		auto c = e.coordinates;
		x1 = c[0][0]; y1 = c[0][1];
		x2 = c[1][0]; y2 = c[1][1];
		x3 = c[2][0]; y3 = c[2][1];
		_j = (x2 - x1) * (y3 - y1) - (x3 - x1) * (y2 - y1);
		binops.insert({ 1, basis_1 });

	}


	double get_jacobian_determian() {
		return _j;
	}

	double basis_1(double x, double y) {
		return -x_reference(x, y) - y_reference(x, y) + 1;
	}
	double basis_1_derivative_x(double x, double y) {
		return -(y3 - y1) / _j - (y1 - y2) / _j;
	}
	double basis_1_derivative_y(double x, double y) {
		return -(x1 - x3) / _j - (x2 - x1) / _j;
	}

	double basis_2(double x, double y) {
		return x_reference(x, y);
	}
	double basis_2_derivative_x(double x, double y) {
		return (y3 - y1) / _j;
	}
	double basis_2_derivative_y(double x, double y) {
		return (x1 - x3) / _j;
	}

	double basis_3(double x, double y) {
		auto z = y_reference(x, y);
		return z;
	}
	double basis_3_derivative_x(double x, double y) {
		return (y1 - y2) / _j;
	}
	double basis_3_derivative_y(double x, double y) {
		return (x2 - x1) / _j;
	}


private:

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