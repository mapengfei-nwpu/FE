#pragma once
#include <functional>
#include <map>
#include <string>
#include <vector>
#include "MeshElement.h"

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
		/// this part is ugly.
		/// "bind" can be used to replace these lines.
		std::map<std::string, basis_type> basis_map; 
		basis_map.erase(basis_map.begin(), basis_map.end());
		basis_map["original"] = &BasisFunction::basis_1;
		basis_map["dx"]		  = &BasisFunction::basis_1_dx;
		basis_map["dy"]       = &BasisFunction::basis_1_dy;
		basis.push_back(basis_map);
		basis_map.erase(basis_map.begin(), basis_map.end());
		basis_map["original"] = &BasisFunction::basis_2;
		basis_map["dx"]       = &BasisFunction::basis_2_dx;
		basis_map["dy"]       = &BasisFunction::basis_2_dy;
		basis.push_back(basis_map);
		basis_map.erase(basis_map.begin(), basis_map.end());
		basis_map["original"] = &BasisFunction::basis_3;
		basis_map["dx"]       = &BasisFunction::basis_3_dx;
		basis_map["dy"]       = &BasisFunction::basis_3_dy;
		basis.push_back(basis_map);
	}


	double get_jacobian_determian() {
		return _j;
	}

	double basis_1(double x, double y) {
		return -x_reference(x, y) - y_reference(x, y) + 1;
	}
	double basis_1_dx(double x, double y) {
		return -(y3 - y1) / _j - (y1 - y2) / _j;
	}
	double basis_1_dy(double x, double y) {
		return -(x1 - x3) / _j - (x2 - x1) / _j;
	}

	double basis_2(double x, double y) {
		return x_reference(x, y);
	}
	double basis_2_dx(double x, double y) {
		return (y3 - y1) / _j;
	}
	double basis_2_dy(double x, double y) {
		return (x1 - x3) / _j;
	}

	double basis_3(double x, double y) {
		auto z = y_reference(x, y);
		return z;
	}
	double basis_3_dx(double x, double y) {
		return (y1 - y2) / _j;
	}
	double basis_3_dy(double x, double y) {
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