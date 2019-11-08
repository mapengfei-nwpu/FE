#include <vector>
#include <iostream>
#include "MeshElement.h"

/// Compute quadrature rule for triangle.
///
/// *Arguments*
///     coordinates (std::vector<Point>)
///         Vertex coordinates for the triangle
///
/// *Returns*
///     std::pair<std::vector<double>, std::vector<double>>
///         A flattened array of quadrature points and a
///         corresponding array of quadrature weights.
//std::pair<std::vector<double>, std::vector<double>> quadrature_weights_nodes() {

//}

class GaussQuadrature {
public:
	///
	GaussQuadrature(std::size_t order = 4, std::size_t geo_dim = 2)
	{
		if (order != 4 || geo_dim != 2) {
			std::cout << "现在只有四阶精度的三角形积分" << std::endl;
		}
	}

	////
	std::array<std::pair<std::array<double, 2>, double>, 4>
		quadrature_weights_nodes(const MeshElement& me) {

		auto affine_x = [&me](double x, double y) {
			auto p = me.coordinates;
			return (p[1][0] - p[0][0]) * x + (p[2][0] - p[0][0]) * y + x;
		};

		auto affine_y = [&me](double x, double y) {
			auto p = me.coordinates;
			return (p[1][1] - p[0][1]) * x + (p[2][1] - p[0][1]) * y + y;
		};

		std::array<std::pair<std::array<double, 2>, double>, 4> result;
		for (std::size_t i = 0; i < 4; ++i) {
			auto x = affine_x(gauss_nodes[i][0], gauss_nodes[i][1]);
			auto y = affine_y(gauss_nodes[i][0], gauss_nodes[i][1]);
			std::array<double, 2> a = { x, y };
			auto z = std::make_pair(a, gauss_weights[i]);
			result[i] = z;
		}
		return result;
	}

	//
	static const double gauss_nodes[4][2];
	static const double gauss_weights[4];
};

const double GaussQuadrature::gauss_nodes[4][2] =
{ { 1.0 / 3.0,  1.0 / 3.0},
  { 3.0 / 5.0,  1.0 / 5.0},
  { 1.0 / 5.0,  1.0 / 5.0},
  { 1.0 / 5.0,  3.0 / 5.0} };

const double GaussQuadrature::gauss_weights[4] =
{ -27.0 / 48.0, 25.0 / 48.0, 25.0 / 48.0, 25.0 / 48.0 };