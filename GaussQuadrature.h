#include <vector>
#include <iostream>
#include "MeshElement.h"
/// Gauss quadrature rule

class GaussQuadrature {
public:
	/// 
	GaussQuadrature(std::size_t order = 4, std::size_t geo_dim = 2)
	{
		if (order != 4 || geo_dim != 2) {
			std::cout << "只能四阶精度的三角形积分" << std::endl;
		}
	}

	/// return an array of local gauss nodes and weights.
	///
	/// 
	/// *Arguments*
	///     element (MeshElement)
	///         The tabulation of which the indices of the element consist
	std::array<std::pair<std::array<double, 2>, double>, 4>
		quadrature_weights_nodes(const MeshElement& me) {

		/// gauss nodes are on the reference triangle
		/// they should be mapped to the local triangle
		auto local_x = [&me](double x, double y) {
			auto p = me.coordinates;
			return (p[1][0] - p[0][0]) * x + (p[2][0] - p[0][0]) * y + p[0][0];
		};

		auto local_y = [&me](double x, double y) {
			auto p = me.coordinates;
			return (p[1][1] - p[0][1]) * x + (p[2][1] - p[0][1]) * y + p[0][1];
		};

		std::array<std::pair<std::array<double, 2>, double>, 4> result;
		for (std::size_t i = 0; i < 4; ++i) {
			auto x = local_x(gauss_nodes[i][0], gauss_nodes[i][1]);
			auto y = local_y(gauss_nodes[i][0], gauss_nodes[i][1]);
			std::array<double, 2> a = { x, y };
			auto z = std::make_pair(a, gauss_weights[i]);
			result[i] = z;
		}
		return result;
	}
/// these should be static variables.
/// but I don't know how to use static variables in class.
	const double gauss_nodes[4][2] =
	{ { 1.0 / 3.0,  1.0 / 3.0},
	  { 3.0 / 5.0,  1.0 / 5.0},
	  { 1.0 / 5.0,  1.0 / 5.0},
	  { 1.0 / 5.0,  3.0 / 5.0} };

	const double gauss_weights[4] =
	{ -27.0 / 96.0, 25.0 / 96.0, 25.0 / 96.0, 25.0 / 96.0 };
};

