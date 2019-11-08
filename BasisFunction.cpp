#include <functional>
#include "Point.h"
#include <boost/multi_array.hpp>

std::function<double(double, double)> phi_x;
auto mod = [](int a, int b) {return a % b;};
void compute_jacobian_determinants() { ; }

auto psi = [](int a) {return a; };
auto psi_derivative_x  = [](int a) {return a; };
auto psi_derivative_y  = [](int a) {return a; };
auto psi_derivative_xy = [](int a) {return 0; };
auto psi_derivative_xy = [](int a) {return 0; };
auto psi_derivative_xx = [](int a) {return 0; };
auto psi_derivative_yy = [](int a) {return 0; };
int main() {
	
	std::vector<std::array<double, 2>> mesh_nodes(100);
	std::vector<std::array<std::size_t, 3>> mesh_elements(100);
	std::vector<std::array<double, 2>> finite_nodes(100);
	std::vector<std::array<std::size_t, 3>> finite_elements(100);
	//finite element is used to index dof. 
	std::vector<std::array<std::size_t, 4>> edges(100);
	auto p1 = mesh_nodes[mesh_elements[0][0]];
	auto p2 = mesh_nodes[mesh_elements[0][1]];
	auto p3 = mesh_nodes[mesh_elements[0][2]];
	boost::multi_array<double, 2> coordinates(3,3);
	for (std::size_t i = 0; i < 3; ++i) {

	}

	std::array<std::array<double, 2>, 3> triangle = { p1, p2, p3 };
	(finite_elements[0][0], finite_elements[0][0]);
	//刚度矩阵的坐标
	//finite_elements[n][i]对应的是第n个单元的第n个基函数。
	
	LocalElement en(&triangle);
	auto local_stiff = en.quadrature();

	const std::size_t dof_dim = 3;
	const std::size_t geo_dim = 2;

	if (coordinates.shape()[0] != dof_dim or coordinates.shape()[1] != geo_dim)
	{
		boost::multi_array<double, 2>::extent_gen extents;
		coordinates.resize(extents[dof_dim][geo_dim]);
	}

}



class LocalElement {
	LocalElement(std::array<std::array<double,2>,3> &e) {
		// calculate the jacobian determiant
		// calculate the gauss quadrature nodes and weights
		// 
	}
public:
	double quadrature() {
		return 1.0;
	}
};





class MeshElementsCollection :public GeneralElementsCollection<3,4> {
public:
	
};
















