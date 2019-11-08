#include "MeshElement.h"
#include <iostream>
MeshElement::MeshElement(std::vector<std::array<double, 2>> nodes) {
	//1. the size of nodes must be three.
	const std::size_t dof_dim = 3;
	const std::size_t geo_dim = 2;
	if (coordinates.shape()[0] != dof_dim or coordinates.shape()[1] != geo_dim)
	{
		boost::multi_array<double, 2>::extent_gen extents;
		coordinates.resize(extents[dof_dim][geo_dim]);
	}
	for (std::size_t i = 0; i < 3; i++) {
		for (std::size_t j = 0; j < 2; j++) {
			coordinates[i][j] = nodes[i][j]; 
		}
	}
}
void MeshElement::print() {
	std::cout << "Points:" << std::endl;
	for (std::size_t i = 0; i < 3; ++i) {
		std::cout << coordinates[i][0] << "," << coordinates[i][1] << std::endl;
	}
}