#pragma once
#include <vector>
#include <array>
#include <boost/multi_array.hpp>
#include <memory>
#include <iostream>

class FiniteElement {
public:
	/// construct MeshElement with coordinates.
	FiniteElement(const std::vector<std::array<double, 2>> & nodes) {
		//1. the size of nodes must be three.
		const std::size_t dof_dim = 6;
		const std::size_t geo_dim = 2;
		if (coordinates.shape()[0] != dof_dim or coordinates.shape()[1] != geo_dim)
		{
			boost::multi_array<double, 2>::extent_gen extents;
			coordinates.resize(extents[dof_dim][geo_dim]);
		}
		for (std::size_t i = 0; i < 6; i++) {
			for (std::size_t j = 0; j < 2; j++) {
				coordinates[i][j] = nodes[i][j];
			}
		}
	}

	/// print coordinates.
	void print() const {
		std::cout << "Points:" << std::endl;
		for (std::size_t i = 0; i < 6; ++i) {
			std::cout << coordinates[i][0] << "," << coordinates[i][1] << std::endl;
		}
	}

	/// coordinate of three nodes
	boost::multi_array<double, 2> coordinates;

};