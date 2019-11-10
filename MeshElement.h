#pragma once
#include <vector>
#include <array>
#include <boost/multi_array.hpp>
#include <memory>

class MeshElement {
public:
	/// construct MeshElement with coordinates.
	MeshElement(std::vector<std::array<double, 2>>);
	
	/// print coordinates.
	void print();

	/// coordinate of three nodes
	boost::multi_array<double, 2> coordinates;

};