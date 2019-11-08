#pragma once
#include <vector>
#include <array>
#include <boost/multi_array.hpp>
class MeshElement {
public:
	/// three points and every point point contains two coordinates.
	MeshElement(std::vector<std::array<double, 2>>);
	void print();
	boost::multi_array<double, 2> coordinates;
};