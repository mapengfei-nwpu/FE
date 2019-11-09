#pragma once
#include <vector>
#include <array>
#include <boost/multi_array.hpp>
#include "MeshElementsCollection.h"

class MeshElement {
public:
	/// 
	MeshElement(std::vector<std::array<double, 2>>);
	
	/// print coordinates.
	void print();

	/// coordinate of three nodes
	boost::multi_array<double, 2> coordinates;

private:
	std::shared_ptr<MeshElementsCollection> _m;

};