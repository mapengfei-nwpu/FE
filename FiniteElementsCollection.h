#pragma once
#include <array>
#include <iostream>
#include <vector>
#include <memory>
#include <boost/multi_array.hpp>

#include "FiniteElement.h"
#include "MeshElementsCollection.h"

class FiniteElementsCollection :public GeneralElementsCollection<2, 3> {
public:
	/// construct from mesh and order.
	FiniteElementsCollection(MeshElementsCollection mec, std::size_t order = 1) {
		if (order != 1) std::cout << "Finite element order must be one." << std::endl;
		_T = mec._T;
		_P = mec._P;
	}

private:
};

