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


	std::size_t element_dim() { return 3; }
	
	/// the tabulation of local dof and global
	/// dof for element i.  
	std::array<std::size_t, 3> dofmap(std::size_t i) {
		return _T[i];
	}

private:
	// The mesh
	// MeshElementsCollection *_mesh;

	// Topological dimension
	// std::size_t _dim;

	// Local index of entity within topological dimension
	// std::size_t _local_index;

};

