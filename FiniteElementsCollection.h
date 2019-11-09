#pragma once
#include <array>
#include <vector>
#include <memory>
#include <boost/multi_array.hpp>

#include "FiniteElement.h"
#include "MeshElementsCollection.h"

/// I found that GeneralElementCollection cannot use template.
class FiniteElementsCollection :public GeneralElementsCollection<2, 3> {
public:
	FiniteElementsCollection(MeshElementsCollection mec);
	std::size_t element_dim();
	std::size_t element_number();
	std::size_t dof_number() { return element_dim() * element_number(); };
	FiniteElement get_Finite_element(std::size_t i);

private:
	// The mesh
	MeshElementsCollection *_mesh;

	// Topological dimension
	std::size_t _dim;

	// Local index of entity within topological dimension
	std::size_t _local_index;

};

