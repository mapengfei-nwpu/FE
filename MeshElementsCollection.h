#pragma once

#include "GeneralElementsCollection.h"
#include "MeshElement.h"

class MeshElementsCollection :public GeneralElementsCollection<2, 3> {
public:
	MeshElementsCollection();
	MeshElement get_mesh_element(std::size_t i);
};
