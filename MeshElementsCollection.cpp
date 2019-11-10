#include "MeshElementsCollection.h"

MeshElementsCollection::MeshElementsCollection() {

	const std::size_t nx = 64;
	const std::size_t ny = 64;

	const double a = 0.0;
	const double b = 1.0;
	const double c = 0.0;
	const double d = 1.0;

	std::array<double, 2> node;
	std::size_t vertex = 0;
	for (std::size_t iy = 0; iy <= ny; iy++)
	{
		node[1] = c + ((static_cast<double>(iy)) * (d - c) / static_cast<double>(ny));
		for (std::size_t ix = 0; ix <= nx; ix++)
		{
			node[0] = a + ((static_cast<double>(ix)) * (b - a) / static_cast<double>(nx));
			add_node(node);
		}
	}
	std::vector<std::array<std::size_t, 3>> elements(2);
	for (std::size_t iy = 0; iy < ny; iy++)
	{
		for (std::size_t ix = 0; ix < nx; ix++)
		{
			const std::size_t v0 = iy * (nx + 1) + ix;
			const std::size_t v1 = v0 + 1;
			const std::size_t v2 = v0 + (nx + 1);
			const std::size_t v3 = v1 + (nx + 1);

			elements[0][0] = v0; elements[0][1] = v1; elements[0][2] = v2;
			elements[1][0] = v1; elements[1][1] = v2; elements[1][2] = v3;
			add_element(elements[0]);
			add_element(elements[1]);
		}
	}
}

MeshElement MeshElementsCollection::get_mesh_element(std::size_t i) {
	std::vector<std::array<double,2>> nodes;
	auto element = get_element(i);
	for (std::size_t j = 0; j < 3; ++j) {
		nodes.push_back(get_node(element[j]));
	}
	MeshElement mesh_element(nodes);
	return mesh_element;
}