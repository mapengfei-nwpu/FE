#pragma once
#include <dolfin.h>

dolfin::Mesh mesh;

class RegularMesh {
public:
	RegularMesh(std::shared_ptr<const dolfin::Mesh> mesh): 
		_mesh(mesh), _mpi_comm(mesh->mpi_comm()) {
		;
	}

	void seperate_mode(std::vector<dolfin::Point> points, std::vector<size_t> dims) {

	}

	size_t hash(dolfin::Point point) {
		return 0;
	}

	void index_mesh() {
		// The local map local to global
		std::vector<std::array<size_t, 3>> local_map;
		for (dolfin::MeshEntityIterator e(*_mesh, 2); !e.end(); ++e) {
			std::array<size_t, 3> single_map;
			auto center = e->midpoint();
			single_map[0] = _mpi_comm.rank();
			single_map[1] = e->index();
			single_map[2] = hash(center);
			local_map.push_back(single_map);
			e->global_index();
		}
		auto num_cell_global = _mesh->num_entities_global(2);
		global_map.resize(num_cell_global);
		///dolfin::MPI::all_to_all()
	}

	std::vector<std::array<size_t, 2>> global_map;
private:
	// MPI communicator
	dolfin::MPI::Comm _mpi_comm;

	// The mesh
	std::shared_ptr<const dolfin::Mesh> _mesh;


};