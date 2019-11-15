#pragma once
#include <dolfin.h>

dolfin::Mesh mesh;

class RegularMesh {
public:
	RegularMesh(std::shared_ptr<const dolfin::Mesh> mesh): 
		_mesh(mesh), _mpi_comm(mesh->mpi_comm()) {
		;
	}
	double x0, x1, y0, y1;
	size_t nx, ny;
	void seperate_mode(std::vector<dolfin::Point> points, std::vector<size_t> dims) {
		
		nx = dims[0];
		ny = dims[0];

		x0 = points[0].x();
		x1 = points[1].x();
		y0 = points[0].y();
		y1 = points[1].y();
	}

	size_t hash(dolfin::Point point) {

		double x = point.x();
		double y = point.y();
		double dx = (x1 - x0) / static_cast<double>(nx);
		double dy = (y1 - y0) / static_cast<double>(ny);

		return static_cast<size_t>((x - x0) / dx) +
			static_cast<size_t>((y - y0) / dy) * nx;
	}

	void index_mesh() {
		// The local map local to global
		std::vector<std::array<size_t, 3>> local_map;
		for (dolfin::CellIterator e(*_mesh); !e.end(); ++e) {
			std::array<size_t, 3> single_map;
			auto center = e->midpoint();
			single_map[0] = _mpi_comm.rank();
			single_map[1] = hash(center);
			single_map[2] = e->global_index();
			local_map.push_back(single_map);
		}
		auto num_cell_global = _mesh->num_entities_global(2);
		global_map.resize(num_cell_global);
		std::vector<std::vector<std::array<size_t,3>>> mpi_collect;
		///dolfin::MPI::all_to_all()
		for (auto iter = mpi_collect.cbegin(); iter != mpi_collect.cend(); iter++) {
			for (auto jter = iter->begin(); jter != iter->cend(); jter++) {
				auto single_map = *jter;
				global_map[single_map[2]][0] = single_map[0];
				global_map[single_map[2]][1] = single_map[1];
			}
		}
	}

	// The map of global index to hash index for cells.
	std::vector<std::array<size_t, 2>> global_map;

private:
	// MPI communicator
	dolfin::MPI::Comm _mpi_comm;

	// The mesh
	std::shared_ptr<const dolfin::Mesh> _mesh;

};