#include <vector>
#include <stdio.h>

#include <dolfin.h>

#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h> // which can return an object of class vector.
#include <pybind11/stl_bind.h>
PYBIND11_MAKE_OPAQUE(std::vector<double>);
namespace py = pybind11;

class RegularMesh {
public:
	RegularMesh(std::shared_ptr<const dolfin::Mesh> mesh) :
		_mesh(mesh)/*, _mpi_comm(mesh->mpi_comm())*/ {
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
		std::vector<size_t> local_map;
		for (dolfin::CellIterator e(*_mesh); !e.end(); ++e) {
			auto center = e->midpoint();
			local_map.push_back(e->global_index());
			local_map.push_back(dolfin::MPI::rank(_mesh->mpi_comm()));
			local_map.push_back(hash(center));
		}
		// send local map to every peocess.
		std::vector<std::vector<size_t>> mpi_collect(dolfin::MPI::size(_mesh->mpi_comm()));
		dolfin::MPI::all_gather(_mesh->mpi_comm(), local_map, mpi_collect);
		// alloc memory for global map.
		auto num_cell_global = _mesh->num_entities_global(2);
		global_map.resize(num_cell_global);
		for (auto iter = mpi_collect.cbegin(); iter != mpi_collect.cend(); iter++) {
			for (auto jter = iter->begin(); jter != iter->cend();) {
				size_t cell_index = *jter;
				jter++;
				global_map[cell_index][0] = *jter; /// mpi_rank;size_t
				jter++;
				global_map[cell_index][1] = *jter; /// cell hash;
				jter++;
			}
		}
	}
	std::array<size_t, 2> map(size_t i) {
		return global_map[i];
	}
	// The map of global index to hash index for cells.
	std::vector<std::array<size_t, 2>> global_map;

private:
	// MPI communicator
	// dolfin::MPI::Comm _mpi_comm;

	// The mesh
	std::shared_ptr<const dolfin::Mesh> _mesh;

};


class DeltaInterplation {
	/// collect all body vectors on every processes and 
	/// assign them to a function.
	void fun0(dolfin::Function& u, std::vector<double>& body) {
		/// TODO : all reduce body.
		/// dolfin::MPI::all_reduce(body)
		auto size = u.vector()->local_size();
		auto start = u.vector()->local_range().first;
		std::vector<double> values(size);
		for (size_t i = 0; i < size; ++i)
			values[i] = body[i + start];
		u.vector()->set_local(values);
	}
	/// distribute v to vector body on every process.
	void fun1(dolfin::Function& v, std::vector<double>& body.boost::multi_array<double, 2> & body_coordinates) {

		/// there is no need to iterate all mesh.
		auto mesh = v.function_space()->mesh();
		auto dofmap = v.function_space()->dofmap();

		/// TODO : return the reference of the global_map.
		/// std::vector<std::array<size_t,2>> get_global_map(); 
		for (size_t i = 0; i < body.size(); i++) {
			std::array<double, 2> body_coordinate;
			body_coordinate[0] = body_coordinates[i][0];
			body_coordinate[1] = body_coordinates[i][1];
			/// TODO : return the adjacents for given index in mesh.
			/// std::vector<size_t> get_adjacent(dolfin::Mesh& mesh, std::array<double,2> body_coordinate);
			std::vector<size_t> adjacents;
			for (size_t j = 0; j < adjacents.size(); j++) {
				/// TODO : Cell constructor take local index to initial.
				///        so it is waiting to be corrected.
				Cell cell(*mesh, adjacents[j]);
				boost::multi_array<double, 2> coordinates;
				std::vector<double> coordinate_dofs;
				auto _element = v.function_space()->element();
				_element->tabulate_dof_coordinates(coordinates, coordinate_dofs, *cell);
				/// local index of the cell is needed rather than global index.
				auto cell_dofmap = dofmap->cell_dofs(cell.index());
				for (size_t k = 0; k < 3; k++) {
					/// TODO : multiply with delta function.
					/// double delta(double x, double y);
					double delta = 0.001;
					/// four components are needed here: body_coordinate, coordinates, coordinate_dofs, cel_dofmap
					body[cell_dofmap[i]] += delta(body_coordinate, coordinates[k]) * coordinate_dofs[k];
				}
			}
		}
	}
}

PYBIND11_MODULE(SIGNATURE, m)
{
	py::class_<RegularMesh>(m, "RegularMesh")
		.def(py::init<std::shared_ptr<const dolfin::Mesh>>())
		.def("index_mesh", &RegularMesh::index_mesh)
		.def("map", &RegularMesh::map)
		.def("hash", &RegularMesh::hash)
		.def("seperate_mode", &RegularMesh::seperate_mode);
}
