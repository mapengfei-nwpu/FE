#include <iostream>
#include <dolfin.h>
using namespace dolfin;

class BoxAdjacents
{
public:
	BoxAdjacents(std::array<dolfin::Point, 2> points, std::vector<size_t> dims, CellType::Type cell_type)
	{

		// toplogy dimesion
		top_dim = dims.size();

		// generate mesh
		if (cell_type == CellType::Type::hexahedron)
		{
			dolfin_assert(top_dims == 3);
			auto _mesh = BoxMesh::create(points, { dims[0], dims[1], dims[2] }, cell_type);
			mesh_ptr = std::make_shared<Mesh>(_mesh);
		}
		else if (cell_type == CellType::Type::quadrilateral)
		{
			dolfin_assert(top_dims == 2);
			auto _mesh = RectangleMesh::create(points, { dims[0], dims[1] }, cell_type);
			mesh_ptr = std::make_shared<Mesh>(_mesh);
		}
		else
		{
			dolfin_error("BoxMesh.h",
				"generate uniform mesh",
				"Wrong cell type '%d'", cell_type);
		}

		// check again.
		if (top_dim != 2 && top_dim != 3)
			dolfin_error("the size of dims must be 2 & 3.", ".", ".");

		nx = dims[0];
		ny = dims[1];
		nz = 0;

		x0 = points[0].x();
		x1 = points[1].x();
		y0 = points[0].y();
		y1 = points[1].y();
		z0 = 0.0;
		z1 = 0.0;

		if (top_dim == 3)
		{
			nz = dims[2];
			z0 = points[0].z();
			z1 = points[1].z();
		}
		index_mesh();
	}
	size_t hash(dolfin::Point point)
	{
		if (top_dim != 2 && top_dim != 3)
			dolfin_error("the size of dims must be 2 and 3.", ".", ".");

		double x = point.x();
		double y = point.y();
		double z = point.z();

		double dx = (x1 - x0) / static_cast<double>(nx);
		double dy = (y1 - y0) / static_cast<double>(ny);
		double dz = 0.0;

		size_t i = static_cast<size_t>((x - x0) / dx);
		size_t j = static_cast<size_t>((y - y0) / dy);
		size_t k = 0;

		if (top_dim == 3)
		{
			dz = (z1 - z0) / static_cast<double>(nz);
			k = static_cast<size_t>((z - z0) / dz);
		}

		return k * nx * ny + j * nx + i;
	}
	bool check()
	{
		for (CellIterator cell(*mesh_ptr); !cell.end(); ++cell)
		{
			auto point = cell->midpoint();
			std::cout << "global index: "
				<< cell->global_index()
				<< " or "
				<< hash(point)
				<< ". local index: "
				<< cell->index()
				<< " or "
				<< map(hash(point))[1]
				<< ". mpi rank: "
				<< map(hash(point))[0]
				<< ". number of adjacents: "
				// << get_adjacents(point).size()
				<< std::endl;
			if (cell->global_index() != hash(cell->midpoint()))
			{
				return false;
			}
		}
		return true;
	}

	// point -> global index of adjacent cells -> local index of adjacent cells
	std::vector<size_t> get_adjacents(Point point)
	{
		/// only finished for 2d occassion.
		std::vector<size_t> adjacents;
		size_t global_index = hash(point);
		std::vector<size_t> a = {
			global_index - 2 * nx - 2, 	global_index - 2 * nx - 1, 	global_index - 2 * nx, 	global_index - 2 * nx + 1, 	global_index - 2 * nx + 2,
			global_index - nx - 2, 		global_index - nx - 1, 		global_index - nx, 		global_index - nx + 1, 		global_index - nx + 2,
			global_index - 2, 			global_index - 1, 			global_index, 			global_index + 1, 			global_index + 2,
			global_index + nx - 2, 		global_index + nx - 1, 		global_index + nx, 		global_index + nx + 1, 		global_index + nx + 2,
			global_index + 2 * nx - 2, 	global_index + 2 * nx - 1, 	global_index + 2 * nx, 	global_index + 2 * nx + 1, 	global_index + 2 * nx + 2 };
		for (size_t i = 0; i < a.size(); i++)
		{
			if (a[i] >= 0 && a[i] < global_map.size())
			{

				if (map(a[i])[0] == dolfin::MPI::rank(mesh_ptr->mpi_comm()))
				{
					adjacents.push_back(map(a[i])[1]);
				}
			}
		}
		return adjacents;
	}

	// global index and cell center
	std::array<size_t, 2> map(size_t i)
	{
		return global_map[i];
	}
	std::shared_ptr<Mesh> mesh()
	{
		return mesh_ptr;
	}

	std::vector<double> side_length()
	{
		std::vector<double> side_lengths;
		side_lengths.push_back((x1 - x0) / nx);
		side_lengths.push_back((y1 - y0) / ny);
		return side_lengths;
	}

private:
	void index_mesh()
	{
		// The local map local to global
		std::vector<size_t> local_map;
		for (dolfin::CellIterator e(*mesh_ptr); !e.end(); ++e)
		{
			auto center = e->midpoint();
			local_map.push_back(e->global_index());
			local_map.push_back(dolfin::MPI::rank(mesh_ptr->mpi_comm()));
			local_map.push_back(e->index());
		}
		// send local map to every peocess.
		std::vector<std::vector<size_t>> mpi_collect(dolfin::MPI::size(mesh_ptr->mpi_comm()));
		dolfin::MPI::all_gather(mesh_ptr->mpi_comm(), local_map, mpi_collect);
		// alloc memory for global map.
		auto num_cell_global = mesh_ptr->num_entities_global(2);
		global_map.resize(num_cell_global);
		for (auto iter = mpi_collect.cbegin(); iter != mpi_collect.cend(); iter++)
		{
			for (auto jter = iter->begin(); jter != iter->cend();)
			{
				size_t cell_index = *jter;
				jter++;
				global_map[cell_index][0] = *jter; /// mpi_rank;size_t
				jter++;
				global_map[cell_index][1] = *jter; /// cell hash;
				jter++;
			}
		}
	}

	double x0, x1, y0, y1, z0, z1;
	size_t nx, ny, nz;
	size_t top_dim;
	// The map of global index to hash index for cells.
	std::vector<std::array<size_t, 2>> global_map;
	std::shared_ptr<Mesh> mesh_ptr;
};
