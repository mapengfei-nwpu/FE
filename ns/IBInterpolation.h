#include <dolfin.h>
#include <numeric>
#include "BoxAdjacents.h"
using namespace dolfin;
std::vector<std::array<double, 2>> get_global_dof_coordinates(const Function &f)
{
	/// some shorcut
	auto mesh = f.function_space()->mesh();
	auto mpi_comm = mesh->mpi_comm();
	auto mpi_size = dolfin::MPI::size(mpi_comm);

	/// get local coordinate_dofs
	auto local_dof_coordinates = f.function_space()->tabulate_dof_coordinates();

	/// collect local coordinate_dofs on every process
	std::vector<std::vector<double>> mpi_dof_coordinates(mpi_size);
	dolfin::MPI::all_gather(mpi_comm, local_dof_coordinates, mpi_dof_coordinates);
	std::vector<double> dof_coordinates_long;

	/// unwrap mpi_dof_coordinates.
	for (size_t i = 0; i < mpi_dof_coordinates.size(); i++)
	{
		for (size_t j = 0; j < mpi_dof_coordinates[i].size(); j++)
		{
			dof_coordinates_long.push_back(mpi_dof_coordinates[i][j]);
		}
	}
	std::vector<std::array<double, 2>> dof_coordinates(dof_coordinates_long.size() / 4);
	for (size_t i = 0; i < dof_coordinates_long.size(); i += 4)
	{
		dof_coordinates[i / 4][0] = dof_coordinates_long[i];
		dof_coordinates[i / 4][1] = dof_coordinates_long[i + 1];
	}

	/// whatch the type of this function return.
	/// it could be changed to "std::vector<std::array<double,3>>" if necessary.
	return dof_coordinates;
}

std::vector<double> get_global_dofs(const Function &f)
{
	/// some shorcut
	auto mesh = f.function_space()->mesh();
	auto mpi_comm = mesh->mpi_comm();
	auto mpi_size = dolfin::MPI::size(mpi_comm);

	/// get local values.
	std::vector<double> local_values;
	f.vector()->get_local(local_values);

	/// collect local values on every process.
	std::vector<std::vector<double>> mpi_values(dolfin::MPI::size(mpi_comm));
	dolfin::MPI::all_gather(mpi_comm, local_values, mpi_values);
	std::vector<double> values;

	/// unwrap mpi_values.
	for (size_t i = 0; i < mpi_values.size(); i++)
	{
		for (size_t j = 0; j < mpi_values[i].size(); j++)
		{
			values.push_back(mpi_values[i][j]);
		}
	}
	return values;
}

class DeltaInterplation
{
public:
	/// information about mesh structure.
	BoxAdjacents &um;
	std::vector<double> side_lengths;

	/// construct function.
	DeltaInterplation(BoxAdjacents &uniform_mesh) : um(uniform_mesh)
	{
		side_lengths = um.side_length();
	}

	/// Assign the solid displacement with the velocity of fluid.
	void fluid_to_solid(Function &fluid, Function &solid)
	{
		/// calculate global dof coordinates and dofs.
		auto dof_coordinates = get_global_dof_coordinates(solid);
		std::vector<double> values(dof_coordinates.size() * 2);
		fluid_to_solid_raw(fluid, values, dof_coordinates);

		/// collect all body vectors on every processes and
		/// assign them to a function.
		auto size = solid.vector()->local_size();
		auto offset = solid.vector()->local_range().first;
		std::vector<double> local_values(size);
		for (size_t i = 0; i < size; ++i)
			local_values[i] = values[i + offset];
		solid.vector()->set_local(local_values);

		/// Finalize assembly of tensor.
		/// this step is quite important.
		solid.vector()->apply("insert");
	}

	void fluid_to_solid_raw(Function &v, std::vector<double> &body,
							std::vector<std::array<double, 2>> &body_coordinates)
	{
		/// the meshes of v and um should be the same.
		/// TODO : compare two meshes
		dolfin_assert(v.value_size() == body.size() / body_coordinates.size());

		/// smart shortcut
		auto rank = dolfin::MPI::rank(v.function_space()->mesh()->mpi_comm());
		auto mesh = v.function_space()->mesh();		// pointer to a mesh
		auto dofmap = v.function_space()->dofmap(); // pointer to a dofmap

		std::cout << "value size: " << v.value_size() << std::endl;
		std::cout << "value size: " << body.size() << std::endl;

		/// iterate every body coordinate.
		for (size_t i = 0; i < body.size() / v.value_size(); i++)
		{
			/// initialize body to zero.
			for (size_t l = 0; l < v.value_size(); l++)
			{
				body[i * v.value_size() + l] = 0.0;
			}

			Point body_coordinate(body_coordinates[i][0], body_coordinates[i][1]);
			auto adjacents = um.get_adjacents(body_coordinate);
			for (size_t j = 0; j < adjacents.size(); j++)
			{

				/// Cell constructor take local index to initial
				Cell cell(*mesh, adjacents[j]);
				boost::multi_array<double, 2> coordinates; /// coordinates of the cell
				std::vector<double> coordinate_dofs;	   /// coordinate_dofs of the cell
				cell.get_coordinate_dofs(coordinate_dofs);
				auto _element = v.function_space()->element(); /// element of the function space
				/// get coordinates and coordinate_dofs
				_element->tabulate_dof_coordinates(coordinates, coordinate_dofs, cell);
				/// local index of cell is needed rather than global index.
				auto cell_dofmap = dofmap->cell_dofs(cell.index());
				for (size_t k = 0; k < cell_dofmap.size() / v.value_size(); k++)
				{
					Point on_cell(coordinates[k][0], coordinates[k][1]);
					double d = delta(body_coordinate, on_cell);
					for (size_t l = 0; l < v.value_size(); l++)
					{
						body[i * v.value_size() + l] += d * (*(v.vector()))[cell_dofmap[k] + l] * side_lengths[0] * side_lengths[1];
					}

					///////////////////////////// WATCH OUT! /////////////////////////////////
					///                                                                    ///
					///  vector in a Function is initialized with a dofmap. this dofmap    ///
					///  contains many informations especially the layout of the vector.   ///
					///  Besides, dofmap tells which is the ghost entries of a element.    ///
					///  see Function::init_vector()                                       ///
					///  and https://fenicsproject.discourse.group/t/it-seems-that         ///
					///  -local-size-didnt-return-real-local-size-of-a-vector/1929         ///
					///                                                                    ///
					//////////////////////////////////////////////////////////////////////////

				} // end loop inside the cell
			}
			// end adjacents loop
		}
		////////////////////////// gather body on one processor //////////////////////////
		//////////////////  TODO : this part can use MPI_reduce directly  ////////////////
		std::vector<std::vector<double>> mpi_collect(dolfin::MPI::size(mesh->mpi_comm()));
		dolfin::MPI::all_gather(mesh->mpi_comm(), body, mpi_collect);
		for (size_t i = 0; i < body.size(); i++)
		{
			body[i] = 0;
			for (size_t j = 0; j < mpi_collect.size(); j++)
			{
				body[i] += mpi_collect[j][i];
			}
		}
		//////////////////////////////////////////////////////////////////////////////////
	}

	/// Assign the solid displacement with the velocity of fluid.
	void solid_to_fluid(Function &fluid, Function &solid)
	{
		/// calculate global dof coordinates and dofs of solid.
		auto solid_dof_coordinates = get_global_dof_coordinates(solid);
		auto solid_values = get_global_dofs(solid);

		/// interpolates solid values into fluid mesh.
		/// the returned fluid_values is the dofs of fluid function.
		auto fluid_values = solid_to_fluid_raw(fluid, solid_values, solid_dof_coordinates);

		/// assemble fluid_values into a function.
		auto local_size = fluid.vector()->local_size();
		auto offset = fluid.vector()->local_range().first;
		std::vector<double> local_values(local_size);
		for (size_t i = 0; i < local_size; ++i)
		{
			local_values[i] = fluid_values[i + offset];
		}
		fluid.vector()->set_local(local_values);
		fluid.vector()->apply("insert");
	}

	///
	std::vector<double> solid_to_fluid_raw(Function &fluid, std::vector<double> &solid_values,
										   std::vector<std::array<double, 2>> &solid_coordinates)
	{
		/// smart shortcut
		std::cout << "interpolate solid to fluid now." << std::endl;
		auto rank = dolfin::MPI::rank(fluid.function_space()->mesh()->mpi_comm());
		auto mesh = fluid.function_space()->mesh();		// pointer to a mesh
		auto dofmap = fluid.function_space()->dofmap(); // pointer to a dofmap

		/// Get local to global dofmap
		std::vector<size_t> local_to_global;
		dofmap->tabulate_local_to_global_dofs(local_to_global);

		auto offset = fluid.vector()->local_range().first;

		/// get the element of function space
		auto element = fluid.function_space()->element();
		auto value_size = fluid.value_size();
		auto global_fluid_size = fluid.function_space()->dim();

		/// initial local fluid values.
		std::vector<double> local_fluid_values(global_fluid_size);

		/// iterate every body coordinate.
		for (size_t i = 0; i < solid_values.size() / value_size; i++)
		{

			/// get indices of adjacent cells on fluid mesh.
			Point solid_point(solid_coordinates[i][0], solid_coordinates[i][1]);
			auto adjacents = um.get_adjacents(solid_point);

			/// iterate adjacent cells.
			for (size_t j = 0; j < adjacents.size(); j++)
			{
				/// 1. get dofmap of the cell.
				/// 2. get the coordinates of the cell.

				/// step 1
				/// get the cell
				Cell cell(*mesh, adjacents[j]);
				/// get the coordinate_dofs of the cell
				std::vector<double> coordinate_dofs;
				cell.get_coordinate_dofs(coordinate_dofs);
				/// get the coordinates of the cell
				boost::multi_array<double, 2> coordinates;
				element->tabulate_dof_coordinates(coordinates, coordinate_dofs, cell);

				/// step 2
				/// local index of cell is needed rather than global index.
				auto cell_dofmap = dofmap->cell_dofs(cell.index());

				/// iterate coordinates on cell.
				for (size_t k = 0; k < cell_dofmap.size() / value_size; k++)
				{
					Point cell_point(coordinates[k][0], coordinates[k][1]);
					double param = delta(solid_point, cell_point, side_lengths[0]);
					for (size_t l = 0; l < value_size; l++)
					{
						local_fluid_values[local_to_global[cell_dofmap[k] + l]] += solid_values[i * value_size + l] * param * side_lengths[0] * side_lengths[1];
					}

				} // end loop inside the cell
			}
			// end adjacents loop
		}

		//////////////////  TODO : this part can use MPI_reduce directly  //////////////////////
		std::vector<double> fluid_values(global_fluid_size);
		std::vector<std::vector<double>> mpi_collect(dolfin::MPI::size(mesh->mpi_comm()));
		dolfin::MPI::all_gather(mesh->mpi_comm(), local_fluid_values, mpi_collect);
		for (size_t i = 0; i < fluid_values.size(); i++)
		{
			for (size_t j = 0; j < mpi_collect.size(); j++)
			{
				fluid_values[i] += mpi_collect[j][i];
			}
		}
		return fluid_values;
	}

	////////////////////////////////////////////
	// thses methods must not be modified!!  ///
	////////////////////////////////////////////
	double phi(double r)
	{
		r = fabs(r);
		if (r > 2)
			return 0;
		else
			return 0.25 * (1 + cos(FENICS_PI * r * 0.5));
	}
	double delta(Point p0, Point p1, double h = 0.0625)
	{
		double ret = 1.0;
		for (unsigned i = 0; i < 2; ++i)
		{
			double dx = p0.coordinates()[i] - p1.coordinates()[i];
			ret *= 1. / h * phi(dx / h);
		}
		return ret;
	}
	////////////////////////////////////////////
	// thses methods must not be modified!!  ///
	////////////////////////////////////////////
};