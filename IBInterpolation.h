#include <dolfin.h>
using namespace dolfin;
class DeltaInterplation
{
	/// collect all body vectors on every processes and
	/// assign them to a function.
	void fun0(Function& u, std::vector<double>& body)
	{
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
	void fun1(Function& v, BoxAdjacents& um, std::vector<double>& body,
		std::vector<std::array<double, 2>> & body_coordinates)
	{
		/// the meshes of v and um should be the same.
		/// TODO : compare these two meshed.

		/// smart shortcut
		auto rank = dolfin::MPI::rank(v.function_space()->mesh()->mpi_comm());
		auto mesh = v.function_space()->mesh();     // pointer to a mesh
		auto dofmap = v.function_space()->dofmap(); // pointer to a dofmap

		/// iterate every body coordinate.
		for (size_t i = 0; i < body.size(); i++)
		{
			Point body_coordinate(body_coordinates[i][0], body_coordinates[i][0]);
			auto adjacents = um.get_adjacents(body_coordinate);
			for (size_t j = 0; j < adjacents.size(); j++)
			{
				/// Cell constructor take local index to initial.
				Cell cell(*mesh, adjacents[j]);
				boost::multi_array<double, 2> coordinates; /// coordinates of the cell
				std::vector<double> coordinate_dofs;       /// coordinate_dofs of the cell
				cell.get_coordinate_dofs(coordinate_dofs);
				auto _element = v.function_space()->element(); /// element of the function space
				/// get coordinates and coordinate_dofs
				_element->tabulate_dof_coordinates(coordinates, coordinate_dofs, cell);
				/// local index of the cell is needed rather than global index.
				auto cell_dofmap = dofmap->cell_dofs(cell.index());
				for (size_t k = 0; k < cell_dofmap.size(); k++)
				{
					/// multiply with delta function.
					/// four components are needed here: body_coordinate,
					/// coordinates, coordinate_dofs, cel_dofmap
					Point on_cell(coordinates[k][0], coordinates[k][1]);
					body[cell_dofmap[i]] += delta(body_coordinate, on_cell) * coordinate_dofs[k];
				} // end loop inside the cell
			}
			// end adjacents loop
		}
		// end body cycle

	}
	// end

	double phi(double r)
	{
		r = fabs(r);
		if (r > 2) return 0;
		else return 0.25 * (1 + cos(FENICS_PI * r * 0.5));
	}
	double delta(Point p0, Point p1, double h = 0.02)
	{
		double ret = 1.0;
		for (unsigned i = 0; i < 2; ++i)
		{
			double dx = p0.coordinates()[i] - p1.coordinates()[i];
			ret *= 1. / h * phi(dx / h);
		}
		return ret;
	}
};