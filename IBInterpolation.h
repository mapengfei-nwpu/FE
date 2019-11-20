#include <dolfin.h>
#include <numeric>
using namespace dolfin;
class DeltaInterplation
{
public:
	/// collect all body vectors on every processes and
	/// assign them to a function.
	static void fun0(Function& u, std::vector<double>& body)
	{
		auto size = u.vector()->local_size();
		auto start = u.vector()->local_range().first;
		std::vector<double> values(size);
		for (size_t i = 0; i < size; ++i)
			values[i] = body[i + start];
		u.vector()->set_local(values);
	}

	static void fun1(Function& v, BoxAdjacents& um, std::vector<double>& body,
		std::vector<std::array<double, 2>> & body_coordinates)
	{
		/// the meshes of v and um should be the same.
		/// TODO : compare two meshes
		std::cout << v.vector()->local_size() << std::endl;
		/// smart shortcut
		auto rank = dolfin::MPI::rank(v.function_space()->mesh()->mpi_comm());
		auto mesh = v.function_space()->mesh();		// pointer to a mesh
		auto dofmap = v.function_space()->dofmap(); // pointer to a dofmap
		size_t n = 0;
		/// iterate every body coordinate.
		for (size_t i = 0; i < body.size(); i++)
		{
			body[i] = 0.0;
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
				for (size_t k = 0; k < cell_dofmap.size(); k++)
				{
					Point on_cell(coordinates[k][0], coordinates[k][1]);
					body[i] += delta(body_coordinate, on_cell) * (*(v.vector()))[cell_dofmap[k]] / 16.0 / 16.0;

					//////////////////////////// WATCH OUT!//////////////////////////////////
					///                                                                   ///
					/// vector in a Function is initialized with a dofmap. this dofmap    ///
					/// contains many informations especially the layout of the vector.   ///
					/// Besides, dofmap tells which is the ghost entries of a element.    ///
					/// see Function::init_vector()                                       ///
					/// and https://fenicsproject.discourse.group/t/it-seems-that         ///
					/// -local-size-didnt-return-real-local-size-of-a-vector/1929         ///    
					///                                                                   ///
					/////////////////////////////////////////////////////////////////////////
					/*
					if (cell_dofmap[k] > v.vector()->local_size() && rank == 0 )
					{
						std::cout << "index: " << cell_dofmap[k] << "vector size: " << v.vector()->local_size() << std::endl;
					}*/

				} // end loop inside the cell
			}
			// end adjacents loop
		}
		// end body cycle
		/* output body before mpi all_gather.
		for (size_t i = 0; i < body.size(); ++i)
		{
			if (rank == 0)
				std::cout << body[i] << std::endl;
		}*/
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
		/* output after all_gather.
		for (size_t i = 0; i < body.size(); ++i)
		{
			if (rank == 0)
				std::cout << body[i] << std::endl;
		}
		*/
	}

	////////////////////////////////////////////
	// thses methods must not be modified!!  ///
	////////////////////////////////////////////
	static double phi(double r)
	{
		r = fabs(r);
		if (r > 2)
			return 0;
		else
			return 0.25 * (1 + cos(FENICS_PI * r * 0.5));
	}
	static double delta(Point p0, Point p1, double h = 0.125)
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