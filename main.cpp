#include <iostream>
#include <dolfin.h>
#include "Poisson.h"
#include "BoxAdjacents.h"
#include "IBInterpolation.h"
using namespace dolfin;

class xplusy : public Expression
{
public:
	xplusy() : Expression(2) {}
	void eval(Array<double>& values, const Array<double>& x) const
	{
		values[0] = x[0];
		values[1] = x[1];
	}
};

std::vector<std::array<double, 2>> get_global_dof_coordinates(const Function & f) {
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

std::vector<double> get_global_dofs(const Function& f) {
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

int main()
{
	Point p0(0, 0, 0);
	Point p1(1, 1, 0);
	BoxAdjacents ba({ p0, p1 }, { 16, 16 }, CellType::Type::quadrilateral);
	auto V = std::make_shared<Poisson::FunctionSpace>(ba.mesh());
	auto g = std::make_shared<xplusy>();
	Function v(V);
	v.interpolate(*g);

	Point p2(0.25, 0.25, 0);
	Point p3(0.75, 0.75, 0);
	BoxAdjacents bb({ p2, p3 }, { 16, 16 }, CellType::Type::quadrilateral);
	auto U = std::make_shared<Poisson::FunctionSpace>(bb.mesh());
	Function u(U);


	auto dof_coordinates = get_global_dof_coordinates(u);

	std::vector<double> values(dof_coordinates.size() * 2);
	DeltaInterplation::fun1(v, ba, values, dof_coordinates);
	DeltaInterplation::fun0(u, values);

	File file("xplusy.pvd");
	file << u;
}
