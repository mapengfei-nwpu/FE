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
	u.interpolate(*g);
	std::vector<double> body(U->dim());
	std::vector<std::array<double, 2>> body_coordinates(body.size() / u.value_size());

	auto local_dof_coordinates = U->tabulate_dof_coordinates();
	std::vector<std::vector<double>> mpi_dof_coordinates(dolfin::MPI::size(bb.mesh()->mpi_comm()));
	dolfin::MPI::all_gather(bb.mesh()->mpi_comm(), local_dof_coordinates, mpi_dof_coordinates);
	std::vector<double> dof_coordinates;

	for (size_t i = 0; i < mpi_dof_coordinates.size(); i++)
	{
		for (size_t j = 0; j < mpi_dof_coordinates[i].size(); j++)
		{
			dof_coordinates.push_back(mpi_dof_coordinates[i][j]);
		}
	}

	std::vector<double> local_values;
	u.vector()->get_local(local_values);
	std::vector<std::vector<double>> mpi_values(dolfin::MPI::size(bb.mesh()->mpi_comm()));
	dolfin::MPI::all_gather(bb.mesh()->mpi_comm(), local_values, mpi_values);
	std::vector<double> values;

	for (size_t i = 0; i < mpi_values.size(); i++)
	{
		for (size_t j = 0; j < mpi_values[i].size(); j++)
		{
			values.push_back(mpi_values[i][j]);
		}
	}

	if (dolfin::MPI::rank(bb.mesh()->mpi_comm()) == 0)
	{
		std::cout << "mpi rank: " << dolfin::MPI::rank(bb.mesh()->mpi_comm()) << std::endl;
		std::cout << "coor number: " << dof_coordinates.size() << std::endl;
		std::cout << "dof  number: " << values.size() << std::endl;

		for (size_t i = 0; i < values.size(); i++)
		{
			std::cout
				<< " NO. "
				<< i
				<< " . coordinate: "
				<< dof_coordinates[2 * i]
				<< " , "
				<< dof_coordinates[2 * i + 1]
				<< " . dof: "
				<< values[i]
				<< std::endl;
		}
		for (size_t i = 0; i < values.size(); i++)
		{
			if (i % 2 == 0)
			{
				std::cout
					<< dof_coordinates[2 * i] - values[i]
					<< std::endl;
			}
			else
			{
				std::cout
					<< dof_coordinates[2 * i + 1] - values[i]
					<< std::endl;
			}
		}
	}

	///DeltaInterplation::fun1(v, ba, body, body_coordinates);

	File file("xplusy.pvd");
	file << u;
}
