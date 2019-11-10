
#include <Eigen/Dense>
#include "MeshElementsCollection.h"
#include "FiniteElementsCollection.h"
#include "GaussQuadrature.h"
#include "BasisFunction.h"
#include "MeshElement.h"
class Assembler {
public:
	Assembler(MeshElementsCollection mesh, FiniteElementsCollection functionspace) :
		mesh(mesh), functionspace(functionspace) {
		;
	}

	/// construct function.
	Eigen::MatrixXd localAssembler(std::size_t index) {

		MeshElement mesh_element = mesh.get_mesh_element(index);

		GaussQuadrature quadrature;
		auto quadrature_coefficient = quadrature.quadrature_weights_nodes(mesh_element);

		BasisFunction b(mesh_element);
		std::size_t local_dim = 3;// fc.element_dimension();

		Eigen::MatrixXd a(local_dim, local_dim);

		/// I want to wrap basis functions into an array of "function object"
		/// Then, we can use for to loop over.
		/// a(i,j) = \int\psi^{(i)}_x\psi^{(j)}_x + \psi^{(i)}_y\psi^{(j)}_ydxdy
		auto basis = b.basis;
		for (std::size_t i = 0; i < local_dim; ++i) {
			auto psi_i_x = basis[i]["dx"];
			auto psi_i_y = basis[i]["dy"];
			for (std::size_t j = 0; j < local_dim; ++j) {
				auto psi_j_x = basis[j]["dx"];
				auto psi_j_y = basis[j]["dy"];
				/// construct a lambda function to be integralled.
				auto f = [&](double x, double y) {
					return (b.*psi_i_x)(x, y) * (b.*psi_j_x)(x, y) +
						(b.*psi_i_y)(x, y) * (b.*psi_j_y)(x, y);
				};
				double sum = 0;
				a(i, j) = 0;
				/// calculate x_k, y_k, w_k for gauss quadrature.
				/// a(i,j) = \sum w_k(\psi^{(i)}_x(x_k,y_k)\psi^{(j)}_x(x_k,y_k)
				///        + \psi^{(i)}_y(x_k,y_k)\psi^{(j)}_y(x_k,y_k))
				for (std::size_t k = 0; k < quadrature_coefficient.size(); ++k) {
					sum += quadrature_coefficient[k].second * f(quadrature_coefficient[k].first[0], quadrature_coefficient[k].first[1]);
				}
				a(i, j) = sum;
			}
		}
		return a;
	}
	Eigen::MatrixXd globalAssembler() {
		std::size_t dof_dim = functionspace.global_dimension();
		std::size_t element_number = functionspace.element_number();
		std::size_t element_dimension = functionspace.element_dimension();
		Eigen::MatrixXd A(dof_dim, dof_dim);
		
		/// iterate over all elements
		for (std::size_t i = 0; i < element_number; i++)
		{
			auto dofmap = functionspace.get_element(i);
			auto local_A = localAssembler(i);

			for (std::size_t j = 0; j < element_dimension; j++)
				for (std::size_t k = 0; k < element_dimension; k++)
					A(dofmap[i], dofmap[j]) = local_A(i, j);
		}
		return A;
	}
	MeshElementsCollection mesh;
	FiniteElementsCollection functionspace;
};