
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include "MeshElementsCollection.h"
#include "FiniteElementsCollection.h"
#include "GaussQuadrature.h"
#include "BasisFunction.h"
#include "MeshElement.h"
class Assembler {
public:

	/// construct function.
	Assembler(MeshElementsCollection mesh, FiniteElementsCollection functionspace) :
		mesh(mesh), functionspace(functionspace) {
		;
	}


	Eigen::MatrixXd matrixLocalAssembler(std::size_t index) {

		MeshElement mesh_element = mesh.get_mesh_element(index);

		GaussQuadrature quadrature;
		auto quadrature_coefficient = quadrature.quadrature_weights_nodes(mesh_element);

		BasisFunction b(mesh_element);
		std::size_t local_dim = functionspace.element_dimension();

		Eigen::MatrixXd a(local_dim, local_dim);

		/// I want to wrap basis functions into an array of "function object"
		/// Then, we can use for to loop over.
		/// a(i,j) = \int\psi^{(i)}_x\psi^{(j)}_x + \psi^{(i)}_y\psi^{(j)}_ydxdy
		auto basis = b.basis;
		auto jacob = std::abs(b.get_jacobian_determian());
		for (std::size_t i = 0; i < local_dim; ++i) {
			auto psi_i_x = basis[i]["dx"];
			auto psi_i_y = basis[i]["dy"];
			auto psi_i = basis[i]["original"];
			for (std::size_t j = 0; j < local_dim; ++j) {
				auto psi_j_x = basis[j]["dx"];
				auto psi_j_y = basis[j]["dy"];
				auto psi_j = basis[j]["original"];
				/// construct a lambda function to be integralled.
				auto f = [&](double x, double y) {
					return - (b.*psi_i_x)(x, y) * (b.*psi_j_x)(x, y) 
						- (b.*psi_i_y)(x, y) * (b.*psi_j_y)(x, y) 
						+ 25.0*(b.*psi_i)(x,y)*(b.*psi_j)(x,y);
				};
				double sum = 0;
				a(i, j) = 0;
				/// calculate x_k, y_k, w_k for gauss quadrature.
				/// a(i,j) = \sum w_k(\psi^{(i)}_x(x_k,y_k)\psi^{(j)}_x(x_k,y_k)
				///        + \psi^{(i)}_y(x_k,y_k)\psi^{(j)}_y(x_k,y_k))
				for (std::size_t k = 0; k < quadrature_coefficient.size(); ++k) {
					sum += quadrature_coefficient[k].second * f(quadrature_coefficient[k].first[0], quadrature_coefficient[k].first[1]);
				}
				a(i, j) = jacob*sum;
			}
		}
		return a;
	}
	Eigen::SparseMatrix<double> matrixGlobalAssembler() {
		std::size_t global_dim = functionspace.global_dimension();
		std::size_t element_number = functionspace.element_number();
		std::size_t element_dimension = functionspace.element_dimension();

		Eigen::SparseMatrix<double> A_sparse(global_dim, global_dim);         // default is column major
		A_sparse.reserve(Eigen::VectorXi::Constant(global_dim, 6));


		/// iterate over all elements
		for (std::size_t i = 0; i < element_number; i++)
		{
			auto dofmap = functionspace.get_element(i);
			auto local_A = matrixLocalAssembler(i);

			for (std::size_t j = 0; j < element_dimension; j++)
				for (std::size_t k = 0; k < element_dimension; k++)
				{
					A_sparse.coeffRef(dofmap[j], dofmap[k]) += local_A(j, k);
				}
					
		}
		A_sparse.makeCompressed();
		return A_sparse;
	}

	Eigen::VectorXd vectorLocalAssembler(std::size_t index) {

		MeshElement mesh_element = mesh.get_mesh_element(index);

		GaussQuadrature quadrature;
		auto quadrature_coefficient = quadrature.quadrature_weights_nodes(mesh_element);

		BasisFunction b(mesh_element);
		std::size_t local_dim = functionspace.element_dimension();

		Eigen::VectorXd rhs(local_dim);

		/// I want to wrap basis functions into an array of "function object"
		/// Then, we can use for to loop over.
		/// a(i,j) = \int\psi^{(i)}_x\psi^{(j)}_x + \psi^{(i)}_y\psi^{(j)}_ydxdy
		auto basis = b.basis;
		auto jacob = std::abs(b.get_jacobian_determian());
		for (std::size_t i = 0; i < local_dim; ++i) {
			auto psi_i = basis[i]["original"];

			/// construct a lambda function to be integralled.
			auto f = [&](double x, double y) {
				return (b.*psi_i)(x, y) * source(x, y);
			};
			double sum = 0;
			/// calculate x_k, y_k, w_k for gauss quadrature.
			for (std::size_t k = 0; k < quadrature_coefficient.size(); ++k) {
				sum += quadrature_coefficient[k].second * f(quadrature_coefficient[k].first[0], quadrature_coefficient[k].first[1]);
			}
			rhs(i) = (jacob * sum);
			
		}

		return rhs;
	}

	Eigen::VectorXd vectorGlobalAssembler() {

		std::size_t global_dim = functionspace.global_dimension();
		std::size_t element_number = functionspace.element_number();
		std::size_t element_dimension = functionspace.element_dimension();

		/// generate vector of the right side term.
		Eigen::VectorXd F = Eigen::VectorXd::Zero(global_dim);

		/// iterate over all elements
		for (std::size_t i = 0; i < element_number; i++)
		{
			auto dofmap = functionspace.get_element(i);
			auto local_F = vectorLocalAssembler(i);

			for (std::size_t j = 0; j < element_dimension; j++)
				F(dofmap[j]) += local_F(j);
		}
		return F;
	}

	double source(double x, double y) {
		return -1.0;
	}

	MeshElementsCollection mesh;
	FiniteElementsCollection functionspace;
};