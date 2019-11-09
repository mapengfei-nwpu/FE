
#include <Eigen/Dense>
#include "MeshElementsCollection.h"
#include "GaussQuadrature.h"
#include "BasisFunction.h"
#include "MeshElement.h"
class LocalAssembler {
public:
	/// construct function.
	LocalAssembler(std::size_t index) {
		MeshElementsCollection mc;
		MeshElement me = mc.get_mesh_element(index);

		//FiniteElementsCollection fc;
		//FiniteElement fe = mc.get_finite_element(index);

		GaussQuadrature gq;
		auto egq = gq.quadrature_weights_nodes(me);

		BasisFunction b(me);
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

				/// calculate x_k, y_k, w_k for gauss quadrature.
				/// a(i,j) = \sum w_k(\psi^{(i)}_x(x_k,y_k)\psi^{(j)}_x(x_k,y_k)
				///        + \psi^{(i)}_y(x_k,y_k)\psi^{(j)}_y(x_k,y_k))
				for (std::size_t k = 0; k < egq.size(); ++k) {
					sum += egq[k].second * f(egq[k].first[0], egq[k].first[1]);
				}
				a(i, j) = sum;
			}
		}
	}


};