#include <Eigen/Dense>
#include "GaussQuadrature.h"
#include "BasisFunction.h"
#include "MeshElement.h"
class LocalAssembler {
public:
	/// construct function.
	LocalAssembler(std::size_t index) {
		MeshElementsCollection mc;
		MeshElement me = mc.get_mesh_element(index);
		
		FiniteElementsCollection fc;
		FiniteElement fe = mc.get_finite_element(index);
		
		GaussQuadrature gq;
		auto egq = gq.quadrature_weights_nodes(me);

		BasisFunction b(me);
		std::size_t local_dim = fc.element_dimension();

		Eigen::MatrixXd a(local_dim, local_dim);
		
		/// I want to wrap basis functions into an array of "function object"
		/// Then, we can use for to loop over.
		/// a(i,j) = \int\psi^{(i)}_x\psi^{(j)}_x + \psi^{(i)}_y\psi^{(j)}_ydxdy
		for (std::size_t i = 0; i < local_dim; ++i) {
			for (std::size_t j = 0; j < local_dim; ++j) {

				/// construct a lambda function to be integralled.
				auto f = [](double x, double y) {return x + y; };
				double sum = 0;

				/// calculate x_k, y_k, w_k for gauss quadrature.
				/// a(i,j) = \sum w_k(\psi^{(i)}_x(x_k,y_k)\psi^{(j)}_x(x_k,y_k)
				///        + \psi^{(i)}_y(x_k,y_k)\psi^{(j)}_y(x_k,y_k))
				for (std::size_t k = 0; k < egq.size(); ++k){
					sum += egq[k].second * f(egq[k].first[0], egq[k].first[1]);
				}
				a(i, j) = sum;
			}
		}
	}

private:

	
}