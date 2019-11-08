#include <functional>
#include <boost/multi_array.hpp>

#include "BasisFunction.h"

std::function<double(double, double)> phi_x;
auto mod = [](int a, int b) {return a % b;};
void compute_jacobian_determinants() { ; }

auto psi = [](int a) {return a; };
auto psi_derivative_x  = [](int a) {return a; };
auto psi_derivative_y  = [](int a) {return a; };

int mammin() {

	return 0;
}


















