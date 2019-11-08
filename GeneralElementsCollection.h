#pragma once
#include <vector>
template<std::size_t N, std::size_t M>
class GeneralElementsCollection {
public:
	const std::size_t geometric_dimension() { return N; }
	const std::size_t fenitespace_dimension() { return M; }
private:
	std::vector<std::array<double, N>> _p;
	std::vector<std::array<std::size_t, M>> _T;

};