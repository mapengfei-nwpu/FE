#pragma once
#include <array>
#include <vector>
template<std::size_t N, std::size_t M>
class GeneralElementsCollection {
public:
	const std::size_t geometric_dimension() { return N; }
	const std::size_t finitespace_dimension() { return M; }
	void add_node(std::array<double, N> node) { _P.push_back(node); }
	std::array<double,N> get_node(double i) { return _P[i]; };
	std::array<std::size_t, M> get_element(std::size_t i) { return _T[i]; };
	void add_element(std::array<std::size_t, M> element) { _T.push_back(element); }
	void print() {
		std::cout << "Points:" << std::endl;
		for (std::size_t i = 0; i < _P.size(); ++i) {
			std::cout << _P[i][0] << "," << _P[i][1] << std::endl;
		}
		std::cout << "Elements:" << std::endl;
		for (std::size_t i = 0; i < _T.size(); ++i) {
			std::cout << _T[i][0] << "," << _T[i][1] << "," << _T[i][2] << std::endl;
		}
	}
private:
	std::vector<std::array<double, N>> _P;
	std::vector<std::array<std::size_t, M>> _T;
};