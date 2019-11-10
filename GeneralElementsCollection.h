#pragma once
#include <array>
#include <vector>
template<std::size_t N, std::size_t M>
class GeneralElementsCollection {
public:

	/// return geometric dimesnion. it decides how many
	/// axis every node has.
	const std::size_t geometric_dimension() { return N; }

	/// return the dimension of element.
	/// when it is not P1 element, the dimension
	/// of finite element is higher than
	/// the dimension of mesh element.
	const std::size_t element_dimension() { return M; }
	
	/// return element dimension which is the number of
	/// nodes every element contain
	const std::size_t element_number() { return _T.size(); }

	/// for finite element collection, it returns the number
	/// nodes. for mesh element collection, it returns the 
	/// total number of degree of freedom.
	std::size_t global_dimension() { return element_dimension() * element_number(); }

	/// Get a node
	///
	/// *Arguments*
	///     i (std::size_t)
	///         the index of the node you want to get.
	std::array<double, N> get_node(std::size_t i) { return _P[i]; };

	/// Get an element.
	///
	/// *Arguments*
	///     i (std::size_t)
	///         the index of the element you want to get.
	std::array<std::size_t, M> get_element(std::size_t i) { return _T[i]; };

	/// add a node to nodes collection.
	///
	/// *Arguments*
	///     node (std::array<double, N>)
	void add_node(std::array<double, N> node) { _P.push_back(node); }

	/// add an element to the tabulation. Here, "element" means the
	/// connectivity or tabulation of an element, which contians indices
	/// of all nodes on the element. no matter mesh nodes or dof nodes,
	/// finite element or mesh element.
	///
	/// *Arguments*
	///     element (std::array<std::size_t, M>)
	///         The tabulation of which the indices of the element consist
	void add_element(std::array<std::size_t, M> element) { _T.push_back(element); }

	/// Print the points and connectivity of the mesh.
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

	/// The collection which contains the coordinates of all points.
	std::vector<std::array<double, N>> _P;

	/// The tabulation for every element.
	std::vector<std::array<std::size_t, M>> _T;

};