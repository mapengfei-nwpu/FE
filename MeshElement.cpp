#include <vector>
#include <array>
#include <boost/multi_array.hpp>
class MeshElement {
	/// three points and every point point contains two coordinates.
	boost::multi_array<double, 2> coordinates;
};
