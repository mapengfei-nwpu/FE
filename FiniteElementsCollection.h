#pragma once
#include "GeneralElementsCollection.h"
#include "FiniteElement.h"
#include <array>
#include <vector>
#include <boost/multi_array.hpp>
class FiniteElementsCollection :public GeneralElementsCollection<2, 3> {
public:

	// FiniteElement get_finite_element(std::size_t i);
};

