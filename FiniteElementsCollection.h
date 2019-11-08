#pragma once

class FiniteElementsCollection :public GeneralElementsCollection<3, 4> {
public:
	FiniteElement get_finite_element(std::size_t i);
};

