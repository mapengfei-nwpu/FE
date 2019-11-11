#pragma once
#include <array>
#include <map>
#include <iostream>
#include <vector>
#include <memory>
#include <boost/multi_array.hpp>

#include "FiniteElement.h"
#include "MeshElementsCollection.h"

class FiniteElementsCollection :public GeneralElementsCollection<2, 6> {
public:
	/// construct from mesh and order.
	FiniteElementsCollection(MeshElementsCollection mec, std::size_t order = 1) {
		if (order == 2) {
			/// generate the index of finite element nodes.
			std::size_t index = mec._P.size();
			std::map<std::size_t, std::map<std::size_t, std::size_t>> midpoint;
			for (size_t i = 0; i < mec._T.size(); i++)
			{
				_T[i][_T[0].size() - 1], _T[i][0];
				for (size_t j = 1; j < _T[0].size(); j++)
				{
					std::size_t start, end;
					if (_T[i][j - 1] > _T[i][j]) {
						start = _T[i][j];
						end = _T[i][j - 1];
					}
					else {
						end = _T[i][j];
						start = _T[i][j - 1];
					}
					midpoint[start][end] = index++;
				}
			}
			/// generate new _T;
			for (size_t i = 0; i < mec._T.size(); i++)
			{
				std::array<std::size_t, 6> t;
				std::size_t start, end;
				t[0] = mec._T[i][0];
				/// start = 
				/// end   =
				t[1] = midpoint[start][end];
				t[2] = mec._T[i][1];
				/// start = 
				/// end   =
				t[3] = midpoint[start][end];
				t[4] = mec._T[i][2];
				/// start = 
				/// end   =
				t[5] = midpoint[start][end];
				_T.push_back(t);
			}
		}
		else {
			std::cout << "Finite element order must be one." << std::endl;
		}
	}
private:
};