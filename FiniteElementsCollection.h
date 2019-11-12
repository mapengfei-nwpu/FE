#pragma once
#include <array>
#include <map>
#include <iostream>
#include <vector>
#include <memory>
#include <boost/multi_array.hpp>

#include "FiniteElement.h"
#include "MeshElementsCollection.h"

class Edge {
public:
	Edge(size_t start, size_t end) : start(start), end(end) {
		if (end < start) { std::swap(end, start); }
	}
	size_t start;
	size_t end;
};

class FiniteElementsCollection :public GeneralElementsCollection<2, 6> {
public:
	/// construct from mesh and order.
	FiniteElementsCollection(MeshElementsCollection mesh, size_t order = 2) {
		if (order != 2) {
			std::cout << "Finite element order must be two." << std::endl;
		}
		mesh.print();
		/// generate index for midpoints on every edge.
		size_t mesh_element_size = mesh._T[0].size();
		size_t mesh_point_number = mesh._P.size();
		size_t finite_point_number = mesh_point_number;
		std::map<size_t, std::map<size_t, size_t>> midpoint;
		for (size_t i = 0; i < mesh._T.size(); i++)
		{
			Edge edge(mesh._T[i][mesh_element_size - 1], mesh._T[i][0]);
			if (midpoint[edge.start][edge.end] == 0) {
				midpoint[edge.start][edge.end] = finite_point_number++;
			}
			for (size_t j = 1; j < mesh_element_size; j++)
			{
				Edge edge(mesh._T[i][j - 1], mesh._T[i][j]);
				if (midpoint[edge.start][edge.end] == 0) {
					midpoint[edge.start][edge.end] = finite_point_number++;
				}
			}
		}
		for (auto i = midpoint.begin(); i != midpoint.end(); ++i) {
			auto secondpoint = i->second;
			for (auto j = secondpoint.begin(); j != secondpoint.end(); ++j) {
				std::cout << i->first << "  " << j->first << " " << j->second << std::endl;
			}
		}
		/// generate new _T which is the tabulation of finite elements.
		for (size_t i = 0; i < mesh._T.size(); i++)
		{
			std::array<size_t, 6> t;
			t[0] = mesh._T[i][0];
			Edge edge01(mesh._T[i][0], mesh._T[i][1]);
			t[1] = midpoint[edge01.start][edge01.end];
			t[2] = mesh._T[i][1];
			Edge edge12(mesh._T[i][1], mesh._T[i][2]);
			t[3] = midpoint[edge12.start][edge12.end];
			t[4] = mesh._T[i][2];
			Edge edge20(mesh._T[i][2], mesh._T[i][0]);
			t[5] = midpoint[edge20.start][edge20.end];
			_T.push_back(t);
		}

		/// generate new _P for the finite elements
		_P.resize(finite_point_number);
		for (size_t i = 0; i < mesh_point_number; i++)
		{
			_P[i] = mesh._P[i];
		}
		for (auto i = midpoint.begin(); i != midpoint.end(); ++i) {
			auto secondpoint = i->second;
			for (auto j = secondpoint.begin(); j != secondpoint.end(); ++j) {
				_P[j->second][0] = (_P[i->first][0] + _P[j->first][0]) / 2.0;
				_P[j->second][1] = (_P[i->first][1] + _P[j->first][1]) / 2.0;
			}
		}
		print();
	}

};