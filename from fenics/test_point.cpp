#include <iostream>
#include "Point.h"
int test_point() {
	dolfin::Point a(1, 2, 3);
	std::cout << a << std::endl;
	return 0;
}