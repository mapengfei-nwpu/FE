#include <iostream>
#include "MeshElementsCollection.h"
int main() {
	/********** test 1 ************
	MeshElementsCollection a;
	a.print();
	*******************************/
	///* test 2
	MeshElementsCollection a;
	auto b = a.get_mesh_element(0);
	a.print();
	b.print();
	//*/

	std::cout << "hello world!" << std::endl;
}