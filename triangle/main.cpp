#include <dolfin.h>
#include "IBMesh.h"
int main()
{
    Point point0(0.0, 0.0);
    Point point1(1.0, 1.0);
    IBMesh ib_mesh({point0, point1}, {32, 32});
    ib_mesh.check();


}