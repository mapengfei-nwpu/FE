#include <iostream>
#include <dolfin.h>
using namespace dolfin;

class BoxAdjacents
{
public:
    BoxAdjacents(std::array<dolfin::Point, 2> points, std::vector<size_t> dims, CellType::Type cell_type)
    {

        // toplogy dimesion
        top_dim = dims.size();

        // generate mesh
        if (cell_type == CellType::Type::hexahedron)
        {
            dolfin_assert(top_dims == 3);
            _mesh = BoxMesh::create(points, {dims[0], dims[1], dims[2]}, cell_type);
        }
        else if (cell_type == CellType::Type::quadrilateral)
        {
            dolfin_assert(top_dims == 2);
            _mesh = RectangleMesh::create(points, {dims[0], dims[1]}, cell_type);
        }
        else
        {
            dolfin_error("BoxMesh.h",
                         "generate uniform mesh",
                         "Wrong cell type '%d'", cell_type);
        }

        // check again.
        if (top_dim != 2 && top_dim != 3)
            dolfin_error("the size of dims must be 2 & 3.", ".", ".");

        nx = dims[0];
        ny = dims[1];
        nz = 0;

        x0 = points[0].x();
        x1 = points[1].x();
        y0 = points[0].y();
        y1 = points[1].y();
        z0 = 0.0;
        z1 = 0.0;

        if (top_dim == 3)
        {
            nz = dims[2];
            z0 = points[0].z();
            z1 = points[1].z();
        }
    }
    size_t hash(dolfin::Point point)
    {
        if (top_dim != 2 && top_dim != 3)
            dolfin_error("the size of dims must be 2 and 3.", ".", ".");

        double x = point.x();
        double y = point.y();
        double z = point.z();

        double dx = (x1 - x0) / static_cast<double>(nx);
        double dy = (y1 - y0) / static_cast<double>(ny);
        double dz = 0.0;

        size_t i = static_cast<size_t>((x - x0) / dx);
        size_t j = static_cast<size_t>((y - y0) / dy);
        size_t k = 0;

        if (top_dim == 3)
        {
            dz = (z1 - z0) / static_cast<double>(nz);
            k = static_cast<size_t>((z - z0) / dz);
        }

        return k * nx * ny + j * nx + i;
    }
    bool check()
    {
        for (CellIterator cell(_mesh); !cell.end(); ++cell)
        {
            std::cout << "global:" << cell->global_index() << ",hash:";
            std::cout << hash(cell->midpoint()) << std::endl;
            if (cell->global_index() != hash(cell->midpoint()))
            {
                dolfin_error("mesh is not consistent.", " ", " ");
                return false;
            }
        }
        return true;
    }

    double x0, x1, y0, y1, z0, z1;
    size_t nx, ny, nz;
    size_t top_dim;
    Mesh _mesh;
};

int main()
{

    Point p0(0, 0, 0);
    Point p1(1, 1, 1);
    BoxAdjacents ba({p0, p1}, {8, 8}, CellType::Type::quadrilateral);
    //BoxAdjacents ba({p0,p1}, {8,8,8}, CellType::Type::hexahedron);
    ba.check();
}
