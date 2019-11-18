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
        index_mesh();
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
            auto point = cell->midpoint();
            std::cout << "global index: "
                      << cell->global_index()
                      << " or "
                      << hash(point)
                      << ". local index: "
                      << cell->index()
                      << " or "
                      << map(hash(point))[1]
                      << ". mpi rank: "
                      << map(hash(point))[0]
                      << ". number of adjacents: "
                      // << get_adjacents(point).size()
                      << std::endl;
            if (cell->global_index() != hash(cell->midpoint()))
            {
                return false;
            }
        }
        return true;
    }

    // point -> global index of adjacent cells -> local index of adjacent cells
    std::vector<size_t> get_adjacents(Point point)
    {
        /// only finished for 2d occassion.
        std::vector<size_t> adjacents;
        size_t global_index = hash(point);
        size_t a[9] = {global_index - nx - 1, global_index - nx, global_index - nx + 1,
                       global_index - 1, global_index, global_index + 1,
                       global_index + nx - 1, global_index + nx, global_index + nx + 1};
        for (size_t i = 0; i < 9; i++)
        {
            if (global_index >= 0 && global_index <= global_map.size())
            {
                if (map(a[i])[0] == dolfin::MPI::rank(_mesh.mpi_comm()))
                {
                    adjacents.push_back(map(a[i])[0]);
                }
            }
        }
        return adjacents;
    }

    // global index and cell center
    std::array<size_t, 2> map(size_t i)
    {
        return global_map[i];
    }

private:
    void index_mesh()
    {
        // The local map local to global
        std::vector<size_t> local_map;
        for (dolfin::CellIterator e(_mesh); !e.end(); ++e)
        {
            auto center = e->midpoint();
            local_map.push_back(e->global_index());
            local_map.push_back(dolfin::MPI::rank(_mesh.mpi_comm()));
            local_map.push_back(e->index());
        }
        // send local map to every peocess.
        std::vector<std::vector<size_t>> mpi_collect(dolfin::MPI::size(_mesh.mpi_comm()));
        dolfin::MPI::all_gather(_mesh.mpi_comm(), local_map, mpi_collect);
        // alloc memory for global map.
        auto num_cell_global = _mesh.num_entities_global(2);
        global_map.resize(num_cell_global);
        for (auto iter = mpi_collect.cbegin(); iter != mpi_collect.cend(); iter++)
        {
            for (auto jter = iter->begin(); jter != iter->cend();)
            {
                size_t cell_index = *jter;
                jter++;
                global_map[cell_index][0] = *jter; /// mpi_rank;size_t
                jter++;
                global_map[cell_index][1] = *jter; /// cell hash;
                jter++;
            }
        }
    }

    double x0, x1, y0, y1, z0, z1;
    size_t nx, ny, nz;
    size_t top_dim;
    // The map of global index to hash index for cells.
    std::vector<std::array<size_t, 2>> global_map;
    Mesh _mesh;
};

class DeltaInterplation
{
    /// collect all body vectors on every processes and
    /// assign them to a function.
    void fun0(dolfin::Function &u, std::vector<double> &body)
    {
        /// TODO : all reduce body.
        /// dolfin::MPI::all_reduce(body)
        auto size = u.vector()->local_size();
        auto start = u.vector()->local_range().first;
        std::vector<double> values(size);
        for (size_t i = 0; i < size; ++i)
            values[i] = body[i + start];
        u.vector()->set_local(values);
    }

    /// distribute v to vector body on every process.
    void fun1(dolfin::Function &v, BoxAdjacents &um, std::vector<double> &body,
              std::vector<std::array<double, 2>> &body_coordinates)
    {
        /// the meshes of v and um should be the same.
        /// TODO : compare these two meshed.

        /// smart shortcut
        auto rank = dolfin::MPI::rank(v.function_space()->mesh()->mpi_comm());
        auto mesh = v.function_space()->mesh();     // pointer to a mesh
        auto dofmap = v.function_space()->dofmap(); // pointer to a dofmap

        /// iterate every body coordinate.
        for (size_t i = 0; i < body.size(); i++)
        {
            Point body_coordinate(body_coordinates[i][0], body_coordinates[i][0]);
            auto adjacents = um.get_adjacents(body_coordinate);
            for (size_t j = 0; j < adjacents.size(); j++)
            {
                /// Cell constructor take local index to initial.
                Cell cell(*mesh, adjacents[j]);
                boost::multi_array<double, 2> coordinates; /// coordinates of the cell
                std::vector<double> coordinate_dofs;       /// coordinate_dofs of the cell
                cell.get_coordinate_dofs(coordinate_dofs);
                auto _element = v.function_space()->element(); /// element of the function space
                /// get coordinates and coordinate_dofs
                _element->tabulate_dof_coordinates(coordinates, coordinate_dofs, cell);
                /// local index of the cell is needed rather than global index.
                auto cell_dofmap = dofmap->cell_dofs(cell.index());
                for (size_t k = 0; k < cell_dofmap.size(); k++)
                {
                    /// multiply with delta function.
                    /// four components are needed here: body_coordinate,
                    /// coordinates, coordinate_dofs, cel_dofmap
                    Point on_cell(coordinates[k][0], coordinates[k][1]);
                    body[cell_dofmap[i]] += delta(body_coordinate, on_cell) * coordinate_dofs[k];
                } // end loop inside the cell
            }
            // end adjacents loop
        }
        // end body cycle

    }
    // end

    double phi(double r)
    {
        r = fabs(r);
        if (r > 2) return 0;
        else return 0.25 * (1 + cos(FENICS_PI * r * 0.5));
    }
    double delta(Point p0, Point p1, double h = 0.02)
    {
        double ret = 1.0;
        for (unsigned i = 0; i < 2; ++i)
        {
            double dx = p0.coordinates()[i] - p1.coordinates()[i];
            ret *= 1. / h * phi(dx / h);
        }
        return ret;
    }
};

int main()
{

    Point p0(0, 0, 0);
    Point p1(1, 1, 0);
    BoxAdjacents ba({p0, p1}, {8, 8}, CellType::Type::quadrilateral);
    std::cout << "index:"
              << ba.hash(Point(0.5, 0.5))
              << "adjacent number"
              << ba.get_adjacents(Point(0.5, 0.5)).size()
              << std::endl;
    auto a = 1.0;

    // BoxAdjacents ba({p0,p1}, {8,8,8}, CellType::Type::hexahedron);
    // ba.check();
}
