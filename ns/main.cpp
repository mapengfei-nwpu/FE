
// Begin demo
#include "Circle.h"
#include "BoxAdjacents.h"
#include "IBInterpolation.h"

#include <dolfin.h>
#include "TentativeVelocity.h"
#include "PressureUpdate.h"
#include "VelocityUpdate.h"


using namespace dolfin;

// Define noslip domain
class NoslipDomain : public SubDomain
{
  bool inside(const Array<double> &x, bool on_boundary) const
  {
    return near(x[1], 0) || near(x[1], 0.41);
  }
};

// Define inflow domain
class InflowDomain : public SubDomain
{
  bool inside(const Array<double> &x, bool on_boundary) const
  {
    return near(x[0], 0);
  }
};

// Define inflow domain
class OutflowDomain : public SubDomain
{
  bool inside(const Array<double> &x, bool on_boundary) const
  {
    return near(x[0], 2.2);
  }
};

// Define pressure boundary value at inflow
class InflowVelocity : public Expression
{
public:
  // Constructor
  InflowVelocity() : Expression(2) {}

  // Evaluate pressure at inflow
  void eval(Array<double> &values, const Array<double> &x) const
  {
    values[0] = 4.0 * 1.5 * x[1] * (0.41 - x[1]) / pow(0.41, 2);
    values[1] = 0.0;
  }
};

int main()
{
  // Create chanel mesh
  Point point0(0, 0, 0);
  Point point1(2.2, 0.41, 0);
  BoxAdjacents ba({point0, point1}, {220, 41}, CellType::Type::quadrilateral);
  DeltaInterplation interpolation(ba);

  // Create circle mesh
  auto circle = std::make_shared<Mesh>("../circle.xml.gz");
  auto U = std::make_shared<Circle::FunctionSpace>(circle);
  auto body_velocity = std::make_shared<Function>(U);
  auto body_force = std::make_shared<Function>(U);

  // Create function spaces
  auto V = std::make_shared<VelocityUpdate::FunctionSpace>(ba.mesh());
  auto Q = std::make_shared<PressureUpdate::FunctionSpace>(ba.mesh());

  // Set parameter values
  double dt = 0.0005;
  double T = 5;

  // Define values for boundary conditions
  auto v_in = std::make_shared<InflowVelocity>();
  auto zero = std::make_shared<Constant>(0.0);
  auto zero_vector = std::make_shared<Constant>(0.0, 0.0);

  // Define subdomains for boundary conditions
  auto noslip_domain = std::make_shared<NoslipDomain>();
  auto inflow_domain = std::make_shared<InflowDomain>();
  auto outflow_domain = std::make_shared<OutflowDomain>();

  // Define boundary conditions
  DirichletBC noslip(V, zero_vector, noslip_domain);
  DirichletBC inflow(V, v_in, inflow_domain);
  DirichletBC outflow(Q, zero, outflow_domain);
  std::vector<DirichletBC *> bcu = {{&noslip, &inflow}};
  std::vector<DirichletBC *> bcp = {&outflow};

  // Create functions
  auto u0 = std::make_shared<Function>(V);
  auto u1 = std::make_shared<Function>(V);
  auto p1 = std::make_shared<Function>(Q);

  // Create coefficients
  auto k = std::make_shared<Constant>(dt);
  auto f = std::make_shared<Function>(V);
  auto fluid_force = std::make_shared<Function>(V);

  // Create forms
  TentativeVelocity::BilinearForm a1(V, V);
  TentativeVelocity::LinearForm L1(V);
  PressureUpdate::BilinearForm a2(Q, Q);
  PressureUpdate::LinearForm L2(Q);
  VelocityUpdate::BilinearForm a3(V, V);
  VelocityUpdate::LinearForm L3(V);

  // Set coefficients
  a1.k = k;
  L1.k = k;
  L1.u0 = u0;
  L1.f = f;
  L2.k = k;
  L2.u1 = u1;
  L3.k = k;
  L3.u1 = u1;
  L3.p1 = p1;

  // Assemble matrices
  Matrix A1, A2, A3;
  assemble(A1, a1);
  assemble(A2, a2);
  assemble(A3, a3);

  // Create vectors
  Vector b1, b2, b3;

  // Use amg preconditioner if available
  const std::string prec(has_krylov_solver_preconditioner("amg") ? "amg" : "default");

  // Create files for storing solution
  File ufile("results/velocity.pvd");
  File pfile("results/pressure.pvd");
  File ffile("results/force.pvd");

  // Time-stepping
  double t = dt;
  while (t < T + DOLFIN_EPS)
  {
    // Interpolate velocity to solid.
    interpolation.fluid_to_solid(*u1, *body_velocity);

    // calculate body force.
    *body_force = FunctionAXPY(body_velocity,-1);

    // interpolate force into fluid.
    interpolation.solid_to_fluid(*f, *body_force);

    // assign force term
    L1.f = f;

    // Compute tentative velocity step
    begin("Computing tentative velocity");
    assemble(b1, L1);
    for (std::size_t i = 0; i < bcu.size(); i++)
      bcu[i]->apply(A1, b1);
    solve(A1, *u1->vector(), b1, "gmres", "default");
    end();

    // Pressure correction
    begin("Computing pressure correction");
    assemble(b2, L2);
    for (std::size_t i = 0; i < bcp.size(); i++)
    {
      bcp[i]->apply(A2, b2);
      bcp[i]->apply(*p1->vector());
    }
    solve(A2, *p1->vector(), b2, "bicgstab", prec);
    end();

    // Velocity correction
    begin("Computing velocity correction");
    assemble(b3, L3);
    for (std::size_t i = 0; i < bcu.size(); i++)
      bcu[i]->apply(A3, b3);
    solve(A3, *u1->vector(), b3, "gmres", "default");
    end();

    // Save to file
    ufile << *u1;
    pfile << *p1;
    ffile << *f;

    // Move to next time step
    *u0 = *u1;
    t += dt;
    cout << "t = " << t << endl;
  }

  return 0;
}
