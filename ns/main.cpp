
// Begin demo
#include <dolfin/geometry/SimplexQuadrature.h>

#include "Circle.h"
#include "BoxAdjacents.h"
#include "IBInterpolation.h"

#include <dolfin.h>
#include "TentativeVelocity.h"
#include "PressureUpdate.h"
#include "VelocityUpdate.h"
#include "ElasticStructure.h"

using namespace dolfin;

// Define noslip domain
class NoslipDomain : public SubDomain
{
  bool inside(const Array<double> &x, bool on_boundary) const
  {
    return near(x[0], 0) || near(x[0], 1.0) || near(x[1], 0);
  }
};

// Define inflow domain
class InflowDomain : public SubDomain
{
  bool inside(const Array<double> &x, bool on_boundary) const
  {
    return near(x[1], 1);
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
    values[0] = 1.0;
    values[1] = 0.0;
  }
};

int main()
{
  // Create chanel mesh
  Point point0(0, 0, 0);
  Point point1(1.0, 1.0, 0);
  BoxAdjacents ba({point0, point1}, {64, 64}, CellType::Type::quadrilateral);
  DeltaInterplation interpolation(ba);

  // Create circle mesh
  auto circle = std::make_shared<Mesh>("../circle.xml.gz");
  auto U = std::make_shared<ElasticStructure::FunctionSpace>(circle);
  auto body_velocity = std::make_shared<Function>(U);
  auto body_force = std::make_shared<Function>(U);
  auto body_disp = std::make_shared<Function>(U);

  // Create function spaces
  auto V = std::make_shared<VelocityUpdate::FunctionSpace>(ba.mesh());
  auto Q = std::make_shared<PressureUpdate::FunctionSpace>(ba.mesh());

  // Set parameter values
  double dt = 0.0001;
  double T = 10;

  // Define values for boundary conditions
  auto v_in = std::make_shared<InflowVelocity>();
  auto zero = std::make_shared<Constant>(0.0);
  auto zero_vector = std::make_shared<Constant>(0.0, 0.0);

  // Define subdomains for boundary conditions
  auto noslip_domain = std::make_shared<NoslipDomain>();
  auto inflow_domain = std::make_shared<InflowDomain>();

  // Define boundary conditions
  DirichletBC noslip(V, zero_vector, noslip_domain);
  DirichletBC inflow(V, v_in, inflow_domain);
  std::vector<DirichletBC *> bcu = {{&inflow, &noslip}};
  std::vector<DirichletBC *> bcp = {};

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
  ElasticStructure::BilinearForm a4(U, U);
  ElasticStructure::LinearForm L4(U);
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
  L4.u = body_velocity;

  // Assemble matrices
  Matrix A1, A2, A3, A4;
  assemble(A1, a1);
  assemble(A2, a2);
  assemble(A3, a3);
  assemble(A4, a4);

  // Create vectors
  Vector b1, b2, b3, b4;

  // Create files for storing solution
  File ufile("results/velocity.pvd");
  File pfile("results/pressure.pvd");
  File ffile("results/force.pvd");
  File bfile("results/body.pvd");

  // Time-stepping
  double t = dt;
  while (t < T + DOLFIN_EPS)
  {
    // Interpolate velocity to solid.
    //interpolation.fluid_to_solid(*u1, *body_velocity);

    // calculate body force.
    //*body_force = FunctionAXPY(body_velocity,-1);

    // interpolate force into fluid.
    //interpolation.solid_to_fluid(*f, *body_force);

    // assign force term
    //L1.f = f;

    // Compute tentative velocity step
    begin("Computing tentative velocity");
    assemble(b1, L1);
    for (std::size_t i = 0; i < bcu.size(); i++)
      bcu[i]->apply(A1, b1);
    solve(A1, *u1->vector(), b1, "bicgstab", "hypre_amg");
    end();

    // Pressure correction
    begin("Computing pressure correction");
    assemble(b2, L2);
    for (std::size_t i = 0; i < bcp.size(); i++)
    {
      bcp[i]->apply(A2, b2);
      bcp[i]->apply(*p1->vector());
    }
    solve(A2, *p1->vector(), b2, "bicgstab", "hypre_amg");
    end();

    // Velocity correction
    begin("Computing velocity correction");
    assemble(b3, L3);
    for (std::size_t i = 0; i < bcu.size(); i++)
      bcu[i]->apply(A3, b3);
    solve(A3, *u1->vector(), b3, "cg", "sor");
    end();

    // Velocity correction
    begin("Computing elastic force");
    interpolation.fluid_to_solid(*u1, *body_velocity);
    auto temp_disp = std::make_shared<Function>(U);
    *temp_disp = FunctionAXPY(body_velocity, dt);
    ALE::move(*circle, *temp_disp);
    *temp_disp = FunctionAXPY(body_velocity, dt)+body_disp;
    *body_disp = *temp_disp;
    L4.u = body_disp;
    assemble(b4, L4);
    solve(A4, *body_force->vector(), b4, "cg", "sor");
    interpolation.solid_to_fluid(*f, *body_force);
    L1.f = f;
    end();

    // Save to file
    ufile << *u1;
    pfile << *p1;
    ffile << *f;
    bfile << *body_force;

    // Move to next time step
    *u0 = *u1;
    t += dt;
    cout << "t = " << t << endl;
  }

  return 0;
}
