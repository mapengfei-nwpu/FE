#include <cmath>
#include <sstream>
#include "constants.h"
#include "Point.h"


/**********************warning********************/
//additional functions start
/**********************warning********************/
bool dolfin::near(double x, double x0, double eps){
  return x0 - eps <= x && x <= x0 + eps;
}
#include<stdexcept>
#include<iostream>
void dolfin::dolfin_assert(bool a){
    if (!a)throw std::runtime_error("It is wrong!haha~");
    else std::cout<<"OK"<<std::endl;
}
/**********************warning********************/
//additional functions end
/**********************warning********************/


using namespace dolfin;
const Point Point::cross(const Point& p) const
{
  Point q;

  q._x[0] = _x[1]*p._x[2] - _x[2]*p._x[1];
  q._x[1] = _x[2]*p._x[0] - _x[0]*p._x[2];
  q._x[2] = _x[0]*p._x[1] - _x[1]*p._x[0];

  return q;
}
//-----------------------------------------------------------------------------
double Point::squared_distance(const Point& p) const
{
  const double dx = p._x[0] - _x[0];
  const double dy = p._x[1] - _x[1];
  const double dz = p._x[2] - _x[2];

  return dx*dx + dy*dy + dz*dz;
}
//-----------------------------------------------------------------------------
double Point::dot(const Point& p) const
{
  return _x[0]*p._x[0] + _x[1]*p._x[1] + _x[2]*p._x[2];
}
//-----------------------------------------------------------------------------
Point Point::rotate(const Point& k, double theta) const
{
  dolfin_assert(near(k.norm(), 1.0));

  const Point& v = *this;
  const double cosTheta = cos(theta);
  const double sinTheta = sin(theta);

  //Rodriques' rotation formula
  return v*cosTheta + k.cross(v)*sinTheta + k*k.dot(v)*(1-cosTheta);
}
//-----------------------------------------------------------------------------
std::string Point::str(bool verbose) const
{
  std::stringstream s;
  s << "<Point x = " << x() << " y = " << y() << " z = " << z() << ">";
  return s.str();
}
//-----------------------------------------------------------------------------
