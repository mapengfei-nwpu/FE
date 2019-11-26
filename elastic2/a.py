from fenics import *

# Scaled variables
L = 1; W = 0.2
mu = 1
rho = 1
delta = W/L
gamma = 0.4*delta**2
beta = 1.25
lambda_ = beta
g = gamma

# Create mesh and define function space
mesh = RectangleMesh(Point(0, 0), Point(L, W), 10, 3)
V = VectorFunctionSpace(mesh, 'P', 1)

# Define boundary condition
tol = 1E-14

def clamped_boundary(x, on_boundary):
    return on_boundary and x[0] < tol

bc = DirichletBC(V, Constant((0, 0)), clamped_boundary)

# Define strain and stress

def epsilon(u):
    return 0.5*(nabla_grad(u) + nabla_grad(u).T)
    #return sym(nabla_grad(u))

def sigma(u):
    return lambda_*div(u)*Identity(d) + 2*mu*epsilon(u)

# Define variational problem
u = TrialFunction(V)
d = u.geometric_dimension()  # space dimension
v = TestFunction(V)
f = Constant((0,  -rho*g))
T = Constant((0,  0))
a = inner(sigma(u), epsilon(v))*dx
L = dot(f, v)*dx + dot(T, v)*ds

# Compute solution
u = Function(V)
solve(a == L, u, bc)
move(mesh,u)









mu = 1
beta = 1.25
lambda_ = beta

# define function element
V = VectorElement("P", triangle, 1)

# Define strain and stress
def epsilon(u):
    return 0.5*(nabla_grad(u) + nabla_grad(u).T)

def sigma(u):
    return lambda_*div(u)*Identity(2) + 2*mu*epsilon(u)

# Define variational problem
u = TrialFunction(V)
v = TestFunction(V)

# body force and tranction force
f = Coefficient(V)

Q = Element("P", triangle, 1)
p_f = Coefficient(Q)
u_f = Coefficient(V)
n = FacetNormal(triangle)
mu_f = 0.001
T = - p_f*n + 2*mu_f*epsilon(u_f)*n

a = inner(sigma(u), epsilon(v))*dx
L = dot(f, v)*dx + dot(T, v)*ds