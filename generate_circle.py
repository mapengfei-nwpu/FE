from fenics import *
from mshr import *

domain = Circle(Point(0.2, 0.2), 0.05)
mesh = generate_mesh(domain, 64)

File("circle.xml.gz") << mesh
