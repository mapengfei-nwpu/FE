from fenics import *
from mshr import *

domain = Circle(Point(0.6, 0.5), 0.2)
mesh = generate_mesh(domain, 6)

File("circle.xml.gz") << mesh

print(mesh.hmax())
