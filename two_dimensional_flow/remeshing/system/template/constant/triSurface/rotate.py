import sys
from stl import mesh

angle = float(sys.argv[1])

mesh = mesh.Mesh.from_file('inner.stl')
mesh.rotate([0, 0, 1], angle)
mesh.save('inner.stl')

