import fenics as fnx
import mshr
import numpy as np
import sympy as sym
import itertools
from matplotlib import pyplot as plt

p = fnx.Expression("A / (2*pi*sigma) * exp(-1/2 * (pow(x[0] - x0, 2) + pow(x[1] - y0, 2)))",
                   degree=1, A=1, sigma=1, x0=1, y0=1)

domain = mshr.Circle(fnx.Point(0, 0), 1)
mesh = mshr.generate_mesh(domain, 64)

beta = 8
R0 = 0.6
p = fnx.Expression("4*exp(-pow(beta, 2)*(pow(x[0], 2) + pow(x[1] - R0, 2)))", degree=1, beta=beta, R0=R0)

V = fnx.FunctionSpace(mesh, 'Lagrange', 1)
u_D = fnx.Constant(0)
def boundary(x, on_mesh_boundary):
    return on_mesh_boundary
bc = fnx.DirichletBC(V, u_D, boundary)

w = fnx.TrialFunction(V)
v = fnx.TestFunction(V)
a = fnx.dot(fnx.grad(w), fnx.grad(v))*fnx.dx
L = p*v*fnx.dx

w = fnx.Function(V)
fnx.solve(a == L, w, bc)

vtkfile = fnx.File("deflection.pvd")
vtkfile << w
vtkfile_p = fnx.File("deflection_pressure.pvd")
# p = fnx.interpolate(p, V)
p = fnx.project(p, V)
vtkfile_p << p

# Create a curve plot along x=0.
tol = 1e-3
nodes = np.linspace(-1 + tol, 1 - tol, 300)
points = [(0, y) for y in nodes]
plt.plot(nodes, 50*np.array(list(map(w, points))), "b--")
plt.plot(nodes, np.array(list(map(p, points))), "k")
plt.xlabel("$y$")
plt.legend(["Deflection x 50", "Load"], loc="upper left")
plt.savefig("deflection_curve.pdf")
