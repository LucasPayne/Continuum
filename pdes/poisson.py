import fenics as fnx
import numpy as np
import sympy as sym
import itertools

def poisson_problem(tess):
    # Create mesh and define function space.
    mesh = fnx.UnitSquareMesh(tess, tess)
    V = fnx.FunctionSpace(mesh, 'Lagrange', 1)
    # Lagrange/P: Lagrange family of elements.
    # DP: Discontinuous Galerkin elements.
    
    # Define boundary condition.
    #     Formula is in C++ syntax.
    #     u_D is the true solution, so it should be sampled on the mesh at a higher degree for comparison.
    #     Since this is a quadratic, setting degree=2 gives exact sampling.
    # Expression objects have built-in access to
    # standard math.h trig functions, exp, log, pow, sqrt, ceil, floor, fabs, etc.
    # Instead of 'M_PI', 'pi' is used.
    u_D = fnx.Expression("1 + x[0]*x[0] + 2*x[1]*x[1]", degree=2)
    # on_mesh_boundary is w/r/t the physical boundary of the mesh. Since the physical boundary is the boundary used for the PDE,
    # this function just passes on_mesh_boundary through.
    # boundary() is called for each discrete point in the mesh, so boundary conditions can also be in the interior of the mesh.
    def boundary(x, on_mesh_boundary):
        return on_mesh_boundary
    # Alternatively
    # tol = 1e-14
    # def boundary(x):
    #     return any(fnx.near(x[i], j, tol) for i,j in itertools.product(range(2), repeat=2))
    
    bc = fnx.DirichletBC(V, u_D, boundary)
    
    # Define the variational problem.
    u = fnx.TrialFunction(V)
    v = fnx.TestFunction(V)
    f = fnx.Constant(-6)
    # Equivalent to
    # f = fnx.Expression("-6", degree=0)
    a = fnx.dot(fnx.grad(u), fnx.grad(v))*fnx.dx
    L = f*v*fnx.dx
    
    # Compute solution.
    u = fnx.Function(V)
    fnx.solve(a == L, u, bc)
    
    # Plot solution and mesh.
    fnx.plot(u)
    fnx.plot(mesh)
    
    error_L2 = fnx.errornorm(u_D, u, "L2")
    vertex_values_u_D = u_D.compute_vertex_values(mesh)
    vertex_values_u = u.compute_vertex_values(mesh)
    error_max = np.max(np.abs(vertex_values_u_D - vertex_values_u))
    
    print("tess =", tess)
    print("    error_L2 =", error_L2)
    print("    error_max =", error_max)
    return u

# for tess in range(5, 20+1):
#     poisson_problem(tess)

u = poisson_problem(10)
# Export for visualization in ParaView.
vtkfile = fnx.File("poisson_solution.pvd")
vtkfile << u
