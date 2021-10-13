import sympy as sym

# Helper function, convert barycentric index to an index into
# a flat array.
def barycentric_index_to_linear(i,j,k):
    assert(i+j+k == 2)
    if (i==2):
        return 0
    if (j==2):
        return 1
    if (k==2):
        return 2
    if (i==1 and j==1):
        return 3
    if (j==1 and k==1):
        return 4
    if (k==1 and i==1):
        return 5

# Define algebraic symbols for the variables.
x,y,z = sym.symbols("x y z")

# The reference triangle has vertices v,vp,vpp.
# We define K1,K2,K3 by
#     K1 = (v - vpp)^perp,
#     K2 = (vp - v)^perp,
#     K3 = (vpp - vp)^perp.
# K1 and K2 are the columns of the matrix K (which is the inverse transpose Jacobian without the constant term).
# while K3 = -K1-K2.
K1,K2,K3 = sym.symbols("K1 K2 K3")

# Define the quadratic basis functions given in barycentric coordinates (x,y,z).
# The order with respect to barycentric indices is 200,020,002, 110,011,101.
nodal_basis_functions_barycentric = [
    x*(1 - 2*y - 2*z),
    y*(1 - 2*z - 2*x),
    z*(1 - 2*x - 2*y),
    4*x*y,
    4*y*z,
    4*z*x
]
# Map the nodal basis functions to the reference triangle by substituting z with 1-x-y.
nodal_basis_functions = [p.subs(z,1-x-y).expand() for p in nodal_basis_functions_barycentric]
# Compute the gradients of the nodal basis functions on the reference triangle.
nodal_basis_gradients = [(sym.diff(p, x), sym.diff(p, y)) for p in nodal_basis_functions]


# Compute the integral
def integrate_gradients(a,b,c, i,j,k):
    grad1 = nodal_basis_gradients[barycentric_index_to_linear(a,b,c)]
    grad2 = nodal_basis_gradients[barycentric_index_to_linear(i,j,k)]
    f = (K1*grad1[0] + K2*grad1[1])*(K1*grad2[0] + K2*grad2[1])

    f_dy = sym.integrate(f, (y, 0,1-x))
    f_dy_dx = sym.integrate(f_dy, (x, 0,1))

    print("{}{}{},{}{}{}:".format(a,b,c,i,j,k), f_dy_dx.simplify().subs(-K1-K2,K3))


def grads():
    print("vertex integrals")
    integrate_gradients(0,0,2,  0,0,2)
    integrate_gradients(0,0,2,  2,0,0)
    integrate_gradients(0,0,2,  0,2,0)
    
    integrate_gradients(0,0,2,  1,0,1)
    integrate_gradients(0,0,2,  0,1,1)
    integrate_gradients(0,0,2,  1,1,0)
    
    print("midpoint integrals")
    integrate_gradients(1,1,0,  1,1,0)
    integrate_gradients(1,1,0,  0,1,1)
    integrate_gradients(1,1,0,  1,0,1)
    
    integrate_gradients(1,1,0,  2,0,0)
    integrate_gradients(1,1,0,  0,2,0)
    integrate_gradients(1,1,0,  0,0,2)


# Source term integration
def integrate_scalars(a,b,c, i,j,k):
    f1 = nodal_basis_functions[barycentric_index_to_linear(a,b,c)]
    f2 = nodal_basis_functions[barycentric_index_to_linear(i,j,k)]
    f = f1*f2
    f_dy = sym.integrate(f, (y, 0,1-x))
    f_dy_dx = sym.integrate(f_dy, (x, 0,1))
    print("{}{}{},{}{}{}:".format(a,b,c,i,j,k), f_dy_dx.simplify().subs(-K1-K2,K3))


def scalars():
    print("source term vertex integrals")
    integrate_scalars(0,0,2,  0,0,2)
    integrate_scalars(0,0,2,  2,0,0)
    integrate_scalars(0,0,2,  0,2,0)
    
    integrate_scalars(0,0,2,  1,0,1)
    integrate_scalars(0,0,2,  0,1,1)
    integrate_scalars(0,0,2,  1,1,0)

    print("source term midpoint integrals")
    integrate_scalars(1,1,0,  0,0,2)
    integrate_scalars(1,1,0,  2,0,0)
    integrate_scalars(1,1,0,  0,2,0)
    
    integrate_scalars(1,1,0,  1,0,1)
    integrate_scalars(1,1,0,  0,1,1)
    integrate_scalars(1,1,0,  1,1,0)


# grads()
# scalars()

# Computation for deriving values for Stokes flow assembly.
pressure_basis_functions = [
    1-x-y,
    x,
    y
] # v, v', v''

def integrate_scalar(psi_p_index, phi_u_i, phi_u_j, phi_u_k):
    for velocity_component in [x,y]:
        psi_p = pressure_basis_functions[psi_p_index]
        phi_u = nodal_basis_functions[barycentric_index_to_linear(phi_u_i, phi_u_j, phi_u_k)]
        diff_phi_u = sym.diff(phi_u, velocity_component)
        f = psi_p * diff_phi_u
        f_dy = sym.integrate(f, (y, 0,1-x))
        f_dy_dx = sym.integrate(f_dy, (x, 0,1))
        print("{}, v{}, {}{}{}: {}".format(velocity_component, '\''*psi_p_index, phi_u_i,phi_u_j,phi_u_k, f_dy_dx))

print("vertex")
integrate_scalar(0, 0,0,2)
integrate_scalar(0, 2,0,0)
integrate_scalar(0, 0,2,0)
print("midpoint")
integrate_scalar(0, 1,1,0)
integrate_scalar(0, 0,1,1)
integrate_scalar(0, 1,0,1)
