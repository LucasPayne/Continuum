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

def linear_index_to_barycentric(index):
    assert(index >= 0 and index <= 5)
    if index == 0:
        return (2,0,0)
    if index == 1:
        return (0,2,0)
    if index == 2:
        return (0,0,2)
    if index == 3:
        return (1,1,0)
    if index == 4:
        return (0,1,1)
    if index == 5:
        return (1,0,1)


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

pressure_basis_functions = [
    1-x-y,
    x,
    y
] # v, v', v''

# # Compute the integral
# def integrate_gradients(a,b,c, i,j,k):
#     grad1 = nodal_basis_gradients[barycentric_index_to_linear(a,b,c)]
#     grad2 = nodal_basis_gradients[barycentric_index_to_linear(i,j,k)]
#     f = (K1*grad1[0] + K2*grad1[1])*(K1*grad2[0] + K2*grad2[1])
# 
#     f_dy = sym.integrate(f, (y, 0,1-x))
#     f_dy_dx = sym.integrate(f_dy, (x, 0,1))
# 
#     print("{}{}{},{}{}{}:".format(a,b,c,i,j,k), f_dy_dx.simplify().subs(-K1-K2,K3))
# 
# 
# def grads():
#     print("vertex integrals")
#     integrate_gradients(0,0,2,  0,0,2)
#     integrate_gradients(0,0,2,  2,0,0)
#     integrate_gradients(0,0,2,  0,2,0)
#     
#     integrate_gradients(0,0,2,  1,0,1)
#     integrate_gradients(0,0,2,  0,1,1)
#     integrate_gradients(0,0,2,  1,1,0)
#     
#     print("midpoint integrals")
#     integrate_gradients(1,1,0,  1,1,0)
#     integrate_gradients(1,1,0,  0,1,1)
#     integrate_gradients(1,1,0,  1,0,1)
#     
#     integrate_gradients(1,1,0,  2,0,0)
#     integrate_gradients(1,1,0,  0,2,0)
#     integrate_gradients(1,1,0,  0,0,2)
# 
# 
# # Source term integration
# def integrate_scalars(a,b,c, i,j,k):
#     f1 = nodal_basis_functions[barycentric_index_to_linear(a,b,c)]
#     f2 = nodal_basis_functions[barycentric_index_to_linear(i,j,k)]
#     f = f1*f2
#     f_dy = sym.integrate(f, (y, 0,1-x))
#     f_dy_dx = sym.integrate(f_dy, (x, 0,1))
#     print("{}{}{},{}{}{}:".format(a,b,c,i,j,k), f_dy_dx.simplify().subs(-K1-K2,K3))
# 
# 
# def scalars():
#     print("source term vertex integrals")
#     integrate_scalars(0,0,2,  0,0,2)
#     integrate_scalars(0,0,2,  2,0,0)
#     integrate_scalars(0,0,2,  0,2,0)
#     
#     integrate_scalars(0,0,2,  1,0,1)
#     integrate_scalars(0,0,2,  0,1,1)
#     integrate_scalars(0,0,2,  1,1,0)
# 
#     print("source term midpoint integrals")
#     integrate_scalars(1,1,0,  0,0,2)
#     integrate_scalars(1,1,0,  2,0,0)
#     integrate_scalars(1,1,0,  0,2,0)
#     
#     integrate_scalars(1,1,0,  1,0,1)
#     integrate_scalars(1,1,0,  0,1,1)
#     integrate_scalars(1,1,0,  1,1,0)
# 
# 
# # grads()
# # scalars()
# 
# # Computation for deriving values for Stokes flow assembly.
# 
# # def integrate_scalar(psi_p_index, phi_u_i, phi_u_j, phi_u_k):
# #     for velocity_component in [x,y]:
# #         psi_p = pressure_basis_functions[psi_p_index]
# #         phi_u = nodal_basis_functions[barycentric_index_to_linear(phi_u_i, phi_u_j, phi_u_k)]
# #         diff_phi_u = sym.diff(phi_u, velocity_component)
# #         f = psi_p * diff_phi_u
# #         f_dy = sym.integrate(f, (y, 0,1-x))
# #         f_dy_dx = sym.integrate(f_dy, (x, 0,1))
# #         print("{}, v{}, {}{}{}: {}".format(velocity_component, '\''*psi_p_index, phi_u_i,phi_u_j,phi_u_k, f_dy_dx))
# # 
# # print("vertex")
# # integrate_scalar(0, 0,0,2)
# # integrate_scalar(0, 2,0,0)
# # integrate_scalar(0, 0,2,0)
# # print("midpoint")
# # integrate_scalar(0, 1,1,0)
# # integrate_scalar(0, 0,1,1)
# # integrate_scalar(0, 1,0,1)
# 
# 
# 
# # def integrate_scalar(a,b,c,  phi_p_index):
# #     for velocity_component in [x,y]:
# #         psi_u = nodal_basis_functions[barycentric_index_to_linear(a, b, c)]
# #         phi_p = pressure_basis_functions[phi_p_index]
# #         diff_phi_p = sym.diff(phi_p, velocity_component)
# #         f = diff_phi_p * psi_u
# #         f_dy = sym.integrate(f, (y, 0,1-x))
# #         f_dy_dx = sym.integrate(f_dy, (x, 0,1))
# #         print("{}, {}{}{}, {}: {}".format(velocity_component, a,b,c, phi_p_index, f_dy_dx))
# # 
# # print("vertex")
# # integrate_scalar(0,0,2, 0)
# # integrate_scalar(0,0,2, 1)
# # integrate_scalar(0,0,2, 2)
# # print("midpoint")
# # integrate_scalar(1,1,0, 0)
# # integrate_scalar(1,1,0, 1)
# # integrate_scalar(1,1,0, 2)
# 
# def integrate_scalar(a,b,c):
#     f = nodal_basis_functions[barycentric_index_to_linear(a, b, c)]
#     f_dy = sym.integrate(f, (y, 0,1-x))
#     f_dy_dx = sym.integrate(f_dy, (x, 0,1))
#     print("{}{}{}, {}".format(a,b,c, f_dy_dx))
# 
# print("vertex")
# integrate_scalar(0,0,2)
# integrate_scalar(2,0,0)
# integrate_scalar(0,2,0)
# print("midpoint")
# integrate_scalar(1,1,0)
# integrate_scalar(0,1,1)
# integrate_scalar(1,0,1)
# 
# 
# print("gramian")
# def gramian_integrate(i,j):
#     phi_i = pressure_basis_functions[i]
#     psi_j = pressure_basis_functions[j]
#     f = phi_i * psi_j
#     f_dy = sym.integrate(f, (y, 0,1-x))
#     f_dy_dx = sym.integrate(f_dy, (x, 0,1))
#     print("{},{}: {}".format(i,j,f_dy_dx))
# gramian_integrate(0,0)
# gramian_integrate(0,1)
# gramian_integrate(0,2)
# 
# 
# 
# print("grads")
# def grads(a,b,c):
#     phi_u = nodal_basis_functions[barycentric_index_to_linear(a, b, c)]
#     # print("{}{}{}: (X: {}, Y: {})".format(a,b,c, sym.diff(phi_u, x), sym.diff(phi_u, y)))
#     psi_p = pressure_basis_functions[0]
#     f = (K1*sym.diff(phi_u, x) + K2*sym.diff(phi_u, y)) * psi_p
#     f_dy = sym.integrate(f, (y, 0,1-x))
#     f_dy_dx = sym.integrate(f_dy, (x, 0,1))
#     print("{}{}{}: {}".format(a,b,c, f_dy_dx))
# grads(0,0,2)
# grads(0,2,0)
# grads(2,0,0)
# grads(0,1,1)
# grads(1,1,0)
# grads(1,0,1)
#     
# 
# print("vertex grads u to u")
# def grads(i,j,k, a,b,c):
#     psi_u = nodal_basis_functions[barycentric_index_to_linear(a,b,c)]
#     phi_u = nodal_basis_functions[barycentric_index_to_linear(i,j,k)]
#     f = (K1*sym.diff(phi_u, x) + K2*sym.diff(phi_u, y)) * psi_u
#     f_dy = sym.integrate(f, (y, 0,1-x))
#     f_dy_dx = sym.integrate(f_dy, (x, 0,1))
#     print("{}{}{}, {}{}{}: {}".format(i,j,k, a,b,c, f_dy_dx))
# grads(0,0,2,  0,0,2)
# grads(0,0,2,  0,2,0)
# grads(0,0,2,  2,0,0)
# grads(0,0,2,  0,1,1)
# grads(0,0,2,  1,1,0)
# grads(0,0,2,  1,0,1)
#     
# 
# print("midpoint grads u to u")
# grads(1,1,0,  0,0,2)
# grads(1,1,0,  0,2,0)
# grads(1,1,0,  2,0,0)
# grads(1,1,0,  0,1,1)
# grads(1,1,0,  1,1,0)
# grads(1,1,0,  1,0,1)
# 
# 
# def scalar_velocity_gramian(i,j,k, a,b,c):
#     phi_u = nodal_basis_functions[barycentric_index_to_linear(i,j,k)]
#     psi_u = nodal_basis_functions[barycentric_index_to_linear(a,b,c)]
#     f = phi_u * psi_u
#     f_dy = sym.integrate(f, (y, 0,1-x))
#     f_dy_dx = sym.integrate(f_dy, (x, 0,1))
#     print("{}{}{}, {}{}{}: {}".format(i,j,k, a,b,c, f_dy_dx))
# print("scalar velocity gramian")
# print("vertex")
# scalar_velocity_gramian(0,0,2,  0,0,2)
# scalar_velocity_gramian(0,0,2,  0,2,0)
# scalar_velocity_gramian(0,0,2,  2,0,0)
# scalar_velocity_gramian(0,0,2,  0,1,1)
# scalar_velocity_gramian(0,0,2,  1,1,0)
# scalar_velocity_gramian(0,0,2,  1,0,1)
# print("midpoint")
# scalar_velocity_gramian(1,1,0,  0,0,2)
# scalar_velocity_gramian(1,1,0,  0,2,0)
# scalar_velocity_gramian(1,1,0,  2,0,0)
# scalar_velocity_gramian(1,1,0,  0,1,1)
# scalar_velocity_gramian(1,1,0,  1,1,0)
# scalar_velocity_gramian(1,1,0,  1,0,1)
# 
# 
# def vel(i,j,k, a,b,c):
#     psi_u = nodal_basis_functions[barycentric_index_to_linear(i,j,k)]
#     phi_u = nodal_basis_functions[barycentric_index_to_linear(a,b,c)]
#     f = -phi_u * (K1*sym.diff(psi_u, x) + K2*sym.diff(psi_u, y))
#     f_dy = sym.integrate(f, (y, 0,1-x))
#     f_dy_dx = sym.integrate(f_dy, (x, 0,1))
#     print("{}{}{}, {}{}{}: {}".format(i,j,k, a,b,c, f_dy_dx))
# print("------------------------------------------------------------")
# print("vertex")
# vel(0,0,2,  0,0,2)
# vel(0,0,2,  0,2,0)
# vel(0,0,2,  2,0,0)
# vel(0,0,2,  0,1,1)
# vel(0,0,2,  1,1,0)
# vel(0,0,2,  1,0,1)
# print("midpoint")
# vel(1,1,0,  0,0,2)
# vel(1,1,0,  0,2,0)
# vel(1,1,0,  2,0,0)
# vel(1,1,0,  0,1,1)
# vel(1,1,0,  1,1,0)
# vel(1,1,0,  1,0,1)
# 
# 
# def div(i,j,k, a,b,c):
#     psi_u = nodal_basis_functions[barycentric_index_to_linear(i,j,k)]
#     phi_u = nodal_basis_functions[barycentric_index_to_linear(a,b,c)]
#     f = (K1*sym.diff(phi_u, x) + K2*sym.diff(phi_u, y)) * psi_u
#     # f = -(K1*sym.diff(psi_u, x) + K2*sym.diff(psi_u, y)) * phi_u
#     f_dy = sym.integrate(f, (y, 0,1-x))
#     f_dy_dx = sym.integrate(f_dy, (x, 0,1))
#     print("{}{}{}, {}{}{}: {}".format(i,j,k, a,b,c, f_dy_dx))
# print("------------------------------------------------------------")
# print("vertex")
# div(0,0,2,  0,0,2)
# div(0,0,2,  0,2,0)
# div(0,0,2,  2,0,0)
# div(0,0,2,  0,1,1)
# div(0,0,2,  1,1,0)
# div(0,0,2,  1,0,1)
# print("midpoint")
# div(1,1,0,  0,0,2)
# div(1,1,0,  0,2,0)
# div(1,1,0,  2,0,0)
# div(1,1,0,  0,1,1)
# div(1,1,0,  1,1,0)
# div(1,1,0,  1,0,1)
# 
# 
# # def div2(i,j,k, a,b,c):
# #     psi_u = nodal_basis_functions[barycentric_index_to_linear(i,j,k)]
# #     phi_u = nodal_basis_functions[barycentric_index_to_linear(a,b,c)]
# #     f = -phi_u * (K1*sym.diff(psi_u, x) + K2*sym.diff(psi_u, y))
# #     f_dy = sym.integrate(f, (y, 0,1-x))
# #     f_dy_dx = sym.integrate(f_dy, (x, 0,1))
# #     print("{}{}{}, {}{}{}: {}".format(i,j,k, a,b,c, f_dy_dx))
# # print("div2")
# # print("------------------------------------------------------------")
# # print("vertex")
# # div2(0,0,2,  0,0,2)
# # div2(0,0,2,  0,2,0)
# # div2(0,0,2,  2,0,0)
# # div2(0,0,2,  0,1,1)
# # div2(0,0,2,  1,1,0)
# # div2(0,0,2,  1,0,1)
# # print("midpoint")
# # div2(1,1,0,  0,0,2)
# # div2(1,1,0,  0,2,0)
# # div2(1,1,0,  2,0,0)
# # div2(1,1,0,  0,1,1)
# # div2(1,1,0,  1,1,0)
# # div2(1,1,0,  1,0,1)

def gramian_integrate(a,b,c, n):
    psi_u = nodal_basis_functions[barycentric_index_to_linear(a,b,c)]
    phi_u = nodal_basis_functions[n]
    f = phi_u * psi_u
    f_dy = sym.integrate(f, (y, 0,1-x))
    f_dy_dx = sym.integrate(f_dy, (x, 0,1))
    i,j,k = linear_index_to_barycentric(n)
    print("{},{},{} {},{},{}: {}".format(a,b,c, i,j,k,  f_dy_dx))
print("vertex")
for n in range(6):
    gramian_integrate(0,0,2, n)
print("midpoint")
for n in range(6):
    gramian_integrate(1,1,0, n)


print("vertex grads u to u")
def grads(i,j,k, a,b,c):
    psi_u = nodal_basis_functions[barycentric_index_to_linear(a,b,c)]
    phi_u = nodal_basis_functions[barycentric_index_to_linear(i,j,k)]
    f = (K1*sym.diff(phi_u, x) + K2*sym.diff(phi_u, y)) * (K1*sym.diff(psi_u, x) + K2*sym.diff(psi_u, y))
    f_dy = sym.integrate(f, (y, 0,1-x))
    f_dy_dx = sym.integrate(f_dy, (x, 0,1))
    print("{}{}{}, {}{}{}: {}".format(i,j,k, a,b,c, f_dy_dx))
grads(0,0,2,  0,0,2)
grads(0,0,2,  0,2,0)
grads(0,0,2,  2,0,0)
grads(0,0,2,  0,1,1)
grads(0,0,2,  1,1,0)
grads(0,0,2,  1,0,1)
    

print("midpoint grads u to u")
grads(1,1,0,  0,0,2)
grads(1,1,0,  0,2,0)
grads(1,1,0,  2,0,0)
grads(1,1,0,  0,1,1)
grads(1,1,0,  1,1,0)
grads(1,1,0,  1,0,1)



def integrate(a,b,c, n):
    psi_u = nodal_basis_functions[barycentric_index_to_linear(a,b,c)]
    phi_p = pressure_basis_functions[n]
    f = -phi_p * (K1*sym.diff(psi_u, x) + K2*sym.diff(psi_u, y))
    f_dy = sym.integrate(f, (y, 0,1-x))
    f_dy_dx = sym.integrate(f_dy, (x, 0,1))
    print("{},{},{} {}: {}".format(a,b,c, n,  f_dy_dx))
print("vertex")
for n in range(3):
    integrate(0,0,2, n)
print("midpoint")
for n in range(3):
    integrate(1,1,0, n)


def integrate(n, m):
    psi_p = pressure_basis_functions[n]
    phi_u = nodal_basis_functions[m]
    f = -psi_p * (K1*sym.diff(phi_u, x) + K2*sym.diff(phi_u, y))
    f_dy = sym.integrate(f, (y, 0,1-x))
    f_dy_dx = sym.integrate(f_dy, (x, 0,1))
    i,j,k = linear_index_to_barycentric(m)
    print("{} {},{},{}: {}".format(n, i,j,k, f_dy_dx))
print("vertex")
for m in range(6):
    integrate(0, m)
