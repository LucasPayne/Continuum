import fenics as fnx
import sympy as sym
import numpy as np
import scipy as sci
import scipy.integrate
import scipy.linalg
import itertools
from matplotlib import pyplot as plt
from printing import print_matrix
from math import sqrt
    

def least_squares(basis_functions, domain_limits, f, symbolic=False):
    """
        basis_functions: List of functions of x.
        domain_limits: Tuple (a,b) giving the domain of approximation.
        f: Function of x to approximate with basis functions.
    """
    n = len(basis_functions)
    A = sym.zeros(n, n)
    x = sym.symbols('x')
    for i,j in itertools.combinations_with_replacement(range(n), 2):
        g = basis_functions[i] * basis_functions[j]
        if symbolic:
            A[i,j] = sym.integrate(g, (x, *domain_limits))
        else:
            integrand = sym.lambdify([x], g)
            val, err = sci.integrate.quadrature(integrand, *domain_limits)
            A[i, j] = val
        A[j,i] = A[i, j]
    print_matrix(A)

    if symbolic:
        b = sym.Matrix([sym.integrate(f*g, (x, *domain_limits)) for g in basis_functions])
    else:
        b = sym.Matrix([sci.integrate.quadrature(sym.lambdify([x], f*g), *domain_limits)[0] for g in basis_functions])
    print_matrix(b)

    if symbolic:
        c = A.LUsolve(b)
    else:
        c = np.linalg.solve(np.array(A).astype(np.float32), np.array(b).astype(np.float32))
    print_matrix(sym.Matrix(c))

    return sum(c[i,0]*basis_functions[i] for i in range(n))

x = sym.symbols('x')
f = sym.sin(sym.pi * x) + sym.cos(sym.pi * 4 * x)
domain_limits = (0, 10)
fp = least_squares([1] + [trig(x*i*2*sym.pi/10) for i in range(1, 5+1) for trig in [sym.sin, sym.cos]], (0,10), f, False)
nodes = np.linspace(*domain_limits, 1000)
ffunc, fpfunc = (sym.lambdify([x], func, modules="numpy") for func in [f, fp])

error_L2 = sqrt(sci.integrate.quadrature(lambda x: (ffunc(x) - fpfunc(x))**2, *domain_limits)[0])
print("error_L2:", error_L2)

plt.plot(nodes, ffunc(nodes), "k")
plt.plot(nodes, fpfunc(nodes), "b")
plt.show()
