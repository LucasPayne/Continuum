from matplotlib import pyplot as plt
from math import pi,sin,cos,asin,acos,exp
import numpy as np
import itertools
import sympy as sym

def vf(x, y):
    # return (x**2-y, x**2+y)
    return (np.sin(x)**2 + np.sin(y)**2, -x**2 + np.exp(y-x**2))
def vf_div(x, y):
    # return 2*x + 1
    return 2*np.cos(x)*np.sin(x) + np.exp(y-x**2)

def plot_vf(f,f_div, xa,xb, ya,yb):
    X,Y = np.meshgrid(np.linspace(xa, xb, 100),
                      np.linspace(ya, yb, 100))
    P,Q = f(X,Y)
    plt.contour(X,Y,f_div(X,Y))
    plt.streamplot(X,Y, P,Q)
            

plot_vf(vf, vf_div, -1,1, -1,1)
plt.show()

# x,y,z = sym.symbols("x y z", real=True)
# F1 = sym.Function("F_1")(x,y,z)
# F2 = sym.Function("F_2")(x,y,z)
# F3 = sym.Function("F_3")(x,y,z)
# F = sym.Matrix([F1, F2, F3])
# J = F.jacobian([x,y,z])
# N = J - J.transpose()
# a1,a2,a3 = sym.symbols("a_1 a_2 a_3")
# b1,b2,b3 = sym.symbols("b_1 b_2 b_3")
# a = sym.Matrix([a1,a2,a3])
# b = sym.Matrix([b1,b2,b3])
# val = (a.transpose() * N * b).expand()[0]
# lprint("", val)
# ldone()

