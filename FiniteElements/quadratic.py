import sympy as sym

def basis_index(i,j,k):
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


x,y,z = sym.symbols("x y z")
K1,K2,K3 = sym.symbols("K1 K2 K3")
ps = [
    x*(1 - 2*y - 2*z),
    y*(1 - 2*z - 2*x),
    z*(1 - 2*x - 2*y),
    4*x*y,
    4*y*z,
    4*z*x
]
ps_ref = [p.subs(z,1-x-y).expand() for p in ps]

ps_ref_gradients = []
for p in ps_ref:
    ddx = sym.diff(p, x)
    ddy = sym.diff(p, y)
    ps_ref_gradients.append(sym.Matrix([ddx, ddy]))
for g in ps_ref_gradients:
    print(g)

def comp_int(f):
    fi = sym.integrate(f, (y, 0,1-x))
    return sym.integrate(fi, (x, 0,1))

def grad_int(a,b,c, i,j,k):
    g1 = ps_ref_gradients[basis_index(a,b,c)]
    g2 = ps_ref_gradients[basis_index(i,j,k)]
    f = (K1*g1[0] + K2*g1[1])*(K1*g2[0] + K2*g2[1])
    inte = comp_int(f)
    print(inte.simplify().subs(-K1-K2,K3))

grad_int(0,0,2,  0,0,2)
grad_int(0,0,2,  2,0,0)
grad_int(0,0,2,  0,2,0)

grad_int(0,0,2,  1,0,1)
grad_int(0,0,2,  0,1,1)
grad_int(0,0,2,  1,1,0)

print("midpoint integrals")
grad_int(1,1,0,  1,1,0)
grad_int(1,1,0,  0,1,1)
grad_int(1,1,0,  1,0,1)

grad_int(1,1,0,  2,0,0)
grad_int(1,1,0,  0,2,0)
grad_int(1,1,0,  0,0,2)

# print(comp_int(ps_ref[basis_index(1,1,0)]**2))
# print(comp_int(ps_ref[basis_index(1,1,0)] * ps_ref[basis_index(2,0,0)]))
# print(comp_int(ps_ref[basis_index(1,1,0)] * ps_ref[basis_index(0,2,0)]))
# print(comp_int(ps_ref[basis_index(1,1,0)] * ps_ref[basis_index(0,1,1)]))
# print(comp_int(ps_ref[basis_index(1,1,0)] * ps_ref[basis_index(1,0,1)]))
