import sympy as sym

x,y = sym.symbols("x y")

f = 0.5*sym.cos(8*x)*y + sym.exp(-x*x)*0.5*y*y+1.3;

lap = -(sym.diff(f, x,x) + sym.diff(f, y,y))
print(lap)

