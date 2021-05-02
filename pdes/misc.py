from fenics import *


# Parameters are initialized by keyword arguments.
f = Expression("exp(-kappa*pow(pi, 2)*t)*sin(pi*k*x[0])", degree=2,
               kappa=1, t=0, k=4)
print(f)
# Supposedly the generated code allows parameters instead of hardcoding
# constants.
f.t += 0.1
print(f)
