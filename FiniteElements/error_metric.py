import numpy as np
from matplotlib import pyplot as plt
import sys

assert(len(sys.argv) == 2)

filename = sys.argv[1]

f = open(filename, "r")
lines = f.readlines()
hs = []
errors = []
for line in lines:
    vals = [float(val) for val in line.strip().split()]
    mesh_N = vals[0]
    hs.append(vals[1])
    errors.append(vals[2])
hs = np.array(hs)
errors = np.array(errors)

plt.loglog(hs, errors, color="k", linewidth=2)
plt.loglog(hs, hs, color="b")
plt.loglog(hs, hs**2, color="r")
plt.scatter(hs, errors)
plt.legend([r"$\log(error)$", r"$\log(h)$", r"$\log(h^2)$"])
plt.xlabel(r"$\log(h)$")
plt.title("Log-log plot $L^2$ error of the linear triangulation FVM for Poisson's equation.")
plt.show()
