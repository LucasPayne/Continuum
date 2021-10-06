import numpy as np
from matplotlib import pyplot as plt
import sys
from pathlib import Path
import os
import subprocess


fig = plt.gcf()
fig.set_size_inches(8, 3)

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


plt.loglog(hs, errors, color="k", linewidth=1)



plt.tight_layout(pad=3)
style = {"linestyle": "--", "linewidth": 0.8}
plt.loglog(hs, (errors[-1]/hs[-1]) * hs, color="b", **style)
plt.loglog(hs, (errors[-1]/(hs[-1]**2)) * hs**2, color="r", **style)
plt.scatter(hs, errors)
plt.legend([r"$\log(error)$", r"$\log(r)$ + const", r"$\log(r^2)$ + const"])
plt.xlabel(r"$r$ (labels), $\log(r)$ (scale)")
plt.title("Log-log plot of the $L^2$ error for varying $r$.")

figures_directory = Path("/home/lucas/Continuum/document/figures")
fig_path = figures_directory / Path("graphs") / Path(Path(filename).name).with_suffix(".png")
plt.savefig(fig_path, dpi=300)
subprocess.call(["convert", fig_path, "-trim", fig_path])
