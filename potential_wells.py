from matplotlib import pyplot as plt
import sympy as sym
import numpy as np
from math import sin, cos, sqrt


def draw_phase_curves(f):
    extent = 3
    N = 300
    X,Y = np.mgrid[-extent:extent:1j*N, -extent:extent:1j*N]

    fig, axs = plt.subplots(1, 2)
    nodes = np.linspace(-extent, extent, N)
    axs[0].plot(nodes, f(nodes))
    axs[0].set_ylim(-extent, extent)
    axs[0].plot([-extent, extent], [0,0], "k")
    axs[1].contour(X,Y, f(X) + Y**2/2, 30)
    axs[1].plot([-extent, extent], [0,0], "k")
    axs[1].plot([0,0], [-extent, extent], "k")
    plt.show()


draw_phase_curves(lambda x: 4*x**3)
draw_phase_curves(lambda x: x**3 - x**2)
draw_phase_curves(lambda x: -np.sin(5*x))

