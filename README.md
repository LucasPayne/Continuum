## Continuum
------------------------
My undergraduate honours project report, completed November 2021.
The purpose of this project was to learn some continuum mechanics and fluid mechanics, and to implement numerical finite element
and finite volume solvers. The main files are

- report.pdf: An overview of concepts and the solvers (aimed as being a sort of tutorial and discussion).
- poster.pdf: A more brief poster presentation of the project (as an A1 wall poster).

The code is not very general, and is somewhat haphazard, as after a certain point I just wanted each application to create specific figures.
Each executable is an interactive 3D program which demonstrates a solver. The figures in the report and poster have been created
by flying a camera around and snapshots of the framebuffer. Here is a [video](https://www.youtube.com/watch?v=m5He2x3gNz8) showing some simulation results.

## Dependencies
-----------------------
For the integration loops I am using my SurfaceMesh data structure in [mesh_processing](https://github.com/LucasPayne/mesh_processing).
I am using [cg_sandbox](https://github.com/LucasPayne/cg_sandbox), which implements some basic concepts (entities, behaviours, a loop, etc.) which make
it simple to start implementing the actual demos.

# External dependencies
-----------------------
- [Eigen3 (C++ header-only linear algebra library.](https://gitlab.com/libeigen/eigen)
- [triangle (Delaunay triangulator. C source code, built in-tree)](http://www.cs.cmu.edu/~quake/triangle.html)
- [GLFW (OpenGL context, window, and system input)](https://www.glfw.org/)
