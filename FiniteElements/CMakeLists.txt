cmake_minimum_required(VERSION 3.1)
project(FiniteElements)

# cg_sandbox
set(cg_path /home/lucas/computer_graphics/cg_sandbox)
add_subdirectory(${cg_path} ${cg_path}/build EXCLUDE_FROM_ALL)
set(cg_include_path ${cg_path}/include ${cg_path}/lib)

# mesh_processing
set(mp_path /home/lucas/computer_graphics/mesh_processing)
set(mp_include_path ${mp_path}/include)

# Eigen3
set(eigen3_include_path /usr/include/eigen3)

# Triangle
set(triangle_path ${CMAKE_CURRENT_SOURCE_DIR}/dependencies/triangle)
add_compile_options(-DANSI_DECLARATORS -DTRILIBRARY) # So the header is included correctly.
# Note: REAL and VOID types must be specified in the application.
add_subdirectory(${triangle_path})

message(STATUS "triangle: " ${triangle_path})
message(STATUS "cg_sandbox include: " ${cg_include_path})
message(STATUS "mesh_processing include: " ${mp_include_path})

include_directories(
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${cg_include_path}
    ${mp_include_path}
    ${eigen3_include_path}
    ${triangle_path}
)

add_library(common INTERFACE)
target_compile_definitions(common INTERFACE
    DATA="${CMAKE_CURRENT_SOURCE_DIR}/data/"
    MODELS="${CMAKE_CURRENT_SOURCE_DIR}/models/"
    SHADERS="${CMAKE_CURRENT_SOURCE_DIR}/shaders/"
)
set(CMAKE_CXX_FLAGS "-g")
target_link_libraries(common INTERFACE
    triangle
    cg_sandbox
)

add_executable(LaplaceSolver LaplaceSolver.cpp)
target_link_libraries(LaplaceSolver PRIVATE common)

add_executable(HeatSolver HeatSolver.cpp)
target_link_libraries(HeatSolver PRIVATE common)

add_executable(discretization_figure discretization_figure.cpp)
target_link_libraries(discretization_figure PRIVATE common)

add_executable(FiniteVolumes FiniteVolumes.cpp)
target_link_libraries(FiniteVolumes PRIVATE common)

add_executable(quadratic_figure quadratic_figure.cpp)
target_link_libraries(quadratic_figure PRIVATE common)

add_executable(basis_quadratic_figure basis_quadratic_figure.cpp)
target_link_libraries(basis_quadratic_figure PRIVATE common)

add_executable(QuadraticLaplaceSolver
    QuadraticLaplaceSolver.cpp
    core.cpp
)
target_link_libraries(QuadraticLaplaceSolver PRIVATE common)

add_executable(LidDrivenCavity
    LidDrivenCavity/LidDrivenCavity.cpp
)
target_link_libraries(LidDrivenCavity PRIVATE common)

add_executable(WeaklyIncompressible
    WeaklyIncompressible/WeaklyIncompressible.cpp
    core.cpp
)
target_link_libraries(WeaklyIncompressible PRIVATE common)


add_executable(NavierStokes
    NavierStokes/main.cpp

    NavierStokes/demo.cpp
    CameraController.cpp
    core.cpp

    NavierStokes/NavierStokesSolver.cpp
    NavierStokes/add_nonlinear_velocity_residual.cpp
    NavierStokes/compute_linear_term_matrix_top_left.cpp
    NavierStokes/compute_linear_term_matrix_bottom_left.cpp
)
target_link_libraries(NavierStokes PRIVATE common)

add_executable(SurfaceNavierStokes
    SurfaceNavierStokes/main.cpp

    SurfaceNavierStokes/demo.cpp
    CameraController.cpp
    core.cpp

    SurfaceNavierStokes/SurfaceNavierStokesSolver.cpp
    SurfaceNavierStokes/compute_velocity_block_coefficients.cpp
    SurfaceNavierStokes/compute_pressure_block_coefficients.cpp
    SurfaceNavierStokes/compute_centripetal_block_coefficients.cpp
    SurfaceNavierStokes/explicit_advection.cpp

    SurfaceNavierStokes/make_sparsity_image.cpp
)
target_link_libraries(SurfaceNavierStokes PRIVATE common)

add_executable(ReactionDiffusion
    ReactionDiffusion/main.cpp

    ReactionDiffusion/demo.cpp
    ReactionDiffusion/CameraController.cpp

    ReactionDiffusion/ReactionDiffusion.cpp
)
target_link_libraries(ReactionDiffusion PRIVATE common)


add_executable(FlowVisualization
    FlowVisualization/main.cpp
    CameraController.cpp
    core.cpp
    Filmer.cpp
)
target_link_libraries(FlowVisualization PRIVATE common)
