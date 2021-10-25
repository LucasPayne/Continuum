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

add_executable(QuadraticLaplaceSolver QuadraticLaplaceSolver.cpp)
target_link_libraries(QuadraticLaplaceSolver PRIVATE common)

add_executable(LidDrivenCavity
    LidDrivenCavity/LidDrivenCavity.cpp
)
target_link_libraries(LidDrivenCavity PRIVATE common)

add_executable(WeaklyIncompressible
    WeaklyIncompressible/WeaklyIncompressible.cpp
)
target_link_libraries(WeaklyIncompressible PRIVATE common)


add_executable(NavierStokes
    NavierStokes/main.cpp

    NavierStokes/demo.cpp
    NavierStokes/mesh_generators.cpp
    NavierStokes/CameraController.cpp

    NavierStokes/NavierStokesSolver.cpp
    NavierStokes/compute_velocity_residual.cpp
    NavierStokes/compute_pressure_residual.cpp
    NavierStokes/compute_gateaux_matrix_top_left.cpp
    NavierStokes/compute_gateaux_matrix_bottom_left.cpp
)
target_link_libraries(NavierStokes PRIVATE common)
