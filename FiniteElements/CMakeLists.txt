cmake_minimum_required(VERSION 3.1)
project(FiniteElements)


set(cg_path /home/lucas/computer_graphics/cg_sandbox)
add_subdirectory(${cg_path} ${cg_path}/build EXCLUDE_FROM_ALL)
set(cg_include_path ${cg_path}/include ${cg_path}/lib)

set(mp_path /home/lucas/computer_graphics/mesh_processing)
set(mp_include_path ${mp_path}/include)
set(mp_library_path ${mp_path}/build/libmesh_processing.a)

set(eigen3_include_path /usr/include/eigen3)

message("cg_sandbox include: " ${cg_include_path})
message("mesh_processing include: " ${mp_include_path})
message("mesh_processing library: " ${mp_library_path})

include_directories(
    ${cg_include_path}
    ${mp_include_path}
    ${eigen3_include_path}
)

add_executable(FiniteElements
    main.cpp
)

target_link_libraries(FiniteElements
    cg_sandbox
)
