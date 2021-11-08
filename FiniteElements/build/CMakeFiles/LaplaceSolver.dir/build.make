# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.16

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:


#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:


# Remove some rules from gmake that .SUFFIXES does not remove.
SUFFIXES =

.SUFFIXES: .hpux_make_needs_suffix_list


# Suppress display of executed commands.
$(VERBOSE).SILENT:


# A target that is always out of date.
cmake_force:

.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /usr/bin/cmake

# The command to remove a file.
RM = /usr/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/lucas/Continuum/FiniteElements

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/lucas/Continuum/FiniteElements/build

# Include any dependencies generated for this target.
include CMakeFiles/LaplaceSolver.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/LaplaceSolver.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/LaplaceSolver.dir/flags.make

CMakeFiles/LaplaceSolver.dir/LaplaceSolver.cpp.o: CMakeFiles/LaplaceSolver.dir/flags.make
CMakeFiles/LaplaceSolver.dir/LaplaceSolver.cpp.o: ../LaplaceSolver.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/lucas/Continuum/FiniteElements/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/LaplaceSolver.dir/LaplaceSolver.cpp.o"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/LaplaceSolver.dir/LaplaceSolver.cpp.o -c /home/lucas/Continuum/FiniteElements/LaplaceSolver.cpp

CMakeFiles/LaplaceSolver.dir/LaplaceSolver.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/LaplaceSolver.dir/LaplaceSolver.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/lucas/Continuum/FiniteElements/LaplaceSolver.cpp > CMakeFiles/LaplaceSolver.dir/LaplaceSolver.cpp.i

CMakeFiles/LaplaceSolver.dir/LaplaceSolver.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/LaplaceSolver.dir/LaplaceSolver.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/lucas/Continuum/FiniteElements/LaplaceSolver.cpp -o CMakeFiles/LaplaceSolver.dir/LaplaceSolver.cpp.s

# Object files for target LaplaceSolver
LaplaceSolver_OBJECTS = \
"CMakeFiles/LaplaceSolver.dir/LaplaceSolver.cpp.o"

# External object files for target LaplaceSolver
LaplaceSolver_EXTERNAL_OBJECTS =

LaplaceSolver: CMakeFiles/LaplaceSolver.dir/LaplaceSolver.cpp.o
LaplaceSolver: CMakeFiles/LaplaceSolver.dir/build.make
LaplaceSolver: dependencies/triangle/libtriangle.a
LaplaceSolver: /home/lucas/computer_graphics/cg_sandbox/build/src/libcg_sandbox.a
LaplaceSolver: /home/lucas/computer_graphics/cg_sandbox/build/lib/opengl_utilities/libopengl_utilities.a
LaplaceSolver: /home/lucas/computer_graphics/cg_sandbox/build/lib/interactive_graphics_context/libinteractive_graphics_context.a
LaplaceSolver: /home/lucas/computer_graphics/cg_sandbox/build/lib/reflector/libreflector.a
LaplaceSolver: /home/lucas/computer_graphics/cg_sandbox/build/lib/table/libtable.a
LaplaceSolver: /home/lucas/computer_graphics/cg_sandbox/build/lib/mathematics/libmathematics.a
LaplaceSolver: /home/lucas/computer_graphics/cg_sandbox/build/lib/mathematics/spatial_algebra/libspatial_algebra.a
LaplaceSolver: /home/lucas/computer_graphics/mesh_processing/build/libmesh_processing.a
LaplaceSolver: /home/lucas/computer_graphics/cg_sandbox/build/dependencies/glfw/src/libglfw3.a
LaplaceSolver: /usr/lib/x86_64-linux-gnu/librt.so
LaplaceSolver: /usr/lib/x86_64-linux-gnu/libm.so
LaplaceSolver: /home/lucas/computer_graphics/assimp/build/bin/libassimp.so
LaplaceSolver: CMakeFiles/LaplaceSolver.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/lucas/Continuum/FiniteElements/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable LaplaceSolver"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/LaplaceSolver.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/LaplaceSolver.dir/build: LaplaceSolver

.PHONY : CMakeFiles/LaplaceSolver.dir/build

CMakeFiles/LaplaceSolver.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/LaplaceSolver.dir/cmake_clean.cmake
.PHONY : CMakeFiles/LaplaceSolver.dir/clean

CMakeFiles/LaplaceSolver.dir/depend:
	cd /home/lucas/Continuum/FiniteElements/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/lucas/Continuum/FiniteElements /home/lucas/Continuum/FiniteElements /home/lucas/Continuum/FiniteElements/build /home/lucas/Continuum/FiniteElements/build /home/lucas/Continuum/FiniteElements/build/CMakeFiles/LaplaceSolver.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/LaplaceSolver.dir/depend

