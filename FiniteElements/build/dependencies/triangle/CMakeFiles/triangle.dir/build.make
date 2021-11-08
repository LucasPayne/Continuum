# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.20

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:

#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:

# Disable VCS-based implicit rules.
% : %,v

# Disable VCS-based implicit rules.
% : RCS/%

# Disable VCS-based implicit rules.
% : RCS/%,v

# Disable VCS-based implicit rules.
% : SCCS/s.%

# Disable VCS-based implicit rules.
% : s.%

.SUFFIXES: .hpux_make_needs_suffix_list

# Command-line flag to silence nested $(MAKE).
$(VERBOSE)MAKESILENT = -s

#Suppress display of executed commands.
$(VERBOSE).SILENT:

# A target that is always out of date.
cmake_force:
.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /usr/local/bin/cmake

# The command to remove a file.
RM = /usr/local/bin/cmake -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/lucas/Continuum/FiniteElements

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/lucas/Continuum/FiniteElements/build

# Include any dependencies generated for this target.
include dependencies/triangle/CMakeFiles/triangle.dir/depend.make
# Include any dependencies generated by the compiler for this target.
include dependencies/triangle/CMakeFiles/triangle.dir/compiler_depend.make

# Include the progress variables for this target.
include dependencies/triangle/CMakeFiles/triangle.dir/progress.make

# Include the compile flags for this target's objects.
include dependencies/triangle/CMakeFiles/triangle.dir/flags.make

dependencies/triangle/CMakeFiles/triangle.dir/triangle/triangle.c.o: dependencies/triangle/CMakeFiles/triangle.dir/flags.make
dependencies/triangle/CMakeFiles/triangle.dir/triangle/triangle.c.o: ../dependencies/triangle/triangle/triangle.c
dependencies/triangle/CMakeFiles/triangle.dir/triangle/triangle.c.o: dependencies/triangle/CMakeFiles/triangle.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/lucas/Continuum/FiniteElements/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building C object dependencies/triangle/CMakeFiles/triangle.dir/triangle/triangle.c.o"
	cd /home/lucas/Continuum/FiniteElements/build/dependencies/triangle && /usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -MD -MT dependencies/triangle/CMakeFiles/triangle.dir/triangle/triangle.c.o -MF CMakeFiles/triangle.dir/triangle/triangle.c.o.d -o CMakeFiles/triangle.dir/triangle/triangle.c.o -c /home/lucas/Continuum/FiniteElements/dependencies/triangle/triangle/triangle.c

dependencies/triangle/CMakeFiles/triangle.dir/triangle/triangle.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/triangle.dir/triangle/triangle.c.i"
	cd /home/lucas/Continuum/FiniteElements/build/dependencies/triangle && /usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /home/lucas/Continuum/FiniteElements/dependencies/triangle/triangle/triangle.c > CMakeFiles/triangle.dir/triangle/triangle.c.i

dependencies/triangle/CMakeFiles/triangle.dir/triangle/triangle.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/triangle.dir/triangle/triangle.c.s"
	cd /home/lucas/Continuum/FiniteElements/build/dependencies/triangle && /usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /home/lucas/Continuum/FiniteElements/dependencies/triangle/triangle/triangle.c -o CMakeFiles/triangle.dir/triangle/triangle.c.s

# Object files for target triangle
triangle_OBJECTS = \
"CMakeFiles/triangle.dir/triangle/triangle.c.o"

# External object files for target triangle
triangle_EXTERNAL_OBJECTS =

dependencies/triangle/libtriangle.a: dependencies/triangle/CMakeFiles/triangle.dir/triangle/triangle.c.o
dependencies/triangle/libtriangle.a: dependencies/triangle/CMakeFiles/triangle.dir/build.make
dependencies/triangle/libtriangle.a: dependencies/triangle/CMakeFiles/triangle.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/lucas/Continuum/FiniteElements/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking C static library libtriangle.a"
	cd /home/lucas/Continuum/FiniteElements/build/dependencies/triangle && $(CMAKE_COMMAND) -P CMakeFiles/triangle.dir/cmake_clean_target.cmake
	cd /home/lucas/Continuum/FiniteElements/build/dependencies/triangle && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/triangle.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
dependencies/triangle/CMakeFiles/triangle.dir/build: dependencies/triangle/libtriangle.a
.PHONY : dependencies/triangle/CMakeFiles/triangle.dir/build

dependencies/triangle/CMakeFiles/triangle.dir/clean:
	cd /home/lucas/Continuum/FiniteElements/build/dependencies/triangle && $(CMAKE_COMMAND) -P CMakeFiles/triangle.dir/cmake_clean.cmake
.PHONY : dependencies/triangle/CMakeFiles/triangle.dir/clean

dependencies/triangle/CMakeFiles/triangle.dir/depend:
	cd /home/lucas/Continuum/FiniteElements/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/lucas/Continuum/FiniteElements /home/lucas/Continuum/FiniteElements/dependencies/triangle /home/lucas/Continuum/FiniteElements/build /home/lucas/Continuum/FiniteElements/build/dependencies/triangle /home/lucas/Continuum/FiniteElements/build/dependencies/triangle/CMakeFiles/triangle.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : dependencies/triangle/CMakeFiles/triangle.dir/depend

