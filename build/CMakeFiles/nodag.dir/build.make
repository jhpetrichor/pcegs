# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.22

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
CMAKE_COMMAND = /usr/bin/cmake

# The command to remove a file.
RM = /usr/bin/cmake -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/jh/JHPC

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/jh/JHPC/build

# Include any dependencies generated for this target.
include CMakeFiles/nodag.dir/depend.make
# Include any dependencies generated by the compiler for this target.
include CMakeFiles/nodag.dir/compiler_depend.make

# Include the progress variables for this target.
include CMakeFiles/nodag.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/nodag.dir/flags.make

CMakeFiles/nodag.dir/src/no_dag.cpp.o: CMakeFiles/nodag.dir/flags.make
CMakeFiles/nodag.dir/src/no_dag.cpp.o: ../src/no_dag.cpp
CMakeFiles/nodag.dir/src/no_dag.cpp.o: CMakeFiles/nodag.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/jh/JHPC/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/nodag.dir/src/no_dag.cpp.o"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/nodag.dir/src/no_dag.cpp.o -MF CMakeFiles/nodag.dir/src/no_dag.cpp.o.d -o CMakeFiles/nodag.dir/src/no_dag.cpp.o -c /home/jh/JHPC/src/no_dag.cpp

CMakeFiles/nodag.dir/src/no_dag.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/nodag.dir/src/no_dag.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/jh/JHPC/src/no_dag.cpp > CMakeFiles/nodag.dir/src/no_dag.cpp.i

CMakeFiles/nodag.dir/src/no_dag.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/nodag.dir/src/no_dag.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/jh/JHPC/src/no_dag.cpp -o CMakeFiles/nodag.dir/src/no_dag.cpp.s

# Object files for target nodag
nodag_OBJECTS = \
"CMakeFiles/nodag.dir/src/no_dag.cpp.o"

# External object files for target nodag
nodag_EXTERNAL_OBJECTS =

../bin/nodag: CMakeFiles/nodag.dir/src/no_dag.cpp.o
../bin/nodag: CMakeFiles/nodag.dir/build.make
../bin/nodag: ../lib/libcpdp.so
../bin/nodag: CMakeFiles/nodag.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/jh/JHPC/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable ../bin/nodag"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/nodag.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/nodag.dir/build: ../bin/nodag
.PHONY : CMakeFiles/nodag.dir/build

CMakeFiles/nodag.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/nodag.dir/cmake_clean.cmake
.PHONY : CMakeFiles/nodag.dir/clean

CMakeFiles/nodag.dir/depend:
	cd /home/jh/JHPC/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/jh/JHPC /home/jh/JHPC /home/jh/JHPC/build /home/jh/JHPC/build /home/jh/JHPC/build/CMakeFiles/nodag.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/nodag.dir/depend

