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
include CMakeFiles/pcegs.dir/depend.make
# Include any dependencies generated by the compiler for this target.
include CMakeFiles/pcegs.dir/compiler_depend.make

# Include the progress variables for this target.
include CMakeFiles/pcegs.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/pcegs.dir/flags.make

CMakeFiles/pcegs.dir/src/pcegs.cpp.o: CMakeFiles/pcegs.dir/flags.make
CMakeFiles/pcegs.dir/src/pcegs.cpp.o: ../src/pcegs.cpp
CMakeFiles/pcegs.dir/src/pcegs.cpp.o: CMakeFiles/pcegs.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/jh/JHPC/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/pcegs.dir/src/pcegs.cpp.o"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/pcegs.dir/src/pcegs.cpp.o -MF CMakeFiles/pcegs.dir/src/pcegs.cpp.o.d -o CMakeFiles/pcegs.dir/src/pcegs.cpp.o -c /home/jh/JHPC/src/pcegs.cpp

CMakeFiles/pcegs.dir/src/pcegs.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/pcegs.dir/src/pcegs.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/jh/JHPC/src/pcegs.cpp > CMakeFiles/pcegs.dir/src/pcegs.cpp.i

CMakeFiles/pcegs.dir/src/pcegs.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/pcegs.dir/src/pcegs.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/jh/JHPC/src/pcegs.cpp -o CMakeFiles/pcegs.dir/src/pcegs.cpp.s

# Object files for target pcegs
pcegs_OBJECTS = \
"CMakeFiles/pcegs.dir/src/pcegs.cpp.o"

# External object files for target pcegs
pcegs_EXTERNAL_OBJECTS =

../bin/pcegs: CMakeFiles/pcegs.dir/src/pcegs.cpp.o
../bin/pcegs: CMakeFiles/pcegs.dir/build.make
../bin/pcegs: ../lib/libcpdp.so
../bin/pcegs: CMakeFiles/pcegs.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/jh/JHPC/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable ../bin/pcegs"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/pcegs.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/pcegs.dir/build: ../bin/pcegs
.PHONY : CMakeFiles/pcegs.dir/build

CMakeFiles/pcegs.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/pcegs.dir/cmake_clean.cmake
.PHONY : CMakeFiles/pcegs.dir/clean

CMakeFiles/pcegs.dir/depend:
	cd /home/jh/JHPC/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/jh/JHPC /home/jh/JHPC /home/jh/JHPC/build /home/jh/JHPC/build /home/jh/JHPC/build/CMakeFiles/pcegs.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/pcegs.dir/depend
