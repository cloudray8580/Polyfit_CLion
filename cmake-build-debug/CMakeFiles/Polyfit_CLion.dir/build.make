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
CMAKE_SOURCE_DIR = /mnt/d/Polyfit_CLion

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /mnt/d/Polyfit_CLion/cmake-build-debug

# Include any dependencies generated for this target.
include CMakeFiles/Polyfit_CLion.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/Polyfit_CLion.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/Polyfit_CLion.dir/flags.make

CMakeFiles/Polyfit_CLion.dir/main.cpp.o: CMakeFiles/Polyfit_CLion.dir/flags.make
CMakeFiles/Polyfit_CLion.dir/main.cpp.o: ../main.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/mnt/d/Polyfit_CLion/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/Polyfit_CLion.dir/main.cpp.o"
	/usr/bin/g++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/Polyfit_CLion.dir/main.cpp.o -c /mnt/d/Polyfit_CLion/main.cpp

CMakeFiles/Polyfit_CLion.dir/main.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/Polyfit_CLion.dir/main.cpp.i"
	/usr/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /mnt/d/Polyfit_CLion/main.cpp > CMakeFiles/Polyfit_CLion.dir/main.cpp.i

CMakeFiles/Polyfit_CLion.dir/main.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/Polyfit_CLion.dir/main.cpp.s"
	/usr/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /mnt/d/Polyfit_CLion/main.cpp -o CMakeFiles/Polyfit_CLion.dir/main.cpp.s

# Object files for target Polyfit_CLion
Polyfit_CLion_OBJECTS = \
"CMakeFiles/Polyfit_CLion.dir/main.cpp.o"

# External object files for target Polyfit_CLion
Polyfit_CLion_EXTERNAL_OBJECTS =

Polyfit_CLion: CMakeFiles/Polyfit_CLion.dir/main.cpp.o
Polyfit_CLion: CMakeFiles/Polyfit_CLion.dir/build.make
Polyfit_CLion: /mnt/d/CPLEXLinux/cplex/lib/x86-64_linux/static_pic/libilocplex.a
Polyfit_CLion: /mnt/d/CPLEXLinux/concert/lib/x86-64_linux/static_pic/libconcert.a
Polyfit_CLion: /mnt/d/CPLEXLinux/cplex/lib/x86-64_linux/static_pic/libcplex.a
Polyfit_CLion: /mnt/d/xgboost/lib/libxgboost.so
Polyfit_CLion: CMakeFiles/Polyfit_CLion.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/mnt/d/Polyfit_CLion/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable Polyfit_CLion"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/Polyfit_CLion.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/Polyfit_CLion.dir/build: Polyfit_CLion

.PHONY : CMakeFiles/Polyfit_CLion.dir/build

CMakeFiles/Polyfit_CLion.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/Polyfit_CLion.dir/cmake_clean.cmake
.PHONY : CMakeFiles/Polyfit_CLion.dir/clean

CMakeFiles/Polyfit_CLion.dir/depend:
	cd /mnt/d/Polyfit_CLion/cmake-build-debug && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /mnt/d/Polyfit_CLion /mnt/d/Polyfit_CLion /mnt/d/Polyfit_CLion/cmake-build-debug /mnt/d/Polyfit_CLion/cmake-build-debug /mnt/d/Polyfit_CLion/cmake-build-debug/CMakeFiles/Polyfit_CLion.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/Polyfit_CLion.dir/depend

