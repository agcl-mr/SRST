# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.5

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
CMAKE_COMMAND = /home/raman/Downloads/cmake-3.5.2-Linux-x86_64/bin/cmake

# The command to remove a file.
RM = /home/raman/Downloads/cmake-3.5.2-Linux-x86_64/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = "/media/raman/My Files/SRST with PP"

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = "/media/raman/My Files/SRST with PP"

# Include any dependencies generated for this target.
include CMakeFiles/srst.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/srst.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/srst.dir/flags.make

CMakeFiles/srst.dir/srst.cpp.o: CMakeFiles/srst.dir/flags.make
CMakeFiles/srst.dir/srst.cpp.o: srst.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir="/media/raman/My Files/SRST with PP/CMakeFiles" --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/srst.dir/srst.cpp.o"
	/usr/bin/c++   $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/srst.dir/srst.cpp.o -c "/media/raman/My Files/SRST with PP/srst.cpp"

CMakeFiles/srst.dir/srst.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/srst.dir/srst.cpp.i"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E "/media/raman/My Files/SRST with PP/srst.cpp" > CMakeFiles/srst.dir/srst.cpp.i

CMakeFiles/srst.dir/srst.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/srst.dir/srst.cpp.s"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S "/media/raman/My Files/SRST with PP/srst.cpp" -o CMakeFiles/srst.dir/srst.cpp.s

CMakeFiles/srst.dir/srst.cpp.o.requires:

.PHONY : CMakeFiles/srst.dir/srst.cpp.o.requires

CMakeFiles/srst.dir/srst.cpp.o.provides: CMakeFiles/srst.dir/srst.cpp.o.requires
	$(MAKE) -f CMakeFiles/srst.dir/build.make CMakeFiles/srst.dir/srst.cpp.o.provides.build
.PHONY : CMakeFiles/srst.dir/srst.cpp.o.provides

CMakeFiles/srst.dir/srst.cpp.o.provides.build: CMakeFiles/srst.dir/srst.cpp.o


# Object files for target srst
srst_OBJECTS = \
"CMakeFiles/srst.dir/srst.cpp.o"

# External object files for target srst
srst_EXTERNAL_OBJECTS =

srst: CMakeFiles/srst.dir/srst.cpp.o
srst: CMakeFiles/srst.dir/build.make
srst: /usr/lib/x86_64-linux-gnu/libmpfr.so
srst: /usr/lib/x86_64-linux-gnu/libgmp.so
srst: /usr/lib/x86_64-linux-gnu/libCGAL_Core.so.11.0.1
srst: /usr/lib/x86_64-linux-gnu/libCGAL.so.11.0.1
srst: /usr/lib/x86_64-linux-gnu/libboost_thread.so
srst: /usr/lib/x86_64-linux-gnu/libboost_system.so
srst: /usr/lib/x86_64-linux-gnu/libpthread.so
srst: /usr/lib/x86_64-linux-gnu/libCGAL_Core.so.11.0.1
srst: /usr/lib/x86_64-linux-gnu/libCGAL.so.11.0.1
srst: /usr/lib/x86_64-linux-gnu/libboost_thread.so
srst: /usr/lib/x86_64-linux-gnu/libboost_system.so
srst: /usr/lib/x86_64-linux-gnu/libpthread.so
srst: CMakeFiles/srst.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir="/media/raman/My Files/SRST with PP/CMakeFiles" --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable srst"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/srst.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/srst.dir/build: srst

.PHONY : CMakeFiles/srst.dir/build

CMakeFiles/srst.dir/requires: CMakeFiles/srst.dir/srst.cpp.o.requires

.PHONY : CMakeFiles/srst.dir/requires

CMakeFiles/srst.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/srst.dir/cmake_clean.cmake
.PHONY : CMakeFiles/srst.dir/clean

CMakeFiles/srst.dir/depend:
	cd "/media/raman/My Files/SRST with PP" && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" "/media/raman/My Files/SRST with PP" "/media/raman/My Files/SRST with PP" "/media/raman/My Files/SRST with PP" "/media/raman/My Files/SRST with PP" "/media/raman/My Files/SRST with PP/CMakeFiles/srst.dir/DependInfo.cmake" --color=$(COLOR)
.PHONY : CMakeFiles/srst.dir/depend

