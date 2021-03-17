# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.10

# Default target executed when no arguments are given to make.
default_target: all

.PHONY : default_target

# Allow only one "make -f Makefile2" at a time, but pass parallelism.
.NOTPARALLEL:


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
CMAKE_SOURCE_DIR = /home/ruben/ne_ooc/ne_ooc

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/ruben/ne_ooc/ne_ooc

#=============================================================================
# Targets provided globally by CMake.

# Special rule for the target rebuild_cache
rebuild_cache:
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --cyan "Running CMake to regenerate build system..."
	/usr/bin/cmake -H$(CMAKE_SOURCE_DIR) -B$(CMAKE_BINARY_DIR)
.PHONY : rebuild_cache

# Special rule for the target rebuild_cache
rebuild_cache/fast: rebuild_cache

.PHONY : rebuild_cache/fast

# Special rule for the target edit_cache
edit_cache:
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --cyan "No interactive CMake dialog available..."
	/usr/bin/cmake -E echo No\ interactive\ CMake\ dialog\ available.
.PHONY : edit_cache

# Special rule for the target edit_cache
edit_cache/fast: edit_cache

.PHONY : edit_cache/fast

# The main all target
all: cmake_check_build_system
	$(CMAKE_COMMAND) -E cmake_progress_start /home/ruben/ne_ooc/ne_ooc/CMakeFiles /home/ruben/ne_ooc/ne_ooc/CMakeFiles/progress.marks
	$(MAKE) -f CMakeFiles/Makefile2 all
	$(CMAKE_COMMAND) -E cmake_progress_start /home/ruben/ne_ooc/ne_ooc/CMakeFiles 0
.PHONY : all

# The main clean target
clean:
	$(MAKE) -f CMakeFiles/Makefile2 clean
.PHONY : clean

# The main clean target
clean/fast: clean

.PHONY : clean/fast

# Prepare targets for installation.
preinstall: all
	$(MAKE) -f CMakeFiles/Makefile2 preinstall
.PHONY : preinstall

# Prepare targets for installation.
preinstall/fast:
	$(MAKE) -f CMakeFiles/Makefile2 preinstall
.PHONY : preinstall/fast

# clear depends
depend:
	$(CMAKE_COMMAND) -H$(CMAKE_SOURCE_DIR) -B$(CMAKE_BINARY_DIR) --check-build-system CMakeFiles/Makefile.cmake 1
.PHONY : depend

#=============================================================================
# Target rules for targets named main

# Build rule for target.
main: cmake_check_build_system
	$(MAKE) -f CMakeFiles/Makefile2 main
.PHONY : main

# fast build rule for target.
main/fast:
	$(MAKE) -f CMakeFiles/main.dir/build.make CMakeFiles/main.dir/build
.PHONY : main/fast

src/Vertex2EdgePart.o: src/Vertex2EdgePart.cpp.o

.PHONY : src/Vertex2EdgePart.o

# target to build an object file
src/Vertex2EdgePart.cpp.o:
	$(MAKE) -f CMakeFiles/main.dir/build.make CMakeFiles/main.dir/src/Vertex2EdgePart.cpp.o
.PHONY : src/Vertex2EdgePart.cpp.o

src/Vertex2EdgePart.i: src/Vertex2EdgePart.cpp.i

.PHONY : src/Vertex2EdgePart.i

# target to preprocess a source file
src/Vertex2EdgePart.cpp.i:
	$(MAKE) -f CMakeFiles/main.dir/build.make CMakeFiles/main.dir/src/Vertex2EdgePart.cpp.i
.PHONY : src/Vertex2EdgePart.cpp.i

src/Vertex2EdgePart.s: src/Vertex2EdgePart.cpp.s

.PHONY : src/Vertex2EdgePart.s

# target to generate assembly for a file
src/Vertex2EdgePart.cpp.s:
	$(MAKE) -f CMakeFiles/main.dir/build.make CMakeFiles/main.dir/src/Vertex2EdgePart.cpp.s
.PHONY : src/Vertex2EdgePart.cpp.s

src/conversions.o: src/conversions.cpp.o

.PHONY : src/conversions.o

# target to build an object file
src/conversions.cpp.o:
	$(MAKE) -f CMakeFiles/main.dir/build.make CMakeFiles/main.dir/src/conversions.cpp.o
.PHONY : src/conversions.cpp.o

src/conversions.i: src/conversions.cpp.i

.PHONY : src/conversions.i

# target to preprocess a source file
src/conversions.cpp.i:
	$(MAKE) -f CMakeFiles/main.dir/build.make CMakeFiles/main.dir/src/conversions.cpp.i
.PHONY : src/conversions.cpp.i

src/conversions.s: src/conversions.cpp.s

.PHONY : src/conversions.s

# target to generate assembly for a file
src/conversions.cpp.s:
	$(MAKE) -f CMakeFiles/main.dir/build.make CMakeFiles/main.dir/src/conversions.cpp.s
.PHONY : src/conversions.cpp.s

src/graph.o: src/graph.cpp.o

.PHONY : src/graph.o

# target to build an object file
src/graph.cpp.o:
	$(MAKE) -f CMakeFiles/main.dir/build.make CMakeFiles/main.dir/src/graph.cpp.o
.PHONY : src/graph.cpp.o

src/graph.i: src/graph.cpp.i

.PHONY : src/graph.i

# target to preprocess a source file
src/graph.cpp.i:
	$(MAKE) -f CMakeFiles/main.dir/build.make CMakeFiles/main.dir/src/graph.cpp.i
.PHONY : src/graph.cpp.i

src/graph.s: src/graph.cpp.s

.PHONY : src/graph.s

# target to generate assembly for a file
src/graph.cpp.s:
	$(MAKE) -f CMakeFiles/main.dir/build.make CMakeFiles/main.dir/src/graph.cpp.s
.PHONY : src/graph.cpp.s

src/hep_partitioner.o: src/hep_partitioner.cpp.o

.PHONY : src/hep_partitioner.o

# target to build an object file
src/hep_partitioner.cpp.o:
	$(MAKE) -f CMakeFiles/main.dir/build.make CMakeFiles/main.dir/src/hep_partitioner.cpp.o
.PHONY : src/hep_partitioner.cpp.o

src/hep_partitioner.i: src/hep_partitioner.cpp.i

.PHONY : src/hep_partitioner.i

# target to preprocess a source file
src/hep_partitioner.cpp.i:
	$(MAKE) -f CMakeFiles/main.dir/build.make CMakeFiles/main.dir/src/hep_partitioner.cpp.i
.PHONY : src/hep_partitioner.cpp.i

src/hep_partitioner.s: src/hep_partitioner.cpp.s

.PHONY : src/hep_partitioner.s

# target to generate assembly for a file
src/hep_partitioner.cpp.s:
	$(MAKE) -f CMakeFiles/main.dir/build.make CMakeFiles/main.dir/src/hep_partitioner.cpp.s
.PHONY : src/hep_partitioner.cpp.s

src/main.o: src/main.cpp.o

.PHONY : src/main.o

# target to build an object file
src/main.cpp.o:
	$(MAKE) -f CMakeFiles/main.dir/build.make CMakeFiles/main.dir/src/main.cpp.o
.PHONY : src/main.cpp.o

src/main.i: src/main.cpp.i

.PHONY : src/main.i

# target to preprocess a source file
src/main.cpp.i:
	$(MAKE) -f CMakeFiles/main.dir/build.make CMakeFiles/main.dir/src/main.cpp.i
.PHONY : src/main.cpp.i

src/main.s: src/main.cpp.s

.PHONY : src/main.s

# target to generate assembly for a file
src/main.cpp.s:
	$(MAKE) -f CMakeFiles/main.dir/build.make CMakeFiles/main.dir/src/main.cpp.s
.PHONY : src/main.cpp.s

src/ne_graph.o: src/ne_graph.cpp.o

.PHONY : src/ne_graph.o

# target to build an object file
src/ne_graph.cpp.o:
	$(MAKE) -f CMakeFiles/main.dir/build.make CMakeFiles/main.dir/src/ne_graph.cpp.o
.PHONY : src/ne_graph.cpp.o

src/ne_graph.i: src/ne_graph.cpp.i

.PHONY : src/ne_graph.i

# target to preprocess a source file
src/ne_graph.cpp.i:
	$(MAKE) -f CMakeFiles/main.dir/build.make CMakeFiles/main.dir/src/ne_graph.cpp.i
.PHONY : src/ne_graph.cpp.i

src/ne_graph.s: src/ne_graph.cpp.s

.PHONY : src/ne_graph.s

# target to generate assembly for a file
src/ne_graph.cpp.s:
	$(MAKE) -f CMakeFiles/main.dir/build.make CMakeFiles/main.dir/src/ne_graph.cpp.s
.PHONY : src/ne_graph.cpp.s

src/ne_partitioner.o: src/ne_partitioner.cpp.o

.PHONY : src/ne_partitioner.o

# target to build an object file
src/ne_partitioner.cpp.o:
	$(MAKE) -f CMakeFiles/main.dir/build.make CMakeFiles/main.dir/src/ne_partitioner.cpp.o
.PHONY : src/ne_partitioner.cpp.o

src/ne_partitioner.i: src/ne_partitioner.cpp.i

.PHONY : src/ne_partitioner.i

# target to preprocess a source file
src/ne_partitioner.cpp.i:
	$(MAKE) -f CMakeFiles/main.dir/build.make CMakeFiles/main.dir/src/ne_partitioner.cpp.i
.PHONY : src/ne_partitioner.cpp.i

src/ne_partitioner.s: src/ne_partitioner.cpp.s

.PHONY : src/ne_partitioner.s

# target to generate assembly for a file
src/ne_partitioner.cpp.s:
	$(MAKE) -f CMakeFiles/main.dir/build.make CMakeFiles/main.dir/src/ne_partitioner.cpp.s
.PHONY : src/ne_partitioner.cpp.s

src/util.o: src/util.cpp.o

.PHONY : src/util.o

# target to build an object file
src/util.cpp.o:
	$(MAKE) -f CMakeFiles/main.dir/build.make CMakeFiles/main.dir/src/util.cpp.o
.PHONY : src/util.cpp.o

src/util.i: src/util.cpp.i

.PHONY : src/util.i

# target to preprocess a source file
src/util.cpp.i:
	$(MAKE) -f CMakeFiles/main.dir/build.make CMakeFiles/main.dir/src/util.cpp.i
.PHONY : src/util.cpp.i

src/util.s: src/util.cpp.s

.PHONY : src/util.s

# target to generate assembly for a file
src/util.cpp.s:
	$(MAKE) -f CMakeFiles/main.dir/build.make CMakeFiles/main.dir/src/util.cpp.s
.PHONY : src/util.cpp.s

# Help Target
help:
	@echo "The following are some of the valid targets for this Makefile:"
	@echo "... all (the default if no target is provided)"
	@echo "... clean"
	@echo "... depend"
	@echo "... rebuild_cache"
	@echo "... main"
	@echo "... edit_cache"
	@echo "... src/Vertex2EdgePart.o"
	@echo "... src/Vertex2EdgePart.i"
	@echo "... src/Vertex2EdgePart.s"
	@echo "... src/conversions.o"
	@echo "... src/conversions.i"
	@echo "... src/conversions.s"
	@echo "... src/graph.o"
	@echo "... src/graph.i"
	@echo "... src/graph.s"
	@echo "... src/hep_partitioner.o"
	@echo "... src/hep_partitioner.i"
	@echo "... src/hep_partitioner.s"
	@echo "... src/main.o"
	@echo "... src/main.i"
	@echo "... src/main.s"
	@echo "... src/ne_graph.o"
	@echo "... src/ne_graph.i"
	@echo "... src/ne_graph.s"
	@echo "... src/ne_partitioner.o"
	@echo "... src/ne_partitioner.i"
	@echo "... src/ne_partitioner.s"
	@echo "... src/util.o"
	@echo "... src/util.i"
	@echo "... src/util.s"
.PHONY : help



#=============================================================================
# Special targets to cleanup operation of make.

# Special rule to run CMake to check the build system integrity.
# No rule that depends on this can have commands that come from listfiles
# because they might be regenerated.
cmake_check_build_system:
	$(CMAKE_COMMAND) -H$(CMAKE_SOURCE_DIR) -B$(CMAKE_BINARY_DIR) --check-build-system CMakeFiles/Makefile.cmake 0
.PHONY : cmake_check_build_system
