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
CMAKE_COMMAND = /opt/anaconda3/envs/altar/bin/cmake

# The command to remove a file.
RM = /opt/anaconda3/envs/altar/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /Users/bato/tools/altar

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /Users/bato/tools/altar/build

# Include any dependencies generated for this target.
include models/reverso/CMakeFiles/reversomodule.dir/depend.make

# Include the progress variables for this target.
include models/reverso/CMakeFiles/reversomodule.dir/progress.make

# Include the compile flags for this target's objects.
include models/reverso/CMakeFiles/reversomodule.dir/flags.make

models/reverso/CMakeFiles/reversomodule.dir/ext/reverso/reverso.cc.o: models/reverso/CMakeFiles/reversomodule.dir/flags.make
models/reverso/CMakeFiles/reversomodule.dir/ext/reverso/reverso.cc.o: ../models/reverso/ext/reverso/reverso.cc
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/bato/tools/altar/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object models/reverso/CMakeFiles/reversomodule.dir/ext/reverso/reverso.cc.o"
	cd /Users/bato/tools/altar/build/models/reverso && /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/reversomodule.dir/ext/reverso/reverso.cc.o -c /Users/bato/tools/altar/models/reverso/ext/reverso/reverso.cc

models/reverso/CMakeFiles/reversomodule.dir/ext/reverso/reverso.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/reversomodule.dir/ext/reverso/reverso.cc.i"
	cd /Users/bato/tools/altar/build/models/reverso && /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/bato/tools/altar/models/reverso/ext/reverso/reverso.cc > CMakeFiles/reversomodule.dir/ext/reverso/reverso.cc.i

models/reverso/CMakeFiles/reversomodule.dir/ext/reverso/reverso.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/reversomodule.dir/ext/reverso/reverso.cc.s"
	cd /Users/bato/tools/altar/build/models/reverso && /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/bato/tools/altar/models/reverso/ext/reverso/reverso.cc -o CMakeFiles/reversomodule.dir/ext/reverso/reverso.cc.s

models/reverso/CMakeFiles/reversomodule.dir/ext/reverso/exceptions.cc.o: models/reverso/CMakeFiles/reversomodule.dir/flags.make
models/reverso/CMakeFiles/reversomodule.dir/ext/reverso/exceptions.cc.o: ../models/reverso/ext/reverso/exceptions.cc
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/bato/tools/altar/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building CXX object models/reverso/CMakeFiles/reversomodule.dir/ext/reverso/exceptions.cc.o"
	cd /Users/bato/tools/altar/build/models/reverso && /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/reversomodule.dir/ext/reverso/exceptions.cc.o -c /Users/bato/tools/altar/models/reverso/ext/reverso/exceptions.cc

models/reverso/CMakeFiles/reversomodule.dir/ext/reverso/exceptions.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/reversomodule.dir/ext/reverso/exceptions.cc.i"
	cd /Users/bato/tools/altar/build/models/reverso && /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/bato/tools/altar/models/reverso/ext/reverso/exceptions.cc > CMakeFiles/reversomodule.dir/ext/reverso/exceptions.cc.i

models/reverso/CMakeFiles/reversomodule.dir/ext/reverso/exceptions.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/reversomodule.dir/ext/reverso/exceptions.cc.s"
	cd /Users/bato/tools/altar/build/models/reverso && /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/bato/tools/altar/models/reverso/ext/reverso/exceptions.cc -o CMakeFiles/reversomodule.dir/ext/reverso/exceptions.cc.s

models/reverso/CMakeFiles/reversomodule.dir/ext/reverso/source.cc.o: models/reverso/CMakeFiles/reversomodule.dir/flags.make
models/reverso/CMakeFiles/reversomodule.dir/ext/reverso/source.cc.o: ../models/reverso/ext/reverso/source.cc
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/bato/tools/altar/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Building CXX object models/reverso/CMakeFiles/reversomodule.dir/ext/reverso/source.cc.o"
	cd /Users/bato/tools/altar/build/models/reverso && /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/reversomodule.dir/ext/reverso/source.cc.o -c /Users/bato/tools/altar/models/reverso/ext/reverso/source.cc

models/reverso/CMakeFiles/reversomodule.dir/ext/reverso/source.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/reversomodule.dir/ext/reverso/source.cc.i"
	cd /Users/bato/tools/altar/build/models/reverso && /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/bato/tools/altar/models/reverso/ext/reverso/source.cc > CMakeFiles/reversomodule.dir/ext/reverso/source.cc.i

models/reverso/CMakeFiles/reversomodule.dir/ext/reverso/source.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/reversomodule.dir/ext/reverso/source.cc.s"
	cd /Users/bato/tools/altar/build/models/reverso && /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/bato/tools/altar/models/reverso/ext/reverso/source.cc -o CMakeFiles/reversomodule.dir/ext/reverso/source.cc.s

# Object files for target reversomodule
reversomodule_OBJECTS = \
"CMakeFiles/reversomodule.dir/ext/reverso/reverso.cc.o" \
"CMakeFiles/reversomodule.dir/ext/reverso/exceptions.cc.o" \
"CMakeFiles/reversomodule.dir/ext/reverso/source.cc.o"

# External object files for target reversomodule
reversomodule_EXTERNAL_OBJECTS =

models/reverso/reverso.cpython-38-darwin.so: models/reverso/CMakeFiles/reversomodule.dir/ext/reverso/reverso.cc.o
models/reverso/reverso.cpython-38-darwin.so: models/reverso/CMakeFiles/reversomodule.dir/ext/reverso/exceptions.cc.o
models/reverso/reverso.cpython-38-darwin.so: models/reverso/CMakeFiles/reversomodule.dir/ext/reverso/source.cc.o
models/reverso/reverso.cpython-38-darwin.so: models/reverso/CMakeFiles/reversomodule.dir/build.make
models/reverso/reverso.cpython-38-darwin.so: models/reverso/libreverso.dylib
models/reverso/reverso.cpython-38-darwin.so: altar/libaltar.dylib
models/reverso/reverso.cpython-38-darwin.so: models/reverso/CMakeFiles/reversomodule.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/Users/bato/tools/altar/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_4) "Linking CXX shared module reverso.cpython-38-darwin.so"
	cd /Users/bato/tools/altar/build/models/reverso && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/reversomodule.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
models/reverso/CMakeFiles/reversomodule.dir/build: models/reverso/reverso.cpython-38-darwin.so

.PHONY : models/reverso/CMakeFiles/reversomodule.dir/build

models/reverso/CMakeFiles/reversomodule.dir/clean:
	cd /Users/bato/tools/altar/build/models/reverso && $(CMAKE_COMMAND) -P CMakeFiles/reversomodule.dir/cmake_clean.cmake
.PHONY : models/reverso/CMakeFiles/reversomodule.dir/clean

models/reverso/CMakeFiles/reversomodule.dir/depend:
	cd /Users/bato/tools/altar/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /Users/bato/tools/altar /Users/bato/tools/altar/models/reverso /Users/bato/tools/altar/build /Users/bato/tools/altar/build/models/reverso /Users/bato/tools/altar/build/models/reverso/CMakeFiles/reversomodule.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : models/reverso/CMakeFiles/reversomodule.dir/depend
