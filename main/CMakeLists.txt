# This is the cmake file for SSS framework on the Delta robot
cmake_minimum_required(VERSION 3.24.0)
project(main)
add_executable(main main.cpp)

if(MINGW OR CYGWIN)
	add_definitions(-O3)
endif()

find_package(Eigen3 3.3 REQUIRED NO_MODULE)
target_link_libraries(main PUBLIC Eigen3::Eigen)

list(PREPEND CMAKE_MODULE_PATH ../libigl)

# Libigl
include(libigl)

# Enable the target igl::glfw
igl_include(glfw)

# Link igl (and the glfw module) to your project
target_link_libraries(main PUBLIC igl::glfw)

if (NOT TARGET SSS)
	add_subdirectory("../SSS" SSS)
endif()
target_link_libraries(main PUBLIC SSS)
target_include_directories(main PUBLIC "../SSS")

if (NOT TARGET ReadSSSCommand)
	add_subdirectory("../ReadSSSCommand" ReadSSSCommand)
endif()
target_link_libraries(main PUBLIC ReadSSSCommand)
target_include_directories(main PUBLIC "../ReadSSSCommand")

