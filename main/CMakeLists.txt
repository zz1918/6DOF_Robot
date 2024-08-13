# This is the cmake file for SE3 framework on the Delta robot
cmake_minimum_required(VERSION 3.24.0)
project(main)
add_executable(main main.cpp)
find_package(Eigen3 3.3 REQUIRED NO_MODULE)
target_link_libraries(main PUBLIC Eigen3::Eigen)

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
