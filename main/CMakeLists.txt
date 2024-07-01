# This is the cmake file for SE3 framework on the Delta robot
cmake_minimum_required(VERSION 3.24.0)
project(SE3)
add_executable(SE3 main.cpp)
find_package(Eigen3 3.3 REQUIRED NO_MODULE)
target_link_libraries(SE3 PUBLIC Eigen3::Eigen)

if (NOT TARGET SO3)
	add_subdirectory("../SO3" SO3)
endif()
target_link_libraries(SE3 PUBLIC SO3)
target_include_directories(SE3 PUBLIC "../SO3")

if (NOT TARGET BoxTree)
	add_subdirectory("../BoxTree" BoxTree)
endif()
target_link_libraries(SE3 PUBLIC BoxTree)
target_include_directories(SE3 PUBLIC "../BoxTree")

if (NOT TARGET WtFp)
	add_subdirectory("../WtFp" WtFp)
endif()
target_link_libraries(SE3 PUBLIC WtFp)
target_include_directories(SE3 PUBLIC "../WtFp")

