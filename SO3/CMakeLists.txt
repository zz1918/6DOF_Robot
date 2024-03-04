cmake_minimum_required(VERSION 3.24.0)
add_library(SO3 SO3.h SO3.cpp)
find_package (Eigen3 3.3 REQUIRED NO_MODULE)
target_link_libraries(SO3 PUBLIC Eigen3::Eigen)