set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_library(Solid Solid.h Solid.cpp)
find_package (Eigen3 3.3 REQUIRED NO_MODULE)
target_link_libraries(Solid PUBLIC Eigen3::Eigen)