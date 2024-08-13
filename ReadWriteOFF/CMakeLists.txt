add_library(ReadWriteOFF ReadWriteOFF.h ReadWriteOFF.cpp)
find_package (Eigen3 3.3 REQUIRED NO_MODULE)
target_link_libraries(ReadWriteOFF PUBLIC Eigen3::Eigen)