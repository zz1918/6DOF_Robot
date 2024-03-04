add_library(Interval interval.h interval.cpp)
find_package (Eigen3 3.3 REQUIRED NO_MODULE)
target_link_libraries(Interval PUBLIC Eigen3::Eigen)