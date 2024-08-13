add_library(ReadSSSCommand ReadSSSCommand.h ReadSSSCommand.cpp)
find_package (Eigen3 3.3 REQUIRED NO_MODULE)
target_link_libraries(ReadSSSCommand PUBLIC Eigen3::Eigen)