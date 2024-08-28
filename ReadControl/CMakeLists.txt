add_library(ReadControl ReadControl.h ReadControl.cpp)
find_package(Eigen3 3.3 REQUIRED NO_MODULE)
target_link_libraries(ReadControl PUBLIC Eigen3::Eigen)

if (NOT TARGET FeatureSet)
	add_subdirectory("../FeatureSet" FeatureSet)
endif()
target_link_libraries(ReadControl PUBLIC FeatureSet)
target_include_directories(ReadControl PUBLIC "../FeatureSet")

if (NOT TARGET ReadSSSCommand)
	add_subdirectory("../ReadSSSCommand" ReadSSSCommand)
endif()
target_link_libraries(ReadControl PUBLIC ReadSSSCommand)
target_include_directories(ReadControl PUBLIC "../ReadSSSCommand")

if (NOT TARGET JSON)
	add_subdirectory("../JSON" JSON)
endif()
target_link_libraries(ReadControl PUBLIC JSON)
target_include_directories(ReadControl PUBLIC "../JSON")

if (NOT TARGET FileCheck)
	add_subdirectory("../FileCheck" FileCheck)
endif()
target_link_libraries(ReadControl PUBLIC FileCheck)
target_include_directories(ReadControl PUBLIC "../FileCheck")


