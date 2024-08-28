add_library(IntroControl IntroControl.h IntroControl.cpp)
find_package(Eigen3 3.3 REQUIRED NO_MODULE)
target_link_libraries(IntroControl PUBLIC Eigen3::Eigen)

if (NOT TARGET FeatureSet)
	add_subdirectory("../FeatureSet" FeatureSet)
endif()
target_link_libraries(IntroControl PUBLIC FeatureSet)
target_include_directories(IntroControl PUBLIC "../FeatureSet")

if (NOT TARGET ReadSSSCommand)
	add_subdirectory("../ReadSSSCommand" ReadSSSCommand)
endif()
target_link_libraries(IntroControl PUBLIC ReadSSSCommand)
target_include_directories(IntroControl PUBLIC "../ReadSSSCommand")

