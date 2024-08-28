add_library(RunControl RunControl.h RunControl.cpp)

find_package(Eigen3 3.3 REQUIRED NO_MODULE)
target_link_libraries(RunControl PUBLIC Eigen3::Eigen)

if (NOT TARGET FeatureSet)
	add_subdirectory("../FeatureSet" FeatureSet)
endif()
target_link_libraries(RunControl PUBLIC FeatureSet)
target_include_directories(RunControl PUBLIC "../FeatureSet")

if (NOT TARGET ReadSSSCommand)
	add_subdirectory("../ReadSSSCommand" ReadSSSCommand)
endif()
target_link_libraries(RunControl PUBLIC ReadSSSCommand)
target_include_directories(RunControl PUBLIC "../ReadSSSCommand")
if (NOT TARGET IntroControl)
	add_subdirectory("../IntroControl" IntroControl)
endif()
target_link_libraries(RunControl PUBLIC IntroControl)
target_include_directories(RunControl PUBLIC "../IntroControl")

if (NOT TARGET ViewerControl)
	add_subdirectory("../ViewerControl" ViewerControl)
endif()
target_link_libraries(RunControl PUBLIC ViewerControl)
target_include_directories(RunControl PUBLIC "../ViewerControl")

if (NOT TARGET SSS)
	add_subdirectory("../SSS" SSS)
endif()
target_link_libraries(RunControl PUBLIC SSS)
target_include_directories(RunControl PUBLIC "../SSS")


