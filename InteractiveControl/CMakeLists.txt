add_library(InteractiveControl InteractiveControl.h InteractiveControl.cpp)

find_package(Eigen3 3.3 REQUIRED NO_MODULE)
target_link_libraries(InteractiveControl PUBLIC Eigen3::Eigen)

if (NOT TARGET FeatureSet)
	add_subdirectory("../FeatureSet" FeatureSet)
endif()
target_link_libraries(InteractiveControl PUBLIC FeatureSet)
target_include_directories(InteractiveControl PUBLIC "../FeatureSet")

if (NOT TARGET ReadSSSCommand)
	add_subdirectory("../ReadSSSCommand" ReadSSSCommand)
endif()
target_link_libraries(InteractiveControl PUBLIC ReadSSSCommand)
target_include_directories(InteractiveControl PUBLIC "../ReadSSSCommand")


if (NOT TARGET IntroControl)
	add_subdirectory("../IntroControl" IntroControl)
endif()
target_link_libraries(InteractiveControl PUBLIC IntroControl)
target_include_directories(InteractiveControl PUBLIC "../IntroControl")

if (NOT TARGET RunControl)
	add_subdirectory("../RunControl" RunControl)
endif()
target_link_libraries(InteractiveControl PUBLIC RunControl)
target_include_directories(InteractiveControl PUBLIC "../RunControl")

if (NOT TARGET FileCheck)
	add_subdirectory("../FileCheck" FileCheck)
endif()
target_link_libraries(InteractiveControl PUBLIC FileCheck)
target_include_directories(InteractiveControl PUBLIC "../FileCheck")