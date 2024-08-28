add_library(ViewerControl ViewerControl.h ViewerControl.cpp)

find_package(Eigen3 3.3 REQUIRED NO_MODULE)
target_link_libraries(ViewerControl PUBLIC Eigen3::Eigen)

list(PREPEND CMAKE_MODULE_PATH ../libigl)

# Libigl
include(libigl)

# Enable the target igl::glfw
igl_include(glfw)

# Link igl (and the glfw module) to your project
target_link_libraries(ViewerControl PUBLIC igl::glfw)

if (NOT TARGET FeatureSet)
	add_subdirectory("../FeatureSet" FeatureSet)
endif()
target_link_libraries(ViewerControl PUBLIC FeatureSet)
target_include_directories(ViewerControl PUBLIC "../FeatureSet")

if (NOT TARGET IntroControl)
	add_subdirectory("../IntroControl" IntroControl)
endif()
target_link_libraries(ViewerControl PUBLIC IntroControl)
target_include_directories(ViewerControl PUBLIC "../IntroControl")


