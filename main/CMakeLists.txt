# This is the cmake file for SSS framework on the Delta robot
cmake_minimum_required(VERSION 3.24.0)
project(main)
add_executable(main main.cpp)

if(MINGW OR CYGWIN)
	add_definitions(-O3)
endif()

find_package(Eigen3 3.3 REQUIRED NO_MODULE)
target_link_libraries(main PUBLIC Eigen3::Eigen)

if (NOT TARGET FeatureSet)
	add_subdirectory("../FeatureSet" FeatureSet)
endif()
target_link_libraries(main PUBLIC FeatureSet)
target_include_directories(main PUBLIC "../FeatureSet")


if (NOT TARGET ReadSSSCommand)
	add_subdirectory("../ReadSSSCommand" ReadSSSCommand)
endif()
target_link_libraries(main PUBLIC ReadSSSCommand)
target_include_directories(main PUBLIC "../ReadSSSCommand")

if (NOT TARGET ReadControl)
	add_subdirectory("../ReadControl" ReadControl)
endif()
target_link_libraries(main PUBLIC ReadControl)
target_include_directories(main PUBLIC "../ReadControl")


if (NOT TARGET ViewerControl)
	add_subdirectory("../ViewerControl" ViewerControl)
endif()
target_link_libraries(main PUBLIC ViewerControl)
target_include_directories(main PUBLIC "../ViewerControl")

if (NOT TARGET RunControl)
	add_subdirectory("../RunControl" RunControl)
endif()
target_link_libraries(main PUBLIC RunControl)
target_include_directories(main PUBLIC "../RunControl")

if (NOT TARGET InteractiveControl)
	add_subdirectory("../InteractiveControl" InteractiveControl)
endif()
target_link_libraries(main PUBLIC InteractiveControl)
target_include_directories(main PUBLIC "../InteractiveControl")

if (NOT TARGET DefMesh)
	add_subdirectory("../DefMesh" DefMesh)
endif()
target_link_libraries(main PUBLIC DefMesh)
target_include_directories(main PUBLIC "../DefMesh")


