add_library(SSS SSS.h SSS.cpp)
find_package (Eigen3 3.3 REQUIRED NO_MODULE)
target_link_libraries(SSS PUBLIC Eigen3::Eigen)

if (NOT TARGET SE3Box)
	add_subdirectory("../SE3Box" SE3Box)
endif()
target_link_libraries(SSS PUBLIC SE3Box)
target_include_directories(SSS PUBLIC "../SE3Box")

if (NOT TARGET WtFp)
	add_subdirectory("../WtFp" WtFp)
endif()
target_link_libraries(SSS PUBLIC WtFp)
target_include_directories(SSS PUBLIC "../WtFp")

if (NOT TARGET Graph)
	add_subdirectory("../Graph" Graph)
endif()
target_link_libraries(SSS PUBLIC Graph)
target_include_directories(SSS PUBLIC "../Graph")

if (NOT TARGET FeatureSet)
	add_subdirectory("../FeatureSet" FeatureSet)
endif()
target_link_libraries(SSS PUBLIC FeatureSet)
target_include_directories(SSS PUBLIC "../FeatureSet")

if (NOT TARGET ReadSSSCommand)
	add_subdirectory("../ReadSSSCommand" ReadSSSCommand)
endif()
target_link_libraries(SSS PUBLIC ReadSSSCommand)
target_include_directories(SSS PUBLIC "../ReadSSSCommand")


