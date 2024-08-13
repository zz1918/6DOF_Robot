add_library(FeatureSet FeatureSet.h FeatureSet.cpp)
find_package (Eigen3 3.3 REQUIRED NO_MODULE)
target_link_libraries(FeatureSet PUBLIC Eigen3::Eigen)

if (NOT TARGET WtFp)
	add_subdirectory("../WtFp" WtFp)
endif()
target_link_libraries(FeatureSet PUBLIC WtFp)
target_include_directories(FeatureSet PUBLIC "../WtFp")

if (NOT TARGET ReadWriteOFF)
	add_subdirectory("../ReadWriteOFF" ReadWriteOFF)
endif()
target_link_libraries(FeatureSet PUBLIC ReadWriteOFF)
target_include_directories(FeatureSet PUBLIC "../ReadWriteOFF")

