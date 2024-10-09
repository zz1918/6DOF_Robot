add_library(DefMesh DefMesh.h DefMesh.cpp)
find_package (Eigen3 3.3 REQUIRED NO_MODULE)
target_link_libraries(DefMesh PUBLIC Eigen3::Eigen)


if (NOT TARGET SO3)
	add_subdirectory("../SO3" SO3)
endif()
target_link_libraries(DefMesh PUBLIC SO3)
target_include_directories(DefMesh PUBLIC "../SO3")

if (NOT TARGET ReadWriteOFF)
	add_subdirectory("../ReadWriteOFF" ReadWriteOFF)
endif()
target_link_libraries(DefMesh PUBLIC ReadWriteOFF)
target_include_directories(DefMesh PUBLIC "../ReadWriteOFF")

if (NOT TARGET ReadControl)
	add_subdirectory("../ReadControl" ReadControl)
endif()
target_link_libraries(DefMesh PUBLIC ReadControl)
target_include_directories(DefMesh PUBLIC "../ReadControl")