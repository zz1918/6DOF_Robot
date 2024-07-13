add_library(Graph Graph.h Graph.cpp)

if (NOT TARGET Bimap)
	add_subdirectory("../Bimap" Bimap)
endif()
target_link_libraries(Graph PUBLIC Bimap)
target_include_directories(Graph PUBLIC "../Bimap")