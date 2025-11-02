add_library(Indicator Indicator.h Indicator.cpp)

if (NOT TARGET BinStr)
	add_subdirectory("../BinStr" BinStr)
endif()
target_link_libraries(SE3Box PUBLIC BinStr)
target_include_directories(SE3Box PUBLIC "../BinStr")

