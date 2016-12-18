cmake_minimum_required(VERSION 3.0.2)
project(campp)
set(CMAKE_MODULE_PATH "/usr/local/lib/cmake/;${CMAKE_MODULE_PATH}")
SET(CMAKE_CXX_FLAGS "-std=c++11")
find_package(OpenCV REQUIRED)
find_package(raspicam REQUIRED)
file(GLOB LIB "lib/*.cpp")

set( NAME_SRC
    src/debug.cpp
)

set( NAME_DIFF
    src/diff_test.cpp
)

set (NAME_NEEDLE
	src/find_needle.cpp
)

set( NAME_HEADERS
     include/center.hpp
     include/cv_algo.hpp
     include/geometry.hpp
     include/needle.hpp
     include/perimeter.hpp
     include/transform.hpp
)

INCLUDE_DIRECTORIES( ${CMAKE_CURRENT_SOURCE_DIR}/include )
link_directories( ${CMAKE_BINARY_DIR}/bin)
set(EXECUTABLE_OUTPUT_PATH ${CMAKE_BINARY_DIR}/bin)
add_executable(campp ${NAME_SRC} ${LIB} ${NAME_HEADERS} )
add_executable(cdiff ${NAME_DIFF} ${LIB} ${NAME_HEADERS} )
add_executable(cneedle ${NAME_NEEDLE} ${LIB} ${NAME_HEADERS})
target_link_libraries(campp ${OpenCV_LIBS} )
target_link_libraries(cdiff ${OpenCV_LIBS} ${raspicam_CV_LIBS})
target_link_librarieS(cneedle ${OpenCV_LIBS})
