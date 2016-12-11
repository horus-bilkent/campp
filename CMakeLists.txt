cmake_minimum_required(VERSION 3.0.2)
project(campp)
find_package(OpenCV REQUIRED)

file(GLOB LIB "lib/*.cpp")

set( NAME_SRC
    src/debug.cpp
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

target_link_libraries(campp ${OpenCV_LIBS} )