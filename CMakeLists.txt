cmake_minimum_required(VERSION 3.0.2)
project(campp)
set(CMAKE_MODULE_PATH "/usr/local/lib/cmake/;${CMAKE_MODULE_PATH}")
set(CMAKE_CXX_FLAGS "-std=c++11")
option(Raspberry OFF)
find_package(OpenCV REQUIRED)
file(GLOB LIB "lib/*.cpp")


set( NAME_HEADERS
     include/center.hpp
     include/cv_algo.hpp
     include/geometry.hpp
     include/needle.hpp
     include/perimeter.hpp
     include/transform.hpp
)

if (Raspberry)
set (NAME_NEEDLE
	src/find_needle.cpp
)

else()
set( NAME_SRC
    src/debug.cpp
)

set (NAME_PREPROCESS
	src/client_preprocess.cpp
)

endif()


INCLUDE_DIRECTORIES( ${CMAKE_CURRENT_SOURCE_DIR}/include )
link_directories( ${CMAKE_BINARY_DIR}/bin)
set(EXECUTABLE_OUTPUT_PATH ${CMAKE_BINARY_DIR}/bin)

if (Raspberry)
add_executable(cneedle ${NAME_NEEDLE} ${LIB} ${NAME_HEADERS})
target_link_libraries(cneedle ${OpenCV_LIBS})
else()
# add_executable(campp ${NAME_SRC} ${LIB} ${NAME_HEADERS} )
add_executable(cpreprocess ${NAME_PREPROCESS} ${LIB} ${NAME_HEADERS})
# target_link_libraries(campp ${OpenCV_LIBS})
target_link_libraries(cpreprocess ${OpenCV_LIBS})
endif()