cmake_minimum_required(VERSION 3.1)
project(ind-krig)

find_package(Threads REQUIRED)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -Wall")

set(eigen_SOURCE_DIR vendor/Eigen)
include_directories(${eigen_SOURCE_DIR})

include_directories(src)

set(kriging_SOURCE
	src/model/model.h
	src/model/model.cpp
	src/FileParser/fileparser.h
	src/FileParser/fileparser.cpp
	src/kriging/kriging.h
	src/kriging/kriging.cpp
	src/main.cpp
	)

add_executable(ind-krig ${kriging_SOURCE})
target_link_libraries(ind-krig ${CMAKE_THREAD_LIBS_INIT})
