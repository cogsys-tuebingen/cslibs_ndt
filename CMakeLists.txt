cmake_minimum_required(VERSION 2.8.3)
project(ndt)

## Find catkin macros and libraries
## if COMPONENTS list like find_package(catkin REQUIRED COMPONENTS xyz)
## is used, also find other catkin packages
find_package(catkin REQUIRED COMPONENTS
  roscpp
  sensor_msgs
  tf
)

## OpenMP
find_package(OpenMP REQUIRED)
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")

## Eigen
find_package(Eigen3 REQUIRED)

## OpenCV
find_package(OpenCV QUIET)

## C++11
include(CheckCXXCompilerFlag)
CHECK_CXX_COMPILER_FLAG("-std=c++11" COMPILER_SUPPORTS_CXX11)
CHECK_CXX_COMPILER_FLAG("-std=c++0x" COMPILER_SUPPORTS_CXX0X)
if(COMPILER_SUPPORTS_CXX11)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
elseif(COMPILER_SUPPORTS_CXX0X)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++0x")
else()
        message(STATUS "The compiler ${CMAKE_CXX_COMPILER} has no C++11 support. Please use a different C++ compiler.")
endif()

add_definitions(-W -Wall -Wno-unused-parameter -fno-strict-aliasing -Wno-unused-function -Wno-deprecated-register)
#add_definitions("-Ofast -march=native -ffast-math")
#add_definitions(-g -rdynamic)

catkin_package(
#  INCLUDE_DIRS include
#  LIBRARIES ndt
#  CATKIN_DEPENDS roscpp sensor_msgs tf
#  DEPENDS system_lib
)

include_directories(
    ${catkin_INCLUDE_DIRS}
    ${Eigen3_INCLUDE_DIRS}
    src
)

add_library(ndt
    src/data/pointcloud.hpp
    src/data/laserscan.hpp

    src/math/distribution.hpp

    src/ndt/rolling_distribution.hpp
    src/ndt/ndt_grid.cpp
    src/ndt/ndt_multi_grid.cpp
    src/ndt/types.hpp
)

target_link_libraries(ndt
    ${catkin_LIBRARIES}
)

add_executable(test_distribution
    src/tests/test_distribution.cpp
)

target_link_libraries(test_distribution
    ndt
)

#add_executable(ndt_scan_matching
#    src/nodes/scan_matching.cpp
#)

#target_link_libraries(ndt_scan_matching
#    ${catkin_LIBRARIES}
#    ndt
#)

if(${OpenCV_FOUND})
add_executable(demo_grid_visualization
    src/optional/demo_grid_visualization.cpp
    src/optional/visualize.hpp
)
target_link_libraries(demo_grid_visualization
    ndt
    ${OpenCV_LIBRARIES}
)


#add_executable(ndt_scan_visualization
#    src/nodes/visualization_node.cpp
#)

#target_link_libraries(ndt_scan_visualization
#    ${catkin_LIBRARIES}
#    ndt
#    ${OpenCV_LIBRARIES}
#)
endif()
