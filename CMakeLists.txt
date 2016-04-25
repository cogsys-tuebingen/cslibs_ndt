cmake_minimum_required(VERSION 2.8.3)
project(ndt)

## Find catkin macros and libraries
## if COMPONENTS list like find_package(catkin REQUIRED COMPONENTS xyz)
## is used, also find other catkin packages
find_package(catkin REQUIRED COMPONENTS
  roscpp
  sensor_msgs
  tf
  pcl_conversions
  nav_msgs
)

## OpenMP
find_package(OpenMP REQUIRED)
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")

## Eigen
find_package(Eigen3 REQUIRED)

## OpenCV
find_package(OpenCV REQUIRED)

## PCL
find_package(PCL REQUIRED)

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
  INCLUDE_DIRS include
#  LIBRARIES ndt
#  CATKIN_DEPENDS roscpp sensor_msgs tf
#  DEPENDS system_lib
)

include_directories(
    ${catkin_INCLUDE_DIRS}
    ${Eigen3_INCLUDE_DIRS}
    ${PCL_INCLUDE_DIRS}
    ${OpenCV_INCLUDE_DIRS}
    include
)

add_library(ndt
    include/ndt/data/pointcloud.hpp
    include/ndt/data/laserscan.hpp

    include/ndt/math/distribution.hpp
    include/ndt/math/sgn.hpp

    include/ndt/conversion/convert.hpp

    include/ndt/matching/matcher.hpp
    include/ndt/matching/multi_grid_matcher_2D.hpp
    include/ndt/matching/single_grid_matcher_2D.hpp

    include/ndt/grid/grid.hpp
    include/ndt/grid/multi_grid.hpp
    include/ndt/grid/mask.hpp

    include/ndt/visualization/visualize.hpp
)
set_target_properties(ndt PROPERTIES LINKER_LANGUAGE CXX)

target_link_libraries(ndt
    ${catkin_LIBRARIES}
)

add_executable(ndt_test_distribution
    src/tests/test_distribution.cpp
)

target_link_libraries(ndt_test_distribution
    ndt
)

add_executable(ndt_test_mask
    src/tests/test_mask.cpp
)

add_executable(ndt_test_transformation
    src/tests/test_transformation.cpp
)

add_executable(ndt_test_solver
    src/tests/test_solver.cpp
)

add_executable(ndt_single_grid_scan_matching_node
    src/nodes/single_grid_scan_matching_node.cpp
)

target_link_libraries(ndt_single_grid_scan_matching_node
    ${catkin_LIBRARIES}
    ndt
    ${OpenCV_LIBRARIES}
)


add_executable(ndt_multi_grid_scan_matching_node
    src/nodes/multi_grid_scan_matching_node.cpp
)

target_link_libraries(ndt_multi_grid_scan_matching_node
    ${catkin_LIBRARIES}
    ndt
    ${OpenCV_LIBRARIES}
)

add_executable(ndt_hierarchic_scan_matching_node
    src/nodes/hierarchic_scan_matching_node.cpp
)

target_link_libraries(ndt_hierarchic_scan_matching_node
    ${catkin_LIBRARIES}
    ${PCL_LIBRIARIES}
    ndt
    ${OpenCV_LIBRARIES}
)

add_executable(ndt_demo_grid_visualization
    src/tests/optional/demo_grid_visualization.cpp
)
target_link_libraries(ndt_demo_grid_visualization
    ndt
    ${OpenCV_LIBRARIES}
)

add_executable(ndt_demo_matcher
    src/tests/optional/demo_matcher.cpp
)

target_link_libraries(ndt_demo_matcher
    ndt
    ${OpenCV_LIBRARIES}
)

add_executable(ndt_demo_matcher_pcl
    src/tests/optional/demo_matcher_pcl.cpp
)

target_link_libraries(ndt_demo_matcher_pcl
    ndt
    ${OpenCV_LIBRARIES}
    ${PCL_LIBRARIES}
)

add_executable(ndt_visualization_node
    src/nodes/visualization_node.cpp
)

target_link_libraries(ndt_visualization_node
    ${catkin_LIBRARIES}
    ndt
    ${OpenCV_LIBRARIES}
)
