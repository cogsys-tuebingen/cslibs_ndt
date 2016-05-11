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
  kdtree
)

set(CMAKE_BUILD_TYPE Debug)

## OpenMP ##########################################################################
find_package(OpenMP REQUIRED)
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")

## Eigen ###########################################################################
find_package(Eigen3 REQUIRED)

## OpenCV ##########################################################################
find_package(OpenCV REQUIRED)

## MRPT ############################################################################
find_package(MRPT QUIET)
if(${MRPT_FOUND})
    message("Building with MRPT support!")
    add_definitions("-DMRPT")
    set(MRPT_DEPENDENT_SRCS
        include/ndt/conversion/convert_mrpt.hpp
    )
endif()


## C++11 ###########################################################################
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

## FLAGS ############################################################################
add_definitions(-W -Wall -Wno-unused-parameter
                -fno-strict-aliasing -Wno-unused-function
                -Wno-deprecated-register)
#add_definitions(-mavx)
add_definitions(-msse2)
add_definitions(-msse4)
add_definitions(-Ofast -march=native -ffast-math)
#add_definitions(-g -rdynamic)

## CATKIN ###########################################################################
catkin_package(
  INCLUDE_DIRS include
  LIBRARIES ndt
  CATKIN_DEPENDS kdtree
#  DEPENDS system_lib
)

## INCLUDES #########################################################################
include_directories(
    ${catkin_INCLUDE_DIRS}
    ${Eigen3_INCLUDE_DIRS}
    ${OpenCV_INCLUDE_DIRS}
    ${MRPT_INCLUDE_DIRS}
    include
)

## NDT ##############################################################################
add_library(ndt
    include/ndt/data/pointcloud.hpp
    include/ndt/data/laserscan.hpp

    include/ndt/math/distribution.hpp
    include/ndt/math/sgn.hpp
    include/ndt/math/angle.hpp
    include/ndt/math/euclidian.hpp

    include/ndt/conversion/convert_ros.hpp

    include/ndt/matching/matcher.hpp
    include/ndt/matching/multi_grid_matcher_2D.hpp
#    include/ndt/matching/multi_matcher.hpp
    include/ndt/matching/single_grid_matcher_2D.hpp
    include/ndt/matching/kdtree_matcher_2D.hpp

    include/ndt/grid/grid.hpp
    include/ndt/grid/multi_grid.hpp
    include/ndt/grid/mask.hpp

    include/ndt/tree/kdtree.hpp

    include/ndt/visualization/grid.hpp
    include/ndt/visualization/multi_grid.hpp
    include/ndt/visualization/kdtree.hpp
    include/ndt/visualization/points.hpp

    ${MRPT_DEPENDENT_SRCS}
)

set_target_properties(ndt PROPERTIES LINKER_LANGUAGE CXX)

target_link_libraries(ndt
    ${catkin_LIBRARIES}
    yaml-cpp
)

## TESTS ############################################################################
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

## NODES ############################################################################
add_executable(ndt_single_grid_scan_matching_node
    src/nodes/single_grid_scan_matching_node.cpp
)

target_link_libraries(ndt_single_grid_scan_matching_node
    ${catkin_LIBRARIES}
    ${OpenCV_LIBRARIES}
    ndt
)

add_executable(ndt_multi_grid_scan_matching_node
    src/nodes/multi_grid_scan_matching_node.cpp
)

target_link_libraries(ndt_multi_grid_scan_matching_node
    ${catkin_LIBRARIES}
    ${OpenCV_LIBRARIES}
    ndt
)

add_executable(ndt_hierarchic_scan_matching_node
    src/nodes/hierarchic_scan_matching_node.cpp
)

target_link_libraries(ndt_hierarchic_scan_matching_node
    ${catkin_LIBRARIES}
    ${OpenCV_LIBRARIES}
    ndt
)

add_executable(ndt_kdtree_scan_matching_node
    src/nodes/kdtree_scan_matching_node.cpp
)

target_link_libraries(ndt_kdtree_scan_matching_node
    ${catkin_LIBRARIES}
    ${OpenCV_LIBRARIES}
    ndt
)

add_executable(ndt_kdtree_cluster_node
    src/nodes/kdtree_cluster_node.cpp
)

target_link_libraries(ndt_kdtree_cluster_node
    ${catkin_LIBRARIES}
    ${OpenCV_LIBRARIES}
    ndt
)

## TOOLS ############################################################################
#add_executable(ndt_scan_matcher
#    src/tools/scan_matcher.cpp
#)

#target_link_libraries(ndt_scan_matcher
#    ${OpenCV_LIBRARIES}
#    ndt
#)

## DEMOS ############################################################################
add_executable(ndt_demo_grid_visualization
    src/tests/optional/demo_grid_visualization.cpp
)
target_link_libraries(ndt_demo_grid_visualization
    ${OpenCV_LIBRARIES}
    ndt
)

add_executable(ndt_demo_matcher
    src/tests/optional/demo_matcher.cpp
)

target_link_libraries(ndt_demo_matcher
    ${OpenCV_LIBRARIES}
    ndt
)
