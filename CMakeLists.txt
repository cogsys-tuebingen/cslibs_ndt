cmake_minimum_required(VERSION 2.8.3)
project(cslibs_ndt)

list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

include(cmake/cslibs_ndt_enable_c++11.cmake)
include(cmake/cslibs_ndt_extras.cmake)
include(cmake/cslibs_ndt_openmp.cmake)
include(cmake/cslibs_ndt_show_headers.cmake)
include(cmake/cslibs_ndt_add_unit_test_gtest.cmake)

## Find catkin macros and libraries
## if COMPONENTS list like find_package(catkin REQUIRED COMPONENTS xyz)
## is used, also find other catkin packages
find_package(catkin REQUIRED COMPONENTS
  roscpp
  sensor_msgs
  tf
  pcl_conversions
  nav_msgs
  cslibs_indexed_storage
)

#set(CMAKE_BUILD_TYPE Debug)
set(BUILD_TESTS True)
set(BUILD_NODES False)
set(BUILD_TOOLS False)


## OpenCV ##########################################################################
find_package(OpenCV QUIET)
if(${OpenCV_FOUND})
set(OPENCV_DEPENDENT_SRCS
    include/ndt/visualization/grid.hpp
    include/ndt/visualization/multi_grid.hpp
   # include/ndt/visualization/kdtree.hpp
    include/ndt/visualization/points.hpp
)
endif()
## MRPT ############################################################################
find_package(MRPT QUIET)
if(${MRPT_FOUND})
    message("Building with MRPT support!")
    add_definitions("-DMRPT")
    set(MRPT_DEPENDENT_SRCS
        include/ndt/conversion/convert_mrpt.hpp
    )
endif()

## CATKIN ###########################################################################
catkin_package(
  INCLUDE_DIRS include
  LIBRARIES 
  CATKIN_DEPENDS cslibs_indexed_storage
#  DEPENDS system_lib
)

## INCLUDES #########################################################################
include_directories(
    ${catkin_INCLUDE_DIRS}
    ${EIGEN3_INCLUDE_DIR}
    ${OpenCV_INCLUDE_DIRS}
    ${MRPT_INCLUDE_DIRS}
    include
)

cslibs_ndt_show_headers()

#if(${BUILD_TESTS} AND ${OpenCV_FOUND})
### TESTS ############################################################################
#add_executable(ndt_test_distribution
#    src/tests/test_distribution.cpp
#)

#add_executable(ndt_test_mask
#    src/tests/test_mask.cpp
#)

#add_executable(ndt_test_transformation
#    src/tests/test_transformation.cpp
#)

#add_executable(ndt_test_solver
#    src/tests/test_solver.cpp
#)


### DEMOS ############################################################################
#add_executable(ndt_demo_grid_visualization
#    src/tests/demo_grid_visualization.cpp
#)
#target_link_libraries(ndt_demo_grid_visualization
#    ${OpenCV_LIBRARIES}
#    yaml-cpp
#)

#add_executable(ndt_demo_matcher
#    src/tests/demo_matcher.cpp
#)

#target_link_libraries(ndt_demo_matcher
#    ${OpenCV_LIBRARIES}
#    yaml-cpp
#)
#endif()

#if(${BUILD_NODES} AND ${OpenCV_FOUND})
### NODES ############################################################################
#add_executable(ndt_single_grid_scan_matching_node
#    src/nodes/single_grid_scan_matching_node.cpp
#)

#target_link_libraries(ndt_single_grid_scan_matching_node
#    ${catkin_LIBRARIES}
#    ${OpenCV_LIBRARIES}
#    yaml-cpp
#)

#add_executable(ndt_multi_grid_scan_matching_node
#    src/nodes/multi_grid_scan_matching_node.cpp
#)

#target_link_libraries(ndt_multi_grid_scan_matching_node
#    ${catkin_LIBRARIES}
#    ${OpenCV_LIBRARIES}
#    yaml-cpp
#)

#add_executable(ndt_hierarchic_scan_matching_node
#    src/nodes/hierarchic_scan_matching_node.cpp
#)

#target_link_libraries(ndt_hierarchic_scan_matching_node
#    ${catkin_LIBRARIES}
#    ${OpenCV_LIBRARIES}
#    yaml-cpp
#)

#add_executable(ndt_kdtree_scan_matching_node
#    src/nodes/kdtree_scan_matching_node.cpp
#)

#target_link_libraries(ndt_kdtree_scan_matching_node
#    ${catkin_LIBRARIES}
#    ${OpenCV_LIBRARIES}
#    yaml-cpp
#)

#add_executable(ndt_kdtree_cluster_node
#    src/nodes/kdtree_cluster_node.cpp
#)

#target_link_libraries(ndt_kdtree_cluster_node
#    ${catkin_LIBRARIES}
#    ${OpenCV_LIBRARIES}
#    yaml-cpp
#)
#endif()

#if(${BUILD_TOOLS} AND ${OpenCV_FOUND})
### TOOLS ############################################################################
#add_executable(ndt_scan_matcher
#    src/tools/scan_matcher.cpp
#)

#target_link_libraries(ndt_scan_matcher
#    ${OpenCV_LIBRARIES}
#    yaml-cpp
#)
#endif()
