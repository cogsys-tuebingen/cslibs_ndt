cmake_minimum_required(VERSION 2.8.3)
project(cslibs_ndt_2d CXX)

list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

include(cmake/cslibs_ndt_2d_enable_c++11.cmake)
include(cmake/cslibs_ndt_2d_extras.cmake)
include(cmake/cslibs_ndt_2d_show_headers.cmake)
include(cmake/cslibs_ndt_2d_add_unit_test_gtest.cmake)

find_package(catkin REQUIRED COMPONENTS
    cslibs_ndt
    cslibs_math_2d
    cslibs_math_ros
    cslibs_time
    cslibs_gridmaps
    roscpp
    tf
)

find_package(Boost COMPONENTS filesystem)

catkin_package(
  INCLUDE_DIRS   include
  CATKIN_DEPENDS cslibs_ndt 
                 cslibs_math_2d
                 cslibs_time 
                 cslibs_gridmaps
  CFG_EXTRAS     cslibs_ndt_2d_enable_c++11.cmake
                 cslibs_ndt_2d_extras.cmake
                 cslibs_ndt_2d_show_headers.cmake
                 cslibs_ndt_2d_add_unit_test_gtest.cmake
)

include_directories(
    include
    ${catkin_INCLUDE_DIRS}
    ${Boost_INCLUDE_DIRS}
)

cslibs_ndt_2d_show_headers()

cslibs_ndt_2d_add_unit_test_gtest(${PROJECT_NAME}_test_serialization
    SRCS test/serialization.cpp
)
target_link_libraries(${PROJECT_NAME}_test_serialization
    ${Boost_LIBRARIES}
    yaml-cpp
)

install(DIRECTORY include/${PROJECT_NAME}/
        DESTINATION ${CATKIN_PACKAGE_INCLUDE_DESTINATION})


add_executable(${PROJECT_NAME}_map_loader
    src/ndt_map_loader.cpp
)
target_link_libraries(${PROJECT_NAME}_map_loader
    ${catkin_LIBRARIES}
    yaml-cpp
)
