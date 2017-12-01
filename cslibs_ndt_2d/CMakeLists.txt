cmake_minimum_required(VERSION 2.8.3)
project(cslibs_ndt_2d)

list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

find_package(catkin REQUIRED COMPONENTS
    cslibs_ndt
    cslibs_math_2d
)

catkin_package(
  INCLUDE_DIRS include
  LIBRARIES
  CATKIN_DEPENDS cslibs_ndt cslibs_math_2d
)

include_directories(
    include
    ${catkin_INCLUDE_DIRS}
)

cslibs_ndt_2d_show_headers()

install(DIRECTORY include/${PROJECT_NAME}/
        DESTINATION ${CATKIN_PACKAGE_INCLUDE_DESTINATION})
