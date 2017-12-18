cmake_minimum_required(VERSION 2.8.3)
project(cslibs_ndt_2d)

list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

find_package(catkin REQUIRED COMPONENTS
    cslibs_ndt
    cslibs_math_2d
    cslibs_math_ros
    cslibs_time
    roscpp
    tf
)

catkin_package(
  INCLUDE_DIRS include
  LIBRARIES
  CATKIN_DEPENDS cslibs_ndt cslibs_math_2d cslibs_time
)

include_directories(
    include
    ${catkin_INCLUDE_DIRS}
)

cslibs_ndt_2d_show_headers()

add_executable(${PROJECT_NAME}_matching_node
    src/nodes/matching_node.cpp
)

target_link_libraries(${PROJECT_NAME}_matching_node
    ${catkin_LIBRARIES}
)

install(DIRECTORY include/${PROJECT_NAME}/
        DESTINATION ${CATKIN_PACKAGE_INCLUDE_DESTINATION})
