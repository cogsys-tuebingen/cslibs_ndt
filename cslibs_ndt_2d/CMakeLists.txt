cmake_minimum_required(VERSION 2.8.3)
project(cslibs_ndt_2d CXX)

if(NOT ${CMAKE_BUILD_TYPE} STREQUAL Debug)
    set(TARGET_COMPILE_OPTIONS
        -Wall -Werror -Ofast -ffast-math -rdynamic
    )
    # enabling -march=native my result in problems on some machines,
    # my try -mtune=native

    message(STATUS "[${PROJECT_NAME}]: Compiling with optimization!")
endif()

find_package(catkin REQUIRED COMPONENTS
    cslibs_ndt
    cslibs_math_2d
    cslibs_math_ros
    cslibs_time
    cslibs_gridmaps
    cslibs_indexed_storage
    roscpp
    tf
)

find_package(Boost COMPONENTS filesystem)
find_package(yaml-cpp REQUIRED)

catkin_package(
  INCLUDE_DIRS
    include/
  CATKIN_DEPENDS
    cslibs_ndt
    cslibs_math_2d
    cslibs_time
    cslibs_gridmaps
)


set(TARGET_INCLUDE_DIRS
    include/
    ${catkin_INCLUDE_DIRS}
    ${Boost_INCLUDE_DIRS}
)

cslibs_ndt_2d_add_unit_test_gtest(${PROJECT_NAME}_test_serialization
    INCLUDE_DIRS
        ${TARGET_INCLUDE_DIRS}
    SOURCE_FILES
        test/serialization.cpp
    LINK_LIBRARIES
        ${Boost_LIBRARIES}
        ${YAML_CPP_LIBRARIES}
    COMPILE_OPTIONS
        ${TARGET_COMPILE_OPTIONS}
)

add_executable(${PROJECT_NAME}_map_loader
    src/ndt_map_loader.cpp
)

target_include_directories(${PROJECT_NAME}_map_loader
    PRIVATE
        ${TARGET_INCLUDE_DIRS}
)

target_compile_options(${PROJECT_NAME}_map_loader
    PRIVATE
        ${TARGET_COMPILE_OPTIONS}
)

target_link_libraries(${PROJECT_NAME}_map_loader
    PRIVATE
        ${catkin_LIBRARIES}
        ${YAML_CPP_LIBRARIES}
)


install(DIRECTORY include/${PROJECT_NAME}/
        DESTINATION ${CATKIN_PACKAGE_INCLUDE_DESTINATION})
