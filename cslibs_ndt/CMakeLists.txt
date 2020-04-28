cmake_minimum_required(VERSION 2.8.3)
project(cslibs_ndt CXX)

find_package(catkin REQUIRED COMPONENTS
    cslibs_math
    cslibs_math_2d
    cslibs_math_3d
    cslibs_indexed_storage
    cslibs_utility
)

catkin_package(
  INCLUDE_DIRS include
  LIBRARIES
  CATKIN_DEPENDS
    cslibs_math
    cslibs_math_2d
    cslibs_math_3d
    cslibs_utility
    cslibs_indexed_storage
)

if(NOT ${CMAKE_BUILD_TYPE} STREQUAL Debug)
    set(TARGET_COMPILE_OPTIONS
        -Wall -Werror -Ofast -ffast-math
    )
    # link options  -rdynamic
    # enabling -march=native my result in problems on some machines,
    # my try -mtune=native

    message(STATUS "[${PROJECT_NAME}]: Compiling with optimization!")
endif()


set(TARGET_INCLUDE_DIRS
    include
    ${catkin_INCLUDE_DIRS}
)

cslibs_ndt_add_unit_test_gtest(${PROJECT_NAME}_test_generate_indices
    INCLUDE_DIRS
        ${TARGET_INCLUDE_DIRS}
    SOURCE_FILES
        test/test_generate_indices.cpp
    COMPILE_OPTIONS
        ${TARGET_COMPILE_OPTIONS}
)
cslibs_ndt_add_unit_test_gtest(${PROJECT_NAME}_test_merge_and
    INCLUDE_DIRS
        ${TARGET_INCLUDE_DIRS}
    SOURCE_FILES
        test/test_merge_and.cpp
    COMPILE_OPTIONS
        ${TARGET_COMPILE_OPTIONS}
)

install(DIRECTORY include/${PROJECT_NAME}/
        DESTINATION ${CATKIN_PACKAGE_INCLUDE_DESTINATION})
