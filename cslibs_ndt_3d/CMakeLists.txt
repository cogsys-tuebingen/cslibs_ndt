cmake_minimum_required(VERSION 2.8.3)
project(cslibs_ndt_3d)

set(CSLIBS_NDT_USE_OMP False)

list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

find_package(catkin REQUIRED COMPONENTS
    cslibs_ndt
    cslibs_math_3d
    cslibs_math_ros
    cslibs_time
    cslibs_gridmaps
    roscpp
    tf
    message_generation
)
find_package(Boost COMPONENTS filesystem)
find_package(PCL REQUIRED COMPONENTS common filters segmentation io QUIET)
find_package(rviz QUIET)

add_message_files(
    FILES
    Distribution.msg
    DistributionArray.msg
)
generate_messages(
    DEPENDENCIES
    std_msgs
)

catkin_package(
  INCLUDE_DIRS include
  LIBRARIES
  CATKIN_DEPENDS cslibs_ndt cslibs_math_3d cslibs_time cslibs_gridmaps
)

include_directories(
    include
    ${catkin_INCLUDE_DIRS}
    ${Boost_INCLUDE_DIRS}
)

cslibs_ndt_3d_show_headers()

cslibs_ndt_add_unit_test_gtest(${PROJECT_NAME}_test_serialization
    SRCS test/serialization.cpp
)
target_link_libraries(${PROJECT_NAME}_test_serialization
    ${Boost_LIBRARIES}
    yaml-cpp
)

install(DIRECTORY include/${PROJECT_NAME}/
        DESTINATION ${CATKIN_PACKAGE_INCLUDE_DESTINATION})


if(rviz_FOUND)
    set(CMAKE_AUTOMOC ON)
    include_directories(${rviz_INCLUDE_DIRS})

    if(rviz_QT_VERSION VERSION_LESS "5")
        message(STATUS "Using Qt4 based on the rviz_QT_VERSION: ${rviz_QT_VERSION}")
        find_package(Qt4 ${rviz_QT_VERSION} EXACT REQUIRED QtCore QtGui)
        ## pull in all required include dirs, define QT_LIBRARIES, etc.
        include(${QT_USE_FILE})
    else()
        message(STATUS "Using Qt5 based on the rviz_QT_VERSION: ${rviz_QT_VERSION}")
        find_package(Qt5 ${rviz_QT_VERSION} EXACT REQUIRED Core Widgets)
        ## make target_link_libraries(${QT_LIBRARIES}) pull in all required dependencies
        set(QT_LIBRARIES Qt5::Widgets)
    endif()


    add_library(${PROJECT_NAME}_rviz
        src/rviz/ndt_visual.cpp
        src/rviz/ndt_display.cpp
        src/rviz/ndt_ellipsoid.cpp
        src/rviz/ndt_mesh.cpp
    )
    target_link_libraries(${PROJECT_NAME}_rviz
        ${QT_LIBRARIES}
        ${catkin_LIBRARIES}
    )
    add_dependencies(${PROJECT_NAME}_rviz ${${PROJECT_NAME}_EXPORTED_TARGETS})

    install(TARGETS ${PROJECT_NAME}_rviz
        ARCHIVE DESTINATION ${CATKIN_PACKAGE_LIB_DESTINATION}
        LIBRARY DESTINATION ${CATKIN_PACKAGE_LIB_DESTINATION}
        RUNTIME DESTINATION ${CATKIN_PACKAGE_BIN_DESTINATION}
    )
    install(FILES plugin_description.xml
            DESTINATION ${CATKIN_PACKAGE_SHARE_DESTINATION})
    install(DIRECTORY ogre_media
            DESTINATION ${CATKIN_PACKAGE_SHARE_DESTINATION})
    install(DIRECTORY icons/
            DESTINATION ${CATKIN_PACKAGE_SHARE_DESTINATION}/icons)
endif()


if(PCL_FOUND)
    include_directories(
        include
        ${catkin_INCLUDE_DIRS}
        ${PCL_INCLUDE_DIRS}
    )
    link_directories(${PCL_LIBRARY_DIRS})
    add_executable(${PROJECT_NAME}_map_loader
        src/nodes/ndt_map_loader.cpp
    )
    target_link_libraries(${PROJECT_NAME}_map_loader
        ${catkin_LIBRARIES}
        ${PCL_LIBRARIES}
        yaml-cpp
    )
    add_dependencies(${PROJECT_NAME}_map_loader ${${PROJECT_NAME}_EXPORTED_TARGETS})
endif()
