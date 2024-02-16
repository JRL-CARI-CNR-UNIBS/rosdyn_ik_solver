cmake_minimum_required(VERSION 3.5)
project(rosdyn_ik_solver)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

find_package(Eigen3 3.3 REQUIRED NO_MODULE)
find_package(ik_solver_core REQUIRED)

set(ROS_VERSION $ENV{ROS_VERSION})
if(ROS_VERSION STREQUAL "")
  message(FATAL_ERROR "No supported version of ROS found!")
endif()
if(NOT((ROS_VERSION EQUAL 1) OR (ROS_VERSION EQUAL 2)))
  message(FATAL_ERROR "No supported version of ROS found!")
endif()

if(ROS_VERSION EQUAL 1)
  message("Building for ROS_VERSION == 1")
  find_package(catkin REQUIRED COMPONENTS
    roscpp
    rosdyn_core
    ik_solver
    pluginlib
  )

  catkin_package(
    INCLUDE_DIRS include
    LIBRARIES ${PROJECT_NAME}
    CATKIN_DEPENDS roscpp rosdyn_core ik_solver pluginlib
  #  DEPENDS system_lib
  )

  include_directories(
    include
    ${catkin_INCLUDE_DIRS}
    ${EIGEN3_INCLUDE_DIRS}
  )
elseif(ROS_VERSION EQUAL 2)
  message("Building for ROS_VERSION == 2")
  find_package(ament_cmake REQUIRED)
  find_package(rclcpp REQUIRED)
  find_package(ik_solver REQUIRED)
  find_package(rdyn_core REQUIRED)
  find_package(eigen_matrix_utils REQUIRED) # Until fix dependencies export on rdyn_core
  find_package(kinematics_filters REQUIRED) # Until fix dependencies export on rdyn_core
  find_package(pluginlib REQUIRED)
endif()

#############
## TARGETS ##
#############

add_library(${PROJECT_NAME} SHARED
  src/${PROJECT_NAME}/rosdyn_ik_solver.cpp
)

##########
## LINK ##
##########

if(ROS_VERSION EQUAL 1)
  add_dependencies(${PROJECT_NAME} ${${PROJECT_NAME}_EXPORTED_TARGETS} ${catkin_EXPORTED_TARGETS})

elseif(ROS_VERSION EQUAL 2)

  ament_target_dependencies(${PROJECT_NAME} PUBLIC
    rclcpp
    ik_solver
    pluginlib)
  target_link_libraries(${PROJECT_NAME} PUBLIC
    rdyn_core
    ${EIGEN3_LIBRARIES}
    ${eigen_matrix_utils_LIBRARIES}
    ${kinematics_filters_LIBRARIES})
  target_include_directories(${PROJECT_NAME} PUBLIC
    "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>"
    "$<INSTALL_INTERFACE:include/${PROJECT_NAME}>"
    )

  pluginlib_export_plugin_description_file(ik_solver plugin_ros2.xml)
  install(TARGETS ${PROJECT_NAME}
    EXPORT export_${PROJECT_NAME}
    LIBRARY DESTINATION lib
    ARCHIVE DESTINATION lib
    RUNTIME DESTINATION bin)
  ament_export_targets(export_${PROJECT_NAME})
  ament_export_dependencies(${PROJECT_NAME})
  ament_package()
endif()
