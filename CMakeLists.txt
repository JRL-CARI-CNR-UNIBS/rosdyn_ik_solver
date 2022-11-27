cmake_minimum_required(VERSION 3.0.2)
project(rosdyn_ik_solver)


if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

find_package(Eigen3 3.3 REQUIRED NO_MODULE)

find_package(catkin REQUIRED COMPONENTS
  eigen_conversions
  ik_solver_msgs
  roscpp
  rosdyn_core
  tf_conversions
  ik_solver
)

catkin_package(
  INCLUDE_DIRS include
  LIBRARIES rosdyn_ik_solver
  CATKIN_DEPENDS eigen_conversions ik_solver_msgs roscpp rosdyn_core tf_conversions ik_solver
#  DEPENDS system_lib
)

include_directories(
  include
  ${catkin_INCLUDE_DIRS}
  ${EIGEN3_INCLUDE_DIRS}
)

 add_library(${PROJECT_NAME}
   src/${PROJECT_NAME}/rosdyn_ik_solver.cpp
 )
 add_dependencies(${PROJECT_NAME} ${${PROJECT_NAME}_EXPORTED_TARGETS} ${catkin_EXPORTED_TARGETS})

# add_executable(${PROJECT_NAME}_node src/node.cpp)
# set_target_properties(${PROJECT_NAME}_node PROPERTIES OUTPUT_NAME rosdyn_ik PREFIX "")
# add_dependencies(${PROJECT_NAME}_node ${${PROJECT_NAME}_EXPORTED_TARGETS} ${catkin_EXPORTED_TARGETS})

# target_link_libraries(${PROJECT_NAME}_node
#   ${PROJECT_NAME}
#   ${catkin_LIBRARIES}
# )
