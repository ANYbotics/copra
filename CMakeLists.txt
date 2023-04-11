# BSD 2-Clause License
#
# Copyright (c) 2012-2019 CNRS-UM LIRMM, CNRS-AIST JRL
# Copyright (c) 2020 ANYbotics AG
# All rights reserved.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
#
# 1. Redistributions of source code must retain the above copyright notice, this
# list of conditions and the following disclaimer.
#
# 2. Redistributions in binary form must reproduce the above copyright notice,
# this list of conditions and the following disclaimer in the documentation
# and/or other materials provided with the distribution.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
# DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
# FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
# DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
# SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
# CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
# OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
# OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

cmake_minimum_required(VERSION 3.5.1)
project(copra)

set(CMAKE_CXX_STANDARD 14)
add_compile_options(-Wall -Wextra -Wpedantic)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

set(CATKIN_PACKAGE_DEPENDENCIES
  eigen-quadprog
  eigen-qld
  qp_oases
)

find_package(catkin REQUIRED COMPONENTS
  ${CATKIN_PACKAGE_DEPENDENCIES}
)

find_package(Eigen3 REQUIRED)

###################################
## catkin specific configuration ##
###################################

catkin_package(
  INCLUDE_DIRS
    include
    ${EIGEN3_INCLUDE_DIR} # Eigen does not set variables according to the CMake naming convention.
  LIBRARIES
    ${PROJECT_NAME}
  CATKIN_DEPENDS
    ${CATKIN_PACKAGE_DEPENDENCIES}
)

###########
## Build ##
###########

add_library(${PROJECT_NAME}
  # src/solvers/GUROBISolver.cpp  # commercial solver
  # src/solvers/LSSOLSolver.cpp  # also commercial
  # src/solvers/OSQPSolver.cpp  # we can add this interface by importing https://github.com/jrl-umi3218/eigen-osqp
  src/AutoSpan.cpp
  src/LMPC.cpp
  src/SolverInterface.cpp
  src/constraints.cpp
  src/costFunctions.cpp
  src/debugUtils.cpp
  src/solvers/QLDSolver.cpp
  src/solvers/qpOASESSolver.cpp
  src/solvers/QuadProgSolver.cpp
  src/solvers/utils.cpp
  src/systems/System.cpp
  src/systems/TimeInvariantSystem.cpp
  src/systems/TimeVariantSystem.cpp
)

add_dependencies(${PROJECT_NAME}
  ${catkin_EXPORTED_TARGETS}
)

target_compile_definitions(${PROJECT_NAME} PUBLIC
  # EIGEN_GUROBI_FOUND  # uncomment if adding Gurobi
  # EIGEN_LSSOL_FOUND  # uncomment if adding LSSOL
  # EIGEN_OSQP_FOUND  # uncomment if adding eigen-osqp
  EIGEN_QLD_FOUND
)

target_include_directories(${PROJECT_NAME} PRIVATE
  include
)

target_include_directories(${PROJECT_NAME} SYSTEM PUBLIC
  ${catkin_INCLUDE_DIRS}
  ${EIGEN3_INCLUDE_DIR}
)

target_link_libraries(${PROJECT_NAME}
  ${catkin_LIBRARIES}
)

#############
## Install ##
#############

install(
  TARGETS
    ${PROJECT_NAME}
  ARCHIVE DESTINATION ${CATKIN_PACKAGE_LIB_DESTINATION}
  LIBRARY DESTINATION ${CATKIN_PACKAGE_LIB_DESTINATION}
  RUNTIME DESTINATION ${CATKIN_PACKAGE_BIN_DESTINATION}
)

install(
  DIRECTORY
    include/${PROJECT_NAME}/
  DESTINATION ${CATKIN_PACKAGE_INCLUDE_DESTINATION}
)

##########
## Test ##
##########

if(CATKIN_ENABLE_TESTING)
  find_package(catkin REQUIRED
    COMPONENTS
      ${CATKIN_PACKAGE_DEPENDENCIES}
  )

  catkin_add_gtest(test_${PROJECT_NAME}
    test/TestSolvers.cpp
    test/TestTimeInvariantSystems.cpp
    test/TestTimeVariantSystems.cpp
    test/main.cpp
    test/tools.cpp
  )
endif()

if(TARGET test_${PROJECT_NAME})
  target_link_libraries(test_${PROJECT_NAME}
    ${PROJECT_NAME}
  )

  target_compile_definitions(test_${PROJECT_NAME} PUBLIC
    # EIGEN_GUROBI_FOUND  # uncomment if adding Gurobi
    # EIGEN_LSSOL_FOUND  # uncomment if adding LSSOL
    # EIGEN_OSQP_FOUND  # uncomment if adding eigen-osqp
    EIGEN_QLD_FOUND
  )

  target_include_directories(test_${PROJECT_NAME} PRIVATE
    include
    )

  target_include_directories(test_${PROJECT_NAME} SYSTEM PUBLIC
    ${catkin_INCLUDE_DIRS}
    ${EIGEN3_INCLUDE_DIR}
  )

  find_package(cmake_code_coverage QUIET)
  if(cmake_code_coverage_FOUND)
    add_gtest_coverage(
      TEST_BUILD_TARGETS test_${PROJECT_NAME}
    )
  endif(cmake_code_coverage_FOUND)
endif()
