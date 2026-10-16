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

cmake_minimum_required(VERSION 3.16.3)
project(copra)

set(CMAKE_CXX_STANDARD 17)
add_compile_options(-Wall -Wextra -Wpedantic)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

set(PACKAGE_DEPENDENCIES
  qp_oases
  Eigen3
  )

find_package(ament_cmake REQUIRED)

foreach(PKG ${PACKAGE_DEPENDENCIES})
  find_package(${PKG} REQUIRED)
endforeach()

add_library(${PROJECT_NAME} SHARED
  src/AutoSpan.cpp
  src/LMPC.cpp
  src/SolverInterface.cpp
  src/constraints.cpp
  src/costFunctions.cpp
  src/debugUtils.cpp
  src/solvers/qpOASESSolver.cpp
  src/solvers/utils.cpp
  src/systems/System.cpp
  src/systems/TimeInvariantSystem.cpp
  src/systems/TimeVariantSystem.cpp
  )

target_include_directories(${PROJECT_NAME}
  PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include/${PROJECT_NAME}>
  )

target_link_libraries(${PROJECT_NAME}
  Eigen3::Eigen
  )
  
ament_target_dependencies(${PROJECT_NAME} ${PACKAGE_DEPENDENCIES})

#############
## Install ##
#############

install(DIRECTORY include/
  DESTINATION include/${PROJECT_NAME}
  )
  
install(TARGETS ${PROJECT_NAME}
  EXPORT export_${PROJECT_NAME}
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib
  RUNTIME DESTINATION bin
)

ament_export_targets(export_${PROJECT_NAME} HAS_LIBRARY_TARGET)
ament_export_dependencies(${PACKAGE_DEPENDENCIES})

##########
## Test ##
##########

if(BUILD_TESTING)

  find_package(ament_cmake_gtest REQUIRED)

  ament_add_gtest(test_${PROJECT_NAME}
    test/TestSolvers.cpp
    test/TestTimeInvariantSystems.cpp
    test/TestTimeVariantSystems.cpp
    test/main.cpp
    test/tools.cpp
    )

  target_link_libraries(test_${PROJECT_NAME}
    ${PROJECT_NAME}
  )



endif()

ament_package()
