# Copra

[![License](https://img.shields.io/badge/License-BSD%202--Clause-green.svg)](https://opensource.org/licenses/BSD-2-Clause)
[![Documentation](https://img.shields.io/badge/doxygen-online-brightgreen?logo=read-the-docs&style=flat)](http://jrl-umi3218.github.io/copra/doxygen/HEAD/index.html)

Copra (**Co**ntrol & **pr**eview **a**lgorithms) is a C++ library for convex model predictive control of linear time-invariant (LTI) and linear time-variant (LTV) systems.

## Installation

### Debian Package

Assuming you have configured the ANYbotics PPA, you can install this package by:

    sudo apt install ros-${ROS_DISTRO}-copra

### Building from Source

#### Dependencies

- [CMake](cmake.org) >= 3.5
- [Eigen](http://eigen.tuxfamily.org/index.php?title=Main_Page) >= 3.2
- [Robot Operating System (ROS)](http://wiki.ros.org)
- [g++](https://gcc.gnu.org/)

#### Optional

* [GUROBI](http://www.gurobi.com/) >= 4.0: optional QP solver
* [OSQP](https://osqp.org/): optional QP solver
* [eigen-gurobi](https://github.com/vsamy/eigen-gurobi): bindings for GUROBI
* [eigen-osqp](https://github.com/jrl-umi3218/eigen-osqp.git): optional QP solver

#### Building with Catkin

To build from source, clone the latest version from this repository into your catkin workspace and compile the package using:

    cd ~/catkin_ws/src
    git clone git@github.com:ANYbotics/copra.git
    catkin build copra

### Unit Tests

Run unit tests by:

    catkin run_tests copra

## Documentation

Doxygen documentation is available [online](http://jrl-umi3218.github.io/copra/doxygen/HEAD/index.html)

You can also check out unit tests and the [bipedal locomotion example](https://vsamy.github.io/en/blog/copra-example-cpp).

## Changelog

This package was forked from [jrl-umi3218](https://github.com/jrl-umi3218/copra). The following changes were made:

- Extended to time-variant systems
- Catkin packaging for ROS
- Unit tests converted from doctest to Google Test
- Introduce simple trajectory and control tracking cost functions
- Introduce a method to reset all costs

The library was originally developed by [Vincent Samy](https://github.com/vsamy). It is licensed under the [BSD-2-Clause](https://opensource.org/licenses/BSD-2-Clause). Its default quadratic programming solver (eigen-quadprog) is licensed under the LGPLv2.
