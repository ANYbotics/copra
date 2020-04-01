/*
 * Copyright 2016-2019 CNRS-UM LIRMM, CNRS-AIST JRL
 * Copyright 2020 ANYbotics AG
 */

#pragma once

#include <Eigen/Core>

#include <gtest/gtest.h>

/*
 * Point-mass under gravity with external force input:
 *
 * p(k + 1) = p(k) + v(k) * T + (u / m + g) * T^2 / 2
 * v(k + 1) = v(k) + (u / m + g) * T
 *
 * The system is constrained by u(k) <= u_max and v(k) <= 0.
 *
 * The LTI system is then turned into an LTV one by adding time-variant offsets to the matrices A, B and vector c of its dynamics.
 *
 * \note The final point of the trajectory should be [val, 0] where val can be any value less than 0;
 */
class SmallTimeVariantSystem : public ::testing::Test {
 protected:
  SmallTimeVariantSystem()
      : T(0.005)
      , mass(5)
      , nbStep(10)
      , A(2, 2)
      , B(2, 1)
      , M(2, 2)
      , N(1, 1)
      , A_offset(2, 2)
      , B_offset(2, 1)
      , c(2)
      , uLower(1)
      , uUpper(1)
      , xLower(2)
      , xUpper(2)
      , x0(2)
      , xd(2)
      , ud(1)
      , wx(2)
      , wu(1)
      , c_offset(2)
  {
    // System
    A << 1, T, 0, 1;
    B << 0.5 * T * T / mass, T / mass;
    c << (-9.81 / 2.) * T * T, -9.81 * T;
    x0 << 0, -1.5;
    wx << 10, 10000;
    wu << 1e-4;
    A_offset << +0.05, -0.1 * T, 0., -0.2;
    B_offset << -0.05 * B(0), +0.02 * B(1);
    c_offset << 0.1 * c(0), -0.1 * c(1);

    // Cost
    M << 1, 0, 0, 1;
    N << 1;
    xd << 0, -1;
    ud << 2;

    // Control bound
    uLower.setConstant(-std::numeric_limits<double>::infinity());
    uUpper.setConstant(200); // force can't be superior to that

    // Trajectory bound
    xLower.setConstant(-std::numeric_limits<double>::infinity());
    xUpper(0) = std::numeric_limits<double>::infinity(); // no position limit
    xUpper(1) = 0; // velocity has to be negative

    // Expected solution, update if you change (A, B, c) or their offsets
    expectedTrajectory.resize(2 * (nbStep + 1));
    expectedTrajectory << 0, -1.5, -0.00768, -1.572, -0.0156781, -1.60926, -0.0239068, -1.61074, -0.0322835, -1.57793, -0.0407411, -1.51474,
        -0.049238, -1.42719, -0.0577663, -1.32296, -0.0663596, -1.21074, -0.0750978, -1.09973, -0.0841117, -0.999249;
    expectedControl.resize(nbStep);
    expectedControl << -22.952, -23.6299, -24.936, -26.9786, -29.9306, -34.0551, -39.7467, -47.5976, -58.5042, -73.8445;
  }

  Eigen::MatrixXd smallStateMatrix(int step) {
    assert(step < nbStep);
    return A + step * A_offset / (nbStep - 1);
  }

  Eigen::MatrixXd smallInputMatrix(int step) {
    assert(step < nbStep);
    return B + step * B_offset / (nbStep - 1);
  }

  Eigen::VectorXd smallDriftVector(int step) {
    assert(step < nbStep);
    return c + step * c_offset / (nbStep - 1);
  }

 protected:
  double T, mass;
  int nbStep;
  int stateDim = 2;
  int inputDim = 1;
  Eigen::MatrixXd A, B, M, N;
  Eigen::MatrixXd A_offset, B_offset;
  Eigen::VectorXd c, uLower, uUpper, xLower, xUpper, x0, xd, ud, wx, wu;
  Eigen::VectorXd c_offset;
  Eigen::VectorXd expectedTrajectory;
  Eigen::VectorXd expectedControl;
};
