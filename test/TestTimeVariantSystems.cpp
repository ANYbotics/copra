/*
 * Copyright 2020 ANYbotics AG
 */

#include <memory>
#include <vector>

#include <gtest/gtest.h>

#include "copra/LMPC.h"
#include "copra/constraints.h"
#include "copra/costFunctions.h"
#include "copra/systems/TimeInvariantSystem.h"
#include "copra/systems/TimeVariantSystem.h"

#include "time_invariant_systems.h"
#include "time_variant_systems.h"
#include "tools.h"

TEST_F(BoundedSystem, checkLTVMatrices) {  // NOLINT
  auto xCost = std::make_shared<copra::TargetCost>(M, xd);
  auto uCost = std::make_shared<copra::ControlCost>(N, ud);
  auto trajConstr = std::make_shared<copra::TrajectoryBoundConstraint>(xLower, xUpper);
  auto contConstr = std::make_shared<copra::ControlBoundConstraint>(uLower, uUpper);
  xCost->weights(wx);
  uCost->weights(wu);

  int stateDim = static_cast<int>(A.cols());
  int inputDim = static_cast<int>(B.cols());

  auto ltvSystem = std::make_shared<copra::TimeVariantSystem>();
  ltvSystem->reset(stateDim, inputDim, x0, nbStep);
  ltvSystem->setStateMatrix([this](int) { return A; });
  ltvSystem->setInputMatrix([this](int) { return B; });
  ltvSystem->setDriftVector([this](int) { return c; });

  for (int k = 0; k < nbStep; k++) {
    auto A_ltv = ltvSystem->getStateMatrix(k);
    auto B_ltv = ltvSystem->getInputMatrix(k);
    auto c_ltv = ltvSystem->getDriftVector(k);
    ASSERT_LE((A_ltv - A).norm(), 1e-8);
    ASSERT_LE((B_ltv - B).norm(), 1e-8);
    ASSERT_LE((c_ltv - c).norm(), 1e-8);
  }
}

TEST_F(BoundedSystem, checkLTVSolutionToLTISystem) {  // NOLINT
  auto xCost = std::make_shared<copra::TargetCost>(M, xd);
  auto uCost = std::make_shared<copra::ControlCost>(N, ud);
  xCost->weights(wx);
  uCost->weights(wu);

  auto ltiSystem = std::make_shared<copra::TimeInvariantSystem>();
  ltiSystem->reset(A, B, c, x0, nbStep);

  auto ltvSystem = std::make_shared<copra::TimeVariantSystem>();
  int stateDim = static_cast<int>(A.cols());
  int inputDim = static_cast<int>(B.cols());
  ltvSystem->reset(stateDim, inputDim, x0, nbStep);
  ltvSystem->setStateMatrix([this](int) { return A; });
  ltvSystem->setInputMatrix([this](int) { return B; });
  ltvSystem->setDriftVector([this](int) { return c; });

  auto ltiController = copra::LMPC(ltiSystem);
  auto ltiTrajConstr = std::make_shared<copra::TrajectoryBoundConstraint>(xLower, xUpper);
  auto ltiControlConstr = std::make_shared<copra::ControlBoundConstraint>(uLower, uUpper);
  ltiController.addCost(xCost);
  ltiController.addCost(uCost);
  ltiController.addConstraint(ltiTrajConstr);
  ltiController.addConstraint(ltiControlConstr);
  ltiController.selectQPSolver(copra::SolverFlag::DEFAULT);

  auto ltvController = copra::LMPC(ltvSystem);
  auto ltvTrajConstr = std::make_shared<copra::TrajectoryBoundConstraint>(xLower, xUpper);
  auto ltvControlConstr = std::make_shared<copra::ControlBoundConstraint>(uLower, uUpper);
  ltvController.addCost(xCost);
  ltvController.addCost(uCost);
  ltvController.addConstraint(ltvTrajConstr);
  ltvController.addConstraint(ltvControlConstr);
  ltvController.selectQPSolver(copra::SolverFlag::DEFAULT);

  ASSERT_TRUE(ltiController.solve());
  ASSERT_TRUE(ltvController.solve());
  ASSERT_LE((ltiController.trajectory() - ltvController.trajectory()).norm(), 1e-10);
  ASSERT_LE((ltiController.control() - ltvController.control()).norm(), 1e-10);
}

TEST_F(SmallTimeVariantSystem, solveSmallTimeVariantSystem) {  // NOLINT
  auto ps = std::make_shared<copra::TimeVariantSystem>();
  ps->reset(stateDim, inputDim, x0, nbStep);
  ps->setStateMatrix([this](int step) { return smallStateMatrix(step); });
  ps->setInputMatrix([this](int step) { return smallInputMatrix(step); });
  ps->setDriftVector([this](int step) { return smallDriftVector(step); });
  ps->update();

  auto controller = copra::LMPC(ps);
  auto xCost = std::make_shared<copra::TargetCost>(M, xd);
  auto uCost = std::make_shared<copra::ControlCost>(N, ud);
  auto trajConstr = std::make_shared<copra::TrajectoryBoundConstraint>(xLower, xUpper);
  auto contConstr = std::make_shared<copra::ControlBoundConstraint>(uLower, uUpper);
  xCost->weights(wx);
  uCost->weights(wu);

  controller.addCost(xCost);
  controller.addCost(uCost);
  controller.addConstraint(trajConstr);
  controller.addConstraint(contConstr);

  auto pcCheck = [&](const std::string& /* solverName */, copra::SolverFlag sFlag, std::unique_ptr<copra::SolverInterface>&& solver = nullptr) {
    if (solver) {
      controller.useSolver(std::move(solver));
    } else {
      controller.selectQPSolver(sFlag);
    }
    ASSERT_TRUE(controller.solve());

    Eigen::VectorXd trajectory = controller.trajectory();
    ASSERT_LE((trajectory - expectedTrajectory).norm(), 1e-4);

    auto trajLen = trajectory.rows() / 2;
    Eigen::VectorXd posTraj(trajLen);
    Eigen::VectorXd velTraj(trajLen);
    for (auto i = 0; i < trajLen; ++i) {
      posTraj(i) = trajectory(2 * i);
      velTraj(i) = trajectory(2 * i + 1);
    }
    Eigen::VectorXd control = controller.control();
    ASSERT_LE((control - expectedControl).norm(), 2e-4);

    // Check system dynamics on full solution
    for (int k = 0; 2 * (k + 1) + 1 < trajLen; k++) {
      Eigen::MatrixXd A_k = smallStateMatrix(k);
      Eigen::MatrixXd B_k = smallInputMatrix(k);
      Eigen::MatrixXd c_k = smallDriftVector(k);
      Eigen::VectorXd x_k = trajectory.segment<2>(2 * k);
      Eigen::VectorXd x_next = trajectory.segment<2>(2 * (k + 1));
      ASSERT_LE((x_next - A_k * x_k - B_k * control(k) - c_k).norm(), 1e-10);
    }

    // Check that terminal condition is fulfilled
    EXPECT_LE(std::abs(xd(1) - velTraj.tail(1)(0)), 0.001);

    // Check constraints
    ASSERT_LE(velTraj.maxCoeff(), xUpper(1) + 1e-6);
    ASSERT_LE(control.maxCoeff(), uUpper(0) + 1e-6);
  };

  for (auto s : tools::Solvers) {
    std::unique_ptr<copra::SolverInterface> solver(nullptr);
#ifdef EIGEN_LSSOL_FOUND
    if (s.second == copra::SolverFlag::LSSOL) {
            solver = copra::solverFactory(copra::SolverFlag::LSSOL);
            solver->SI_maxIter(200);
        }
#endif
    pcCheck(s.first, s.second, std::move(solver));
  }
}
