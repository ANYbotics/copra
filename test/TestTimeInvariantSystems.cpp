/*
 * Copyright 2016-2019 CNRS-UM LIRMM, CNRS-AIST JRL
 * Copyright 2020 ANYbotics AG
 */

#include <algorithm>
#include <memory>
#include <numeric>
#include <vector>

#include <Eigen/Core>

#include <gtest/gtest.h>

#include "copra/LMPC.h"
#include "copra/constraints.h"
#include "copra/costFunctions.h"
#include "copra/solvers/all.h"
#include "copra/systems/TimeInvariantSystem.h"

#include "time_invariant_systems.h"
#include "tools.h"

/********************************************************************************************************
 *                               Check Bound constraint                                                 *
 ********************************************************************************************************/

TEST_F(BoundedSystem, mpcTargetCostWithBoundConstraints) {  // NOLINT
    auto ps = std::make_shared<copra::TimeInvariantSystem>();
    ps->reset(A, B, c, x0, nbStep);
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

        Eigen::VectorXd fullTraj = controller.trajectory();
        auto trajLen = fullTraj.rows() / 2;
        Eigen::VectorXd posTraj(trajLen);
        Eigen::VectorXd velTraj(trajLen);
        for (auto i = 0; i < trajLen; ++i) {
            posTraj(i) = fullTraj(2 * i);
            velTraj(i) = fullTraj(2 * i + 1);
        }
        Eigen::VectorXd control = controller.control();

        // Check result
        EXPECT_LE(std::abs(xd(1) - velTraj.tail(1)(0)), 0.001);

        // Check constraints
        ASSERT_LE(posTraj.maxCoeff(), x0(0));
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

TEST_F(BoundedSystem, mpcTrajectoryCostWithBoundConstraints) {  // NOLINT
    auto ps = std::make_shared<copra::TimeInvariantSystem>();
    ps->reset(A, B, c, x0, nbStep);
    auto controller = copra::LMPC(ps);
    auto xCost = std::make_shared<copra::TrajectoryCost>(M, xd);
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
        if (solver)
            controller.useSolver(std::move(solver));
        else
            controller.selectQPSolver(sFlag);

#ifdef EIGEN_OSQP_FOUND
    // Increase precision
    if (sFlag == copra::SolverFlag::OSQP) {
        auto& bs = static_cast<copra::OSQPSolver&>(controller.solver()).baseSolver();
        bs.scalingIter(0);
        bs.absConvergenceTol(1e-6);
        bs.relConvergenceTol(1e-6);
        bs.primalInfeasibilityTol(1e-7);
        bs.dualInfeasibilityTol(1e-7);
    }
#endif

        ASSERT_TRUE(controller.solve());

        Eigen::VectorXd fullTraj = controller.trajectory();
        auto trajLen = fullTraj.rows() / 2;
        Eigen::VectorXd posTraj(trajLen);
        Eigen::VectorXd velTraj(trajLen);
        for (auto i = 0; i < trajLen; ++i) {
            posTraj(i) = fullTraj(2 * i);
            velTraj(i) = fullTraj(2 * i + 1);
        }
        Eigen::VectorXd control = controller.control();

        // Check system dynamics on full solution
        for (int k = 0; 2 * (k + 1) + 1 < trajLen; k++) {
            Eigen::VectorXd x_k = fullTraj.segment<2>(2 * k);
            Eigen::VectorXd x_next = fullTraj.segment<2>(2 * (k + 1));
            ASSERT_LE((x_next - A * x_k - B * control(k) - c).norm(), 1e-10);
        }

        // Check result
        EXPECT_LE(std::abs(xd(1) - velTraj.tail(1)(0)), 0.001);

        // Check constraints
        ASSERT_LE(posTraj.maxCoeff(), x0(0));  // actually not an explicit constraint
        ASSERT_LE(velTraj.maxCoeff(), xUpper(1) + 1e-6);
        ASSERT_LE(control.maxCoeff(), uUpper(0) + 1e-6); // QuadProg allows to exceeds the constrain of a small amount.
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

TEST_F(BoundedSystem, mpcMixedCostWithBoundConstraints) {  // NOLINT
    auto ps = std::make_shared<copra::TimeInvariantSystem>();
    ps->reset(A, B, c, x0, nbStep);
    auto controller = copra::LMPC(ps);
    auto xCost = std::make_shared<copra::MixedCost>(M, Eigen::MatrixXd::Zero(2, 1), xd); // min(||X - Xt||^2)
    auto uCost = std::make_shared<copra::MixedCost>(Eigen::MatrixXd::Zero(1, 2), N, ud); // min(||U - Ut||^2)
    auto trajConstr = std::make_shared<copra::TrajectoryBoundConstraint>(xLower, xUpper);
    auto contConstr = std::make_shared<copra::ControlBoundConstraint>(uLower, uUpper);
    xCost->weights(wx);
    uCost->weights(wu);

    controller.addCost(xCost);
    controller.addCost(uCost);
    controller.addConstraint(trajConstr);
    controller.addConstraint(contConstr);

    ASSERT_TRUE(controller.solve());

    Eigen::VectorXd fullTraj = controller.trajectory();
    auto trajLen = fullTraj.rows() / 2;
    Eigen::VectorXd posTraj(trajLen);
    Eigen::VectorXd velTraj(trajLen);
    for (auto i = 0; i < trajLen; ++i) {
        posTraj(i) = fullTraj(2 * i);
        velTraj(i) = fullTraj(2 * i + 1);
    }
    Eigen::VectorXd control = controller.control();

    // Check system dynamics on full solution
    for (int k = 0; 2 * (k + 1) + 1 < trajLen; k++) {
        Eigen::VectorXd x_k = fullTraj.segment<2>(2 * k);
        Eigen::VectorXd x_next = fullTraj.segment<2>(2 * (k + 1));
        ASSERT_LE((x_next - A * x_k - B * control(k) - c).norm(), 1e-10);
    }

    // Check result
    EXPECT_LE(std::abs(xd(1) - velTraj.tail(3)(0)), 0.001); // Check X_{N-1} for mixed cost because X_N is not evaluated.

    // Check constrains
    ASSERT_LE(posTraj.maxCoeff(), x0(0));
    ASSERT_LE(velTraj.maxCoeff(), xUpper(1) + 1e-6);
    ASSERT_LE(control.maxCoeff(), uUpper(0) + 1e-6); // QuadProg allows to exceeds the constrain of a small amount.
}

/********************************************************************************************************
 *                                      Check Simple cost                                               *
 ********************************************************************************************************/

TEST_F(BoundedSystem, mpcSimpleTrajectoryCostEquivalenceTest) {  // NOLINT
  Eigen::VectorXd trajectory1 ;
  Eigen::VectorXd trajectory2;
  Eigen::VectorXd control1;
  Eigen::VectorXd control2;

  {
    auto ps = std::make_shared<copra::TimeInvariantSystem>();
    ps->reset(A, B, c, x0, nbStep);
    N.setIdentity();
    auto controller = copra::LMPC(ps);
    auto xCost = std::make_shared<copra::TrajectoryCost>(M, xd);
    auto uCost = std::make_shared<copra::ControlCost>(N, ud);
    auto trajConstr = std::make_shared<copra::TrajectoryBoundConstraint>(xLower, xUpper);
    auto contConstr = std::make_shared<copra::ControlBoundConstraint>(uLower, uUpper);
    xCost->weights(wx);
    uCost->weights(wu);

    controller.addCost(xCost);
    controller.addCost(uCost);
    controller.addConstraint(trajConstr);
    controller.addConstraint(contConstr);
    ASSERT_TRUE(controller.solve());
    trajectory1 = controller.trajectory();
    control1 = controller.trajectory();
  }

  {
    auto ps = std::make_shared<copra::TimeInvariantSystem>();
    ps->reset(A, B, c, x0, nbStep);
    N.setIdentity();
    auto controller = copra::LMPC(ps);
    auto xCost = std::make_shared<copra::TrajectoryCost>(M, xd);
    auto uCost = std::make_shared<copra::SimpleControlCost>(N, ud);
    auto trajConstr = std::make_shared<copra::TrajectoryBoundConstraint>(xLower, xUpper);
    auto contConstr = std::make_shared<copra::ControlBoundConstraint>(uLower, uUpper);
    xCost->weights(wx);
    uCost->weights(wu);

    controller.addCost(xCost);
    controller.addCost(uCost);
    controller.addConstraint(trajConstr);
    controller.addConstraint(contConstr);
    ASSERT_TRUE(controller.solve());
    trajectory2 = controller.trajectory();
    control2 = controller.trajectory();
  }

  ASSERT_TRUE(trajectory1.isApprox(trajectory2));
  ASSERT_TRUE(control1.isApprox(control2));

}

TEST_F(BoundedSystem, mpcSimpleControlCostEquivalenceTest) {  // NOLINT
  Eigen::VectorXd trajectory1 ;
  Eigen::VectorXd trajectory2;
  Eigen::VectorXd control1;
  Eigen::VectorXd control2;

  {
    auto ps = std::make_shared<copra::TimeInvariantSystem>();
    ps->reset(A, B, c, x0, nbStep);
    N.setIdentity();
    auto controller = copra::LMPC(ps);
    auto xCost = std::make_shared<copra::TrajectoryCost>(M, xd);
    auto uCost = std::make_shared<copra::ControlCost>(N, ud);
    auto trajConstr = std::make_shared<copra::TrajectoryBoundConstraint>(xLower, xUpper);
    auto contConstr = std::make_shared<copra::ControlBoundConstraint>(uLower, uUpper);
    xCost->weights(wx);
    uCost->weights(wu);

    controller.addCost(xCost);
    controller.addCost(uCost);
    controller.addConstraint(trajConstr);
    controller.addConstraint(contConstr);
    ASSERT_TRUE(controller.solve());
    trajectory1 = controller.trajectory();
    control1 = controller.trajectory();
  }

  {
    auto ps = std::make_shared<copra::TimeInvariantSystem>();
    ps->reset(A, B, c, x0, nbStep);
    N.setIdentity();
    auto controller = copra::LMPC(ps);
    auto xCost = std::make_shared<copra::SimpleTrajectoryCost>(M, xd);
    auto uCost = std::make_shared<copra::ControlCost>(N, ud);
    auto trajConstr = std::make_shared<copra::TrajectoryBoundConstraint>(xLower, xUpper);
    auto contConstr = std::make_shared<copra::ControlBoundConstraint>(uLower, uUpper);
    xCost->weights(wx);
    uCost->weights(wu);

    controller.addCost(xCost);
    controller.addCost(uCost);
    controller.addConstraint(trajConstr);
    controller.addConstraint(contConstr);
    ASSERT_TRUE(controller.solve());
    trajectory2 = controller.trajectory();
    control2 = controller.trajectory();
  }

  ASSERT_TRUE(trajectory1.isApprox(trajectory2));
  ASSERT_TRUE(control1.isApprox(control2));

}

/********************************************************************************************************
 *                            Check inequality constraint                                               *
 ********************************************************************************************************/

TEST_F(IneqSystem, mpcTargetCostWithInequalityConstraints) {  // NOLINT
    auto ps = std::make_shared<copra::TimeInvariantSystem>();
    ps->reset(A, B, c, x0, nbStep);
    auto controller = copra::LMPC(ps);
    auto xCost = std::make_shared<copra::TargetCost>(M, xd);
    auto uCost = std::make_shared<copra::ControlCost>(N, ud);
    auto trajConstr = std::make_shared<copra::TrajectoryConstraint>(E, p);
    auto contConstr = std::make_shared<copra::ControlConstraint>(G, h);
    xCost->weights(wx);
    uCost->weights(wu);

    controller.addCost(xCost);
    controller.addCost(uCost);
    controller.addConstraint(trajConstr);
    controller.addConstraint(contConstr);

    auto pcCheck = [&](const std::string& /* solverName */, copra::SolverFlag sFlag, std::unique_ptr<copra::SolverInterface>&& solver = nullptr) {
        if (solver)
            controller.useSolver(std::move(solver));
        else
            controller.selectQPSolver(sFlag);

#ifdef EIGEN_OSQP_FOUND
    // Increase precision
    if (sFlag == copra::SolverFlag::OSQP) {
        auto& bs = static_cast<copra::OSQPSolver&>(controller.solver()).baseSolver();
        bs.scalingIter(0);
        bs.absConvergenceTol(1e-6);
        bs.relConvergenceTol(1e-6);
        bs.primalInfeasibilityTol(1e-7);
        bs.dualInfeasibilityTol(1e-7);
        bs.relaxationParam(1);
        bs.maxIter(6000);
    }
#endif

        ASSERT_TRUE(controller.solve());

        Eigen::VectorXd fullTraj = controller.trajectory();
        auto trajLen = fullTraj.rows() / 2;
        Eigen::VectorXd posTraj(trajLen);
        Eigen::VectorXd velTraj(trajLen);
        for (auto i = 0; i < trajLen; ++i) {
            posTraj(i) = fullTraj(2 * i);
            velTraj(i) = fullTraj(2 * i + 1);
        }
        Eigen::VectorXd control = controller.control();

        // Check result
        EXPECT_LE(std::abs(xd(1) - velTraj.tail(1)(0)), 0.001);

        // Check constrains
        ASSERT_LE(posTraj.maxCoeff(), x0(0));
        ASSERT_LE(velTraj.maxCoeff(), p(0) + 1e-6);
        ASSERT_LE(control.maxCoeff(), h(0) + 1e-6);
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

TEST_F(IneqSystem, mpcTrajectoryCostWithInequalityConstraints) {  // NOLINT
    auto ps = std::make_shared<copra::TimeInvariantSystem>();
    ps->reset(A, B, c, x0, nbStep);
    auto controller = copra::LMPC(ps);
    auto xCost = std::make_shared<copra::TrajectoryCost>(M, xd);
    auto uCost = std::make_shared<copra::ControlCost>(N, ud);
    auto trajConstr = std::make_shared<copra::TrajectoryConstraint>(E, p);
    auto contConstr = std::make_shared<copra::ControlConstraint>(G, h);
    xCost->weights(wx);
    uCost->weights(wu);

    controller.addCost(xCost);
    controller.addCost(uCost);
    controller.addConstraint(trajConstr);
    controller.addConstraint(contConstr);

    auto pcCheck = [&](const std::string& /* solverName */, copra::SolverFlag sFlag, std::unique_ptr<copra::SolverInterface>&& solver = nullptr) {
        if (solver)
            controller.useSolver(std::move(solver));
        else
            controller.selectQPSolver(sFlag);

#ifdef EIGEN_OSQP_FOUND
    // Increase precision
    if (sFlag == copra::SolverFlag::OSQP) {
        auto& bs = static_cast<copra::OSQPSolver&>(controller.solver()).baseSolver();
        bs.scalingIter(0);
        bs.absConvergenceTol(1e-6);
        bs.relConvergenceTol(1e-6);
        bs.primalInfeasibilityTol(1e-7);
        bs.dualInfeasibilityTol(1e-7);
    }
#endif

        ASSERT_TRUE(controller.solve());

        Eigen::VectorXd fullTraj = controller.trajectory();
        auto trajLen = fullTraj.rows() / 2;
        Eigen::VectorXd posTraj(trajLen);
        Eigen::VectorXd velTraj(trajLen);
        for (auto i = 0; i < trajLen; ++i) {
            posTraj(i) = fullTraj(2 * i);
            velTraj(i) = fullTraj(2 * i + 1);
        }
        Eigen::VectorXd control = controller.control();

        // Check result
        EXPECT_LE(std::abs(xd(1) - velTraj.tail(1)(0)), 0.001);

        // Check constrains
        ASSERT_LE(posTraj.maxCoeff(), x0(0));
        ASSERT_LE(velTraj.maxCoeff(), p(0) + 1e-6);
        ASSERT_LE(control.maxCoeff(), h(0) + 1e-6); // QuadProg allows to exceeds the constrain of a small amount.
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

TEST_F(IneqSystem, mpcMixedCostWithInequalityConstraints) {  // NOLINT
    auto ps = std::make_shared<copra::TimeInvariantSystem>();
    ps->reset(A, B, c, x0, nbStep);
    auto controller = copra::LMPC(ps);
    auto xCost = std::make_shared<copra::MixedCost>(M, Eigen::MatrixXd::Zero(2, 1), xd);
    auto uCost = std::make_shared<copra::MixedCost>(Eigen::MatrixXd::Zero(1, 2), N, ud);
    auto trajConstr = std::make_shared<copra::TrajectoryConstraint>(E, p);
    auto contConstr = std::make_shared<copra::ControlConstraint>(G, h);
    xCost->weights(wx);
    uCost->weights(wu);

    controller.addCost(xCost);
    controller.addCost(uCost);
    controller.addConstraint(trajConstr);
    controller.addConstraint(contConstr);

    ASSERT_TRUE(controller.solve());

    Eigen::VectorXd fullTraj = controller.trajectory();
    auto trajLen = fullTraj.rows() / 2;
    Eigen::VectorXd posTraj(trajLen);
    Eigen::VectorXd velTraj(trajLen);
    for (auto i = 0; i < trajLen; ++i) {
        posTraj(i) = fullTraj(2 * i);
        velTraj(i) = fullTraj(2 * i + 1);
    }
    Eigen::VectorXd control = controller.control();

    // Check result
    EXPECT_LE(std::abs(xd(1) - velTraj.tail(3)(0)), 0.001); // Check X_{N-1} for mixed cost because X_N is not evaluated.

    // Check constrains
    ASSERT_LE(posTraj.maxCoeff(), x0(0));
    ASSERT_LE(velTraj.maxCoeff(), p(0) + 1e-6);
    ASSERT_LE(control.maxCoeff(), h(0) + 1e-6); // QuadProg allows to exceeds the constrain of a small amount.
}

/********************************************************************************************************
 *                               Check Mixed constraint                                                 *
 ********************************************************************************************************/

TEST_F(MixedSystem, mpcTargetCostWithMixedConstraints) {  // NOLINT
    auto ps = std::make_shared<copra::TimeInvariantSystem>();
    ps->reset(A, B, c, x0, nbStep);
    auto controller = copra::LMPC(ps);
    auto xCost = std::make_shared<copra::TargetCost>(M, xd);
    auto uCost = std::make_shared<copra::ControlCost>(N, ud);
    auto mixedConstr = std::make_shared<copra::MixedConstraint>(E, G, p);
    xCost->weights(wx);
    uCost->weights(wu);

    controller.addCost(xCost);
    controller.addCost(uCost);
    controller.addConstraint(mixedConstr);

    auto pcCheck = [&](const std::string& /* solverName */, copra::SolverFlag sFlag, std::unique_ptr<copra::SolverInterface>&& solver = nullptr) {
        if (solver)
            controller.useSolver(std::move(solver));
        else
            controller.selectQPSolver(sFlag);

        ASSERT_TRUE(controller.solve());

        Eigen::VectorXd fullTraj = controller.trajectory();
        auto trajLen = fullTraj.rows() / 2;
        Eigen::VectorXd posTraj(trajLen);
        Eigen::VectorXd velTraj(trajLen);
        for (auto i = 0; i < trajLen; ++i) {
            posTraj(i) = fullTraj(2 * i);
            velTraj(i) = fullTraj(2 * i + 1);
        }
        Eigen::VectorXd control = controller.control();

        // Check result
        EXPECT_LE(std::abs(xd(1) - velTraj.tail(1)(0)), 0.001);

        // Check constrains
        ASSERT_LE(posTraj.maxCoeff(), x0(0));
        for (int i = 0; i < nbStep; ++i) {
            auto res = E * fullTraj.segment(i * E.cols(), E.cols()) + G * control.segment(i * G.cols(), G.cols());
            if (!(res(0) <= p(0) + 1e-6)) {
                ASSERT_TRUE(false) << "Mixed constraint violated!";
            }
        }
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

TEST_F(MixedSystem, mpcTrajectoryCostWithMixedConstraints) {  // NOLINT
    auto ps = std::make_shared<copra::TimeInvariantSystem>();
    ps->reset(A, B, c, x0, nbStep);
    auto controller = copra::LMPC(ps);
    auto xCost = std::make_shared<copra::TrajectoryCost>(M, xd);
    auto uCost = std::make_shared<copra::ControlCost>(N, ud);
    auto mixedConstr = std::make_shared<copra::MixedConstraint>(E, G, p);
    xCost->weights(wx);
    uCost->weights(wu);

    controller.addCost(xCost);
    controller.addCost(uCost);
    controller.addConstraint(mixedConstr);

    auto pcCheck = [&](const std::string& /* solverName */, copra::SolverFlag sFlag, std::unique_ptr<copra::SolverInterface>&& solver = nullptr) {
        if (solver)
            controller.useSolver(std::move(solver));
        else
            controller.selectQPSolver(sFlag);

        ASSERT_TRUE(controller.solve());

        Eigen::VectorXd fullTraj = controller.trajectory();
        auto trajLen = fullTraj.rows() / 2;
        Eigen::VectorXd posTraj(trajLen);
        Eigen::VectorXd velTraj(trajLen);
        for (auto i = 0; i < trajLen; ++i) {
            posTraj(i) = fullTraj(2 * i);
            velTraj(i) = fullTraj(2 * i + 1);
        }
        Eigen::VectorXd control = controller.control();

        // Check result
        EXPECT_LE(std::abs(xd(1) - velTraj.tail(1)(0)), 0.001);

        // Check constrains
        ASSERT_LE(posTraj.maxCoeff(), x0(0));
        for (int i = 0; i < nbStep; ++i) {
            auto res = E * fullTraj.segment(i * E.cols(), E.cols()) + G * control.segment(i * G.cols(), G.cols());
            if (!(res(0) <= p(0) + 1e-6)) {
                ASSERT_TRUE(false) << "Mixed constraint violated!";
            }
        }
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

TEST_F(MixedSystem, mpcMixedCostWithMixedConstraints) {  // NOLINT
    auto ps = std::make_shared<copra::TimeInvariantSystem>();
    ps->reset(A, B, c, x0, nbStep);
    auto controller = copra::LMPC(ps);
    auto xCost = std::make_shared<copra::MixedCost>(M, Eigen::MatrixXd::Zero(2, 1), xd);
    auto uCost = std::make_shared<copra::MixedCost>(Eigen::MatrixXd::Zero(1, 2), N, ud);
    auto mixedConstr = std::make_shared<copra::MixedConstraint>(E, G, p);
    xCost->weights(wx);
    uCost->weights(wu);

    controller.addCost(xCost);
    controller.addCost(uCost);
    controller.addConstraint(mixedConstr);

    ASSERT_TRUE(controller.solve());

    Eigen::VectorXd fullTraj = controller.trajectory();
    auto trajLen = fullTraj.rows() / 2;
    Eigen::VectorXd posTraj(trajLen);
    Eigen::VectorXd velTraj(trajLen);
    for (auto i = 0; i < trajLen; ++i) {
        posTraj(i) = fullTraj(2 * i);
        velTraj(i) = fullTraj(2 * i + 1);
    }
    Eigen::VectorXd control = controller.control();

    // Check result
    EXPECT_LE(std::abs(xd(1) - velTraj.tail(3)(0)), 0.001); // Check X_{N-1} for mixed cost because X_N is not evaluated.

    // Check constrains
    ASSERT_LE(posTraj.maxCoeff(), x0(0));
    for (int i = 0; i < nbStep; ++i) {
        auto res = E * fullTraj.segment(i * E.cols(), E.cols()) + G * control.segment(i * G.cols(), G.cols());
        if (!(res(0) <= p(0) + 1e-6)) {
            ASSERT_TRUE(false) << "Mixed constraint violated!";
        }
    }
}

/********************************************************************************************************
 *                              Check Equality constraint                                               *
 ********************************************************************************************************/

TEST_F(EqSystem, mpcTargetCostWithEqualityConstraints) {  // NOLINT
    auto ps = std::make_shared<copra::TimeInvariantSystem>();
    ps->reset(A, B, c, x0, nbStep);
    auto controller = copra::LMPC(ps);
    auto xCost = std::make_shared<copra::TargetCost>(M, xd);
    auto uCost = std::make_shared<copra::ControlCost>(N, ud);
    auto trajConstr = std::make_shared<copra::TrajectoryConstraint>(E, p, false);
    xCost->weights(wx);
    uCost->weights(wu);

    controller.addCost(xCost);
    controller.addCost(uCost);
    controller.addConstraint(trajConstr);

    auto pcCheck = [&](const std::string& /* solverName */, copra::SolverFlag sFlag, std::unique_ptr<copra::SolverInterface>&& solver = nullptr) {
        if (solver)
            controller.useSolver(std::move(solver));
        else
            controller.selectQPSolver(sFlag);

#ifdef EIGEN_OSQP_FOUND
    // Increase precision
    if (sFlag == copra::SolverFlag::OSQP) {
        auto& bs = static_cast<copra::OSQPSolver&>(controller.solver()).baseSolver();
        bs.scalingIter(0);
        bs.absConvergenceTol(1e-7);
        bs.relConvergenceTol(1e-7);
        bs.primalInfeasibilityTol(1e-8);
        bs.dualInfeasibilityTol(1e-8);
        bs.relaxationParam(0.01);
    }
#endif

        ASSERT_TRUE(controller.solve());

        Eigen::VectorXd fullTraj = controller.trajectory();
        auto trajLen = fullTraj.rows() / 2;
        Eigen::VectorXd posTraj(trajLen);
        Eigen::VectorXd velTraj(trajLen);
        for (auto i = 0; i < trajLen; ++i) {
            posTraj(i) = fullTraj(2 * i);
            velTraj(i) = fullTraj(2 * i + 1);
        }
        Eigen::VectorXd control = controller.control();

        // Check result
        EXPECT_LE(std::abs(xd(1) - velTraj.tail(1)(0)), 0.001);

        // Check constrains
        ASSERT_LE(posTraj.maxCoeff(), x0(0) + 1e-6);
#ifdef EIGEN_OSQP_FOUND
        if (sFlag == copra::SolverFlag::OSQP)
            ASSERT_LE(velTraj.maxCoeff(), p(0) + 1e-4); // I could not get a better precision...
        else
#endif
        ASSERT_LE(velTraj.maxCoeff(), p(0) + 1e-6);
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

TEST_F(EqSystem, mpcTrajectoryCostWithEqualityConstraints) {  // NOLINT
    auto ps = std::make_shared<copra::TimeInvariantSystem>();
    ps->reset(A, B, c, x0, nbStep);
    auto controller = copra::LMPC(ps);
    auto xCost = std::make_shared<copra::TrajectoryCost>(M, xd);
    auto uCost = std::make_shared<copra::ControlCost>(N, ud);
    auto trajConstr = std::make_shared<copra::TrajectoryConstraint>(E, p, false);
    xCost->weights(wx);
    uCost->weights(wu);

    controller.addCost(xCost);
    controller.addCost(uCost);
    controller.addConstraint(trajConstr);

    auto pcCheck = [&](const std::string& /* solverName */, copra::SolverFlag sFlag, std::unique_ptr<copra::SolverInterface>&& solver = nullptr) {
        if (solver)
            controller.useSolver(std::move(solver));
        else
            controller.selectQPSolver(sFlag);

        ASSERT_TRUE(controller.solve());

        Eigen::VectorXd fullTraj = controller.trajectory();
        auto trajLen = fullTraj.rows() / 2;
        Eigen::VectorXd posTraj(trajLen);
        Eigen::VectorXd velTraj(trajLen);
        for (auto i = 0; i < trajLen; ++i) {
            posTraj(i) = fullTraj(2 * i);
            velTraj(i) = fullTraj(2 * i + 1);
        }
        Eigen::VectorXd control = controller.control();

        // Check result
        EXPECT_LE(std::abs(xd(1) - velTraj.tail(1)(0)), 0.001);

        // Check constrains
        ASSERT_LE(posTraj.maxCoeff(), x0(0) + 1e-6);
#ifdef EIGEN_OSQP_FOUND
        if (sFlag == copra::SolverFlag::OSQP)
            ASSERT_LE(velTraj.maxCoeff(), p(0) + 1e-4); // I could not get a better precision...
        else
#endif
        ASSERT_LE(velTraj.maxCoeff(), p(0) + 1e-6);
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

TEST_F(EqSystem, mpcMixedCostWithEqualityConstraints) {  // NOLINT
    auto ps = std::make_shared<copra::TimeInvariantSystem>();
    ps->reset(A, B, c, x0, nbStep);
    auto controller = copra::LMPC(ps);
    auto xCost = std::make_shared<copra::MixedCost>(M, Eigen::MatrixXd::Zero(2, 1), xd);
    auto uCost = std::make_shared<copra::MixedCost>(Eigen::MatrixXd::Zero(1, 2), N, ud);
    auto trajConstr = std::make_shared<copra::TrajectoryConstraint>(E, p, false);
    xCost->weights(wx);
    uCost->weights(wu);

    controller.addCost(xCost);
    controller.addCost(uCost);
    controller.addConstraint(trajConstr);

    ASSERT_TRUE(controller.solve());

    Eigen::VectorXd fullTraj = controller.trajectory();
    auto trajLen = fullTraj.rows() / 2;
    Eigen::VectorXd posTraj(trajLen);
    Eigen::VectorXd velTraj(trajLen);
    for (auto i = 0; i < trajLen; ++i) {
        posTraj(i) = fullTraj(2 * i);
        velTraj(i) = fullTraj(2 * i + 1);
    }
    Eigen::VectorXd control = controller.control();

    // Check result
    EXPECT_LE(std::abs(xd(1) - velTraj.tail(3)(0)), 0.001); // Check X_{N-1} for mixed cost because X_N is not evaluated.

    // Check constrains
    ASSERT_LE(posTraj.maxCoeff(), x0(0) + 1e-6);
    ASSERT_LE(velTraj.maxCoeff(), p(0) + 1e-6);
}

/********************************************************************************************************
 *                                   Check Autospan                                                     *
 ********************************************************************************************************/

TEST_F(BoundedSystem, checkAutospanAndWholeMatrixOnBoundConstraint) {  // NOLINT
    auto ps = std::make_shared<copra::TimeInvariantSystem>();
    ps->reset(A, B, c, x0, nbStep);
    auto controller = copra::LMPC(ps);
    int nbXStep = nbStep + 1;

    auto checkSpan = [&](const Eigen::VectorXd& xLower, const Eigen::VectorXd& xUpper, const Eigen::VectorXd& uLower, const Eigen::VectorXd& uUpper) {
        auto trajConstr = std::make_shared<copra::TrajectoryBoundConstraint>(xLower, xUpper);
        trajConstr->autoSpan();

        auto contConstr = std::make_shared<copra::ControlBoundConstraint>(uLower, uUpper);
        contConstr->autoSpan();

        ASSERT_NO_THROW(controller.addConstraint(trajConstr));
        ASSERT_NO_THROW(controller.addConstraint(contConstr));
    };

    auto fullxLower = tools::spanVector(xLower, nbXStep);
    auto fullxUpper = tools::spanVector(xUpper, nbXStep);
    auto fulluLower = tools::spanVector(uLower, nbStep);
    auto fulluUpper = tools::spanVector(uUpper, nbStep);

    checkSpan(xLower, xUpper, uLower, uUpper);
    checkSpan(fullxLower, xUpper, fulluLower, uUpper);
    checkSpan(xLower, fullxUpper, uLower, fulluUpper);
    checkSpan(fullxLower, fullxUpper, fulluLower, fulluUpper);
}

TEST_F(IneqSystem, checkAutospanAndWholeMatrixOnInequalityConstraint) {  // NOLINT
    auto ps = std::make_shared<copra::TimeInvariantSystem>();
    ps->reset(A, B, c, x0, nbStep);
    auto controller = copra::LMPC(ps);
    int nbXStep = nbStep + 1;

    auto checkSpan = [&](const Eigen::MatrixXd& E, const Eigen::VectorXd& p, const Eigen::MatrixXd& G, const Eigen::VectorXd& h) {
        auto trajConstr = std::make_shared<copra::TrajectoryConstraint>(E, p);
        trajConstr->autoSpan();

        auto contConstr = std::make_shared<copra::ControlConstraint>(G, h);
        contConstr->autoSpan();

        ASSERT_NO_THROW(controller.addConstraint(trajConstr));
        ASSERT_NO_THROW(controller.addConstraint(contConstr));
    };

    auto fullE = tools::spanMatrix(E, nbXStep);
    auto fullf = tools::spanVector(p, nbXStep);
    auto fullG = tools::spanMatrix(G, nbStep);
    auto fullh = tools::spanVector(h, nbStep);

    checkSpan(E, p, G, h);
    checkSpan(fullE, p, fullG, h);
    checkSpan(E, fullf, G, fullh);
    checkSpan(fullE, fullf, fullG, fullh);
}

TEST_F(MixedSystem, checkAutospanAndWholeMatrixOnMixedConstraint) {  // NOLINT
    auto ps = std::make_shared<copra::TimeInvariantSystem>();
    ps->reset(A, B, c, x0, nbStep);
    auto controller = copra::LMPC(ps);

    auto checkSpan = [&](const Eigen::MatrixXd& E, const Eigen::MatrixXd& G, const Eigen::VectorXd& p) {
        auto mixedConstr = std::make_shared<copra::MixedConstraint>(E, G, p);
        mixedConstr->autoSpan();

        ASSERT_NO_THROW(controller.addConstraint(mixedConstr));
    };

    auto fullE = tools::spanMatrix(E, nbStep, 1);
    auto fullG = tools::spanMatrix(G, nbStep);
    auto fullf = tools::spanVector(p, nbStep);

    checkSpan(E, G, p);
    checkSpan(fullE, G, p);
    checkSpan(fullE, fullG, p);
    checkSpan(fullE, G, fullf);
    checkSpan(E, fullG, p);
    checkSpan(fullE, fullG, p);
    checkSpan(E, fullG, fullf);
    checkSpan(E, G, fullf);
    checkSpan(fullE, G, fullf);
    checkSpan(E, fullG, fullf);
    checkSpan(fullE, fullG, fullf);
}

TEST_F(IneqSystem, checkAutospanAndWholeMatrixOnTrajectoryCost) {  // NOLINT
    auto ps = std::make_shared<copra::TimeInvariantSystem>();
    ps->reset(A, B, c, x0, nbStep);
    auto controller = copra::LMPC(ps);
    int nbXStep = nbStep + 1;

    auto checkSpan = [&](const Eigen::MatrixXd& M, const Eigen::VectorXd& p, const Eigen::VectorXd& weights) {
        auto cost = std::make_shared<copra::TrajectoryCost>(M, p);
        cost->weights(weights);
        cost->autoSpan();

        ASSERT_NO_THROW(controller.addCost(cost));
    };

    auto fullM = tools::spanMatrix(M, nbXStep);
    auto fullxd = tools::spanVector(xd, nbXStep);

    checkSpan(M, xd, wx);
    checkSpan(M, fullxd, wx);
    checkSpan(fullM, xd, wx);
    checkSpan(fullM, fullxd, wx);
}

TEST_F(IneqSystem, checkAutospanAndWholeMatrixOnControlCost) {  // NOLINT
    auto ps = std::make_shared<copra::TimeInvariantSystem>();
    ps->reset(A, B, c, x0, nbStep);
    auto controller = copra::LMPC(ps);

    auto checkSpan = [&](const Eigen::MatrixXd& M, const Eigen::VectorXd& p, const Eigen::VectorXd& weights) {
        auto cost = std::make_shared<copra::ControlCost>(M, p);
        cost->weights(weights);
        cost->autoSpan();

        ASSERT_NO_THROW(controller.addCost(cost));
    };

    auto fullN = tools::spanMatrix(N, nbStep);
    auto fullud = tools::spanVector(ud, nbStep);

    checkSpan(N, ud, wu);
    checkSpan(N, fullud, wu);
    checkSpan(fullN, ud, wu);
    checkSpan(fullN, fullud, wu);
}

TEST_F(IneqSystem, checkAutospanAndWholeMatrixOnMixedCost) {  // NOLINT
    auto ps = std::make_shared<copra::TimeInvariantSystem>();
    ps->reset(A, B, c, x0, nbStep);
    auto controller = copra::LMPC(ps);

    auto checkSpan = [&](const Eigen::MatrixXd& M, const Eigen::MatrixXd& N, const Eigen::VectorXd& p, const Eigen::VectorXd& weights) {
        auto cost = std::make_shared<copra::MixedCost>(M, N, p);
        cost->weights(weights);
        cost->autoSpan();

        ASSERT_NO_THROW(controller.addCost(cost));
    };

    auto MVec = std::vector<Eigen::MatrixXd>();
    MVec.push_back(M);
    MVec.push_back(tools::spanMatrix(M, nbStep, 1));
    auto nnVec = std::vector<Eigen::MatrixXd>();
    nnVec.push_back(Eigen::MatrixXd::Ones(2, 1));
    nnVec.push_back(tools::spanMatrix(Eigen::MatrixXd::Ones(2, 1), nbStep));
    auto xdVec = std::vector<Eigen::VectorXd>();
    xdVec.push_back(xd);
    xdVec.push_back(tools::spanVector(xd, nbStep));

    for (auto& i : MVec)
        for (auto& j : nnVec)
            for (auto& k : xdVec)
                checkSpan(i, j, k, wx);
}

/********************************************************************************************************
 *                                Check Error Messages                                                  *
 ********************************************************************************************************/

TEST_F(IneqSystem, errorHandlerForSystem) {  // NOLINT
    auto ps = std::make_shared<copra::TimeInvariantSystem>();
    ASSERT_THROW(ps->reset(Eigen::MatrixXd::Ones(5, 2), B, c, x0, nbStep), std::domain_error);
    ASSERT_THROW(ps->reset(Eigen::MatrixXd::Ones(2, 5), B, c, x0, nbStep), std::domain_error);
    ASSERT_THROW(ps->reset(A, Eigen::MatrixXd::Ones(5, 1), c, x0, nbStep), std::domain_error);
    ASSERT_THROW(ps->reset(A, B, Eigen::VectorXd::Ones(5), x0, nbStep), std::domain_error);
    ASSERT_THROW(ps->reset(A, B, c, x0, -1), std::domain_error);
}

TEST_F(IneqSystem, errorHandlerForWeigths) {  // NOLINT
    auto ps = std::make_shared<copra::TimeInvariantSystem>();
    ps->reset(A, B, c, x0, nbStep);
    auto controller = copra::LMPC(ps);
    auto cost = std::make_shared<copra::TrajectoryCost>(M, xd);

    ASSERT_NO_THROW(cost->weight(2));
    ASSERT_THROW(cost->weights(Eigen::VectorXd::Ones(5)), std::domain_error);
    ASSERT_NO_THROW(cost->weights(wx));
    ASSERT_NO_THROW(controller.addCost(cost));
    ASSERT_NO_THROW(cost->weights(Eigen::VectorXd::Ones(2)));
}

TEST_F(IneqSystem, errorHandlerForTrajectoryCost) {  // NOLINT
    auto ps = std::make_shared<copra::TimeInvariantSystem>(A, B, c, x0, nbStep);
    auto controller = copra::LMPC(ps);

    auto badCost1 = std::make_shared<copra::TrajectoryCost>(Eigen::MatrixXd::Identity(5, 5), Eigen::VectorXd::Ones(2));
    ASSERT_THROW(controller.addCost(badCost1), std::domain_error);
    auto badCost2 = std::make_shared<copra::TrajectoryCost>(Eigen::MatrixXd::Identity(5, 5), Eigen::VectorXd::Ones(5));
    ASSERT_THROW(controller.addCost(badCost2), std::domain_error);
}

TEST_F(IneqSystem, errorHandlerForTargetCost) {  // NOLINT
    auto ps = std::make_shared<copra::TimeInvariantSystem>(A, B, c, x0, nbStep);
    auto controller = copra::LMPC(ps);

    auto badCost1 = std::make_shared<copra::TargetCost>(Eigen::MatrixXd::Identity(5, 5), Eigen::VectorXd::Ones(2));
    ASSERT_THROW(controller.addCost(badCost1), std::domain_error);
    auto badCost2 = std::make_shared<copra::TargetCost>(Eigen::MatrixXd::Identity(5, 5), Eigen::VectorXd::Ones(5));
    ASSERT_THROW(controller.addCost(badCost2), std::domain_error);
}

TEST_F(IneqSystem, errorHandlerForControlCost) {  // NOLINT
    auto ps = std::make_shared<copra::TimeInvariantSystem>(A, B, c, x0, nbStep);
    auto controller = copra::LMPC(ps);

    auto badCost1 = std::make_shared<copra::ControlCost>(Eigen::MatrixXd::Identity(5, 5), Eigen::VectorXd::Ones(2));
    ASSERT_THROW(controller.addCost(badCost1), std::domain_error);
    auto badCost2 = std::make_shared<copra::ControlCost>(Eigen::MatrixXd::Identity(5, 5), Eigen::VectorXd::Ones(5));
    ASSERT_THROW(controller.addCost(badCost2), std::domain_error);
}

TEST_F(IneqSystem, errorHandlerForMixedCost) {  // NOLINT
    auto ps = std::make_shared<copra::TimeInvariantSystem>(A, B, c, x0, nbStep);
    auto controller = copra::LMPC(ps);

    auto badCost1 = std::make_shared<copra::MixedCost>(Eigen::MatrixXd::Identity(5, 5), Eigen::MatrixXd::Identity(2, 1), Eigen::VectorXd::Ones(2));
    ASSERT_THROW(controller.addCost(badCost1), std::domain_error);
    auto badCost2 = std::make_shared<copra::MixedCost>(Eigen::MatrixXd::Identity(2, 1), Eigen::MatrixXd::Identity(5, 5), Eigen::VectorXd::Ones(2));
    ASSERT_THROW(controller.addCost(badCost2), std::domain_error);
    auto badCost3 = std::make_shared<copra::MixedCost>(Eigen::MatrixXd::Identity(5, 5), Eigen::MatrixXd::Identity(5, 5), Eigen::VectorXd::Ones(5));
    ASSERT_THROW(controller.addCost(badCost3), std::domain_error);
}

TEST_F(IneqSystem, errorHandlerForTrajectoryConstraint) {  // NOLINT
    auto ps = std::make_shared<copra::TimeInvariantSystem>();
    ps->reset(A, B, c, x0, nbStep);
    auto controller = copra::LMPC(ps);

    auto badConstr = std::make_shared<copra::TrajectoryConstraint>(Eigen::MatrixXd::Identity(5, 5), Eigen::VectorXd::Ones(2));
    ASSERT_THROW(controller.addConstraint(badConstr), std::domain_error);
    auto trajConstr = std::make_shared<copra::TrajectoryConstraint>(Eigen::MatrixXd::Identity(5, 5), Eigen::VectorXd::Ones(5));
    ASSERT_THROW(controller.addConstraint(trajConstr), std::domain_error);
}

TEST_F(IneqSystem, errorHandlerForControlConstraint) {  // NOLINT
    auto ps = std::make_shared<copra::TimeInvariantSystem>();
    ps->reset(A, B, c, x0, nbStep);
    auto controller = copra::LMPC(ps);

    auto badConstr1 = std::make_shared<copra::ControlConstraint>(Eigen::MatrixXd::Identity(5, 5), Eigen::VectorXd::Ones(2));
    ASSERT_THROW(controller.addConstraint(badConstr1), std::domain_error);
    auto badConstr2 = std::make_shared<copra::ControlConstraint>(Eigen::MatrixXd::Identity(5, 5), Eigen::VectorXd::Ones(5));
    ASSERT_THROW(controller.addConstraint(badConstr2), std::domain_error);

    auto goodConstr = std::make_shared<copra::ControlConstraint>(G, h);
    controller.addConstraint(goodConstr);
    ASSERT_THROW(controller.addConstraint(goodConstr), std::runtime_error);
}

TEST_F(IneqSystem, errorHandlerForMixedConstraint) {  // NOLINT
    auto ps = std::make_shared<copra::TimeInvariantSystem>();
    ps->reset(A, B, c, x0, nbStep);
    auto controller = copra::LMPC(ps);

    auto badConstr1 = std::make_shared<copra::MixedConstraint>(Eigen::MatrixXd::Identity(5, 5), Eigen::MatrixXd::Identity(2, 1), Eigen::VectorXd::Ones(2));
    ASSERT_THROW(controller.addConstraint(badConstr1), std::domain_error);
    auto badConstr2 = std::make_shared<copra::MixedConstraint>(Eigen::MatrixXd::Identity(2, 1), Eigen::MatrixXd::Identity(5, 5), Eigen::VectorXd::Ones(2));
    ASSERT_THROW(controller.addConstraint(badConstr2), std::domain_error);
    auto badConstr3 = std::make_shared<copra::MixedConstraint>(Eigen::MatrixXd::Identity(5, 5), Eigen::MatrixXd::Identity(5, 5), Eigen::VectorXd::Ones(5));
    ASSERT_THROW(controller.addConstraint(badConstr3), std::domain_error);
}

TEST_F(BoundedSystem, errorHandlerForTrajectoryBoundConstraint) {  // NOLINT
    auto ps = std::make_shared<copra::TimeInvariantSystem>();
    ps->reset(A, B, c, x0, nbStep);
    auto controller = copra::LMPC(ps);

    auto badConstr = std::make_shared<copra::TrajectoryBoundConstraint>(Eigen::VectorXd::Ones(3), Eigen::VectorXd::Ones(2));
    ASSERT_THROW(controller.addConstraint(badConstr), std::domain_error);
    auto tbConstr = std::make_shared<copra::TrajectoryBoundConstraint>(Eigen::VectorXd::Ones(3), Eigen::VectorXd::Ones(3));
    ASSERT_THROW(controller.addConstraint(tbConstr), std::domain_error);
}

TEST_F(BoundedSystem, errorHandlerForControlBoundConstraint) {  // NOLINT
    auto ps = std::make_shared<copra::TimeInvariantSystem>();
    ps->reset(A, B, c, x0, nbStep);
    auto controller = copra::LMPC(ps);

    auto badConstr1 = std::make_shared<copra::ControlBoundConstraint>(Eigen::VectorXd::Ones(3), Eigen::VectorXd::Ones(2));
    ASSERT_THROW(controller.addConstraint(badConstr1), std::domain_error);
    auto badConstr2 = std::make_shared<copra::ControlBoundConstraint>(Eigen::VectorXd::Ones(3), Eigen::VectorXd::Ones(3));
    ASSERT_THROW(controller.addConstraint(badConstr2), std::domain_error);

    auto goodConstr = std::make_shared<copra::ControlBoundConstraint>(uLower, uUpper);
    controller.addConstraint(goodConstr);
    ASSERT_THROW(controller.addConstraint(goodConstr), std::runtime_error);
}

/********************************************************************************************************
 *                               Check remove functions                                                 *
 ********************************************************************************************************/

TEST_F(IneqSystem, removeCostAndConstraint) { // NOLINT
    auto ps = std::make_shared<copra::TimeInvariantSystem>(A, B, c, x0, nbStep);
    auto controller = copra::LMPC(ps);
    {
        auto xCost = std::make_shared<copra::TargetCost>(M, xd);
        auto uCost = std::make_shared<copra::ControlCost>(N, ud);
        auto trajConstr = std::make_shared<copra::TrajectoryConstraint>(E, p);
        auto contConstr = std::make_shared<copra::ControlConstraint>(G, h);

        controller.addCost(xCost);
        controller.addCost(uCost);
        controller.addConstraint(trajConstr);
        controller.addConstraint(contConstr);

        controller.removeCost(xCost);
        controller.removeCost(uCost);
        controller.removeConstraint(trajConstr);
        controller.removeConstraint(contConstr);
    }
    ASSERT_TRUE(controller.solve());
}

TEST_F(SmallSystem, solveSmallSystem) {  // NOLINT
  auto ps = std::make_shared<copra::TimeInvariantSystem>();
  ps->reset(A, B, c, x0, nbStep);
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
      Eigen::VectorXd x_k = trajectory.segment<2>(2 * k);
      Eigen::VectorXd x_next = trajectory.segment<2>(2 * (k + 1));
      ASSERT_LE((x_next - A * x_k - B * control(k) - c).norm(), 1e-10);
    }

    // Check that terminal condition is fulfilled
    EXPECT_LE(std::abs(xd(1) - velTraj.tail(1)(0)), 0.001);

    // Check constraints
    ASSERT_LE(posTraj.maxCoeff(), x0(0));  // actually not an explicit constraint
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
