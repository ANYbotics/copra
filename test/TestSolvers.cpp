/*
 * Copyright 2016-2019 CNRS-UM LIRMM, CNRS-AIST JRL
 * Copyright 2020 ANYbotics AG
 */

#include <iostream>
#include <numeric>

#include <Eigen/Core>

#include <gtest/gtest.h>

#include "copra/solvers/all.h"

#include "time_invariant_systems.h"

TEST_F(Problem, QuadProgTest) {  // NOLINT
    copra::QuadProgDenseSolver qpQuadProg;

    qpQuadProg.SI_problem(nrvars, nreqs, nrineqs);
    ASSERT_TRUE(qpQuadProg.SI_solve(Q, c, Aeq, beq, Aineq, bineq, XL, XU));

    ASSERT_EQ(qpQuadProg.SI_fail(), 0);
}

TEST_F(Problem, qpOASESTest) {  // NOLINT
  copra::qpOASESSolver qpOasesSolver;

  qpOasesSolver.SI_problem(nrvars, nreqs, nrineqs);
  ASSERT_TRUE(qpOasesSolver.SI_solve(Q, c, Aeq, beq, Aineq, bineq, XL, XU));

  ASSERT_EQ(qpOasesSolver.SI_fail(), 0);
}

TEST_F(Problem, QLDOnQuadProgTest) {  // NOLINT
    copra::QLDSolver qpQLD;
    copra::QuadProgDenseSolver qpQuadProg;

    qpQLD.SI_problem(nrvars, nreqs, nrineqs);
    qpQuadProg.SI_problem(nrvars, nreqs, nrineqs);
    ASSERT_TRUE(qpQLD.SI_solve(Q, c, Aeq, beq, Aineq, bineq, XL, XU));
    ASSERT_TRUE(qpQuadProg.SI_solve(Q, c, Aeq, beq, Aineq, bineq, XL, XU));

    Eigen::VectorXd resQLD = qpQLD.SI_result();
    Eigen::VectorXd resQuadProg = qpQuadProg.SI_result();
    EXPECT_TRUE(resQuadProg.isApprox(resQLD));
    ASSERT_EQ(qpQLD.SI_fail(), 0);
    ASSERT_EQ(qpQuadProg.SI_fail(), 0);
}

#ifdef EIGEN_LSSOL_FOUND
TEST_F(Problem, LSSOLOnQuadProgTest) {  // NOLINT
    copra::QuadProgDenseSolver qpQuadProg;
    copra::LSSOLSolver qpLSSOL;

    qpQuadProg.SI_problem(nrvars, nreqs, nrineqs);
    qpLSSOL.SI_problem(nrvars, nreqs, nrineqs);
    ASSERT_TRUE(qpQuadProg.SI_solve(Q, c, Aeq, beq, Aineq, bineq, XL, XU));
    ASSERT_TRUE(qpLSSOL.SI_solve(Q, c, Aeq, beq, Aineq, bineq, XL, XU));

    Eigen::VectorXd resQuadProg = qpQuadProg.SI_result();
    Eigen::VectorXd resLSSOL = qpLSSOL.SI_result();
    EXPECT_TRUE(resLSSOL.isApprox(resQuadProg));
    ASSERT_EQ(qpLSSOL.SI_fail(), 0);
}
#endif

#ifdef EIGEN_GUROBI_FOUND
TEST_F(Problem, GUROBIOnQuadProgTest) {  // NOLINT
    copra::QuadProgDenseSolver qpQuadProg;
    copra::GUROBISolver qpGUROBI;

    qpQuadProg.SI_problem(nrvars, nreqs, nrineqs);
    qpGUROBI.SI_problem(nrvars, nreqs, nrineqs);
    ASSERT_TRUE(qpQuadProg.SI_solve(Q, c, Aeq, beq, Aineq, bineq, XL, XU));
    ASSERT_TRUE(qpGUROBI.SI_solve(Q, c, Aeq, beq, Aineq, bineq, XL, XU));

    Eigen::VectorXd resQuadProg = qpQuadProg.SI_result();
    Eigen::VectorXd resGUROBI = qpGUROBI.SI_result();
    EXPECT_TRUE(resGUROBI.isApprox(resQuadProg, 1e-6));
    ASSERT_EQ(qpGUROBI.SI_fail(), GRB_OPTIMAL);
}
#endif

#ifdef EIGEN_OSQP_FOUND
TEST_F(Problem, OSQPOnQuadProgTest) {  // NOLINT
    copra::QuadProgDenseSolver qpQuadProg;
    copra::OSQPSolver osqp;

    qpQuadProg.SI_problem(nrvars, nreqs, nrineqs);
    osqp.SI_problem(nrvars, nreqs, nrineqs);
    ASSERT_TRUE(qpQuadProg.SI_solve(Q, c, Aeq, beq, Aineq, bineq, XL, XU));
    ASSERT_TRUE(osqp.SI_solve(Q, c, Aeq, beq, Aineq, bineq, XL, XU));

    Eigen::VectorXd resQuadProg = qpQuadProg.SI_result();
    Eigen::VectorXd resOSQP = osqp.SI_result();
    EXPECT_TRUE(resOSQP.isApprox(resQuadProg, 1e-6));
    ASSERT_EQ(osqp.SI_fail(), 1);
}
#endif
