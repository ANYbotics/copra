/*
 * Copyright 2016-2019 CNRS-UM LIRMM, CNRS-AIST JRL
 * Copyright 2020 ANYbotics AG
 */

#include <iostream>

#include "copra/solvers/OSQPSolver.h"

namespace copra
{

/*
 * OSQP
 */

OSQPSolver::OSQPSolver()
{
}

int OSQPSolver::SI_fail() const
{
    return static_cast<int>(solver_.status());
}

void OSQPSolver::SI_inform() const
{
    solver_.inform(std::cout);
}

int OSQPSolver::SI_iter() const
{
    return static_cast<int>(solver_.iter());
}

int OSQPSolver::SI_maxIter() const
{
    return static_cast<int>(solver_.maxIter());
}

void OSQPSolver::SI_maxIter(int maxIter)
{
    solver_.maxIter(maxIter);
}

bool OSQPSolver::SI_warmStart() const
{
    return solver_.warmStart();
}

void OSQPSolver::SI_warmStart(bool w)
{
    solver_.warmStart(w);
}

void OSQPSolver::SI_printLevel(int pl)
{
    solver_.verbose(pl != 0);
}

void OSQPSolver::SI_feasibilityTolerance(double tol)
{
    solver_.absConvergenceTol(tol);
}

const Eigen::VectorXd& OSQPSolver::SI_result() const
{
    return solver_.result();
}

void OSQPSolver::SI_problem(int nrVar, int nrEq, int nrInEq)
{
    solver_.problem(nrVar, nrEq + nrInEq);
}

bool OSQPSolver::SI_solve(const Eigen::MatrixXd& Q, const Eigen::VectorXd& c,
    const Eigen::MatrixXd& Aeq, const Eigen::VectorXd& beq,
    const Eigen::MatrixXd& Aineq, const Eigen::VectorXd& bineq,
    const Eigen::VectorXd& XL, const Eigen::VectorXd& XU)
{
    Eigen::MatrixXd A(Aeq.rows() + Aineq.rows(), Aeq.cols());
    A.topRows(Aeq.rows()) = Aeq;
    A.bottomRows(Aineq.rows()) = Aineq;
    Eigen::VectorXd bl(Aeq.rows() + Aineq.rows());
    Eigen::VectorXd bu(Aeq.rows() + Aineq.rows());
    bl.head(Aeq.rows()) = beq;
    bl.tail(Aineq.rows()).setConstant(-std::numeric_limits<double>::max());
    bu.head(Aeq.rows()) = beq;
    bu.tail(Aineq.rows()) = bineq;

    return solver_.solve(Q, c, A, bl, bu, XL, XU);
}

} // namespace copra
