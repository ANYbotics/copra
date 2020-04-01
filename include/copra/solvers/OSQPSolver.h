/*
 * Copyright 2016-2019 CNRS-UM LIRMM, CNRS-AIST JRL
 * Copyright 2020 ANYbotics AG
 */

#pragma once

#include <Eigen/Core>
#include <eigen-osqp/OSQP.h>

#include "copra/api.h"
#include "copra/SolverInterface.h"

namespace copra {

/**
 * OSQP solver for dense matrix.
 */

// TODO: Enable sparse matrix
class COPRA_DLLAPI OSQPSolver : public SolverInterface {
public:
    /**
       * QLDSolver default constructor
       */
    OSQPSolver();

    /**
     * Get information of eventual fail's solver output as define by the
     * solver documentation.
     * \return 0 The optimality conditions are satisfied.
     * \return 1 The algorithm has been stopped after too many iterations.
     * \return 2 Termination accuracy insufficient to satisfy convergence
     * criterion.
     * \return 3 Internal inconsistency of QL, division by zero.
     * \return 4 Numerical instability prevents successful termination.
     * \return 5 Length of a working array is too short.
     * \return >100 Constraints are inconsistent and fail=100+ICON, where ICON
     * denotes a constraint causing the conflict.
     */
    int SI_fail() const override;

    void SI_inform() const override;
    int SI_iter() const override;
    int SI_maxIter() const override;
    void SI_maxIter(int maxIter) override;
    void SI_printLevel(int pl) override;
    void SI_feasibilityTolerance(double tol) override;
    bool SI_warmStart() const override;
    void SI_warmStart(bool w) override;

    const Eigen::VectorXd& SI_result() const override;
    void SI_problem(int nrVar, int nrEq, int nrInEq) override;


    bool SI_solve(const Eigen::MatrixXd& Q, const Eigen::VectorXd& c,
        const Eigen::MatrixXd& Aeq, const Eigen::VectorXd& beq,
        const Eigen::MatrixXd& Aineq, const Eigen::VectorXd& bineq,
        const Eigen::VectorXd& XL, const Eigen::VectorXd& XU) override;

    Eigen::OSQP& baseSolver() noexcept { return solver_; }

private:
    Eigen::OSQP solver_;
};

} // namespace copra
