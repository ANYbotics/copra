/*
 * Copyright 2016-2019 CNRS-UM LIRMM, CNRS-AIST JRL
 * Copyright 2020 ANYbotics AG
 */

#include "copra/solvers/utils.h"

namespace copra {

std::unique_ptr<SolverInterface> solverFactory(SolverFlag flag)
{
    switch (flag) {
#ifdef EIGEN_LSSOL_FOUND
    case SolverFlag::LSSOL:
        return std::unique_ptr<LSSOLSolver>(new LSSOLSolver());
#endif
#ifdef EIGEN_GUROBI_FOUND
    case SolverFlag::GUROBIDense:
        return std::unique_ptr<GUROBISolver>(new GUROBISolver());
#endif
    case SolverFlag::QLD:
        return std::unique_ptr<QLDSolver>(new QLDSolver());
#ifdef EIGEN_OSQP_FOUND
    case SolverFlag::OSQP:
        return std::unique_ptr<OSQPSolver>(new OSQPSolver());
#endif
    // case SolverFlag::QuadProgSparse:
    //    return std::make_unique<QuadProgSparseSolver>();
    case SolverFlag::QuadProgDense:
        return std::unique_ptr<QuadProgDenseSolver>(new QuadProgDenseSolver());
    case SolverFlag::qpOASES:
    case SolverFlag::DEFAULT:
        return std::unique_ptr<qpOASESSolver>(new qpOASESSolver());
    default:
        throw std::runtime_error("Invalid solver flag");
    }
}

} // namespace copra
