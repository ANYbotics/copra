/*
 * Copyright 2016-2019 CNRS-UM LIRMM, CNRS-AIST JRL
 * Copyright 2020 ANYbotics AG
 */

#pragma once

#include "copra/api.h"

#include "QuadProgSolver.h"
#include "QLDSolver.h"
#include "qpOASESSolver.h"
#ifdef EIGEN_LSSOL_FOUND
#include "LSSOLSolver.h"
#endif
#ifdef EIGEN_GUROBI_FOUND
#include "GUROBISolver.h"
#endif
#ifdef EIGEN_OSQP_FOUND
#include "OSQPSolver.h"
#endif
#include <memory>
#include <utility>

namespace copra {

#ifdef __GNUC__
#pragma GCC diagnostic push
// Work around GCC (< 6) bug see: https://gcc.gnu.org/bugzilla/show_bug.cgi?id=43407
#pragma GCC diagnostic ignored "-Wattributes"
#endif

/**
 * Enum class that handles flag for selecting a qp solver.
 */
enum class COPRA_DLLAPI SolverFlag {
    DEFAULT, /**< Default solver (QuadProgDense solver) */
#ifdef EIGEN_LSSOL_FOUND
    LSSOL, /**< Stanford LSSOL solver */
#endif
#ifdef EIGEN_GUROBI_FOUND
    GUROBIDense, /**< Gurobi quadratic dense solver */
#endif
    QLD, /**< Scilab QLD solver */
#ifdef EIGEN_OSQP_FOUND
    OSQP,
#endif
    QuadProgDense, /**< DenseMatrix version of QuadProg solver */
    // QuadProgSparse
    qpOASES,
};

#ifdef __GNUC__
#pragma GCC diagnostic pop
#endif

/**
 * Helper function to get an unique pointer to a desired solver.
 * \param flag Flag of the solver.
 * \return An unique pointer to the desired solver.
 */
COPRA_DLLAPI std::unique_ptr<SolverInterface> solverFactory(SolverFlag flag);

} // namespace copra
