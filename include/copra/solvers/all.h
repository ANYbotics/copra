/*
 * Copyright 2016-2019 CNRS-UM LIRMM, CNRS-AIST JRL
 * Copyright 2020 ANYbotics AG
 */

#pragma once

#ifdef EIGEN_QUADPROG_FOUND
#include "QuadProgSolver.h"
#endif
#ifdef EIGEN_QLD_FOUND
#include "QLDSolver.h"
#endif
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
