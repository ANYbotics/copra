/*
 * Copyright 2016-2019 CNRS-UM LIRMM, CNRS-AIST JRL
 * Copyright 2020 ANYbotics AG
 */

#pragma once

#include <Eigen/Core>
#include <memory>
#include <qp_oases/qpOASES.hpp>
#include <vector>

#include "copra/SolverInterface.h"
#include "copra/api.h"

namespace copra {

//! qpOASES solver for dense matrix.

class COPRA_DLLAPI qpOASESSolver : public SolverInterface {
 public:
  /*!
   * qpOASES default constructor
   */
  qpOASESSolver();

  /*!
   * Get information of eventual fail's solver output as define by the
   * solver documentation.
   * \return 0 If the solver succeeded in finding a solution to the problem.
   * \return 1 If the solver failed to find a solution to the problem.
   */
  int SI_fail() const override;

  //! Prints the solver status.
  void SI_inform() const override;

  //! Get the maximum number of internal solver iterations.
  int SI_maxIter() const override;

  /*! Set the maximum number of internal solver iterations.
   *
   * @param maxIter Maximum number of working space recalculations.
   */
  void SI_maxIter(int maxIter) override;

  //! Get the maximum allowed CPU time to compute a QP solution
  double maxCpuTime() const;

  /*! Set the maximum allowed CPU time to compute a QP solution
   *
   * @param maxCpuTime Maximum allowed CPU time.
   */
  void maxCpuTime(double maxCpuTime);

  /*! Set the logger level for the solver.
   *
   * @param pl verbosity level: 0 = NONE, 1 = LOW, 2 = MEDIUM and 3 = HIGH verbosity level.
   *
   */
  void SI_printLevel(int pl) override;

  /*! Get the solution of the optimization.
   *
   * @return Vector of solution.
   */
  const Eigen::VectorXd& SI_result() const override;

  /*! Setup the solver. Re-create the problem if the problem dimension has changed.
   *
   * @param nrVar Number of optimization variables.
   * @param nrEq Number of equality constraints.
   * @param nrInEq Number of inequality constraints.
   */
  void SI_problem(int nrVar, int nrEq, int nrInEq) override;

  /*! Call the solver to solve the optimization problem.
   *
   * @note Must be called after SI_problem()
   *
   * @param Q Hessian matrix of the objective function
   * @param c Linear term of the objective function
   * @param Aeq Equality constraint matrix
   * @param beq Equality constraint vector
   * @param Aineq Inequality constraint matrix
   * @param bineq Inequality constraint vector
   * @param XL Lower bound of the optimization variables
   * @param XU Upper bound of the optimization variables
   * @return True, iff optimization was successful.
   */
  bool SI_solve(const Eigen::MatrixXd& Q, const Eigen::VectorXd& c, const Eigen::MatrixXd& Aeq, const Eigen::VectorXd& beq,
                const Eigen::MatrixXd& Aineq, const Eigen::VectorXd& bineq, const Eigen::VectorXd& XL, const Eigen::VectorXd& XU) override;

  //! Get the base solver pointer.
  qpOASES::QProblem& baseSolver() noexcept { return *solver_; }

 private:
  //! qpOASES solver.
  std::unique_ptr<qpOASES::SQProblem> solver_;

  //! Solution vector.
  Eigen::VectorXd solution_;

  //! Number of optimization variables.
  int numVariables_ = 0;

  //! Total number of constraints.
  int numConstraints_ = 0;

  //! Primal solution vector.
  std::vector<qpOASES::real_t> xOpt_;

  //! Dual solution vector.
  std::vector<qpOASES::real_t> yOpt_;

  //! Boolean flag for re-initializing the problem.
  bool doInitWorkspace_ = true;

  //! Maximum number of internal solver iterations.
  int numIterations_ = 500;

  //! Maximum allowed CPU time to compute a QP solution (disabled if zero, which is the default).
  double maxCpuTime_ = 0.0;
};

}  // namespace copra
