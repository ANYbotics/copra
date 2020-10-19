/*
 * Copyright 2016-2019 CNRS-UM LIRMM, CNRS-AIST JRL
 * Copyright 2020 ANYbotics AG
 */

#include "copra/solvers/qpOASESSolver.h"
#include <iostream>

namespace copra {

/*
 * qpOASES
 */

qpOASESSolver::qpOASESSolver() {
  solver_ = std::make_unique<qpOASES::SQProblem>();
}

int qpOASESSolver::SI_fail() const {
  return static_cast<int>(!solver_->isSolved());
}

void qpOASESSolver::SI_inform() const {
  std::cout << "qpOASES solver status is " << solver_->getStatus() << "\n";
}

void qpOASESSolver::SI_printLevel(int pl) {
  qpOASES::PrintLevel printLevel;
  switch (pl) {
    case 0:
      printLevel = qpOASES::PL_NONE;
      break;
    case 1:
      printLevel = qpOASES::PL_LOW;
      break;
    case 2:
      printLevel = qpOASES::PL_MEDIUM;
      break;
    case 3:
    default:
      printLevel = qpOASES::PL_HIGH;
  }
  solver_->setPrintLevel(printLevel);
}

const Eigen::VectorXd& qpOASESSolver::SI_result() const {
  return solution_;
}

void qpOASESSolver::SI_problem(int nrVar, int nrEq, int nrInEq) {
  if (numVariables_ != nrVar || numConstraints_ != (nrEq + nrInEq)) {
    numVariables_ = nrVar;
    numConstraints_ = nrEq + nrInEq;

    xOpt_ = new qpOASES::real_t[numVariables_];
    yOpt_ = new qpOASES::real_t[numVariables_ + numConstraints_];
    solution_.resize(numVariables_, 1);
    doInitWorkspace_ = true;

    qpOASES::Options op;
    op.setToMPC();
    op.printLevel = solver_ != nullptr ? solver_->getPrintLevel() : qpOASES::PL_NONE;
    solver_.reset();
    solver_ = std::make_unique<qpOASES::SQProblem>(numVariables_, numConstraints_, qpOASES::HessianType::HST_POSDEF);
    solver_->setOptions(op);
  }
}

bool qpOASESSolver::SI_solve(const Eigen::MatrixXd& Q, const Eigen::VectorXd& c, const Eigen::MatrixXd& Aeq, const Eigen::VectorXd& beq,
                             const Eigen::MatrixXd& Aineq, const Eigen::VectorXd& bineq, const Eigen::VectorXd& XL,
                             const Eigen::VectorXd& XU) {
  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> Q_(Q);
  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> A_(Aeq.rows() + Aineq.rows(), Aeq.cols());
  const Eigen::VectorXd& g_(c);
  const Eigen::VectorXd& XL_(XL);
  const Eigen::VectorXd& XU_(XU);
  Eigen::VectorXd bl_(Aeq.rows() + Aineq.rows());
  Eigen::VectorXd bu_(Aeq.rows() + Aineq.rows());

  A_.topRows(Aeq.rows()) = Aeq;
  A_.bottomRows(Aineq.rows()) = Aineq;
  bl_.head(Aeq.rows()) = beq;
  bl_.tail(Aineq.rows()).setConstant(-std::numeric_limits<double>::max());
  bu_.head(Aeq.rows()) = beq;
  bu_.tail(Aineq.rows()) = bineq;

  // Maximum number of working set recalculations. Refer to https://github.com/ANYbotics/qpOASES/blob/master/doc/manual.pdf for more info.
  qpOASES::int_t numIterations = 500;
  qpOASES::real_t cpuTime = 0.02;
  int exitFlag = -1;

  if (doInitWorkspace_) {
    exitFlag = solver_->init(Q_.data(), g_.data(), A_.data(), XL_.data(), XU_.data(), bl_.data(), bu_.data(), numIterations);
  } else {
    exitFlag = solver_->hotstart(Q_.data(), g_.data(), A_.data(), XL_.data(), XU_.data(), bl_.data(), bu_.data(), numIterations, &cpuTime);
  }

  if (exitFlag != qpOASES::SUCCESSFUL_RETURN) {
    std::cerr << "qpOASES could not be successfully initialized! Init returned exit flag = " << exitFlag << "\n";
  } else {
    exitFlag = solver_->getPrimalSolution(xOpt_);
    solver_->getPrimalSolution(solution_.data());
    solver_->getDualSolution(yOpt_);

    if (exitFlag != qpOASES::SUCCESSFUL_RETURN) {
      std::cerr << "qpOASES did not solve the problem! Returned exit flag = " << exitFlag << "\n";
    }
  }

  doInitWorkspace_ = (exitFlag != qpOASES::SUCCESSFUL_RETURN);

  return (exitFlag == qpOASES::SUCCESSFUL_RETURN);
}

}  // namespace copra
