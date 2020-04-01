/*
 * Copyright 2016-2019 CNRS-UM LIRMM, CNRS-AIST JRL
 * Copyright 2020 ANYbotics AG
 */

#include "copra/debugUtils.h"
#include "copra/systems/System.h"

namespace copra {

void System::reset(int stateDim, int inputDim, int numberOfSteps)
{
  if (numberOfSteps <= 0) {
    DOMAIN_ERROR_EXCEPTION("The number of steps should be a positive number! ");
  }

  xDim = stateDim;
  uDim = inputDim;

  isUpdated = false;
  nrUStep = numberOfSteps;
  nrXStep = numberOfSteps + 1;
  fullXDim = xDim * nrXStep;
  fullUDim = uDim * nrUStep;

  Phi.resize(fullXDim, xDim);
  Psi.resize(fullXDim, fullUDim);
  xi.resize(fullXDim);

  Phi.setZero();
  Phi.block(0, 0, xDim, xDim) = Eigen::MatrixXd::Identity(xDim, xDim);
  Psi.setZero();
  xi.setZero();
}

} // namespace copra
