/*
 * Copyright 2016-2019 CNRS-UM LIRMM, CNRS-AIST JRL
 * Copyright 2020 ANYbotics AG
 */

#include "copra/debugUtils.h"
#include "copra/systems/TimeVariantSystem.h"

namespace copra {

void TimeVariantSystem::reset(int stateDim, int inputDim, const Eigen::VectorXd& xInit, int numberOfSteps)
{
    System::reset(stateDim, inputDim, numberOfSteps);
    x0 = xInit;
}

void TimeVariantSystem::update() noexcept
{
  /*
   * Loop invariant: x_k = phi_k * x_0 + psi_k * U + xi_k
   */

  if (isUpdated) {
    RUNTIME_ERROR_EXCEPTION("You need to reset() your system before update()-ing it.");
  }

  /*
   * x_1 = A_0 * x_0 + B_0 * u_0 + d_0
   * => phi_1 = A_0
   * => psi_1 = B_0
   * => xi_1 = d_0
   *
   * Note that Phi, Psi and xi have been initially set to zero.
   */
  Phi.block(xDim, 0, xDim, xDim) = getStateMatrix(0);
  Psi.block(xDim, 0, xDim, uDim) = getInputMatrix(0);
  xi.segment(xDim, xDim) = getDriftVector(0);

  for (int k = 1; k < nrXStep - 1; k++) {
    const Eigen::MatrixXd A_k = getStateMatrix(k);
    const Eigen::MatrixXd B_k = getInputMatrix(k);
    const Eigen::MatrixXd d_k = getDriftVector(k);
    /*
     * x_{k+1} = A_k * x_k + B_k * u_k + d_k
     * => phi_{k+1} = A_k * phi_k
     * => psi_{k+1} = A_k * psi_k + [0 0 ... 0 B_k 0 ... 0]
     * => xi_{k+1} = A_k * xi_k + d_k
     */
    Phi.block((k + 1) * xDim, 0, xDim, xDim).noalias() = A_k * Phi.block(k * xDim, 0, xDim, xDim);
    Psi.block((k + 1) * xDim, 0, xDim, fullUDim).noalias() = A_k * Psi.block(k * xDim, 0, xDim, fullUDim);
    Psi.block((k + 1) * xDim, k * uDim, xDim, uDim) += B_k;
    xi.segment((k + 1) * xDim, xDim).noalias() = A_k * xi.segment(k * xDim, xDim) + d_k;
  }

  isUpdated = true;
}

} // namespace copra
