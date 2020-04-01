/*
 * Copyright 2016-2019 CNRS-UM LIRMM, CNRS-AIST JRL
 * Copyright 2020 ANYbotics AG
 */

#include "copra/debugUtils.h"
#include "copra/systems/TimeInvariantSystem.h"

namespace copra {

void TimeInvariantSystem::reset(const Eigen::MatrixXd& state, const Eigen::MatrixXd& control,
    const Eigen::VectorXd& bias, const Eigen::VectorXd& xInit, int numberOfSteps)
{
    if (xInit.rows() != state.rows()) {
        DOMAIN_ERROR_EXCEPTION(throwMsgOnRows("xInit", "state", xInit, state));
    } else if (state.rows() != state.cols()) {
        DOMAIN_ERROR_EXCEPTION(throwMsgOnSquareMat("state", state));
    } else if (xInit.rows() != control.rows()) {
        DOMAIN_ERROR_EXCEPTION(throwMsgOnRows("xInit", "control", xInit, control));
    } else if (xInit.rows() != bias.rows()) {
        DOMAIN_ERROR_EXCEPTION(throwMsgOnRows("xInit", "state", xInit, bias));
    }
    int stateDim = static_cast<int>(state.cols());
    int inputDim = static_cast<int>(control.cols());
    System::reset(stateDim, inputDim, numberOfSteps);
    x0 = xInit;
    A_ = state;
    B_ = control;
    d_ = bias;
}

void TimeInvariantSystem::update() noexcept
{
    /*
     * Loop invariant: x_k = phi_k * x_0 + psi_k * U + xi_k
     */

    if (isUpdated) {
        RUNTIME_ERROR_EXCEPTION("You need to reset() your system before update()-ing it.");
    }

    /*
     * x_1 = A * x_0 + B * u_0 + d
     * => phi_1 = A
     * => psi_1 = B
     * => xi_1 = d
     *
     * Note that Phi, Psi and xi have been initially set to zero.
     */
    Phi.block(xDim, 0, xDim, xDim) = A_;
    Psi.block(xDim, 0, xDim, uDim) = B_;
    xi.segment(xDim, xDim) = d_;

    for (auto k = 2; k < nrXStep; ++k) {
        /*
         * x_{k+1} = A * x_k + B * u_k + d
         * => phi_{k+1} = A * phi_k
         * => psi_{k+1} = A * psi_k + [0 0 ... 0 B 0 ... 0]
         * => xi_{k+1} = A * xi_k + d
         */
        Phi.block(k * xDim, 0, xDim, xDim).noalias() = A_ * Phi.block((k - 1) * xDim, 0, xDim, xDim);
        Psi.block(k * xDim, 0, xDim, uDim).noalias() = A_ * Psi.block((k - 1) * xDim, 0, xDim, uDim);
        for (auto j = 1; j < k; ++j) {
            // the computation of Psi in TimeVariantSystem::update() could be faster than this for loop, try it if performance needed
            Psi.block(k * xDim, j * uDim, xDim, uDim) = Psi.block((k - 1) * xDim, (j - 1) * uDim, xDim, uDim);
        }
        xi.segment(k * xDim, xDim).noalias() = A_ * xi.segment((k - 1) * xDim, xDim) + d_;
    }

    isUpdated = true;
}

} // namespace copra
