/*
 * Copyright 2016-2019 CNRS-UM LIRMM, CNRS-AIST JRL
 * Copyright 2020 ANYbotics AG
 */

#pragma once

#include <Eigen/Core>

#include "copra/api.h"
#include "copra/systems/System.h"

namespace copra {

/*! Linear time-invariant (LTI) system.
 *
 * An LTI system is defined by:
 * \f$x_{k+1} = A x_{k} + B u_{k} + d\f$
 * where \f$x\f$ is the state, \f$u\f$ is the input, \f$A\f$ is the state matrix, \f$B\f$ is the input matrix, and \f$d\f$ is the drift.
 *
 * States and inputs can be grouped into stacked vectors \f$X = [x_0^T x_1^T ... x_N^T]^T\f$ and \f$U = [u_0^T u_1^T ... u_{N-1}^T]^T\f$
 * where \f$N\f$ is the number of steps of the receding horizon. Then, \f$dim(X) = (N + 1) dim(x)\f$ and \f$dim(U) = N dim(u)\f$.
 *
 * Applying the system equation recursively, the stacked state vector can be written as
 * \f$X = \Phi x_{0} + \Psi U + \Xi\f$
 * following the single-shooting method. See for instance "A condensed and sparse qp formulation for predictive control" (Jerez et al.) or
 * its application in "Model preview control in multi-contact motion-application to a humanoid robot" (Audren et al.).
 */
struct COPRA_DLLAPI TimeInvariantSystem : public System {
    /**
     * Default constructor.
     */
    TimeInvariantSystem() = default;

   /**
    * Constructor of the class.
    * \param state The state matrix of the system.
    * \param control The control matrix of the system.
    * \param bias The bias vector of the system.
    * \param xInit The initial state.
    * \param numberOfSteps The number of step to perform.
    * \throw std::domain_error if the dimension of the matrices mismatch.
    */
    TimeInvariantSystem(const Eigen::MatrixXd& state, const Eigen::MatrixXd& control,
        const Eigen::VectorXd& bias, const Eigen::VectorXd& xInit, int numberOfSteps) {
      reset(state, control, bias, xInit, numberOfSteps);
    }

    /*! Reset system.
     *
     * \param state The state matrix of the system.
     * \param control The control matrix of the system.
     * \param bias The bias vector of the system.
     * \param xInit The initial state.
     * \param numberOfSteps The number of step to perform.
     * \throw std::domain_error if the dimension of the matrices mismatch.
     */
    void reset(const Eigen::MatrixXd& state, const Eigen::MatrixXd& control,
        const Eigen::VectorXd& bias, const Eigen::VectorXd& xInit, int numberOfSteps);

    /**
     * \brief Update the system.
     * Fill Phi, Psi, xi in System
     */
    void update() noexcept override;

    /*! Get the state matrix at a given time step.
     *
     * @param step Time step \f$k\f$.
     * @return State matrix \f$A_{k} = A\f$ (LTI system).
     *
     * \note Returns a copy.
     */
    Eigen::MatrixXd getStateMatrix(int /* step */) const override { return A_; }

    /*! Get the input matrix at a given time step.
     *
     * @param step Time step \f$k\f$.
     * @return Input matrix \f$B_{k} = B\f$ (LTI system).
     *
     * \note Returns a copy.
     */
    Eigen::MatrixXd getInputMatrix(int /* step */) const override { return B_; }

    /*! Get the drift vector at a given time step.
     *
     * @param step Time step \f$k\f$.
     * @return Drift vector \f$d_{k} = d\f$ (LTI system).
     *
     * \note Returns a copy.
     */
    Eigen::VectorXd getDriftVector(int /* step */) const override { return d_; }

protected:
    Eigen::MatrixXd A_; /**< Buffer for state matrices */
    Eigen::MatrixXd B_; /**< Buffer for control matrices */
    Eigen::VectorXd d_; /**< Buffer for bias (drift) vectors */
};

} // namespace copra
