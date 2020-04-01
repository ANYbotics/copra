/*
 * Copyright 2016-2019 CNRS-UM LIRMM, CNRS-AIST JRL
 * Copyright 2020 ANYbotics AG
 */

#pragma once

#include <Eigen/Core>

#include "copra/api.h"

namespace copra {

/*! Linear time-variant system.
 *
 * A linear time-variant system is defined by:
 * \f$x_{k+1} = A_{k} x_{k} + B_{k} u_{k} + d_k \f$
 * where \f$x\f$ is the state, \f$u\f$ is the input, \f$A\f$ is the state matrix, \f$B\f$ is the input matrix, and \f$d\f$ is the drift.
 *
 * States and inputs can be grouped into stacked vectors \f$X = [x_0^T x_1^T ... x_N^T]^T\f$ and \f$U = [u_0^T u_1^T ... u_{N-1}^T]^T\f$
 * where \f$N\f$ is the number of steps of the receding horizon. Then, \f$dim(X) = (N + 1) dim(x)\f$ and \f$dim(U) = N dim(u)\f$.
 *
 * Applying the system equation recursively, the stacked state vector can be written as:
 * \f$X = \Phi x_{0} + \Psi U + \Xi\f$.
 * See for instance "A condensed and
 * sparse qp formulation for predictive control" (Jerez et al.) or its application in "Model preview control in multi-contact
 * motion-application to a humanoid robot" (Audren et al.).
 */
struct COPRA_DLLAPI System {
    /**
     * \brief Default constructor.
     */
    System() = default;

    /** \brief Constructor from system dimensions.
     *
     * @param stateDim Dimension of system state.
     * @param inputDim Dimension of system input.
     * @param nbSteps Number of time steps in the receding horizon.
     */
    System(int stateDim, int inputDim, int nbSteps) {
      reset(stateDim, inputDim, nbSteps);
    }

    /** \brief Reset system dimensions.
     *
     * @param stateDim Dimension of system state.
     * @param inputDim Dimension of system input.
     * @param nbSteps Number of time steps in the receding horizon.
     */
    void reset(int stateDim, int inputDim, int numberOfSteps);

    /**
     * \brief Update the system.
     * Fill Phi, Psi, xi in System
     */
    virtual void update() noexcept = 0;

    /**
     * \brief Update the initial state.
     */
    void xInit(const Eigen::VectorXd& xInit) { x0 = xInit; }

    /**! Set initial state.
     *
     * @param state Initial state.
     */
    void setInitialState(const Eigen::VectorXd& state) { xInit(state); }

    /*! Get the state matrix at a given time step.
     *
     * @param step Time step \f$k\f$.
     * @return State matrix \f$A_{k}\f$.
     */
    virtual Eigen::MatrixXd getStateMatrix(int step) const = 0;

    /*! Get the input matrix at a given time step.
     *
     * @param step Time step \f$k\f$.
     * @return Input matrix \f$B_{k}\f$.
     */
    virtual Eigen::MatrixXd getInputMatrix(int step) const = 0;

    /*! Get the drift vector at a given time step.
     *
     * @param step Time step \f$k\f$.
     * @return Drift vector \f$d_{k}\f$.
     */
    virtual Eigen::VectorXd getDriftVector(int step) const = 0;

    bool isUpdated = false; /**< State whether or not the preview system has been updated. This is done when calling the solve function of a mpc. Calling \see system or setting this to false will force a new update*/
    int nrUStep = 0; /**< The number of iteration to perform for U (it is the dimension of \f$U\f$, \f$N\f$). */
    int nrXStep = 0; /**< The number of iteration to perform for X (it is the dimension of \f$X\f$, \f$N+1\f$). */
    int xDim = 0; /**< The dimension of the state vector */
    int uDim = 0; /**< The dimension of the control vector */
    int fullXDim = 0; /**< The full dimension of the state vector (xDim*nbStep) */
    int fullUDim = 0; /**< The full dimension of the control vector (uDim*nbStep) */
    Eigen::VectorXd x0; /**< The initial state */
    Eigen::MatrixXd Phi; /**< The full (after nrStep) state matrix */
    Eigen::MatrixXd Psi; /**< The full (after nrStep) control matrix */
    Eigen::VectorXd xi; /**< The full (after nrStep) bias vector */
};

} // namespace copra
