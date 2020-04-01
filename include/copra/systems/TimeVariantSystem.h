/*
 * Copyright 2016-2019 CNRS-UM LIRMM, CNRS-AIST JRL
 * Copyright 2020 ANYbotics AG
 */

#pragma once

#include <Eigen/Core>

#include "copra/api.h"
#include "copra/systems/System.h"

namespace copra {

/*! Linear time-variant (LTV) system.
 *
 * An LTV system is defined by:
 * \f$x_{k+1} = A_{k} x_{k} + B_{k} u_{k} + d_k \f$
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
struct COPRA_DLLAPI TimeVariantSystem : public System {
    /**
     * Default constructor.
     */
    TimeVariantSystem() = default;

    /*! Create system from system dimensions and initial state.
     *
     * @param stateDim State dimension.
     * @param inputDim Input dimension.
     * @param xInit Initial state.
     * @param nbSteps Number of receding-horizon steps.
     */
    TimeVariantSystem(int stateDim, int inputDim, const Eigen::VectorXd& xInit, int nbSteps) {
      reset(stateDim, inputDim, xInit, nbSteps);
    }

    /*! Reset system from system dimensions and initial state.
     *
     * @param stateDim State dimension.
     * @param inputDim Input dimension.
     * @param xInit Initial state.
     * @param nbSteps Number of receding-horizon steps.
     */
    void reset(int stateDim, int inputDim, const Eigen::VectorXd& xInit, int numberOfSteps);

    /**
     * \brief Update the system.
     * Fill Phi, Psi, xi in System
     */
    void update() noexcept override;

    /*! Get the state matrix at a given time step.
       *
       * @param step Time step \f$k\f$.
       * @return State matrix \f$A_{k}\f$.
       */
    Eigen::MatrixXd getStateMatrix(int step) const override { return A_(step); }

    /*! Set the state matrix for all time steps.
     *
     * @param A Lambda function such that A(k) returns the state matrix \f$A_{k}\f$ at step \f$k\f$.
     */
    void setStateMatrix(std::function<Eigen::MatrixXd(int)> A) { A_ = std::move(A); }

    /*! Get the input matrix at a given time step.
     *
     * @param step Time step \f$k\f$.
     * @return Input matrix \f$B_{k}\f$.
     */
    Eigen::MatrixXd getInputMatrix(int step) const override { return B_(step); }

    /*! Set the input matrix for all time steps.
     *
     * @param B Lambda function such that B(k) returns the state matrix \f$B_{k}\f$ at step \f$k\f$.
     */
    void setInputMatrix(std::function<Eigen::MatrixXd(int)> B) { B_ = std::move(B); }

    /*! Get the drift vector at a given time step.
     *
     * @param step Time step \f$k\f$.
     * @return Drift vector \f$d_{k}\f$.
     */
    Eigen::VectorXd getDriftVector(int step) const override { return d_(step); }

    /*! Set the drift vector for all time steps.
     *
     * @param d Lambda function such that d(k) returns the drift vector \f$d_{k}\f$ at step \f$k\f$.
     */
    void setDriftVector(std::function<Eigen::VectorXd(int)> d) { d_ = std::move(d); }

protected:
    //! Function A(k) returning the state matrix at receding-horizon time k
    std::function<Eigen::MatrixXd(int)> A_;

    //! Function B(k) returning the input matrix at receding-horizon time k
    std::function<Eigen::MatrixXd(int)> B_;

    //! Function d(k) returning the drift vector at receding-horizon time k
    std::function<Eigen::VectorXd(int)> d_;
};

} // namespace copra
