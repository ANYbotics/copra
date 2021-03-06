/*
 * Copyright 2016-2019 CNRS-UM LIRMM, CNRS-AIST JRL
 * Copyright 2020 ANYbotics AG
 */

#include "copra/costFunctions.h"
#include "copra/AutoSpan.h"

namespace copra {

/*************************************************************************************************
 *                                        Cost Function                                          *
 *************************************************************************************************/

CostFunction::CostFunction(std::string&& name)
    : name_(std::move(name))
    , fullSizeEntry_(false)
{
}

void CostFunction::autoSpan()
{
}

void CostFunction::initializeCost(const System& ps)
{
    Q_.resize(ps.fullUDim, ps.fullUDim);
    c_.resize(ps.fullUDim);
}

/*************************************************************************************************
 *                                  Trajectory Cost Function                                     *
 *************************************************************************************************/

void TrajectoryCost::autoSpan()
{
    auto max_dim = std::max(M_.rows(), std::max(weights_.rows(), p_.rows()));
    AutoSpan::spanMatrix(M_, max_dim);
    AutoSpan::spanVector(p_, max_dim);
    AutoSpan::spanVector(weights_, max_dim);
}

void TrajectoryCost::initializeCost(const System& ps)
{
    CostFunction::initializeCost(ps);
    if (M_.rows() != p_.rows())
        DOMAIN_ERROR_EXCEPTION(throwMsgOnRowsAskAutoSpan("M", "p", M_, p_));

    if (M_.cols() == ps.xDim) {
        Q_.setZero();
        c_.setZero();
    } else if (M_.cols() == ps.fullXDim) {
        fullSizeEntry_ = true;
    } else {
        DOMAIN_ERROR_EXCEPTION(throwMsgOnColsOnPSXDim("M", M_, &ps));
    }
}

void TrajectoryCost::update(const System& ps)
{
    if (fullSizeEntry_) {
        Eigen::MatrixXd tmp{ M_ * ps.Psi };
        Q_ = tmp.transpose() * weights_.asDiagonal() * tmp;
        c_ = (M_ * (ps.Phi * ps.x0 + ps.xi) - p_).transpose() * weights_.asDiagonal() * tmp;
    } else {
        for (int i = 0; i < ps.nrXStep; ++i) { // Can be optimized. Lot of sums of zero here
            Eigen::MatrixXd tmp{ M_ * ps.Psi.block(i * ps.xDim, 0, ps.xDim, ps.fullUDim) };
            Q_.noalias() += tmp.transpose() * weights_.asDiagonal() * tmp;
            c_.noalias() += (M_ * (ps.Phi.block(i * ps.xDim, 0, ps.xDim, ps.xDim) * ps.x0 + ps.xi.segment(i * ps.xDim, ps.xDim)) - p_).transpose() * weights_.asDiagonal() * tmp;
        }
    }
}

/*************************************************************************************************
 *                           Simple Trajectory Cost Function                                     *
 *************************************************************************************************/

void SimpleTrajectoryCost::autoSpan()
{
    auto max_dim = std::max(weights_.rows(), p_.rows());
    AutoSpan::spanVector(p_, max_dim);
    AutoSpan::spanVector(weights_, max_dim);
}

void SimpleTrajectoryCost::initializeCost(const System& ps)
{
    CostFunction::initializeCost(ps);

    if (p_.rows() == ps.xDim) {
        Q_.setZero();
        c_.setZero();
    } else if (p_.rows() == ps.fullXDim) {
        fullSizeEntry_ = true;
    } else {
        DOMAIN_ERROR_EXCEPTION(throwMsgOnRowsOnPSXDim("p", p_, &ps));
    }
}

void SimpleTrajectoryCost::update(const System& ps)
{
    if (fullSizeEntry_) {
        Q_.noalias() = ps.Psi.transpose() * weights_.asDiagonal() * ps.Psi;
        c_.noalias() = ((ps.Phi * ps.x0 + ps.xi) - p_).cwiseProduct(weights_).transpose() * ps.Psi;
    } else {
        for (int i = 0; i < ps.nrXStep; ++i) { // Can be optimized. Lot of sums of zero here
            Eigen::MatrixXd tmp{ps.Psi.block(i * ps.xDim, 0, ps.xDim, ps.fullUDim)};
            Q_.noalias() += tmp.transpose() * weights_.asDiagonal() * tmp;
            c_.noalias() += ((ps.Phi.block(i * ps.xDim, 0, ps.xDim, ps.xDim) * ps.x0 + ps.xi.segment(i * ps.xDim, ps.xDim)) - p_).transpose() * weights_.asDiagonal() * tmp;
        }
    }
}

/*************************************************************************************************
 *                                    Target Cost Function                                       *
 *************************************************************************************************/

void TargetCost::initializeCost(const System& ps)
{
    CostFunction::initializeCost(ps);
    if (M_.rows() != p_.rows())
        DOMAIN_ERROR_EXCEPTION(throwMsgOnRowsAskAutoSpan("M", "p", M_, p_));

    if (M_.cols() != ps.xDim)
        DOMAIN_ERROR_EXCEPTION(throwMsgOnRowsOnPSxDim("M", M_, &ps));
}

void TargetCost::update(const System& ps)
{
    Eigen::MatrixXd tmp{ M_ * ps.Psi.bottomRows(ps.xDim) };
    Q_.noalias() = tmp.transpose() * weights_.asDiagonal() * tmp;
    c_.noalias() = (M_ * (ps.Phi.bottomRows(ps.xDim) * ps.x0 + ps.xi.bottomRows(ps.xDim)) - p_).transpose() * weights_.asDiagonal() * tmp;
}

/*************************************************************************************************
 *                                   Control Cost Function                                       *
 *************************************************************************************************/

void ControlCost::autoSpan()
{
    auto max_dim = std::max(N_.rows(), std::max(weights_.rows(), p_.rows()));
    AutoSpan::spanMatrix(N_, max_dim);
    AutoSpan::spanVector(p_, max_dim);
    AutoSpan::spanVector(weights_, max_dim);
}

void ControlCost::initializeCost(const System& ps)
{
    CostFunction::initializeCost(ps);
    if (N_.rows() != p_.rows())
        DOMAIN_ERROR_EXCEPTION(throwMsgOnRowsAskAutoSpan("N", "p", N_, p_));

    if (N_.cols() == ps.uDim)
        Q_.setZero();
    else if (N_.cols() == ps.fullUDim)
        fullSizeEntry_ = true;
    else
        DOMAIN_ERROR_EXCEPTION(throwMsgOnColsOnPSUDim("N", N_, &ps));
}

void ControlCost::update(const System& ps)
{
    if (fullSizeEntry_) {
        Q_.noalias() = N_.transpose() * weights_.asDiagonal() * N_;
        c_.noalias() = -p_.transpose() * weights_.asDiagonal() * N_;
    } else {
        Eigen::MatrixXd mat{ N_.transpose() * weights_.asDiagonal() * N_ };
        Eigen::VectorXd vec{ -p_.transpose() * weights_.asDiagonal() * N_ };
        for (int i = 0; i < ps.nrUStep; ++i) {
            Q_.block(i * ps.uDim, i * ps.uDim, ps.uDim, ps.uDim) = mat;
            c_.segment(i * ps.uDim, ps.uDim) = vec;
        }
    }
}

/*************************************************************************************************
 *                            Simple Control Cost Function                                       *
 *************************************************************************************************/

void SimpleControlCost::autoSpan()
{
    auto max_dim = std::max(weights_.rows(), p_.rows());
    AutoSpan::spanVector(p_, max_dim);
    AutoSpan::spanVector(weights_, max_dim);
}

void SimpleControlCost::initializeCost(const System& ps)
{
    CostFunction::initializeCost(ps);

    if (p_.rows() == ps.uDim)
        Q_.setZero();
    else if (p_.rows() == ps.fullUDim)
        fullSizeEntry_ = true;
    else
        DOMAIN_ERROR_EXCEPTION(throwMsgOnRowsOnPSUDim("p", p_, &ps));
}

void SimpleControlCost::update(const System& ps)
{
    if (fullSizeEntry_) {
        Q_.noalias() = weights_.asDiagonal().toDenseMatrix();
        c_.noalias() = -p_.cwiseProduct(weights_);
    } else {
        Eigen::MatrixXd mat{weights_.asDiagonal().toDenseMatrix()};
        Eigen::VectorXd vec{-p_.cwiseProduct(weights_)};
        for (int i = 0; i < ps.nrUStep; ++i) {
            Q_.block(i * ps.uDim, i * ps.uDim, ps.uDim, ps.uDim) = mat;
            c_.segment(i * ps.uDim, ps.uDim) = vec;
        }
    }
}

/*************************************************************************************************
 *                                    Mixed Cost Function                                        *
 *************************************************************************************************/

void MixedCost::autoSpan()
{
    auto max_dim = std::max(M_.rows(), std::max(N_.rows(), std::max(weights_.rows(), p_.rows())));
    AutoSpan::spanMatrix(M_, max_dim, 1); // This is tricky. Has X and U are not the same dimensions, we need to handle it.
    AutoSpan::spanMatrix(N_, max_dim);
    AutoSpan::spanVector(p_, max_dim);
    AutoSpan::spanVector(weights_, max_dim);
}

void MixedCost::initializeCost(const System& ps)
{
    CostFunction::initializeCost(ps);
    if (M_.rows() != p_.rows())
        DOMAIN_ERROR_EXCEPTION(throwMsgOnRowsAskAutoSpan("M", "p", M_, p_));
    if (N_.rows() != p_.rows())
        DOMAIN_ERROR_EXCEPTION(throwMsgOnRowsAskAutoSpan("N", "p", N_, p_));

    if (M_.cols() == ps.xDim && N_.cols() == ps.uDim) {
        Q_.setZero();
        c_.setZero();
    } else if (M_.cols() == ps.fullXDim && N_.cols() == ps.fullUDim) {
        fullSizeEntry_ = true;
    } else {
        DOMAIN_ERROR_EXCEPTION(throwMsgOnColsOnPSXUDim("M", "N", M_, N_, &ps));
    }
}

void MixedCost::update(const System& ps)
{
    if (fullSizeEntry_) {
        Eigen::MatrixXd tmp = M_ * ps.Psi + N_;
        Q_.noalias() = tmp.transpose() * weights_.asDiagonal() * tmp;
        c_.noalias() = (M_ * (ps.Phi * ps.x0 + ps.xi) - p_).transpose() * weights_.asDiagonal() * tmp;
    } else {
        for (int i = 0; i < ps.nrUStep; ++i) { // Can be optimized. Lot of sums of zero here
            Eigen::MatrixXd tmp{ M_ * ps.Psi.block(i * ps.xDim, 0, ps.xDim, ps.fullUDim) };
            tmp.block(0, i * ps.uDim, M_.rows(), ps.uDim).noalias() += N_;
            Q_.noalias() += tmp.transpose() * weights_.asDiagonal() * tmp;
            c_.noalias() += (M_ * (ps.Phi.block(i * ps.xDim, 0, ps.xDim, ps.xDim) * ps.x0 + ps.xi.segment(i * ps.xDim, ps.xDim)) - p_).transpose() * weights_.asDiagonal() * tmp;
        }
    }
}

} // namespace copra
