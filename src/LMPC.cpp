/*
 * Copyright 2016-2019 CNRS-UM LIRMM, CNRS-AIST JRL
 * Copyright 2020 ANYbotics AG
 */

#include <algorithm>
#include <exception>
#include <numeric>

#include "copra/LMPC.h"
#include "copra/constraints.h"
#include "copra/costFunctions.h"
#include "copra/systems/System.h"

namespace copra {

LMPC::Constraints::Constraints()
    : nrEqConstr(0)
    , nrIneqConstr(0)
    , spConstr()
    , spEqConstr()
    , spIneqConstr()
    , spBoundConstr()
{
}

void LMPC::Constraints::updateNr()
{
    nrEqConstr = 0;
    nrIneqConstr = 0;

    auto countConstr = [](auto& spc, int& nreq) {
        for (auto& sp : spc)
            nreq += sp->nrConstr();
    };

    countConstr(spEqConstr, nrEqConstr);
    countConstr(spIneqConstr, nrIneqConstr);
}

void LMPC::Constraints::clear()
{
    nrEqConstr = 0;
    nrIneqConstr = 0;
    spConstr.clear();
    spEqConstr.clear();
    spIneqConstr.clear();
    spBoundConstr.clear();
}

/*************************************************************************************************
 *                                             LMPC                                               *
 *************************************************************************************************/

LMPC::LMPC(SolverFlag sFlag)
    : system_(nullptr)
    ,
      solver_(solverFactory(sFlag))
    , constraints_()
    , Q_()
    , Aineq_()
    , Aeq_()
    , c_()
    , bineq_()
    , beq_()
    , lb_()
    , ub_()
    , solveTime_()
    , solveAndBuildTime_()
{
}

LMPC::LMPC(const std::shared_ptr<System>& ps, SolverFlag sFlag)
    : system_(ps)
    ,
      solver_(solverFactory(sFlag))
    , constraints_()
    , Q_(system_->fullUDim, system_->fullUDim)
    , Aineq_(0, system_->fullUDim)
    , Aeq_(0, system_->fullUDim)
    , c_(system_->fullUDim)
    , bineq_()
    , beq_()
    , lb_(system_->fullUDim)
    , ub_(system_->fullUDim)
    , solveTime_()
    , solveAndBuildTime_()
{
    lb_.setConstant(-std::numeric_limits<double>::max());
    ub_.setConstant(std::numeric_limits<double>::max());
}

void LMPC::selectQPSolver(SolverFlag flag)
{
  solver_ = solverFactory(flag);
}

void LMPC::useSolver(std::unique_ptr<SolverInterface>&& solver)
{
  solver_ = std::move(solver);
}

void LMPC::initializeController(const std::shared_ptr<System>& ps)
{
  system_ = ps;
    clearConstraintMatrices();

    Q_.resize(system_->fullUDim, system_->fullUDim);
    c_.resize(system_->fullUDim);
}

bool LMPC::solve()
{
    using namespace std::chrono;
    auto sabTime = high_resolution_clock::now();

    updateSystem();
    makeQPForm();
    solver_->SI_problem(system_->fullUDim, constraints_.nrEqConstr, constraints_.nrIneqConstr);
    auto sTime = high_resolution_clock::now();
    bool success = solver_->SI_solve(Q_, c_, Aeq_, beq_, Aineq_, bineq_, lb_, ub_);
    solveTime_ = duration_cast<duration<double>>(high_resolution_clock::now() - sTime);

    checkDeleteCostsAndConstraints();

    solveAndBuildTime_ = duration_cast<duration<double>>(high_resolution_clock::now() - sabTime);
    return success;
}

void LMPC::inform() const noexcept
{
    solver_->SI_inform();
}

const Eigen::VectorXd& LMPC::control() const noexcept
{
    return solver_->SI_result();
}

Eigen::VectorXd LMPC::trajectory() const noexcept
{
    return system_->Phi * system_->x0 + system_->Psi * control() + system_->xi;
}

double LMPC::solveTime() const noexcept
{
    return solveTime_.count();
}

double LMPC::solveAndBuildTime() const noexcept
{
    return solveAndBuildTime_.count();
}

void LMPC::addCost(const std::shared_ptr<CostFunction>& costFun)
{
    costFun->initializeCost(*system_);
    spCost_.emplace_back(costFun);
}

void LMPC::addConstraint(const std::shared_ptr<Constraint>& constr)
{
    constr->initializeConstraint(*system_);
    addConstraintByType(constr);
}

void LMPC::resetConstraints() noexcept
{
    constraints_.clear();
}

void LMPC::resetCosts() noexcept {
  spCost_.clear();
}

void LMPC::removeCost(const std::shared_ptr<CostFunction>& costFun)
{
    auto rCost = std::find(spCost_.begin(), spCost_.end(), costFun);
    if (rCost != spCost_.end())
        spCost_.erase(rCost);
}

void LMPC::removeConstraint(const std::shared_ptr<Constraint>& constr)
{
    auto rConstr = std::find(constraints_.spConstr.begin(), constraints_.spConstr.end(), constr);
    if (rConstr != constraints_.spConstr.end()) {
        constraints_.spConstr.erase(rConstr);
        switch (constr->constraintType()) {
        case ConstraintFlag::Constraint:
        case ConstraintFlag::EqualityConstraint:
            constraints_.spEqConstr.erase(std::find(constraints_.spEqConstr.begin(), constraints_.spEqConstr.end(), constr));
            break;
        case ConstraintFlag::InequalityConstraint:
            constraints_.spIneqConstr.erase(std::find(constraints_.spIneqConstr.begin(), constraints_.spIneqConstr.end(), constr));
            break;
        case ConstraintFlag::BoundConstraint:
            constraints_.spBoundConstr.erase(std::find(constraints_.spBoundConstr.begin(), constraints_.spBoundConstr.end(), constr));
            break;
        }
    }
}

/*
 *  Protected methods
 */

void LMPC::addConstraintByType(const std::shared_ptr<Constraint>& constr)
{
    switch (constr->constraintType()) {
    case ConstraintFlag::EqualityConstraint: {
        constraints_.nrEqConstr += constr->nrConstr();
        // DownCasting to std::shared_ptr<EqIneqConstraint>
        // This is a safe operation since we know that the object is a derived class of a EqIneqConstraint
        constraints_.spEqConstr.emplace_back(std::static_pointer_cast<EqIneqConstraint>(constr));
    } break;
    case ConstraintFlag::InequalityConstraint: {
        constraints_.nrIneqConstr += constr->nrConstr();
        // DownCasting to std::shared_ptr<EqIneqConstraint>
        // This is a safe operation since we know that the object is a derived class of a EqIneqConstraint
        constraints_.spIneqConstr.emplace_back(std::static_pointer_cast<EqIneqConstraint>(constr));
    } break;
    case ConstraintFlag::BoundConstraint: {
        // DownCasting to std::shared_ptr<ControlBoundConstraint>
        // This is a safe operation since we know that the object is a ControlBoundConstraint
        constraints_.spBoundConstr.emplace_back(std::static_pointer_cast<ControlBoundConstraint>(constr));
    } break;
    default:
        return;
    }
    constraints_.spConstr.emplace_back(constr);
}

void LMPC::clearConstraintMatrices()
{
    assert(system_ != nullptr);

    Aineq_.resize(0, system_->fullUDim);
    Aeq_.resize(0, system_->fullUDim);
    bineq_.resize(0);
    beq_.resize(0);
    lb_.resize(system_->fullUDim);
    ub_.resize(system_->fullUDim);
    lb_.setConstant(-std::numeric_limits<double>::max());
    ub_.setConstant(std::numeric_limits<double>::max());
}

void LMPC::updateSystem()
{
    // Reset the QP variables
    Q_.setIdentity();
    Q_ *= 1e-6; // Ensure that Q is positive definite if no cost functions are used
    c_.setZero();

    // Update the system
    if (!system_->isUpdated) system_->update();

    // Update the constraints
    constraints_.updateNr();
    Aeq_.resize(constraints_.nrEqConstr, system_->fullUDim);
    beq_.resize(constraints_.nrEqConstr);
    Aineq_.resize(constraints_.nrIneqConstr, system_->fullUDim);
    bineq_.resize(constraints_.nrIneqConstr);

    for (auto& cstr : constraints_.spConstr)
        cstr->update(*system_);

    // Update the costs
    for (auto& sp : spCost_)
        sp->update(*system_);
}

void LMPC::makeQPForm()
{
    for (auto& cost : spCost_) {
        Q_ += cost->Q();
        c_ += cost->c();
    }

    int nrLines = 0;
    // Get Equality constraints
    for (auto& cstr : constraints_.spEqConstr) {
        Aeq_.block(nrLines, 0, cstr->nrConstr(), system_->fullUDim) = cstr->A();
        beq_.segment(nrLines, cstr->nrConstr()) = cstr->b();
        nrLines += cstr->nrConstr();
    }

    // Get Inequality constraints
    nrLines = 0;
    for (auto& cstr : constraints_.spIneqConstr) {
        Aineq_.block(nrLines, 0, cstr->nrConstr(), system_->fullUDim) = cstr->A();
        bineq_.segment(nrLines, cstr->nrConstr()) = cstr->b();
        nrLines += cstr->nrConstr();
    }

    // Get Bound constraints
    nrLines = 0;
    for (auto& cstr : constraints_.spBoundConstr) {
        lb_.segment(nrLines, cstr->nrConstr()) = cstr->lower();
        ub_.segment(nrLines, cstr->nrConstr()) = cstr->upper();
        nrLines += cstr->nrConstr();
    }
}

void LMPC::checkDeleteCostsAndConstraints()
{
    auto check = [](auto& sp, bool useWarn = false, int delLimit = 2) {
        for (auto it = sp.begin(); it != sp.end();) {
            if ((*it).use_count() <= delLimit) {
                CONSTRAINT_DELETION_WARN(useWarn, "%s%s%s", "A '", (*it)->name().c_str(),
                    "' has been destroyed.\nIt has been removed from the controller\n");
                it = sp.erase(it);
            } else {
                ++it;
            }
        }
    };

    check(constraints_.spConstr, true);
    check(constraints_.spEqConstr);
    check(constraints_.spIneqConstr);
    check(constraints_.spBoundConstr);
    check(spCost_, true, 1);
}

} // namespace copra
