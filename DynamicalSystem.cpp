/*
 * File:   DynamicalSystem.cpp
 * Author: Richard Mace
 *
 * March 2020
 *
 */

#include "DynamicalSystem.h"
#include <iostream>
#include <cmath>

// Constructor
//

DynamicalSystem::DynamicalSystem(std::size_t dimension) 
    : dimension_(dimension)
    , stepper_accuracy_(kDefaultAccuracyGoal_) {
    state_.resize(2 * dimension_);
}

// Return the dimension of the dynamical system.
//
// The dimension is equivalent to the number of degrees of freedom of the 
// system. The order of the ordinary differential equation system corresponding
// to a DynamicalSystem is 2 * dimension.
//

std::size_t DynamicalSystem::getDimension() const {
    return dimension_;
}

// Return a constant reference to the system state vector

const std::vector<double>& DynamicalSystem::state() const {
    return state_;
}

// Integrate the system differential equations from time t1 to time t2.

void DynamicalSystem::step(double t1, double t2) {
    ODEIntegrate(t1, t2, stepper_accuracy_);
}

// Set integration accuracy.

void DynamicalSystem::setAccuracy(double accuracy) {
    stepper_accuracy_ = accuracy;
}

// Get integration accuracy.

double DynamicalSystem::getAccuracy() const {
    return stepper_accuracy_;
}

///////////////////////////////////////////////////////////////////////////////
// private helper/utility member functions                                   //
///////////////////////////////////////////////////////////////////////////////

// RungeKuttaCashKarp:
//
// Given a vector currDerivative containing the time derivative of the system 
// state vector (state_) at time t, executes a low-level Runge-Kutta-Cash-Karp 
// step that reads the current member state_ (vector) of the dynamical system 
// at time t and returns in the vector argument newState the new state of the
// system at time t + deltat. The member variable state_ is NOT changed. An 
// estimate for the error in the solution is returned in argument stateError.
//
// This function uses the pure virtual member function stateDerivative(), which
// must be overriden in the derived class. The function 
// std::vector<double> stateDerivative(double t, std::vector<double>& state) 
// must return the time derivative at time t of an ARBITRARY given state vector, 
// state (at time t). The stateDerivative function effectively defines the 
// first order ordinary differential equation system governing the dynamics of
// the dynamical system (derived from abstract base class DynamicalSystem) being 
// modelled.
//
// REFERENCE: Numerical Recipes in C, p. 713.
//

void DynamicalSystem::RungeKuttaCashKarp(
                                    double t, 
                                    double deltat,
                                    const std::vector<double>& currDerivative,
                                    std::vector<double>& newState,
                                    std::vector<double>& stateError) {
    
    const double A2 = 0.2, A3 = 0.3, A4 = 0.6, A5 = 1.0,
                 A6 = 0.875, B21 = 0.2, B31 = 3.0 / 40.0,
                 B32 = 9.0 / 40.0, B41 = 0.3, B42 = -0.9,
                 B43 = 1.2, B51 = -11.0 / 54.0, B52 = 2.5,
                 B53 = -70.0 / 27.0, B54 = 35.0 / 27.0,
                 B61 = 1631.0 / 55296.0, B62 = 175.0 / 512.0,
                 B63 = 575.0 / 13824.0, B64 = 44275.0 / 110592.0,
                 B65 = 253.0 / 4096.0,
                 C1 = 37.0 / 378.0, C3 = 250.0 / 621.0,
                 C4 = 125.0 / 594.0, C6 = 512.0 / 1771.0,
                 DC1 = C1 - 2825.0 / 27648.0,
                 DC3 = C3 - 18575.0 / 48384.0,
                 DC4 = C4 - 13525.0 / 55296.0,
                 DC5 = -277.0 / 14336.0, DC6 = C6 - 0.25;
    
    std::size_t order = 2 * dimension_;
    
    // First step. We already have the derivs in currDerivative, so needn't
    // call stateDerivative here...
    std::vector<double> tempState(order);    
    for (std::size_t i = 0; i < order; ++i) {
        tempState[i] = state_[i] + B21 * deltat * currDerivative[i];
    }
    
    // Second step
    std::vector<double> ak2 = stateDerivative(t + A2 * deltat, tempState);
    for (std::size_t i = 0; i < order; ++i) {
        tempState[i] = state_[i] + deltat * (B31 * currDerivative[i] 
                + B32 * ak2[i]);
    }

    // Third step
    std::vector<double> ak3 = stateDerivative(t + A3 * deltat, tempState);
    for (std::size_t i = 0; i < order; ++i) {
        tempState[i] = state_[i] + deltat * (B41 * currDerivative[i] 
                + B42 * ak2[i] + B43 * ak3[i]);
    }

    // Fourth step
    std::vector<double> ak4 = stateDerivative(t + A4 * deltat, tempState);
    for (std::size_t i = 0; i < order; ++i) {
        tempState[i] = state_[i] + deltat * (B51 * currDerivative[i] 
                + B52 * ak2[i] + B53 * ak3[i] + B54 * ak4[i]);
    }

    // Fifth step
    std::vector<double> ak5 = stateDerivative(t + A5 * deltat, tempState);
    for (std::size_t i = 0; i < order; ++i) {
        tempState[i] = state_[i] + deltat * (B61 * currDerivative[i] 
                + B62 * ak2[i] + B63 * ak3[i] + B64 * ak4[i] + B65 * ak5[i]);
    }

    // Evaluate the new state
    std::vector<double> ak6 = stateDerivative(t + A6 * deltat, tempState);
    for (std::size_t i = 0; i < order; ++i) {
        newState[i] = state_[i] + deltat * (C1 * currDerivative[i] + C3 * ak3[i] 
                + C4 * ak4[i] + C6 * ak6[i]);
    }

    // The error is estimated as the difference between the
    // fourth and fifth order Runge-Kutta methods
    for (std::size_t i = 0; i < order; ++i) {
        stateError[i] = deltat * (DC1 * currDerivative[i] + DC3 * ak3[i] 
                + DC4 * ak4[i] + DC5 * ak5[i] + DC6 * ak6[i]);
    }
}


// Quality controlled Runge-Kutta-Cash-Karp step.
//
// Used by ODEIntegrate to integrate the ODE system defined in virtual member
// function stateDerivative(...) from time t to t + Delta, where Delta is 
// determined by accuracy and convergence conditions. The actual value of Delta
// used is returned in variable deltaDid, and a guess for Delta to be used in 
// the next integration step is returned in deltaNext. 
//
// This function updates the member variable state_, representing the current
// system state vector.

void DynamicalSystem::RungeKuttaQualityControlled(
                                    double t, 
                                    double deltaTry,
                                    double epsAcc, 
                                    const std::vector<double>& stateDerivative,
                                    const std::vector<double>& stateScale,
                                    double& deltaDid,
                                    double& deltaNext) {
    
    const double k_safety     = 0.9;
    const double k_pow_grow   = -0.2;
    const double k_pow_shrink = -0.25;
    const double k_error_cond = 1.89e-4;

    std::size_t order = 2 * dimension_;
    
    std::vector<double> stateError(order);
    std::vector<double> newState(order);
    
    double delta = deltaTry;
    double errorMax = 0.0;
    
    do {
        RungeKuttaCashKarp(t, delta, stateDerivative, newState, stateError);

        // Evaluate the error
        errorMax = 0.0;
        for (std::size_t i = 0; i < order; ++i) {
            errorMax = std::max(errorMax, fabs(stateError[i] / stateScale[i]));
        }

        errorMax /= epsAcc;

        if (errorMax > 1.0) {
            double deltaTemp = k_safety * delta * pow(errorMax, k_pow_shrink);
            delta = copysign(std::max(fabs(deltaTemp), 0.1 * fabs(delta)), delta);
            if (t + delta == t) {
                std::cerr <<"RungeKuttaQualityControlled: stepsize underflow!"
                          << '\n';
            }
        }
    } while (errorMax > 1.0);

    // We have a solution of the required accuracy. Set next stepsize.
    if (errorMax > k_error_cond) {
        deltaNext = k_safety * delta * pow(errorMax, k_pow_grow);
    } 
    else {
        deltaNext = 0.5 * delta;
    }

    // Transfer the new state vector to the member state_.
    deltaDid = delta;
    for (std::size_t i = 0; i < order; ++i) {
        state_[i] = newState[i];
    }
}


// High-level adaptive ode integrator. Integrates the ODE system defined in 
// virtual member function stateDerivative(...) from time t1 to time t2, where 
// usually t2 > t1. The member variable state_ representing the current system 
// state vector is updated by this function (through calls to
// RungeKuttaQualityControlled(...)).

void DynamicalSystem::ODEIntegrate(double t1, double t2, double accuracy) {

    const double k_tiny = 1.0e-30;
    
    std::size_t order = 2 * dimension_;

    double delta = 0.1;
    double t = t1;
    double deltaDid;
    double deltaNext;
    
    std::vector<double> stateDot(order);
    std::vector<double> stateScale(order);
    
    do {
        stateDot = stateDerivative(t, state_);
        
        for (std::size_t i = 0; i < order; ++i) {
            stateScale[i] = fabs(state_[i]) + fabs(delta * stateDot[i]) + k_tiny;
        }
        
        RungeKuttaQualityControlled(t, delta, accuracy, stateDot, stateScale, 
                deltaDid, deltaNext);
        
        t += deltaDid;
        delta = deltaNext;

        // make sure we don't overshoot
        if (t + delta > t2) {
            delta = t2 - t;
        }
    } while (t < t2);
}