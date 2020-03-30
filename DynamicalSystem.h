/* 
 * File:   DynamicalSystem.h
 * Author: Richard Mace
 *
 * Created on 03 March 2020, 10:40
 */

#ifndef DYNAMICALSYSTEM_H
#define DYNAMICALSYSTEM_H

#include <cstddef>
#include <vector>

class DynamicalSystem {
public:
    explicit DynamicalSystem(std::size_t);
    virtual ~DynamicalSystem() = default;
    
    std::size_t getDimension() const;
    const std::vector<double>& state() const;
    virtual std::vector<double> 
        stateDerivative(double, std::vector<double>&) = 0;
    virtual void step(double, double);
    void setAccuracy(double);
    double getAccuracy() const;
    
protected:
    static constexpr double kDefaultAccuracyGoal_ = 1.0e-8; // integration 
                                                            // accuracy goal
    
    std::size_t         dimension_; // dimension (deg. of freedom) of system
    std::vector<double> state_;     // state vector (q_1,..,q_n,p_1,..,p_n)
    double              stepper_accuracy_; // integration accuracy
    
    // helper functions for system ODE integration
    void RungeKuttaCashKarp(double t, double deltat, 
                            const std::vector<double>& currDerivative, 
                            std::vector<double>& newState, 
                            std::vector<double>& stateError);
    
    void RungeKuttaQualityControlled(double t, 
                                     double deltaTry,
                                     double epsAcc, 
                                     const std::vector<double>& stateDerivative,
                                     const std::vector<double>& stateScale,
                                     double& deltaDid,
                                     double& deltaNext);
    
    void ODEIntegrate(double t1, double t2, double accuracy);
};

#endif /* DYNAMICALSYSTEM_H */

