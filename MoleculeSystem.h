/* 
 * File:   MoleculeSystem.h
 * Author: Richard Mace
 *
 * Created on 03 March 2020, 10:48
 */

#ifndef MOLECULESYSTEM_H
#define MOLECULESYSTEM_H

#include "DynamicalSystem.h"

class MoleculeSystem : public DynamicalSystem {
public:        
    struct Extents {
        double xmin;
        double xmax;
        double ymin;
        double ymax;
        double zmin;
        double zmax;
    };
    
    enum class LoadType {SOLID, LIQUID, GAS};
    
    MoleculeSystem(std::size_t, std::size_t, const Extents&, LoadType);
    ~MoleculeSystem();
    
    virtual std::vector<double> 
        stateDerivative(double, std::vector<double>&) override;
    void step(double t1, double t2) override;
    void setGravitationalAcceleration(double);
    double getGravitationalAcceleration() const;
    std::size_t getNumber() const; 
    std::size_t getDegreesOfFreedom() const;
    Extents getExtents() const;
    void setExtents(const MoleculeSystem::Extents&);
    void scaleVelocities(double factor);
    void immerseCoolBath();
    void immerseHotBath();
    
private:
    static constexpr double kKineticMax_                = 10.0;
    static constexpr double kKToverEpsilonGas_          = 2.0;
    static constexpr double kKToverEpsilonLiquid_       = 0.8;
    static constexpr double kKToverEpsilonSolid_        = 1.0e-5;
    
    std::size_t numMolecules_;        // number of molecules
    std::size_t numDegreesOfFreedom_; // number of DoF per molecule
    Extents     extents_;        // container extents in x-, y-, z-directions
    double      gravitationalAcceleration_;  // external gravitational field g
    
    // private helper functions
    void initRandom(double);
    void initLatticeHCP(double);
    void applyBoundaryConditions();
};

#endif /* MOLECULESYSTEM_H */

