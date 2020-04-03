/*
 * File: MoleculeSystem.cpp
 * Author: Richard Mace
 *
 * March 2020
 */

#include "MoleculeSystem.h"
#include <cmath>
#include <random>
#include <iostream>

/*
 * NOTE: With dimension_ = N * numDegreesOfFreedom_, where N is the number of 
 * molecules, and numDegreesOfFreedom_ = 3 per molecule, the state_ vector 
 * layout is as follows:
 *
 * state_[0]                  = x_1,  
 * state_[1]                  = y_1,  
 * state_[2]                  = z_1,
 * .
 * .
 * state_[dimension_ - 3]     = x_N    [N.B. dimension_ - 3 = 3(N - 1)]
 * state_[dimension_ - 2]     = y_N
 * state_[dimension_ - 1]     = z_N
 * state_[dimension_ + 0]     = vx_1, 
 * state_[dimension_ + 1]     = vy_1, 
 * state_[dimension_ + 2]     = vz_1,
 * .
 * .
 * state_[2 * dimension_ - 3] = vx_N   [N.B. 2 * dimension_ - 3 = 3N + 3(N - 1)]
 * state_[2 * dimension_ - 2] = vy_N
 * state_[2 * dimension_ - 1] = vz_N
 * 
 * In other words, indices [0 ... 3N - 1] contain particle position vector
 * coordinates, while indices [3N ... 6N - 1] contain particle velocity vector
 * coordinates. By labelling particles with index i running from 0 to N - 1, 
 * we can find the position vector (xi, yi, zi) of particle i via the following:
 * 
 *   offset_i = numDegreesOfFreedom_ * i;
 *   xi = state_[offset_i]; 
 *   yi = state_[offset_i+1]; 
 *   zi = state_[offset_i+2]; 
 *
 * The velocity of particle i, (vxi, vyi, vzi) is located via (recall that 
 * dimension_ = N * numDegreesOfFreedom_):
 * 
 *   offset_vi = dimension_ + numDegreesOfFreedom_ * i;
 *   vxi = state_[offset_vi]; 
 *   vyi = state_[offset_vi+1]; 
 *   vzi = state_[offset_vi+2];
 * 
 * In this code numDegreesOfFreedom_ = 3.
 * 
 */

///////////////////////////////////////////////////////////////////////////////
// non-member utility functions                                              //
///////////////////////////////////////////////////////////////////////////////

// gaussian_deviate()
//
// Return a Gaussian deviate with zero mean and unit standard deviation.
//

static double gaussian_deviate() {
    static std::random_device random; 
    static std::mt19937 mersenne_twister(random()); // seed twister engine with
                                                    // random() 
    static std::normal_distribution<double> gaussian(0.0, 1.0);
    
    return gaussian(mersenne_twister);
}

// uniform deviate()
//
// Return a uniform deviate in the range [0.0,1.0).
//
static double uniform_deviate() {
    static std::random_device random; 
    static std::mt19937 mersenne_twister(random()); // seed twister engine with
                                                    // random()
    
    static std::uniform_real_distribution<double> uniform(0.0, 1.0);
    
    return uniform(mersenne_twister);

}

///////////////////////////////////////////////////////////////////////////////
// Member functions                                                          //
///////////////////////////////////////////////////////////////////////////////

MoleculeSystem::MoleculeSystem(std::size_t numMolecules, std::size_t numDOF, 
        const Extents& extents, LoadType loadType) 

    : DynamicalSystem(numMolecules * numDOF)
    , numMolecules_(numMolecules)
    , numDegreesOfFreedom_(numDOF)
    , extents_(extents)
    , gravitationalAcceleration_(0.0) {
        
    if (loadType == LoadType::SOLID) {
        initLatticeHCP(kKToverEpsilonSolid_);
    }
    else if (loadType == LoadType::LIQUID) {
        initLatticeHCP(kKToverEpsilonLiquid_);
    }
    else if (loadType == LoadType::GAS) {
        initRandom(kKToverEpsilonGas_);
    }
    else {
        std::cerr << "Invalid load type. Defaulting to SOLID." << '\n';
        initLatticeHCP(kKToverEpsilonSolid_);
    }
}


MoleculeSystem::~MoleculeSystem() {
    // empty
    std::cout << "MoleculeSystem: Destructor." << '\n';
}

// Return the number of molecules.

std::size_t MoleculeSystem::getNumber() const {
    return numMolecules_;
}

// Return the number of degrees of freedom per molecule.

std::size_t MoleculeSystem::getDegreesOfFreedom() const {
    return numDegreesOfFreedom_;
}

// Return the container extents.

MoleculeSystem::Extents MoleculeSystem::getExtents() const {
    return extents_;
}

// Set the container extents.

void MoleculeSystem::setExtents(const MoleculeSystem::Extents& extents) {
    extents_ = extents;
}

// Set the magnitude of the gravitational acceleration.

void MoleculeSystem::setGravitationalAcceleration(double g) {
    gravitationalAcceleration_ = fabs(g);
}

// Get the magnitude of the gravitational acceleration.

double MoleculeSystem::getGravitationalAcceleration() const {
    return gravitationalAcceleration_;
}

// Step the system in time from t1 to t2.

void MoleculeSystem::step(double t1, double t2) {
    DynamicalSystem::step(t1, t2);
    //DynamicalSystem::EulerStep(t1, t2);
    applyBoundaryConditions();
}

// pairInteraction(...)
//
// Non-member function that calculates the force, force_ij, on particle i due 
// to particle j.
//
// Particle i has position vector (xi, yi, zi) = (pos_i[0], pos_i[1], pos_i[2]) 
// and particle j has position vector (xj, yj, zj) = (pos_j[0], pos_j[1], 
// pos_j[2]).
// 
// The force on particle i due to j is given by the derivative of the 
// interaction potential function. In this case, the potential function 
// governing the internal forces between pairs of atoms (point molecules) is the 
// Lennard-Jones 6-12 potential written in symmetric form as follows:
//
//     V(r) = 4.0 * epsilon [(a/r)^12 - (a/r)^6]
//
// where r is the distance between the pair of atoms and a is their diameter. 
// The parameter epsilon is the depth of the potential well.
//

static void pairInteraction(double* const force_ij, const double* const pos_i, 
                     const double* const pos_j) {
    double rjix = pos_i[0] - pos_j[0]; 
    double rjiy = pos_i[1] - pos_j[1]; 
    double rjiz = pos_i[2] - pos_j[2];
    
    double r2 = rjix * rjix + rjiy * rjiy + rjiz * rjiz; 
    double r6 = r2 * r2 * r2;
    double r8 = r6 * r2;
                
    // Compute the force on particle i due to particle j. This comes from the 
    // derivative of the Lennard-Jones potential, Eq. (3.5) of Flowers and 
    // Mendoza, in normalised form.
    double flaw = -(1.0 - 2.0 / r6) / r8;
    double fjix = rjix * flaw;
    double fjiy = rjiy * flaw;
    double fjiz = rjiz * flaw;
                 
    // Add it to the running total for force on particle i (due to particle j).
    force_ij[0] += fjix;
    force_ij[1] += fjiy;
    force_ij[2] += fjiz;
}

// stateDerivative:
//
// Evaluate the time derivative of the given state vector thisState at given 
// time t (in units of 1 / omega_E -- see below). 
//
// The interaction potential governing the internal forces between pairs of 
// atoms (point molecules) is the Lennard-Jones 6-12 potential (see function
// pairInteraction above).
//
// The equations of motion (Newton's equations) have been written in 
// dimensionless units as follows. Lengths have been normalised by a, the 
// diameter of an atom (point molecule). Times have been normalised by 
// 1 / omega_E, where omega_E is the Einstein angular frequency. The Einstein
// (angular) frequency is a measure of the frequency of oscillation of atoms 
// trapped in the potential well and is defined by (c.f. Flowers and Mendoza, 
// Properties of Matter, p. 46, Eq. (3.20)):
//
//            omega_E = sqrt[24 * epsilon / (m * a^2)]
// 
// where epsilon is the depth of the Lennard-Jones potential (see above), m is 
// the mass of an atom (molecule), and a is its diameter. 
//
// With the layout of the state vector (state_) described above, the elements 
// of its time derivative, stateDot, have the following meanings:
// stateDot[0] ... stateDot[dimension_ - 1] correspond to velocity components, 
// whereas stateDot[dimension_] ... stateDot[2 * dimension_ - 1] correspond to 
// acceleration components. In other words, the elements stateDot[dimension_] 
// ... stateDot[2 * dimension_ - 1] are determined via Newton's equations for
// the system of particles.
// 

std::vector<double> MoleculeSystem::stateDerivative(
                                            double t, 
                                            std::vector<double>& thisState) {
    // Components of total (sum of) external force(s).
    double f_external_x = 0.0;
    double f_external_y = 0.0;
    double f_external_z = -1.0 * gravitationalAcceleration_;
    
    std::vector<double> stateDot(2 * dimension_);
    
    for (std::size_t i = 0; i < numMolecules_; ++i) {
        std::size_t pos_offset_i = i * numDegreesOfFreedom_;
        
        // vxi, vyi, vzi
        stateDot[pos_offset_i + 0] = thisState[pos_offset_i + 0 + dimension_]; 
        stateDot[pos_offset_i + 1] = thisState[pos_offset_i + 1 + dimension_]; 
        stateDot[pos_offset_i + 2] = thisState[pos_offset_i + 2 + dimension_];
        
        // Calculate the acceleration.
        double force_i[3] = {0.0, 0.0, 0.0};
        
        for (std::size_t j = 0; j < numMolecules_; ++j) {
            if (j != i) {  
                std::size_t pos_offset_j = j * numDegreesOfFreedom_;
                pairInteraction(force_i, &thisState[pos_offset_i], 
                        &thisState[pos_offset_j]);
            }
        } 

        // axi, ayi, azi
        stateDot[pos_offset_i + dimension_ + 0] = force_i[0] + f_external_x;
        stateDot[pos_offset_i + dimension_ + 1] = force_i[1] + f_external_y;
        stateDot[pos_offset_i + dimension_ + 2] = force_i[2] + f_external_z;
    }
    
    return stateDot;
}

///////////////////////////////////////////////////////////////////////////////
// private member functions                                                  //
///////////////////////////////////////////////////////////////////////////////

void MoleculeSystem::initRandom(double thermalEnergy) {
    double lx = extents_.xmax - extents_.xmin;
    double ly = extents_.ymax - extents_.ymin;
    double lz = extents_.zmax - extents_.zmin;
    
    // Initialise thermal speed such that the ratio of the thermal energy, kT, 
    // to the depth of the Lennard-Jones potential well, epsilon, is given by
    // thermalEnergy.
    double vtherm = sqrt(thermalEnergy / 24.0);
    
    double minx = extents_.xmin;
    double miny = extents_.ymin;
    double minz = extents_.zmin;
    
    double* ptrMolecule = state_.data();
    
    for (std::size_t i = 0; i < numMolecules_; ++i) {
        // Set initial x, y and z coords of each molecule randomly inside  box.
        ptrMolecule[0] = minx + uniform_deviate() * lx;
        ptrMolecule[1] = miny + uniform_deviate() * ly;
        ptrMolecule[2] = minz + uniform_deviate() * lz;
        
        // Set initial vx, vy and vz of each molecule using Gaussian deviates 
        // with zero mean and standard deviation equal to the (normalised) 
        // thermal speed. Gives an approximately Maxwellian distribution initial 
        // loading, with thermal energy kT (in units of epsilon) given as above.
        ptrMolecule[dimension_ + 0] = vtherm * gaussian_deviate();
        ptrMolecule[dimension_ + 1] = vtherm * gaussian_deviate();
        ptrMolecule[dimension_ + 2] = vtherm * gaussian_deviate();
        
        // next molecule address
        ptrMolecule += numDegreesOfFreedom_;
    } 
}

// initLatticeHCP():
//
// Load particles on a simple hexagonal close packed (HCP) lattice. 
//
// Algorithm is a bit crude at the moment. Unless numMolecules_ is a "nice" 
// number, the last crystal face may be incomplete. Needs improvement.
//

void MoleculeSystem::initLatticeHCP(double thermalEnergy) {
    double xmin = std::numeric_limits<double>::max();
    double xmax = std::numeric_limits<double>::lowest();
    double ymin = xmin;
    double ymax = xmax;
    double zmin = xmin;
    double zmax = xmax;
    
    double d = 1.12246204831; // 2^(1/6) = equilibrium intermolecular spacing
    double root2 = sqrt(2.0);
    double root3 = sqrt(3.0);
    
    int numx = static_cast<int>(
                   ceil(pow(static_cast<double>(numMolecules_), 1.0/3.0))
                            );
    int numy = numx;
    
    // Calculate the thermal speed of a molecule given the ratio of 
    // thermal energy to depth of Lennard-Jones potential well (epsilon), i.e.,
    // the value of the thermalEnergy argument.
    double vtherm = sqrt(thermalEnergy / 24.0);
    
    for (std::size_t nptcl = 0; nptcl < numMolecules_; ++nptcl) {
        int k = static_cast<int>(nptcl) / (numx * numy);
        int nplane = static_cast<int>(nptcl) % (numx * numy);
       
        int j = nplane / numx;
        int i = nplane % numx;
       
        double x = i * d + j * (d / 2.0) + k * (d / 2.0);
        double y = j * (root3 * d / 2.0) + k * (d / (2.0 * root3));
        double z = k * (root2 * d / root3);
        
        std::size_t offset = nptcl * numDegreesOfFreedom_;
                
        state_[offset + 0] = x;
        state_[offset + 1] = y;
        state_[offset + 2] = z;
                
        if (x < xmin) {
            xmin = x;
        }
        if (x > xmax) {
            xmax = x;
        }
        if (y < ymin) { 
            ymin = y;
        }
        if (y > ymax) {
            ymax = y;
        }
        if (z < zmin) {
            zmin = z;
        }
        if (z > zmax) {
            zmax = z;
        }
        
        // vx, vy and vz are given by deviates from a gaussian (normal) 
        // distribution with zero mean and unit standard deviation, multiplied
        // by the thermal speed.
        state_[offset + 0 + dimension_] = vtherm * gaussian_deviate();
        state_[offset + 1 + dimension_] = vtherm * gaussian_deviate();
        state_[offset + 2 + dimension_] = vtherm * gaussian_deviate();
    }

    // Centre lattice at the origin of our fixed coordinate system.
    
    for (std::size_t i = 0; i < numMolecules_; ++i) {
        std::size_t offset = i * numDegreesOfFreedom_;
        state_[offset + 0] -= (xmax - xmin) / 2.0;
        state_[offset + 1] -= (ymax - ymin) / 2.0;
        state_[offset + 2] -= (zmax - zmin) / 2.0;
    }
}

/*
 * applyBoundaryConditions:
 *
 * Applies simple specular reflection of molecules at the rigid container 
 * walls. A bit crude at the moment. Could be more careful to conserve energy.
 *
 */

void MoleculeSystem::applyBoundaryConditions() {
    double minx = extents_.xmin;
    double maxx = extents_.xmax;
    double miny = extents_.ymin;
    double maxy = extents_.ymax;
    double minz = extents_.zmin;
    double maxz = extents_.zmax;
    
    for (std::size_t i = 0; i < numMolecules_; ++i) {
        std::size_t n = i * numDegreesOfFreedom_;
        double x = state_[n + 0];
        double y = state_[n + 1];
        double z = state_[n + 2];
        
        double vx = state_[n + 0 + dimension_];
        double vy = state_[n + 1 + dimension_];
        double vz = state_[n + 2 + dimension_];
        
        // Left hand side.
        if (x < minx) {
            x = minx + (minx - x);
            vx = -vx;
        }
        
        // Right hand side.
        if (x > maxx) {
            x = maxx - (x - maxx);
            vx = -vx;
        }
        
        // Bottom. 
        if (y < miny) {
            y = miny + (miny - y);
            vy = -vy;
        }
        
        // Top.
        if (y > maxy) {
            y = maxy - (y - maxy);
            vy = -vy;
        }
        
        // Near. 
        if (z < minz) {
            z = minz + (minz - z);
            vz = -vz; 
        }
        
        // Far.
        if (z > maxz) {
            z = maxz - (z - maxz);
            vz = -vz; 
        }
        
        // Restore values to the state vector.
        state_[n + 0] = x;
        state_[n + 1] = y;
        state_[n + 2] = z;
        state_[n + 0 + dimension_] = vx;
        state_[n + 1 + dimension_] = vy;
        state_[n + 2 + dimension_] = vz;
    }
}

/*
 * scaleVelocities:
 *
 * Scales the velocities (momenta) by multiplicative factor, factor. Used to
 * simulate heating or cooling of the molecule system.
 */

void MoleculeSystem::scaleVelocities(double factor) {
    for (std::size_t i = 0; i < numMolecules_; ++i) {
        std::size_t offset = i * numDegreesOfFreedom_ + dimension_;
        double vx = state_[offset + 0];
        double vy = state_[offset + 1];
        double vz = state_[offset + 2];
        
        // In our chosen units, the kinetic energy (in units of epsilon -- the 
        // depth of the Lennard-Jones potential) is equal to 12 * v^2, where
        // v is the speed in normalised units.
        double kineticEnergy = 12 * (vx * vx + vy * vy + vz * vz);
        
        if (factor < 1.0 || kineticEnergy < kKineticMax_) {
            state_[offset + 0] *= factor;
            state_[offset + 1] *= factor;
            state_[offset + 2] *= factor;
        }
    }
}

// immserseCoolBath()
//
// Simulate the effects of cooling the system by immersion in a bath (solvent)
// of cooler molecules. Very crude at the moment, but effective.
//
// The thermal energy determining the thermal speed used to perturb the initial
// particle velocity is (arbitrarily) 120.0 * (KT/epsilon)_solid.

void MoleculeSystem::immerseCoolBath() {
    constexpr double vtherm = sqrt(120.0 * kKToverEpsilonSolid_ / 24.0);
    
    for (std::size_t i = 0; i < numMolecules_; ++i) {
        std::size_t offset = i * numDegreesOfFreedom_ + dimension_;
        double vx = state_[offset + 0];
        double vy = state_[offset + 1];
        double vz = state_[offset + 2];
        
        state_[offset + 0] = 0.98 * vx + vtherm * gaussian_deviate();
        state_[offset + 1] = 0.98 * vy + vtherm * gaussian_deviate();
        state_[offset + 2] = 0.98 * vz + vtherm * gaussian_deviate();
    }
}

// immerseHotBath()
//
// Simulate the effects of heating the system by immersion in a bath (solvent)
// of hotter molecules. Very crude at the moment, but effective.
//
// The thermal energy determining the thermal speed used to perturb the initial
// particle velocity is (arbitrarily) 2000.0 * (KT/epsilon)_solid.

void MoleculeSystem::immerseHotBath() {
    constexpr double vtherm = sqrt(2000.0 * kKToverEpsilonSolid_ / 24.0);
    
    for (std::size_t i = 0; i < numMolecules_; ++i) {
        std::size_t offset = i * numDegreesOfFreedom_ + dimension_;
        double vx = state_[offset + 0];
        double vy = state_[offset + 1];
        double vz = state_[offset + 2];
        
        // In our chosen units, the kinetic energy (in units of epsilon -- the 
        // depth of the Lennard-Jones potential) is equal to 12 * v^2, where
        // v is the speed in normalised units.
        double kineticEnergy = 12 * (vx * vx + vy * vy + vz * vz);
        
        if (kineticEnergy < kKineticMax_) {
            state_[offset + 0] += vtherm * gaussian_deviate();
            state_[offset + 1] += vtherm * gaussian_deviate();
            state_[offset + 2] += vtherm * gaussian_deviate();
        }
    }
}