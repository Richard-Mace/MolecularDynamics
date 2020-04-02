/*
 * File: MoleculeAnimation.cpp
 * Author: Richard Mace
 *
 * March 2020
 */

#include "MoleculeAnimation.h"
#include <iostream>
#include <stdexcept>


/*
 * MoleculeAnimation:
 *
 * Constructor for the MoleculeAnimation class. 
 *
 * NOTE: the deallocation of the MoleculeSystem, MoleculeVisualisation and
 *        GraphicsWindow objects dynamically allocated here is handled (smartly)
 *        by the Animation class destructor.
 */

MoleculeAnimation::MoleculeAnimation(std::size_t numMolecules, 
                                     const MoleculeSystem::Extents& extents,
                                     MoleculeSystem::LoadType loadType,
                                     GraphicsWindow::Type windowType) 

    : Animation(new MoleculeSystem(numMolecules, kDoF_, extents, loadType),
                nullptr, 
                new GraphicsWindow) {
    
    if (!ptrWindow_->isInitialised()) {
        throw std::runtime_error{"MoleculeAnimation: could not init window"};
    }
    
    ptrWindow_->setTitle("Molecular Dynamics Animation");
    
    bool success = true;
    if (windowType == GraphicsWindow::Type::WINDOWED) {
        success = ptrWindow_->open(kDefaultWidth_, kDefaultHeight_);
    }
    else {
        success = ptrWindow_->open();
    }
    if (!success) {
        throw std::runtime_error{"MoleculeAnimation: could not open window"};
    }
    
    ptrSystem_->setAccuracy(kDefaultAccuracyGoal_);
    
    // Create a molecule visualisation for the system. This must be done
    // only after a GraphicsWindow has been successfully opened.
    ptrVisualisation_ = new MoleculeVisualisation(ptrWindow_, kDefaultNearClip_, 
            kDefaultFarClip_);
}

MoleculeAnimation::~MoleculeAnimation() {
    std::cout << "MoleculeAnimation: Destructor." << '\n';
}

// step()
//
// The main sequence of calls in run(). This overrides Animation::step() to 
// allow a scaling of time (in seconds). Since time has been normalised such 
// that t' = omega_E t, it follows that, with t' measured in seconds, 2 * M_PI 
// seconds would correspond to one oscillation period. This is too long for a 
// visually appealing animation. To accommodate this we scale the molecular 
// times by a factor of several times 2 * M_PI. With a factor of 4.0 * 
// (2.0 * M_PI), there are approximately four lattice oscillation periods per 
// second, which is pleasing visually.

void MoleculeAnimation::step(double t1, double t2) {
    double molecTimeScale = 8.0 * M_PI;
    
    ptrVisualisation_->display(ptrSystem_);
    ptrSystem_->step(molecTimeScale * t1, molecTimeScale * t2);
    ptrVisualisation_->step(t1, t2);
}


/*
 * handleKeys:
 *
 * The keyboard event processing handler.
 *
 */
bool MoleculeAnimation::handleKeys(const GraphicsWindow::Event& event) {
    
    if (event.type == SDL_KEYDOWN) {
        switch (event.key.keysym.sym) {
            case SDLK_ESCAPE:
                return false;
                break;
                
            case SDLK_UP:
                heat();
                break;
                
            case SDLK_DOWN:
                cool();
                break;
                
            case SDLK_LEFT:
                contract();
                break;
                
            case SDLK_RIGHT:
                expand();
                break;
                
            case SDLK_g:
                toggleGravity();
                break;
                
            case SDLK_p:
                toggleCameraMotion();
                break;
                
            case SDLK_z:
                toggleCameraZoom();
                break;
                
            case SDLK_r:
                toggleCameraDirection();
                break;
                
            case SDLK_c:
                changeColour();
                break;
                
            default:
                break;
        }
    }
    
    return true;
}

/*
 * handleMouse:
 *
 * The mouse event processing handler. Presently not used.
 *
 */

bool MoleculeAnimation::handleMouse(const GraphicsWindow::Event& event) {
    return true;
}


///////////////////////////////////////////////////////////////////////////////
// Molecule System (physical) effects                                        //
///////////////////////////////////////////////////////////////////////////////

/*
 * expand()
 *
 * Expand the system container (indefinitely) without volume checks.
 * 
 */

void MoleculeAnimation::expand() {
    MoleculeSystem::Extents extents = 
                        dynamic_cast<MoleculeSystem*>(ptrSystem_)->getExtents();
    double speed = 0.1;
    
    extents.xmin -= speed * MoleculeVisualisation::kMoleculeRadius_;
    extents.xmax += speed * MoleculeVisualisation::kMoleculeRadius_;
    extents.ymin -= speed * MoleculeVisualisation::kMoleculeRadius_;
    extents.ymax += speed * MoleculeVisualisation::kMoleculeRadius_;
    extents.zmin -= speed * MoleculeVisualisation::kMoleculeRadius_;
    extents.zmax += speed * MoleculeVisualisation::kMoleculeRadius_;
    
    dynamic_cast<MoleculeSystem*>(ptrSystem_)->setExtents(extents);
}

/*
 * contract()
 *
 * Contract the system container until the volume equals the number of molecules
 * times the volume of a cube with side length equal to the molecule diameter.
 *
 * TODO: Small volumes create stress on the integration routines, whose accuracy
 * has been minimised for speed. Think about how better to deal with low volumes
 * and high energy density. Perhaps limit the contraction by keeping track of
 * total energy.   
 */
void MoleculeAnimation::contract() {
    
    MoleculeSystem* ptrMolySys = dynamic_cast<MoleculeSystem*>(ptrSystem_);
    
    MoleculeSystem::Extents extents = ptrMolySys->getExtents();
    
    double volume = (extents.xmax - extents.xmin) 
                    * (extents.ymax - extents.ymin)
                    * (extents.zmax - extents.zmin);
    
    double r = MoleculeVisualisation::kMoleculeRadius_;
    
    double num = static_cast<double>(ptrMolySys->getNumber());
    
    // Minimum volume corresponds to each molecule confined to a volume
    // equal to 2.0 * the volume of a cube of side length equal to 2 * r.
    double volume_min = 2.0 * 8.0 * r * r * r * num;
    
    // Don't allow unlimited squashing.
    if (volume <= volume_min) {
        return;
    } 
    
    double speed = -0.1;
    
    extents.xmin -= speed * MoleculeVisualisation::kMoleculeRadius_;
    extents.xmax += speed * MoleculeVisualisation::kMoleculeRadius_;
    extents.ymin -= speed * MoleculeVisualisation::kMoleculeRadius_;
    extents.ymax += speed * MoleculeVisualisation::kMoleculeRadius_;
    extents.zmin -= speed * MoleculeVisualisation::kMoleculeRadius_;
    extents.zmax += speed * MoleculeVisualisation::kMoleculeRadius_;
    
    ptrMolySys->setExtents(extents);
}

/*
 * toggleGravity()
 *
 * Alternatively turn on and turn off the gravitational field of an external
 * (planetary, say) system.
 */

void MoleculeAnimation::toggleGravity() {
    MoleculeSystem* ptrMolySys = dynamic_cast<MoleculeSystem*>(ptrSystem_);
    
    if (ptrMolySys->getGravitationalAcceleration() == 0.0) {
        ptrMolySys->setGravitationalAcceleration(0.004);
    }
    else {
        ptrMolySys->setGravitationalAcceleration(0.0);
    }
}

/*
 * heat()
 *
 * Simulate the effects of heating the molecule system, simply by scaling
 * the velocity of each molecule by a pre-determined factor.
 * 
 * This is a bit artificial and should be replaced by a process in which each
 * molecule's velocity is perturbed by an amount determined by deviates from a 
 * 3D Gaussian distribution with a given temperature.
 */

void MoleculeAnimation::heat() {
    //dynamic_cast<MoleculeSystem*>(ptrSystem_)->scaleVelocities(1.05);
    dynamic_cast<MoleculeSystem*>(ptrSystem_)->immerseHotBath();
}

/*
 * cool()
 *
 * Simulate the effects of cooling the molecule system, simply by scaling 
 * the velocity of each molecule by a pre-determined factor. See comments for
 * member function heat().
 */

void MoleculeAnimation::cool() {
    //dynamic_cast<MoleculeSystem*>(ptrSystem_)->scaleVelocities(0.95);
    dynamic_cast<MoleculeSystem*>(ptrSystem_)->immerseCoolBath();
}


///////////////////////////////////////////////////////////////////////////////
// Visualisation effects                                                     //
///////////////////////////////////////////////////////////////////////////////

/*
 * toggleCameraMotion:
 *
 * Alternatively stop or start the camera orbit.
 */

void MoleculeAnimation::toggleCameraMotion() {
    dynamic_cast<MoleculeVisualisation*>(ptrVisualisation_)->toggleCameraMotion();
}

/*
 * toggleCameraDirection()
 *
 * Reverses the direction of motion of the camera in its orbit.
 *
 */

void MoleculeAnimation::toggleCameraDirection() {
    dynamic_cast<MoleculeVisualisation*>(ptrVisualisation_)->toggleCameraDirection();
}

/*
 * toggleCameraZoom()
 *
 * Alternatively zooms in and out of the scene.
 *
 */

void MoleculeAnimation::toggleCameraZoom() {
    dynamic_cast<MoleculeVisualisation*>(ptrVisualisation_)->toggleCameraZoom();
}

/*
 * changeColour()
 *
 * Cycles through a stack of pre-set materials (colours).
 *
 */

void MoleculeAnimation::changeColour() {
    MoleculeVisualisation *ptrMolyVis = 
        dynamic_cast<MoleculeVisualisation*>(ptrVisualisation_);
    
    std::size_t num = ptrMolyVis->getNumMaterials();
    std::size_t currIndex = ptrMolyVis->getMaterialIndex();
    std::size_t newIndex = (currIndex + 1) % num;
    ptrMolyVis->setMaterialIndex(newIndex);
}