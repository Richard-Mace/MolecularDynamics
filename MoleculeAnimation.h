/* 
 * File:   MoleculeAnimation.h
 * Author: Richard Mace
 *
 * Created on 17 March 2020, 15:02
 */

#ifndef MOLECULEANIMATION_H
#define MOLECULEANIMATION_H

#include "Animation.h"
#include "MoleculeSystem.h"
#include "MoleculeVisualisation.h"

class MoleculeAnimation : public Animation {
public:        
    MoleculeAnimation(std::size_t, const MoleculeSystem::Extents&, 
            MoleculeSystem::LoadType, GraphicsWindow::Type);
    
    ~MoleculeAnimation();

    virtual bool handleKeys(const GraphicsWindow::Event&) override;
    virtual bool handleMouse(const GraphicsWindow::Event&) override;
    
    // Clients of MoleculeSystem
    void expand();
    void contract();
    void heat();
    void cool();
    void toggleGravity();
    
    // Clients of MoleculeVisualisation
    void toggleCameraMotion();
    void toggleCameraDirection();
    void toggleCameraZoom();
    void changeColour();
    
private:
    static const std::size_t kDoF_                 = 3;
    static const std::size_t kDefaultWidth_        = 1920;
    static const std::size_t kDefaultHeight_       = 1080; 
    static constexpr double  kDefaultNearClip_     = 1.0;
    static constexpr double  kDefaultFarClip_      = 200.0;
    static constexpr double  kDefaultAccuracyGoal_ = 5.0e-5;
};

#endif /* MOLECULEANIMATION_H */

