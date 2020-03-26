/* 
 * File:   MoleculeVisualisation.h
 * Author: Richard Mace
 *
 * Created on 03 March 2020, 11:10
 */

#ifndef MOLECULEVISUALISATION_H
#define MOLECULEVISUALISATION_H

#include "GraphicsWindow.h"
#include "Visualisation.h"
#include "MoleculeSystem.h"

class MoleculeVisualisation : public Visualisation {
public:
    struct Light {
        GLfloat ambient_[4];
        GLfloat diffuse_[4];
        GLfloat specular_[4];
        GLfloat position_[4];
    };
    
    struct Material {
        GLfloat ambient_[4];   
        GLfloat diffuse_[4];   
        GLfloat specular_[4];  
        GLfloat shininess_;    
    };
    
    MoleculeVisualisation(GraphicsWindow*, double, double);
    virtual ~MoleculeVisualisation() override;
    
    void setFocalPoint(double, double, double);
    void getFocalPoint(double*, double*, double*) const;
    void setAmbientLight(GLfloat, GLfloat, GLfloat);
    void addLight(const Light&);
    void setMaterialIndex(std::size_t);
    std::size_t getMaterialIndex() const;
    std::size_t getNumMaterials() const;
    virtual void display(const DynamicalSystem*) const override;
    virtual void step(double, double) override;
    void toggleCameraMotion();
    void toggleCameraZoom();
    void toggleCameraDirection();
    
    static constexpr double kMoleculeRadius_ = 0.5; // Radius -- don't change. 
    
private:
    static const int        kSphereSegments_ = 20;  // For sphere tessellation.
    static const int        kNumLines_       = 60;  // For the surface plane.
    static constexpr double kZoomFactor_     = 0.3; // Simulates camera zoom.
    
    bool                    cameraPaused_;
    double                  cameraDirection_;
    double                  cameraZoomFactor_;
    double                  cameraPosition_[3];
    double                  focalPoint_[3]; 
    std::vector<Material>   materials_;
    std::size_t             materialIndex_;
    unsigned int            numLights_;  
    
    // utility functions
    void draw_box(const MoleculeSystem::Extents&) const;
    void orbit_camera(double, double);
    void init_material_buffer();
};

#endif /* MOLECULEVISUALISATION_H */

