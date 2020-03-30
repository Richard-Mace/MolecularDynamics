/*
 * File:   MoleculeVisualisation.cpp
 * Author: Richard Mace
 *
 * March 2020
 *
 */

#include "GraphicsWindow.h"
#include "MoleculeVisualisation.h"
#include <iostream>
#include <GL/gl.h>
#include <GL/glu.h>
                    
MoleculeVisualisation::MoleculeVisualisation(GraphicsWindow* ptrRenderTarget,
                                             double nearClip, double farClip) 

    : Visualisation(ptrRenderTarget, nearClip, farClip)
    , cameraPaused_(false)
    , cameraDirection_(1.0)
    , cameraZoomFactor_(1.0)
    , cameraPosition_{0.0, 18.0, 75.0}
    , focalPoint_{0.0, 0.0, 0.0}
    , materialIndex_(0)
    , numLights_(0) {
    
    setAmbientLight(0.4f, 0.4f, 0.4f);
    
    Light default_light { {0.5f, 0.5f, 0.5f, 1.0f},
                          {0.6f, 0.6f, 0.6f, 1.0f},
                          {0.6f, 0.6f, 0.6f, 1.0f},
                          {0.0f, 40.0f, 20.0f, 1.0f} };
                        
    addLight(default_light);
    
    Light second_light { {0.5f, 0.5f, 0.5f, 1.0f},
                         {0.6f, 0.6f, 0.6f, 1.0f},
                         {0.6f, 0.6f, 0.6f, 1.0f},
                         {0.0f, 40.0f, -20.0f, 1.0f} };
                         
    addLight(second_light);
    
    init_material_buffer();
}

    
MoleculeVisualisation::~MoleculeVisualisation() {    
    numLights_ = 0;
    
    std::cout << "MoleculeVisualisation: Destructor." << '\n';
}

void MoleculeVisualisation::setFocalPoint(double x, double y, double z) {
    focalPoint_[0] = x;
    focalPoint_[1] = y;
    focalPoint_[2] = z;
}

void MoleculeVisualisation::getFocalPoint(double* x, double* y, 
                                          double* z) const {
    *x = focalPoint_[0];
    *y = focalPoint_[1];
    *z = focalPoint_[2];
}

void MoleculeVisualisation::setAmbientLight(GLfloat r, GLfloat g, GLfloat b) {
    GLfloat ambient[] = {r, g, b, 1.0f};
    glLightModelfv(GL_LIGHT_MODEL_AMBIENT, ambient);
}

void MoleculeVisualisation::addLight(const Light& light) {
    bool fail = false;
    unsigned int flag = 0;
    
    switch (numLights_) {
        case 0:
            flag = GL_LIGHT0;
            break;
            
        case 1: 
            flag = GL_LIGHT1;
            break;
        
        case 2:
            flag = GL_LIGHT2;
            break;
            
        case 3:
            flag = GL_LIGHT3;
            break;
            
        case 4:
            flag = GL_LIGHT4;
            break;
            
        case 5:
            flag = GL_LIGHT5;
            break;
            
        case 6:
            flag = GL_LIGHT6;
            break;
            
        case 7:
            flag = GL_LIGHT7;
            break;
            
        default:
            std::cerr <<"MoleculeVisualisation: Cannot add light. Too many lights" 
                      << std::endl;
            fail = true;
            break;
    }
    
    if (!fail) {
        glLightfv(flag, GL_AMBIENT,  light.ambient_);
        glLightfv(flag, GL_DIFFUSE,  light.diffuse_);
        glLightfv(flag, GL_SPECULAR, light.specular_);
        glLightfv(flag, GL_POSITION, light.position_);
        glEnable(flag);
        
        ++numLights_;
    }
}

void MoleculeVisualisation::setMaterialIndex(std::size_t index) {
    if (index <= materials_.size()) {
        materialIndex_ = index;
    }
}

std::size_t MoleculeVisualisation::getMaterialIndex() const {
    return materialIndex_;
}

std::size_t MoleculeVisualisation::getNumMaterials() const {
    return materials_.size();
}

void MoleculeVisualisation::display(const DynamicalSystem* ptrSystem) const {
    
    // Clear (fill) the frame and depth (z-) buffers.
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    
    glLoadIdentity();
    
    // Follow a molecule for fun. A bit jiggly. Experiment with it.
    //const double* ptrState = ptrSystem->getStatePtr();
    //setFocalPoint(ptrState[0], ptrState[1], ptrState[2]);
    
    // Set up eye (coordinate system) transformations.
    gluLookAt(cameraPosition_[0], cameraPosition_[1], cameraPosition_[2], 
              focalPoint_[0], focalPoint_[1], focalPoint_[2], 0.0, 1.0, 0.0);
    
    // Draw the confining box.
    MoleculeSystem::Extents extents = 
            dynamic_cast<const MoleculeSystem*>(ptrSystem)->getExtents();
    
    draw_box(extents);
    
    // Set the material properties for all molecules.
    glMaterialfv(GL_FRONT, GL_SPECULAR, materials_[materialIndex_].specular_);
    glMaterialfv(GL_FRONT, GL_AMBIENT,  materials_[materialIndex_].ambient_);
    glMaterialfv(GL_FRONT, GL_DIFFUSE,  materials_[materialIndex_].diffuse_);
    glMaterialf(GL_FRONT, GL_SHININESS, materials_[materialIndex_].shininess_);
    
    const size_t degreesOfFreedom = 
        dynamic_cast<const MoleculeSystem*>(ptrSystem)->getDegreesOfFreedom();
    
    GLUquadric* ptr_sphere = gluNewQuadric();
    
    const size_t numMolecules = 
        dynamic_cast<const MoleculeSystem*>(ptrSystem)->getNumber();
    
    for (std::size_t i = 0; i < numMolecules; ++i) {
        std::size_t offset = i * degreesOfFreedom;
        GLfloat x = static_cast<GLfloat>(ptrSystem->state()[offset + 0]);
        GLfloat y = static_cast<GLfloat>(ptrSystem->state()[offset + 1]);
        GLfloat z = static_cast<GLfloat>(ptrSystem->state()[offset + 2]);
        
        // Position and draw molecule.
        glPushMatrix();
        glTranslatef(y, z, x);
        gluSphere(ptr_sphere, kMoleculeRadius_, kSphereSegments_, 
                kSphereSegments_);
        glPopMatrix();
    }
    
    gluDeleteQuadric(ptr_sphere);
    
    // Swap backbuffer with frontbuffer.
    ptrRenderTarget_->swapBuffers();
}

void MoleculeVisualisation::step(double t1, double t2) {
    orbit_camera(t1, t2);
}

void MoleculeVisualisation::toggleCameraMotion() {
    cameraPaused_ = !cameraPaused_;
}

void MoleculeVisualisation::toggleCameraZoom() {
    if (cameraZoomFactor_ == 1.0) {
        cameraZoomFactor_ = kZoomFactor_;
    }
    else {
        cameraZoomFactor_ = 1.0;
    }
}

void MoleculeVisualisation::toggleCameraDirection() {
    if (cameraDirection_ == 1.0) {
        cameraDirection_ = -1.0;
    }
    else {
        cameraDirection_ = 1.0;
    }
}

///////////////////////////////////////////////////////////////////////////////
// private utility functions                                                 //
///////////////////////////////////////////////////////////////////////////////

void MoleculeVisualisation::draw_box(
                                const MoleculeSystem::Extents& extents) const {
    
    GLfloat spec1[] = {0.9f, 0.9f, 0.9f, 1.0f};
    GLfloat amb1[] = {1.0f, 1.0f, 1.0f, 1.0f};
    GLfloat diff1[] = {1.0f, 1.0f, 1.0f, 1.0f};
    GLfloat shin1[] = {32.0f};
    
    glMaterialfv(GL_FRONT, GL_SPECULAR, spec1);
    glMaterialfv(GL_FRONT, GL_SHININESS, shin1);
    glMaterialfv(GL_FRONT, GL_AMBIENT, amb1);
    glMaterialfv(GL_FRONT, GL_DIFFUSE, diff1);
    
    // Account for the finite size of the "point" molecules
    GLfloat xmin = static_cast<GLfloat>(extents.xmin) - kMoleculeRadius_;
    GLfloat xmax = static_cast<GLfloat>(extents.xmax) + kMoleculeRadius_;
    GLfloat ymin = static_cast<GLfloat>(extents.ymin) - kMoleculeRadius_;
    GLfloat ymax = static_cast<GLfloat>(extents.ymax) + kMoleculeRadius_;
    GLfloat zmin = static_cast<GLfloat>(extents.zmin) - kMoleculeRadius_;
    GLfloat zmax = static_cast<GLfloat>(extents.zmax) + kMoleculeRadius_;
    
    glBegin(GL_LINES);
    
    glVertex3f(xmin, ymin, zmin);
    glVertex3f(xmax, ymin, zmin);
    glVertex3f(xmax, ymin, zmin);
    glVertex3f(xmax, ymax, zmin);
    glVertex3f(xmax, ymax, zmin);
    glVertex3f(xmin, ymax, zmin);
    glVertex3f(xmin, ymax, zmin);
    glVertex3f(xmin, ymin, zmin);
    
    glVertex3f(xmin, ymin, zmax);
    glVertex3f(xmax, ymin, zmax);
    glVertex3f(xmax, ymin, zmax);
    glVertex3f(xmax, ymax, zmax);
    glVertex3f(xmax, ymax, zmax);
    glVertex3f(xmin, ymax, zmax);
    glVertex3f(xmin, ymax, zmax);
    glVertex3f(xmin, ymin, zmax);
    
    glVertex3f(xmin, ymin, zmin);
    glVertex3f(xmin, ymin, zmax);
    
    glVertex3f(xmax, ymin, zmin);
    glVertex3f(xmax, ymin, zmax);
    
    glVertex3f(xmax, ymax, zmin);
    glVertex3f(xmax, ymax, zmax);
    
    glVertex3f(xmin, ymax, zmin);
    glVertex3f(xmin, ymax, zmax);
    
    // draw floor
    GLfloat deltax = 20.5f / 4.0f;
    GLfloat deltay = 20.5f / 4.0f;
   
    GLfloat col_green[] {0.0f, 0.25f, 0.0f, 1.0f};
    
    glMaterialfv(GL_FRONT, GL_AMBIENT_AND_DIFFUSE, col_green);
    for (int i = 0; i < kNumLines_; ++i) {
        glVertex3f(i * deltax, zmin - 0.1f, -kNumLines_ * deltay);
        glVertex3f(i * deltax, zmin - 0.1f,  kNumLines_ * deltay);
        
        glVertex3f(-kNumLines_ * deltax, zmin - 0.1f, i * deltay);
        glVertex3f(kNumLines_ * deltax, zmin - 0.1f,  i * deltay);
    }
    
    for (int i = 1; i < kNumLines_; ++i) {
        glVertex3f(-i * deltax, zmin - 0.1f, -kNumLines_ * deltay);
        glVertex3f(-i * deltax, zmin - 0.1f,  kNumLines_ * deltay);
        
        glVertex3f(-kNumLines_ * deltax, zmin - 0.1f, -i * deltay);
        glVertex3f(kNumLines_ * deltax, zmin - 0.1f,  -i * deltay);
    }
     
    glEnd();
}

void MoleculeVisualisation::orbit_camera(double t1, double t2) {

    const double k_cam_max_r = 50.0;
    const double k_cam_min_r = 30.0;
    const double k_cam_max_z = 30.0;
    const double k_slow = 0.5;
    
    static double cam_time = 0.0;

    
    double delta_t;
    if (cameraPaused_) {
        delta_t = 0.0;
    }
    else {
        delta_t = t2 - t1;
    }
    
    if (cameraDirection_ == -1.0) {
        delta_t *= -1.0;
    }

    cam_time += delta_t;

    double r = k_cam_max_r + (k_cam_max_r - k_cam_min_r) * 
        cos(0.1 * cam_time * k_slow);
    cameraPosition_[0] = cameraZoomFactor_ * r * cos(0.6 * cam_time * k_slow);
    cameraPosition_[2] = cameraZoomFactor_ * r * sin(0.3 * cam_time * k_slow);
    cameraPosition_[1] = cameraZoomFactor_ * k_cam_max_z * 
            (sin(0.5 * cam_time * k_slow) + 1.02);
}

/*
 * init_material_buffer:
 *
 * Defines a few pre-set materials for the user to cycle through, if desired.
 *
 */

void MoleculeVisualisation::init_material_buffer() {
    
    Material def_mat { {0.5f, 0.1995f, 0.0745f, 1.0f},
                       {0.75164f, 0.60648f, 0.22648f, 1.0f},
                       {0.9f, 0.9f, 0.9f, 1.0f},
                        32.0f };
                        
    materials_.push_back(def_mat);
    
    Material emerald { {0.0215f, 0.1745f, 0.0215f, 1.0f},
                       {0.07568f, 0.61424f, 0.07568f, 1.0f},
                       {0.633f, 0.727811f, 0.633f, 1.0f}, 
                        76.8f };
                       
    materials_.push_back(emerald);
    
    Material obsidian { {0.05375f, 0.05f, 0.06625f, 1.0f},
                        {0.18275f, 0.17f, 0.22525f, 1.0f},
                        {0.332741f, 0.328634f, 0.346435f, 1.0f}, 
                         38.4f };
                         
    materials_.push_back(obsidian);
    
    Material ruby { {0.1745f, 0.01175f, 0.01175f, 1.0f},
                    {0.61424f, 0.04136f, 0.04136f, 1.0f},
                    {0.727811f, 0.626959f, 0.626959f, 1.0f},
                     76.8f };
                     
    materials_.push_back(ruby);
    
    Material chrome { { 0.25f, 0.25f, 0.25f, 1.0f}, 
                      {0.4f, 0.4f, 0.4f, 1.0f},
                      {0.774597f, 0.774597f, 0.774597f, 1.0f},
                       76.8f };
                      
    materials_.push_back(chrome);
    
    Material cyan_plastic{ {0.0f, 0.1f, 0.06f, 1.0f},
                           {0.0f, 0.50980392f, 0.50980392f, 1.0f},
                           {0.50196078f, 0.50196078f, 0.50196078f, 1.0f},
                            32.0f };
                            
    materials_.push_back(cyan_plastic);    
    
    Material copper { {0.19125f, 0.0735f, 0.0225f, 1.0f},
                      {0.7038f, 0.27048f, 0.0828f, 1.0f},
                      {0.25677f, 0.137622f, 0.086014f, 1.0f}, 
                       12.8f };
                       
    materials_.push_back(copper);   
}