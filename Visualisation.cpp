/*
 * File:   Visualisation.cpp
 * Author: Richard Mace
 *
 * March 2020
 *
 */

#include "Visualisation.h"
#include <GL/glu.h>
#include <iostream>
#include <stdexcept>


Visualisation::Visualisation(GraphicsWindow* ptrRenderTarget, double nearClip, 
                             double farClip) 

    : ptrRenderTarget_(ptrRenderTarget)
    , nearClip_(nearClip)
    , farClip_(farClip) {
    
    if (ptrRenderTarget_ != nullptr) {
        bool status = init_opengl();
    
        if (status == false) {
            throw std::runtime_error{"Visualisation: could not initialise OpenGL"};
        }
    }
}

Visualisation::~Visualisation() {
    ptrRenderTarget_ = nullptr;
    nearClip_ = farClip_ = 0.0;
}

GraphicsWindow* Visualisation::getRenderTarget() const {
    return ptrRenderTarget_;
}

void Visualisation::setRenderTarget(GraphicsWindow* ptrRenderTarget) {
    if (ptrRenderTarget != nullptr) {
        ptrRenderTarget_ = ptrRenderTarget;
        bool status = init_opengl();
        if (status == false) {
            throw std::runtime_error{"setRenderTarget: could not initialise OpenGL!\n"};
        }
    }
}

///////////////////////////////////////////////////////////////////////////////
// protected member functions                                                //
///////////////////////////////////////////////////////////////////////////////

bool Visualisation::init_opengl() {
    glClearColor(0.0f, 0.0f, 0.0f, 0.0f);
    glClearDepth(1.0f);
    glEnable(GL_DEPTH_TEST);
    glDepthFunc(GL_LESS);
    glHint(GL_PERSPECTIVE_CORRECTION_HINT, GL_NICEST);
    
    glViewport(0, 0, static_cast<GLsizei>(ptrRenderTarget_->widthPixels()), 
            static_cast<GLsizei>(ptrRenderTarget_->heightPixels()));
    
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    gluPerspective(45.0f, 
                   static_cast<GLfloat>(ptrRenderTarget_->widthPixels()) / 
                       static_cast<GLfloat>(ptrRenderTarget_->heightPixels()), 
                   static_cast<GLfloat>(nearClip_), 
                   static_cast<GLfloat>(farClip_));
    
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    
    glEnable(GL_LIGHTING);
    
    return true;
}
