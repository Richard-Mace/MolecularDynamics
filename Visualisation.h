/* 
 * File:   Visualisation.h
 * Author: Richard Mace
 *
 * Created on 03 March 2020, 11:06
 */

#ifndef VISUALISATION_H
#define VISUALISATION_H

#include "GraphicsWindow.h"
#include "DynamicalSystem.h"

class Visualisation {
public:
    Visualisation(GraphicsWindow*, double, double);
    virtual ~Visualisation();

    virtual void display(const DynamicalSystem*) const = 0;
    virtual void step(double, double) = 0;
    virtual void setRenderTarget(GraphicsWindow*);
    virtual GraphicsWindow* getRenderTarget() const;

protected:
    GraphicsWindow* ptrRenderTarget_;
    double          nearClip_;
    double          farClip_;
    
    // utility function
    bool init_opengl();
};

#endif /* VISUALISATION_H */

