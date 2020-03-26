/* 
 * File:   Animation.h
 * Author: Richard Mace
 *
 * Created on 16 March 2020, 11:07
 */

#ifndef ANIMATION_H
#define ANIMATION_H

#include "DynamicalSystem.h"
#include "Visualisation.h"

class Animation {
public:
    Animation(DynamicalSystem*, Visualisation*, GraphicsWindow*);
    virtual ~Animation();
    
    void setInitialTime(double);
    double getInitialTime() const;
    virtual void step(double, double);
    virtual void run();
    virtual bool handleKeys(const GraphicsWindow::Event&) = 0;
    virtual bool handleMouse(const GraphicsWindow::Event&) = 0;
    
protected:
    DynamicalSystem* ptrSystem_;
    Visualisation*   ptrVisualisation_;
    GraphicsWindow*  ptrWindow_;
    double           timeInitial_;
};

#endif /* ANIMATION_H */

