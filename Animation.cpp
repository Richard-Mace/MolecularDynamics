/*
 * File: Animation.cpp
 * Author: Richard Mace
 * March 2020
 */

#include "Animation.h"
#include <chrono>
#include <iostream>

Animation::Animation(DynamicalSystem* system, Visualisation* visualiser, 
        GraphicsWindow* window) 

    : ptrSystem_(system)
    , ptrVisualisation_(visualiser)
    , ptrWindow_(window)
    , timeInitial_(0.0) {
    
    // empty
}

// ~Animation:
//
// Destructor for the Animation class. 
//
// NOTE: The Animation object is responsible for the DELETION of the 
// GraphicsWindow object, the Visualisation object and the DynamicalSystem
// object. In effect is behaves like a smart triple pointer!
//

Animation::~Animation() {
    delete ptrWindow_;
    delete ptrVisualisation_;
    delete ptrSystem_;
}

/*
 * setInitialTime
 *
 * Sets the initial time, t_0,  for the animation.
 * 
 */

void Animation::setInitialTime(double time) {
    timeInitial_ = time;
}

/*
 * getInitialTime()
 *
 * Returns the initial time, t_0.
 * 
 */

double Animation::getInitialTime() const {
    return timeInitial_;
}

// step()
//
// The main sequence of calls in run(). Can be overridden in derived classes
// to provide additional functionality.

void Animation::step(double t1, double t2) {
    ptrVisualisation_->display(ptrSystem_);
    ptrSystem_->step(t1, t2);
    ptrVisualisation_->step(t1, t2);
}

/*
 * run()
 *
 * The main Animation loop.
 * 
 * Will return if the user closes the window, or if either handleKeys() or 
 * handleMouse() returns false.
 * 
 */

void Animation::run() {
    
    std::chrono::time_point<std::chrono::steady_clock> timeNow;
    std::chrono::time_point<std::chrono::steady_clock> timeStart;
    
    timeStart = timeNow = std::chrono::steady_clock::now();
    
    double timeLast = timeInitial_;
    
    while (true) {
    
        GraphicsWindow::Event event;
        while (ptrWindow_->pollEvent(&event)) {
            switch (event.type) {
                case SDL_QUIT:
                    return;                  // exit run loop
                    break;
                    
                case SDL_KEYDOWN:
                case SDL_KEYUP:
                    if (!handleKeys(event)) { // return false to exit run loop
                        return;
                    }
                    break;
                    
                case SDL_MOUSEMOTION:
                    if (!handleMouse(event)) { // return false to exit run loop
                        return;
                    }
                    break;
                
                default:
                    break;
            }
        }
        
        std::chrono::duration<double> elapsed_seconds = timeNow - timeStart;
        double time = elapsed_seconds.count() + timeInitial_;
        
        step(timeLast, time);
        
        timeLast = time;
        timeNow = std::chrono::steady_clock::now();
    }
}