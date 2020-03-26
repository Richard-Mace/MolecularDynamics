/*
 * File: GraphicsWindow.cpp
 * Author: Richard Mace
 *
 * March 2020
 */

#include "GraphicsWindow.h"
#include <iostream>


GraphicsWindow::GraphicsWindow() 

    : initialised_(false)
    , pCanvas_(nullptr)
    , glContext_(nullptr)
    , widthPixels_(0)
    , heightPixels_(0) 
    , windowTitle_("Graphics Window") {
    
    if (SDL_Init(SDL_INIT_VIDEO) < 0) {
        std::cerr << "Could not initialise SDL: " << SDL_GetError() << '\n';
    } 
    else {
        initialised_ = true;
    }
}

GraphicsWindow::~GraphicsWindow() {
    if (pCanvas_) {
        std::cerr << "Destroying SDL window." << std::endl;
        SDL_DestroyWindow(pCanvas_);
    }
    
    initialised_ = false;
    
    std::cerr << "GraphicsWindow::Destructor." << std::endl;
    
    SDL_Quit();
}

bool GraphicsWindow::isInitialised() const {
    return initialised_;
}

void GraphicsWindow::setTitle(const std::string& title) {
    windowTitle_ = title;
    
    if (pCanvas_ != nullptr) {
        SDL_SetWindowTitle(pCanvas_, windowTitle_.c_str());
    }
}

std::string GraphicsWindow::getTitle() const {
    return windowTitle_;
}

// Open fullscreen graphics window.

bool GraphicsWindow::open() {
    return openHelper(0, 0, Type::FULLSCREEN);
}

// Open windowed graphics window.

bool GraphicsWindow::open(std::size_t width, std::size_t height) {
    return openHelper(width, height, Type::WINDOWED);
}

// Close the graphics window.

void GraphicsWindow::close() {
    if (!initialised_) {
        return;
    }
    
    std::cerr << "Closing graphics window" << std::endl;
    SDL_DestroyWindow(pCanvas_);
    
    pCanvas_      = nullptr;
    glContext_    = nullptr;
    widthPixels_  = 0;
    heightPixels_ = 0;
}

std::size_t GraphicsWindow::widthPixels() const {
    return widthPixels_;
}

std::size_t GraphicsWindow::heightPixels() const {
    return heightPixels_;
}

bool GraphicsWindow::pollEvent(Event* ptrEvent) {
    return SDL_PollEvent(ptrEvent);
}

void GraphicsWindow::swapBuffers() {
    SDL_GL_SwapWindow(pCanvas_);
}

///////////////////////////////////////////////////////////////////////////////
// private utility/helper functions                                          //
///////////////////////////////////////////////////////////////////////////////

bool GraphicsWindow::openHelper(std::size_t width, std::size_t height, 
                                Type type) {
    if (!initialised_) {
        return false;
    }
    
    // Tell SDL to use OpenGL 2.1, for now.
    SDL_GL_SetAttribute(SDL_GL_CONTEXT_MAJOR_VERSION, 2);
    SDL_GL_SetAttribute(SDL_GL_CONTEXT_MINOR_VERSION, 1);
    
    unsigned int flags = SDL_WINDOW_OPENGL | SDL_WINDOW_SHOWN;
    
    if (type == Type::FULLSCREEN) {
        flags |= SDL_WINDOW_FULLSCREEN_DESKTOP;
    }
    // Create window.
    Canvas* pCanvas = SDL_CreateWindow(
                                windowTitle_.c_str(), 
                                SDL_WINDOWPOS_UNDEFINED,
                                SDL_WINDOWPOS_UNDEFINED, 
                                static_cast<int>(width), 
                                static_cast<int>(height), 
                                flags
                                );
    
    if (pCanvas == nullptr) {
        std::cerr << "Could not create window: " << SDL_GetError() << '\n';
        
        return false;
    }

    // Create OpenGL context.
    GLcontext glContext = SDL_GL_CreateContext(pCanvas);
    if (glContext == nullptr) {
        std::cerr << "Could not create OpenGL context: " << SDL_GetError() 
                  << '\n';
        
        return false;
    }
    
    // Sychronise. Use vsync.
    if (SDL_GL_SetSwapInterval(1) < 0) {
        std::cerr << "WARNING: unable to VSYNC: " << SDL_GetError() << '\n';
        
        return false;
    }
    
    SDL_ShowCursor(SDL_DISABLE);
    

    // Get screen dimensions in case of fullscreen usage. SDL requires int's.
    int w;
    int h;
    SDL_GL_GetDrawableSize(pCanvas, &w, &h);
    
    // All successful. Update member variables.
    pCanvas_      = pCanvas;
    glContext_    = glContext;
    widthPixels_  = static_cast<std::size_t>(w);
    heightPixels_ = static_cast<std::size_t>(h);
    
    return true;
}
