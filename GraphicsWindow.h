/* 
 * File:   GraphicsWindow.h
 * Author: Richard Mace
 *
 * Created on 04 March 2020, 09:58
 */

#ifndef GRAPHICSWINDOW_H
#define GRAPHICSWINDOW_H

#include <SDL2/SDL.h>
#include <SDL2/SDL_opengl.h>
#include <string>

class GraphicsWindow {
public:
    using Event     = SDL_Event;
    using Canvas    = SDL_Window;
    using GLcontext = SDL_GLContext;
    
    enum class Type {WINDOWED, FULLSCREEN};

    GraphicsWindow();
    ~GraphicsWindow();
    
    bool isInitialised() const;
    bool open(); // fullscreen
    bool open(std::size_t, std::size_t); // windowed
    void close();
    
    std::size_t widthPixels() const;
    std::size_t heightPixels() const;
    void setTitle(const std::string&);
    std::string getTitle() const;
    bool pollEvent(Event*);
    void swapBuffers();
    
private:
    bool        initialised_;
    Canvas*     pCanvas_;
    GLcontext   glContext_;
    std::size_t widthPixels_;
    std::size_t heightPixels_;
    std::string windowTitle_;

    // helper functions
    bool openHelper(std::size_t, std::size_t, Type);
};

#endif /* GRAPHICSWINDOW_H */

