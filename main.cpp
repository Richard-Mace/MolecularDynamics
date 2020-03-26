/* 
 * File:   main.cpp
 * Author: Richard Mace
 *
 * Created on 03 March 2020, 10:39
 */

#include "GraphicsWindow.h"
#include "MoleculeSystem.h"
#include "MoleculeVisualisation.h"
#include "MoleculeAnimation.h"


int main() {
    MoleculeSystem::Extents  extents {-20.0, 20.0, -20.0, 20.0, -20.0, 20.0};
    MoleculeSystem::LoadType loadType   {MoleculeSystem::LoadType::SOLID};
    GraphicsWindow::Type     windowType {GraphicsWindow::Type::FULLSCREEN};

    MoleculeAnimation molecules(100, extents, loadType, windowType);
    
    molecules.run();
    
    return 0;
}

