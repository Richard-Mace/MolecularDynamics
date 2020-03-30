# MolecularDynamics
Simple molecular dynamics simulation and animation of a point molecule system using the Lennard-Jones interaction potential
```
V(r) = 4.0 * epsilon * [(a/r)^12 - (a/r)^6]
```
between pairs of molecules. In the above, *r* is the distance between molecules, *a* is their diameter and *epsilon* is the 
depth of the potential well (a trapping depth). As "points", each molecule has only three (translational) degrees of freedom. 
However, by removing heat from the system, the user will occasionally observe diatomic and higher molecule formation from the
single spheres, with these compound molecules exhibiting *angular* and *linear* momenta.

The animation is an example of the use of a physics animation framework that I am developing, which provides three abstract 
base classes: *DynamicalSystem*, *Visualisation* and *Animation* from which clients can derive concrete classes for animating
dynamical (physical) systems in 3D. The code also provides a *GraphicsWindow* object, which is a thin wrapper around SDL2.

The file 'main.cpp' is a simple example of usage. There are 3 types of molecule loading: SOLID, LIQUID and GAS, with average 
particle energies provided to give a representative configuration. The molecule system is visualised in 3D using hard spheres
contained in a box whose volume can be controlled by the user (see *Execution and control,* below). The camera continually orbits
the box in an interesting orbital path which can be paused and reversed by the user. The user can also simulate
adding and removing heat from the system. An external planetary gravitational field can be toggled on and off using the *g* key. 


# Build instructions
## Linux
You will need OpenGL, SDL2, and GLU libraries installed. Building is very easy using GNU g++. At the command line, simply execute:
```
$ g++ -O2 Animation.cpp DynamicalSystem.cpp GraphicsWindow.cpp MoleculeAnimation.cpp 
MoleculeSystem.cpp MoleculeVisualisation.cpp Visualisation.cpp main.cpp 
 -o moleculardynamics -lGL -lGLU -lSDL2
```
# Execution and control
Then run the executable with
```
$ ./moleculardynamics
```

## Keys
The following keys control physical effects on the system, as well as visualisation effects.

- **Up Arrow:** add heat to the system (causes melting and vapourisation)
- **Down Arrow:** remove heat from the system (causes freezing and condensation, and compound molecule formation)
- **Arrow left:** decreases the container size (slowly). Hold down for large effects.
- **Arrow right:** increases the container size (slowly). Hold down for large effects.
- **g:** toggle external gravitational field. (Watch your crystal crash to the floor.)
- **p:** alternatively pause and resume camera motion.
- **r:** reverse camera orbit direction.
- **z:** toggle camera zoom.
- **c:** change molecule colour.
