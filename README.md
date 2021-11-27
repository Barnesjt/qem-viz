# QEM Visualization

A VS2019 project for simplifying (decimating) 3d images, encoded as OBJ files. This uses a QEM algorithm to simplify the geometry of the input files. Error ellipsoids for each vertex can also be toggled

# Repository Contents

| Path | Contents |
| ---- | -------- |
| [data/geometry/](data/geometry/) | Input obj files  |
| [Debug/](Debug/) | Preloaded with required DLLs to run the project |
| [learnply/](learnply/) | Visual Studio project |
| [libraries/](libraries/)  | Glut and Glew libraries needed to run the project |
| [README.md](README.md) | This file |

# Running the Project

To run the project, launch learnply/learnply.sln with Visual Studio 2019. The project should be setup to build without the need for additional configuration.

Input OBJ files should be added to data/geometry/ prior to launching the executable.

# Controls

## Keyboard

| Key | Action |
| --- | ------ |
|  t  | Cycles through OBJ files |
|  q  | Begin mesh simplification (requires console input for the target number of faces) |
|  e  | Toggle error ellipsoid visualization |
|  r  | Reset view |

## Mouse

| Left click and drag | Pan view |
| Right click and drag | Rotate view |
| Scroll wheel | Zoom view |

# Code Attributions

- Reading all files from a directory:
	 - Adapted from https://stackoverflow.com/a/20847429
- Natural sorting (natural_sort.hpp) (for sorting file names):
	 - https://github.com/scopeInfinity/NaturalSort
- DoRasterString (on-screen information):
	 - CS450 Sample Code http://web.engr.oregonstate.edu/~mjb/cs550/
- Code to read OBJ file (adapted):
   - CS450 http://web.engr.oregonstate.edu/~mjb/cs550/loadobjfile.cpp
- Learnply base project:
   - Oregon State University CS553 (Scientific Visualization)

# References

Garland, Michael & Heckbert, Paul. (1997). Surface Simplification Using Quadric Error Metrics. Proceedings of the ACM SIGGRAPH Conference on Computer Graphics. 1997. 10.1145/258734.258849. 
