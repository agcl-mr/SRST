# SRST 

Surface Reconstruction using Solitary Triangle, is a program for reconstructing a 3D model retaining the solitary triangles.
Copyright (C) 2017 Subhasree Methirumangalath, Amal Dev Parakkat, Shyam Sundar Kannan, , and Ramanathan Muthuganapathy, Advanced Geometric Computing Lab, Department of Engineering Design, IIT Madras, Chennai, India.

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

Surface Reconstruction using Solitary Triangle is proposed in the paper - Reconstruction using a simple triangle removal approach by Subhasree Methirumangalath, Amal Dev Parakkat, Shyam Sundar Kannan, Ramanathan Muthuganapathy:accepted as a Technical Brief in SIGGRAPH Asia 2017
The algorithm uses Delaunay triangualtion to reconstructs a 3D point cloud in to a mesh by retaining solitary triangles.

Compilation:
The code needs CGAL (Computational Geometry Algorithm Library) for compilation. 
Refer to https://www.cgal.org/ for CGAL installation and compilation process.

Usage: 
The program needs the input and output file names as arguments
The program takes .xyz file as an input and creates a .stl file as an output.

$ ./srst <input .xyz file name> <output .stl file name>

Example
$ ./srst horse.xyx horse.stl

Parameter: 
The algorithm needs a parameter to obtain ideal results.
The parameter is used to specify the circumradius range for the triangles to be retained apart from solitary triangles.
Enter the parameter when prompted. (Note: The parameter should always be > 0)

Hint: The parameter can be tuned by visually inspecting the result. 
If the result many extra triangles decrease the value of the parameter and vice versa.

Refer to Sample Input folder for input samples.

