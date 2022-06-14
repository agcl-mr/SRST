# SRST (Surface Reconstruction using Solitary Triangle)
CGAL implementation of the paper- Reconstruction using a simple triangle removal approach by Subhasree Methirumangalath, Amal Dev Parakkat, Shyam Sundar Kannan, Ramanathan Muthuganapathy: Published as aTechnical Brief in SIGGRAPH Asia 2017
The algorithm uses Delaunay triangualtion to reconstructs a 3D point cloud in to a mesh by retaining solitary triangles

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

