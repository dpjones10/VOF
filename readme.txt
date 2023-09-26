-----------------------------------------------------------------------
|                           V O F T o o l s                           |
|                                                                     |
|	 A package of FORTRAN subroutines with analytical and         |
|                                                                     |  
|      geometrical tools for 2D/3D VOF methods in general convex      |
|                                                                     |
|                     grids  and Cartesian geometry                   |
|                                                                     |
|                                                                     |
|                        (Version 3.2, June 2019)                     |
|                                                                     |
|				                                      |
|             Copyright (C) 2019 J. Lopez and J. Hernandez            |
|                                                                     |
|                                                                     |
-----------------------------------------------------------------------

-----------------------------------------------------------------------

VOFTools is a collection of FORTRAN subroutines with analytical and 
geometrical tools for 2D/3D Volume of Fluid (VOF) methods in
general convex grids and Cartesian geometry.

This library is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as
published by the Free Software Foundation, either version 3 of the
License, or (at your option) any later version.
                                                                     
This library is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
General Public License for more details.
                                                                     
You should have received a copy of the GNU General Public License 
along with this library.  If not, see <http://www.gnu.org/licenses/>.

Publications reporting results obtained using VOFTools should make
appropriate citation to the software and papers that describe it:

[1] J. Lopez and J. Hernandez, Analytical and geometrical tools for 3D
volume of fluid methods in general grids, Journal of Computational
Physics 227 (2008) 5939-5948. 

[2] J. Lopez, P. Gomez, J. Hernandez, F. Faura, A two-grid adaptive
volume of fluid approach for dendritic solidification, Computers &
Fluids 86 (2013) 326-342.

[3] J. Lopez, J. Hernandez, P. Gomez, F. Faura, A new volume 
conservation enforcement method for PLIC reconstruction in general 
convex grids, Journal of Computational Physics 316 (2016) 338-359.

[4] J. Lopez, J. Hernandez, P. Gomez, F. Faura, VOFTools - A software
package of calculation tools for volume of fluid methods using general
convex grids, Computer Physics Communications 223 (2018) 45-54.

[5] J. Lopez, J. Hernandez, P. Gomez, C. Zanzi, R. Zamora, VOFTools 3.2 
- A software package of calculation tools for VOF methods with an added 
capability to initialize the liquid volume fraction in general convex
grids, submitted to Computer Physics Communications.

For more information contact joaquin.lopez@upct.es 

-----------------------------------------------------------------------
                       V O F T o o l s  3.2
-----------------------------------------------------------------------

In this directory you will find the following files: 

user-manual.pdf ---> instructions for the library installation and test
                     programs execution and description of the meaning
		     of routine arguments

voftools.f      ---> source code of the analytical and geometrical tools

uservoftools.f  ---> source code of user-defined functions.

mesh.f          ---> definitions of different cell types

dim.h           ---> array dimensions for FORTRAN codes

dimc.h          ---> array dimensions for C codes

cvoftools.h     ---> declaration of the subroutines of the VOFTools 
                     library to be used in C codes

cuservoftools.h ---> declaration of the user-defined functions included in
                     the uservoftools.f file to be used in C codes

cmesh.h         ---> declaration of the subroutines used to define the 
                     cell geometries considered in the test programs 
                     to be used in C codes

test2d.f        ---> 2D test program in FORTRAN

test2d.c        ---> 2D test program in C

test3d.f        ---> 3D test program in FORTRAN

test3d.c        ---> 3D test program in C

Makefile        ---> constructs the VOFTools library and makes executable
                     files for test programs in C and Fortran

change.log      ---> list of notable changes made to the VOFTools package

COPYING         ---> copy of the GNU General Public License, Version 3

-----------------------------------------------------------------------
                            INSTALLATION
-----------------------------------------------------------------------
To install the VOFTools library and execute the test programs supplied
in FORTRAN and C, perform the following steps:

1. Decompress the downloaded package in the working directory:

   tar -zxvf voftools-3.2.tgz

2. Go to the VOFTools directory:

   cd voftools-3.2

3. Edit the Makefile file and choose the compilers (COMPILER = 1 
for GNU compilers and 2 for Intel compilers). The user can introduce 
other compilers by setting variables CC, F77 and LIBS.

4. Build the VOFTools library (libvoftools.a):

   make

5. Type

   make all

to compile the four versions of the test program (test2d_c, test2d_f,
test3d_c and test3d_f).

Optionally, the built VOFTools library can be moved to a search path 
(for example, /usr/lib/) and the test program (testprogram.f) can be 
compiled as follows:

   sudo cp libvoftools.a /usr/lib
   sudo chmod 777 /usr/lib/libvoftools.a
   ifort -o testprogram testprogram.f -lvoftools

6. Execute the test program/s.

If the user needs to modify the values of parameters ns or nv in the dim.h
file for the FORTRAN version of the code or in the dimc.h file for the C
version of the code, for example in order to use cells with a number of
faces or vertices higher than that initially specified in the above files,
the VOFTools library must be recompiled. To recompile the VOFTools library
in the same directory, type first

   make clean

and then proceed as indicated above.
