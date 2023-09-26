#---------------------------------------------------------------------#
#---------------------------------------------------------------------#
#                              Makefile                               #
#---------------------------------------------------------------------#
#            Copyright (C) 2019 J. Lopez and J. Hernandez             #
#---------------------------------------------------------------------#
#      Makefile is used to construct the VOFTools libraries           #
#      (libvoftools.a) and compile the test programs in C and         #
#      FORTRAN. Type:                                                 #
#      1) 'make' to construct the VOFTools libraries                  #
#      2) 'make fortran_3d' for the 3D test program in FORTRAN        #
#      3) 'make c_3d' for the 3D test program in C                    #
#      4) 'make fortran_2d' for the 2D test program in FORTRAN        #
#      5) 'make c_2d' for the 2D test program in C                    #
#      6) 'make all' for all the test programs in FORTRAN and C       #
#---------------------------------------------------------------------#
# This file is part of VOFTools.                                      #
#                                                                     #
# VOFTools is free software: you can redistribute it and/or           #
# modify it under the terms of the GNU General Public License         #
# as published by the Free Software Foundation, either version 3 of   #
# the License, or (at your option) any later version.                 #
#                                                                     #
# VOFTools is distributed in the hope that it will be useful,         #
# but WITHOUT ANY WARRANTY; without even the implied warranty of      #
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the       #
# GNU General Public License for more details.                        #
#                                                                     #
# You should have received a copy of the GNU General Public License   #
# along with VOFTools. If not, see <http://www.gnu.org/licenses/>.    #
#---------------------------------------------------------------------# 
#Choose your compilers:
#----------------------
# 1 for GNU compilers
# 2 for Intel compilers
#----------------------
COMPILER=1

ifeq "$(COMPILER)" "1"
CC=gcc
F77=gfortran
LIBS=-lm -lgfortran
endif

ifeq "$(COMPILER)" "2"
CC=icc
F77=ifort
LIBS=-lifcore
endif

OPT1=-O2 -o
OPT2=-O2 -c

libvoftools.a: voftools.o
	ar rv libvoftools.a voftools.o

fortran_3d: test3d.f mesh.o uservoftools.o libvoftools.a
	$(F77) $(OPT1) test3d_f test3d.f mesh.o uservoftools.o libvoftools.a

c_3d: test3d.c mesh.o uservoftools.o libvoftools.a
	$(CC) $(OPT1) test3d_c test3d.c mesh.o uservoftools.o libvoftools.a $(LIBS)

fortran_2d: test2d.f mesh.o uservoftools.o libvoftools.a
	$(F77) $(OPT1) test2d_f test2d.f mesh.o uservoftools.o libvoftools.a

c_2d: test2d.c mesh.o uservoftools.o libvoftools.a
	$(CC) $(OPT1) test2d_c test2d.c mesh.o uservoftools.o libvoftools.a $(LIBS)

all: fortran_3d c_3d fortran_2d c_2d

voftools.o: voftools.f
	$(F77) $(OPT2) voftools.f

mesh.o: mesh.f
	$(F77) $(OPT2) mesh.f

uservoftools.o: uservoftools.f
	$(F77) $(OPT2) uservoftools.f

clean:
	rm -f *.o
