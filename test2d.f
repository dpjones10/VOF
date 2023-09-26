c---------------------------------------------------------------------c
c---------------------------------------------------------------------c
c                              TEST2D                                 c
c---------------------------------------------------------------------c
c            Copyright (C) 2019 J. Lopez and J. Hernandez             c
c---------------------------------------------------------------------c
c      Test program to solve the local volume enforcement problem     c
c      and to calculate the distance from a given point to the        c
c      interfacial segment                                            c
c---------------------------------------------------------------------c
c This file is part of VOFTools.                                      c
c                                                                     c
c VOFTools is free software: you can redistribute it and/or           c
c modify it under the terms of the GNU General Public License         c
c as published by the Free Software Foundation, either version 3 of   c
c the License, or (at your option) any later version.                 c
c                                                                     c
c VOFTools is distributed in the hope that it will be useful,         c
c but WITHOUT ANY WARRANTY; without even the implied warranty of      c
c MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the       c
c GNU General Public License for more details.                        c
c                                                                     c
c You should have received a copy of the GNU General Public License   c
c along with VOFTools. If not, see <http://www.gnu.org/licenses/>.    c
c---------------------------------------------------------------------c 
      INCLUDE "dim.h"
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
c* Original polygon
      DIMENSION IPV(NV),VERTP(NV,2)
c* Work polygon 0
      DIMENSION IPV0(NV),VERTP0(NV,2)
c* Segment
      DIMENSION X(2),Y(2)
c* External function to define the interface shape      
      EXTERNAL FUNC2D1, FUNC2D2
      WRITE(6,*)'-----------------------------------------------------'
      WRITE(6,*)'|---------------------------------------------------|'
      WRITE(6,*)'|          TEST PROGRAM IN 2D OF VOFTools           |'
      WRITE(6,*)'|                                                   |'
      WRITE(6,*)'|             (Version 3.2, June 2019)              |'
      WRITE(6,*)'|                                                   |'
      WRITE(6,*)'|                       by                          |'
      WRITE(6,*)'|                                                   |'
      WRITE(6,*)'|            J. Lopez and J. Hernandez              |'
      WRITE(6,*)'|---------------------------------------------------|'
      WRITE(6,*)'-----------------------------------------------------'
      WRITE(6,*)''
c* Material body shape to be initialized:
c----------------------------------------
c. ISHAPE  =  1, circle with radious 0.325 centered at (0.5,0.5)
c.         =  2, ellipse with semi-major axis 0.6, semi-minor axis 0.2
c.               and centered at (0.5,0.5)      
      ISHAPE=1
c* Subdivision number in the volume fraction cell initialization:
      NC=10
c* Tolerance for the initialization procedure:
      TOL=10.0
c* Choose the desired cell type:
c------------------------------
c* Square mesh
      CALL SQUAREMESH(IPV,NTP,NTV,VERTP)
c* Hexagonal mesh
c      CALL HEXAGOMESH(IPV,NTP,NTV,VERTP)
c* Triangular mesh
c      CALL TRIANGLEMESH(IPV,NTP,NTV,VERTP)
c* Quadrangular mesh
c      CALL QUADRANGLEMESH(IPV,NTP,NTV,VERTP)
c* Pentagonal mesh
c      CALL PENTAGONMESH(IPV,NTP,NTV,VERTP)
c* Irregular hexagonal mesh
c      CALL HEXAGONMESH(IPV,NTP,NTV,VERTP)
c* Calculate the area VT of the cell:
c------------------------------------
      CALL TOOLV2D(IPV,NTV,VERTP,VT)
      WRITE(6,*)'-----------------------------------------------------'
      WRITE(6,*)'|---------------------------------------------------|'
      WRITE(6,*)'|         OUTPUT OF THE TOOLV2D SUBROUTINE          |'
      WRITE(6,*)'|---------------------------------------------------|'
      WRITE(6,*)'-----------------------------------------------------'
      WRITE(6,*)''
      WRITE(6,*)'Area of the selected cell:',VT
      WRITE(6,*)''
c* Choose the desired volume fraction of liquid
      F=0.5
c* Calculate the volume of liquid in the cell
      V=F*VT
c* Choose the desired unit-length normal vector of the truncation plane
      XNC=0.70710678
      YNC=0.70710678
c* Solve the local volume enforcement problem:
c--------------------------------------------
      CALL ENFORV2D(C,IPV,NTP,NTV,V,VT,VERTP,XNC,YNC)
c.. For rectangular cells, like that defined in the SQUAREMESH subroutine,
c.. it is more efficient to use the analytical method of Scardovelli and 
c.. Zaleski [Journal of Computational Physics, 164 (2000) 228-237], which
c.. was proposed to be used specifically for this type of cells. To use
c.. this method, uncomment the following three lines.
c      DX=DABS(VERTP(1,1)-VERTP(2,1))    !x-side length
c      DY=DABS(VERTP(2,2)-VERTP(3,2))    !y-side length
c      CALL ENFORV2DSZ(C,DX,DY,V,VERTP,XNC,YNC)
      WRITE(6,*)'-----------------------------------------------------'
      WRITE(6,*)'|---------------------------------------------------|'
      WRITE(6,*)'|        OUTPUT OF THE ENFORV2D SUBROUTINE          |'
      WRITE(6,*)'|---------------------------------------------------|'
      WRITE(6,*)'-----------------------------------------------------'
      WRITE(6,*)''
      WRITE(6,*)'Solution of the problem:',C
      WRITE(6,*)''
c* Calculate the distance from a given point P to the interfacial segment:
c------------------------------------------------------------------------
c* copy the original polygon that defines the cell to the working
c* polygon 0
      CALL CPPOL2D(IPV,IPV0,NTP,NTP0,NTV,NTV0,VERTP,VERTP0)
c* intersection between the cell and the line defined as X·NC+C=0
      CALL INTE2D(C,ICONTN,ICONTP,IPV0,NTP0,NTV0,VERTP0,XNC,YNC)
c* choose the desired point P from which the distance is calculated
      XP=0.0
      YP=0.0
c* interfacial segment defined as the last edge of the truncated polygon 0
      X(1)=VERTP0(NTP0-1,1)
      Y(1)=VERTP0(NTP0-1,2)
      X(2)=VERTP0(NTP0,1)
      Y(2)=VERTP0(NTP0,2)
c* calculate the distance
      CALL DIST2D(D,X,Y,XP,YP)
      WRITE(6,*)'-----------------------------------------------------'
      WRITE(6,*)'|---------------------------------------------------|'
      WRITE(6,*)'|         OUTPUT OF THE DIST2D SUBROUTINE           |'
      WRITE(6,*)'|---------------------------------------------------|'
      WRITE(6,*)'-----------------------------------------------------'
      WRITE(6,*)''
      WRITE(6,*)'Distance from P to the interfacial segment:',D
      WRITE(6,*)''
      WRITE(6,*)'-----------------------------------------------------'
c* Initialize the material volume fraction in the cell:
c-----------------------------------------------------      
      IF(ISHAPE.EQ.1) THEN
         CALL INITF2D(FUNC2D1,IPV,NC,NTP,NTV,TOL,VERTP,VF)
      ELSEIF(ISHAPE.EQ.2) THEN
         CALL INITF2D(FUNC2D2,IPV,NC,NTP,NTV,TOL,VERTP,VF)
      END IF
      WRITE(6,*)'-----------------------------------------------------'
      WRITE(6,*)'|---------------------------------------------------|'
      WRITE(6,*)'|         OUTPUT OF THE INITF3D SUBROUTINE          |'
      WRITE(6,*)'|---------------------------------------------------|'
      WRITE(6,*)'-----------------------------------------------------'
      WRITE(6,*)''
      WRITE(6,*)'Material volume fraction in the selected cell:',VF
      WRITE(6,*)''
      WRITE(6,*)'-----------------------------------------------------'
      END
c-------------------------- END OF TEST2D ----------------------------c
c---------------------------------------------------------------------c
