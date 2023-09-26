c---------------------------------------------------------------------c
c---------------------------------------------------------------------c
c                              TEST3D                                 c
c---------------------------------------------------------------------c
c            Copyright (C) 2019 J. Lopez and J. Hernandez             c
c---------------------------------------------------------------------c
c      Test program to solve the local volume enforcement problem     c
c      and to calculate the distance from a given point to the        c
c      interfacial polygon                                            c 
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
      PROGRAM TEST3D
      INCLUDE "dim.h"
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
c* Original polyhedron
      DIMENSION CS(NS),IPV(NS,NV),NIPV(NS),VERTP(NV,3),XNS(NS),YNS(NS),
     -     ZNS(NS)
c* Working polyhedron 0
      DIMENSION CS0(NS),IPV0(NS,NV),NIPV0(NS),VERTP0(NV,3),XNS0(NS),
     -     YNS0(NS),ZNS0(NS)
c* Interfacial polygon
      DIMENSION X(NV),Y(NV),Z(NV)
c* External function to define the interface shape      
      EXTERNAL FUNC3D1,FUNC3D2
      WRITE(6,*)'-----------------------------------------------------'
      WRITE(6,*)'|---------------------------------------------------|'
      WRITE(6,*)'|          TEST PROGRAM IN 3D OF VOFTools           |'
      WRITE(6,*)'|                                                   |'
      WRITE(6,*)'|            (Version 3.2, June 2019)               |'
      WRITE(6,*)'|                                                   |'
      WRITE(6,*)'|                       by                          |'
      WRITE(6,*)'|                                                   |'
      WRITE(6,*)'|            J. Lopez and J. Hernandez              |'
      WRITE(6,*)'|---------------------------------------------------|'
      WRITE(6,*)'-----------------------------------------------------'
      WRITE(6,*)''
c* Material body shape to be initialized:
c----------------------------------------
c. ISHAPE  =  1, sphere with radious 0.6 centered at (0.5,0.5,0.5)
c.         =  2, Torus with major radius 2/3, minor radius 1/3 and
c.               centered at (0.5,0.5,0.5)
      ISHAPE=1
c* Subdivision number in the volume fraction cell initialization:
      NC=10
c* Tolerance for the initialization procedure:
      TOL=10.0
c* Choose the desired cell type:
c------------------------------
c* Cubic mesh
      CALL CUBICMESH(IPV,NIPV,NTP,NTS,NTV,VERTP,XNS,YNS,ZNS)
c* General hexahedrical mesh
c      CALL HEXAHEMESH(IPV,NIPV,NTP,NTS,NTV,VERTP,XNS,YNS,ZNS)
c* Tetrahedrical mesh
c      CALL TETRAMESH(IPV,NIPV,NTP,NTS,NTV,VERTP,XNS,YNS,ZNS)
c* Dodecahedral mesh
c      CALL DODECAMESH(IPV,NIPV,NTP,NTS,NTV,VERTP,XNS,YNS,ZNS)
c* Icosahedral mesh
c      CALL ICOSAMESH(IPV,NIPV,NTP,NTS,NTV,VERTP,XNS,YNS,ZNS)
c* Complex mesh
c      CALL COMPLEXMESH(IPV,NIPV,NTP,NTS,NTV,VERTP,XNS,YNS,ZNS)
c* Calculate the volume VT of the cell:
c-------------------------------------
      CALL TOOLV3D(IPV,NIPV,NTS,VERTP,VT,XNS,YNS,ZNS)
      WRITE(6,*)'-----------------------------------------------------'
      WRITE(6,*)'|---------------------------------------------------|'
      WRITE(6,*)'|         OUTPUT OF THE TOOLV3D SUBROUTINE          |'
      WRITE(6,*)'|---------------------------------------------------|'
      WRITE(6,*)'-----------------------------------------------------'
      WRITE(6,*)''
      WRITE(6,*)'Volume of the selected cell:',VT
      WRITE(6,*)''
c* choose the desired volume fraction of liquid
      F=0.5
c* calculate the volume of liquid in the cell
      V=F*VT
c* choose the desired unit-length normal vector of the truncation plane
      XNC=0.57735027
      YNC=0.57735027
      ZNC=0.57735027
c* Solve the local volume enforcement problem:
c--------------------------------------------
      WRITE(6,*)'-----------------------------------------------------'
      WRITE(6,*)'|---------------------------------------------------|'
      WRITE(6,*)'|        OUTPUT OF THE ENFORV3D SUBROUTINE          |'
      WRITE(6,*)'|---------------------------------------------------|'
      WRITE(6,*)'-----------------------------------------------------'
      WRITE(6,*)''
      CALL ENFORV3D(C,IPV,NIPV,NTP,NTS,NTV,V,VT,VERTP,XNC,XNS,YNC,YNS,
     -     ZNC,ZNS)
c.. For rectangular parallelepiped cells, like that defined in the 
c.. CUBICMESH subroutine,it is more efficient to use the analytical method 
c.. of Scardovelli and Zaleski [Journal of Computational Physics, 164 
c.. (2000) 228-237], which was proposed to be used specifically for this 
c.. type of cells. To use this method, uncomment the following four lines.
c      DX=DABS(VERTP(1,1)-VERTP(5,1))     !x-side length
c      DY=DABS(VERTP(1,2)-VERTP(4,2))     !y-side length
c      DZ=DABS(VERTP(2,3)-VERTP(1,3))     !z-side length
c      CALL ENFORV3DSZ(C,DX,DY,DZ,V,VERTP,XNC,YNC,ZNC)
      WRITE(6,*)'Solution for C:',C
      WRITE(6,*)''
      WRITE(6,*)'-----------------------------------------------------'
c* Calculate the distance from a given point P to the interfacial polygon:
c------------------------------------------------------------------------
c* copy the original polyhedron that defines the cell to the working
c* polyhedron 0
      CALL CPPOL3D(CS0,CS,IPV0,IPV,NIPV0,NIPV,NTP0,NTP,NTS0,NTS,NTV0,
     -     NTV,VERTP0,VERTP,XNS0,XNS,YNS0,YNS,ZNS0,ZNS)
c* intersection between the cell and the plane defined as X·NC+C=0  
      CALL INTE3D(C,ICONTN,ICONTP,IPV0,NIPV0,NTP0,NTS0,NTV0,VERTP0,XNC,
     -     XNS0,YNC,YNS0,ZNC,ZNS0)
c* choose the desired point P from which the distance is calculated
      XP=0.0
      YP=0.0
      ZP=0.0
c* interfacial polygon defined as the last face of the truncated polyhedron 0
      N=NIPV0(NTS0)
      DO I=1,N
         IP=IPV0(NTS0,I)
         X(I)=VERTP0(IP,1)
         Y(I)=VERTP0(IP,2)
         Z(I)=VERTP0(IP,3)
      END DO
c* calculate the distance
      CALL DIST3D(D,N,X,Y,Z,XP,YP,ZP)
      WRITE(6,*)'-----------------------------------------------------'
      WRITE(6,*)'|---------------------------------------------------|'
      WRITE(6,*)'|         OUTPUT OF THE DIST3D SUBROUTINE           |'
      WRITE(6,*)'|---------------------------------------------------|'
      WRITE(6,*)'-----------------------------------------------------'
      WRITE(6,*)''
      WRITE(6,*)'Distance from P to the interfacial polygon:',D
      WRITE(6,*)''
      WRITE(6,*)'-----------------------------------------------------'
c* Initialize the material volume fraction in the cell:
c-----------------------------------------------------      
      IF(ISHAPE.EQ.1) THEN
         CALL INITF3D(FUNC3D1,IPV,NC,NIPV,NTP,NTS,NTV,TOL,VERTP,VF,
     -        XNS,YNS,ZNS)
      ELSE IF(ISHAPE.EQ.2) THEN
         CALL INITF3D(FUNC3D2,IPV,NC,NIPV,NTP,NTS,NTV,TOL,VERTP,VF,
     -        XNS,YNS,ZNS)
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
      END PROGRAM TEST3D
c-------------------------- END OF TEST3D ----------------------------c
c---------------------------------------------------------------------c

