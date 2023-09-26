//---------------------------------------------------------------------//
//---------------------------------------------------------------------//
//                           TEST2D in C                               //
//---------------------------------------------------------------------//
//            Copyright (C) 2019 J. Lopez and J. Hernandez             //
//---------------------------------------------------------------------//
//      Test program in C to solve the local volume enforcement        //
//      problem and to calculate the distance from a given point to    //
//      the interfacial polygon                                        // 
//---------------------------------------------------------------------//
// This file is part of VOFTools.                                      //
//                                                                     //
// VOFTools is free software: you can redistribute it and/or           //
// modify it under the terms of the GNU General Public License         //
// as published by the Free Software Foundation, either version 3 of   //
// the License, or (at your option) any later version.                 //
//                                                                     //
// VOFTools is distributed in the hope that it will be useful,         //
// but WITHOUT ANY WARRANTY; without even the implied warranty of      //
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the       //
// GNU General Public License for more details.                        //
//                                                                     //
// You should have received a copy of the GNU General Public License   //
// along with VOFTools. If not, see <http://www.gnu.org/licenses/>.    //
//---------------------------------------------------------------------//
#include <stdio.h>
#include "dimc.h"
#include "cvoftools.h"
#include "cmesh.h"
#include "cuservoftools.h"

int main()
{
  //* Original polyhedron
  double vertp[nv*2];
  int ipv[nv];
  //* working polyhedron 0
  double vertp0[nv*2];
  int ipv0[nv];
  //* working polyhedron 1
  double vertp1[nv*2];
  int ipv1[nv];
  //* working polyhedron 2
  double vertp2[nv*2];
  int ipv2[nv];
  //* Interfacial polygon
  double x[2],y[2];
  //* other.......
  int i, icontn, icontp, ip, ishape, n, nc, ntp, ntp0, ntv, ntv0;
  double c,d,dx,dy,f,tol,v,vf,vt,xnc,xp,ync,yp;

  printf("-----------------------------------------------------\n");
  printf("|---------------------------------------------------|\n");
  printf("|          TEST PROGRAM IN 2D OF VOFTools           |\n");
  printf("|                                                   |\n");
  printf("|            (Version 3.2, June 2019)               |\n");
  printf("|                                                   |\n");
  printf("|                       by                          |\n");
  printf("|                                                   |\n");
  printf("|            J. Lopez and J. Hernandez              |\n");
  printf("|---------------------------------------------------|\n");
  printf("-----------------------------------------------------\n");
  printf("\n");
  //* Material body shape to be initialized:
  //----------------------------------------
  //. ISHAPE  =  1, circle with radious 0.325 centered at (0.5,0.5)
  //.         =  2, ellipse with semi-major axis 0.6, semi-minor axis 0.2
  //.               and centered at (0.5,0.5)      
  ishape=1;
  //* Subdivision number in the volume fraction cell initialization:
  nc=10;
  //* Tolerance for the initialization procedure:
  tol=10.0;
  //* Choose the desired cell type:
  //------------------------------
  //* Square mesh
  squaremesh_(ipv,&ntp,&ntv,vertp);
  //* Regular hexagonal mesh
  //hexagomesh_(ipv,&ntp,&ntv,vertp);
  //* Triangular mesh
  //trianglemesh_(ipv,&ntp,&ntv,vertp);
  //* Quadrangular mesh
  //quadranglemesh_(ipv,&ntp,&ntv,vertp);
  //* Pentagonal mesh
  //pentagonmesh_(ipv,&ntp,&ntv,vertp);
  //* Irregular hexagonal mesh
  //hexagonmesh_(ipv,&ntp,&ntv,vertp);
  //* Calculate the area VT of the cell:
  //------------------------------------
  toolv2d_(ipv,&ntv,vertp,&vt);
  printf("-----------------------------------------------------\n");
  printf("|---------------------------------------------------|\n");
  printf("|         OUTPUT OF THE TOOLV2D SUBROUTINE          |\n");
  printf("|---------------------------------------------------|\n");
  printf("-----------------------------------------------------\n");
  printf("\n");
  printf("Area of the selected cell:%f\n",vt);
  printf("\n");
  printf("-----------------------------------------------------\n");
  //* choose the desired volume fraction of liquid
  f=0.5;
  //* calculate the volume of liquid in the cell
  v=f*vt;
  //* choose the desired unit-length normal vector of the truncation plane
  xnc=0.70710678;
  ync=0.70710678;
  //* Solve the local volume enforcement problem:
  //--------------------------------------------
  printf("-----------------------------------------------------\n");
  printf("|---------------------------------------------------|\n");
  printf("|        OUTPUT OF THE ENFORV2D SUBROUTINE          |\n");
  printf("|---------------------------------------------------|\n");
  printf("-----------------------------------------------------\n");
  printf("\n");
  //enforv2d_(&c,ipv,&ntp,&ntv,&v,&vt,vertp,&xnc,&ync);
  //* For rectangular cells, like that defined in the squaremesh_ subroutine,
  //* it is more efficient to use the analytical method of Scardovelli and 
  //* Zaleski [Journal of Computational Physics, 164 (2000) 228-237], which
  //* was proposed to be used specifically for this type of cells. To use 
  //* this method, uncomment the following three lines.
  dx=vertp[0*nv+2-1]-vertp[0*nv+1-1];    // x-side length
  dy=vertp[1*nv+3-1]-vertp[1*nv+2-1];    // y-side length
  enforv2dsz_(&c,&dx,&dy,&v,vertp,&xnc,&ync);
  printf("Solution of the problem:%f\n",c);
  printf("\n");
  printf("-----------------------------------------------------\n");
  //* Calculate the distance from a given point P to the interfacial segment:
  //------------------------------------------------------------------------
  //* copy the original polygon that defines the cell to the working
  //* polygon 0
  cppol2d_(ipv,ipv0,&ntp,&ntp0,&ntv,&ntv0,vertp,vertp0);
  //* intersection between the cell and the line defined as X·NC+C=0
  inte2d_(&c,&icontn,&icontp,ipv0,&ntp0,&ntv0,vertp0,&xnc,&ync);
  //* choose the desired point P from which the distance is calculated
  xp=0.0;
  yp=0.0;
  //* interfacial segment defined as the last edge of the truncated polygon 0
  ip=ntp0-1;
  x[0]=vertp0[0*nv+ip-1];
  y[0]=vertp0[1*nv+ip-1];
  ip=ntp0;
  x[1]=vertp0[0*nv+ip-1];
  y[1]=vertp0[1*nv+ip-1];
  //* calculate the distance
  dist2d_(&d,x,y,&xp,&yp);
  printf("-----------------------------------------------------\n");
  printf("|---------------------------------------------------|\n");
  printf("|         OUTPUT OF THE DIST2D SUBROUTINE           |\n");
  printf("|---------------------------------------------------|\n");
  printf("-----------------------------------------------------\n");
  printf("\n");
  printf("Distance from P to the interfacial segment:%f\n",d);
  printf("\n");
  printf("-----------------------------------------------------\n");
  //* Initialize the material volume fraction in the cell:
  //-----------------------------------------------------      
  if(ishape==1){
    initf2d_(func2d1_,ipv,&nc,&ntp,&ntv,&tol,vertp,&vf);
  }
  else if(ishape==2){
    initf2d_(func2d2_,ipv,&nc,&ntp,&ntv,&tol,vertp,&vf);
  }
  printf("-----------------------------------------------------\n");
  printf("|---------------------------------------------------|\n");
  printf("|         OUTPUT OF THE INITF2D SUBROUTINE          |\n");
  printf("|---------------------------------------------------|\n");
  printf("-----------------------------------------------------\n");
  printf("\n");
  printf("Material volume fraction in the selected cell:%f\n",vf);
  printf("\n");
  printf("-----------------------------------------------------\n");
} 
