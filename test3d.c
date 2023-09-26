//---------------------------------------------------------------------//
//---------------------------------------------------------------------//
//                              TEST3D                                 //
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
  double cs[ns],vertp[nv*3],xns[ns],yns[ns],zns[ns];
  int ipv[ns*nv],nipv[ns];
  //* working polyhedron 0
  double cs0[ns],vertp0[nv*3],xns0[ns],yns0[ns],zns0[ns];
  int ipv0[ns*nv],nipv0[ns];
  //* working polyhedron 1
  double cs1[ns],vertp1[nv*3],xns1[ns],yns1[ns],zns1[ns];
  int ipv1[ns*nv],nipv1[ns];
  //* working polyhedron 2
  double cs2[ns],vertp2[nv*3],xns2[ns],yns2[ns],zns2[ns];
  int ipv2[ns*nv],nipv2[ns];
  //* Interfacial polygon
  double x[nv],y[nv],z[nv];
  //* other.......
  int i, icontn, icontp, ip, ishape, n, nc, ntp, ntp0, nts, nts0, ntv, ntv0;
  double c,d,dx,dy,dz,f,tol,v,vf,vt,xnc,xp,ync,yp,znc,zp;

  printf("-----------------------------------------------------\n");
  printf("|---------------------------------------------------|\n");
  printf("|          TEST PROGRAM IN 3D OF VOFTools           |\n");
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
  //. ISHAPE  =  1, sphere with radious 0.6 centered at (0.5,0.5,0.5)
  //.         =  2, Torus with major radius 2/3, minor radius 1/3 and
  //.               centered at (0.5,0.5,0.5)
  ishape=1;
  //* Subdivision number in the volume fraction cell initialization:
  nc=10;
  //* Tolerance for the initialization procedure:
  tol=10.0;
  //* Choose the desired cell type:
  //------------------------------
  //* Cubic mesh
  cubicmesh_(ipv,nipv,&ntp,&nts,&ntv,vertp,xns,yns,zns);
  //* general hexahedrical mesh
  //hexahemesh_(ipv,nipv,&ntp,&nts,&ntv,vertp,xns,yns,zns);
  //* Tetrahedrical mesh
  //tetramesh_(ipv,nipv,&ntp,&nts,&ntv,vertp,xns,yns,zns);
  //* dodecahedral mesh
  //dodecamesh_(ipv,nipv,&ntp,&nts,&ntv,vertp,xns,yns,zns);
  //* icosahedral mesh
  //icosamesh_(ipv,nipv,&ntp,&nts,&ntv,vertp,xns,yns,zns);
  //* complex mesh
  //complexmesh_(ipv,nipv,&ntp,&nts,&ntv,vertp,xns,yns,zns);
  //* Calculate the volume VT of the cell:
  //-------------------------------------
  toolv3d_(ipv,nipv,&nts,vertp,&vt,xns,yns,zns);
  printf("-----------------------------------------------------\n");
  printf("|---------------------------------------------------|\n");
  printf("|         OUTPUT OF THE TOOLV3D SUBROUTINE          |\n");
  printf("|---------------------------------------------------|\n");
  printf("-----------------------------------------------------\n");
  printf("\n");
  printf("Volume of the selected cell:%f\n",vt);
  printf("\n");
  printf("-----------------------------------------------------\n");
  //* choose the desired volume fraction of liquid
  f=0.5;
  //* calculate the volume of liquid in the cell
  v=f*vt;
  //* choose the desired unit-length normal vector of the truncation plane
  xnc=0.57735027;
  ync=0.57735027;
  znc=0.57735027;
  //* Solve the local volume enforcement problem:
  //--------------------------------------------
  printf("-----------------------------------------------------\n");
  printf("|---------------------------------------------------|\n");
  printf("|        OUTPUT OF THE ENFORV3D SUBROUTINE          |\n");
  printf("|---------------------------------------------------|\n");
  printf("-----------------------------------------------------\n");
  printf("\n");
  enforv3d_(&c,ipv,nipv,&ntp,&nts,&ntv,&v,&vt,vertp,&xnc,xns,&ync,yns,
	    &znc,zns);
  //* For rectangular parallelepiped cells, like that defined in the 
  //* CUBICMESH subroutine,it is more efficient to use the analytical method 
  //* of Scardovelli and Zaleski [Journal of Computational Physics, 164 
  //* (2000) 228-237], which was proposed to be used specifically for this 
  //* type of cells. To use this method, uncomment the following four lines.
  //dx=vertp[0*nv+1-1]-vertp[0*nv+5-1];     // x-side length
  //dy=vertp[1*nv+4-1]-vertp[1*nv+1-1];     // y-side length
  //dz=vertp[2*nv+1-1]-vertp[2*nv+2-1];     // z-side length
  //enforv3dsz_(&c,&dx,&dy,&dz,&v,vertp,&xnc,&ync,&znc);
  printf("Solution for c:%f\n",c);
  printf("\n");
  printf("-----------------------------------------------------\n");
  //* Calculate the distance from a given point P to the interfacial polygon:
  //------------------------------------------------------------------------
  //* copy the original polyhedron that defines the cell to the working
  //* polyhedron 0
  cppol3d_(cs0,cs,ipv0,ipv,nipv0,nipv,&ntp0,&ntp,&nts0,&nts,&ntv0,
	   &ntv,vertp0,vertp,xns0,xns,yns0,yns,zns0,zns);
  //c* intersection between the cell and the plane defined as X·NC+C=0  
  inte3d_(&c,&icontn,&icontp,ipv0,nipv0,&ntp0,&nts0,&ntv0,vertp0,&xnc,
	 xns0,&ync,yns0,&znc,zns0);
  //* choose the desired point P from which the distance is calculated
  xp=0.0;
  yp=0.0;
  zp=0.0;
  //* interfacial polygon defined as the last face of the truncated polyhedron 0
  n=nipv0[nts0-1];
  for(i=0; i<n; ++i){
    ip=ipv0[i*ns+nts0-1];
    x[i]=vertp0[0*nv+ip-1];
    y[i]=vertp0[1*nv+ip-1];
    z[i]=vertp0[2*nv+ip-1];
  }
  //* calculate the distance
  dist3d_(&d,&n,x,y,z,&xp,&yp,&zp);
  printf("-----------------------------------------------------\n");
  printf("|---------------------------------------------------|\n");
  printf("|         OUTPUT OF THE DIST3D SUBROUTINE           |\n");
  printf("|---------------------------------------------------|\n");
  printf("-----------------------------------------------------\n");
  printf("\n");
  printf("Distance from P to the interfacial polygon:%f\n",d);
  printf("\n");
  printf("-----------------------------------------------------\n");
  //* Initialize the material volume fraction in the cell:
  //-----------------------------------------------------      
  if(ishape==1){
    initf3d_(func3d1_,ipv,&nc,nipv,&ntp,&nts,&ntv,&tol,vertp,&vf,xns,yns,zns);
  }
  else if(ishape==2){
    initf3d_(func3d2_,ipv,&nc,nipv,&ntp,&nts,&ntv,&tol,vertp,&vf,xns,yns,zns);
  }
  printf("-----------------------------------------------------\n");
  printf("|---------------------------------------------------|\n");
  printf("|         OUTPUT OF THE INITF3D SUBROUTINE          |\n");
  printf("|---------------------------------------------------|\n");
  printf("-----------------------------------------------------\n");
  printf("\n");
  printf("Material volume fraction in the selected cell:%f\n",vf);
  printf("\n");
  printf("-----------------------------------------------------\n");
} 
