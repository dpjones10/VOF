//---------------------------------------------------------------------//
//---------------------------------------------------------------------//
//                             cmesh.h                                 //
//---------------------------------------------------------------------//
//            Copyright (C) 2016 J. Lopez and J. Hernandez             //
//---------------------------------------------------------------------//
//      Declaration of the subroutines used to define the cell         //
//      geometries considered in the tests                             //
//---------------------------------------------------------------------//
// This file is part of VOFTools.                                      //
//                                                                     //
// VOFTools is free software: you can redistribute it and/or           //
// modify it under the terms of the GNU General Public License as      //
// published by the Free Software Foundation, either version 3 of      //
// the License, or (at your option) any later version.                 //
//                                                                     //
// VOFTools is distributed in the hope that it will be useful,         //
// but WITHOUT ANY WARRANTY; without even the implied warranty of      //
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the       //
// GNU General Public License for more details.                        //
//                                                                     //
// You should have received a copy of the GNU General Public License   //
// along with VOFTools.  If not, see <http://www.gnu.org/licenses/>.   //
//---------------------------------------------------------------------//

// 3D cell geometries

extern void cubicmesh_(int *ipv,int *nipv,int *ntp,int *nts,int *ntv,
		       double *vertp,double *xns,double *yns,double *zns);

extern void hexahemesh_(int *ipv,int *nipv,int *ntp,int *nts,int *ntv,
		       double *vertp,double *xns,double *yns,double *zns);

extern void tetramesh_(int *ipv,int *nipv,int *ntp,int *nts,int *ntv,
		       double *vertp,double *xns,double *yns,double *zns);

extern void dodecamesh_(int *ipv,int *nipv,int *ntp,int *nts,int *ntv,
		       double *vertp,double *xns,double *yns,double *zns);

extern void icosamesh_(int *ipv,int *nipv,int *ntp,int *nts,int *ntv,
		       double *vertp,double *xns,double *yns,double *zns);

extern void complexmesh_(int *ipv,int *nipv,int *ntp,int *nts,int *ntv,
		       double *vertp,double *xns,double *yns,double *zns);

// 2D cell geometries

extern void squaremesh_(int *ipv,int *ntp,int *ntv,double *vertp);

extern void hexagomesh_(int *ipv,int *ntp,int *ntv,double *vertp);

extern void trianglemesh_(int *ipv,int *ntp,int *ntv,double *vertp);

extern void quadranglemesh_(int *ipv,int *ntp,int *ntv,double *vertp);

extern void pentagonmesh_(int *ipv,int *ntp,int *ntv,double *vertp);

extern void hexagonmesh_(int *ipv,int *ntp,int *ntv,double *vertp);
