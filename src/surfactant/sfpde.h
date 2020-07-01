/// \file
/// \brief Including file: Some setting about PDE para, and some global function
/// \author LSEC: Song Lu
/*
 * This file is part of DROPS.
 *
 * DROPS is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * DROPS is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with DROPS. If not, see <http://www.gnu.org/licenses/>.
 *
 *
 * Copyright 2009 LNM/SC RWTH Aachen, Germany
*/
//#pragma once
#include "misc/container.h"
#include "geom/simplex.h"
#include "num/discretize.h"

#ifndef DROPS_SFPDE_H
#define DROPS_SFPDE_H
//namespace DROPS
//{
double xyz_rhs (const DROPS::Point3DCL& p, double);
double laplace_beltrami_xyz_sol (const DROPS::Point3DCL& p, double);
void lsFun(double x, double y, double z, double *value);
void lsGrad(double x, double y, double z, double *grad);
//}

void vecMinus(double a[3],double b[3],double (&result)[3]);
void crossMul(double a[3],double b[3],double (&p)[3]);
double dotP3(double a[3],double b[3]);
double getBaryCoord(double tetra[4][3],int i,double x,double y,double z);
DROPS::BaryCoordCL getBaryCoords(double tetra[4][3],double x,double y,double z);
void GetTet2DArr(const DROPS::TetraCL& t,double tet[4][3]);

extern double tet[4][3];
extern int iG;
extern int jG;
extern int orderG;
extern double gradTri[4][3];//store gradients of shape functions
//extern DROPS::LocalP2CL<> localP2Set[10];
#endif
