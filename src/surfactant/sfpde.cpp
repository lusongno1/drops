/// \file
/// \brief Some settings about PDE para, and some global function
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
#include "sfpde.h"
/*
// test case 1
//define right hand side and true solution
//my test case，f = 3*(x+y+z) for problem -\Delta u + u = f
double xyz_rhs (const DROPS::Point3DCL& p, double)
{

    return 3*(p[0]+p[1]+p[2]);//p.norm();
}
//my test case u=x+y+z
double laplace_beltrami_xyz_sol (const DROPS::Point3DCL& p, double)
{
    return (p[0]+p[1]+p[2]);//p.norm();
}
*/

// test case 2
double xyz_rhs (const DROPS::Point3DCL& p, double)
{

    return 1;
}
double laplace_beltrami_xyz_sol (const DROPS::Point3DCL& p, double)
{
    return 1;
}

//define level set funcion and its gradient
//unit ball zero level set
void lsFun(double x, double y, double z, double *value)
{
    *value = x * x + y * y + z * z - 1.0;
}

void lsGrad(double x, double y, double z, double *grad)
/* the gradient of the level set function */
{
    grad[0] = x + x;
    grad[1] = y + y;
    grad[2] = z + z;
}

//**********************some useful funtion********************/

void vecMinus(double a[3],double b[3],double (&result)[3])
{
    for(int i=0; i<3; i++)
    {
        result[i] = a[i]-b[i];
    }
}

void crossMul(double a[3],double b[3],double (&p)[3])
{
    p[0] = a[1]*b[2] - a[2]*b[1];
    p[1] = a[2]*b[0] - a[0]*b[2];
    p[2] = a[0]*b[1] - a[1]*b[0];
}

double dotP3(double a[3],double b[3])
{

    return a[0]*b[0]+a[1]*b[1]+a[2]*b[2];
}


double getBaryCoord(double tetra[4][3],int i,double x,double y,double z)
{
    double pValue = 0;
    double v0[3] = {tetra[i][0],tetra[i][1],tetra[i][2]};
    int idx = 0;
    double vGround[3][3];
    for(int j; j<4; j++)
    {
        if(j==i)
            continue;
        for(int k=0; k<3; k++)
            vGround[idx][k] = tetra[j][k];
        idx++;
    }
    double vec1[3];
    double vec2[3];
    double n[3];
    vecMinus(vGround[1],vGround[0],vec1);
    vecMinus(vGround[2],vGround[0],vec2);
    crossMul(vec1,vec2,n);
    double n_norm = std::sqrt(n[0]*n[0]+n[1]*n[1]+n[2]*n[2]);
    for(int j=0; j<3; j++)
        n[j] /=n_norm;
    double vecV[3] = {v0[0] - vGround[0][0],v0[1] - vGround[0][1],v0[2] - vGround[0][2]};
    double vecX[3] = {x - vGround[0][0],y - vGround[0][1],z - vGround[0][2]};
    double valXYZ = dotP3(vecX,n);
    double valV = dotP3(vecV,n);
    //assert((valV>=0&&valXYZ>=0)||(valV<=0&&valXYZ<=0));
    pValue = valXYZ/valV;
    //assert()
    return pValue;
}



