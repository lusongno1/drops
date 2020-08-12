/// \file
/// \brief Some settings about PDE paras, and some global functions
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

//// test case 1
////define right hand side and true solution
////my test case，f = 3*(x+y+z) for problem -\Delta u + u = f
////then u = f/3
double xyz_rhs (const DROPS::Point3DCL& p, double)
{

    return 3*(p[0]+p[1]+p[2]);//p.norm();
}
//my test case u=x+y+z
double laplace_beltrami_xyz_sol (const DROPS::Point3DCL& p, double)
{
    return (p[0]+p[1]+p[2]);//p.norm();
}

DROPS::Point3DCL laplace_beltrami_xyz_sol_grad (const DROPS::Point3DCL& p, double)
{
    DROPS::Point3DCL tmp{1,1,1};
    return tmp;
}


// test case 2
//double xyz_rhs (const DROPS::Point3DCL& p, double)
//{
//    return 1;
//}
//double laplace_beltrami_xyz_sol (const DROPS::Point3DCL& p, double)
//{
//    return 1;
//}
//
//DROPS::Point3DCL laplace_beltrami_xyz_sol_grad (const DROPS::Point3DCL& p, double)
//{
//    DROPS::Point3DCL tmp{0,0,0};
//    return tmp;
//}

// test case 3
//define right hand side and true solution
//u = a*|x|^2/(12+|x|^2)*(3x1^2x2-x2^3)
//f = a*(3x1^2x2-x2^3)
//double a(1.0);
//double xyz_rhs (const DROPS::Point3DCL& p, double)
//{
//
//    return a/std::pow( p.norm(), 3.)*(3.*p[0]*p[0]*p[1]-p[1]*p[1]*p[1]);
//}
//double laplace_beltrami_xyz_sol (const DROPS::Point3DCL& p, double)
//{
//    return (p.norm_sq()/(12.+p.norm_sq()))*xyz_rhs(p,0.);
//}
//DROPS::Point3DCL laplace_beltrami_xyz_sol_grad (const DROPS::Point3DCL& p, double)
//{
// //   DROPS::Point3DCL tmp{6*a/13*p[0]*p[1],-3*a/13*p[1]*p[1],0};
//    DROPS::Point3DCL tmp= 3./std::pow( p.norm(), 3)
//    *( DROPS::MakePoint3D(2.*p[0]*p[1], p[0]*p[0] - p[1]*p[1], 0.) -
//      (3.*p[0]*p[0]*p[1] - std::pow(p[1], 3))/p.norm_sq()*p);
//    return tmp;// This equals tmp - inner_prod( p/p.norm(), tmp)*p/p.norm().
//}

//test case 4
//define right hand side and true solution
//u = x*y*z
//f = x*y*z - 12*x*y*z*(x^2 + y^2 + z^2 - 2)
//double xyz_rhs (const DROPS::Point3DCL& p, double)
//{
//
//    return p[0]*p[1]*p[2]-12.*p[0]*p[1]*p[2]*(p.norm_sq()-2);
//}
//double laplace_beltrami_xyz_sol (const DROPS::Point3DCL& p, double)
//{
//    return p[0]*p[1]*p[2];
//}
//
//DROPS::Point3DCL laplace_beltrami_xyz_sol_grad (const DROPS::Point3DCL& p, double)
//{
//    DROPS::Point3DCL tmp{p[1]*p[2],p[0]*p[2],p[0]*p[1]};
//    return tmp;
//}

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
//template <typename T>
//double dotP3(T a,T b)
//{
//
//    return a[0]*b[0]+a[1]*b[1]+a[2]*b[2];
//}


double getBaryCoord(double tetra[4][3],int i,double x,double y,double z)
{
    double pValue = 0;
    double v0[3] = {tetra[i][0],tetra[i][1],tetra[i][2]};
    int idx = 0;
    double vGround[3][3];
    for(int j=0; j<4; j++)
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

DROPS::BaryCoordCL getBaryCoords(double tetra[4][3],double x,double y,double z)
{
    double BaryCoordArr[4];
    for(int i=0; i<4; i++)
    {
        BaryCoordArr[i] = getBaryCoord(tet,i, x,y,z);
    }
    const DROPS::BaryCoordCL& BaryCoord{BaryCoordArr[0],BaryCoordArr[1],BaryCoordArr[2],BaryCoordArr[3]};
    return BaryCoord;
}


void GetTet2DArr(const DROPS::TetraCL& t,double tet[4][3])
{
    for (int i= 0; i < 4; ++i)
    {
        auto vtx = t.GetVertex(i);
        auto coord = vtx->GetCoord();
        for(int j=0; j<3; j++)
        {

            tet[i][j] = coord[j];

        }
    }

}


void getSfNormalVec(double x,double y,double z,double (&n)[3])
{
    double ls_grad[3];
    lsGrad(x,y,z,ls_grad);
    double ls_grad_norm = std::sqrt(ls_grad[0]*ls_grad[0]+ls_grad[1]*ls_grad[1]+ls_grad[2]*ls_grad[2]);
    for(int i=0; i<3; i++)
    {
        n[i] = ls_grad[i]/ls_grad_norm;
    }

}

//template <typename T>
void getSurfaceGradient(DROPS::Point3DCL v,double n[3],double (&sf_grad)[3])
{
    double proj_norm = 0;
    for(int i=0; i<3; i++)
        proj_norm += v[i]*n[i];
    for(int i=0; i<3; i++)
    {
        sf_grad[i] = v[i] - proj_norm*n[i];
    }
}

void getSurfaceGradient(DROPS::Point3DCL v,double n[3],DROPS::SVectorCL<3> &sf_grad)
{
    double proj_norm = 0;
    for(int i=0; i<3; i++)
        proj_norm += v[i]*n[i];
    for(int i=0; i<3; i++)
    {
        sf_grad[i] = v[i] - proj_norm*n[i];
    }
}


void ouput_valarray(std::valarray<double> v)
{
    std::cout<<"begin output valarray:"<<std::endl;
    for(int i=0; i<v.size(); i++)
    {
        std::cout<<v[i]<<" ";
    }
    std::cout<<std::endl;
}
void cout2txt(double a)
{
    std::ofstream mycout("./debug.txt",std::ios_base::app);
    mycout<<a<<std::endl;
    mycout.close();
}

int nc = 1;
void coutTet(const DROPS::TetraCL& t)
{
    std::cout<<std::endl<<nc++<<":"<<std::endl;;
    for (int i= 0; i < 4; ++i)
    {
        auto vtx = t.GetVertex(i);
        auto coord = vtx->GetCoord();
        for(int j=0; j<3; j++)
        {

            auto tmp = coord[j];
            std::cout<<tmp<<" ";

        }
        std::cout<<std::endl;
    }

}

double tet[4][3];
int iG;
int jG;
int orderG = 10;
double gradTri[4][3];
//DROPS::LocalP2CL<> localP2Set[10];


