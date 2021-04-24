/// \file
/// \brief Some settings about PDE paras, and some global variables, functions, and classes
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
#include <cmath>
#include "misc/funcmap.h"
#include "num/accumulator.h"
#include "geom/principallattice.h"
#include "num/lattice-eval.h"
#include "geom/subtriangulation.h"
using namespace DROPS;
//**********************************************set global variables*******************************************************/
double tet[4][3];
int iG;
int jG;
int orderG = 10;
double gradTri[4][3];

//***************************************************define test case*******************************************************/
static DROPS::RegisterScalarFunction regsca_level_set_function_drops( "LevelSetFunDrops", level_set_function_drops);
static DROPS::RegisterScalarFunction regsca_xyz_rhs( "xyzRhs", xyz_rhs);
static DROPS::RegisterScalarFunction regsca_laplace_beltrami_xyz_sol( "LaplaceBeltramixyzSol", laplace_beltrami_xyz_sol);


// test case 1
//define right hand side and true solution
//my test case，f = 3*(x+y+z) for problem -\Delta u + u = f
//then u = f/3
#if 0
double xyz_rhs (const DROPS::Point3DCL& p, double)
{

    return 3*(p[0]+p[1]+p[2]);//p.norm();
    //return 3*(p[0]+p[1]+p[2])/p.norm();
}
//my test case u=x+y+z
double laplace_beltrami_xyz_sol (const DROPS::Point3DCL& p, double)
{
    return (p[0]+p[1]+p[2]);//p.norm();
    //return (p[0]+p[1]+p[2])/p.norm();
}

DROPS::Point3DCL laplace_beltrami_xyz_sol_grad (const DROPS::Point3DCL& p, double)
{
    DROPS::Point3DCL tmp{1,1,1};
    return tmp;
}
#endif


// test case 2, constant
#if 0
double xyz_rhs (const DROPS::Point3DCL& p, double)
{
    return 1;
}
double laplace_beltrami_xyz_sol (const DROPS::Point3DCL& p, double)
{
    return 1;
}

DROPS::Point3DCL laplace_beltrami_xyz_sol_grad (const DROPS::Point3DCL& p, double)
{
    DROPS::Point3DCL tmp{0,0,0};
    return tmp;
}
#endif


// test case 3, a general test case
//define right hand side and true solution
//u = a*|x|^2/(12+|x|^2)*(3x1^2x2-x2^3)
//f = a*(3x1^2x2-x2^3)
#if 1
double level_set_function_drops(const DROPS::Point3DCL& p, double)//directly modified in routine
{

    double x = p[0],y=p[1],z=p[2];
    return x * x + y * y + z * z - 1.0;
    //return p.norm()-1.0;
}


void lsFun(double x, double y, double z, double *value)
{
    *value = x * x + y * y + z * z - 1.0;
}


void lsGrad(double x, double y, double z, double *grad)
///* the gradient of the level set function */
{
    grad[0] = x + x;
    grad[1] = y + y;
    grad[2] = z + z;
}
double a(1.0);
double xyz_rhs (const DROPS::Point3DCL& p, double)
{
    return 3.0*p[0]*p[0]*p[1]-p[1]*p[1]*p[1];
    //return a/std::pow( p.norm(), 3.)*(3.*p[0]*p[0]*p[1]-p[1]*p[1]*p[1]);
}
double laplace_beltrami_xyz_sol (const DROPS::Point3DCL& p, double)
{
    return (p.norm_sq()/(12.+p.norm_sq()))*xyz_rhs(p,0.);
}
DROPS::Point3DCL laplace_beltrami_xyz_sol_grad (const DROPS::Point3DCL& p, double)
{
//   DROPS::Point3DCL tmp{6*a/13*p[0]*p[1],-3*a/13*p[1]*p[1],0};
    DROPS::Point3DCL tmp= 3./std::pow( p.norm(), 3)
                          *( DROPS::MakePoint3D(2.*p[0]*p[1], p[0]*p[0] - p[1]*p[1], 0.) -
                             (3.*p[0]*p[0]*p[1] - std::pow(p[1], 3))/p.norm_sq()*p);
    return tmp;// This equals tmp - inner_prod( p/p.norm(), tmp)*p/p.norm().
}

#endif


//test case 4
//define right hand side and true solution
//u = x*y*z
//f = x*y*z - 12*x*y*z*(x^2 + y^2 + z^2 - 2)
#if 0
double xyz_rhs (const DROPS::Point3DCL& p, double)
{

    return p[0]*p[1]*p[2]-12.*p[0]*p[1]*p[2]*(p.norm_sq()-2);
}
double laplace_beltrami_xyz_sol (const DROPS::Point3DCL& p, double)
{
    return p[0]*p[1]*p[2];
}


DROPS::Point3DCL laplace_beltrami_xyz_sol_grad (const DROPS::Point3DCL& p, double)
{
    DROPS::Point3DCL tmp{p[1]*p[2],p[0]*p[2],p[0]*p[1]};
    return tmp;
}

double level_set_function_drops (const DROPS::Point3DCL& p, double)
{
    DROPS::Point3DCL RadDrop(1,1,1);
    DROPS::Point3DCL PosDrop(0,0,0);
    DROPS::Point3DCL x( p - PosDrop);
    //double value=0;
    //lsFun(x[0], x[1], x[2], &value);
    return x.norm() - RadDrop[0];
    //return value;
}
static DROPS::RegisterScalarFunction regsca_sphere_dist_lset( "LevelSetFunDrops", level_set_function_drops);

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
#endif

//test for the gyroid, heart-shape and torus, just have a thy
#if 0
double level_set_function_drops(const DROPS::Point3DCL& p, double)//directly modified in routine
{
    //std::cout<<M_PI<<std::endl;
    //double phi = cos(M_PI*p[0])*sin(M_PI*p[1])+
    //cos(M_PI*p[1])*sin(M_PI*p[2])+cos(M_PI*p[2])*sin(M_PI*p[0]);
    double x = p[0],y=p[1],z=p[2];
    //double phi = 2*pow(x-0.5,2)-8*pow(y-0.5,3)-16*pow(z-0.5,4)-1/50;
    //double phi = pow((pow(x,2)+(9/4)*pow(y,2)+pow(z,2)-1),3)-pow(x,2)*pow(z,3)-(9/80)*pow(y,2)*pow(z,3);
    double phi = sqrt(pow(sqrt(x*x+y*y)-1.0,2.0)+z*z)-3.0/5.0;;
    return phi;
//    return std::sqrt( std::pow( RadTorus[0] - std::sqrt(p[0]*p[0] + p[1]*p[1]), 2) + std::pow( p[2], 2)) - RadTorus[1];
}
static DROPS::RegisterScalarFunction regsca_sphere_dist_lset( "LevelSetFunDrops", level_set_function_drops);

//void lsFun(double x, double y, double z, double *value)
//{
//    *value = x * x + y * y + z * z - 1.0;
//}


//void lsGrad(double x, double y, double z, double *grad)
///* the gradient of the level set function */
//{
//    grad[0] = x + x;
//    grad[1] = y + y;
//    grad[2] = z + z;
//}

void lsFun(double x, double y, double z, double *value)
{
    double R = 1;
    double r = 0.6;
    //double phi = pow((pow(x,2)+(9/4)*pow(y,2)+pow(z,2)-1),3)-pow(x,2)*pow(z,3)-(9/80)*pow(y,2)*pow(z,3);
    double phi = sqrt(pow(sqrt(x*x+y*y)-1.0,2.0)+z*z)-3.0/5.0;;
    *value = phi;//2*pow(x-0.5,2)-8*pow(y-0.5,3)-16*pow(z-0.5,4)-1/50;
    //cos(M_PI*x)*sin(M_PI*y)+cos(M_PI*y)*sin(M_PI*z)+cos(M_PI*z)*sin(M_PI*x);
}
//
void lsGrad(double x, double y, double z, double *grad)
/* the gradient of the level set function */
{
    grad[0] = x*1.0/sqrt(x*x+y*y)*1.0/sqrt(pow(sqrt(x*x+y*y)-1.0,2.0)+z*z)*(sqrt(x*x+y*y)-1.0);
    grad[1] = y*1.0/sqrt(x*x+y*y)*1.0/sqrt(pow(sqrt(x*x+y*y)-1.0,2.0)+z*z)*(sqrt(x*x+y*y)-1.0);
    grad[2] = z*1.0/sqrt(pow(sqrt(x*x+y*y)-1.0,2.0)+z*z);

    //grad[0] = 6*x*pow(pow(x,2) + (9*pow(y,z))/4 + pow(z,2) - 1,2) - 2*x*pow(z,3);
    //grad[1] =      (27*y*pow(pow(x,2) + (9*pow(y,z))/4 + pow(z,2) - 1,2))/2 - (9*y*pow(z,3))/40;
    //grad[2] = 6*z*pow(pow(x,2) + (9*pow(y,z))/4 + pow(z,2) - 1,2) - (27*pow(y,z)*pow(z,2))/80 - 3*pow(x,2)*pow(z,2);
    //grad[0] =   4*x - 2;
    //grad[1]  = -24*pow(y - 1/2,2);
    //grad[2] = -64*pow(z - 1/2,3);
    // grad[0] = M_PI*cos(M_PI*x)*cos(M_PI*z) - M_PI*sin(M_PI*x)*sin(M_PI*y);
    // grad[1] = M_PI*cos(M_PI*x)*cos(M_PI*y) - M_PI*sin(M_PI*y)*sin(M_PI*z);
    // grad[2] = M_PI*cos(M_PI*y)*cos(M_PI*z) - M_PI*sin(M_PI*x)*sin(M_PI*z);
}

double xyz_rhs (const DROPS::Point3DCL& p, double)
{
    double x = p[0],y=p[1],z=p[2];
    //return pow(x,2)*sin(y)*exp(z);

    //return p[0]*p[1]*p[2]-12.*p[0]*p[1]*p[2]*(p.norm_sq()-2);
    //double tmp = x+y+z;
    //return std::sqrt( p[2]*p[2] + std::pow( std::sqrt( p[0]*p[0] + p[1]*p[1]) - 1.0, 2));
    return x*sin(y)*exp(z);


    //return tmp;
}

double laplace_beltrami_xyz_sol (const DROPS::Point3DCL& p, double)
{
    return 1.0;
}

DROPS::Point3DCL laplace_beltrami_xyz_sol_grad (const DROPS::Point3DCL& p, double)
{
    DROPS::Point3DCL tmp{1.0,1.0,1.0};
    return tmp;
}
#endif



//****************************************************some useful funtions ***************************************************************/

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

//****************************************************some useful classed ***************************************************************/
//DROPS::LocalP2CL<> localP2Set[10];
//std::vector<double,4> v4type;
using v3 = DROPS::SVectorCL<3>;
using v4 = DROPS::SVectorCL<4>;
using v43 = DROPS::SMatrixCL<4,3>;
class TetrahedronFECL
{
private:
    v43 coordinates;//(4,std::vector<double>(3));
    v4 baryCoordTmp;
public:

    TetrahedronFECL(v43 coordinates):coordinates(coordinates) {};
    void getBaryCoord(double x,double y,double z)
    {
        for(DROPS::Uint i=0; i<4; i++)
        {
            double pValue = 0;
            double v0[3] = {coordinates(i,0),coordinates(i,1),coordinates(i,2)};
            int idx = 0;
            double vGround[3][3];
            for(int j=0; j<4; j++)
            {
                if(j==i)
                    continue;
                for(int k=0; k<3; k++)
                    vGround[idx][k] = coordinates(j,k);
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
            baryCoordTmp[i] = valXYZ/valV;
            //assert()
        }
    };
};

class TetrahedronP3FECL:public TetrahedronFECL
{
};


//**************************************************************************
// Class:   FE_P3CL                                                        *
// Purpose: Shape functions and their gradients for piecewise quadratic,   *
//          continuous finite elements on the reference tetrahedron        *
//          The number of the H-functions refers to the number of the      *
//          (mid-) vertex in the tetrahedron as defined in topo.h, where   *
//          the degree of freedom is located.                              *
//**************************************************************************
class FE_P3CL
{
  private:
    static const double _D2H[20][3][3];

  public:
    // default ctor, copy-ctor, assignment-op, dtor

    static const Uint NumDoFC= 20;

    // restriction of the shape functions to reference edge
    static double H0(double v1) { return 1. +v1*(2.*v1 -3.); }
    static double H1(double v1) { return v1*(2.*v1 -1.); }
    static double H2(double v1) { return 4.*v1*(1. -v1); }

    // restriction of the shape functions to reference face
    static double H0(double v1, double v2) { const double sum= v1 + v2; return 1. +sum*(2.*sum -3.); }
    static double H1(double v1, double)    { return v1*(2.*v1 -1.); }
    static double H2(double, double v2)    { return v2*(2.*v2 -1.); }
    static double H3(double v1, double v2) { return 4.*v1*( 1. -(v1 + v2) ); }
    static double H4(double v1, double v2) { return 4.*v2*( 1. -(v1 + v2) ); }
    static double H5(double v1, double v2) { return 4.*v1*v2; }

    // restriction of the shape functions to reference tetrahedron
    static double H0(double v1, double v2, double v3) { const double sum= v1 + v2 + v3; return 1. +sum*(2.*sum -3.); }
    static double H1(double v1, double, double)       { return v1*(2.*v1 -1.); }
    static double H2(double, double v2, double)       { return v2*(2.*v2 -1.); }
    static double H3(double, double, double v3)       { return v3*(2.*v3 -1.); }
    static double H4(double v1, double v2, double v3) { return 4.*v1*( 1. -(v1 + v2 + v3) ); }
    static double H5(double v1, double v2, double v3) { return 4.*v2*( 1. -(v1 + v2 + v3) ); }
    static double H6(double v1, double v2, double)    { return 4.*v1*v2; }
    static double H7(double v1, double v2, double v3) { return 4.*v3*( 1. -(v1 + v2 + v3) ); }
    static double H8(double v1, double, double v3)    { return 4.*v1*v3; }
    static double H9(double, double v2, double v3)    { return 4.*v2*v3; }
    static inline double H (Uint dof, double v1, double v2, double v3);

    // restriction of the shape functions to reference tetrahedron, barycentric coordinates
    static inline double H0(const BaryCoordCL& p) { return p[0]*(2.*p[0] - 1.); }
    static inline double H1(const BaryCoordCL& p) { return p[1]*(2.*p[1] - 1.); }
    static inline double H2(const BaryCoordCL& p) { return p[2]*(2.*p[2] - 1.); }
    static inline double H3(const BaryCoordCL& p) { return p[3]*(2.*p[3] - 1.); }
    static inline double H4(const BaryCoordCL& p) { return 4.*p[0]*p[1]; }
    static inline double H5(const BaryCoordCL& p) { return 4.*p[0]*p[2]; }
    static inline double H6(const BaryCoordCL& p) { return 4.*p[1]*p[2]; }
    static inline double H7(const BaryCoordCL& p) { return 4.*p[0]*p[3]; }
    static inline double H8(const BaryCoordCL& p)    { return 4.*p[1]*p[3]; }
    static inline double H9(const BaryCoordCL& p)    { return 4.*p[2]*p[3]; }
    static inline double H (Uint dof, const BaryCoordCL& p);

    template <class Cont>
      static inline typename ValueHelperCL<Cont>::value_type
      val(const Cont& c, const BaryCoordCL& p) {
          return LinearCombinationCL<Cont, typename ValueHelperCL<Cont>::value_type>::do_it( c, H0( p), H1( p), H2( p), H3( p), H4( p),  H5( p), H6( p), H7( p), H8( p), H9( p));
      }
    template <class Cont>
      static inline typename ValueHelperCL<Cont>::value_type
      val(const Cont& c, double v1, double v2, double v3) {
          return c[0] * H0( v1, v2, v3) + c[1] * H1( v1, v2, v3) + c[2] * H2( v1, v2, v3) + c[3] * H3( v1, v2, v3)
               + c[4] * H4( v1, v2, v3) + c[5] * H5( v1, v2, v3) + c[6] * H6( v1, v2, v3) + c[7] * H7( v1, v2, v3)
               + c[8] * H8( v1, v2, v3) + c[9] * H9( v1, v2, v3);
      }

    // pt[0]...pt[numpt-1] are coordinates where the shape-functions are evaluated.
    // v is an array of 20 valarrays. They are resized to have numpt components.
    // v[i] contains H_i( pt[0])...H_i( pt[numpt-1])
    static void ApplyAll(Uint numpt, const BaryCoordCL* const pt, std::valarray<double>* v);

    // gradients of the shape functions on the reference tetrahedron.
    // To obtain the gradient on tetra T: See comments in FE_P1CL.
    static inline SVectorCL<3> DH0Ref(double, double, double);
    static inline SVectorCL<3> DH1Ref(double, double, double);
    static inline SVectorCL<3> DH2Ref(double, double, double);
    static inline SVectorCL<3> DH3Ref(double, double, double);
    static inline SVectorCL<3> DH4Ref(double, double, double);
    static inline SVectorCL<3> DH5Ref(double, double, double);
    static inline SVectorCL<3> DH6Ref(double, double, double);
    static inline SVectorCL<3> DH7Ref(double, double, double);
    static inline SVectorCL<3> DH8Ref(double, double, double);
    static inline SVectorCL<3> DH9Ref(double, double, double);

    // DHREef(i, ...) == DHiRef(...)
    static inline SVectorCL<3> DHRef(Uint dof, double v1, double v2, double v3);

    // D2HRef(i) == second derivative of H_i; this is constant
    // diff( H_d(x), k,i) = sum( M_ij*M_kl*D2HRef(d, j, l), j,l=0..2)
    static inline double D2HRef(Uint dof, Uint r, Uint s)
        { return _D2H[dof][r][s]; }

    // Laplace(H_d) on a tetrahedron T; M:= transpose(inverse(A)), where A is the matrix
    // of the affine transformation that maps the reference tetrahedron onto T.
    static inline double Laplace(Uint dof, const SMatrixCL<3,3>& M);

    // The barycentric coordinates of the dofs.
    static const BaryCoordCL bary_coord[20];
};

//**************************************************************************
// Class:   FE_P3CL                                                        *
//**************************************************************************
inline double
FE_P3CL::H(Uint dof, double v1, double v2, double v3)
{
    switch(dof) {
      case 0: return H0( v1, v2, v3);
      case 1: return H1( v1, v2, v3);
      case 2: return H2( v1, v2, v3);
      case 3: return H3( v1, v2, v3);
      case 4: return H4( v1, v2, v3);
      case 5: return H5( v1, v2, v3);
      case 6: return H6( v1, v2, v3);
      case 7: return H7( v1, v2, v3);
      case 8: return H8( v1, v2, v3);
      case 9: return H9( v1, v2, v3);
      default: throw DROPSErrCL("FE_P3CL::H: Invalid shape function.");
    };
}

inline double
FE_P3CL::H(Uint dof, const BaryCoordCL& p)
{
    switch(dof) {
      case 0: return H0( p);
      case 1: return H1( p);
      case 2: return H2( p);
      case 3: return H3( p);
      case 4: return H4( p);
      case 5: return H5( p);
      case 6: return H6( p);
      case 7: return H7( p);
      case 8: return H8( p);
      case 9: return H9( p);
      default: throw DROPSErrCL("FE_P3CL::H: Invalid shape function.");
    };
}

inline SVectorCL<3>
FE_P3CL::DH0Ref(double v1, double v2, double v3)
{
    SVectorCL<3> ret(4.*(v1 + v2 + v3) -3.);
    return ret;
}
inline SVectorCL<3>
FE_P3CL::DH1Ref(double v1, double, double)
{
    SVectorCL<3> ret;
    ret[0]= 4.*v1 -1.; ret[1]= ret[2]= 0.;
    return ret;
}
inline SVectorCL<3>
FE_P3CL::DH2Ref(double, double v2, double)
{
    SVectorCL<3> ret;
    ret[0]= ret[2]= 0.; ret[1]= 4*v2 -1.;
    return ret;
}
inline SVectorCL<3>
FE_P3CL::DH3Ref(double, double, double v3)
{
    SVectorCL<3> ret;
    ret[0]= 0.; ret[1]= 0.; ret[2]= 4.*v3 -1.;
    return ret;
}
inline SVectorCL<3>
FE_P3CL::DH4Ref(double v1, double v2, double v3)
{
    SVectorCL<3> ret;
    ret[0]= 4.*( 1. - (2.*v1 + v2 + v3) ); ret[1]= ret[2]= -4.*v1;
    return ret;
}
inline SVectorCL<3>
FE_P3CL::DH5Ref(double v1, double v2, double v3)
{
    SVectorCL<3> ret;
    ret[0]= ret[2]= -4.*v2; ret[1]= 4.*( 1. -(2.*v2 + v1 + v3) );
    return ret;
}
inline SVectorCL<3>
FE_P3CL::DH6Ref(double v1, double v2, double)
{
    SVectorCL<3> ret;
    ret[0]= 4.*v2; ret[1]= 4.*v1; ret[2]= 0.;
    return ret;
}
inline SVectorCL<3>
FE_P3CL::DH7Ref(double v1, double v2, double v3)
{
    SVectorCL<3> ret;
    ret[0]= ret[1]= -4.*v3; ret[2]= 4.*( 1. -(2.*v3 + v1 + v2) );
    return ret;
}
inline SVectorCL<3>
FE_P3CL::DH8Ref(double v1, double, double v3)
{
    SVectorCL<3> ret;
    ret[0]= 4.*v3; ret[1]= 0.; ret[2]= 4.*v1;
    return ret;
}
inline SVectorCL<3>
FE_P3CL::DH9Ref(double, double v2, double v3)
{
    SVectorCL<3> ret;
    ret[0]= 0.; ret[1]= 4.*v3; ret[2]= 4.*v2;
    return ret;
}
inline SVectorCL<3>
FE_P3CL::DHRef(Uint dof, double v1, double v2, double v3)
{
    switch (dof)
    {
      case 0: return DH0Ref(v1, v2, v3);
      case 1: return DH1Ref(v1, v2, v3);
      case 2: return DH2Ref(v1, v2, v3);
      case 3: return DH3Ref(v1, v2, v3);
      case 4: return DH4Ref(v1, v2, v3);
      case 5: return DH5Ref(v1, v2, v3);
      case 6: return DH6Ref(v1, v2, v3);
      case 7: return DH7Ref(v1, v2, v3);
      case 8: return DH8Ref(v1, v2, v3);
      case 9: return DH9Ref(v1, v2, v3);
      default: throw DROPSErrCL("FE_P3CL::DHRef: Invalid shape function.");
    };
}

inline double
FE_P3CL::Laplace(Uint dof, const SMatrixCL<3,3>& M)
{
    double ret= 0.;
    for (Uint i=0; i<3; ++i)
      for(Uint j=0; j<3; ++j)
        for (Uint k=0; k<3; ++k)
        {
            ret+= M(i,j)*M(i,k)*D2HRef(dof, j, k);
        }
    return ret;
}


template<class T= double>
class LocalP3CL: public GridFunctionCL<T>
{
  public:
    typedef GridFunctionCL<T> base_type;
    typedef typename base_type::value_type value_type;
    typedef typename base_type::instat_fun_ptr instat_fun_ptr;
    typedef FE_P3CL FETYPE;//FE_P2CL

  protected:
    typedef LocalP3CL<T> self_;

  public:
    LocalP3CL() : base_type( value_type(), FE_P3CL::NumDoFC) {}
    LocalP3CL(const value_type& t): base_type( t, FE_P3CL::NumDoFC) {}
    // Initialize from a given function
    LocalP3CL(const TetraCL&, instat_fun_ptr , double= 0.0);
    // Initialize from VecDescCL and boundary-data
    template<class BndDataT>
      LocalP3CL(const TetraCL&, const VecDescCL&, const BndDataT&);
    // Initialize from PiEvalCL
    template <class P3FunT>
      LocalP3CL(const TetraCL&, const P3FunT&);
    // Initialize from LocalP1CL
    LocalP3CL(const LocalP1CL<T>&);

DROPS_DEFINE_VALARRAY_DERIVATIVE(LocalP3CL, T, base_type)

    // These "assignment-operators" correspond to the constructors
    // with multiple arguments
    inline self_&
    assign(const TetraCL&, instat_fun_ptr, double= 0.0);
    template<class BndDataT>
      inline self_&
      assign(const TetraCL&, const VecDescCL&, const BndDataT&);
    template<class BndDataT>
      inline self_&
      assign_on_tetra(const TetraCL&, const VecDescCL&, const BndDataT&);
    template <class P3FunT>
      inline self_&
      assign(const TetraCL&, const P3FunT&);
    inline self_&
    assign(const LocalP1CL<T>&);

    // pointwise evaluation in barycentric coordinates
    inline value_type operator()(const BaryCoordCL&) const;
};
//**************************************************************************
// RestrictP3: Stores the DoF-values of a P3-function corresponding to vd  *
//     and bnd for tetrahedron s in the container c.                       *
// Precondition: vd is a VecDescCL for a P3-function on level l, bnd is a  *
//     BndDataCL and s a tetrahedron on a level <= l. c is a container     *
//     (component access with []) that can hold at least 20 values of f's  *
//     return type.                                                        *
// Postcondition: c contains the value of f in the 20 DoF in the order used*
//     by FE_P3CL.                                                         *
//**************************************************************************
//template <class VecDescT, class BndDataT, class Cont>
//void RestrictP3(const TetraCL& s, const VecDescT& vd, const BndDataT& bnd, Cont& c);//RestrictP2
template <class VecDescT, class BndDataT, class Cont>
void RestrictP3(const TetraCL& s, const VecDescT& vd, const BndDataT& bnd, Cont& c)
{
    const Uint slvl= s.GetLevel();
    const Uint flvl= vd.GetLevel();
    Assert( slvl<=flvl, DROPSErrCL("RestrictP3: Tetra is on a finer level"
            "than the function."), ~0);

    typedef typename VecDescT::DataType VecT;
    typedef DoFHelperCL< typename BndDataT::bnd_type, VecT> DoFT;
    const VecT& v= vd.Data;
    const Uint idx= vd.RowIdx->GetIdx();
    for (Uint i= 0; i<NumVertsC; ++i)
        c[i]= !bnd.IsOnDirBnd( *s.GetVertex( i))
                ? DoFT::get( v, s.GetVertex( i)->Unknowns( idx))
                : bnd.GetDirBndValue( *s.GetVertex( i), vd.t);
    for (Uint i= 0; i<NumEdgesC; ++i) {
        const EdgeCL& e= *s.GetEdge( i);
        c[i+NumVertsC]= (slvl < flvl && e.IsRefined())
            ? ( !bnd.IsOnDirBnd( *e.GetMidVertex())
                ? DoFT::get( v, e.GetMidVertex()->Unknowns( idx))
                : bnd.GetDirBndValue( *e.GetMidVertex(), vd.t))
            : ( !bnd.IsOnDirBnd( e)
                ? DoFT::get( v, e.Unknowns( idx))
                : bnd.GetDirBndValue( e, vd.t));
    }
}

template <class P3FuncT, class Cont>
void RestrictP3(const TetraCL& s, const P3FuncT& f, Cont& c)
{
    RestrictP3( s, *f.GetSolution(), *f.GetBndData(), c);
}

//**************************************************************************
// Class:   LocalP3CL                                                      *
//**************************************************************************
template<class T>
  inline LocalP3CL<T>&
  LocalP3CL<T>::assign(const TetraCL& s, instat_fun_ptr f, double t)
{
    for (Uint i= 0; i< NumVertsC; ++i)
        (*this)[i]= f( s.GetVertex( i)->GetCoord(), t);
    for (Uint i= 0; i< NumEdgesC; ++i)
        (*this)[i+NumVertsC]= f( GetBaryCenter( *s.GetEdge( i)), t);
    return *this;
}

template<class T>
  template<class BndDataT>
    inline LocalP3CL<T>&
    LocalP3CL<T>::assign(const TetraCL& s,
        const VecDescCL& vd, const BndDataT& bnd)
{
    typedef VecDescCL::DataType VecT;
    typedef DoFHelperCL<value_type, VecT> DoFT;
    const VecT& v= vd.Data;
    if (vd.RowIdx->IsDG())
    { // This is just for P3
        Uint idx_num = vd.RowIdx->GetIdx();
        Uint first = s.Unknowns(idx_num);
        for (int i = 0; i < 20; ++i)
        {
            (*this)[i] = DoFT::get( vd.Data, first++);
        }
        return *this;
    }
    const Uint tlvl= s.GetLevel();
    const Uint vlvl= vd.GetLevel();
    const Uint idx= vd.RowIdx->GetIdx();
    if (tlvl == vlvl) {
        for (Uint i= 0; i< NumVertsC; ++i)
            (*this)[i]= !bnd.IsOnDirBnd( *s.GetVertex( i))
                ? DoFT::get( v, s.GetVertex( i)->Unknowns( idx))
                : bnd.GetDirBndValue( *s.GetVertex( i), vd.t);
        for (Uint i= 0; i< NumEdgesC; ++i)
            (*this)[i+NumVertsC]= !bnd.IsOnDirBnd( *s.GetEdge( i))
                ? DoFT::get( v, s.GetEdge( i)->Unknowns( idx))
                : bnd.GetDirBndValue( *s.GetEdge( i), vd.t);
    }
    else {
        if (tlvl < vlvl) RestrictP3( s, vd, bnd, *this);//RestrictP2
        else throw DROPSErrCL( "LocalP3CL::Assign: Prolongation not implemented.\n");
    }
    return *this;
}

template<class T>
  template<class BndDataT>
    inline LocalP3CL<T>&
    LocalP3CL<T>::assign_on_tetra(const TetraCL& s,
        const VecDescCL& vd, const BndDataT& bnd)
{
    typedef VecDescCL::DataType VecT;
    typedef DoFHelperCL<value_type, VecT> DoFT;
    const VecT& v= vd.Data;
    const Uint idx= vd.RowIdx->GetIdx();

    //const Uint tlvl= s.GetLevel();
    //const Uint flvl= vd.GetLevel();

    if (vd.RowIdx->IsDG())
    { // This is just for P3
        Uint first = s.Unknowns(idx);
        for (int i = 0; i < 20; ++i)
        {
            (*this)[i] = DoFT::get( v, first++);
        }
    }
    else
    {
        for (Uint i= 0; i< NumVertsC; ++i)
            (*this)[i]= !bnd.IsOnDirBnd( *s.GetVertex( i))
                ? DoFT::get( v, s.GetVertex( i)->Unknowns( idx))
                : bnd.GetDirBndValue( *s.GetVertex( i), vd.t);
        for (Uint i= 0; i< NumEdgesC; ++i)
            (*this)[i+NumVertsC]= !bnd.IsOnDirBnd( *s.GetEdge( i))
                ? DoFT::get( v, s.GetEdge( i)->Unknowns( idx))
                : bnd.GetDirBndValue( *s.GetEdge( i), vd.t);
    }
    return *this;
}

template<class T>
  template<class P3FunT>
    inline LocalP3CL<T>&
    LocalP3CL<T>::assign(const TetraCL& s, const P3FunT& f)
{
    typedef VecDescCL::DataType VecT;
    typedef DoFHelperCL<value_type, VecT> DoFT;
    const VecDescCL& vd = *(f.GetSolution());
    const Uint tlvl= s.GetLevel();
    const Uint flvl= vd.GetLevel();
    //const Uint idx= vd.RowIdx->GetIdx();

    if (vd.RowIdx->IsDG())
    { // This is just for P3
        if (tlvl != flvl)
            throw DROPSErrCL( "LocalP3CL::Assign: Prolongation not implemented.\n");
        Uint idx_num = vd.RowIdx->GetIdx();
        Uint first = s.Unknowns(idx_num);
        for (int i = 0; i < 20; ++i)
        {
            (*this)[i] = DoFT::get( vd.Data, first++);
        }
        return *this;
    }
    else

    {
        const Uint tlvl= s.GetLevel();
        const Uint flvl= f.GetLevel();
        if (tlvl == flvl)
            f.GetDoF( s, *this);
        else
            if (tlvl < flvl) RestrictP3( s, f, *this);
            else throw DROPSErrCL( "LocalP3CL::Assign: Prolongation not implemented.\n");
    }
    return *this;
}

template<class T>
  inline LocalP3CL<T>&
  LocalP3CL<T>::assign(const LocalP1CL<T>& p1)
{
    for (size_t i= 0; i < 4; ++i)
        (*this)[i]= p1[i];
    for (size_t i= 0; i < 6; ++i)
        (*this)[i + 4]= 0.5*(p1[VertOfEdge( i, 0)] + p1[VertOfEdge( i, 1)]);
    return *this;
}


template<class T>
  LocalP3CL<T>::LocalP3CL(const TetraCL& s, instat_fun_ptr f , double t)
  : base_type( value_type(), FE_P3CL::NumDoFC)
{
    this->assign( s, f, t);
}

template<class T>
  template <class P3FunT>
    LocalP3CL<T>::LocalP3CL(const TetraCL& s, const P3FunT& f)
    : base_type( value_type(), FE_P3CL::NumDoFC)
{
    this->assign( s, f);
}

template<class T>
  template<class BndDataT>
    LocalP3CL<T>::LocalP3CL(const TetraCL& s,
        const VecDescCL& vd, const BndDataT& bnd)
    : base_type( value_type(), FE_P3CL::NumDoFC)
{
    this->assign( s, vd, bnd);
}

template<class T>
  LocalP3CL<T>::LocalP3CL (const LocalP1CL<T>& p1)
    : base_type( value_type(), FE_P3CL::NumDoFC)
{
    this->assign( p1);
}

template<class T>
  inline typename LocalP3CL<T>::value_type
  LocalP3CL<T>::operator() (const BaryCoordCL& p) const
{
    return FE_P3CL::val( *this, p);
}

template<class T>
void ExtendP1onChild( const LocalP3CL<T>& isoP3, int child, LocalP3CL<T>& P1onParent)
{
    static const int childOrder[8]= { 0, 1, 7, 6, 2, 4, 5, 3};
    const Uint ch= childOrder[child];
    // children ordered such that
    // A) children 0,...,3 located at corners of parent, parent vertex ch is also a child vertex
    // B) children 4,...,7 located inside parent forming an octahedron, one face is on parent face ch-4
    // Note: data can be looked up in topo.h/cpp

    // first step: compute values on parent's vertices
    if (ch<4) // case A
    {
        const Uint pv= ch; // parent vertex which is also vertex of child

        P1onParent[pv]= isoP3[pv]; // copy value for vertex pv
        // 1D interpolation along the 3 edges starting at parent vertex pv
        for (Uint v=0; v<4; ++v)
            if (v != pv)
                P1onParent[v]= 2*isoP3[EdgeByVert(v,pv)+4] - isoP3[pv];
    }
    else // case B
    {
        const Uint pf= ch - 4; // parent face containing one of the child's faces
        // parent face contains 3 vertices of child, remaining fourth vertex is located on parent's edge pe (= edge 1, 4, 1, 4).
        const Uint ev= (pf + 2) % 4; // parent vertex in face pf which is located on parent edge pe (= vertex 2, 3, 0, 1)

        // 2D interpolation on parent face
        for (Uint i=0; i<3; ++i) {
            const Uint fv= VertOfFace(pf,i);
            P1onParent[fv]= 0.;
            for (Uint j=0; j<3; ++j) {
                const Uint fe= EdgeOfFace(pf,j);
                const T val= isoP3[fe+4];
                if (fv == VertOfEdge(fe,0))
                    P1onParent[fv]+= val;
                else if (fv == VertOfEdge(fe,1))
                    P1onParent[fv]+= val;
                else // vert fv is opposite to edge fe in face pf
                    P1onParent[fv]-= val;
             }
        }

        // 1D interpolation along edge pe whose midvertex is remaining 4th vertex of child
        P1onParent[OppVert(pf)]= 2*isoP3[EdgeByVert(OppVert(pf),ev)+4] - P1onParent[ev];
    }

    // second step: linear interpolation of edge values
    for (Uint e=0; e<6; ++e)
        P1onParent[e+4]= 0.5*(P1onParent[VertOfEdge(e,0)] + P1onParent[VertOfEdge(e,1)]);
}


class P3DiscCL
{
  public:
    // gradients on reference tetra
    static void GetGradientsOnRef( LocalP1CL<Point3DCL> GRef[20]);
    static void GetGradientsOnRef( Quad2CL<Point3DCL> GRef[20]);
    static void GetGradientsOnRef( Quad5CL<Point3DCL> GRef[20]);
    // The 2nd arg points to 3 vertices of the triangle
    static void GetGradientsOnRef( Quad5_2DCL<Point3DCL> GRef[20], const BaryCoordCL* const);
    // p3[i] contains a LocalP3CL-object that is initialized with FE_P3CL::Hi
    static void GetP3Basis( LocalP3CL<> p3[20]);
    // p3[i] contains a Quad5_2DCL-object that is initialized with FE_P3CL::Hi
    static void GetP3Basis( Quad5_2DCL<> p3[20], const BaryCoordCL* const p);
    // compute gradients
    static void GetGradients( LocalP1CL<Point3DCL> G[20], const LocalP1CL<Point3DCL> GRef[20], const SMatrixCL<3,3> &T)
    { for (int i=0; i<20; ++i) for (int j=0; j<4; ++j) G[i][j]= T*GRef[i][j]; }
    static void GetGradients( Quad2CL<Point3DCL> G[20], Quad2CL<Point3DCL> GRef[20], const SMatrixCL<3,3> &T)
    { for (int i=0; i<20; ++i) for (int j=0; j<5; ++j) G[i][j]= T*GRef[i][j]; }
    static void GetGradients( Quad5CL<Point3DCL> G[20], Quad5CL<Point3DCL> GRef[20], const SMatrixCL<3,3> &T)
    { for (int i=0; i<20; ++i) for (int j=0; j<Quad5DataCL::NumNodesC; ++j) G[i][j]= T*GRef[i][j]; }
    static void GetGradients( Quad5_2DCL<Point3DCL> G[20], Quad5_2DCL<Point3DCL> GRef[20], const SMatrixCL<3,3> &T)
    { for (int i=0; i<20; ++i) for (int j=0; j<Quad5_2DDataCL::NumNodesC; ++j) G[i][j]= T*GRef[i][j]; }
    static void GetGradient( Quad2CL<Point3DCL> &G, Quad2CL<Point3DCL> &GRef, const SMatrixCL<3,3> &T)
    { for (int j=0; j<5; ++j) G[j]= T*GRef[j]; }
    static void GetGradient( Quad5CL<Point3DCL> &G, Quad5CL<Point3DCL> &GRef, const SMatrixCL<3,3> &T)
    { for (int j=0; j<Quad5DataCL::NumNodesC; ++j) G[j]= T*GRef[j]; }
    /// compute gradient of a function provided as LocalP3CL<double> object
    template<class GradT>
    static void GetFuncGradient( GradT& gradF, const LocalP3CL<>& F, const GradT G[20])
    { gradF= F[0]*G[0]; for (int i=1; i<20; ++i) gradF+= F[i]*G[i]; }
    // Compute the Hessians H[d]= M*Href[d]*M^T
    static void GetHessians (SMatrixCL<3,3> H[20], const SMatrixCL<3,3>& M) {
        for (Uint d= 0; d< 20; ++d) {
            std::memset( &H[d], 0, 3*3*sizeof( double));
            for (Uint i= 0; i < 3; ++i)
                for (Uint j= 0; j < 3; ++j)
                    for (Uint k= 0; k < 3; ++k)
                        for (Uint l= 0; l < 3; ++l)
                            H[d](i,j)+= M(i, k)*M(j, l)*FE_P3CL::D2HRef( d, k, l);
        }
    }
    // cubatur formula for int f(x)*phi_i dx, exact up to degree 1
    static inline SVectorCL<3> Quad( const TetraCL& tetra, instat_vector_fun_ptr, Uint, double= 0.0);
    // cubatur formula for int f(x)*phi_i dx, exact up to degree 2
    template<class valT>
    static inline valT Quad( valT f[20], int i);
    // returns int phi_i phi_j dx
    static inline double GetMass( int i, int j);
    // returns int phi_i dx
    static inline double GetLumpedMass( int i) { return i<4 ? -1./120. : 1./30.; }
};


class LocalP3GradientCL
{
  private:
    const VecDescCL& f_;
    const BndDataCL<>& fbnd_;

    SMatrixCL<3,3> M;
    LocalP1CL<Point3DCL> GradRefLP1[20],
                         GradLP1[20],
                         p1grad;
    LocalP3CL<> p3;//LocalP2CL
    LocalP3CL<Point3DCL> p3grad;

  public:
    typedef Point3DCL value_type;
    static const int num_components= 3;

    LocalP3GradientCL (const VecDescCL& f, const BndDataCL<>& fbnd)
        : f_( f), fbnd_( fbnd)  { P3DiscCL::GetGradientsOnRef( GradRefLP1); }//P2DiscCL

    void set_tetra (const TetraCL* t);
    value_type&       operator[] (size_t i)       { return p3grad[i]; }
    const value_type& operator[] (size_t i) const { return p3grad[i]; }
    bool invalid_p (size_t /*i*/) const { return false; }
    void finalize_accumulation () const {}
};

void P3DiscCL::GetGradientsOnRef( LocalP1CL<Point3DCL> GRef[20])
{
    for (int i= 0; i < 20; ++i)
    {
        GRef[i][0]= FE_P3CL::DHRef( i, 0,0,0);//FE_P2CL
        GRef[i][1]= FE_P3CL::DHRef( i, 1,0,0);
        GRef[i][2]= FE_P3CL::DHRef( i, 0,1,0);
        GRef[i][3]= FE_P3CL::DHRef( i, 0,0,1);
    }
}

void P3DiscCL::GetGradientsOnRef( Quad2CL<Point3DCL> GRef[20])
{
    for (int i=0; i<20; ++i)
    {
        GRef[i][0]= FE_P3CL::DHRef( i, 0,0,0);
        GRef[i][1]= FE_P3CL::DHRef( i, 1,0,0);
        GRef[i][2]= FE_P3CL::DHRef( i, 0,1,0);
        GRef[i][3]= FE_P3CL::DHRef( i, 0,0,1);
        GRef[i][4]= FE_P3CL::DHRef( i, 0.25,0.25,0.25);
    }
}

void P3DiscCL::GetGradientsOnRef( Quad5CL<Point3DCL> GRef[20])
{
    for (int i=0; i<20; ++i)
        for (int j=0; j<Quad5DataCL::NumNodesC; ++j)
        {
            const BaryCoordCL& Node= Quad5DataCL::Node[j];
            GRef[i][j]= FE_P3CL::DHRef( i, Node[1], Node[2], Node[3]);
        }
}

void P3DiscCL::GetGradientsOnRef( Quad5_2DCL<Point3DCL> GRef[20],
    const BaryCoordCL* const p)
{
    BaryCoordCL NodeInTetra[Quad5_2DDataCL::NumNodesC];
    Quad5_2DCL<Point3DCL>::SetInterface( p, NodeInTetra);
    for (int i= 0; i < 20; ++i)
        for (int j= 0; j < Quad5_2DDataCL::NumNodesC; ++j) {
            const BaryCoordCL& Node= NodeInTetra[j];
            GRef[i][j]= FE_P3CL::DHRef( i, Node[1], Node[2], Node[3]);
        }
}

void P3DiscCL::GetP3Basis( LocalP3CL<> p3[20])
{
    for (int i= 0; i < 20; ++i) {
        for (int j= 0; j < 20; ++j)
            p3[i][j]= 0;
        p3[i][i]= 1;
    }
}

void P3DiscCL::GetP3Basis( Quad5_2DCL<> p3[20], const BaryCoordCL* const p)
{
    BaryCoordCL NodeInTetra[Quad5_2DDataCL::NumNodesC];
    Quad5_2DCL<>::SetInterface( p, NodeInTetra);
    for (int j= 0; j < Quad5_2DDataCL::NumNodesC; ++j) {
        const BaryCoordCL& Node= NodeInTetra[j];
        p3[0][j]= FE_P3CL::H0( Node);
        p3[1][j]= FE_P3CL::H1( Node);
        p3[2][j]= FE_P3CL::H2( Node);
        p3[3][j]= FE_P3CL::H3( Node);
        p3[4][j]= FE_P3CL::H4( Node);
        p3[5][j]= FE_P3CL::H5( Node);
        p3[6][j]= FE_P3CL::H6( Node);
        p3[7][j]= FE_P3CL::H7( Node);
        p3[8][j]= FE_P3CL::H8( Node);
        p3[9][j]= FE_P3CL::H9( Node);
    }
}

void LocalP3GradientCL::set_tetra (const TetraCL* t)
{
    double det; // dummy
    GetTrafoTr( M, det, *t);
    P3DiscCL::GetGradients( GradLP1, GradRefLP1, M);
    p3.assign( *t, f_, fbnd_);
    p1grad= Point3DCL();
    for (Uint i= 0; i < 20; ++i)
        p1grad+= p3[i]*GradLP1[i];
    for (Uint i= 0; i < 4; ++i)
        p3grad[i]= p1grad[i];
    for (Uint i= 0; i < 6; ++i)
        p3grad[i+4]= 0.5*(p1grad[VertOfEdge( i, 0)] + p1grad[VertOfEdge( i, 1)]);

}

/// \brief Collect indices of unknowns, boundary-segments and boundary
///     conditions on a tetrahedron.
///
/// This is convenient for discretisation of operators in the Setup-routines.
class LocalNumbP3CL
{
  public:
    /// \brief Field of unknown-indices; NoIdx, iff the degree of freedom lies
    /// on a boundary without unknowns. (Formerly called Numb.)
    IdxT     num   [20];
    /// \brief On boundaries, the number of the relevant BndSegDataCL-object
    /// in the corresponding BndDataCL-object, else NoBndC.
    BndIdxT  bndnum[20];
    /// \brief The relevant BndCondT, NoBC in the interior dofs.
    BndCondT bc    [20];

    /// \brief The default constructor leaves everything uninitialized.
    LocalNumbP3CL() {}
    /// \brief Read indices, boundary-segment numbers and boundary conditions
    /// from a tetrahedron and a BndDataCL-like object.
    template<class BndDataT>
      LocalNumbP3CL(const TetraCL&, const IdxDescCL&, const BndDataT&);

    /// \brief Read indices only
    /// from a tetrahedron.
    LocalNumbP3CL(const TetraCL&, const IdxDescCL&);

    /// \brief Read indices, boundary-segment numbers and boundary conditions
    ///     from a tetrahedron and a BndDataCL-like object.
    template<class BndDataT>
      void
      assign(const TetraCL& s, const IdxDescCL& idx, const BndDataT& bnd);

    /// \brief Compute the indices only.
    /// Only num is set up.
    void assign_indices_only (const TetraCL& s, const IdxDescCL& idx);

    /// \brief True, iff index i has a dof associated with it.
    bool WithUnknowns(IdxT i) const { return num[i] != NoIdx; }
};
void
LocalNumbP3CL::assign_indices_only (const TetraCL& s, const IdxDescCL& idx)
{
    const Uint sys= idx.GetIdx();
    if (!idx.IsDG())
    {
        for (Uint i= 0; i < 4; ++i)
            num[i]= s.GetVertex( i)->Unknowns.Exist( sys) ? s.GetVertex( i)->Unknowns( sys) : NoIdx;
        for(Uint i= 0; i < 6; ++i)
            num[i+4]= s.GetEdge( i)->Unknowns.Exist( sys) ? s.GetEdge( i)->Unknowns( sys)   : NoIdx;
    }
    else
    {
        Uint first = s.Unknowns(sys);
        for (int i = 0; i < 20; ++i)
            num[i] = first++;
    }
}

LocalNumbP3CL::LocalNumbP3CL(const TetraCL& s, const IdxDescCL& idx)
/// \param s The tet, from which index-numbers are read.
/// \param idx The IdxDescCL-object to be used.
{
    this->assign_indices_only( s, idx);
}

#if 1
template <typename LocalP3T>
class OswaldProjectionP3AccuCL : public TetraAccumulatorCL
{
  private:
    LocalP3T loc_;
    std::valarray<double>* n_,
                         * n_invalid_;
    bool check_averaging_;
    VecDescCL& avg_;

    LocalNumbP3CL numg;//LocalNumbP2CL

    const VecDescCL*   ls;      // a P3-level-set function
    const BndDataCL<>* lsetbnd; // boundary data for the level set function
    std::valarray<double> ls_loc;
    const PrincipalLatticeCL* lat;

    OswaldProjectionP3AccuCL& set_n (std::valarray<double>* n) { // The clones must refer to the n_ of thread 0.
        n_= n;
        return *this;
    }
    OswaldProjectionP3AccuCL& set_n_invalid (std::valarray<double>* n) { // The clones must refer to the n_invalid of thread 0.
        n_invalid_= n;
        return *this;
    }

  public:
    OswaldProjectionP3AccuCL (LocalP3T loc, VecDescCL& avg)
        : loc_( loc), n_( 0), n_invalid_( 0), check_averaging_( false), avg_( avg), ls( 0), lsetbnd( 0), lat( 0) {}

    OswaldProjectionP3AccuCL& set_check_averaging (bool b= true) {
        check_averaging_= b;
        return *this;
    }

    OswaldProjectionP3AccuCL& set_level_set_function (const VecDescCL* lsarg, const BndDataCL<>* lsetbndarg, const PrincipalLatticeCL* latarg) {
        ls= lsarg;
        lsetbnd= lsetbndarg;
        lat= latarg;
        ls_loc.resize ( lat ? lat->vertex_size () : 0);
        return *this;
    }

    virtual void begin_accumulation   () {
        n_= new std::valarray<double>( avg_.Data.size()/loc_.num_components);
        if (check_averaging_)
            n_invalid_= new std::valarray<double>( avg_.Data.size()/loc_.num_components);
    }
    virtual void finalize_accumulation() {
        loc_.finalize_accumulation ();
        if (check_averaging_)
            for (size_t i= 0; i < n_->size (); ++i)
                if (n_[0][i] == 0 && n_invalid_[0][i] > 0)
                    std::cerr << "OswaldProjectionP3AccuCL::finalize_accumulation: No local value for " << i << "; invalid_p: " << n_invalid_[0][i] << ".\n";
        delete n_;
        delete n_invalid_;
    }

    virtual void visit (const TetraCL& t) {
        if (ls != 0) {
            LocalP3CL<> locp3_ls( t, *ls, *lsetbnd);
            evaluate_on_vertexes( locp3_ls, *lat, Addr( ls_loc));
            if (equal_signs( ls_loc))
                return;
        }
        loc_.set_tetra( &t);
        numg.assign_indices_only( t, *avg_.RowIdx);
        for (Uint i= 0; i < 20; ++i) {
            if (!numg.WithUnknowns( i))
                continue;
            const IdxT dof= numg.num[i];
        if (loc_.invalid_p (i)) {
            if (check_averaging_)
                ++n_invalid_[0][dof/loc_.num_components];
            continue;
        }
            double& n= n_[0][dof/loc_.num_components]; // This assumes that the local gradient is in the components dof/3..dof/3 + 2.
            n+= 1.;
            typedef typename LocalP3T::value_type value_type;
            const value_type& oldavg= DoFHelperCL<value_type, VectorCL>::get( avg_.Data, dof);
            DoFHelperCL<value_type, VectorCL>::set( avg_.Data, dof, ((n - 1.)/n)*oldavg + (1./n)*loc_[i]);
        }
    }

    virtual TetraAccumulatorCL* clone (int /*clone_id*/) {
        OswaldProjectionP3AccuCL* p= new OswaldProjectionP3AccuCL( loc_, avg_);
        p->set_n( n_)
          .set_n_invalid (n_invalid_)
          .set_check_averaging (check_averaging_)
          .set_level_set_function (ls, lsetbnd, lat);
       return p;
    }
};
#endif

void averaging_P3_gradient_recovery (const MultiGridCL& mg, const VecDescCL& f, const BndDataCL<>& fbnd, VecDescCL& grad)
{
    LocalP3GradientCL loc( f, fbnd);
    OswaldProjectionP3AccuCL<LocalP3GradientCL> accu( loc, grad);//OswaldProjectionP2AccuCL
    TetraAccumulatorTupleCL accus;
    accus.push_back( &accu);
    accumulate( accus, mg, f.RowIdx->TriangLevel(), f.RowIdx->GetBndInfo());
};




