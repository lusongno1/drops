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
//namespace DROPS
//{
double xyz_rhs (const DROPS::Point3DCL& p, double)
{
    //my test case，f = 3*(x+y+z) for problem -\Delta u + u = f
    return 3*(p[0]+p[1]+p[2]);//p.norm();
}
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
//}
