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
#pragma once
#include "misc/container.h"
//namespace DROPS
//{
double xyz_rhs (const DROPS::Point3DCL& p, double);
void lsFun(double x, double y, double z, double *value);
void lsGrad(double x, double y, double z, double *grad);
//}
