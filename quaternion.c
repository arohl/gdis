/*
Copyright (C) 2003 by Craig Andrew James Fisher
Copyright (C) 2000 by Sean David Fleming

sean@ivec.org

This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation; either version 2
of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.

The GNU GPL can also be found at http://www.gnu.org
*/

#include <math.h>
#include <stdio.h>
#include "gdis.h"
#include "matrix.h"
#include "interface.h"

/*****************************/
/* quaternion multiplication */
/*****************************/
void quat_mult(gdouble *q1, gdouble *q2)
{
gdouble q[4];

q[0] = q1[0]*q2[0] - q1[1]*q2[1] - q1[2]*q2[2] - q1[3]*q2[3];
q[1] = q1[0]*q2[1] + q1[1]*q2[0] + q1[2]*q2[3] - q1[3]*q2[2];
q[2] = q1[0]*q2[2] - q1[1]*q2[3] + q1[2]*q2[0] + q1[3]*q2[1];
q[3] = q1[0]*q2[3] + q1[1]*q2[2] - q1[2]*q2[1] + q1[3]*q2[0];

/* place the result in q1 */
ARR4SET(q1, q);
}

/***********************************************/
/* convert a (rotation) matrix to a quaternion */
/***********************************************/
void quat_convert_matrix(gdouble *mat, gdouble *q)
{
gdouble s, trace;

trace = 1.0 + mat[0] + mat[4] + mat[8];

if (trace < FRACTION_TOLERANCE)
  printf("WARNING: ");

s = 2.0 * sqrt(trace);
q[0] = (mat[7] - mat[5]) / s;
q[1] = (mat[2] - mat[6]) / s;
q[2] = (mat[3] - mat[1]) / s;
q[3] = 0.25 * s;
}

/***********************************************************/
/* concatenate a quaternion rotation (v = axis, a = angle) */
/***********************************************************/
void quat_concat(gdouble *q, gdouble *v, gdouble a)
{
gdouble ha, sa, qr[4];

/* create and apply the desired quaternion rotation */
ha = 0.5*a;
sa = sin(ha);
VEC4SET(qr, cos(ha), v[0]*sa, v[1]*sa, v[2]*sa);
quat_mult(q, qr);
}

/*********************************/
/* standard quaternion rotations */
/*********************************/
void quat_concat_euler(gdouble *q, gint type, gdouble a)
{
gdouble v[3];

switch (type)
  {
  case PITCH:
    VEC3SET(v, 1.0, 0.0, 0.0);
    break;
  case ROLL:
    VEC3SET(v, 0.0, 1.0, 0.0);
    break;
  default:
  case YAW:
    VEC3SET(v, 0.0, 0.0, 1.0);
    break;
  }

/* general axis rotation */
quat_concat(q, v, a);
}

/*********************************************/
/* convert quaternion into a rotation matrix */
/*********************************************/
void quat_matrix(gdouble *mat, gdouble *quat)
{
gdouble q0, q1, q2, q3;
gdouble a01, a02, a03, a12, a13, a23;

q0 = quat[0];  q1 = quat[1]; q2 = quat[2];  q3 = quat[3];

a01 = 2.0*q0*q1;  a02 = 2.0*q0*q2;  a03 = 2.0*q0*q3;
                  a12 = 2.0*q1*q2;  a13 = 2.0*q1*q3;
                                    a23 = 2.0*q2*q3;
mat[1] = a12 - a03;
mat[2] = a13 + a02;
mat[3] = a12 + a03;
mat[5] = a23 - a01;
mat[6] = a13 - a02;
mat[7] = a23 + a01;

q0 = q0*q0;  q1 = q1*q1;  q2 = q2*q2;  q3 = q3*q3;

mat[0] = q0 + q1 - q2 - q3;
mat[4] = q0 - q1 + q2 - q3;
mat[8] = q0 - q1 - q2 + q3;
}

/*****************************************************************/
/* transform vector using a quaternion generated rotation matrix */
/*****************************************************************/
void quat_rotate(gdouble *vec, gdouble *quat)
{
gdouble rot[9];
 
quat_matrix(rot, quat);
vecmat(rot, vec);
}

