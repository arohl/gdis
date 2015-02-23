/*
Copyright (C) 2003 by Sean David Fleming

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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "gdis.h"
#include "coords.h"
#include "edit.h"
#include "matrix.h"
#include "quaternion.h"
#include "model.h"
#include "morph.h"
#include "opengl.h"
#include "render.h"
#include "select.h"
#include "space.h"
#include "spatial.h"
#include "interface.h"

/* main structure */
extern struct sysenv_pak sysenv;

/**********************************************************/
/* adjust a rotation matrix wrt curent camera position */
/**********************************************************/
void matrix_camera_transform(gdouble *rot, struct model_pak *model)
{
gdouble q[9], qi[9], mat[9];

g_assert(model != NULL);

/* build component matrices */
quat_matrix(qi, camera_q(model->camera));
quat_matrix(q, camera_q(model->camera));
matrix_invert(qi);

/* build transformation */
memcpy(mat, model->latmat, 9*sizeof(gdouble));
matmat(qi, mat);
matmat(rot, mat);
matmat(q, mat);
matmat(model->ilatmat, mat);

memcpy(rot, mat, 9*sizeof(gdouble));
}

/**********************************************************/
/* construct a rotation matrix wrt curent camera position */
/**********************************************************/
void matrix_relative_rotation(gdouble *m, gdouble a, gint type, struct model_pak *model)
{
g_assert(model != NULL);

/* get rotation and transform according to camera */
matrix_rotation(m, a, type);
matrix_camera_transform(m, model);
}

/*******************************/
/* multiply two (3x3) matrices */
/*******************************/
void matmat(gdouble *mat1, gdouble *mat2)
/* mat2 -> mat1*mat2 */
{
gint i;
gdouble tmp[9], *res, *ptr1, *ptr2;

/* init */
matrix_transpose(mat2);
res = tmp;
ptr1 = mat1;

/* mult - loop over rows of mat 1*/
for (i=3 ; i-- ; )
  {
  ptr2 = mat2;

  *res = *(ptr1++) * (*(ptr2++));
  *res += *(ptr1++) * (*(ptr2++));
  *(res++) += *(ptr1--) * (*(ptr2++));
  ptr1--;

  *res = *(ptr1++) * (*(ptr2++));
  *res += *(ptr1++) * (*(ptr2++));
  *(res++) += *(ptr1--) * (*(ptr2++));
  ptr1--;

  *res = *(ptr1++) * (*(ptr2++));
  *res += *(ptr1++) * (*(ptr2++));
  *(res++) += *(ptr1++) * (*ptr2);
  }

/* copy result to mat2 */
memcpy(mat2,tmp,9*sizeof(gdouble));
}

/*******************************/
/* multiply two (4x4) matrices */
/*******************************/
void mat4mat(gdouble *mat1, gdouble *mat2)
/* mat2 -> mat1*mat2 */
{
gint i;
gdouble tmp[16], *res, *ptr1, *ptr2;

/* init */
matrix_transpose_44(mat2);
res = tmp;
ptr1 = mat1;

/* mult - loop over rows of mat 1*/
for (i=4 ; i-- ; )
  {
  ptr2 = mat2;

  *res = *(ptr1++) * (*(ptr2++));
  *res += *(ptr1++) * (*(ptr2++));
  *res += *(ptr1++) * (*(ptr2++));
  *(res++) += *(ptr1--) * (*(ptr2++));
  ptr1--;
  ptr1--;

  *res = *(ptr1++) * (*(ptr2++));
  *res += *(ptr1++) * (*(ptr2++));
  *res += *(ptr1++) * (*(ptr2++));
  *(res++) += *(ptr1--) * (*(ptr2++));
  ptr1--;
  ptr1--;

  *res = *(ptr1++) * (*(ptr2++));
  *res += *(ptr1++) * (*(ptr2++));
  *res += *(ptr1++) * (*(ptr2++));
  *(res++) += *(ptr1--) * (*(ptr2++));
  ptr1--;
  ptr1--;

  *res = *(ptr1++) * (*(ptr2++));
  *res += *(ptr1++) * (*(ptr2++));
  *res += *(ptr1++) * (*(ptr2++));
  *(res++) += *(ptr1++) * (*ptr2);
  }

/* copy result to mat2 */
memcpy(mat2,tmp,16*sizeof(gdouble));

return;
}

/****************************************/
/* transform a 3D vector by matrix mult */
/****************************************/
void vecmat(gdouble *mat, gdouble *vec)
{
gdouble x, y, z;
gdouble *mptr=mat, *vptr=vec;

/* init */
x = *vptr++;
y = *vptr++;
z = *vptr--;
/* mult */
*(--vptr) *= *mptr++;
*vptr += *mptr++ * y;
*vptr++ += *mptr++ * z;

*vptr = *mptr++ * x;
*vptr += *mptr++ * y;
*vptr++ += *mptr++ * z;

*vptr = *mptr++ * x;
*vptr += *mptr++ * y;
*vptr += *mptr * z;
}

/****************************************/
/* transform a 4D vector by matrix mult */
/****************************************/
void vec4mat(gdouble *mat, gdouble *vec)
{
gdouble x, y, z, t;
gdouble *mptr=mat, *vptr=vec;

/* init */
x = *vptr++;
y = *vptr++;
z = *vptr++;
t = *vptr--;
vptr--;
vptr--;

/* mult */
*vptr *= *mptr++;
*vptr += *mptr++ * y;
*vptr += *mptr++ * z;
*vptr++ += *mptr++ * t;

*vptr = *mptr++ * x;
*vptr += *mptr++ * y;
*vptr += *mptr++ * z;
*vptr++ += *mptr++ * t;

*vptr = *mptr++ * x;
*vptr += *mptr++ * y;
*vptr += *mptr++ * z;
*vptr++ += *mptr++ * t;

*vptr = *mptr++ * x;
*vptr += *mptr++ * y;
*vptr += *mptr++ * z;
*vptr += *mptr * t;
}

/************************************************************/
/* transform a 3D vector by matrix mult (transposed matrix) */
/************************************************************/
void vectmat(gdouble *mat, gdouble *vec)
{
gdouble x, y, z;
gdouble *mptr=mat, *vptr=vec;

/* init */
x = *vptr++;
y = *vptr++;
z = *vptr--;

/* mult */
*(--vptr) *= *mptr++;
*(++vptr) *= *mptr++;
*(++vptr) *= *mptr++;

vptr--;

*(--vptr) += *mptr++ * y;
*(++vptr) += *mptr++ * y;
*(++vptr) += *mptr++ * y;

vptr--;

*(--vptr) += *mptr++ * z;
*(++vptr) += *mptr++ * z;
*(++vptr) += *mptr * z;
}

/*****************/
/* transpose 4x4 */
/*****************/
void matrix_transpose_44(gdouble *mat)
{
gdouble tmp[16], *ptr=mat;

memcpy(tmp, mat, 16*sizeof(gdouble));
ptr++; 
*ptr++ = tmp[4]; 
*ptr++ = tmp[8]; 
*ptr++ = tmp[12]; 
*ptr++ = tmp[1]; 
ptr++;
*ptr++ = tmp[9]; 
*ptr++ = tmp[13]; 
*ptr++ = tmp[2]; 
*ptr++ = tmp[6]; 
ptr++;
*ptr++ = tmp[14]; 
*ptr++ = tmp[3]; 
*ptr++ = tmp[7]; 
*ptr   = tmp[11]; 
}

/*****************/
/* transpose 3x3 */
/*****************/
void matrix_transpose(gdouble *mat)
{
gdouble tmp[9], *ptr=mat;

memcpy(tmp, mat, 9*sizeof(gdouble));
ptr++; 
*ptr++ = tmp[3]; 
*ptr++ = tmp[6]; 
*ptr++ = tmp[1]; 
ptr++;
*ptr++ = tmp[7]; 
*ptr++ = tmp[2]; 
*ptr   = tmp[5]; 
}

void make_cofmat(gdouble *cof, gdouble *mat)
{
cof[0] = mat[4]*mat[8] - mat[5]*mat[7];
cof[1] = mat[5]*mat[6] - mat[3]*mat[8];
cof[2] = mat[3]*mat[7] - mat[4]*mat[6];
cof[3] = mat[2]*mat[7] - mat[1]*mat[8];
cof[4] = mat[0]*mat[8] - mat[2]*mat[6];
cof[5] = mat[1]*mat[6] - mat[0]*mat[7];
cof[6] = mat[1]*mat[5] - mat[2]*mat[4];
cof[7] = mat[2]*mat[3] - mat[0]*mat[5];
cof[8] = mat[0]*mat[4] - mat[1]*mat[3];
}

/********************************/
/* new matrix inversion routine */
/********************************/
#define DEBUG_INVMAT 0
gint matrix_invert(gdouble *mat)
{
gdouble d, id;
gdouble cof[9];

/* determinant */
d = matrix_determinant(mat);

#if DEBUG_INVMAT
P3MAT("input:", mat);
printf("determinant = %f\n", d);
#endif

/* FIXME - what is a good value to use here??? */
if (fabs(d) < 0.0001)
  {
  gui_text_show(ERROR, "Bad lattice matrix.");
  return(1);
  }
id = 1.0/d;

/* compute co-factor matrix */
make_cofmat(cof, mat);

matrix_transpose(cof);

VEC3MUL(&cof[0], id);
VEC3MUL(&cof[3], id);
VEC3MUL(&cof[6], id);

#if DEBUG_INVMAT
P3MAT("inverse matrix: ", cof);
#endif

matmat(cof, mat);

/* TODO - checks */
#if DEBUG_INVMAT
P3MAT("Should be I: ", mat);
#endif

memcpy(mat,cof,9*sizeof(gdouble));

return(0);
}

/****************************/
/* general vector magnitude */
/****************************/
gdouble magnitude(gdouble *vec, gint dim)
{
gint i;
gdouble sum=0.0;

for (i=0 ; i<dim ; i++)
  sum += (*(vec+i)) * (*(vec+i));

return(sqrt(sum));
}

/********************************/
/* general vector normalization */
/********************************/
gint normalize(gdouble *vec, gint dim)
{
gint i;
gdouble len, *ptr=vec;

len = magnitude(vec, dim);
if (len < FRACTION_TOLERANCE)
  return(1);

for (i=0 ; i<dim ; i++)
  *ptr++ /= len;

return(0);
}

/*****************************/
/* vector intersection angle */
/*****************************/
gdouble via(gdouble *vec1, gdouble *vec2, gint dim)
{
gint i;
gdouble lenprod, dot=0.0;
gdouble *ptr1=vec1, *ptr2=vec2;
gdouble cosa=0.0;

for (i=0 ; i<dim ; i++)
  dot += (*ptr1++) * (*ptr2++);

lenprod = magnitude(vec1, dim) * magnitude(vec2, dim);

/* get cos of the angle */
if (lenprod > FRACTION_TOLERANCE)
  cosa = dot / lenprod;
 
/* enforce range */
cosa = CLAMP(cosa, -1.0, 1.0);

return(acos(cosa));
}

/*************************/
/* cross product for 3x3 */
/*************************/
void crossprod(gdouble *res, gdouble *vec1, gdouble *vec2)
{
*(res+0) = *(vec1+1) * *(vec2+2) - *(vec1+2) * *(vec2+1);
*(res+1) = *(vec1+2) * *(vec2)   - *(vec1)   * *(vec2+2);
*(res+2) = *(vec1)   * *(vec2+1) - *(vec1+1) * *(vec2);
}

/***************************************************/
/* calculate normal of a plane defined by 3 points */
/***************************************************/
void calc_norm(gdouble *res, gdouble *a, gdouble *b, gdouble *c)
{
gdouble vec1[3], vec2[3];

/* a-b */
ARR3SET(vec1, b);
ARR3SUB(vec1, a);
/* a-c */
ARR3SET(vec2, c);
ARR3SUB(vec2, a);
/* normal */
crossprod(res, vec1, vec2);
}

/*********************************************/
/* calculate the distance between two points */
/*********************************************/
gdouble calc_sep(gdouble *a, gdouble *b)
{
gdouble x[3];

ARR3SET(x, a);
ARR3SUB(x, b);

return(VEC3MAG(x));
}

/*****************************/
/* project vector onto plane */
/*****************************/
void proj_vop(gdouble *res, gdouble *vec, gdouble *plane_norm)
{
gdouble len;
gdouble tmp[3];

/* project vector onto normal */
/* NB: assumes plane_norm is normalized */
ARR3SET(tmp, plane_norm);
ARR3MUL(tmp, vec);
len = tmp[0] + tmp[1] + tmp[2];

/* get vector to plane */
ARR3SET(tmp, plane_norm);
VEC3MUL(tmp, len);

/* subtract to get projected vector */
ARR3SET(res, vec);
ARR3SUB(res, tmp);
}

/***********************/
/* determinant for 3x3 */
/***********************/
gdouble matrix_determinant(gdouble *mat)
{
gdouble vec1[3], vec2[3], vec3[3], tmp[3];

/* make vectors from columns */
VEC3SET(vec1, mat[0], mat[3], mat[6]);
VEC3SET(vec2, mat[1], mat[4], mat[7]);
VEC3SET(vec3, mat[2], mat[5], mat[8]);
/* get cross product of 2 & 3 */
crossprod(tmp, vec2, vec3);
/* return dot prod with 1 */
ARR3MUL(vec1, tmp);
return(vec1[0]+vec1[1]+vec1[2]);
}

/**********************************************/
/* compute volume, assuming 3D lattice matrix */
/**********************************************/
gdouble calc_volume(gdouble *latmat)
{
gdouble vec[3], tmat[9];

/* NB: latmat has cell vectors in columns */
memcpy(tmat, latmat, 9*sizeof(gdouble));
matrix_transpose(tmat);

/* compute volume */
crossprod(vec, &tmat[0], &tmat[3]); 
ARR3MUL(vec, &tmat[6]);

return(fabs(vec[0]+vec[1]+vec[2]));
}

/****************************************************/
/* compute surface area, assuming 2D lattice matrix */
/****************************************************/
gdouble calc_area(gdouble *latmat)
{
gdouble a, b, c, tmat[9];

/* NB: latmat has cell vectors in columns */
memcpy(tmat, latmat, 9*sizeof(gdouble));
matrix_transpose(tmat);

/* get lengths and angle between */
a = VEC3MAG(&tmat[0]);
b = VEC3MAG(&tmat[3]);
c = via(&tmat[0], &tmat[3], 3);

/* return volume */
return(a * b * sin(c));
}

/*****************************/
/* initialize a 3x3 identity */
/*****************************/
void matrix_identity(gdouble *mat)
{
VEC3SET(&mat[0], 1.0, 0.0, 0.0);
VEC3SET(&mat[3], 0.0, 1.0, 0.0);
VEC3SET(&mat[6], 0.0, 0.0, 1.0);
}

/************************************************************/
/* assign a transformation matrix that aligns a vector to z */
/************************************************************/
#define DEBUG_Z_ALIGNMENT 0
void matrix_z_alignment(gdouble *mat, gdouble *vec)
{
gdouble len, proj, iproj, n[3], mop[9];

/* default return matrix */
matrix_identity(mat);
ARR3SET(n, vec);

/* check */
len = VEC3MAG(n);
if (len < FRACTION_TOLERANCE)
  {
  printf("matrix_z_alignment(): bad orientation vector.\n");
  return;
  }

/* normalize */
VEC3MUL(n, 1.0/len);
proj = sqrt(n[0]*n[0] + n[1]*n[1]);
iproj = 1.0 / proj;

/* rot z -> yz plane */
if (proj > FRACTION_TOLERANCE)
  {
  VEC3SET(&mop[0], n[1]*iproj, -n[0]*iproj, 0.0);
  VEC3SET(&mop[3], n[0]*iproj,  n[1]*iproj, 0.0);
  VEC3SET(&mop[6], 0.0, 0.0, 1.0);
  matmat(mop, mat);
  }

/* rot x -> colinear with z */
VEC3SET(&mop[0], 1.0, 0.0, 0.0);
VEC3SET(&mop[3], 0.0, n[2], -proj);
VEC3SET(&mop[6], 0.0, proj, n[2]);
matmat(mop, mat);

#if DEBUG_Z_ALIGNMENT
ARR3SET(n, vec);
vecmat(mat, n);
P3MAT("mat: ", mat);
P3VEC("v0: ", vec);
P3VEC("v1: ", n);
#endif
}

/********************************************************/
/* arbitrary alignment, maps v1 to be co-linear with v2 */
/********************************************************/
#define DEBUG_V_ALIGNMENT 0
void matrix_v_alignment(gdouble *m1, gdouble *v1, gdouble *v2)
{
gdouble m2[9];

/* map v1 to z, then to v2 via inverse of v2's z alignment */
matrix_z_alignment(m1, v1);
matrix_z_alignment(m2, v2);
matrix_invert(m2);
matmat(m2, m1);
}

/****************************************************/
/* construct rotation matrix about arbitrary vector */
/****************************************************/
void matrix_v_rotation(gdouble *mat, gdouble *v, gdouble angle)
{
gdouble m1[9], m2[9];

/* align with z */
matrix_z_alignment(m1, v);
memcpy(mat, m1, 9*sizeof(gdouble));

/* rotation operation (about z) */
matrix_rotation(m2, angle, YAW);
matmat(m2, mat);

/* inverse z alignment */
matrix_invert(m1);
matmat(m1, mat);
}

/****************************************************/
/* construct rotation matrix about arbitrary vector */
/****************************************************/
void matrix_v_reflection(gdouble *mat, gdouble *v)
{
gdouble m1[9], m2[9];

/* align with z */
matrix_z_alignment(m1, v);
memcpy(mat, m1, 9*sizeof(gdouble));

/* reflection about z */
matrix_identity(m2);
m2[8] = -1.0;
matmat(m2, mat);

/* inverse z alignment */
matrix_invert(m1);
matmat(m1, mat);
}

/*******************************/
/* construct a rotation matrix */
/*******************************/
void matrix_rotation(gdouble *rot, gdouble da, gint type)
{
gdouble cosa, sina;

/* precalc */
cosa = cos(da);
sina = sin(da);

/* dependant matrix elements */
switch(type)
  {
  case PITCH:
    VEC3SET(&rot[0], 1.0,  0.0, 0.0);
    VEC3SET(&rot[3], 0.0, cosa,-sina);
    VEC3SET(&rot[6], 0.0, sina, cosa);
    break;

  case YAW:
    VEC3SET(&rot[0], cosa, sina, 0.0);
    VEC3SET(&rot[3], -sina, cosa, 0.0);
    VEC3SET(&rot[6],  0.0,  0.0, 1.0);
    break;

  case ROLL:
    VEC3SET(&rot[0], cosa, 0.0,-sina);
    VEC3SET(&rot[3],  0.0, 1.0,  0.0);
    VEC3SET(&rot[6], sina, 0.0, cosa);
    break;

  case INITIAL:
  default:
    VEC3SET(&rot[0], 1.0, 0.0, 0.0);
    VEC3SET(&rot[3], 0.0, 1.0, 0.0);
    VEC3SET(&rot[6], 0.0, 0.0, 1.0);
    break;
  }
}

/**************************************************************/
/* get reciprocal lattice vectors from direct lattice vectors */
/**************************************************************/
#define DEBUG_RLAT 0
void matrix_reciprocal_init(struct model_pak *data)
{
gdouble norm;
gdouble a[3], b[3], c[3];
gdouble axb[3], bxc[3], cxa[3];

#if DEBUG_RLAT
P3MAT("direct lattice matrix:",data->latmat);
#endif

/* extract direct lattice vectors */
a[0] = data->latmat[0];
a[1] = data->latmat[3];
a[2] = data->latmat[6];
b[0] = data->latmat[1];
b[1] = data->latmat[4];
b[2] = data->latmat[7];
c[0] = data->latmat[2];
c[1] = data->latmat[5];
c[2] = data->latmat[8];

/* calc intermediate products */
crossprod(bxc, b, c);
crossprod(cxa, c, a);
crossprod(axb, a, b);
norm = 1.0/(a[0]*bxc[0] + a[1]*bxc[1] + a[2]*bxc[2]);

/* create reciprocal lattice vectors */
ARR3SET(a, bxc);
ARR3SET(b, cxa);
ARR3SET(c, axb);
VEC3MUL(a, norm);
VEC3MUL(b, norm);
VEC3MUL(c, norm);

/* store */
data->rlatmat[0] = a[0];
data->rlatmat[3] = a[1];
data->rlatmat[6] = a[2];
data->rlatmat[1] = b[0];
data->rlatmat[4] = b[1];
data->rlatmat[7] = b[2];
data->rlatmat[2] = c[0];
data->rlatmat[5] = c[1];
data->rlatmat[8] = c[2];

#if DEBUG_RLAT
P3MAT("reciprocal lattice matrix:",data->rlatmat);
#endif
}

/******************************************/
/* construct the unit translation vectors */
/******************************************/
#define DEBUG_XLAT 0
void matrix_lattice_init(struct model_pak *data)
{
gdouble n[3], v1[3], v2[3], v3[3];
gdouble b1, b2, b3, c1, c2, c3;

/* use a supplied latmat, rather than generating from pbc's */
/* NB: should be in gdis style colum vector format (gulp/marvin are in rows) */
if (data->construct_pbc)
  {
#if DEBUG_XLAT
printf("constructing pbc...\n");
#endif
/* get lattice vector lengths */
  VEC3SET(v1, data->latmat[0], data->latmat[3], data->latmat[6]);
  VEC3SET(v2, data->latmat[1], data->latmat[4], data->latmat[7]);
  VEC3SET(v3, data->latmat[2], data->latmat[5], data->latmat[8]);
/* set lengths */
  data->pbc[0] = VEC3MAG(v1);
  data->pbc[1] = VEC3MAG(v2);
  data->pbc[2] = VEC3MAG(v3);
/* get cell angles */
  data->pbc[3] = via(v2,v3,3);
  data->pbc[4] = via(v1,v3,3);
  data->pbc[5] = via(v1,v2,3);

/* NEW - handle cell angle signs */
  crossprod(n, v2 ,v3);
  ARR3MUL(n, v3);
  if (n[0]+n[1]+n[2] < 0.0)
    data->pbc[3] *= -1.0;

  crossprod(n, v1 ,v3);
  ARR3MUL(n, v3);
  if (n[0]+n[1]+n[2] < 0.0)
    data->pbc[4] *= -1.0;

  crossprod(n, v1 ,v2);
  ARR3MUL(n, v3);
  if (n[0]+n[1]+n[2] < 0.0)
    data->pbc[5] *= -1.0;
  }
else
  {
#if DEBUG_XLAT
printf("constructing latmat...\n");
#endif
/* construct lattice matrix from the unit cell lengths & angles */
/* this basically works by using the cosine rule in conjunction */
/* with a few geometric constraints (eg normalized vectors) */

/* FIXME - correction needed for 1D case as well? */
if (data->periodic == 2)
  {
/* lattice code requires non existent c parameter to be 1 */
  data->pbc[2] = 1.0;
  }

/* compute the translation vector for b */
  b1 = cos(data->pbc[5]);
  b2 = sin(data->pbc[5]);
  b3 = 0.0;              /* constrain a,b to remain on the x,y plane */

  if (b2 == 0.0)
    {
    printf("matrix_lattice_init(): bad cell parameters.\n");
    b2 = 1.0;
    }

/* compute the translation vector for c */
  c1 = cos(data->pbc[4]);
  c2 = (2.0*cos(data->pbc[3]) + b1*b1 + b2*b2 - 2.0*b1*c1 - 1.0)/(2.0*b2);
  c3 = sqrt(1.0 - c1*c1 - c2*c2);

/* assign in rows to make it easier to scale */
/* x & a are assumed co-linear */
  VEC3SET(&data->latmat[0], data->pbc[0], 0.0, 0.0);
  VEC3SET(&data->latmat[3], b1, b2, b3);
  VEC3SET(&data->latmat[6], c1, c2, c3);

/* scale b & c vectors up */
  VEC3MUL(&data->latmat[3], data->pbc[1]);
  VEC3MUL(&data->latmat[6], data->pbc[2]);

/* get vectors in cols */
  matrix_transpose(data->latmat);
  }

/* update dependents */
make_axes(data);
make_cell(data);
memcpy(data->ilatmat, data->latmat, 9*sizeof(gdouble));
matrix_invert(data->ilatmat);
matrix_reciprocal_init(data);
switch(data->periodic)
  {
  case 3:
    data->volume = calc_volume(data->latmat);
    break;
  case 2:
    data->area = calc_area(data->latmat);
    break;
  }

/* NEW - store angles in degrees as well as radians */
ARR3SET(data->cell_angles, &data->pbc[3]);
VEC3MUL(data->cell_angles, R2D);

#if DEBUG_XLAT
printf("cell dimensions: [%6.2f %6.2f %6.2f] (%6.2f %6.2f %6.2f)\n",
       data->pbc[0], data->pbc[1], data->pbc[2],
       data->cell_angles[0], data->cell_angles[1], data->cell_angles[2]);
P3MAT("lattice matrix:",data->latmat);
#endif
return;
}

/****************************/
/* build new lattice matrix */
/****************************/
/* v - new orientation (may need projection if model is 2D) */
/* ix - which lattice vector to replace (0-2) */
#define DEBUG_NEW_LATTICE 0
void matrix_lattice_new(gdouble *latmat, struct model_pak *model)
{
gint i;
gint ia, ib, ic, ma, mb, mc;
gdouble tmat[9], mat1[9], mat2[9];
gdouble vec[3];
GSList *list;
struct core_pak *core;
struct shel_pak *shel;
struct model_pak *dest;

/* checks */
g_assert(latmat != NULL);
g_assert(model != NULL);

#if DEBUG_NEW_LATTICE
P3MAT("transformation matrix: :", latmat);
P3MAT("old latmat:", model->latmat);
#endif

/* get the transformation matrix */
memcpy(tmat, latmat, 9*sizeof(gdouble));

/* transpose the source lattice matrix */
memcpy(mat1, model->latmat, 9*sizeof(gdouble));
matrix_transpose(mat1);

/* get the new lattice matrix */
matmat(tmat, mat1);
matrix_transpose(mat1);

#if DEBUG_NEW_LATTICE
P3MAT("new latmat:", mat1);
#endif

memcpy(mat2, mat1, 9*sizeof(gdouble));
if (matrix_invert(mat2))
  {
  gui_text_show(ERROR, "Two or more lattice vectors are not independant.");
  return;
  }

/* transformation -> new model */
dest = model_new();
dest->periodic = model->periodic;
memcpy(dest->latmat, mat1, 9*sizeof(gdouble));
gulp_data_copy(model, dest);
gulp_extra_copy(model, dest);

/* setup repeats required to fill the new cell */
ma=mb=mc=1;

/* a direction repeat */
for (i=0 ; i<3 ; i++)
  {
/* NB: add tolerance, since precision can be a problem */
  ia = (gint) (FRACTION_TOLERANCE + fabs(tmat[i]));
  if (ia > ma)
    ma = ia;
  }
/* b direction repeat */
for (i=3 ; i<6 ; i++)
  {
/* NB: add tolerance, since precision can be a problem */
  ib = (gint) (FRACTION_TOLERANCE + fabs(tmat[i]));
  if (ib > mb)
    mb = ib;
  }
/* c direction repeat */
for (i=6 ; i<9 ; i++)
  {
/* NB: add tolerance, since precision can be a problem */
  ic = (gint) (FRACTION_TOLERANCE + fabs(tmat[i]));
  if (ic > mc)
    mc = ic;
  }

/* TODO - check if these values are very large - warn user */
/* TDOO - implement a correct y/n/continue? popup */
#if DEBUG_NEW_LATTICE
printf("ma=%d, mb=%d, mc=%d\n", ma, mb, mc);
#endif

/* modify bulk energy */
dest->gulp.sbulkenergy *= ma*mb*mc;

/* full loop required to cover one cell in the transformed coordinates */
for (ic=0 ; ic<mc ; ic++)
  {
  for (ib=0 ; ib<mb ; ib++)
    {
    for (ia=0 ; ia<ma ; ia++)
      {
/* current source cell translation */
      VEC3SET(vec, ia, ib, ic);

/* create CARTESIAN core/shell images */
      for (list = model->cores ; list ; list=g_slist_next(list))
        {
        core = dup_core(list->data);
        core->primary = TRUE;
        core->primary_core = NULL;
        core->orig = TRUE;
        ARR3ADD(core->x, vec);
        vecmat(model->latmat, core->x);
        dest->cores = g_slist_prepend(dest->cores, core);
        if (core->shell)
          {
          shel = core->shell;
          shel->primary = TRUE;
          shel->primary_shell = NULL;
          shel->orig = TRUE;
          ARR3ADD(shel->x, vec);
          vecmat(model->latmat, shel->x);
          dest->shels = g_slist_prepend(dest->shels, shel);
          }
        }

      }
    }
  }

/* init new model */
dest->fractional = FALSE;
dest->construct_pbc = TRUE;
model_prep(dest);
tree_model_add(dest);
redraw_canvas(SINGLE);
}

/************************/
/* distinguish identity */
/************************/
gint matrix_is_identity(gdouble *mat)
{
if (*(mat+0) != 1.0)
  return(FALSE);
if (*(mat+4) != 1.0)
  return(FALSE);
if (*(mat+8) != 1.0)
  return(FALSE);

if (*(mat+1) != 0.0)
  return(FALSE);
if (*(mat+2) != 0.0)
  return(FALSE);
if (*(mat+3) != 0.0)
  return(FALSE);

if (*(mat+5) != 0.0)
  return(FALSE);
if (*(mat+6) != 0.0)
  return(FALSE);
if (*(mat+7) != 0.0)
  return(FALSE);

return(TRUE);
}

/*************************/
/* distinguish inversion */
/*************************/
gint matrix_is_inversion(gdouble *mat)
{
if (*(mat+0) != -1.0)
  return(FALSE);
if (*(mat+4) != -1.0)
  return(FALSE);
if (*(mat+8) != -1.0)
  return(FALSE);

if (*(mat+1) != 0.0)
  return(FALSE);
if (*(mat+2) != 0.0)
  return(FALSE);
if (*(mat+3) != 0.0)
  return(FALSE);

if (*(mat+5) != 0.0)
  return(FALSE);
if (*(mat+6) != 0.0)
  return(FALSE);
if (*(mat+7) != 0.0)
  return(FALSE);

return(TRUE);
}

/**************************************/
/* get the order of a symmetry matrix */
/**************************************/
gint matrix_order(gdouble *mat)
{
gint i;
gdouble m1[9], m2[9];

ARR3SET(&m1[0], (mat+0));
ARR3SET(&m1[3], (mat+3));
ARR3SET(&m1[6], (mat+6));
ARR3SET(&m2[0], (mat+0));
ARR3SET(&m2[3], (mat+3));
ARR3SET(&m2[6], (mat+6));

/* keep applying until get the identity */
i=1;
while(i<17)
  {
  if (matrix_is_identity(m2))
    return(i);
  matmat(m1, m2);
  i++;
  }

/* not a symmetry operator */
return(0);
}

/************************************************/
/* return the order if a matrix is a z rotation */
/************************************************/
/* FIXME - only does 2 fold (for defect setup) */
gint matrix_is_z_rotation(gdouble *mat)
{
if (*(mat+0) != -1.0)
  return(0);
if (*(mat+1) != 0.0)
  return(0);
if (*(mat+2) != 0.0)
  return(0);

if (*(mat+3) != 0.0)
  return(0);
if (*(mat+4) != -1.0)
  return(0);
if (*(mat+5) != 0.0)
  return(0);

if (*(mat+6) != 0.0)
  return(0);
if (*(mat+7) != 0.0)
  return(0);
if (*(mat+8) != 1.0)
  return(0);

return(2);
}
