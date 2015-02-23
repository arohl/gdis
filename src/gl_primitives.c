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
#include <math.h>

#include "gdis.h"
#include "coords.h"
#include "matrix.h"
#include "numeric.h"
#include "quaternion.h"
#include "opengl.h"
#include "render.h"
#include "interface.h"
#ifdef __APPLE__
#include <OpenGL/gl.h>
#else
#include <GL/gl.h>
#endif

/* externals */
extern struct sysenv_pak sysenv;

/***************************/
/* free an array of points */
/***************************/
void gl_free_points(struct point_pak *p)
{
g_return_if_fail(p->num_points != 0);

if (p->num_points)
  {
  g_free(p->x[0]);
  g_free(p->x[1]);
  g_free(p->x[2]);
  p->num_points = 0;
  }
}

/***********************************/
/* more efficient point generation */
/***********************************/
#define DEBUG_GLO_INIT_SPHERE 0
void gl_init_sphere(struct point_pak *p, struct model_pak *model)
{
gint i, j, n;
gdouble layers, divisions, alpha, beta, phi, theta, sp, cp, st, ct;
gdouble x[3], y[3], mat[9];
struct camera_pak *camera;

g_assert(model != NULL);
g_assert(model->camera != NULL);

/* number of layers - azimuthal sweep */
layers = sysenv.render.sphere_quality+1.0;
alpha = 0.5*PI/layers;

/* alloc */
n=0;
for (i=0 ; i<(gint) 1.0+layers ; i++)
  n += i;
n *= 6;
n++;

p->x[0] = g_malloc(n * sizeof(gdouble));
p->x[1] = g_malloc(n * sizeof(gdouble));
p->x[2] = g_malloc(n * sizeof(gdouble));
p->num_points = n;
p->type = sysenv.render.sphere_quality;

/* TODO - make it use camera_init() as this is just the initial v->x vector */
VEC3SET(x, 0.0, -1.0, 0.0);

/* orient so hemisphere is facing the camera */
camera = model->camera;
if (camera->mode == LOCKED)
  quat_rotate(x, camera->q);
else
  {
/* in the case of a free camera - compute alignment transformation */
  ARR3SET(y, camera->v);
  VEC3MUL(y, -1.0);
  matrix_v_alignment(mat, x, y);
  vecmat(mat, x);
  }

*(p->x[0]+0) = x[0];
*(p->x[1]+0) = x[1];
*(p->x[2]+0) = x[2];

n=1;

for (i=1 ; i<(gint) 1.0+layers ; i++)
  {
#if DEBUG_GLO_INIT_SPHERE
printf("[%d/%d]\n", i, (gint) layers);
#endif

  phi = i*alpha;
  sp = tbl_sin(phi);
  cp = tbl_cos(phi);

/* divisions at this level - hexagonal sweep */
  divisions = i*6.0;
  beta = 2.0*PI/divisions;
  for (j=0 ; j<(gint) divisions ; j++)
    {
    theta = j*beta;
    st = tbl_sin(theta);
    ct = tbl_cos(theta);

/* NEW - quat camera */
    VEC3SET(x, st*sp, -cp, ct*sp);

/* FIXME - FREE camera mode case */
    if (camera->mode == LOCKED)
      quat_rotate(x, camera->q);
    else
      vecmat(mat, x);

    *(p->x[0]+n) = x[0]; 
    *(p->x[1]+n) = x[1];
    *(p->x[2]+n) = x[2];
    n++;

#if DEBUG_GLO_INIT_SPHERE
printf("%7.4f %7.4f %7.4f\n", ct*sp, st*sp, -cp);
#endif
    }
  }
g_assert(p->num_points == n);
}

/**********************/
/* draw half a sphere */
/**********************/
void gl_draw_hemisphere(struct point_pak *p, gdouble *x, gdouble rad)
{
gint i, j, div, m, n, mstart, nstart;
gdouble vec[3];

g_return_if_fail(p->num_points != 0);

/* top cap - split drawing to minimize vertices */
glBegin(GL_TRIANGLE_FAN);
vec[0] = *(p->x[0]);
vec[1] = *(p->x[1]);
vec[2] = *(p->x[2]);

glNormal3dv(vec);
VEC3MUL(vec, rad);

ARR3ADD(vec, x);
glVertex3dv(vec);
for (i=6 ; i-- ; )
  {
  vec[0] = *(p->x[0]+i+1);
  vec[1] = *(p->x[1]+i+1);
  vec[2] = *(p->x[2]+i+1);

  glNormal3dv(vec);
  VEC3MUL(vec, rad);
  ARR3ADD(vec, x);
  glVertex3dv(vec);
  }
vec[0] = *(p->x[0]+6);
vec[1] = *(p->x[1]+6);
vec[2] = *(p->x[2]+6);

glNormal3dv(vec);
VEC3MUL(vec, rad);
ARR3ADD(vec, x);
glVertex3dv(vec);
glEnd();

/* remaining layers */
n=7;
mstart = 1;
for (i=2 ; i<p->type+2 ; i++)
  {
  nstart = n;
  div = i*6.0;
  m = mstart;

#if DEBUG_GLO_DRAW_SPHERE
printf("%d : [m,n start = %d,%d]\n", i, mstart, nstart);
#endif

/* strip loop (theta) */
  glBegin(GL_TRIANGLE_STRIP);
  for (j=0 ; j<div ; j++)
    {
#if DEBUG_GLO_DRAW_SPHERE
printf("> : m,n = %d,%d\n", m, n);
#endif

/* lower */
    vec[0] = *(p->x[0]+n);
    vec[1] = *(p->x[1]+n);
    vec[2] = *(p->x[2]+n);

    glNormal3dv(vec);
    VEC3MUL(vec, rad);
    ARR3ADD(vec, x);
    glVertex3dv(vec);
    n++;

/* upper */
    vec[0] = *(p->x[0]+m);
    vec[1] = *(p->x[1]+m);
    vec[2] = *(p->x[2]+m);

    glNormal3dv(vec);
    VEC3MUL(vec, rad);
    ARR3ADD(vec, x);
    glVertex3dv(vec);

/* staggered doubling */
    m++;
    if ((j+1) % i == 0)
      m--;
/* cycle around */
    if (m >= nstart)
      {
      m -= nstart;
      m += mstart;
      }
    }

/* finish off the triangle layer */
  vec[0] = *(p->x[0]+nstart);
  vec[1] = *(p->x[1]+nstart);
  vec[2] = *(p->x[2]+nstart);

  glNormal3dv(vec);
  VEC3MUL(vec, rad);
  ARR3ADD(vec, x);
  glVertex3dv(vec);

  vec[0] = *(p->x[0]+mstart);
  vec[1] = *(p->x[1]+mstart);
  vec[2] = *(p->x[2]+mstart);

  glNormal3dv(vec);
  VEC3MUL(vec, rad);
  ARR3ADD(vec, x);
  glVertex3dv(vec);

/* finished current layer */
  glEnd();

  mstart = nstart;
  }
#if DEBUG_GLO_DRAW_SPHERE
printf("n : %d\n", n);
#endif
}

/*********************************/
/* more efficient sphere drawing */
/*********************************/
#define DEBUG_GLO_DRAW_SPHERE 0
void gl_draw_sphere(struct point_pak *p, gdouble *x, gdouble rad)
{
gl_draw_hemisphere(p, x, rad);
}

/***************************/
/* create digitized circle */
/***************************/
void gl_init_circle(struct point_pak *circle, gint segments, struct model_pak *model)
{
gint i;
gdouble da, theta;
gdouble n[3], x[3], y[3], mat[9];
struct camera_pak *camera;

g_assert(model != NULL);
g_assert(model->camera != NULL);
g_return_if_fail(segments != 0);

camera = model->camera;

/* allocate for digitized circle */
circle->x[0] = g_malloc((segments+1) * sizeof(gdouble));
circle->x[1] = g_malloc((segments+1) * sizeof(gdouble));
circle->x[2] = g_malloc((segments+1) * sizeof(gdouble));
circle->num_points = segments + 1;

/* first point */
VEC3SET(x, 1.0, 0.0, 0.0);
if (camera->mode == LOCKED)
  quat_rotate(x, camera->q);
else
  {
/* in the case of a free camera - compute alignment transformation */
  VEC3SET(n, 0.0, -1.0, 0.0);
  ARR3SET(y, camera->v);
  matrix_v_alignment(mat, n, y);
  vecmat(mat, x);
  }

*(circle->x[0]) = x[0];
*(circle->x[1]) = x[1];
*(circle->x[2]) = x[2];

/* sweep in x,z */
da = 2.0 * PI / (gdouble) segments;
for (i=1 ; i<segments ; i++) 
  {
  theta = (gdouble) i * da;
/* store point on rim (doubles as the normal) */
  x[0] = tbl_cos(theta);
  x[1] = 0.0;
  x[2] = tbl_sin(theta);

  if (camera->mode == LOCKED)
    quat_rotate(x, camera->q);
  else
    vecmat(mat, x);

  *(circle->x[0]+i) = x[0];
  *(circle->x[1]+i) = x[1];
  *(circle->x[2]+i) = x[2];
  }

/* duplicate first point, so the circle is closed */
*(circle->x[0]+segments) = *(circle->x[0]);
*(circle->x[1]+segments) = *(circle->x[1]);
*(circle->x[2]+segments) = *(circle->x[2]);
}

/*******************************/
/* draw a pre-digitized circle */
/*******************************/
void gl_draw_circle(struct point_pak *circle, gdouble *x, gdouble rad)
{
gint i;
gdouble vec[3];

g_return_if_fail(circle->num_points != 0);

glBegin(GL_POLYGON);
for (i=0 ; i<circle->num_points ; i++) 
  {
/* get points */
  vec[0] = *(circle->x[0]+i);
  vec[1] = *(circle->x[1]+i);
  vec[2] = *(circle->x[2]+i);
  VEC3MUL(vec, rad);
  ARR3ADD(vec, x);
  glVertex3dv(vec);
  }
glEnd();
}

/*****************************/
/* draw a pre-digitized ring */
/*****************************/
void gl_draw_ring(struct point_pak *circle, gdouble *x, gdouble rad1, gdouble rad2)
{
gint i;
gdouble vec[3];

g_return_if_fail(circle->num_points != 0);

glBegin(GL_QUAD_STRIP);
for (i=0 ; i<circle->num_points ; i++) 
  {
/* inner point */
  vec[0] = *(circle->x[0]+i);
  vec[1] = *(circle->x[1]+i);
  vec[2] = *(circle->x[2]+i);
  VEC3MUL(vec, rad2);
  ARR3ADD(vec, x);
  glVertex3dv(vec);

/* outer point */
  vec[0] = *(circle->x[0]+i);
  vec[1] = *(circle->x[1]+i);
  vec[2] = *(circle->x[2]+i);
  VEC3MUL(vec, rad1);
  ARR3ADD(vec, x);
  glVertex3dv(vec);
  }
glEnd();
}

/******************************************/
/* draw a cylinder between points va & vb */
/******************************************/
/* NB: used for drawing vectors */
void gl_draw_cylinder(gdouble *va, gdouble *vb, gdouble rad, guint segments)
{
gint i;
gdouble theta, st, ct;
gdouble v1[3], v2[3], v[3], v12[3], p[3], q[3];

g_return_if_fail(segments != 0);

/* force ordering of v1 & v2 to get correct quad normals */
if (va[2] > vb[2])
  {
  ARR3SET(v1, vb);
  ARR3SET(v2, va);
  }
else
  {
  ARR3SET(v1, va);
  ARR3SET(v2, vb);
  }

/* vector from v1 to v2 */
ARR3SET(v12, v2);
ARR3SUB(v12, v1);

/* make a guess at a vector with some orthogonal component to v12 */
VEC3SET(p, 1.0, 1.0, 1.0);
crossprod(q, p, v12);

/* our guess was bad - so fix it */
if (VEC3MAGSQ(q) < FRACTION_TOLERANCE)
  {
  VEC3SET(p, 0.0, 1.0, 0.0);
  crossprod(q, p, v12);
  }
crossprod(p, v12, q);

/* p and q are orthonormal and in the plane of the cylinder's cross section */
normalize(p, 3);
normalize(q, 3);

/* build the cylinder from rectangular segments */
/* TODO - see if vertex arrays give any speedup here */
glBegin(GL_QUAD_STRIP);
for (i=segments+1 ; i-- ; )
  {
/* sweep out a circle */
  theta = (gdouble) i * 2.0 * PI / (gdouble) segments;
  st = tbl_sin(theta);
  ct = tbl_cos(theta);
/* construct normal */
  v12[0] = ct * p[0] + st * q[0];
  v12[1] = ct * p[1] + st * q[1];
  v12[2] = ct * p[2] + st * q[2];
/* set the normal for the two subseqent points */
  glNormal3dv(v12);
/* get the vector from centre to rim */
  VEC3MUL(v12, rad);
/* point on disk 1 */
  ARR3SET(v, v2);
  ARR3ADD(v, v12);
  glVertex3dv(v);
/* point on disk 2 */
  ARR3SET(v, v1);
  ARR3ADD(v, v12);
  glVertex3dv(v);
  }
glEnd();
}

/*************************************/
/* draw a cone defined by two points */
/*************************************/
void draw_cone(gdouble *va, gdouble *vb, gdouble base, gint segments)
{
gint i;
gdouble theta, st, ct, da;
gdouble v[3], p[3], q[3], v12[3];

g_return_if_fail(segments != 0);

/* vector from va to vb */
ARR3SET(v12, vb);
ARR3SUB(v12, va);

/* create normal vectors p and q, co-planar with the cylinder's cross-sectional disk */

/* make a guess at a vector with some orthogonal component to v12 */
VEC3SET(p, 1.0, 1.0, 1.0);
crossprod(q, p, v12);

/* if our guess was bad - fix it */
if (VEC3MAGSQ(q) < FRACTION_TOLERANCE)
  {
  VEC3SET(p, 0.0, 1.0, 0.0);
  crossprod(q, p, v12);
  }
crossprod(p, v12, q);

normalize(p, 3);
normalize(q, 3);

/* build the cone */
da = 2.0 * PI / (gdouble) segments;
glBegin(GL_TRIANGLE_FAN);
glVertex3dv(vb);
for (i=segments+1 ; i-- ; )
  {
/* sweep out a circle */
  theta = (gdouble) i * da;
  st = tbl_sin(theta);
  ct = tbl_cos(theta);
/* construct normal */
  v12[0] = ct * p[0] + st * q[0];
  v12[1] = ct * p[1] + st * q[1];
  v12[2] = ct * p[2] + st * q[2];
  glNormal3dv(v12);
/* point 1 on disk 1 */
  ARR3SET(v, va);
  v[0] += base*v12[0];
  v[1] += base*v12[1];
  v[2] += base*v12[2];
  glVertex3dv(v);
  }
glEnd();
}

/*********************************************/
/* draw a vector cylinder between two points */
/*********************************************/
void draw_vector(gdouble *va, gdouble *vb, gdouble size)
{
gdouble vec[3];

gl_draw_cylinder(va, vb, size, 5);

ARR3SET(vec, vb);
ARR3SUB(vec, va);
normalize(vec, 3);
VEC3MUL(vec, 4.0*size);
ARR3ADD(vec, vb);

draw_cone(vb, vec, 3*size, 9);
}

/*********************************************************/
/* draw a circular arc from p1 to p2, with p0 the centre */
/*********************************************************/
void draw_arc(gdouble *p0, gdouble *p1, gdouble *p2)
{
gint i, num_points;
gdouble angle, theta, st, ct;
gdouble v[3], v1[3], v2[3], len1, len2;

/* get vector to arc start */
ARR3SET(v1, p1);
ARR3SUB(v1, p0);

/* get vector to arc end */
ARR3SET(v2, p2);
ARR3SUB(v2, p0);

/* scale back a bit */
VEC3MUL(v1, 0.5);
VEC3MUL(v2, 0.5);

/* get lengths */
len1 = VEC3MAG(v1);
len2 = VEC3MAG(v2);

/* make them both equal to the smaller */
if (len2 > len1)
  {
  VEC3MUL(v2, len1/len2);
  }
else
  {
  VEC3MUL(v1, len2/len1);
  }

angle = via(v1, v2, 3);

num_points = 32;

glBegin(GL_LINE_STRIP);
/* sweep from start to end */
for (i=0 ; i<=num_points ; i++)
  {
/* sweep out a circle */
  theta = (gdouble) i * 0.5 * PI / (gdouble) num_points;
  st = tbl_sin(theta);
  ct = tbl_cos(theta);
/* construct normal */
  ARR3SET(v, p0);
  v[0] += ct * v1[0] + st * v2[0];
  v[1] += ct * v1[1] + st * v2[1];
  v[2] += ct * v1[2] + st * v2[2];
  glVertex3dv(v);
  }
glEnd();
}

