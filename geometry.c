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
#include "matrix.h"
#include "morph.h"
#include "spatial.h"
#include "surface.h"
#include "model.h"
#include "interface.h"

/* main structure */
extern struct sysenv_pak sysenv;

/****************************/
/* project vector x along v */
/****************************/
void vector_v_project(gdouble *x, gdouble *v)
{
gdouble len, n[3], y[3];

ARR3SET(n, v);
normalize(n, 3);

ARR3SET(y, n);
ARR3MUL(y, x);
len = y[0] + y[1] + y[2];

ARR3SET(x, n);
VEC3MUL(x, len);
}

/**************************************/
/* compute vector from point to plane */
/**************************************/
/* vector/point (v,x) are 3D, plane (n) is 4D */
void vector_point2plane(gdouble *v, gdouble *x, gdouble *n)
{
gdouble dp, y[3];

/* get point on the plane */
ARR3SET(y, n);
VEC3MUL(y, n[3]);
/* form vector to input point */
ARR3SUB(y, x);

/* dotprod with normal */
ARR3MUL(y, n);
dp = y[0] + y[1] + y[2];

/* scale the normal */
ARR3SET(v, n);
VEC3MUL(v, dp);
}

/***********************************************************/
/* check if point is inside the halfspace defined by plane */
/***********************************************************/
/* x is 3D n is 4D */
gint point_outside(gdouble *x, gdouble *n)
{
gdouble v[3];

vector_point2plane(v, x, n);

if (via(v, n, 3) > 0.5*G_PI)
  return(TRUE);
return(FALSE);
}

/************************************/
/* compute intersection of 3 planes */
/************************************/
/* n1 n2 n3 are 4D - the plane normals (3) plus the distance to the plane (1) */
/* ie n[0].x + n[1].y + n[2].z = n[3] */
#define DEBUG_PLANE_INTERSECT 0
gint plane_intersect(gdouble *x, gdouble *n1, gdouble *n2, gdouble *n3)
{
gdouble d, dp1, dp2, dp3;
gdouble x1[3], x2[3], x3[3];
gdouble a1[3], a2[3], a3[3];
gdouble mat[9];

/* compute determinant of normal matrix */
mat[0] = n1[0];
mat[1] = n2[0];
mat[2] = n3[0];
mat[3] = n1[1];
mat[4] = n2[1];
mat[5] = n3[1];
mat[6] = n1[2];
mat[7] = n2[2];
mat[8] = n3[2];

d = matrix_determinant(mat);

#if DEBUG_PLANE_INTERSECT
P3MAT("matrix: ", mat);
#endif

if (fabs(d) < FRACTION_TOLERANCE*FRACTION_TOLERANCE)
  {
/*
printf("planes do not intersect.\n");
*/
  return(FALSE);
  }

/* compute actual points on the planes */
ARR3SET(x1, n1);
ARR3SET(x2, n2);
ARR3SET(x3, n3);
VEC3MUL(x1, n1[3]);
VEC3MUL(x2, n2[3]);
VEC3MUL(x3, n3[3]);

/* compute intersection point */
crossprod(a1, n2, n3);
crossprod(a2, n3, n1);
crossprod(a3, n1, n2);

ARR3MUL(x1, n1);
dp1 = x1[0] + x1[1] + x1[2];
ARR3MUL(x2, n2);
dp2 = x2[0] + x2[1] + x2[2];
ARR3MUL(x3, n3);
dp3 = x3[0] + x3[1] + x3[2];

VEC3MUL(a1, dp1);
VEC3MUL(a2, dp2);
VEC3MUL(a3, dp3);

ARR3ADD(a1, a2);
ARR3ADD(a1, a3);

VEC3MUL(a1, 1.0/d);

#if DEBUG_PLANE_INTERSECT
P3VEC("intersection: ", a1);
#endif

ARR3SET(x, a1);

return(TRUE);
}

/***********************************/
/* initialize the halfspace matrix */
/***********************************/
/* m is the number of planes, returns the number of valid planes */
#define DEBUG_MORPH_SETUP 0
gint morph_setup(gint m, gint **hkl, gdouble **mat, struct model_pak *model)
{
gint i, n;
gdouble max, value;
GSList *list;
struct plane_pak *plane;

/* feed planes into matrix */
max = 0.0;
i=0;
for (list=model->planes ; list ; list=g_slist_next(list))
  {
  plane = list->data;

/* enforce update of best energy value */
  update_plane_energy(plane, model);

/* NB - skip invalid values (eg non -ve attachment energy / zero length planes) */
  switch (model->morph_type)
    {
    case DHKL:
      value = 1.0/fabs(plane->dhkl);
      break;

    case EQUIL_UN:
      if (plane->esurf[0] <= 0.0)
        continue;
      value = plane->esurf[0];
      break;

    case EQUIL_RE:
      if (plane->esurf[1] <= 0.0)
        continue;
      value = plane->esurf[1];
      break;

    case GROWTH_UN:
      if (plane->eatt[0] >= 0.0)
        continue;
      value = -plane->eatt[0];
      break;

    case GROWTH_RE:
      if (plane->eatt[1] >= 0.0)
        continue;
      value = -plane->eatt[1];
      break;

    case MORPH_BBPA:
      if (!plane->bbpa)
        continue;
      value = (gdouble) plane->bbpa;
      break;

    default:
      gui_text_show(ERROR, "Bad morphology type.\n");
      return(0);
    }

/* NEW */
  if (model->morph_type == MORPH_BBPA)
    {
    printf("(%2d %2d %2d)", plane->index[0], plane->index[1], plane->index[2]);
    printf(" : cleaved = %.1f", plane->bbpa*plane->area);
    printf(" : area = %f", plane->area);
    printf(" : bbpa = %f", plane->bbpa);
    printf("\n");
    }

  if (fabs(value) > max)
    max = fabs(value);

/* fill out the inequality */
  mat[i] = g_malloc(4*sizeof(gdouble));
  hkl[i] = g_malloc(3*sizeof(gint));

  hkl[i][0] = plane->index[0];
  hkl[i][1] = plane->index[1];
  hkl[i][2] = plane->index[2];

  mat[i][0] = -plane->x[0];
  mat[i][1] = -plane->x[1];
  mat[i][2] = -plane->x[2];
  mat[i][3] = value;
  i++;

  g_assert(i <= m);
  }
n = i;

/* scale as morph is overlayed on unit cells (if present) */
/*
for (i=n ; i-- ; )
  {
  mat[i][3] *= model->rmax;
  mat[i][3] /= max;
  }
*/

if (!n)
  {
  gui_text_show(ERROR, "No valid planes found.\n");
  return(0);
  }

#if DEBUG_MORPH_SETUP
printf("input inequalities: %d/%d\n", n, m);
for (i=0 ; i<n ; i++)
  {
  P4VEC("", mat[i]);
  }
printf("\n\n");
#endif

return(n);
}

/*************************************************************/
/* build a polyhedron from a list of halfspace intersections */
/*************************************************************/
#define MORPH_BUILD 0
gint morph_build(struct model_pak *model)
{
gint i, j, k, l, m, flag, *absent, **hkl;
gdouble x[3], **mat;
GSList *vlist, *list1, *list2;
struct plane_pak *plane;
struct vertex_pak *v1, *v2;
struct spatial_pak *spatial;

g_assert(model != NULL);

spatial_destroy_by_type(SPATIAL_MORPHOLOGY, model);

/* initialize */
m = g_slist_length(model->planes);

/* CHECK - g_malloc0 initialize to NULLs for a pointer array? */
mat = g_malloc0(m * sizeof(gdouble *));
hkl = g_malloc0(m * sizeof(gint *));

m = morph_setup(m, hkl, mat, model);
absent = g_malloc0(m * sizeof(gint));

/* compute plane intersections */
for (i=0 ; i<m ; i++)
  {
  vlist = NULL;

/* match ith plane against all (unique) j,k combinations */
  for (j=0 ; j<m-1 ; j++)
    {
    if (i==j || absent[j])
      continue;

    for (k=j+1 ; k<m ; k++)
      {
      if (k==i || absent[k])
        continue;

/* get intersection point */
      if (plane_intersect(x, mat[i], mat[j], mat[k]))
        {
        flag = 1;
        for (l=m ; l-- ; )
          {
          if (l==i || l==j || l==k || absent[l])
            continue;

          if (point_outside(x, mat[l]))
            {
            flag = 0;
            break;
            }
          }
        if (flag)
          {
          v1 = g_malloc(sizeof(struct vertex_pak));
          vlist = g_slist_prepend(vlist, v1);
          ARR3SET(v1->x, x);
          ARR3SET(v1->n, mat[i]);
          vecmat(model->ilatmat, v1->x);
          vecmat(model->ilatmat, v1->n);
          v1->adj = NULL;
          v1->adj = g_slist_prepend(v1->adj, GINT_TO_POINTER(j));
          v1->adj = g_slist_prepend(v1->adj, GINT_TO_POINTER(k));
          }
        }
      }
    }

/* eliminate redundent vertices */
/* FIXME - this can leave gaps by deleting very skinny polygons */
/* FIXME - BUT removing this code results in errors in the vertex ordering */
  list1 = vlist;
  while (list1)
    {
    v1 = list1->data;
    list2 = g_slist_next(list1);
    while (list2)
      {
      v2 = list2->data;
      list2 = g_slist_next(list2);
      ARR3SET(x, v1->x);
      ARR3SUB(x, v2->x);
      if (VEC3MAGSQ(x) < FRACTION_TOLERANCE)
        vlist = g_slist_remove(vlist, v2);
      }
    list1 = g_slist_next(list1);
    }

/* build a polygon if we have at least 3 points */
  if (g_slist_length(vlist) > 2)
    {
#if MORPH_BUILD
printf("facet %d: %d vertices\n", i, g_slist_length(vlist));
for (list1=vlist ; list1 ; list1=g_slist_next(list1))
  {
  v1 = list1->data;
  P3VEC("x: ", v1->x);
  }
#endif

/* this should do generate the correct vertex ordering for the polygon */
    spatial = spatial_build_facet(mat[i], vlist, model);
    free_slist(vlist);

/* generate hkl label */
    plane = g_slist_nth_data(model->planes, i);
    if (plane)
      spatial->label = g_strdup_printf("(%d%d%d)", hkl[i][0], hkl[i][1], hkl[i][2]);
    }
  else
    absent[i] = TRUE;
  }

/* cleanup */
for (i=m ; i-- ; )
  {
  g_free(hkl[i]);
  g_free(mat[i]);
  }
g_free(hkl);
g_free(mat);
g_free(absent);

return(0);
}

