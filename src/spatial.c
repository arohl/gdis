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
#include "geometry.h"
#include "matrix.h"
#include "morph.h"
#include "numeric.h"
#include "spatial.h"
#include "select.h"
#include "interface.h"
#include "gui_shorts.h"
#ifdef __APPLE__
#include <OpenGL/gl.h>
#else
#include "opengl.h"
#include <GL/gl.h>
#endif

/* main structure */
extern struct sysenv_pak sysenv;

/* handle all the vector/plane creation, intersections etc. */
/* add all the vector normalization/crossprod etc. stuff as well */

/**************************/
/* spatial free primitive */
/**************************/
void spatial_free(gpointer data)
{
GSList *list;
struct spatial_pak *spatial=data;

g_assert(spatial != NULL);

for (list=spatial->list ; list ; list=g_slist_next(list))
  {
  struct vec_pak *v = list->data;
  g_free(v->data);
  g_free(v);
  }

g_slist_free(spatial->list);
g_free(spatial->data);
g_free(spatial);
}

/*****************************/
/* destroy spatial primitive */
/*****************************/
void spatial_destroy(gpointer data, struct model_pak *model)
{
struct spatial_pak *spatial=data;

if (spatial && model)
  {
  model->spatial = g_slist_remove(model->spatial, spatial);
  spatial_free(spatial);
  }
}

/*********************************/
/* destroy all spatial primitive */
/*********************************/
void spatial_destroy_all(struct model_pak *model)
{
GSList *list;
struct spatial_pak *spatial;

g_assert(model != NULL);

list = model->spatial;
while (list)
  {
  spatial = list->data;
  list = g_slist_next(list);
  spatial_free(spatial);
  }

g_slist_free(model->spatial);
model->spatial = NULL;
}

/**************************/
/* destroy by type number */
/**************************/
void spatial_destroy_by_type(gint type, struct model_pak *model)
{
GSList *list;
struct spatial_pak *spatial;

list = model->spatial;
while (list)
  {
  spatial = list->data;
  list = g_slist_next(list);
  if (type == spatial->type)
    spatial_destroy(spatial, model);
  }
}

/***********************************/
/* destroy spatial by id primitive */
/***********************************/
void spatial_destroy_by_label(const gchar *label, struct model_pak *model)
{
gint a, b;
GSList *list;
struct spatial_pak *spatial;

g_assert(label != NULL);
g_assert(model != NULL);

a = strlen(label);

list = model->spatial;
while (list)
  {
  spatial = list->data;
  list = g_slist_next(list);

  b = strlen(spatial->label);
  if (a == b)
    if (g_ascii_strncasecmp(label, spatial->label, a) == 0)
      spatial_destroy(spatial, model);
  }
}

/************************************/
/* primitives for creating spatials */
/************************************/
gpointer spatial_new(const gchar *label,
                     gint type,
                     gint size,
                     gint periodic,
                     struct model_pak *model)
{
static gint n=0;
struct spatial_pak *spatial;

g_assert(model != NULL);

spatial = g_malloc(sizeof(struct spatial_pak));
model->spatial = g_slist_append(model->spatial, spatial);

if (label)
  spatial->label = g_strdup(label);
else
  {
/* FIXME - better facility for creating a unique name */
  spatial->label = g_strdup_printf("unknown_%d", n++);
  }

spatial->type = type;
spatial->size = size;
spatial->periodic = periodic;
spatial->show_label = FALSE;
spatial->data = NULL;
spatial->list = NULL;

ARR3SET(spatial->c, sysenv.render.fg_colour);

switch (size)
  {
  case 1:
    spatial->method = GL_POINTS;
    spatial->material = SPATIAL_LINE;
    break;

  case 2:
    spatial->method = GL_LINES;
    spatial->material = SPATIAL_LINE;
    break;

  case 3:
    spatial->method = GL_TRIANGLES;
    spatial->material = SPATIAL_SURFACE;
    break;

  case 4:
    spatial->method = GL_QUADS;
    spatial->material = SPATIAL_SURFACE;
    break;

  default:
    spatial->method = GL_POLYGON;
    spatial->material = SPATIAL_SURFACE;
  }

/* special case - composite objects */
switch (type)
  {
  case SPATIAL_VECTOR:
    spatial->method = -1;
    spatial->material = SPATIAL_SOLID;
  }

return(spatial);
}

/************************************************************/
/* associates a general purpose data pointer with a spatial */
/************************************************************/
void spatial_data_set(gpointer data, gpointer ptr)
{
struct spatial_pak *spatial=ptr;

spatial->data = data;
}

/***********************************************************/
/* retrieves a general purpose data pointer from a spatial */
/***********************************************************/
gpointer spatial_data_get(gpointer ptr)
{
struct spatial_pak *spatial=ptr;

return(spatial->data);
}

/*************************************/
/* structure manipulation primitives */
/*************************************/
void spatial_label_show(gpointer data)
{
struct spatial_pak *spatial=data;

spatial->show_label = TRUE;
}

void spatial_label_hide(gpointer data)
{
struct spatial_pak *spatial=data;

spatial->show_label = FALSE;
}

/*****************************/
/* vertex creation primitive */
/*****************************/
gpointer vertex_new(gdouble *x, gdouble *n, gdouble *c)
{
struct vec_pak *v;

v = g_malloc(sizeof(struct vec_pak));

if (x)
  {
  ARR3SET(v->x, x);
  }
else
  {
  VEC3SET(v->x, 0.0, 0.0, 0.0);
  }
if (n)
  {
  ARR3SET(v->n, n);
  }
else
  {
  VEC3SET(v->n, 0.0, 0.0, 1.0);
  }
if (c)
  {
  ARR3SET(v->colour, c);
  }
else
  {
  ARR3SET(v->colour, sysenv.render.fg_colour);
  }

v->data = NULL;

return(v);
}

/**************************/
/* add vertex with colour */
/**************************/
void spatial_vertex_add(gdouble *x, gdouble *colour, gpointer ptr)
{
struct spatial_pak *spatial = ptr;
struct vec_pak *vec;

vec = vertex_new(x, NULL, colour);

spatial->list = g_slist_prepend(spatial->list, vec);
}

/*************************************/
/* add vertex with normal and colour */
/*************************************/
void spatial_vnorm_add(gdouble *x, gdouble *n, gdouble *colour, gpointer ptr)
{
struct spatial_pak *spatial = ptr;
struct vec_pak *vec;

vec = vertex_new(x, n, colour);

spatial->list = g_slist_prepend(spatial->list, vec);
}

/********************************************/
/* add vertex with additional property data */
/********************************************/
void spatial_vdata_add(gdouble *x, gdouble *n, gdouble *colour, gpointer data, gpointer ptr)
{
struct spatial_pak *spatial = ptr;
struct vec_pak *vec;

vec = vertex_new(x, n, colour);

vec->data = data;

spatial->list = g_slist_prepend(spatial->list, vec);
}

/****************************************************/
/* compute normal for the plane containing 3 points */
/****************************************************/
void compute_normal(gdouble *n, gdouble *c1, gdouble *c2, gdouble *c3)
{
gdouble v0[3], v1[3], v2[3];

ARR3SET(v0, c1);
ARR3SET(v1, c2);
ARR3SET(v2, c3);

/* make vectors 0->1 & 0->2 */
ARR3SUB(v1, v0);
ARR3SUB(v2, v0);

/* get normal */
crossprod(n, v1, v2);
normalize(n, 3);
}

/************************************************************/
/* compute normal for a facet defined by a list of vertices */
/************************************************************/
void compute_facet_normal(gdouble *n, GSList *vlist)
{
gint i, j;
struct vertex_pak *v1, *v2, *v3;

g_assert(g_slist_length(vlist) > 2);

j = g_slist_length(vlist);

printf("[size = %d] -------------------------------\n", j);
for (i=0 ; i<j-2 ; i++)
  {
v1 = g_slist_nth_data(vlist, i);
v2 = g_slist_nth_data(vlist, i+1);
v3 = g_slist_nth_data(vlist, i+2);

compute_normal(n, v1->x, v2->x, v3->x);

P3VEC("n: ", n);
  }


}

/***********************************************/
/* compute best fit plane for a list of atoms */
/**********************************************/
#define DEBUG_COMPUTE_PLANE 0
void compute_plane(gdouble *res, GList *cores)
{
gdouble tmp[3], vec1[3], vec2[3];
struct core_pak *core1, *core2, *core3;

/* FIXME - this only uses the 1st 3 points, in future - an average for n pts */
/* http://www.geocities.com/CollegePark/Campus/1278/JBook/LSQ_Plane.html */
/* shouldn't matter when this routine is properly implemented, as we'll go */
/* through the atom_list in the usual fashoin and terminate when it's NULL, */
/* not when its data ptr is NULL */ 
g_assert(cores != NULL);
g_assert(g_list_length(cores) > 2);

core1 = g_list_nth_data(cores, 0);
core2 = g_list_nth_data(cores, 1);
core3 = g_list_nth_data(cores, 2);

#if DEBUG_COMPUTE_PLANE
printf("computing normal from atoms: %p, %p, %p\n", core1, core2, core3);
#endif

ARR3SET(tmp, core1->x);
ARR3SET(vec1, core2->x);
ARR3SET(vec2, core3->x);

/* make vectors 0->1 & 0->2 */
ARR3SUB(vec1, tmp);
ARR3SUB(vec2, tmp);

/* get normal */
crossprod(res, vec1, vec2);
}

/******************************************/
/* get vertex from list with common facet */
/******************************************/
gpointer vertex_common_facet(struct vertex_pak *v, GSList *vlist)
{
GSList *list1, *list2;
struct vertex_pak *v1;

for (list1=vlist ; list1 ; list1=g_slist_next(list1))
  {
  v1 = list1->data;
  
  for (list2=v->adj ; list2 ; list2=g_slist_next(list2))
    {
    if (g_slist_find(v1->adj, list2->data))
      return(v1);
    }
  }
return(NULL);
}

/*****************************************/
/* compute the angle between two vectors */
/*****************************************/
/* TODO - replace all via() calls with this */
gdouble vector_angle(gdouble *vec1, gdouble *vec2, gint dim)
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

/**********************************************************/
/* order a set of input coplanar vertices to be clockwise */
/**********************************************************/
/* n is the surface normal */
#define DEBUG_SPATIAL_ORDER_FACET 0
GSList *spatial_order_facet(gdouble *norm, GSList *source, struct model_pak *model)
{
gint new, old;
gdouble total, angle, min;
gdouble a[3], b[3], c[3], m[3], n[3];
GSList *list, *pool, *sorted;
struct vertex_pak *v1, *v2, *v3;

g_assert (model != NULL);

/* convert input miller to cartesian space normal */
ARR3SET(n, norm);
vecmat(model->latmat, n);
normalize(n, 3);

/* compute midpoint (acts as orientation vector) */
total = 0.0;
VEC3SET(m, 0.0, 0.0, 0.0);
for (list=source ; list ; list=g_slist_next(list))
  {
  v1 = list->data;

  ARR3ADD(m, v1->x);
  total++;
  }
VEC3MUL(m, 1.0/total);

/* vertex pool */
pool = g_slist_copy(source);
old = g_slist_length(source);
v1 = pool->data;
pool = g_slist_remove(pool, v1);

/* build a sorted list */
sorted = NULL;
sorted = g_slist_append(sorted, v1);
new = 1;

/* loop until we assign all source vertices */
while (new < old)
  {
/* search the current pool for correct normal orientations */
  min = 2.0*G_PI;
  v3 = NULL;

/* vector from midpoint to current reference vertex */
  ARR3SET(a, v1->x);
  ARR3SUB(a, m);
  for (list=pool ; list ; list=g_slist_next(list))
    {
    v2 = list->data;
    ARR3SET(b, v2->x);
    ARR3SUB(b, m);
    crossprod(c, a, b);

/* consider vertex only if it satisfies clockwise condition */
/*
    if (vector_angle(c, n, 3) < 0.1)
*/
    if (vector_angle(c, n, 3) < 0.5*G_PI)
      {
/* search for minimum angle separation */
      angle = vector_angle(a, b, 3);
      if (angle < min)
        {
        min = angle;
        v3 = v2;
        }
      }
    }

/* best vertex is now the reference */
  if (v3)
    {
    sorted = g_slist_prepend(sorted, v3);
    pool = g_slist_remove(pool, v3);
    v1 = v3;
    }
  else
    gui_text_show(ERROR, "Failed to find next vertex.\n");

  new++;
  }

return(sorted);
}

/***************************************************/
/* build a polygon spatial from a list of vertices */
/***************************************************/
#define DEBUG_SPATIAL_BUILD_FACET 0
gpointer spatial_build_facet(gdouble *n, GSList *source, struct model_pak *model)
{
GSList *list, *vlist;
gpointer spatial;
struct vertex_pak *v;

g_assert(model != NULL);

/* CURRENT */
vlist = spatial_order_facet(n, source, model);

spatial = spatial_new(NULL, SPATIAL_MORPHOLOGY, 0, FALSE, model);

spatial_label_show(spatial);

for (list=vlist ; list ; list=g_slist_next(list))
  {
  v = list->data;
  spatial_vnorm_add(v->x, v->n, sysenv.render.morph_colour, spatial);
  }

g_slist_free(vlist);

return(spatial);
}

/************************************/
/* create a unit cell bounded plane */
/************************************/
#define DEBUG_SPATIAL_CELL_PLANE 0
void spatial_cell_plane(GSList *cores, struct model_pak *model)
{
gint i, j, k, flag;
gdouble a1[3], a2[3], a3[3];
gdouble x1[3], x2[3], x3[3];
gdouble n1[4], n2[4], n3[4];
gdouble walls[24];
gpointer spatial;
#if DEBUG_SPATIAL_CELL_PLANE
GSList *list;
#endif
GSList *vlist;
struct core_pak *core1, *core2, *core3;
struct vertex_pak *v;

g_assert(cores != NULL);
g_assert(model != NULL);

/* get cell vectors */
a1[0] = model->latmat[0];
a1[1] = model->latmat[3];
a1[2] = model->latmat[6];
a2[0] = model->latmat[1];
a2[1] = model->latmat[4];
a2[2] = model->latmat[7];
a3[0] = model->latmat[2];
a3[1] = model->latmat[5];
a3[2] = model->latmat[8];

/* setup bc walls */
crossprod(&walls[0], a2, a3);
walls[3] = model->pbc[0];
normalize(&walls[0], 3);
ARR3SET(&walls[4], &walls[0]);
VEC3MUL(&walls[4], -1.0);
walls[7] = 0.0;
/* get separation between (origin containing) plane to end of a vector */
vector_point2plane(x1, a1, &walls[4]);
walls[3] = VEC3MAG(x1);

/* setup ca walls */
crossprod(&walls[8], a3, a1);
walls[11] = model->pbc[1];
normalize(&walls[8], 3);
ARR3SET(&walls[12], &walls[8]);
VEC3MUL(&walls[12], -1.0);
walls[15] = 0.0;
/* get separation between (origin containing) plane to end of a vector */
vector_point2plane(x1, a2, &walls[12]);
walls[11] = VEC3MAG(x1);

/* setup ab walls */
crossprod(&walls[16], a1, a2);
walls[19] = model->pbc[2];
normalize(&walls[16], 3);
ARR3SET(&walls[20], &walls[16]);
VEC3MUL(&walls[20], -1.0);
walls[23] = 0.0;
/* get separation between (origin containing) plane to end of a vector */
vector_point2plane(x1, a3, &walls[20]);
walls[19] = VEC3MAG(x1);

#if DEBUG_SPATIAL_CELL_PLANE
P3MAT("latmat: ", model->latmat);
for (i=0 ; i<6 ; i++)
  {
  P4VEC("wall: ", &walls[4*i]); 
  }
#endif

core1 = g_slist_nth_data(cores, 0);
core2 = g_slist_nth_data(cores, 1);
core3 = g_slist_nth_data(cores, 2);

g_assert(core1 != NULL);
g_assert(core2 != NULL);
g_assert(core3 != NULL);

/* get cartesian normal */
ARR3SET(x1, core1->x);
ARR3SET(x2, core2->x);
ARR3SET(x3, core3->x);
vecmat(model->latmat, x1);
vecmat(model->latmat, x2);
vecmat(model->latmat, x3);

calc_norm(n1, x1, x2, x3);
normalize(n1, 3);

/* compute distance (from origin) to plane */
/* NB: sign is important */
ARR3SET(x1, core1->x);
vecmat(model->latmat, x1);
ARR3SET(x2, x1);
ARR3MUL(x2, n1);
n1[3] = x2[0] + x2[1] + x2[2];

#if DEBUG_SPATIAL_CELL_PLANE
P4VEC("n1: ", n1);
#endif

/* compute intersections with cell walls */
vlist = NULL;
for (i=0 ; i<5 ; i++)
  {
  for (j=i+1 ; j<6 ; j++)
    {
    ARR4SET(n2, &walls[4*i]);
    ARR4SET(n3, &walls[4*j]);

    if (plane_intersect(x1, n1, n2, n3))
      {
      flag = 0;

/* eliminate intersections outside the unit cell */
      ARR3SET(x2, x1);
      vecmat(model->ilatmat, x2);

      for (k=3 ; k-- ; )
        {
/* NB: precision errors may be a factor here */
        if (x2[k] > 1.0+FRACTION_TOLERANCE || x2[k] < -FRACTION_TOLERANCE)
          flag++;
        }

      if (!flag)
        {

/*
printf("====================================\n");
P4VEC("n2: ", n2);
P4VEC("n3: ", n3);
P4VEC("x1: ", x1);
*/

/* store the intersection point */
        v = g_malloc(sizeof(struct vertex_pak));
        vlist = g_slist_prepend(vlist, v);
        ARR3SET(v->x, x2);

        v->adj = NULL;

/* compute facet adjacencies */
        for (k=0 ; k<6 ; k++)
          {
/* get distance of point to facets */
          vector_point2plane(x2, x1, &walls[4*k]);
          if (VEC3MAGSQ(x2) < FRACTION_TOLERANCE)
            {
/*
printf(" + (%d)\n", k);
*/
            v->adj = g_slist_prepend(v->adj, GINT_TO_POINTER(k));
            }
          }
        }
      }
    }
  }

k = g_slist_length(vlist);

/* order the vertices (via common facets) */
if (k > 2)
  {
  spatial = spatial_new("plane", SPATIAL_GENERIC, 0, TRUE, model);
  v = vlist->data;
  spatial_vertex_add(v->x, NULL, spatial);

#if DEBUG_SPATIAL_CELL_PLANE
P3VEC("v: ", v->x);
for (list=v->adj ; list ; list=g_slist_next(list))
  {
  printf(" + (%d)\n", GPOINTER_TO_INT(list->data));
  }
#endif

vlist = g_slist_remove(vlist, v);
v = vertex_common_facet(v, vlist);
while (v)
  {
  spatial_vertex_add(v->x, NULL, spatial);

#if DEBUG_SPATIAL_CELL_PLANE
P3VEC("v: ", v->x);
for (list=v->adj ; list ; list=g_slist_next(list))
  {
  printf(" + (%d)\n", GPOINTER_TO_INT(list->data));
  }
#endif

    vlist = g_slist_remove(vlist, v);
    v = vertex_common_facet(v, vlist);
    }
  }
g_slist_free(vlist);
}

/**************************************/
/* spatial point definition primitive */
/**************************************/
void spatial_point_add(struct core_pak *core, struct model_pak *data)
{
gint reset=0;
gchar *label;
gdouble plane[4], x1[3], x2[3], ortho1[3], ortho2[3];
GSList *list;
struct core_pak *core1, *core2;
struct spatial_pak *spatial;
static guint count=0;
static GSList *cores=NULL;

g_assert(core != NULL);
g_assert(data != NULL);

/* add the new atom */
cores = g_slist_append(cores, core);
select_add_core(core, data);
data->state++;

switch( data->mode)
  {
  case DEFINE_VECTOR:
    if (data->state==2)
      {
/* new vector */
      label = g_strdup_printf("vector_%d", count++);
      spatial = spatial_new(label, SPATIAL_VECTOR, 2, FALSE, data);
      g_free(label);

/* get the (lattice space) separation vector */
      core1 = g_slist_nth_data(cores, 0);
      g_assert(core1 != NULL);

      ARR3SET(x1, core1->x);
      ARR3SET(x2, core->x);
      ARR3SUB(x2, core1->x);

/* convert to real space & scale */
      vecmat(data->latmat, x2);
      normalize(x2, 3);
      VEC3MUL(x2, 0.7);
/* convert back to lattice space */
      vecmat(data->ilatmat, x2);
      ARR3ADD(x2, x1);

/* add the vertices */
      spatial_vertex_add(x2, NULL, spatial);
      spatial_vertex_add(x1, NULL, spatial);

/* drawing update */
      coords_compute(data);
      sysenv.refresh_dialog=TRUE;
      redraw_canvas(SINGLE);
      reset++;
      }
    break;

  case DEFINE_PLANE:
    if (data->state == 3)
      {

/* CURRENT */

if (data->periodic == 3)
  {
  spatial_cell_plane(cores, data);
  }
else
  {
/* new plane */
      label = g_strdup_printf("plane_%d", count++);
      spatial = spatial_new(label, SPATIAL_GENERIC, 4, FALSE, data);
      g_free(label);

      core1 = g_slist_nth_data(cores, 0);
      core2 = g_slist_nth_data(cores, 1);
      g_assert(core1 != NULL);
      g_assert(core2 != NULL);

/* compute plane normal */
      calc_norm(plane, core1->rx, core2->rx, core->rx);

/* vector from first to second point (NB: normalized in cartesian space) */
      ARR3SET(ortho1, core2->x);
      ARR3SUB(ortho1, core1->x);
      vecmat(data->latmat, ortho1);
      normalize(ortho1, 3);
      vecmat(data->ilatmat, ortho1);

/* get orthogonal vector */
      crossprod(ortho2, plane, ortho1);
      vecmat(data->latmat, ortho2);
      normalize(ortho2, 3);
      vecmat(data->ilatmat, ortho2);

/* define 4 points for the plane */
      ARR3SET(x1, core1->x);    
      ARR3ADD(x1, ortho1);
      ARR3ADD(x1, ortho2);
      spatial_vertex_add(x1, NULL, spatial);

      ARR3SET(x1, core1->x);    
      ARR3ADD(x1, ortho1);
      ARR3SUB(x1, ortho2);
      spatial_vertex_add(x1, NULL, spatial);

      ARR3SET(x1, core1->x);    
      ARR3SUB(x1, ortho1);
      ARR3SUB(x1, ortho2);
      spatial_vertex_add(x1, NULL, spatial);

      ARR3SET(x1, core1->x);    
      ARR3SUB(x1, ortho1);
      ARR3ADD(x1, ortho2);
      spatial_vertex_add(x1, NULL, spatial);
  }

/* drawing update */
      coords_init(REDO_COORDS, data);
      sysenv.refresh_dialog=TRUE;
      redraw_canvas(SINGLE);
      reset++;
      }
    break;
  }

/* completed an object */
if (reset)
  {
  for (list=cores ; list ; list=g_slist_next(list))
    {
    core1 = list->data;
    select_del_core(core1, data);
    }
  g_slist_free(cores);
  cores=NULL;
  data->state = 0;
  }
}

