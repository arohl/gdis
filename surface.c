/*
Copyright (C) 1993 by David H. Gay and Andrew L. Rohl 

dgay@ricx.royal-institution.ac.uk
andrew@ricx.royal-institution.ac.uk

Modified 2003 by Sean David Fleming

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
#include <string.h>
#include <stdlib.h>
#include <math.h>

#include "gdis.h"
#include "coords.h"
#include "model.h"
#include "edit.h"
#include "matrix.h"
#include "numeric.h"
#include "vector.h"
#include "space.h"
#include "surface.h"
#include "zone.h"
#include "interface.h"

typedef gchar boolean;

/**********************/
/* debugging routines */
/**********************/
void surface_print_planes(GSList *planes)
{
GSList *list;
struct plane_pak *plane;

for (list=planes ; list ; list=g_slist_next(list))
  {
  plane = list->data;

printf("[%2d %2d %2d]", plane->index[0], plane->index[1], plane->index[2]);
printf(" : %f , %f", plane->esurf[0], plane->eatt[0]);
printf("\n");
  }
}

/*********************/
/* duplicate a shift */
/*********************/
gpointer shift_dup(struct shift_pak *src)
{
struct shift_pak *dest;

dest = g_malloc(sizeof(struct shift_pak));
memcpy(dest, src, sizeof(struct shift_pak));
return(dest);
}

/********************************************/
/* allocate for new shift & init for safety */
/********************************************/
gpointer shift_new(gdouble shift)
{
struct shift_pak *sdata;

/* alloc */
sdata = g_malloc(sizeof(struct shift_pak));

sdata->locked = FALSE;
sdata->shift = shift;
sdata->dipole_computed = FALSE;
sdata->dipole = 99.9;
sdata->region[0] = 1;
sdata->region[1] = 1;
sdata->esurf[0] = 0.0;
sdata->esurf[1] = 0.0;
sdata->eatt[0] = 0.0;
sdata->eatt[1] = 0.0;
sdata->bbpa = 0.0;
sdata->gnorm = -1.0;
sdata->procfile = NULL;

return(sdata);
}

/*********************/
/* duplicate a plane */
/*********************/
gpointer plane_dup(struct plane_pak *src)
{
GSList *list;
struct plane_pak *dest;
struct shift_pak *shift;

dest = g_malloc(sizeof(struct plane_pak));

memcpy(dest, src, sizeof(struct plane_pak));
dest->shifts = NULL;
for (list=src->shifts ; list ; list=g_slist_next(list))
  {
  shift = shift_dup(list->data);
  dest->shifts = g_slist_prepend(dest->shifts, shift);
  }
return(dest);
}

/****************************************/
/* determines surface and depth vectors */
/****************************************/
/* TODO - call this to eliminate redundancy in generate_surface() */
#define DEBUG_GENSURF 0
void plane_compute_vectors(struct plane_pak *plane, struct model_pak *src)
{
vector normal;
vector lattice_vector[3];    /* lattice vectors */
vector t_mat[3];             /* transformation matrix */
vector work_lat[3];          /* transformed lattice vectors */
vector rec_work_lat[3];      /* reciprocal transformed lattice vectors */
vector s[3];                 /* space we are trying to fill */
vector tempvec;
vector *sv;
vector a, *v_a=NULL, *v_b=NULL;
boolean s2_found = FALSE;
gint c, i, j, h, k, l, flag;
gint GCD(gint, gint);
gint vector_compare(vector *, vector *);
gdouble inv_denom;
gdouble depth, depth_min;
gdouble gcd, cd, x;
gdouble tmat[9], norm[3], va[3];
gdouble wlat[9];
gdouble sign, vec[3], vec2[3], dn[3];
GSList *list1, *list2, *sv_list1, *sv_list2;

g_assert(src != NULL);

/* find surface normal from the miller index and lattice vectors*/
V_ZERO(normal);

/* setup */
h = plane->m[0];
k = plane->m[1];
l = plane->m[2];

/* acquire lattice vectors */
V_X(lattice_vector[0]) = src->latmat[0];
V_Y(lattice_vector[0]) = src->latmat[3];
V_Z(lattice_vector[0]) = src->latmat[6];
V_X(lattice_vector[1]) = src->latmat[1];
V_Y(lattice_vector[1]) = src->latmat[4];
V_Z(lattice_vector[1]) = src->latmat[7];
V_X(lattice_vector[2]) = src->latmat[2];
V_Y(lattice_vector[2]) = src->latmat[5];
V_Z(lattice_vector[2]) = src->latmat[8];

V_CROSS(tempvec, lattice_vector[1], lattice_vector[2]);
V_QADD(normal, normal, +h*, tempvec);
V_CROSS(tempvec, lattice_vector[2], lattice_vector[0]);
V_QADD(normal, normal, +k*, tempvec);
V_CROSS(tempvec, lattice_vector[0], lattice_vector[1]);
V_QADD(normal, normal, +l*, tempvec);

sv_list1 = NULL;
 
c = (float) GCD(h, k);
cd = c ? c : 1.0;
V_2_ASSIGN(a, = k/cd * ,lattice_vector[0], - h/cd *, lattice_vector[1]);
if (V_MAGSQ(a) > EPSILON) 
  {
  sv = g_malloc(sizeof(vector));
  V_EQUATE(*sv, a);
  sv_list1 = g_slist_prepend(sv_list1, sv);
  }
  
c = (float) GCD(h, l);
cd = c ? c : 1.0;
V_2_ASSIGN(a, = l/cd * ,lattice_vector[0], - h/cd *, lattice_vector[2]);
if (V_MAGSQ(a) > EPSILON) 
  {
  sv = g_malloc(sizeof(vector));
  V_EQUATE(*sv, a);
  sv_list1 = g_slist_prepend(sv_list1, sv);
  }
  
c = (float) GCD(k, l);
cd = c ? c : 1.0;
V_2_ASSIGN(a, = l/cd * ,lattice_vector[1], - k/cd *, lattice_vector[2]);
if (V_MAGSQ(a) > EPSILON) 
  {
  sv = g_malloc(sizeof(vector));
  V_EQUATE(*sv, a);
  sv_list1 = g_slist_prepend(sv_list1, sv);
  }

sv_list2 = NULL;
list1 = sv_list1;
for (i=0 ; i<g_slist_length(sv_list1)-1 ; i++)
  {
  v_a = (vector *) list1->data;
  list1 = list2 = g_slist_next(list1);
  for (j=i+1 ; j<g_slist_length(sv_list1) ; j++)
    {
    v_b = (vector *) list2->data;
    list2 = g_slist_next(list2);


    V_2_ASSIGN(a, = , *v_a, + , *v_b);
    if (V_MAGSQ(a) > EPSILON)
      {
      sv = g_malloc(sizeof(vector));
      V_EQUATE(*sv, a);
      sv_list2 = g_slist_prepend(sv_list2, sv);
      }
    V_2_ASSIGN(a, = , *v_a, - , *v_b);
    if (V_MAGSQ(a) > EPSILON)
      {
      sv = g_malloc(sizeof(vector));
      V_EQUATE(*sv, a);
      sv_list2 = g_slist_prepend(sv_list2, sv);
      }
    }
  }
sv_list1 = g_slist_concat(sv_list1, sv_list2);
sv_list1 = g_slist_sort(sv_list1, (gpointer) vector_compare);

/* set v_a to first (shortest now sorted) vector */
list1 = sv_list1;
v_a = (vector *) list1->data;
list1 = g_slist_next(list1);
/* loop over remaining vectors - find shortest w/ orthogonal component */
while (list1)
  {
  v_b = (vector *) list1->data;
  
  V_CROSS(a, *v_a, *v_b);
  x = V_MAGSQ(a);
  if (x > EPSILON) 
    {
    s2_found = TRUE;
    break;
    }
  list1 = g_slist_next(list1);
  }

/* checks */
if (!s2_found) 
  {
  printf("Failed to find surface vectors.\n");
  return;
  }

g_assert(v_a != NULL);
g_assert(v_b != NULL);

/* calculate transformation matrix */
V_SCALER(t_mat[2], (1.0/V_MAG(normal)), normal);
V_SCALER(t_mat[0], (1.0/V_MAG(*v_a)), *v_a);
V_CROSS(t_mat[1], t_mat[2], t_mat[0]);

ARR3SET(norm, normal.element);
ARR3SET(va, v_a->element);
normalize(norm, 3);
normalize(va, 3);
ARR3SET(&tmat[6], norm);
ARR3SET(&tmat[0], va);
crossprod(&tmat[3], norm, va);

/* calculate transformed lattice vectors */
for (j = 0; j < 3; j++)
  for (i = 0; i < 3; i++)
    V_E(work_lat[j], i) = V_DOT(t_mat[i], lattice_vector[j]);

/* calculate reciprocal transformed lattice vectors */
V_CROSS(tempvec, work_lat[1], work_lat[2]);
inv_denom = 1.0 / V_DOT(tempvec, work_lat[0]);
V_CROSS(tempvec, work_lat[1], work_lat[2]);
V_SCALER(rec_work_lat[0], inv_denom, tempvec);
V_CROSS(tempvec, work_lat[2], work_lat[0]);
V_SCALER(rec_work_lat[1], inv_denom, tempvec);
V_CROSS(tempvec, work_lat[0], work_lat[1]);
V_SCALER(rec_work_lat[2], inv_denom, tempvec);

/* calculate transformed surface vectors */
for (i = 0; i < 3; i++)
  V_E(s[0], i) = V_DOT(t_mat[i], *v_a);
for (i = 0; i < 3; i++)
  V_E(s[1], i) = V_DOT(t_mat[i], *v_b);

/* return transformed surface vectors data in dest->latmat */
  /* NB: latmat[0..8] is a 3x3 matrix so, */
  /*
    | a b 0 |
    | c d 0 |
    | 0 0 0 |

    goes into latmat rows first ie a,b,0,c,d,0,0,0,0
  */
  
/* set s[2] = to depths */
/* depth = 1.0e10; */

/* surface specification */
gcd = GCD(GCD(h, k), GCD(k, l));
VEC3SET(vec, h, k, l);
vecmat(src->rlatmat, vec);
 
/* calculate the Dhkl */
/* NB: we want the correct dspacing wrt hkl, ie DON'T remove the GCD */
/*
dest->surface.dspacing = 1.0/VEC3MAG(vec);
*/

/* actual cut depth */
/* NB: done wrt the FULL depth (ie remove the gcd) */
depth = gcd/VEC3MAG(vec);

/*
dest->surface.depth = depth;
*/
/* round off the shift */
/*
shift = decimal_round(dest->surface.shift, 4);
*/

/* store the work lattice */
for (i=0 ; i<3 ; i++)
  {
  wlat[3*i+0] = V_E(work_lat[0], i);
  wlat[3*i+1] = V_E(work_lat[1], i);
  wlat[3*i+2] = V_E(work_lat[2], i);
  }

/* transfer surface vectors */
VEC3SET(&plane->lattice[0], V_X(s[0]), V_X(s[1]), 0.0);
VEC3SET(&plane->lattice[3], V_Y(s[0]), V_Y(s[1]), 0.0);
VEC3SET(&plane->lattice[6], 0.0, 0.0, 1.0);

/* construct the depth vector */
/* compute the surface normal */
VEC3SET(vec, V_X(s[0]), V_Y(s[0]), 0.0);
VEC3SET(vec2, V_X(s[1]), V_Y(s[1]), 0.0);
crossprod(dn, vec, vec2);

/* search for a vector with the smallest */
/* orthogonal component to the surface vectors */
flag=0;
depth_min=0;
for (i=0 ; i<3 ; i++)
  {
/* ignore zero indices, which yield 0 orthogonal component */
  if (plane->m[i])
    {
/* dotproduct of vector with normal */
    ARR3SET(vec, dn);
    vec[0] *= wlat[i];
    vec[1] *= wlat[i+3];
    vec[2] *= wlat[i+6];
    depth = fabs(vec[0] + vec[1] + vec[2]);

   if (flag)
     {
     if (depth > depth_min)
       continue;
     }
   else
     flag++;

    depth_min = depth;

/* ensure the depth vector is pointing in the +ve z direction */
    if (wlat[i+6] < 0.0)
      sign = -1.0;
    else
      sign = 1.0;

    plane->lattice[2] = sign*wlat[i];
    plane->lattice[5] = sign*wlat[i+3];
    plane->lattice[8] = sign*wlat[i+6];
    }
  }
}

/********************************************/
/* allocate for new plane & init for safety */
/********************************************/
#define DEBUG_CREATE_PLANE 0
gpointer plane_new(gdouble *m, struct model_pak *model)
{
gint h, k, l, n;
gdouble len, vec[3];
struct plane_pak *plane=NULL;

/* checks */
g_return_val_if_fail(model != NULL, NULL);

plane = g_malloc(sizeof(struct plane_pak));

#if DEBUG_CREATE_PLANE
printf("creating: %d %d %d -> ", (gint) m[0], (gint) m[1], (gint) m[2]);
#endif

/* save orig */
h = m[0];
k = m[1];
l = m[2];
n = 2;

/* absence testing function */
if (!model->surface.ignore_symmetry)
  {
  while (surf_sysabs(model, h, k, l) && n<100)
    {
    h = n*m[0];
    k = n*m[1];
    l = n*m[2];
    n++;
    }
  }

if (n == 100)
  {
  printf("WARNING: screwed up systematic absence check.\n");
  h /= 99;
  k /= 99;
  l /= 99;
  n = 2;
  }

#if DEBUG_CREATE_PLANE
printf("%d %d %d  (x %d)\n", h, k, l, n-1);
#endif

/* init */
VEC3SET(plane->m, (gdouble) h, (gdouble) k, (gdouble) l);
VEC3SET(plane->norm, (gdouble) h, (gdouble) k, (gdouble) l);
VEC3SET(plane->index, h, k, l);
plane->shifts = NULL;
plane->vertices = NULL;

/* calc Dhkl */
VEC3SET(vec, h, k, l);
vecmat(model->rlatmat, vec);
plane->dhkl = 1.0/VEC3MAG(vec);

/* morphology prediction */
len = VEC3MAG(vec);
if (len > FRACTION_TOLERANCE)
  {
  VEC3MUL(vec, 1.0/len);
  VEC3MUL(vec, -1.0);
  ARR3SET(plane->x, vec);
  }
else
  {
  gui_text_show(WARNING, "plane of zero length created.\n");
  }

plane->esurf_shift = 0.0;
plane->eatt_shift = 0.0;
plane->bbpa_shift = 0.0;
plane->esurf[0] = 0.0;
plane->esurf[1] = 0.0;
plane->eatt[0] = 0.0;
plane->eatt[1] = 0.0;
plane->bbpa = 0.0;
plane->area = 0.0;
plane->f[0] = 0.0;
plane->f[1] = 0.0;
/* calc? */
plane->multiplicity = 1;
plane->present = TRUE;
plane->visible = FALSE;
plane->primary = TRUE;
plane->parent = NULL;
VEC3SET(plane->rx, 0.0, 0.0, 0.0);

/* NEW */
plane_compute_vectors(plane, model);

return(plane);
}

/************************************************************/
/* determine if a particular plane has already been created */
/************************************************************/
gpointer plane_find(gdouble *hkl, struct model_pak *model)
{
GSList *list;
struct plane_pak *comp, *plane;

/* get corrected (sys absent!) miller index */
comp = plane_new(hkl, model);
g_return_val_if_fail(comp != NULL, NULL);

for (list=model->planes ; list ; list=g_slist_next(list))
  {
  plane = list->data;
  if (facet_equiv(model, comp->index, plane->index))
    return(plane);
  }

g_free(comp);
return(NULL);
}

/***********************************/
/* generate symmetry related faces */
/***********************************/
/* FIXME - handle duplication */
#define DEBUG_SURF_SYMMETRY_GENERATE 0
void surf_symmetry_generate(struct model_pak *model)
{
gint i;
gint *sym_hkl;
gdouble m[3];
GSList *list, *new_planes, *sym_equiv, *sym_list;
struct plane_pak *plane=NULL, *new_plane=NULL;

#if DEBUG_SURF_SYMMETRY_GENERATE
printf("======\n");
printf("BEFORE\n");
printf("======\n");
surface_print_planes(model->planes);
#endif

/* compute symmetry related faces */
new_planes = NULL;
for (list=model->planes ; list ; list=g_slist_next(list))
  {
  plane = list->data;
  sym_equiv = get_facet_equiv(model, plane->index);

#if DEBUG_SURF_SYMMETRY_GENERATE
printf(" *** {%d %d %d}\n", plane->index[0], plane->index[1], plane->index[2]);
#endif

/* check each non-primary face for extinctions */
  i=0;
  sym_list = g_slist_next(sym_equiv);
  while (sym_list)
    {
/* get the (normalized) hkl's */
    sym_hkl = sym_list->data;
    ARR3SET(m, sym_hkl);

    new_plane = plane_new(m, model);
    if (plane)
      {
      new_plane->primary = FALSE;
      new_plane->parent = plane;
      new_plane->esurf[0] = plane->esurf[0];
      new_plane->esurf[1] = plane->esurf[1];
      new_plane->eatt[0] = plane->eatt[0];
      new_plane->eatt[1] = plane->eatt[1];
      new_plane->bbpa = plane->bbpa;
      new_plane->area = plane->area;

#if DEBUG_SURF_SYMMETRY_GENERATE
printf("   > (%d %d %d)\n",
       new_plane->index[0], new_plane->index[1], new_plane->index[2]);
#endif

      new_planes = g_slist_prepend(new_planes, new_plane);
      }
    sym_list = g_slist_next(sym_list);
    i++;
    }
  free_slist(sym_equiv);
  }
/* use the new list with the symmetry adjusted planes */
model->planes = g_slist_concat(model->planes, new_planes);
model->num_planes = g_slist_length(model->planes);

#if DEBUG_SURF_SYMMETRY_GENERATE
printf("======\n");
printf("ADDING\n");
printf("======\n");
surface_print_planes(new_planes);
#endif

}

/***********************/
/* free a single shift */
/***********************/
void shift_free(gpointer data)
{
struct shift_pak *shift = data;

g_assert(shift != NULL);

/* FIXME - this may create a memory leak, but it'll be */
/* the users fault */
if (shift->locked)
  {
  printf("ERROR: shift is locked.\n");
  return;
  }
g_free(shift);
}

/*************************************/
/* free the data in a list of shifts */
/*************************************/
void shift_data_free(GSList *shifts)
{
GSList *list;
struct shift_pak *shift;

/* free a list of shift pak structures */
list = shifts;
while (list)
  {
  shift = list->data;
  list = g_slist_next(list);

  shift_free(shift);
  }
}

/***********************/
/* free a single plane */
/***********************/
void plane_free(gpointer data)
{
struct plane_pak *plane = data;

shift_data_free(plane->shifts);
g_slist_free(plane->shifts);

/* FIXME - is this good enough? */
g_slist_free(plane->vertices);
g_free(plane);
}

/*************************************/
/* free the data in a list of planes */
/*************************************/
void plane_data_free(GSList *planes)
{
GSList *list;
struct plane_pak *plane;

/* free a list of plane pak structures */
list = planes;
while (list)
  {
  plane = list->data;
  list = g_slist_next(list);

  plane_free(plane);
  }
}

/* DEBUG */

void print_shell_offset(gint flag, struct model_pak *model)
{
gdouble x[3];
GSList *list;
struct core_pak *core;
struct shel_pak *shell;

printf("-----------------------------------------------\n");
printf("model: %p\n", model);
P3MAT("latmat:", model->latmat);
printf("-----------------------------------------------\n");

for (list=model->cores ; list ; list=g_slist_next(list))
  {
  core = list->data;

  if (core->shell)
    {
    shell = core->shell;

    ARR3SET(x, shell->x);
    ARR3SUB(x, core->x);
    if (flag)
      vecmat(model->latmat, x);

printf("[r = %f] ", VEC3MAG(x));

    P3VEC(" + ", x);
    }
  }
}

/**********************************************/
/* refresh plane energy data for all surfaces */
/**********************************************/
void plane_refresh_all(struct model_pak *model)
{
/* refresh primary planes */
/* refresh symmetry related planes */
}

/****************************************/
/* uses marvin code to create a surface */
/****************************************/
#define DEBUG_GENSURF 0
gint generate_surface(struct model_pak *src, struct model_pak *dest)
{
vector normal;
vector lattice_vector[3];    /* lattice vectors */
vector t_mat[3];             /* transformation matrix */
vector work_lat[3];          /* transformed lattice vectors */
vector rec_work_lat[3];      /* reciprocal transformed lattice vectors */
vector s[3];                 /* space we are trying to fill */
vector rec_s[2];             /* reciprocal of s[0] and s[1] */
vector tempvec;
vector *sv;
vector a, *v_a=NULL, *v_b=NULL;
boolean s2_found = FALSE;
gint c, i, j, h, k, l, flag;
gint ia, ib, ic;
gint amin, amax, bmin, bmax, cmin, cmax;
gint ob, sb;
gint GCD(gint, gint);
gint vector_compare(vector *, vector *);
gdouble inv_denom;
gdouble shift, depth, depth_1, depth_2, depth_min;
gdouble gcd, cd, x;
gdouble z1_max, z1_min, z2_max, z2_min;
gdouble tempfloat;
gdouble tmat[9], norm[3], va[3];
gdouble wlat[9], temp[9], lpm[9];
gdouble sign, vec[3], vec2[3], dn[3];
GSList *list1, *list2, *sv_list1, *sv_list2;
GSList *list, *clist, *slist, *mlist;
struct core_pak *core;
struct shel_pak *shel;
struct mol_pak *mol;

/* checks */
g_assert(src != NULL);
g_assert(dest != NULL);

/* find surface normal from the miller index and lattice vectors*/
V_ZERO(normal);

/* setup */
h = dest->surface.miller[0];
k = dest->surface.miller[1];
l = dest->surface.miller[2];

/* acquire lattice vectors */
V_X(lattice_vector[0]) = src->latmat[0];
V_Y(lattice_vector[0]) = src->latmat[3];
V_Z(lattice_vector[0]) = src->latmat[6];
V_X(lattice_vector[1]) = src->latmat[1];
V_Y(lattice_vector[1]) = src->latmat[4];
V_Z(lattice_vector[1]) = src->latmat[7];
V_X(lattice_vector[2]) = src->latmat[2];
V_Y(lattice_vector[2]) = src->latmat[5];
V_Z(lattice_vector[2]) = src->latmat[8];

V_CROSS(tempvec, lattice_vector[1], lattice_vector[2]);
V_QADD(normal, normal, +h*, tempvec);
V_CROSS(tempvec, lattice_vector[2], lattice_vector[0]);
V_QADD(normal, normal, +k*, tempvec);
V_CROSS(tempvec, lattice_vector[0], lattice_vector[1]);
V_QADD(normal, normal, +l*, tempvec);

sv_list1 = NULL;
 
c = (float) GCD(h, k);
cd = c ? c : 1.0;
V_2_ASSIGN(a, = k/cd * ,lattice_vector[0], - h/cd *, lattice_vector[1]);
if (V_MAGSQ(a) > EPSILON) 
  {
  sv = g_malloc(sizeof(vector));
  V_EQUATE(*sv, a);
  sv_list1 = g_slist_prepend(sv_list1, sv);
  }
  
c = (float) GCD(h, l);
cd = c ? c : 1.0;
V_2_ASSIGN(a, = l/cd * ,lattice_vector[0], - h/cd *, lattice_vector[2]);
if (V_MAGSQ(a) > EPSILON) 
  {
  sv = g_malloc(sizeof(vector));
  V_EQUATE(*sv, a);
  sv_list1 = g_slist_prepend(sv_list1, sv);
  }
  
c = (float) GCD(k, l);
cd = c ? c : 1.0;
V_2_ASSIGN(a, = l/cd * ,lattice_vector[1], - k/cd *, lattice_vector[2]);
if (V_MAGSQ(a) > EPSILON) 
  {
  sv = g_malloc(sizeof(vector));
  V_EQUATE(*sv, a);
  sv_list1 = g_slist_prepend(sv_list1, sv);
  }

sv_list2 = NULL;
list1 = sv_list1;
for (i=0 ; i<g_slist_length(sv_list1)-1 ; i++)
  {
  v_a = (vector *) list1->data;
  list1 = list2 = g_slist_next(list1);
  for (j=i+1 ; j<g_slist_length(sv_list1) ; j++)
    {
    v_b = (vector *) list2->data;
    list2 = g_slist_next(list2);


    V_2_ASSIGN(a, = , *v_a, + , *v_b);
    if (V_MAGSQ(a) > EPSILON)
      {
      sv = g_malloc(sizeof(vector));
      V_EQUATE(*sv, a);
      sv_list2 = g_slist_prepend(sv_list2, sv);
      }
    V_2_ASSIGN(a, = , *v_a, - , *v_b);
    if (V_MAGSQ(a) > EPSILON)
      {
      sv = g_malloc(sizeof(vector));
      V_EQUATE(*sv, a);
      sv_list2 = g_slist_prepend(sv_list2, sv);
      }
    }
  }
sv_list1 = g_slist_concat(sv_list1, sv_list2);
sv_list1 = g_slist_sort(sv_list1, (gpointer) vector_compare);

/* set v_a to first (shortest now sorted) vector */
list1 = sv_list1;
v_a = (vector *) list1->data;
list1 = g_slist_next(list1);
/* loop over remaining vectors - find shortest w/ orthogonal component */
while (list1)
  {
  v_b = (vector *) list1->data;
  
  V_CROSS(a, *v_a, *v_b);
  x = V_MAGSQ(a);
  if (x > EPSILON) 
    {
    s2_found = TRUE;
    break;
    }
  list1 = g_slist_next(list1);
  }

/* checks */
if (!s2_found) 
  {
  gui_text_show(ERROR, "Failed to find surface vectors.\n");
  return(1);
  }
g_assert(v_a != NULL);
g_assert(v_b != NULL);

/* calculate transformation matrix */
V_SCALER(t_mat[2], (1.0/V_MAG(normal)), normal);
V_SCALER(t_mat[0], (1.0/V_MAG(*v_a)), *v_a);
V_CROSS(t_mat[1], t_mat[2], t_mat[0]);

ARR3SET(norm, normal.element);
ARR3SET(va, v_a->element);
normalize(norm, 3);
normalize(va, 3);
ARR3SET(&tmat[6], norm);
ARR3SET(&tmat[0], va);
crossprod(&tmat[3], norm, va);

/* calculate transformed lattice vectors */
for (j = 0; j < 3; j++)
  for (i = 0; i < 3; i++)
    V_E(work_lat[j], i) = V_DOT(t_mat[i], lattice_vector[j]);

/* calculate reciprocal transformed lattice vectors */
V_CROSS(tempvec, work_lat[1], work_lat[2]);
inv_denom = 1.0 / V_DOT(tempvec, work_lat[0]);
V_CROSS(tempvec, work_lat[1], work_lat[2]);
V_SCALER(rec_work_lat[0], inv_denom, tempvec);
V_CROSS(tempvec, work_lat[2], work_lat[0]);
V_SCALER(rec_work_lat[1], inv_denom, tempvec);
V_CROSS(tempvec, work_lat[0], work_lat[1]);
V_SCALER(rec_work_lat[2], inv_denom, tempvec);

/* calculate transformed surface vectors */
for (i = 0; i < 3; i++)
  V_E(s[0], i) = V_DOT(t_mat[i], *v_a);
for (i = 0; i < 3; i++)
  V_E(s[1], i) = V_DOT(t_mat[i], *v_b);

/* return transformed surface vectors data in dest->latmat */
  /* NB: latmat[0..8] is a 3x3 matrix so, */
  /*
    | a b 0 |
    | c d 0 |
    | 0 0 0 |

    goes into latmat rows first ie a,b,0,c,d,0,0,0,0
  */
  
/* set s[2] = to depths */
/* depth = 1.0e10; */

/* surface specification */
gcd = GCD(GCD(h, k), GCD(k, l));
VEC3SET(vec, h, k, l);
vecmat(src->rlatmat, vec);
 
/* calculate the Dhkl */
/* NB: we want the correct dspacing wrt hkl, ie DON'T remove the GCD */
dest->surface.dspacing = 1.0/VEC3MAG(vec);

/* actual cut depth */
/* NB: done wrt the FULL depth (ie remove the gcd) */
depth = gcd/VEC3MAG(vec);
dest->surface.depth = depth;

/* round off the shift */
shift = decimal_round(dest->surface.shift, 4);

#if DEBUG_GENSURF
printf("\n----------------------------\n");
printf("    hkl: %d %d %d\n", h, k, l);
printf("    gcd: %f \n", gcd);
printf("  shift: %.20f\n", shift);
printf("regions: %d %d\n", (gint) dest->surface.region[0], (gint) dest->surface.region[1]);
printf("   Dhkl: %f\n", dest->surface.dspacing);
printf("  depth: %f\n", depth);
printf("----------------------------\n");
#endif

/* NB: depth_1 is region[0], but depth_2 is region[0]+region[1] */
depth_1 = dest->surface.region[0] * depth;
depth_2 = depth_1 + dest->surface.region[1] * depth;

V_ZERO(s[2]);
V_Z(s[2]) = -(depth_2);

/* transfer surface vectors */
VEC3SET(&dest->latmat[0], V_X(s[0]), V_X(s[1]), 0.0);
VEC3SET(&dest->latmat[3], V_Y(s[0]), V_Y(s[1]), 0.0);
VEC3SET(&dest->latmat[6], 0.0, 0.0, 1.0);

/* calculate reciprocal of two surface vectors */
inv_denom = 1.0 / (V_X(s[0])*V_Y(s[1]) - V_Y(s[0])*V_X(s[1]));
V_X(rec_s[0]) =  V_Y(s[1])*inv_denom;
V_Y(rec_s[0]) = -V_X(s[1])*inv_denom;
V_Z(rec_s[0]) =  0.0;
V_X(rec_s[1]) = -V_Y(s[0])*inv_denom;
V_Y(rec_s[1]) =  V_X(s[0])*inv_denom;
V_Z(rec_s[1]) =  0.0;

z2_min = z1_min = 0.00;
z2_max = V_Z(s[2]);
z1_max = depth_1 * (z2_max < 0 ? -1.0: +1.0);

if (z2_max < z2_min) 
  {
  tempfloat = z2_max;
  z2_max = z2_min;
  z2_min = tempfloat;
  }
if (z1_max < z1_min) 
  {
  tempfloat = z1_max;
  z1_max = z1_min;
  z1_min = tempfloat;
  }

/* display name */
g_free(dest->basename);
dest->basename = g_strdup_printf("%s_%-1d%-1d%-1d_%6.4f",
                 g_strstrip(src->basename), h, k, l, shift);

/* store the work lattice */
for (i=0 ; i<3 ; i++)
  {
  wlat[3*i+0] = V_E(work_lat[0], i);
  wlat[3*i+1] = V_E(work_lat[1], i);
  wlat[3*i+2] = V_E(work_lat[2], i);
  }

#if DEBUG_GENSURF
P3MAT(" src lat: ", src->latmat);
#endif

/* construct the depth vector */
/* compute the surface normal */
VEC3SET(vec, dest->latmat[0], dest->latmat[3], dest->latmat[6]);
VEC3SET(vec2, dest->latmat[1], dest->latmat[4], dest->latmat[7]);
crossprod(dn, vec, vec2);

/* search for a vector with the smallest */
/* orthogonal component to the surface vectors */
flag=0;
depth_min=0;
for (i=0 ; i<3 ; i++)
  {
/* ignore zero indices, which yield 0 orthogonal component */
  if (dest->surface.miller[i])
    {
/* dotproduct of vector with normal */
    ARR3SET(vec, dn);
    vec[0] *= wlat[i];
    vec[1] *= wlat[i+3];
    vec[2] *= wlat[i+6];
    depth = fabs(vec[0] + vec[1] + vec[2]);

#if DEBUG_GENSURF
printf("[i=%d, depth=%f]", i, depth);
#endif

   if (flag)
     {
     if (depth > depth_min)
       continue;
     }
   else
     flag++;

    depth_min = depth;

/* ensure the depth vector is pointing in the +ve z direction */
    if (wlat[i+6] < 0.0)
      sign = -1.0;
    else
      sign = 1.0;

    dest->surface.depth_vec[0] = sign*wlat[i];
    dest->surface.depth_vec[1] = sign*wlat[i+3];
    dest->surface.depth_vec[2] = sign*wlat[i+6];
    dest->latmat[2] = sign*wlat[i];
    dest->latmat[5] = sign*wlat[i+3];
    dest->latmat[8] = sign*wlat[i+6];
    }
  }
#if DEBUG_GENSURF
printf(" : Minimum = %f\n", depth_min);
P3VEC("depth vec:", dest->surface.depth_vec);
#endif

g_assert(flag != 0);

/* generate lattice matrix inverse */
memcpy(dest->ilatmat, dest->latmat, 9*sizeof(gdouble));
matrix_invert(dest->ilatmat);

/* locate the position of the new cell's repeat vectors */
/* in terms of the lattice space of the source model */
/* this gives us the lattice point matrix (lpm) which */
/* is used to determine how many source unit cells */
/* will be required to fill out the new surface cell */
memcpy(lpm, dest->latmat, 9*sizeof(gdouble));
memcpy(temp, tmat, 9*sizeof(gdouble));
matrix_invert(temp);
matmat(temp, lpm);
matmat(src->ilatmat, lpm);

#if DEBUG_GENSURF
P3MAT("tmat: ", tmat);
P3MAT("dest latmat: ", dest->latmat);
P3MAT("dest ilatmat: ", dest->ilatmat);
P3MAT("lpm: ", lpm);
#endif

/* setup repeats required to fill the new cell */
amin=bmin=cmin=0;
amax=bmax=cmax=0;

/* required number of a cells */
for (i=0 ; i<3 ; i++)
  {
  ia = rint(lpm[i]);
  if (ia < amin)
    amin = ia;
  if (ia > amax)
    amax = ia;
  }
/* required number of b cells */
for (i=3 ; i<6 ; i++)
  {
  ib = rint(lpm[i]);
  if (ib < bmin)
    bmin = ib;
  if (ib > bmax)
    bmax = ib;
  }
/* required number of c cells */
for (i=6 ; i<9 ; i++)
  {
  ic = rint(lpm[i]);

  if (ic < cmin)
    cmin = ic;
  if (ic > cmax)
    cmax = ic;
  }

/* one cell buffer */
amax+=2;
bmax+=2;
cmax+=2;
amin--;
bmin--;
cmin--;

#if DEBUG_GENSURF
printf("a: %d - %d\n", amin, amax);
printf("b: %d - %d\n", bmin, bmax);
printf("c: %d - %d\n", cmin, cmax);
#endif

/* compute transformation pipeline */
memcpy(temp, src->latmat, 9*sizeof(gdouble));
matmat(tmat, temp);
matmat(dest->ilatmat, temp);

/* create a single transformed unit cell */
/* full loop required to cover one cell in the transformed coordinates */
for (ic=cmin ; ic<cmax ; ic++)
  {
  for (ib=bmin ; ib<bmax ; ib++)
    {
    for (ia=amin ; ia<amax ; ia++)
      {
/* current source cell translation */
      VEC3SET(vec, ia, ib, ic);

/* loop over all molecules */
      for (mlist=src->moles ; mlist ; mlist=g_slist_next(mlist))
        {
        mol = mlist->data;

/* transform the centroid */
        ARR3SET(va, mol->centroid);
        ARR3ADD(va, vec);
        vecmat(temp, va);
        va[2] += shift;

/* ignore molecules outside the transformed cell */
/* NB: the range must be exactly the smallest numerical value below 1.0 [G_MINDOUBLE, 1.0] */
/* NB: the value of EPS adds a bias to promote molecules to the top of the surface */
#define EPS 0.001
        if (va[0] < G_MINDOUBLE+EPS || va[0] > 1.0+EPS)
          continue;
        if (va[1] < G_MINDOUBLE+EPS || va[1] > 1.0+EPS)
          continue;
        if (va[2] < G_MINDOUBLE+EPS || va[2] > 1.0+EPS)
          continue;

/* loop over all atoms in molecule */
        for (clist=mol->cores ; clist ; clist=g_slist_next(clist))
          {
/* dup_core() will duplicate an attached shell as well */
          core = dup_core(clist->data);
          dest->cores = g_slist_prepend(dest->cores, core);

/* transform */
          ARR3ADD(core->x, vec);
          vecmat(temp, core->x);
          core->x[2] += shift;

/* init flags */
          core->primary = TRUE;
          core->primary_core = NULL;
          core->orig = TRUE;

/* shell */
          if (core->shell)
            {
            shel = core->shell;
            dest->shels = g_slist_prepend(dest->shels, shel);

            ARR3ADD(shel->x, vec);
            vecmat(temp, shel->x);
            shel->x[2] += shift;

/* init flags */
            shel->primary = TRUE;
            shel->primary_shell = NULL;
            shel->orig = TRUE;
            }
          }
        }
      }
    }
  }

dest->periodic = 3;
dest->fractional = FALSE;
dest->axes_type = CARTESIAN;
dest->construct_pbc = TRUE;

/* NB: can't leave the (expensive) delete_duplicate_cores() to the */
/* final model_prep() as that treats the model as 2D periodic */
/* and won't clip the top & bottom of the surface*/
matrix_lattice_init(dest);
zone_init(dest);
delete_duplicate_cores(dest);
shell_make_links(dest);

/* terminate here for debugging */
#define DEBUG_SURFACE_CELL 0
#if DEBUG_SURFACE_CELL
dest->sginfo.spacenum = 1;
dest->fractional = TRUE;
dest->axes_type = CARTESIAN;
dest->construct_pbc = TRUE;
dest->gulp.method = CONV;
model_prep(dest);
return(0);
#endif

/* build surface core and shell lists from transformed cell */
clist = slist = NULL;

/* region 1 cores and shells */
amax = dest->surface.region[0] + 1;
for (i=1 ; i<amax ; i++)
  {
  for (list=dest->cores ; list ; list=g_slist_next(list))
    {
    core = dup_core(list->data);
    clist = g_slist_prepend(clist, core);

/*
printf("[%s] c ", core->label);
P3VEC(" : ", core->x);
*/

    core->x[2] -= (gdouble) i;
    vecmat(dest->latmat, core->x);

    core->region = REGION1A;

    if (core->shell)
      {
      shel = core->shell;
      slist = g_slist_prepend(slist, shel);

/* CURRENT - alternate shell coordinate init */
/* more tolerant of pbc offset core-shell links */
      shel->x[2] -= (gdouble) i;
      vecmat(dest->latmat, shel->x);

/*
printf("[%s] s ", shel->label);
P3VEC(" : ", shel->x);
*/

      shel->region = REGION1A;
      }
    }
  }
/* region 2 cores and shells */
bmax = amax + dest->surface.region[1];
for (i=amax ; i<bmax ; i++)
  {
  for (list=dest->cores ; list ; list=g_slist_next(list))
    {
    core = dup_core(list->data);
    clist = g_slist_prepend(clist, core);
    core->x[2] -= (gdouble) i;
    vecmat(dest->latmat, core->x);
    core->region = REGION2A;
    if (core->shell)
      {
      shel = core->shell;
      slist = g_slist_prepend(slist, shel);
      shel->x[2] -= (gdouble) i;
      vecmat(dest->latmat, shel->x);
      shel->region = REGION2A;
      }
    }
  }

/* replace old structure (transformed cell) with the new surface */
free_core_list(dest);
dest->cores = g_slist_reverse(clist);
dest->shels = g_slist_reverse(slist);

/* adjust depth translation vector to account for region sizes */
dest->latmat[2] *= dest->surface.region[0]+dest->surface.region[1];
dest->latmat[5] *= dest->surface.region[0]+dest->surface.region[1];
dest->latmat[8] *= dest->surface.region[0]+dest->surface.region[1];

/* surface cell init */
dest->periodic = 2;
dest->sginfo.spacenum = 1;
dest->fractional = FALSE;
dest->axes_type = CARTESIAN;
dest->construct_pbc = TRUE;
dest->surface.keep_atom_order = src->surface.keep_atom_order;
dest->gulp.method = CONV;

/* surface cell type */
if (src->surface.true_cell)
  dest->surface.true_cell = TRUE;
else
  {
  dest->surface.true_cell = FALSE;

/* want to force c to be // to z */
/* NB: vectors are in columns */
  dest->latmat[2] = 0.0;
  dest->latmat[5] = 0.0;
  }
free_slist(sv_list1);

/* CURRENT - does this screw up the true_cell crap? */
/*
for (list=dest->cores ; list ; list=g_slist_next(list))
  {
  core = list->data;

  ARR3SET(va, core->x);
  vecmat(dest->latmat, va);

  core->x[2] = va[2];
  }
dest->latmat[2] = 0.0;
dest->latmat[5] = 0.0;
dest->latmat[8] = 1.0;
*/

#if DEBUG_GENSURF
printf("gensurf cores: %d\n", g_slist_length(dest->cores));
printf("gensurf shels: %d\n", g_slist_length(dest->shels));
#endif

/* build mode follows source structure */
/*
dest->build_molecules = src->build_molecules;
*/

/* last init */
model_prep(dest);
calc_emp(dest);

/* FIXME - bugged for unbroken bond chain materials */
/* mark the growth slice */
depth_min = dest->surface.region[0]+dest->surface.region[1];
depth_min *= gcd;
depth_min = -1.0/depth_min;
for (list1=dest->moles ; list1 ; list1=g_slist_next(list1))
  {
  mol = list1->data;
  if (mol->centroid[2] > depth_min)
    {
    for (list2=mol->cores ; list2 ; list2=g_slist_next(list2))
      {
      core = list2->data;
      core->growth = TRUE;
      }
    }
  }

/* NEW - bonding analysis */
ob = src->surface.bonds_full;
sb = g_slist_length(dest->bonds);
c = dest->surface.region[0]+dest->surface.region[1];
dest->surface.bonds_cut = ob*c - sb;

/* printout of bonding analysis */
/*
printf("[%2d %2d %2d] (%f)", h, k, l, shift);
printf(" [src = %d][surf = %d][r1+r2 = %d] : cleaved = %d : ", ob, sb, c, dest->surface.bonds_cut);
printf(" : area = %f", dest->area);
printf(" : cleaved = %d", dest->surface.bonds_cut);
printf("\n");
*/

return(0);
}

/****************************************/
/* construct and initialize the surface */
/****************************************/
#define MAKE_SURFACE_DEBUG 0
gpointer make_surface(struct model_pak *data,
                      struct plane_pak *pdata,
                      struct shift_pak *sdata)
{
gint bmode, r1size;
gdouble sbe;
GSList *list;
struct model_pak *surf;
struct core_pak *core;

/* checks */
g_assert(data != NULL);
if (data->periodic != 3)
  {
  gui_text_show(ERROR, "base model is not 3D periodic.\n");
  return(NULL);
  }
if (!pdata->index[0] && !pdata->index[1] && !pdata->index[2])
  {
  gui_text_show(ERROR, "Don't be silly.\n");
  return(NULL);
  }

#if MAKE_SURFACE_DEBUG
printf("Generating surface for model: %d\n", source);
#endif

/* ignore mols when building surface? */
bmode = data->build_molecules;
if (data->surface.ignore_bonding)
  {
  data->build_molecules = FALSE;
  connect_bonds(data);
  connect_molecules(data);
  }

/* allocate & init for surface data */
surf = model_new();

/* NEW - label it as MARVIN, so it's build mode follows the */
/* source model, rather than the GULP setup data - see model_prep() */
surf->id = MARVIN;

/* transfer appropriate GULP setup data to new model */
gulp_data_copy(data, surf);

surf->gulp.run = E_SINGLE;
surf->gulp.method = CONV;

/* transfer surface data to new model */
/* NB: memcpy would produce undesired effects if pointers were */
/* later added to the surface struct */
ARR3SET(surf->surface.miller, pdata->index);
surf->surface.shift = sdata->shift;
surf->surface.region[0] = sdata->region[0];
surf->surface.region[1] = sdata->region[1];
/* setup for region display */
if (surf->surface.region[0])
  surf->region_empty[REGION1A] = FALSE;
if (surf->surface.region[1])
  surf->region_empty[REGION2A] = FALSE;

surf->build_molecules = bmode;

#if MAKE_SURFACE_DEBUG
printf("Miller (%d,%d,%d)\n",surf->surface.miller[0]
                            ,surf->surface.miller[1]
                            ,surf->surface.miller[2]);
printf("Shift %f\n",surf->surface.shift);
printf("Region sizes (%d,%d)\n",surf->surface.region[0]
                               ,surf->surface.region[1]);
#endif

/* TODO - valid shift/region size loop??? */
generate_surface(data, surf);
gulp_files_init(surf);

/*
coords_init(INIT_COORDS, surf);
*/

/* restore connectivity in source model */
if (data->surface.ignore_bonding)
  {
  data->build_molecules = bmode;
  connect_bonds(data);
  connect_molecules(data);
  }

/* calculate the bulk energy needed for GULP surface calcs */
sbe = data->gulp.energy;
if (fabs(sbe) < FRACTION_TOLERANCE)
  {
  gui_text_show(WARNING, "Suspicious total energy. Has it been properly calculated?\n");
  }

/* calc the number of region 1 atoms */
r1size=0;
for (list=surf->cores ; list ; list=g_slist_next(list))
  {
  core = list->data;
  if (core->region == REGION1A)
    r1size++;
  }

sbe /= data->num_atoms;    /* E per atom in unit cell */
sbe *= r1size;             /* E scaled to region 1 size */
surf->gulp.sbulkenergy = sbe;

/* employ src file? */
tree_model_add(surf);
tree_select_model(surf);

return(surf);
}

/*************************************/
/* molecule centroid z value sorting */
/*************************************/
gint surf_molecules_zsort(gpointer p1, gpointer p2)
{
struct mol_pak *m1 = p1, *m2 = p2;

if (m1->centroid[2] > m2->centroid[2])
  return(-1);

return(1);
}

/********************************************************/
/* scan a plane's shift list for the best energy values */
/********************************************************/
#define DEBUG_UPDATE_PLANE_ENERGY 0
void update_plane_energy(struct plane_pak *plane, struct model_pak *model)
{
gint i;
GSList *list;
struct plane_pak *plane2;
struct shift_pak *sdata;

/* checks */
g_assert(model != NULL);
g_assert(plane != NULL);

#if DEBUG_UPDATE_PLANE_ENERGY
printf("(%d %d %d)\n", plane->index[0], plane->index[1], plane->index[2]);
#endif

/* set best values to those of the 1st shift */
list = plane->shifts;
if (list)
  {
  sdata = list->data;

/* init the plane's best values */
  for (i=0 ; i<2 ; i++)
    {
    plane->esurf[i] = sdata->esurf[i];
    plane->esurf_shift = sdata->shift;
    plane->eatt[i] = sdata->eatt[i];
    plane->eatt_shift = sdata->shift;
    plane->bbpa = sdata->bbpa;
    plane->bbpa_shift = sdata->shift;
    }
  }
else
  return;

/* search the shift list for better values */
while (list != NULL)
  {
  sdata = list->data;
  for (i=0 ; i<2 ; i++)
    {
/* surface energy */
    if (sdata->esurf[i] < plane->esurf[i])
      {
      plane->esurf[i] = sdata->esurf[i];
      plane->esurf_shift = sdata->shift;
      }
/* attachment energy */
    if (sdata->eatt[i] > plane->eatt[i])
      {
      plane->eatt[i] = sdata->eatt[i];
      plane->eatt_shift = sdata->shift;
      }
    }
/* broken bonds */
  if (sdata->bbpa < plane->bbpa)
    {
    plane->bbpa = sdata->bbpa;
    plane->bbpa_shift = sdata->shift;
    }
  list = g_slist_next(list);
  }

#if DEBUG_UPDATE_PLANE_ENERGY
printf(" %f\n", plane->bbpa);
#endif

/* update symmetry related planes */
if (plane->primary)
  {
  for (list=model->planes ; list ; list=g_slist_next(list))
    {
    plane2 = list->data;

    if (plane2->parent == plane)
      {
      plane2->esurf[0] = plane->esurf[0];
      plane2->esurf[1] = plane->esurf[1];
      plane2->eatt[0] = plane->eatt[0];
      plane2->eatt[1] = plane->eatt[1];
      }
    }
  }
}

/*******************************************************************/
/* peturb a planar surface in order to test alternate terminations */
/*******************************************************************/
void surf_shift_explore(struct model_pak *src, struct surface_pak *surface)
{
gint i, j, n;
GSList *list;
struct plane_pak *plane;
struct shift_pak *shift;
struct core_pak *core;
struct mol_pak *mol1, *mol2;
struct model_pak *model;

g_assert(src != NULL);
g_assert(src->periodic == 3);

/* setup plane and shift */
plane = plane_new(surface->miller, src);
shift = shift_new(surface->shift);
shift->region[0] = surface->region[0];
shift->region[1] = surface->region[1];

printf("Perturbing [%d %d %d] : %f : [%d][%d]\n",
        plane->index[0], plane->index[1], plane->index[2],
        shift->shift, shift->region[0], shift->region[1]);

/* CURRENT - create a new model for each alternate (zero dipole) surface found */

/* generate surface */
model = make_surface(src, plane, shift);  
coords_init(INIT_COORDS, model);

/* sort molecule centroids */
model->moles = g_slist_sort(model->moles, (gpointer) surf_molecules_zsort);

printf("initial dipole: %f\n", model->gulp.sdipole);

i = 0;

/* top */
mol1 = g_slist_nth_data(model->moles, i);
n = g_slist_length(model->moles);

mol2 = g_slist_nth_data(model->moles, n-1);

/* nth bottommost */
/*
bottom = g_slist_copy(model->moles);
bottom = g_slist_reverse(bottom);
*/

/* for each topmost item, exchange with a bottomost item and recalculate dipole */

/* move all atoms in selection */
for (list=mol1->cores ; list ; list=g_slist_next(list))
  {
  core = list->data;
  region_move_atom(core, DOWN, model); 
  }

mol2 = g_slist_nth_data(model->moles, n-1);

for (j=0 ; j<3 ; j++)
  {

mol2 = g_slist_nth_data(model->moles, n-j-1);

/* test */
for (list=mol2->cores ; list ; list=g_slist_next(list))
  {
  core = list->data;
  region_move_atom(core, UP, model); 

/*
P3VEC("c1: ", core->x);
*/
  }

/* update */
/* FIXME - something wierd going on here - dipole recalc doesnt seem */
/* to work unless connect_molecules() is called - why??? */
coords_compute(model);
/*
connect_bonds(model);
*/
connect_molecules(model);

/* TODO - don't recalc from scratch - can update based on moves */
calc_emp(model);

printf("pert [%d,%d] new dipole: %f\n", i, j, model->gulp.sdipole);

/* TODO - if zero -> leave this model (after init) and create a new */
/* one for exploring further perturbations */

/* restore */
for (list=mol2->cores ; list ; list=g_slist_next(list))
  {
  core = list->data;
  region_move_atom(core, DOWN, model); 

/*
P3VEC("c2: ", core->x);
*/
  }

coords_compute(model);

  }

gui_refresh(GUI_CANVAS);
}

/******************************************************************************
 * GCD
 *      Greatest common denominator
 ******************************************************************************/
gint GCD(gint p, gint q)
{
#if defined(__MWERKS__)  /* fix for compiler bug! */
  gint  i;  
  
  if (q == 0) {
    i = abs(p);
    return(i);
  }
#else
  if (q == 0)
    return(abs(p));
#endif
  else
    return(GCD(q, p%q));
}

/******************************************************************************
 * vector_compare
 *      compare two vectors returning their difference as an integer.  This
 *      routine is to be used in conjuction with SORT().
 ******************************************************************************/
gint vector_compare(vector *a, vector *b)
{
int i;
double  diff;
  
diff = V_MAGSQ(*a) - V_MAGSQ(*b);
if (diff < -EPSILON)
  return(-1);
else if (diff > EPSILON)
  return(1);
else
  {
/* magnitude of a & b is the same, so sort based on element magnitudes */ 
  for (i=0 ; i<3 ; i++)
    {
    diff = V_E(*a,i) - V_E(*b,i);
    if (diff < -EPSILON)
      return(1);
    else if (diff > EPSILON)
      return(-1);
    }
  }
/* a & b are identical */
return(0);
}

