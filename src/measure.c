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

#include <math.h>
#include <stdio.h>
#include <string.h>

#include "gdis.h"
#include "coords.h"
#include "matrix.h"
#include "measure.h"
#include "select.h"
#include "interface.h"
#include "opengl.h"

/* globals */
extern struct sysenv_pak sysenv;
extern struct elem_pak elements[];

/* data structure */
#define MEASURE_MAX_CORES 4

struct measure_pak
{
gint type;                                /* distance, angle, etc. */
gchar *value;
gpointer core[MEASURE_MAX_CORES];         /* participants */
gint image[MEASURE_MAX_CORES][3];         /* perioidic image offset for each core */
gdouble colour[3];
};

/*****************/
/* debugging aid */
/*****************/
#define DEBUG_DUMP_ALL 0
void measure_dump_all(struct model_pak *model)
{
GSList *list;
struct measure_pak *mp;
struct core_pak *c1, *c2, *c3, *c4;

g_assert(model != NULL);

for (list=model->measure_list ; list ; list=g_slist_next(list))
  {
  mp = list->data;

  switch (mp->type)
    {
    case MEASURE_BOND:
    case MEASURE_INTER:
    case MEASURE_INTRA:
    case MEASURE_DISTANCE:
      c1 = mp->core[0];
      c2 = mp->core[1];
      printf("[%s-%s] : ", c1->atom_label, c2->atom_label);
#if DEBUG_DUMP_ALL
      printf("[%d %d %d] ", mp->image[0][0], mp->image[0][1], mp->image[0][2]);
      printf("[%d %d %d] : ", mp->image[1][0], mp->image[1][1], mp->image[1][2]);
#endif
      printf("%s\n", mp->value);
      break;

    case MEASURE_ANGLE:
      c1 = mp->core[0];
      c2 = mp->core[1];
      c3 = mp->core[2];
      printf("[%s-%s-%s] : %s\n", c1->atom_label,
                                  c2->atom_label,
                                  c3->atom_label,
                                  mp->value);
      break;

    case MEASURE_TORSION:
      c1 = mp->core[0];
      c2 = mp->core[1];
      c3 = mp->core[2];
      c4 = mp->core[3];
      printf("[%s-%s-%s-%s] : %s\n", c1->atom_label,
                                     c2->atom_label,
                                     c3->atom_label,
                                     c4->atom_label,
                                     mp->value);
      break;
    }
  }
}

/*****************************************/
/* select the active models measurements */
/*****************************************/
void measure_select_all(void)
{
GSList *list;
struct model_pak *model;
struct measure_pak *m;

model = sysenv.active_model;
if (!model)
  return;

for (list=model->measure_list ; list ; list=g_slist_next(list))
  {
  m = list->data;
  switch (m->type)
    {
    case MEASURE_INTER:
    case MEASURE_INTRA:
    case MEASURE_DISTANCE:
      select_core(m->core[0], FALSE, model);
      select_core(m->core[1], FALSE, model);
      break;
    case MEASURE_ANGLE:
      select_core(m->core[0], FALSE, model);
      select_core(m->core[1], FALSE, model);
      select_core(m->core[2], FALSE, model);
      break;
    }
  }
redraw_canvas(SINGLE);
}

/*******************************/
/* measurement type extraction */
/*******************************/
gint measure_type_get(gpointer data)
{
struct measure_pak *mp = data;

g_assert(mp != NULL);

return(mp->type);
}

/********************************/
/* measurement value extraction */
/********************************/
gchar *measure_value_get(gpointer data)
{
struct measure_pak *mp = data;

g_assert(mp != NULL);

return(mp->value);
}

/********************************/
/* measurement cores extraction */
/********************************/
GSList *measure_cores_get(gpointer data)
{
GSList *list=NULL;
struct measure_pak *mp = data;

g_assert(mp != NULL);

switch (mp->type)
  {
  case MEASURE_TORSION:
    list = g_slist_prepend(list, mp->core[3]);
  case MEASURE_ANGLE:
    list = g_slist_prepend(list, mp->core[2]);
  default:
    list = g_slist_prepend(list, mp->core[1]);
    list = g_slist_prepend(list, mp->core[0]);
  }
return(list);
}

/************************************/
/* measurement type label creation */
/***********************************/
gchar *measure_type_label_create(gpointer data)
{
gchar *text;
struct measure_pak *mp = data;

g_assert(mp != NULL);

switch(mp->type)
  {
  case MEASURE_BOND:
    text = g_strdup("Bond");
    break;

  case MEASURE_DISTANCE:
    text = g_strdup("Dist");
    break;

  case MEASURE_INTER:
    text = g_strdup("Inter");
    break;

  case MEASURE_INTRA:
    text = g_strdup("Intra");
    break;

  case MEASURE_ANGLE:
    text = g_strdup("Angle");
    break;

  case MEASURE_TORSION:
    text = g_strdup("Torsion");
    break;

  default:
    text = g_strdup("Unknown");
  }
return(text);
}

/*************************************/
/* measurement constituents creation */
/*************************************/
gchar *measure_constituents_create(gpointer data)
{
gchar *text;
struct measure_pak *mp = data;
struct core_pak *core[MEASURE_MAX_CORES];

g_assert(mp != NULL);

switch(mp->type)
  {
  case MEASURE_BOND:
  case MEASURE_DISTANCE:
  case MEASURE_INTER:
  case MEASURE_INTRA:
    core[0] = mp->core[0];
    core[1] = mp->core[1];
g_assert(core[0] != NULL);
g_assert(core[1] != NULL);
    text = g_strdup_printf("%4s : %-4s", core[0]->atom_label,
                                         core[1]->atom_label);
    break;

  case MEASURE_ANGLE:
    core[0] = mp->core[0];
    core[1] = mp->core[1];
    core[2] = mp->core[2];

g_assert(core[0] != NULL);
g_assert(core[1] != NULL);
g_assert(core[2] != NULL);

    text = g_strdup_printf("%4s : %-4s : %-4s", core[0]->atom_label,
                                                core[1]->atom_label,
                                                core[2]->atom_label);
    break;

  case MEASURE_TORSION:
    core[0] = mp->core[0];
    core[1] = mp->core[1];
    core[2] = mp->core[2];
    core[3] = mp->core[3];
g_assert(core[0] != NULL);
g_assert(core[1] != NULL);
g_assert(core[2] != NULL);
g_assert(core[3] != NULL);
    text = g_strdup_printf("%4s : %-4s : %-4s : %-4s", core[0]->atom_label,
                                                       core[1]->atom_label,
                                                       core[2]->atom_label,
                                                       core[3]->atom_label);
    break;

  default:
    text = g_strdup(" ? ");
  }
return(text);
}

/*********************************/
/* measurement colour extraction */
/*********************************/
void measure_colour_get(gdouble *colour, gpointer data)
{
struct measure_pak *mp = data;

ARR3SET(colour, mp->colour);
}

/***************************************************************/
/* extract nth measurement component location (periodic aware) */
/***************************************************************/
/* NB: always returns closest (image aware) positions relative to 0th core */
void measure_coord_get(gdouble *x, gint n, gpointer data, struct model_pak *model)
{
gdouble dx[3];
struct core_pak *core1, *core2;
struct measure_pak *mp = data;

/* checks */
g_assert(mp != NULL);
g_assert(n >= 0);
g_assert(n < MEASURE_MAX_CORES);

core1 = mp->core[0];
ARR3SET(x, core1->rx);
if (n)
  {
  core2 = mp->core[n];
  ARR3SET(dx, core1->x);
  ARR3SUB(dx, core2->x);
  fractional_min(dx, model->periodic);
  vecmat(model->latmat, dx);
  ARR3SUB(x, dx);
  }
}

/*********************************/
/* measurement colour alteration */
/*********************************/
void measure_colour_set(gdouble r, gdouble g, gdouble b, gpointer data)
{
struct measure_pak *mp = data;

VEC3SET(mp->colour, r, g, b);
}

/******************************************/
/* measurement duplicate search primitive */
/******************************************/
#define DEBUG_MEASURE_SEARCH 0
gpointer measure_search(gint type,
                        struct core_pak **c,
                        gdouble *xlat,
                        struct model_pak *model)
{
gint i, j, count, match;
gdouble x[3];
GSList *list;
struct measure_pak *mp;

#if DEBUG_MEASURE_SEARCH
printf("search [%p - %p]", c[0], c[1]);
P3VEC(" : ", &xlat[4]);
#endif

/* search for duplicate */
for (list=model->measure_list ; list ; list=g_slist_next(list))
  {
  mp = list->data;
  if (type == mp->type)
    {
/* setup the required matches */
    switch (mp->type)
      {
      case MEASURE_DISTANCE:
      case MEASURE_INTRA:
      case MEASURE_INTER:
#if DEBUG_MEASURE_SEARCH
printf(" * [%p - %p] [%d %d %d]\n", mp->core[0], mp->core[1],
                                    mp->image[1][0], mp->image[1][1], mp->image[1][2]);
#endif
/* skip, if periodic translation vectors don't match */
        ARR3SET(x, &xlat[4]);
        ARR3SUB(x, &mp->image[1][0]);
        if (VEC3MAGSQ(x) > FRACTION_TOLERANCE)
          continue;

      case MEASURE_BOND:
        match = 2;
        break;

      case MEASURE_ANGLE:
        match = 3;
        break;

      default:
        match = 4;
      }

g_assert(match <= MEASURE_MAX_CORES);

/* count the matches and return they satisfy */
    count = 0;
    for (i=match ; i-- ; )
      for (j=match ; j-- ; )
        if (c[i] == mp->core[j])
          count++;
    if (count == match)
      return(mp);

/* can this happen? eg some constituent cores are the same */
g_assert(count < match);

    }
  }
#if DEBUG_MEASURE_SEARCH
printf(" - not found.\n");
#endif

return(NULL);
}


/**********************************/
/* measurement creation primitive */
/**********************************/
gpointer measure_new(void)
{
gint i;
struct measure_pak *m;

m = g_malloc(sizeof(struct measure_pak));

for (i=MEASURE_MAX_CORES ; i-- ; )
  {
  m->core[i] = NULL;
  VEC3SET(m->image[i], 0, 0, 0);
  }
return(m);
}

/**********************************/
/* measurement deletion primitive */
/**********************************/
void measure_free(gpointer data, struct model_pak *model)
{
struct measure_pak *mp = data;

g_free(mp->value);

model->measure_list = g_slist_remove(model->measure_list, mp);

g_free(mp);
}

/*******************************/
/* global measurement deletion */
/*******************************/
void measure_free_all(struct model_pak *model)
{
gpointer data;
GSList *list;

list = model->measure_list;
while (list)
  {
  data = list->data;
  list = g_slist_next(list);

  measure_free(data, model);
  }
model->measure_list = NULL;
}

/***************************************/
/* test if measurement contains a core */
/***************************************/
gboolean measure_has_core(struct core_pak *core, gpointer data)
{
gint i;
struct measure_pak *m = data;

for (i=MEASURE_MAX_CORES ; i-- ; )
  if (m->core[i] == core)
    return(TRUE);
return(FALSE);
}

/*********************************/
/* absolute distance computation */
/*********************************/
gdouble measure_distance(gdouble *x1, gdouble *x2)
{
gdouble a[3];

ARR3SET(a, x2);
ARR3SUB(a, x1);

return(VEC3MAG(a));
}

/******************************/
/* absolute angle computation */
/******************************/
gdouble measure_angle(gdouble *x1, gdouble *x2, gdouble *x3)
{
gdouble a[3], b[3];

/* compute arm vectors (center atom == x2) */
ARR3SET(a, x1);
ARR3SUB(a, x2);
ARR3SET(b, x3);
ARR3SUB(b, x2);

return(R2D*via(a,b,3)); 
}

/**********************************/
/* dihedral calculation primitive */
/**********************************/
gdouble measure_dihedral(gdouble *x1, gdouble *x2, gdouble *x3, gdouble *x4)
{
gdouble len, sign;
gdouble a[3], b[3], c[3], n[3];

/* compute end atom vectors (from 1-2 axis) */
ARR3SET(a, x1);
ARR3SUB(a, x2);
ARR3SET(b, x4);
ARR3SUB(b, x3);

/* compute normal to the plane in which the dihedral is to be calc'd */
ARR3SET(n, x2);
ARR3SUB(n, x3);
normalize(n, 3);

/* CURRENT - compute the direction sense (sign) */
crossprod(c, a, b);
if (via(c,n,3) < 0.5*PI)
  sign = -1.0;
else
  sign = 1.0;

/* n[0..2] is the normal of the plane */
/* project 1-0 onto the normal (ie dist to plane) */
len = a[0]*n[0] + a[1]*n[1] + a[2]*n[2];
/* subtract vector to plane from 1-0 to get projection on plane */
a[0] -= len*n[0];
a[1] -= len*n[1];
a[2] -= len*n[2];
/* project 2-3 onto the normal */
len = b[0]*n[0] + b[1]*n[1] + b[2]*n[2];
/* subtract vector to plane from 2-3 to get projection on plane */
b[0] -= len*n[0];
b[1] -= len*n[1];
b[2] -= len*n[2];
/* compute angle between projected vectors */
return(sign * R2D*via(a,b,3)); 
}

/*********************************************/
/* recompute a measurement value for a model */
/*********************************************/
gdouble measure_update_single(gpointer data, struct model_pak *model)
{
gdouble value=0.0;
gdouble x1[3], x2[3], x3[3];
struct measure_pak *mp = data;

/* checks */
g_assert(model != NULL);
g_assert(mp != NULL);

switch (mp->type)
  {
  case MEASURE_BOND:
  case MEASURE_DISTANCE:
  case MEASURE_INTER:
  case MEASURE_INTRA:
/* get measurement component coords */
    measure_coord_get(x1, 0, mp, model);
    measure_coord_get(x2, 1, mp, model);
/* update distance value */
    ARR3SUB(x1, x2);
    value = VEC3MAG(x1);
    g_free(mp->value);
    mp->value = g_strdup_printf("%8.3f", value); 
    break;

  case MEASURE_ANGLE:
/* get measurement component coords */
    measure_coord_get(x1, 0, mp, model);
    measure_coord_get(x2, 1, mp, model);
    measure_coord_get(x3, 2, mp, model);
/* update angle value */
    value = measure_angle(x1, x2, x3);
    g_free(mp->value);
    mp->value = g_strdup_printf("%8.3f", value); 
    break;
  }
return(value);
}

/************************************************/
/* recompute all measurement values for a model */
/************************************************/
void measure_update_global(struct model_pak *model)
{
GSList *list;

for (list=model->measure_list ; list ; list=g_slist_next(list))
  {
  measure_update_single(list->data, model);
  }
}

/*********************************************************/
/* extract a list of atoms that satisfy the filter label */
/*********************************************************/
#define DEBUG_MEAS_FILTER 0
GSList *measure_filter(const gchar *label, struct model_pak *model)
{
gint code;
gboolean strict=FALSE;
GSList *list, *output;
struct core_pak *core;

/* checks */
g_assert(model != NULL);

/* determine if filter label is an element or an exact (strict) label, */
/* or neither - and set up initial list accordingly */
code = elem_symbol_test(label);

if (code)
  {
  if (strlen(label) != strlen(elements[code].symbol))
    strict = TRUE;
  output = g_slist_copy(model->cores);
  }
else
  {
  if (g_ascii_strncasecmp(label, "any", 3) == 0)
    output = g_slist_copy(model->cores);
  else
    output = g_slist_copy(model->selection);
  }

#if DEBUG_MEAS_FILTER
printf("label = %s, strict = %d:\n", label, strict);
#endif

/* remove undesirables from the initial list */
list = output;
while (list)
  {
  core = list->data;
  list = g_slist_next(list);

/* remove if hidden/deleted */
  if (core->status & (HIDDEN | DELETED))
    {
    output = g_slist_remove(output, core);
    continue;
    }

/* remove if core doesn't match the filter element */
  if (code)
    {
    if (strict)
      {
      if (g_ascii_strcasecmp(label, core->atom_label) != 0)
        output = g_slist_remove(output, core);
      }
    else
      {
      if (code != core->atom_code)
        output = g_slist_remove(output, core);
      }
    }
  }

#if DEBUG_MEAS_FILTER
for (list=output ; list ; list=g_slist_next(list))
  {
  core = list->data;
  printf(" [%s]", core->atom_label);
  }
printf("\n");
#endif

return(output);
}

/**********************************/
/* measurement creation primitive */
/**********************************/
/* NB: return NULL if exists - so we don't keep grafting */
gpointer measure_distance_create(gint type,
                                 struct core_pak **c,
                                 gint *image,
                                 struct model_pak *model)
{
gint i;
gdouble xlat[8];
struct measure_pak *mp;

/* checks */
g_assert(model != NULL);

/* search for duplicate */
VEC4SET(&xlat[0], 0.0, 0.0, 0.0, 0.0);
ARR3SET(&xlat[4], image);
xlat[7] = 0.0;
mp = measure_search(type, c, xlat, model);

if (!mp)
  {
/* create a new geometry label */
  mp = measure_new();
  model->measure_list = g_slist_append(model->measure_list, mp);
  mp->type = type;
  VEC3SET(mp->colour,0.8,0.8,0.8);
  for (i=2 ; i-- ; )
    mp->core[i] = c[i];

/* store raw perioidic image offsets */
  VEC3SET(&mp->image[0][0], 0, 0, 0);
  ARR3SET(&mp->image[1][0], image);
  mp->value = NULL;

/* compute actual cartesian offsets (also calculates measurement value) */
  measure_update_single(mp, model);

  return(mp);
  }
return(NULL);
}

/*********************************************************/
/* measurement testing primitive (perioidic image aware) */
/*********************************************************/
gpointer measure_bond_test(struct core_pak **core,
                           gdouble r1, gdouble r2,
                           struct model_pak *model)
{
gint t[3], cutoff_test;
gdouble r1sq, r2sq, rsq;
gdouble x1[3], x2[3];
struct measure_pak *mp=NULL;

/* checks */
g_assert(model != NULL);
g_assert(core[0] != NULL);
g_assert(core[1] != NULL);

/* in user select mode - don't want a cutoff or image check */
if (r2 < FRACTION_TOLERANCE)
  cutoff_test = FALSE;
else
  cutoff_test = TRUE;

/* setup cutoffs */
r1sq = r1*r1;
r2sq = r2*r2;

/* get measurement component coords */
ARR3SET(x1, core[0]->x);
ARR3SET(x2, core[1]->x);

/*
printf("[%s - %s]\n", core[0]->label, core[1]->label);
P3VEC("core1:", x1);
P3VEC("core2:", x2);
P3VEC("xlat1:", &xlat[0]);
P3VEC("xlat2:", &xlat[4]);
*/

/* test if distance satisfies the cutoffs */
ARR3SUB(x1, x2);
fractional_min(x1, model->periodic);
vecmat(model->latmat, x1);
rsq = VEC3MAGSQ(x1);

if (cutoff_test)
  if (rsq < r1sq || rsq > r2sq)
     return(NULL);

VEC3SET(t, 0, 0, 0);
mp = measure_distance_create(MEASURE_BOND, core, t, model);

return(mp);
}

/**************************************/
/* measurement calculation primitives */
/**************************************/
#define DEBUG_BOND_SEARCH 0
void measure_bond_search(const gchar **label, gdouble r1, gdouble r2, struct model_pak *model)
{
GSList *item1, *item2, *list1, *list2;
struct core_pak *core[2];
struct bond_pak *bond;

/* checks */
g_assert(model != NULL);

/* get pattern matching entries */
list1 = measure_filter(label[0], model);
list2 = measure_filter(label[1], model);

/* for all valid first atoms */
for (item1=list1 ; item1 ; item1=g_slist_next(item1))
  {
  core[0] = item1->data;

/* search for valid (bonded) second atom */
  for (item2=core[0]->bonds ; item2 ; item2=g_slist_next(item2))
    {
    bond = item2->data;

    if (!model->build_hydrogen)
      if (bond->type == BOND_HBOND)
        continue;

    if (core[0] == bond->atom1)
      core[1] = bond->atom2;
    else
      core[1] = bond->atom1;

    if (g_slist_find(list2, core[1]))
      {
#if DEBUG_BOND_SEARCH
printf("valid: %s - %s\n", core[0]->atom_label, core[1]->atom_label);
#endif
      measure_bond_test(core, r1, r2, model);
      }
    }
  }
coords_compute(model);

measure_dump_all(model);
}

/**************************************************/
/* add measurements that may span periodic images */
/**************************************************/
#define DEBUG_DISTANCE_TEST 0
gpointer measure_distance_test(gint type,
                               struct core_pak **core,
                               gdouble r1,
                               gdouble r2,
                               struct model_pak *model)
{
gint i, cutoff_test;
gint t[3], limit[3];
gdouble r1sq, r2sq, rsq, x[3], xlat[8];
struct measure_pak *mp=NULL;

/* checks */
g_assert(model != NULL);
g_assert(core[0] != NULL);
g_assert(core[0] != NULL);

/* setup */
r1sq = r1*r1;
r2sq = r2*r2;
VEC3SET(limit, 0, 0, 0);

/* in user select mode - don't want a cutoff or image check */
if (r2 < FRACTION_TOLERANCE)
  cutoff_test = FALSE;
else
  {
  cutoff_test = TRUE;
/* set up image limits */
  for (i=0 ; i<model->periodic ; i++)
    limit[i] = 1 + r2/model->pbc[i];
  }

#if DEBUG_DISTANCE_TEST
printf("lim: %d %d %d (%f - %f)\n", limit[0], limit[1], limit[2], r1, r2);
#endif

VEC4SET(&xlat[0], 0.0, 0.0, 0.0, 0.0);
VEC4SET(&xlat[4], 0.0, 0.0, 0.0, 0.0);

/* periodic image scan */
for (t[0] = -limit[0] ; t[0] <= limit[0] ; t[0]++)
  {
  for (t[1] = -limit[1] ; t[1] <= limit[1] ; t[1]++)
    {
    for (t[2] = -limit[2] ; t[2] <= limit[2] ; t[2]++)
      {
#if DEBUG_DISTANCE_TEST
printf("periodic image: %d %d %d\n", t[0], t[1], t[2]);
#endif

/* only do the intermolecular check if it's not a periodic image */
    if (type == MEASURE_INTER)
      if (!(t[0] | t[1] | t[2]))
        {
        if (core[0]->mol == core[1]->mol)
          continue;
        }

/* get separation for this image */
      ARR3SET(x, core[0]->x);
      ARR3SUB(x, core[1]->x);
      ARR3SUB(x, t);
      vecmat(model->latmat, x);
      rsq = VEC3MAGSQ(x);

/* no cutoff check for manual measurements */
      if (cutoff_test)
        if (rsq < r1sq || rsq > r2sq)
          continue;

#if DEBUG_DISTANCE_TEST
printf("+ %d %d %d : %s\n", t[0], t[1], t[2], label);
#endif

      mp = measure_distance_create(type, core, t, model);
      }
    }
  }
return(mp);
}

/**************************************/
/* measurement calculation primitives */
/**************************************/
void measure_distance_search(const gchar **label,
                             gint type,
                             gdouble r1, gdouble r2,
                             struct model_pak *model)
{
GSList *item1, *item2, *list1, *list2;
struct core_pak *core[2];

/* checks */
g_assert(model != NULL);

/* get pattern matching entries */
list1 = measure_filter(label[0], model);
list2 = measure_filter(label[1], model);

/* for all valid first atoms */
for (item1=list1 ; item1 ; item1=g_slist_next(item1))
  {
  core[0] = item1->data;

/* search for valid second atom */
  for (item2=list2 ; item2 ; item2=g_slist_next(item2))
    {
    core[1] = item2->data;

    if (core[0] == core[1])
      continue;

    measure_distance_test(type, core, r1, r2, model);
    }
  }
coords_compute(model);
}

/**********************************/
/* measurement creation primitive */
/**********************************/
/* NB: return NULL if exists - so we don't keep grafting */
gpointer measure_angle_create(gchar *label,
                              struct core_pak **c,
                              gdouble *xlat,
                              struct model_pak *model)
{
gint i;
struct measure_pak *mp;

/* checks */
g_assert(model != NULL);

/* search for duplicate */
mp = measure_search(MEASURE_ANGLE, c, xlat, model);
if (!mp)
  {
/* create a new geometry label */
  mp = measure_new();
  model->measure_list = g_slist_append(model->measure_list, mp);
  mp->type = MEASURE_ANGLE;
  VEC3SET(mp->colour,0.8,0.8,0.8);
  for (i=3 ; i-- ; )
    {
    mp->core[i] = c[i];
    VEC3SET(&mp->image[i][0], 0, 0, 0);
    }
  if (label)
    mp->value = g_strdup(label);
  else
    mp->value = g_strdup("?");
  return(mp);
  }
return(NULL);
}

/***************************************/
/* perioidic image aware angle testing */
/***************************************/
gpointer measure_angle_test(struct core_pak **core,
                            gdouble a1, gdouble a2,
                            struct model_pak *model)
{
gchar *label;
gdouble angle;
gdouble x1[3], x2[3], x3[3];
gdouble xlat[12];
struct measure_pak *mp=NULL;

/* checks */
g_assert(model != NULL);
g_assert(core[0] != NULL);
g_assert(core[1] != NULL);
g_assert(core[2] != NULL);

/* get cartesian coordinates of component cores */
ARR3SET(x1, core[0]->rx);
ARR3SET(x2, core[1]->rx);
ARR3SET(x3, core[2]->rx);

/* create a new angle new measurement if it satisfies the cutoffs */
angle = measure_angle(x1, x2, x3);
if (angle >= a1 && angle <= a2)
  {
  label = g_strdup_printf("%8.2f", angle); 

  mp = measure_angle_create(label, core, xlat, model);

  g_free(label);
  }
return(mp);
}

/**************************************/
/* measurement calculation primitives */
/**************************************/
void measure_angle_search(const gchar **label, gdouble *range, struct model_pak *model)
{
gdouble r, r12_min, r12_max, r23_min, r23_max;
gdouble x[3];
GSList *item1, *item2, *item3, *list1, *list2, *list3;
struct core_pak *core[3];

/* checks */
g_assert(model != NULL);

/* get pattern matching entries */
list1 = measure_filter(label[0], model);
list2 = measure_filter(label[1], model);
list3 = measure_filter(label[2], model);

/* setup range tests */
r12_min = range[0]*range[0];
r12_max = range[1]*range[1];
r23_min = range[2]*range[2];
r23_max = range[3]*range[3];

/* loop over center atoms */
for (item2=list2 ; item2 ; item2=g_slist_next(item2))
  {
  core[1] = item2->data;

/* search for valid arm 1 */
  for (item1=list1 ; item1 ; item1=g_slist_next(item1))
    {
    core[0] = item1->data;
    if (core[1] == core[0])
      continue;

    ARR3SET(x, core[1]->rx);
    ARR3SUB(x, core[0]->rx);
    r = VEC3MAGSQ(x);

    if (r < r12_min || r > r12_max)
      continue;

/* search for valid arm 2 */
    for (item3=list3 ; item3 ; item3=g_slist_next(item3))
      {
      core[2] = item3->data;
      if (core[2] == core[0])
        continue;
      if (core[2] == core[1])
        continue;

      ARR3SET(x, core[1]->rx);
      ARR3SUB(x, core[2]->rx);
      r = VEC3MAGSQ(x);

      if (r < r23_min || r > r23_max)
        continue;

      measure_angle_test(core, range[4], range[5], model);
      }
    }
  }
coords_compute(model);
}

/**************************************/
/* measurement calculation primitives */
/**************************************/
void measure_bangle_search(const gchar **label, gdouble a1, gdouble a2, struct model_pak *model)
{
gint m1, m2;
GSList *item1, *item2, *item3, *list1, *list2, *list3;
struct core_pak *core[3];
struct bond_pak *bond1, *bond2;

/* checks */
g_assert(model != NULL);

/* get pattern matching entries */
list1 = measure_filter(label[0], model);
list2 = measure_filter(label[1], model);
list3 = measure_filter(label[2], model);

/* find valid center atoms */
for (item2=list2 ; item2 ; item2=g_slist_next(item2))
  {
  core[1] = item2->data;

/* search for a valid arm */
  item1 = core[1]->bonds;
  while (item1)
    {
    bond1 = item1->data;
    item1 = g_slist_next(item1);

    if (!model->build_hydrogen)
      if (bond1->type == BOND_HBOND)
        continue;

    if (core[1] == bond1->atom1)
      core[0] = bond1->atom2;
    else
      core[0] = bond1->atom1;

    m1 = 0;
    if (g_slist_find(list1, core[0]))
      m1 |= 1;
    if (g_slist_find(list3, core[0]))
      m1 |= 2;

/* search for second arm atom if found a valid first arm atom */
    if (m1)
      {
      item3 = core[1]->bonds;
      while (item3)
        {
        bond2 = item3->data;
        item3 = g_slist_next(item3);

        if (!model->build_hydrogen)
          if (bond2->type == BOND_HBOND)
            continue;

        if (core[1] == bond2->atom1)
          core[2] = bond2->atom2;
        else
          core[2] = bond2->atom1;

        if (core[0] == core[2])
          continue;

/* search for valid second arm atom */
        m2 = 0;
        if (g_slist_find(list1, core[2]))
          m2 |= 1;
        if (g_slist_find(list3, core[2]))
          m2 |= 2;

/* check for a valid second match */
/*
printf("%d & %d = ", m1, m2);
*/
        if (m2)
          {
          m2 |= m1;
/*
printf(" %d\n", m2);
*/
          switch (m2)
            {
/* exceptions - no match, or both matches come from the same list */
            case 0:
            case 1:
            case 2:
              break;

            default:
#if DEBUG_ANGLE_SEARCH
printf("testing valid candidate: %s - %s - %s\n", core1->label, core2->label, core3->label);
#endif
              measure_angle_test(core, a1, a2, model);
            }
          }
        }
      }
    }
  }
coords_compute(model);
}

/*******************************/
/* routines for geometry calcs */
/*******************************/
gdouble measure_torsion(struct core_pak **core)
{
gdouble len, sign;
gdouble a[3], b[3], c[3], n[3];

/* compute end atom vectors (from 1-2 axis) */
ARR3SET(a, core[0]->rx);
ARR3SUB(a, core[1]->rx);
ARR3SET(b, core[3]->rx);
ARR3SUB(b, core[2]->rx);

/* compute normal to the plane in which the dihedral is to be calc'd */
ARR3SET(n, core[1]->rx);
ARR3SUB(n, core[2]->rx);
normalize(n, 3);

/* CURRENT - compute the direction sense (sign) */
crossprod(c, a, b);
if (via(c,n,3) < 0.5*PI)
  sign = -1.0;
else
  sign = 1.0;

/* n[0..2] is the normal of the plane */
/* project 1-0 onto the normal (ie dist to plane) */
len = a[0]*n[0] + a[1]*n[1] + a[2]*n[2];
/* subtract vector to plane from 1-0 to get projection on plane */
a[0] -= len*n[0];
a[1] -= len*n[1];
a[2] -= len*n[2];
/* project 2-3 onto the normal */
len = b[0]*n[0] + b[1]*n[1] + b[2]*n[2];
/* subtract vector to plane from 2-3 to get projection on plane */
b[0] -= len*n[0];
b[1] -= len*n[1];
b[2] -= len*n[2];
/* compute angle between projected vectors */
return(sign * R2D*via(a,b,3)); 
}

/********************************/
/* tosion measurement searching */
/********************************/
void measure_torsion_search(gchar **label,
                            gdouble a1, gdouble a2,
                            struct model_pak *model)
{
}

/**********************************/
/* measurement creation primitive */
/**********************************/
/* NB: return NULL if exists - so we don't keep grafting */
gpointer measure_torsion_create(const gchar *label,
                                struct core_pak **core,
                                struct model_pak *model)
{
gint i;
gdouble xlat[3];
struct measure_pak *m;

g_assert(core != NULL);
g_assert(model != NULL);

VEC3SET(xlat, 0.0, 0.0, 0.0);
m = measure_search(MEASURE_TORSION, core, xlat, model);
if (!m)
  {
/* create a new geometry label */
  m = measure_new();
  model->measure_list = g_slist_append(model->measure_list, m);
  m->type = MEASURE_TORSION;
  VEC3SET(m->colour,0.8,0.8,0.8);
  for (i=4 ; i-- ; )
    {
    m->core[i] = core[i];
/*
  ARR3SET(&m->xlat[i][0], &xlat[3*i]);
  VEC3SET(&m->image[i][0], 0, 0, 0);
*/
    }
  if (label)
    m->value = g_strdup(label);
  else
    m->value = g_strdup("?");

  return(m);
  }
else
  printf("already exists.\n");

return(NULL);
}

/*************************************/
/* measurement calculation primitive */
/*************************************/
gpointer measure_torsion_test(struct core_pak **core,
                              gdouble a1, gdouble a2,
                              struct model_pak *model)
{
gchar *label;
gdouble angle;
gpointer m = NULL;

g_assert(core != NULL);
g_assert(model != NULL);

angle = measure_torsion(core);

if (angle >= a1 && angle <= a2)
  {
  label = g_strdup_printf("%8.2f", angle); 
  m = measure_torsion_create(label, core, model);
  g_free(label);
  }
return(m);
}

/************************/
/* match pairs of atoms */
/************************/
#define DEBUG_PAIR_MATCH 0
gint pair_match(const gchar *label1, const gchar *label2,
                struct core_pak *core1, struct core_pak *core2)
{
gint a, b, i, j, mask;

#if DEBUG_PAIR_MATCH
printf("[%s,%s] : [%s,%s]\n", label1, label2, core1->label, core2->label);
#endif

/* if a or b = 0 => any i or j is accepted */
/* otherwise it must match either i or j */
/* if a and b are non zero, both i & j must match both a & b */
a = elem_symbol_test(label1);
b = elem_symbol_test(label2);
i = core1->atom_code;
j = core2->atom_code;

/* fill out the mask */
mask = 0;
if (i == a)
  {
/* if input label doesn't match the element symbol length - it means the */
/* user has put in something like H1 - compare this with the atom label */
  if (g_ascii_strcasecmp(label1, elements[i].symbol) != 0)
    {
    if (g_ascii_strcasecmp(core1->atom_label, label1) == 0)
      mask |= 1;
    }
  else
    mask |= 1; 
  }

if (j == a)
  {
  if (g_ascii_strcasecmp(label1, elements[j].symbol) != 0)
    {
    if (g_ascii_strcasecmp(core2->atom_label, label1) == 0)
      mask |= 2;
    }
  else
    mask |= 2; 
  }

if (i == b)
  {
  if (g_ascii_strcasecmp(label2, elements[i].symbol) != 0)
    {
    if (g_ascii_strcasecmp(core1->atom_label, label2) == 0)
      mask |= 4;
    }
  else
    mask |= 4; 
  }

if (j == b)
  {
  if (g_ascii_strcasecmp(label2, elements[j].symbol) != 0)
    {
    if (g_ascii_strcasecmp(core2->atom_label, label2) == 0)
      mask |= 8;
    }
  else
    mask |= 8; 
  }

#if DEBUG_PAIR_MATCH
printf("mask = %d\n", mask);
#endif

/* if both types must match - only two possibilities (a,b) or (b,a) */
/* but we can get further matches (than the required 2) when labels are compared */
if (a && b)
  {
  switch(mask)
    {
/* single valid pair match */
    case 6:
    case 9:
/* valid pair match plus one extra */
    case 7:
    case 11:
    case 13:
    case 14:
/* pair match in all possible combinations */
    case 15:
      break;
/* bad match - exit */
    default:
      return(0);
    }
  }

/* if only one type to match - any match at all will do */
if (a || b)
  if (!mask)
    return(0);

#if DEBUG_PAIR_MATCH
printf("accepted [%d][%d] as a match.\n",i,j);
#endif

return(1);
}
