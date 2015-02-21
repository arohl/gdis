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
#include <string.h>
#include <strings.h>
#include <stdlib.h>
#include <math.h>
#ifdef __APPLE__
#include <OpenGL/gl.h>
#else
#include <GL/gl.h>
#endif
#include "gdis.h"
#include "coords.h"
#include "interface.h"
#include "gui_shorts.h"
#include "matrix.h"
#include "model.h"
#include "spatial.h"
#include "numeric.h"
#include "morph.h"
#include "opengl.h"
#include "select.h"
#include "zone.h"

/* data structures */
extern struct sysenv_pak sysenv;
extern struct elem_pak elements[];

/*******************************/
/* debugging routines for bonds */
/*******************************/
void dump_atom_bonds(struct model_pak *data)
{
GSList *list1, *list2;
struct core_pak *core;
struct bond_pak *bond;

for (list1=data->cores ; list1 ; list1=g_slist_next(list1))
  {
  core = list1->data;
  printf("core [%s][%p], bonds: ", core->atom_label, core);
  for (list2=core->bonds ; list2 ; list2=g_slist_next(list2))
    {
    bond = list2->data;

    switch (bond->type)
      {
      case BOND_HBOND:
        printf("(h) ");
        break;
      default:
        printf("(n) ");
        break;
      }
    }
  printf("\n");
  }
}

/*******************************/
/* debugging routines for bonds */
/*******************************/
void dump_bonds(struct model_pak *data)
{
gint h,n,m,p;
GSList *item;
struct bond_pak *bond;

h=n=m=p=0;
for (item=data->bonds ; item ; item=g_slist_next(item))
  {
  bond = item->data;

  if (bond->type == BOND_HBOND)
    h++;
  else
    {
    printf("  normal bond [%p][%4s-%4s][%p-%p]\n",
            bond, (bond->atom1)->atom_label, (bond->atom2)->atom_label, 
            bond->atom1, bond->atom2);
    n++;
    }
  }

printf("bond total=%d : normal=%d, periodic=%d, merged=%d, hbond=%d\n",
       g_slist_length(data->bonds),n,p,m,h);
}

/*******************************/
/* debugging routines for mols */
/*******************************/
void dump_mols(struct model_pak *model)
{
GSList *list, *list2;
struct mol_pak *mol;
struct core_pak *core;

printf("Found %d molecules.\n", g_slist_length(model->moles));

for (list=model->moles ; list ; list=g_slist_next(list))
  {
  mol = list->data;

  printf("mol [%p] has %d cores: ", mol, g_slist_length(mol->cores));

  for (list2=mol->cores ; list2 ; list2=g_slist_next(list2))
    {
    core = list2->data;

    printf("%s ", core->atom_label);
    }
  printf("\n");
  }
}

/**************************************/
/* get list of atoms bonded to a core */
/**************************************/
/* TODO - implement with a level parameter that gives next-nearest etc. */
GSList *connect_neighbours(struct core_pak *core) 
{
GSList *list, *neighbours=NULL;
struct bond_pak *bond;
struct core_pak *core2;

for (list=core->bonds ; list ; list=g_slist_next(list))
  {
  bond = list->data;
  if (bond->type == BOND_HBOND)
    continue;

/* get bonded cores */
  if (core == bond->atom1)
    core2 = bond->atom2;
  else
    core2 = bond->atom1;

  if (core2->status & HIDDEN)
    continue;

  neighbours = g_slist_prepend(neighbours, core2);
  }

return(neighbours);
}

/*******************************************/
/* determine if a bond is in a split state */
/*******************************************/
gint connect_split(struct bond_pak *bond)
{
gdouble x[3];
struct core_pak *c1, *c2;

g_assert(bond != NULL);

c1 = bond->atom1;
c2 = bond->atom2;

/* compare actual coords to bond midpoint offset */
ARR3SET(x, c2->x);
ARR3SUB(x, c1->x);
ARR3SUB(x, bond->offset);
ARR3SUB(x, bond->offset);

/* if magnitude is 0 - bond is not split */
if (VEC3MAGSQ(x) < FRACTION_TOLERANCE)
  return(FALSE);
return(TRUE);
}

/**************************************************************/
/* check for bonds "ruptured" by a periodic image translation */
/**************************************************************/
gint connect_split_check_all(struct model_pak *model)
{
GSList *list;
struct bond_pak *bond;

g_assert(model != NULL);

/* NB: assumes the bonds are still valid */
for (list=model->bonds ; list ; list=g_slist_next(list))
  {
  bond = list->data;
  if (bond->type == BOND_HBOND)
    continue;

  if (connect_split(bond))
    return(TRUE);
  }
return(FALSE);
}

/**************************************************************/
/* check for bonds "ruptured" by a periodic image translation */
/**************************************************************/
gint connect_split_list(GSList *cores)
{
GSList *list1, *list2;
struct bond_pak *bond;
struct core_pak *core;

for (list1=cores ; list1 ; list1=g_slist_next(list1))
  { 
  core = list1->data;

/* NB: assumes the bonds are still valid */
  for (list2=core->bonds ; list2 ; list2=g_slist_next(list2))
    {
    bond = list2->data;
    if (bond->type == BOND_HBOND)
      continue;

    if (connect_split(bond))
      return(TRUE);
    }
  }
return(FALSE);
}

/****************************/
/* display hydrogen bonding */
/****************************/
void build_hbonds(struct model_pak *model)
{
gint status;
GSList *list;
struct bond_pak *bond;

if (model->build_hydrogen)
  status = NORMAL;
else
  status = HIDDEN;

for (list=model->bonds ; list ; list=g_slist_next(list))
  {
  bond = list->data;

  if (bond->type == BOND_HBOND)
    bond->status = status;
  }
}

/************************************/
/* build zeolite style connectivity */
/************************************/
void build_zeolite(struct model_pak *data)
{
GSList *list, *list1, *list2;
struct bond_pak *bond;
struct core_pak *core, *core1, *core2;

for (list=data->cores ; list ; list=g_slist_next(list))
  {
  core = list->data;

/* search for O */
  if (core->atom_code != 8)
    continue;

  for (list1=core->bonds ; list1 ; list1=g_slist_next(list1))
    {
    core1 = list1->data;

/* search for Al or Si */
    if (core1->atom_code == 13 || core1->atom_code == 14)
      {
/* search for Al or Si */
      list2 = g_slist_next(list1);
      while (list2)
        {
        core2 = list2->data;

        if (core2->atom_code == 13 || core2->atom_code == 14)
          {
/* hide the oxygen */
          core->status |= ZEOL_HIDDEN;
/* replace with one zeolite connection */
          bond = g_malloc(sizeof(struct bond_pak)); 
          bond->type = ZEOLITE;
          bond->atom1 = core1;
          bond->atom2 = core2;
          data->bonds = g_slist_prepend(data->bonds, bond);
          }
        list2 = g_slist_next(list2);
        }
      }
    }
  }
}

/***************************************/
/* build polyhedral style connectivity */
/***************************************/
#define DEBUG_BUILD_POLYHEDRAL 0
void create_polyhedra(struct model_pak *data)
{
gint a, b, c, k, n;
gdouble angle, vec1[3], vec2[3], vec3[3], vec4[3], offset[3];
GSList *list1, *list2, *clist, *mlist, *vlist;
struct spatial_pak *spatial;
struct vec_pak *p1, *p2, *p3, *v1, *v2, *v3;
struct core_pak *core1, *core2;
struct bond_pak *bond;

/* if selection, apply only to those cores, else all cores */
if (data->selection)
  mlist = data->selection;
else
  mlist = data->cores;

/* NEW - one spatial for all the triangle vertices */
spatial = spatial_new("polyhedra", SPATIAL_GENERIC, 3, TRUE, data);

/* iterate over candidate center atoms */
for (list1=mlist ; list1 ; list1=g_slist_next(list1))
  {
  core1 = list1->data;
  ARR3SET(vec1, core1->x);

/* search for bonds */
  n=0;
  clist=NULL;
  for (list2=core1->bonds ; list2 ; list2=g_slist_next(list2))
    {
    bond = list2->data;

    ARR3SET(offset, bond->offset);
    if (core1 == bond->atom1)
      {
      core2 = bond->atom2;
      VEC3MUL(offset, 2.0);
      }
    else
      {
      core2 = bond->atom1;
      VEC3MUL(offset, -2.0);
      }

/* NEW - ignore bonds between same atom types (both candidate center atoms) */
    if (core1->atom_code == core2->atom_code)
      continue;

/* candidate polyhedral vertex */
    p1 = g_malloc(sizeof(struct vec_pak));
    p1->data = NULL;
    ARR3SET(p1->x, core1->x);
    ARR3ADD(p1->x, offset);

    clist = g_slist_prepend(clist, p1);
    n++;
    }

#if DEBUG_BUILD_POLYHEDRAL
printf("This [%s] has %d candidates.\n", core1->atom_label, n);
#endif

/* minimum of 4 points for a closed polyhedron */
  if (n > 3)
    {
    k=0;
    vlist = NULL;
/* replace all triples with a triangle */
    for (a=0 ; a<n-2 ; a++)
    for (b=a+1 ; b<n-1 ; b++)
    for (c=b+1 ; c<n ; c++)
      {
      p1 = g_slist_nth_data(clist, a);
      p2 = g_slist_nth_data(clist, b);
      p3 = g_slist_nth_data(clist, c);

/* central atom's colour */
      ARR3SET(p1->colour, core1->colour);
      VEC3MUL(p1->colour, 1.0/65535.0);
      ARR3SET(p2->colour, core1->colour);
      VEC3MUL(p2->colour, 1.0/65535.0);
      ARR3SET(p3->colour, core1->colour);
      VEC3MUL(p3->colour, 1.0/65535.0);

/* triangle, defined by 3 surrounding bonds */
/* compute midpoints */
      ARR3SET(vec2, p1->x);
      ARR3ADD(vec2, p2->x);
      VEC3MUL(vec2, 0.5);

      ARR3SET(vec3, p1->x);
      ARR3ADD(vec3, p3->x);
      VEC3MUL(vec3, 0.5);

      ARR3SET(vec4, p2->x);
      ARR3ADD(vec4, p3->x);
      VEC3MUL(vec4, 0.5);

/* compare with central atom */
      ARR3SUB(vec2, vec1);
      ARR3SUB(vec3, vec1);
      ARR3SUB(vec4, vec1);
      vecmat(data->latmat, vec2);
      vecmat(data->latmat, vec3);
      vecmat(data->latmat, vec4);

/* if one matches, then this is an internal facet, and will be discarded */
/*
      if (VEC3MAGSQ(vec2) < POSITION_TOLERANCE)
        flag++;
      if (VEC3MAGSQ(vec3) < POSITION_TOLERANCE)
        flag++;
      if (VEC3MAGSQ(vec4) < POSITION_TOLERANCE)
        flag++;
*/
      if (VEC3MAGSQ(vec2) < 0.5)
        continue;
      if (VEC3MAGSQ(vec3) < 0.5)
        continue;
      if (VEC3MAGSQ(vec4) < 0.5)
        continue;

/* TODO - a general routine, given a list of points, and a */
/* reference -> reorder (if nec.) so this is true */
/* compute the clockwise normal */
      ARR3SET(vec2, p2->x);
      ARR3SUB(vec2, p1->x);
      ARR3SET(vec3, p3->x);
      ARR3SUB(vec3, p1->x);
      crossprod(vec4, vec2, vec3);
/* CURRENT */
/* keep the normal */
      normalize(vec4, 3);
      ARR3SET(p1->n, vec4);
      ARR3SET(p2->n, vec4);
      ARR3SET(p3->n, vec4);
/* get a vector from poly center to face */
      ARR3SUB(vec2, p1->x);
      VEC3MUL(vec2, -1.0);
      angle = via(vec2, vec4, 3);

/* copy the vertices (since they may be free'd individually later) */
/* TODO - implement copy/duplicate etc. vector method */
      v1 = g_malloc(sizeof(struct vec_pak));
      v2 = g_malloc(sizeof(struct vec_pak));
      v3 = g_malloc(sizeof(struct vec_pak));
v1->data = NULL;
v2->data = NULL;
v3->data = NULL;
      ARR3SET(v1->x, p1->x);
      ARR3SET(v2->x, p2->x);
      ARR3SET(v3->x, p3->x);
      ARR3SET(v1->n, p1->n);
      ARR3SET(v2->n, p2->n);
      ARR3SET(v3->n, p3->n);
      ARR3SET(v1->colour, p1->colour);
      ARR3SET(v2->colour, p2->colour);
      ARR3SET(v3->colour, p3->colour);

/* enforce an outwards pointing clockwise normal */
/* TODO - prepend for speed */
      if (angle < PI/2.0)
        {
        vlist = g_slist_prepend(vlist, v3);
        vlist = g_slist_prepend(vlist, v2);
        vlist = g_slist_prepend(vlist, v1);
        }
      else
        {
        vlist = g_slist_prepend(vlist, v2);
        vlist = g_slist_prepend(vlist, v3);
        vlist = g_slist_prepend(vlist, v1);
        VEC3MUL(v1->n, -1.0);
        VEC3MUL(v2->n, -1.0);
        VEC3MUL(v3->n, -1.0);
        }
      k++;
      }

/* only add if more than 3 facets found */
/* ie min for a closed polyhedron - not infallible but good enough??? */
    if (k > 3)
      spatial->list = g_slist_concat(spatial->list, vlist);
    else
      free_slist(vlist);
    }

/* free all vertices */
/* TODO - free vertex data */
  g_slist_free(clist);
  }

/* TODO - check if spatial->list is NULL (no polyhedra found) & free if so */

}

/********************/
/* duplicate a bond */
/********************/
struct bond_pak *dup_bond(struct bond_pak *bond)
{
struct bond_pak *copy;

copy = g_malloc(sizeof(struct bond_pak));

memcpy(copy, bond, sizeof(struct bond_pak));

return(copy);
}

/*************************************************/
/* check for and creates if sucessful a new bond */
/*************************************************/
gpointer connect_bond_test(gdouble *x, gdouble r2cut,
                           struct core_pak *c1, struct core_pak *c2,
                           struct model_pak *model)
{
gint type, status, test;
gdouble r2, r[3];
struct bond_pak *bond;

g_assert(c1 != NULL);
g_assert(c2 != NULL);
g_assert(model != NULL);

/* test cartesian separation */
ARR3SET(r, x);
vecmat(model->latmat, r);
r2 = VEC3MAGSQ(r);

if (r2 < r2cut)
  {
/* avoid the too close case - should mainly happened if c1 == c2 (and neither is an image) */
/* NB: c1 == c2 is now allowed in order to trap bonding due to periodic images */
  if (r2 < POSITION_TOLERANCE)
    return(NULL);
/* normal */
  status = NORMAL;
  type = BOND_SINGLE;
  }
else
  {
/* hbond candidate */
  if (model->show_hbonds)
    status = NORMAL;
  else
    status = HIDDEN;

  type = BOND_HBOND;
/* test for H - N, O, F */
  test = 0;
  if (c1->atom_code == 1)
    test = c2->atom_code;
  if (c2->atom_code == 1)
    test = c1->atom_code;
  switch (test)
    {
    case 7:
    case 8:
    case 9:
      if (r2 < 6.25)
        break;
    default:
      return(NULL);
    }
  }
if (type == BOND_HBOND)
  if (!(c1->hydrogen_bond || c2->hydrogen_bond))
    return(NULL);

/* create a new bond */
bond = g_malloc(sizeof(struct bond_pak)); 
bond->type = type;
bond->status = status;
bond->atom1 = c1;
bond->atom2 = c2;

/* bond midpoint fractional offset (relative to core1) */
ARR3SET(bond->offset, x);
VEC3MUL(bond->offset, 0.5);

/* create bond references */
c1->bonds = g_slist_prepend(c1->bonds, bond);
c2->bonds = g_slist_prepend(c2->bonds, bond);
model->bonds = g_slist_prepend(model->bonds, bond);

return(bond);
}

/********************************/
/* update bond list for an atom */
/********************************/
/* NEW - using the new zone scheme, the cell translation bit is */
/* encapsulated in zone lookup */
#define DEBUG_ATOM_BONDS 0
void connect_atom_compute(struct core_pak *core1, struct model_pak *model)
{
gint i;
gdouble r1, r2cut;
gdouble x[3], x1[3], x2[3];
gpointer zone;
GSList *list, *locality;
struct core_pak *core2;

g_assert(core1 != NULL);
g_assert(model != NULL);

/* allow connectivity to be turned off */
if (!model->build_molecules)
  return;

#if DEBUG_ATOM_BONDS
printf("[%s] [%f %f %f] ===================\n",
core1->atom_label, core1->x[0], core1->x[1], core1->x[2]);
#endif

/* build a list of neighbour zones (some may be images) */
zone = zone_get(core1->x, model->zone_array);
locality = zone_area_cores(1, zone, model->zone_array);

/* setup for compare */
r1 = core1->bond_cutoff;
ARR3SET(x1, core1->x);

for (list=locality ; list ; list=g_slist_next(list))
  {
  core2 = list->data;

  if (core2->status & DELETED)
    continue;
/* avoid double counting */
  if (core1 > core2)
    continue;

/* calc bond cutoff*/
  r2cut = r1 + core2->bond_cutoff;
  r2cut += BOND_FUDGE*r2cut;
  r2cut *= r2cut;

/* get the minimum separation */
/* FIXME - this doesn't always give the min CARTESIAN separation */
/* ie can happen when the lattice matrix has large non-diagonal elements */
  ARR3SET(x2, core2->x);
  ARR3SUB(x2, x1);
  fractional_min(x2, model->periodic);

#if DEBUG_ATOM_BONDS
printf(" [%s ?< %f]\n", core2->atom_label, r2cut);
P3VEC(" : ", x2);
#endif

/* the minimum separation test */
/*
  if (connect_bond_test(x2, r2cut, core1, core2, model))
    {
*/
/* CURRENT - hack for the above FIXME */
/* ie do periodic image checks in all cases - should trap most cases */
      connect_bond_test(x2, r2cut, core1, core2, model);

/* NEW - test for multiple (due to periodicity) connections */
/* technically this search is not exhaustive */
/* eg should also check (1,1,1) diagonal offsets */
/* fortunately, this is unlikely to affect most cases */
    for (i=0 ; i<model->periodic ; i++)
      {
/* check +ve direction */
      ARR3SET(x, x2);
      x[i] += 1.0;
      connect_bond_test(x, r2cut, core1, core2, model);
#if DEBUG_ATOM_BONDS
P3VEC(" : ", x);
#endif
/* check -ve direction */
      ARR3SET(x, x2);
      x[i] -= 1.0;
      connect_bond_test(x, r2cut, core1, core2, model);
#if DEBUG_ATOM_BONDS
P3VEC(" : ", x);
#endif
      }
/*
    }
*/
  }
g_slist_free(locality);
}

/***********************************************/
/* modify primary connectivity with user bonds */
/***********************************************/
void connect_merge_user(struct model_pak *data)
{
gint found, match;
gdouble x[3];
GSList *list1, *list2;
struct bond_pak *bond1, *bond2;
struct core_pak *core1, *core2;
   
/* checks */
g_assert(data != NULL);
if (!data->build_molecules)
  return;

/* search & modify if bond already exists */
for (list1=data->ubonds ; list1 ; list1=g_slist_next(list1))
  {
  bond1 = list1->data;
  found = 0;

  for (list2=data->bonds ; list2 ; list2=g_slist_next(list2))
    {
    bond2 = list2->data;

    match = 0;
    if (bond1->atom1 == bond2->atom1 && bond1->atom2 == bond2->atom2)
      match++;
    if (bond1->atom1 == bond2->atom2 && bond1->atom2 == bond2->atom1)
      match++;

/* bond already exists -> modify type */
    if (match)
      {
      if (bond1->status == DELETED)
        bond2->status = DELETED;
      else
        {
        bond2->type = bond1->type;
        bond2->status = NORMAL;
        }
      found++;
      break;
      }
    }

  if (!found)
    {
/* if not found, create a copy for the primary list */
    bond2 = dup_bond(bond1);

/* compute offset */
    core1 = bond2->atom1;
    core2 = bond2->atom2;
    ARR3SET(x, core2->x);
    ARR3SUB(x, core1->x);
    fractional_min(x, data->periodic);
    VEC3MUL(x, 0.5);
    ARR3SET(bond2->offset, x);

/* add list references */
    data->bonds = g_slist_prepend(data->bonds, (gpointer) bond2);
    (bond2->atom1)->bonds = g_slist_prepend((bond2->atom1)->bonds, bond2);
    (bond2->atom2)->bonds = g_slist_prepend((bond2->atom2)->bonds, bond2);
    }
  }
}

/************************************************************/
/* force creation of a bond of specified type between i & j */
/************************************************************/
#define DEBUG_USER_BOND 0
void connect_user_bond(struct core_pak *core1, struct core_pak *core2,
                            gint type, struct model_pak *data)
{
gint count, flag;
struct bond_pak *bdata;
GSList *ubond;

#if DEBUG_USER_BOND
printf("submit bond [%d]: %p - %p : ", type, core1, core2);
#endif
/* checks */
g_assert(core1 != NULL);
g_assert(core2 != NULL);
g_assert(data != NULL);
if (core1 == core2)
  return;

/* search... */
flag=0;
ubond = data->ubonds;
while (ubond)
  {
  bdata = ubond->data;

/* atom match */
  count=0;
  if (core1 == bdata->atom1)
    count++;
  if (core1 == bdata->atom2)
    count++;
  if (core2 == bdata->atom1)
    count++;
  if (core2 == bdata->atom2)
    count++;

/* ...and change if found */
  if (count > 1)
    {
#if DEBUG_USER_BOND
printf(" [modifying]\n");
#endif
    if (type == BOND_DELETE)
      bdata->status = DELETED;
    else
      {
      bdata->type = type; 
      bdata->status = NORMAL;
      }
    flag=1;
    break;
    }
  ubond = g_slist_next(ubond);
  }

/* not found - create */
if (!flag)
  {
#if DEBUG_USER_BOND
printf(" [creating]\n");
#endif
  bdata = g_malloc(sizeof(struct bond_pak));
  bdata->atom1 = core1;
  bdata->atom2 = core2;
  if (type == BOND_DELETE)
    {
    bdata->status = DELETED;
    bdata->type = BOND_SINGLE;
    }
  else
    {
    bdata->status = NORMAL;
    bdata->type = type; 
    }
  data->ubonds = g_slist_prepend(data->ubonds, (gpointer) bdata);
  }

/* update primary bond list */
connect_merge_user(data);
/* update connectivity */
connect_molecules(data);
}

/***********************************/
/* make bonds by clicking on atoms */
/***********************************/
void connect_make_bond(struct core_pak *core, gint type, struct model_pak *model)
{
static struct core_pak *first=NULL;

if (first)
  {
  select_del_core(first, model);
  connect_user_bond(first, core, type, model);
  first = NULL;
  }
else
  {
  first = core;
  select_add_core(core, model);
  }
}

/*******************************************/
/* update connectivity for a list of atoms */
/*******************************************/
#define DEBUG_REDO_BONDS 0
void connect_atom_clear(struct core_pak *core, struct model_pak *data)
{
GSList *list;
struct bond_pak *bond;
struct core_pak *core1, *core2;

/* seek input core in the bond list */
list = core->bonds;
while (list)
  {
  bond = list->data;
  list = g_slist_next(list);

  core1 = bond->atom1;
  core2 = bond->atom2;
/* remove both core references and main list reference, then free it */
  core1->bonds = g_slist_remove(core1->bonds, bond);
  core2->bonds = g_slist_remove(core2->bonds, bond);
  data->bonds = g_slist_remove(data->bonds, bond);
  g_free(bond);
  }
core->bonds=NULL;
}

/****************************/
/* recalc for a single core */
/****************************/
void connect_atom_refresh(struct core_pak *core, struct model_pak *data)
{
connect_atom_clear(core, data);
connect_atom_compute(core, data);
}

/****************************/
/* destroy all connectivity */
/****************************/
/* TODO connect_free() - try loop over all atoms & just g_slist_free bond list */
/* do some timing on wipe_bodns & above to compare */
void wipe_bonds(struct model_pak *data)
{
GSList *list;
struct bond_pak *bond;
struct core_pak *core1, *core2;

list = data->bonds;
while (list)
  {
  bond = list->data;
  list = g_slist_next(list);

  core1 = bond->atom1;
  core2 = bond->atom2;
/* remove both core references and main list reference, then free it */
  core1->bonds = g_slist_remove(core1->bonds, bond);
  core2->bonds = g_slist_remove(core2->bonds, bond);
  data->bonds = g_slist_remove(data->bonds, bond);
  g_free(bond);
  }
data->bonds=NULL;
}

/******************************/
/* redo bond midpoint offsets */
/******************************/
void connect_midpoints(struct model_pak *model)
{
gdouble x[3];
GSList *list;
struct core_pak *core1, *core2;
struct bond_pak *bond;

g_assert(model != NULL);

for (list=model->bonds ; list ; list=g_slist_next(list))
  {
  bond = list->data;

/* compute offset */
  core1 = bond->atom1;
  core2 = bond->atom2;
  ARR3SET(x, core2->x);
  ARR3SUB(x, core1->x);
  fractional_min(x, model->periodic);
  VEC3MUL(x, 0.5);
  ARR3SET(bond->offset, x);
  }
}

/*********/
/* BONDS */
/*********/
#define DEBUG_CALC_BONDS 0
void connect_bonds(struct model_pak *model)
{
GSList *list;
struct core_pak *core;

g_assert(model != NULL);

/* enforce a midpoint recalculation (gulp noautobonds + selection rotate fix) */
connect_midpoints(model);

if (model->anim_fix)
  return;
if (!model->build_molecules)
  return;

g_assert(model->zone_array != NULL);

/* TODO - delete polyhedra */
/* redo atom connectivity */
wipe_bonds(model);

for (list=model->cores ; list ; list=g_slist_next(list))
  {
  core = list->data;

  core->status &= ~ZEOL_HIDDEN;
  if (!(core->status & HIDDEN))
    connect_atom_compute(core, model);
  }

/* user bond modification */
connect_merge_user(model);

/* special building modes */
build_hbonds(model);

if (model->build_zeolite)
  build_zeolite(model);

#if DEBUG_CALC_BONDS
dump_bonds(model);
dump_atom_bonds(model);
#endif
}

/************************************/
/* free molecule list plus its data */
/************************************/
void free_mol_list(struct model_pak *data)
{
GSList *list;
struct mol_pak *mol;

for (list=data->moles ; list ; list=g_slist_next(list))
  {
  mol = list->data;

  g_slist_free(mol->cores);
  g_free(mol);
  }

g_slist_free(data->moles);
data->moles=NULL;
}

/******************************/
/* compute molecule centroids */
/******************************/
#define DEBUG_CALC_MOL_CENT 0
void connect_centroid_compute(struct mol_pak *mol)
{
gdouble scale, n;
GSList *list;
struct core_pak *core;

/* centroid calc & core->mol reference */
VEC3SET(mol->centroid, 0.0, 0.0, 0.0);
n=0.0;
for (list=mol->cores ; list ; list=g_slist_next(list))
  {
  core = list->data;

#if DEBUG_CALC_MOL_CENT
printf("[%4s]  [%10.4f %10.4f %10.4f]\n", core->atom_label, core->x[0], core->x[1], core->x[2]);
#endif

  ARR3ADD(mol->centroid, core->x);
  n += 1.0;
  }

if (n)
  {
  scale = 1.0 / n;
  VEC3MUL(mol->centroid, scale);
  }

#if DEBUG_CALC_MOL_CENT
P3VEC(" centroid: ", mol->centroid);
#endif
}

/*********************************/
/* sort molecule cores primitive */
/*********************************/
gint sort_core_order(gpointer ptr_core1, gpointer ptr_core2)
{
struct core_pak *core1 = ptr_core1, *core2 = ptr_core2;

if (core1->atom_order < core2->atom_order)
  return(-1);
return(1);
}

/*********************************************************************/
/* sort cores in molecules into the same order as the main core list */
/*********************************************************************/
void sort_mol_cores(struct model_pak *model)
{
GSList *list;
struct mol_pak *mol;

for (list=model->moles ; list ; list=g_slist_next(list))
  {
  mol = list->data;

  mol->cores = g_slist_sort(mol->cores, (gpointer) sort_core_order);
  }
}

/**************************/
/* compute molecule lists */
/**************************/
#define DEBUG_CALC_MOLS 0
void connect_molecules(struct model_pak *data)
{
gint n=0;
guint i=0;
gdouble x[3];
GSList *list, *blist, *mlist, *list1, *list2;
struct core_pak *core, *core2;
struct bond_pak *bond;
struct mol_pak *mol;

g_assert(data != NULL);

/* clean */
free_mol_list(data);

/* init */
for (list=data->cores ; list ; list=g_slist_next(list))
  {
  core = list->data;
  core->status &= ~PRUNED;
/* NEW - record core order in list */
  core->atom_order = i++;
  }

/* build mol lists */
mlist=NULL;
for (list=data->cores ; list ; list=g_slist_next(list))
  {
  core = list->data;
  if (core->status & PRUNED)
    continue;

/* init a new mol list */
  if (!mlist)
    {
    mlist = g_slist_append(mlist, core);
    core->status |= PRUNED;
    core->molecule = n;
    }

/* while cores are being added to mlist */
  do
    {
    blist = NULL;
    for (list2=mlist ; list2 ; list2=g_slist_next(list2))
      {
      core = list2->data;

/* get list of atoms bonded to the current core */
      for (list1=core->bonds ; list1 ; list1=g_slist_next(list1))
        {
        bond = list1->data;

/* don't include these types in connectivity tracing */
        if (bond->status == DELETED)
          continue;
        if (bond->type == BOND_HBOND)
          continue;

/* get the attached core */
        ARR3SET(x, bond->offset);
        if (core == bond->atom1)
          {
          core2 = bond->atom2;
          VEC3MUL(x, 2.0);
          }
        else
          {
          core2 = bond->atom1;
          VEC3MUL(x, -2.0);
          }

/* only accumulate bonded, non-pruned atoms */
        if (core2->status & PRUNED)
          continue;

        blist = g_slist_prepend(blist, core2); 
        core2->status |= PRUNED;

/* always attempt to unfragment molecules */
/* CURRENT - fix for fractional_min not giving cartesian minimum */
        ARR3ADD(x, core->x);
        ARR3SUB(x, core2->x);
        ARR3ADD(core2->x, x);
        if (core2->shell)
          {
          ARR3ADD((core2->shell)->x, x);
          }

/* deprec. */
        core2->molecule = n;
        }
      }
/* augment the current listing */
    if (blist)
      mlist = g_slist_concat(mlist, blist);
    }
  while (g_slist_length(blist));
  n++;

/* done - completed mol */
  mol = g_malloc(sizeof(struct mol_pak));
  mol->cores = mlist;
  mlist=NULL;
  data->moles = g_slist_prepend(data->moles, mol);

/* label the cores as belonging to the mol */
  for (list1=mol->cores ; list1 ; list1=g_slist_next(list1))
    {
    core = list1->data;
    core->mol = mol;
    }
  }

#if DEBUG_CALC_MOLS
dump_mols(data);
#endif

/* NEW - fine grained molecule split check */
for (list=data->moles ; list ; list=g_slist_next(list))
  {
  mol = list->data;
  if (connect_split_list(mol->cores))
    {
/* confine/calc centroid */
    coords_confine_cores(mol->cores, data);
    connect_centroid_compute(mol);
    }
  else
    {
/* calc/confine centroid */
    connect_centroid_compute(mol);
    coords_confine_centroid(mol, data);
    }
  }
coords_compute(data);
model_content_refresh(data); /* Added by C.Fisher 2005 */

/* NEW - ensure core ordering in molecules is the same as in the main list */
sort_mol_cores(data);
}

/****************************************/
/* initialize model fragment operations */
/****************************************/
void connect_fragment_init(struct model_pak *model)
{
GSList *list;
struct core_pak *core;

/* clear the pruned flag */
for (list=model->cores ; list ; list=g_slist_next(list))
  {
  core = list->data;
  core->status &= ~PRUNED;
  }
}

/*************************************/
/* acquire a fragment from two atoms */
/*************************************/
GSList *connect_fragment_get(struct core_pak *c1, struct core_pak *c2, struct model_pak *model)
{
GSList *list, *item1, *item2;
GSList *branch_new, *branch_all;
struct core_pak *core1, *core2;

/* get initial branches */
branch_new = exit_branches(c1, c2);
item1 = branch_all = g_slist_prepend(branch_new, c2);

#if DEBUG_ADD_FRAGMENT
printf("seeking branches for: %p -> %p : ", c1, c2);
printf("initial branches: %d\n", g_slist_length(branch_new));
for (item1=branch_new ; item1 ; item1 = g_slist_next(item1))
  {
  core1 = item1->data;
  printf(" - %p (%s)\n", core1, core1->atom_label);
  }
#endif

/* iterate through upper branch level */
while (item1)
  {
  core1 = item1->data;
  item1 = g_slist_next(item1);
  g_assert(core1 != NULL);

/* iterate through lower branch level and aquire next level of branches */
  item2 = branch_new;
  branch_new = NULL;
  while (item2)
    {
    core2 = item2->data;
    g_assert(core2 != NULL);
    item2 = g_slist_next(item2);
    if (core2->status & PRUNED)
      continue;

/* prune this core, so its branches are not re-examined */
    core2->status |= PRUNED;
    list = exit_branches(core1, core2);

#if DEBUG_ADD_FRAGMENT
{
struct core_pak *core3;
GSList *item3;

printf("branches: %p (%s) -> %p (%s)", core1, core1->label, core2, core2->label);
printf(" : found %d.\n", g_slist_length(list));
for (item3=list ; item3 ; item3 = g_slist_next(item3))
  {
  core3 = item3->data;
  printf(" - %p (%s)\n", core3, core3->atom_label);
  }
}
#endif

/* aquire all branches at current level */
    branch_new = g_slist_concat(branch_new, list);
    }
/* test if we arrive back at first core (from a different direction) */
  if (g_slist_find(branch_new, c1))
    {
    printf("Bad fragment\n");
    g_slist_free(branch_new);
    g_slist_free(branch_all);
    return(NULL);
    }
/* accumlate new branches into main list */
  branch_all = g_slist_concat(branch_all, branch_new);
  }
return(branch_all);
}

/***************************************/
/* connectivity update for given model */
/***************************************/
void connect_refresh(struct model_pak *model)
{
connect_bonds(model);
connect_molecules(model);
model_content_refresh(model); /* Added by C.Fisher 2005 */
}

/********************************************/
/* complete connectivity refresh and redraw */
/********************************************/
void connect_refresh_global(void)
{
struct model_pak *model;

model = sysenv.active_model;

coords_compute(model);
connect_bonds(model);
connect_molecules(model);

redraw_canvas(SINGLE);
model_content_refresh(model); /* Added by C.Fisher 2005 */
}

