/*
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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <math.h>

#include "gdis.h"
#include "file.h"
#include "coords.h"
#include "matrix.h"
#include "measure.h"
#include "numeric.h"
#include "zmatrix.h"
#include "zmatrix_pak.h"
#include "parse.h"
#include "interface.h"

/* top level data structure */
extern struct sysenv_pak sysenv;

/******************************/
/* output zmatrix coordinates */
/******************************/
/* NB - this is siesta specific */
void zmat_mol_print(FILE *fp, gpointer data)
{
gint i;
GSList *list;
struct zmat_pak *zmat=data;
struct zval_pak *zval;

g_assert(data != NULL);

if (g_slist_length(zmat->zlines))
  {
  if (zmat->fractional)
    fprintf(fp, "molecule fractional\n");
  else
    fprintf(fp, "molecule\n");

  for (list=zmat->zlines ; list ; list=g_slist_next(list))
    {
    zval = list->data;

    fprintf(fp, " %d  %d %d %d", zval->type, zval->connect[0], zval->connect[1], zval->connect[2]);

    for (i=0 ; i<3 ; i++)
      {
      if (zval->name[i])
        fprintf(fp, "  %s", (gchar *) zval->name[i]);
      else
        fprintf(fp, "  %f", zval->value[i]);
      }

    fprintf(fp, "  %d %d %d", zval->fitting[0], zval->fitting[1], zval->fitting[2]);

    fprintf(fp, "\n");
    }
  }
}

/******************************/
/* primtive key building list */
/******************************/
void zmat_build_keys(gpointer key, gpointer val, gpointer data)
{
GList **list = data;

*list = g_list_append(*list, key);
}

/****************************/
/* output zmatrix variables */
/****************************/
/* NB - this is siesta specific */
void zmat_var_print(FILE *fp, gpointer data)
{
gchar *value;
GList *keys=NULL, *list;
struct zmat_pak *zmat=data;

g_assert(data != NULL);

g_hash_table_foreach(zmat->vars, &zmat_build_keys, &keys);

if (g_list_length(keys))
  {
  fprintf(fp, "variables\n");

  for (list=keys ; list ; list=g_list_next(list))
    {
    value = g_hash_table_lookup(zmat->vars, list->data);

    fprintf(fp, " %s  %s\n", (gchar *) list->data, value);
    }
  }
}

/****************************/
/* output zmatrix constants */
/****************************/
/* NB - this is siesta specific */
void zmat_const_print(FILE *fp, gpointer data)
{
gchar *value;
GList *keys=NULL, *list;
struct zmat_pak *zmat=data;

g_assert(data != NULL);

g_hash_table_foreach(zmat->consts, &zmat_build_keys, &keys);

if (g_list_length(keys))
  {
  fprintf(fp, "constants\n");

  for (list=keys ; list ; list=g_list_next(list))
    {
    value = g_hash_table_lookup(zmat->consts, list->data);

    fprintf(fp, " %s  %s\n", (gchar *) list->data, value);
    }
  }
}

/***********************************/
/* output non molecule coordinates */
/***********************************/
/* NB - this is siesta specific */
void zmat_coord_print(FILE *fp, GSList *species_list, gpointer data)
{
gint i, j;
gdouble x[3];
GSList *list, *clist, *slist;
struct model_pak *model=data;
struct species_pak *species_data;
struct core_pak *core;
struct zmat_pak *zmat;

g_assert(data != NULL);

clist = g_slist_copy(model->cores);

zmat = model->zmatrix;

g_assert(zmat != NULL);

for (list=zmat->zcores ; list ; list=g_slist_next(list))
  clist = g_slist_remove(clist, list->data);

j = 1 + g_slist_length(zmat->zcores);

if (clist)
  {
  if (model->fractional)
    fprintf(fp, "fractional\n");
  else
    fprintf(fp, "cartesian\n");

  for (list=clist ; list ; list=g_slist_next(list))
    {
    core = list->data;

    ARR3SET(x, core->x);

/* NB: want fractional if 3D periodic, otherwise cartesian */
    if (!model->fractional)
      vecmat(model->latmat, x);

/* find corresponding number of this element in the species block */
    i=1;
    for (slist=species_list ; slist ; slist=g_slist_next(slist))
      {
      species_data = slist->data;

      if (g_ascii_strcasecmp(core->atom_label, species_data->label) == 0)
        break;

      i++;
      }
    fprintf(fp," %d  %14.9f  %14.9f  %14.9f    1 1 1   %d\n", i, x[0], x[1], x[2], j++);
    }
  }

g_slist_free(clist);

}

/******************************/
/* zmatrix debugging printout */
/******************************/
void zmat_debug(gpointer data)
{
zmat_mol_print(stdout, data);
zmat_var_print(stdout, data);
}

/**********************/
/* creation primitive */
/**********************/
gpointer zmat_new(void)
{
struct zmat_pak *zmat;

zmat = g_malloc(sizeof(struct zmat_pak));

zmat->fractional = FALSE;
zmat->distance_scale = 1.0;
zmat->angle_scale = 1.0;
zmat->distance_units = g_strdup("Ang");
zmat->angle_units = g_strdup("Rad");
zmat->zlines = NULL;
zmat->zcores = NULL;
zmat->vars = g_hash_table_new(g_str_hash, g_str_equal);
zmat->consts = g_hash_table_new(g_str_hash, g_str_equal);

return(zmat);
}

/*************************/
/* destruction primitive */
/*************************/
void zmat_free(gpointer data)
{
struct zmat_pak *zmat = data;

g_assert(zmat != NULL);

g_free(zmat->distance_units);
g_free(zmat->angle_units);

free_slist(zmat->zlines);
/* NB: zcores is a subset of the model's main core list - don't free its data */
g_slist_free(zmat->zcores);
g_hash_table_destroy(zmat->vars);
g_hash_table_destroy(zmat->consts);
g_free(zmat);
}

/**********************/
/* creation primitive */
/**********************/
gpointer zval_new(void)
{
struct zval_pak *zval;

zval = g_malloc(sizeof(struct zval_pak));

zval->type = -1;

zval->connect[0] = -1;
zval->connect[1] = -1;
zval->connect[2] = -1;

zval->fitting[0] = -1;
zval->fitting[1] = -1;
zval->fitting[2] = -1;

zval->name[0] = NULL;
zval->name[1] = NULL;
zval->name[2] = NULL;

VEC3SET(zval->value, 0.0, 0.0, 0.0);

return(zval);
}

/**********************************/
/* zmatrix distance units/scaling */
/**********************************/
void zmat_distance_units_set(gpointer data, gint type)
{
struct zmat_pak *zmat = data;

g_assert(zmat != NULL);

g_free(zmat->distance_units);

switch (type)
  {
  case ANGSTROM:
    zmat->distance_scale = 1.0;
    zmat->distance_units = g_strdup("Ang");
    break;
  case BOHR:
    zmat->distance_scale = AU2ANG;
    zmat->distance_units = g_strdup("Bohr");
    break;
  }
}

/***********************/
/* retrieval primitive */
/***********************/
const gchar *zmat_distance_units_get(gpointer data)
{
struct zmat_pak *zmat = data;

g_assert(zmat != NULL);

return(zmat->distance_units);
}

/*******************************/
/* zmatrix angle units/scaling */
/*******************************/
void zmat_angle_units_set(gpointer data, gint type)
{
struct zmat_pak *zmat = data;

g_assert(zmat != NULL);

g_free(zmat->angle_units);

switch (type)
  {
  case DEGREES:
    zmat->angle_scale = D2R;
    zmat->angle_units = g_strdup("Deg");
    break;
  case RADIANS:
    zmat->angle_scale = 1.0;
    zmat->angle_units = g_strdup("Rad");
    break;
  }
}

/***********************/
/* retrieval primitive */
/***********************/
const gchar *zmat_angle_units_get(gpointer data)
{
struct zmat_pak *zmat = data;

g_assert(zmat != NULL);

return(zmat->angle_units);
}

/*********************************************/
/* set zmatrix molecule coords as cartesians */
/*********************************************/
void zmat_cartesian_set(gpointer data)
{
struct zmat_pak *zmat = data;

g_assert(zmat != NULL);

zmat->fractional=FALSE;
}

/**********************************************/
/* set zmatrix molecule coords as fractionals */
/**********************************************/
void zmat_fractional_set(gpointer data)
{
struct zmat_pak *zmat = data;

g_assert(zmat != NULL);

zmat->fractional=TRUE;
}

/***********************************/
/* get the number of zmatrix lines */
/***********************************/
gint zmat_entries_get(gpointer data)
{
struct zmat_pak *zmat = data;

g_assert(zmat != NULL);

return(g_slist_length(zmat->zlines));
}

/*********************************/
/* add a zmatrix coordinate line */
/*********************************/
void zmat_core_add(const gchar *line, gpointer data)
{
gint i, j, num_tokens;
gchar **buff;
struct zmat_pak *zmat = data;
struct zval_pak *zval;

buff = tokenize(line, &num_tokens);

g_assert(num_tokens > 8);

zval = zval_new(); 
zmat->zlines = g_slist_append(zmat->zlines, zval);

/* 1st column atom label */
zval->type = str_to_float(*buff);

/* zmatrix connectivity data */
for (i=1 ; i<4 ; i++)
  zval->connect[i-1] = str_to_float(*(buff+i));

/* zmatrix coordinate data */
j=0;
for (i=4 ; i<7 ; i++)
  {
/* if it's a number -> value, else stored as variable name */
  if (str_is_float(*(buff+i)))
    zval->value[j] = str_to_float(*(buff+i));
  else
    zval->name[j] = g_strdup(*(buff+i));
  j++;
  }

/* zmatrix fitting flags */
for (i=7 ; i<10 ; i++)
  zval->fitting[i-7] = str_to_float(*(buff+i));

g_strfreev(buff);
}

/***********************/
/* add a variable line */
/***********************/
void zmat_var_add(const gchar *line, gpointer data)
{
gint num_tokens;
gchar **buff;
struct zmat_pak *zmat = data;

buff = tokenize(line, &num_tokens);

g_assert(num_tokens > 1);

g_hash_table_insert(zmat->vars, g_strdup(*buff), g_strdup(*(buff+1)));

g_strfreev(buff);
}

/***********************/
/* add a constant line */
/***********************/
void zmat_const_add(const gchar *line, gpointer data)
{
gint num_tokens;
gchar **buff;
struct zmat_pak *zmat = data;

buff = tokenize(line, &num_tokens);

g_assert(num_tokens > 1);

g_hash_table_insert(zmat->consts, g_strdup(*buff), g_strdup(*(buff+1)));

g_strfreev(buff);
}

/**********************************************/
/* convert silly SIESTA numbers to atom types */
/**********************************************/
void zmat_type(gpointer data, GSList *species)
{
GSList *list;
struct species_pak *species_data;
struct zmat_pak *zmat = data;
struct zval_pak *zval;

/* checks */
if (!zmat)
  return;
if (!species)
  return;

for (list=zmat->zlines ; list ; list=g_slist_next(list))
  {
  zval = list->data;

  species_data = g_slist_nth_data(species, zval->type - 1);

  if (species_data)
    zval->elem = g_strdup(species_data->label);
  else
    zval->elem = g_strdup("x");
  }
}

/*****************************************************************/
/* look up constant and variable list for numerical substitution */
/*****************************************************************/
gdouble zmat_table_lookup(const gchar *name, gpointer data)
{
const gchar *value;
struct zmat_pak *zmat = data;

g_assert(data != NULL);

value = g_hash_table_lookup(zmat->vars, name);
if (!value)
  value = g_hash_table_lookup(zmat->consts, name);

if (value)
  return(str_to_float(value));

return(0.0);
}

/**********************************************/
/* create a model core list from zmatrix data */
/**********************************************/
#define ZMAT_PROCESS_DEBUG 0
void zmat_process(gpointer data, struct model_pak *model)
{
gint i, n;
gdouble r, a, d;
gdouble v1[3], v2[3], v3[3], m1[9], m2[9];
gpointer tmp;
GSList *list;
struct zmat_pak *zmat = data;
struct zval_pak *zval;
struct core_pak *core[4] = {NULL, NULL, NULL, NULL};

/* checks */
if (!zmat)
  return;

matrix_lattice_init(model);

#if ZMAT_PROCESS_DEBUG
printf("distance scale = %f\n", zmat->distance_scale);
printf("angle scale = %f\n", zmat->angle_scale);
#endif

/* track the core # - 1st 3 items in zmatrix are exceptions */
n = 0;
for (list=zmat->zlines ; list ; list=g_slist_next(list))
  {
  zval = list->data;

/* check for variable names */
  for (i=3 ; i-- ; )
    {
    if (zval->name[i])
      {
/* hash table lookup for variable value */
      zval->value[i] = zmat_table_lookup(zval->name[i], zmat);
      }
    }

/* create the core */
#if ZMAT_PROCESS_DEBUG
printf("[%d = %s] [%d %d %d]", zval->type, zval->elem,
                               zval->connect[0], zval->connect[1], zval->connect[2]);
P3VEC(" x: ", zval->value);

#endif

/* TODO - need to mark zmatrix generated cores as special */
/* probably have another list in the zmat struct that contains */
/* all the cores created below */

  switch (n)
    {
    case 0:
/* TODO - convert to cartesian if fractional and enforce cart model */
      core[0] = new_core(zval->elem, model);
      model->cores = g_slist_prepend(model->cores, core[0]);
      zmat->zcores = g_slist_prepend(zmat->zcores, core[0]);
      ARR3SET(core[0]->x, zval->value);
      if (zmat->fractional)
        vecmat(model->latmat, core[0]->x);
      else
        {
        VEC3MUL(core[0]->x, zmat->distance_scale);
        }
      break;

    case 1:
      core[1] = new_core(zval->elem, model);
      model->cores = g_slist_prepend(model->cores, core[1]);
      zmat->zcores = g_slist_prepend(zmat->zcores, core[1]);

      r = zmat->distance_scale * zval->value[0];
/*
      a = zmat->angle_scale * zval->value[1];
      d = zmat->angle_scale * zval->value[2];
*/

/* SIESTA hack : z-axis angle is 2nd, and theta is 3rd (last) */
      a = zmat->angle_scale * zval->value[2];
      d = zmat->angle_scale * zval->value[1];

      v1[0] = v1[1] = r * sin(d);
      v1[0] *= cos(a);
      v1[1] *= sin(a);
      v1[2] = r * cos(d);

      ARR3SET(core[1]->x, core[0]->x);
      ARR3ADD(core[1]->x, v1);
      break;

    case 2:
/* check the connection order */
      if (zval->connect[0] == 2)
        {
        tmp = core[0];
        core[0] = core[1];
        core[1] = tmp;
        }

      r = zmat->distance_scale * zval->value[0];
      a = zmat->angle_scale * zval->value[1];
      d = zmat->angle_scale * zval->value[2];

      ARR3SET(v2, core[1]->x);
      ARR3SUB(v2, core[0]->x);

/* get rotation axis for bond angle */
      VEC3SET(v3, 0.0, 0.0, 1.0);
      crossprod(v1, v3, v2);

/* rotate bondlength scaled vector into position */
      matrix_v_rotation(m2, v1, a);
      ARR3SET(v3, v2);
      normalize(v3, 3);
      VEC3MUL(v3, r);
      vecmat(m2, v3);

/* rotate to give required dihedral */
      matrix_v_rotation(m1, v2, d);
      vecmat(m1, v3);

/* generate the atom position */
      core[2] = new_core(zval->elem, model);
      model->cores = g_slist_prepend(model->cores, core[2]);
      zmat->zcores = g_slist_prepend(zmat->zcores, core[2]);

      ARR3SET(core[2]->x, core[0]->x);
      ARR3ADD(core[2]->x, v3);
      break;

    default:
/* get core connectivity (NB: prepending cores - hence n - number) */
      core[0] = g_slist_nth_data(zmat->zcores, n-zval->connect[0]); 
      core[1] = g_slist_nth_data(zmat->zcores, n-zval->connect[1]); 
      core[2] = g_slist_nth_data(zmat->zcores, n-zval->connect[2]); 
g_assert(core[0] != NULL);
g_assert(core[1] != NULL);
g_assert(core[2] != NULL);

      r = zmat->distance_scale * zval->value[0];
      a = zmat->angle_scale * zval->value[1];
      d = zmat->angle_scale * zval->value[2];

/* setup vectors */
      ARR3SET(v2, core[1]->x);
      ARR3SUB(v2, core[0]->x);
      ARR3SET(v3, core[2]->x);
      ARR3SUB(v3, core[1]->x);

/* rotate v3 about v2 to give dihedral */
      matrix_v_rotation(m1, v2, d);
      vecmat(m1, v3);

/* get rotation axis and matrix for bond angle */
      crossprod(v1, v3, v2);
      matrix_v_rotation(m2, v1, a);
      normalize(v2, 3);
      VEC3MUL(v2, r);
      vecmat(m2, v2);

/* generate the atom position */
      core[3] = new_core(zval->elem, model);
      model->cores = g_slist_prepend(model->cores, core[3]);
      zmat->zcores = g_slist_prepend(zmat->zcores, core[3]);

      ARR3SET(core[3]->x, core[0]->x);
      ARR3ADD(core[3]->x, v2);

/* TODO (maybe) - some zmatrix constructions implicitly assume */
/* some checking for duplicate atoms & reversing the torsional */
/* angle sense to accomodate this. */

      break;
    }
  n++;
  }

/* zmatrix cores are created in cartesians - revert to fractional if required */
if (model->fractional)
  {
  for (list=zmat->zcores ; list ; list=g_slist_next(list))
    {
    core[0] = list->data;
    vecmat(model->ilatmat, core[0]->x);
    }
  }
}

/*********************************/
/* check if two cores are bonded */
/*********************************/
gint zmat_bond_check(struct core_pak *c1, struct core_pak *c2)
{
GSList *list;
struct bond_pak *bond;

for (list=c1->bonds ; list ; list=g_slist_next(list))
  {
  bond = list->data;
  if (bond->type == BOND_HBOND)
    continue;
  if (bond->atom1 == c2 || bond->atom2 == c2)
    return(TRUE);
  }
return(FALSE);
}

/***************************************************/
/* scan past zmatrix connectivity for valid chains */
/***************************************************/
#define DEBUG_ZMAT_CONNECT_FIND 0
gint zmat_connect_find(gint n, struct core_pak **core, struct zmat_pak *zmatrix)
{
gint m;
GSList *list1, *list2;
struct zval_pak *zval;
struct core_pak *seek;
struct bond_pak *bond;

g_assert(core != NULL);
g_assert(zmatrix != NULL);

for (list1=core[0]->bonds ; list1 ; list1=g_slist_next(list1))
  {
  bond = list1->data;

  if (bond->type == BOND_HBOND)
    continue;

  if (bond->atom1 == core[0])
    seek = bond->atom2;
  else
    seek = bond->atom1;

/* TODO - if atom has been PRUNED - should be gtg (faster than index check) */
  m = g_slist_index(zmatrix->zcores, seek);

#if DEBUG_ZMAT_CONNECT_FIND
printf("atom: %d, seeking chain with head: [%d] %p\n", n, m, seek);
#endif

  if (m < 0)
    {
#if DEBUG_ZMAT_CONNECT_FIND
printf("Ignoring atom that is not yet defined.\n");
#endif
    continue;
    }

/* trivial case - use the head atom's zmatrix line (IF it's not one */
/* of the first 2 special case atoms */
  if (m < n && m > 3)
    {
    zval = g_slist_nth_data(zmatrix->zlines, m);

    g_assert(zval != NULL);

    core[1] = g_slist_nth_data(zmatrix->zcores, m);
    core[2] = g_slist_nth_data(zmatrix->zcores, zval->connect[0]-1);
    core[3] = g_slist_nth_data(zmatrix->zcores, zval->connect[1]-1);

#if DEBUG_ZMAT_CONNECT_FIND
printf("Using trivial case: [%d][%d][%d]\n", 
g_slist_index(zmatrix->zcores, core[1]),
g_slist_index(zmatrix->zcores, core[2]),
g_slist_index(zmatrix->zcores, core[3]));
#endif

    return(TRUE);
    }


  list2 = zmatrix->zlines;
  list2 = g_slist_next(list2);
  list2 = g_slist_next(list2);

  while (list2)
    {
    zval = list2->data;

    if (m == zval->connect[0])
      {
#if DEBUG_ZMAT_CONNECT_FIND
printf("found: [%d][%d][%d]\n", zval->connect[0], zval->connect[1], zval->connect[2]);
#endif

      core[1] = g_slist_nth_data(zmatrix->zcores, zval->connect[0]);
      core[2] = g_slist_nth_data(zmatrix->zcores, zval->connect[1]);
      core[3] = g_slist_nth_data(zmatrix->zcores, zval->connect[2]);
      return(TRUE);
      }
    if (m == zval->connect[2])
      {
#if DEBUG_ZMAT_CONNECT_FIND
printf("found: [%d][%d][%d]\n", zval->connect[0], zval->connect[1], zval->connect[2]);
#endif

      core[1] = g_slist_nth_data(zmatrix->zcores, zval->connect[2]);
      core[2] = g_slist_nth_data(zmatrix->zcores, zval->connect[1]);
      core[3] = g_slist_nth_data(zmatrix->zcores, zval->connect[0]);
      return(TRUE);
      }
    list2 = g_slist_next(list2);
    }
  }
return(FALSE);
}

/**********************************************************************/
/* get the next core bonded to input, that itself is maximally bonded */
/**********************************************************************/
gpointer zmat_connect_next(struct core_pak *core)
{
gint size, max;
GSList *list;
struct bond_pak *bond;
struct core_pak *tmp, *next;

max = 0;
next = NULL;
for (list=core->bonds ; list ; list=g_slist_next(list))
  {
  bond = list->data;

  if (bond->type == BOND_HBOND)
    continue;

  if (bond->atom1 == core)
    tmp = bond->atom2;
  else
    tmp = bond->atom1;

  if (tmp->status & PRUNED)
    continue;

  size = g_slist_length(tmp->bonds);
  if (size > max)
    {
    max = size;
    next = tmp;
    }
  }

return(next);
}

/*****************************************************/
/* sort atoms follow zmatrix compatible connectivity */
/*****************************************************/
GSList *zmat_connect_sort(GSList *molecule)
{
GSList *list, *path, *sorted;
struct core_pak *core, *next;

/* flag all atoms as unpruned */
for (list=molecule ; list ; list=g_slist_next(list))
  {
  core = list->data;
  core->status &= ~PRUNED;
  }

/* CURRENT - try to start connectivity chain with an atom with 1 bond */
next = NULL;
path = NULL;
for (list=molecule ; list ; list=g_slist_next(list))
  {
  core = list->data;

  if (g_slist_length(core->bonds) == 1)
    {
    next = core;
    break;
    }
  }

/* if we can't find single bonded atom - start with 1st in selection & hope */
if (next)
  path = g_slist_append(path, core);
else
  path = g_slist_append(path, molecule->data);

/* trace connectivity - store atoms with more than two bonds */
sorted = NULL;
while (path)
  {
  core = path->data;

  if (!(core->status & PRUNED))
    {
    core->status |= PRUNED;
    sorted = g_slist_append(sorted, core);
    }

  while (core)
    {
    next = zmat_connect_next(core);

    if (next)
      {
      if (g_slist_length(next->bonds) > 2)
        path = g_slist_append(path, next);
 
      next->status |= PRUNED;
      sorted = g_slist_append(sorted, next);
      }

    core = next;
    }

/* current chain exhausted - check current node (path) for another to */
/* expore, if nothing - try to get next node */
   while (path) 
     {
     core = zmat_connect_next(path->data);
     if (core)
       break;
 
     path = g_slist_next(path);
     }
  }

return(sorted);
}

/**************************************************/
/* construct a zmatrix for an input list of atoms */
/**************************************************/
#define DEBUG_ZMAT_BUILD 0
void zmat_build(void)
{
gint i, j, k, n, type;
gdouble r, a, d, x[4][3], v[3];
gdouble  zaxis[3] = {0.0, 0.0, 1.0};
gchar *line;
GSList *list, *species;
struct zmat_pak *zmat;
struct core_pak *core[4] = {NULL, NULL, NULL, NULL};
struct model_pak *model;

model = sysenv.active_model;
if (!model)
  return;

/* CURRENT - using selection as our list of cores to generate a zmatrix from */
if (!model->selection)
  {
  gui_text_show(WARNING, "ZMATRIX: please select a molecule.\n");
  return;
  }

/* destroy old zmatrix */
/* TODO - prompt if non null */
zmat_free(model->zmatrix);
zmat = model->zmatrix = zmat_new();
zmat_angle_units_set(model->zmatrix, DEGREES);

/* setup SIESTA species type */
species = fdf_species_build(model);

/* sort the list so it follows molecular connectivity */
model->selection = zmat_connect_sort(model->selection);

n=0;
for (list=model->selection ; list ; list=g_slist_next(list))
  {
/* current atom/zmatrix line init */
  core[0] = list->data;
  type = fdf_species_index(core[0]->atom_label, species);
  line = NULL;

  zmat->zcores = g_slist_append(zmat->zcores, core[0]);

/* build a ZMATRIX line for processing */
  switch (n)
    {
    case 0:
      if (core[0])
        {
        ARR3SET(x[0], core[0]->x);
        vecmat(model->latmat, x[0]);
        }
      line = g_strdup_printf("%d  0 0 0  %f %f %f  0 0 0\n", type, x[0][0], x[0][1], x[0][2]);
      break;

    case 1:
      if (core[0])
        {
        ARR3SET(x[0], core[0]->x);
        vecmat(model->latmat, x[0]);
        }
      if (core[1])
        {
        ARR3SET(x[1], core[1]->x);
        vecmat(model->latmat, x[1]);
        }

      r = measure_distance(x[0], x[1]);

/* angle with z axis */
      ARR3SET(v, x[0]);
      ARR3SUB(v, x[1]);
      a = R2D * via(v, zaxis, 3);

/* angle between xy projection and x axis */
      d = R2D * angle_x_compute(v[0], v[1]);

      line = g_strdup_printf("%d  1 0 0  %f %f %f 0 0 0\n", type, r, a, d);
      break;

    case 2:
/* coords init */
  for (i=3 ; i-- ; )
    {
    if (core[i])
      {
      ARR3SET(x[i], core[i]->x);
      vecmat(model->latmat, x[i]);
      }
    else
      g_assert_not_reached();
    }

      r = measure_distance(x[0], x[1]);
      a = measure_angle(x[0], x[1], x[2]);

/* create a fake core -> 1 unit displaced in the z direction */
      g_assert(core[3] == NULL);
      core[3] = core_new("x", NULL, model);
      ARR3SET(core[3]->rx, core[2]->rx);
      ARR3ADD(core[3]->rx, zaxis); 
      d = measure_torsion(core);
      core_free(core[3]);

      line = g_strdup_printf("%d  2 1 0  %f %f %f 0 0 0\n", type,r,a,d);
      break;

    default:

/* connectivity test */
      if (!zmat_bond_check(core[0], core[1]))
        {
#if DEBUG_ZMAT_BUILD
printf("[%d] non-connected atoms [%f]\n", n, measure_distance(x[0], x[1]));
#endif
/* need to build a new connectivity chain starting from core[0] */
        core[1] = core[2] = core[3] = NULL;
        if (!zmat_connect_find(n, core, zmat))
          {
          gui_text_show(WARNING, "ZMATRIX: bad connectivity (molecule will be incomplete)\n");
          goto zmat_build_done;
          }
        }

/* coords init */
      for (i=3 ; i-- ; )
        {
        if (core[i])
          {
          ARR3SET(x[i], core[i]->x);
          vecmat(model->latmat, x[i]);
          }
        else
          g_assert_not_reached();
        }

      r = measure_distance(x[0], x[1]);
      a = measure_angle(x[0], x[1], x[2]);
      d = measure_torsion(core);

/* NB: indexing starts from 0, siesta starts from 1 (naturally) */
      i = 1+g_slist_index(zmat->zcores, core[1]);
      j = 1+g_slist_index(zmat->zcores, core[2]);
      k = 1+g_slist_index(zmat->zcores, core[3]);

      line = g_strdup_printf("%d  %d %d %d  %f %f %f 0 0 0\n", type,i,j,k,r,a,d);
    }

/* process a successfully constructed ZMATRIX line */
  if (line)
    {
    zmat_core_add(line, model->zmatrix);
    g_free(line);
    }

/* shuffle */
  core[3] = core[2];
  core[2] = core[1];
  core[1] = core[0];

  n++;
  }

zmat_build_done:

/* do the species typing */
zmat_type(model->zmatrix, species);

free_slist(species);
}

