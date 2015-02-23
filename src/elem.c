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
#include <strings.h>
#include <stdlib.h>
#include <ctype.h>

#include "gdis.h"
#include "coords.h"
#include "edit.h"
#include "file.h"
#include "parse.h"
#include "matrix.h"
#include "numeric.h"
#include "surface.h"
#include "select.h"
#include "interface.h"
#include "gui_shorts.h"
#include "opengl.h"

/* main structures */
extern struct sysenv_pak sysenv;
extern struct elem_pak elements[];

/* TODO - XML & add diffraction data to elements file */
/**************************************/
/* parse a file for gdis element data */
/**************************************/
/* type - should the element data be treated as default or patched */
/* ie the standard library, or modified values from a gdisrc file */
gint read_elem_data(FILE *fp, gint type)
{
gint i, n, key, flag, num_tokens; 
gint *set;
gchar **buff, *str, *symbol, line[LINELEN];
gdouble *sfc;
struct elem_pak elem, elem_default;
GSList *list;

/* checks */
g_return_val_if_fail(fp != NULL, 1);

/* FIXME - num_keys() is returning a value that is too small */
/*
n = num_keys();
*/
n = 1000;
set = g_malloc0(n * sizeof(gint));

while(!fgetline(fp,line))
  {
  list = get_keywords(line);
  if (list) 
    {
    switch(GPOINTER_TO_INT(list->data))
      {
      case GDIS_ELEM_START:
/* init for element */
         elem.shell_charge = 0.0;
/* init the set flags */
        for (i=n ; i-- ; )
          set[i] = 0;
/* keep looping until end of element data, or EOF */
        flag=0;
        while(!flag)
          {
/* get next line & process for keywords */
          if (fgetline(fp, line))
            break;
          list = get_keywords(line);
          buff = tokenize(line, &num_tokens);

/* process 1st keyword, and subsequent data (if any) */
          if (list)
            {
            key = GPOINTER_TO_INT(list->data);
g_assert(key < n);
            switch (key)
              {
/* these shouldn't be modified if not default */
              case SYMBOL:
                g_assert(num_tokens > 1);
                elem.symbol = g_strdup(*(buff+1));
                break; 
              case NAME:
                g_assert(num_tokens > 1);
                elem.name = g_strdup(*(buff+1));
                break; 
              case NUMBER:
                g_assert(num_tokens > 1);
                elem.number = (gint) str_to_float(*(buff+1));
                set[key]++;
                break; 
              case WEIGHT:
                g_assert(num_tokens > 1);
                elem.weight = str_to_float(*(buff+1));
                set[key]++;
                break; 
              case COVALENT:
                g_assert(num_tokens > 1);
                elem.cova = str_to_float(*(buff+1));
                set[key]++;
                break; 
              case VDW:
                g_assert(num_tokens > 1);
                elem.vdw = str_to_float(*(buff+1));
                set[key]++;
                break; 
              case CHARGE:
                g_assert(num_tokens > 1);
                elem.charge = str_to_float(*(buff+1));
                set[key]++;
                break; 
              case COLOUR:
                g_assert(num_tokens > 3);
                elem.colour[0] = str_to_float(*(buff+1));
                elem.colour[1] = str_to_float(*(buff+2));
                elem.colour[2] = str_to_float(*(buff+3));
if (VEC3MAGSQ(elem.colour) > 3.0)
  {
  VEC3MUL(elem.colour, 1.0/65535.0);
  }
                set[key]++;
                break; 
              case GDIS_END:
                flag++;
                break; 
              }
            }
          g_strfreev(buff); 
          }
/* store the element data according to type (ie global/model only etc.) */
        switch(type)
          {
          case DEFAULT:
            if (elem.number < 0 || elem.number > MAX_ELEMENTS-1)
              printf("Error: Element number %d outside allocated bounds.\n", elem.number);
            else
              memcpy(&elements[elem.number], &elem, sizeof(struct elem_pak));
            break;
          case MODIFIED:
/* TODO - scan for symbol as well (elem_seek with data = NULL) */
            if (!set[NUMBER])
              {
              printf("Warning: missing element number.\n");
              elem.number = 0;
              }
/* get default values */
            get_elem_data(elem.number, &elem_default, NULL);
/* fill out specified values */
            if (set[WEIGHT])
              elem_default.weight = elem.weight;
            if (set[COVALENT])
              elem_default.cova = elem.cova;
            if (set[VDW])
              elem_default.vdw = elem.vdw;
            if (set[CHARGE])
              elem_default.charge = elem.charge;
            if (set[COLOUR])
              {
              elem_default.colour[0] = elem.colour[0];
              elem_default.colour[1] = elem.colour[1];
              elem_default.colour[2] = elem.colour[2];
              }
            put_elem_data(&elem_default, NULL);
            break;
/* ignore other types */
          default:
            break;
          }
/* ignore all other keywords */
      default:
        break;
      }
    }

/* scattering factor coeffs */
  if (g_ascii_strncasecmp("\%gdis_sfc", line, 9) == 0)
    {
    while ((str = file_read_line(fp)))
      {
      buff = tokenize(str, &num_tokens);

if (num_tokens != 12)
  {
  /*
  printf("suspicious line: %s\n", str);
  */
  g_strfreev(buff);
  g_free(str);
  break;
  }
      symbol = g_strdup(*(buff+0));
      if (elem_symbol_test(symbol))
        {
        list = NULL;
        for (i=1 ; i<num_tokens ; i++)
          {
          sfc = g_malloc(sizeof(gdouble));
          *sfc = str_to_float(*(buff+i));
          list = g_slist_prepend(list, sfc);
          }
        list = g_slist_reverse(list);

        g_hash_table_insert(sysenv.sfc_table, symbol, list);

	g_strfreev(buff);
        }
      g_free(str);
      }
    }
  }
/* done */
g_free(set);
return(0);
}

/*****************************************************/
/* write the global exceptions to a specified stream */
/*****************************************************/
/* NB: should be suitable for appending (hence *fp) to gdisrc file */
gint write_elem_data(FILE *fp)
{
GSList *list;
struct elem_pak *elem;

/* checks */
g_return_val_if_fail(fp != NULL, 1);

list = sysenv.elements;
while(list != NULL)
  {
/* only write those values that are different to default library */
  elem = list->data;
  if (elem)
    {
    fprintf(fp, "%%gdis_elem\n");
    fprintf(fp, "number: %d\n", elem->number);
    fprintf(fp, "weight: %f\n", elem->weight);
    fprintf(fp, "  cova: %f\n", elem->cova);
    fprintf(fp, "   vdw: %f\n", elem->vdw);
    fprintf(fp, "charge: %f\n", elem->charge);
    fprintf(fp, "colour: %f %f %f\n", elem->colour[0], 
                                      elem->colour[1],
                                      elem->colour[2]);
    fprintf(fp, "%%gdis_end\n");
    }
  list = g_slist_next(list);
  }
return(0);
}

void write_sfc(gpointer key, gpointer val, gpointer data)
{
GSList *list;
FILE *fp = data;

fprintf(fp, "%s", (gchar *) key);
for (list=val ; list ; list=g_slist_next(list))
  fprintf(fp, " %f", *((gdouble *) list->data)); 
fprintf(fp, " \n");
}

void write_sfc_data(FILE *fp)
{
fprintf(fp, "%%gdis_sfc\n");
g_hash_table_foreach(sysenv.sfc_table, &write_sfc, fp);
fprintf(fp, "%%gdis_end\n");
}

/***************************************/
/* find all unique elements in a model */
/***************************************/
#define DEBUG_FIND_UNIQUE 0
GSList *find_unique(gint mode, struct model_pak *model)
{
gint found;
GSList *list, *types, *unique;
struct core_pak *core;

unique = NULL;

/* go through all atoms */
for (list=model->cores ; list ; list=g_slist_next(list))
  {
  core = list->data;

/* determine if current atom type is in the list */
  found=0;
  for (types=unique ; types ; types = g_slist_next(types))
    {
    switch(mode)
      {
      case ELEMENT:
        if (GPOINTER_TO_INT(types->data) == core->atom_code)
          found++;
        break;

      case LABEL_FF:
        if (g_ascii_strcasecmp(types->data, core->atom_type) == 0)
          found++;
        break;

      case LABEL:
      case LABEL_NORMAL:
      case LABEL_GHOST:
        if (g_ascii_strcasecmp(types->data, core->atom_label) == 0)
          found++;
        break;
      }
    if (found)
      break;
    }

/* not in list, so add it */
  if (!found)
    {
    switch(mode)
      {
      case ELEMENT:
        unique = g_slist_prepend(unique, GINT_TO_POINTER(core->atom_code));
        break;
      case LABEL:
        if (core->atom_label)
          unique = g_slist_prepend(unique, g_strdup(core->atom_label));
        break;
      case LABEL_NORMAL:
        if (core->atom_label)
          if (!core->ghost)
            unique = g_slist_prepend(unique, g_strdup(core->atom_label));
        break;
      case LABEL_GHOST:
        if (core->atom_label)
          if (core->ghost)
            unique = g_slist_prepend(unique, g_strdup(core->atom_label));
        break;
      case LABEL_FF:
        if (core->atom_type)
          unique = g_slist_prepend(unique, g_strdup(core->atom_type));
        break;
      }
    }
  }

#if DEBUG_FIND_UNIQUE
printf("Found: ");
for (types=unique ; types ; types = g_slist_next(types))
  {
  switch(mode)
    {
    case ELEMENT:
      printf("[%3d] ", GPOINTER_TO_INT(types->data));
      break;
    case LABEL:
    case LABEL_NORMAL:
    case LABEL_GHOST:
    case LABEL_FF:
      printf("[%s] ", (gchar *) types->data);
      break;
    }
  }
printf("\n");
#endif

/* return reversed list (preserves order) */
return(g_slist_reverse(unique));
}

/***********************/
/* set the atom colour */
/***********************/
void init_atom_colour(struct core_pak *core, struct model_pak *model)
{
struct elem_pak elem_data;

get_elem_data(core->atom_code, &elem_data, model);
ARR3SET(core->colour, elem_data.colour);
VEC3MUL(core->colour, 65535.0);
if (core->ghost)
  core->colour[3] = 0.5;
else
  core->colour[3] = 1.0;
}

/***********************/
/* set the atom charge */
/***********************/
void init_atom_charge(struct core_pak *core, struct model_pak *model)
{
struct elem_pak elem_data;

if (core->lookup_charge)
  {
  get_elem_data(core->atom_code, &elem_data, model);
  core->charge = elem_data.charge;
  }
}

/************************/
/* set the atom charges */
/************************/
void init_model_charges(struct model_pak *model)
{
GSList *list;

for (list=model->cores ; list ; list=g_slist_next(list))
  init_atom_charge(list->data, model);
}

/*******************************************/
/* net charge on atom (ie including shell) */
/*******************************************/
gdouble atom_charge(struct core_pak *core)
{
gdouble q;
struct shel_pak *shell;

q = core->charge;
if (core->shell)
  {
  shell = core->shell;
  q += shell->charge;
  }
return(q);
}

/********************************************/
/* initialize element properties for a core */
/********************************************/
void elem_init(struct core_pak *core, struct model_pak *model)
{
struct elem_pak elem;

g_assert(model != NULL);
g_assert(core != NULL);

/* NB: assumes core->atom_code is the atomic number */

get_elem_data(core->atom_code, &elem, model);
core->atom_code = elem.number;
core->bond_cutoff = elem.cova;
core->charge = elem.charge;
ARR3SET(core->colour, elem.colour);
VEC3MUL(core->colour, 65535.0);
}

/*********************************************************/
/* compute net charge and dipole (if model is a surface) */
/*********************************************************/
#define DEBUG_CALC_EMP 0
void calc_emp(struct model_pak *data)
{
gdouble qsum, x[3], p[3];
GSList *list;
struct core_pak *core;
struct shel_pak *shell;

/* checks */
g_assert(data != NULL);

/* TODO - speed up by getting to a local var the elem data for unique atoms */
/* cores & shell dipole calc */
qsum=0.0;
VEC3SET(p, 0.0, 0.0, 0.0);
for (list=data->cores ; list ; list=g_slist_next(list))
  {
  core = list->data;
/*
  if (core->region != REGION1A)
    continue;
*/
  if (core->status & DELETED)
    continue;

  ARR3SET(x, core->x);
  vecmat(data->latmat, x);

/* dipole */
  VEC3MUL(x, core->charge);
  ARR3ADD(p, x);
/* monopole */
  qsum += core->charge;

/* NB: add shell constribution at the core's z location */
  if (core->shell)
    {
    shell = core->shell;
    ARR3SET(x, shell->x);
    vecmat(data->latmat, x);

/* dipole */
    VEC3MUL(x, shell->charge);
    ARR3ADD(p, x);
/* monopole */
    qsum += shell->charge;
    }
  }

#if DEBUG_CALC_EMP
printf("surface: %f %f %f  (%f)\n", 
data->surface.miller[0], data->surface.miller[1], data->surface.miller[2], 
data->surface.shift);
P3VEC("dipole: ", p);
printf("sum charge: %e\n", qsum);
if (qsum*qsum > FRACTION_TOLERANCE)
  printf("Warning: your model has a net charge = %f\n", qsum);
#endif

/* only surfaces get the z dipole */
if (data->periodic == 2)
  data->gulp.sdipole = p[2];

data->gulp.qsum = qsum;
}

/***************************************************/
/* does a given string represent an atomic number? */
/***************************************************/
gint elem_number_test(const gchar *input)
{
gint n=0;

if (str_is_float(input))
  {
  n = str_to_float(input);
  if (n<0 || n>MAX_ELEMENTS-1)
    n=0; 
  }
return(n);
}

/***************************************************/
/* does a given string represent an element symbol */
/***************************************************/
#define DEBUG_ELEM_TYPE 0
gint elem_symbol_test(const gchar *input)
{
gint i, j, m, n;
gchar *elem1, **buff;

/* checks */
if (input == NULL)
  return(0);
/* zero length string matches something, so force it to return 0 */
if (!strlen(input))
  return(0);

/* duplicate for manipulation */
elem1 = g_strdup(input);

/* remove anything but alphabetic chars */
for (i=0 ; i<strlen(elem1) ; i++)
  if (!g_ascii_isalpha(*(elem1+i)))
    *(elem1+i) = ' ';

/* remove trailing/leading spaces & use only the first token */
buff = tokenize(elem1, &n);
if (n)
  {
  g_free(elem1);
  elem1 = g_strdup(*buff);
  }
g_strfreev(buff);

m = strlen(elem1);

/* catch Deuterium */
if (g_ascii_strncasecmp(elem1, "D", m) == 0)
  *elem1 = 'H';

#if DEBUG_ELEM_TYPE
printf("Looking for [%s]...", elem1);
#endif

/* attempt to match atom type with database */
j=0;
/* FIXME - elim dependence on const */
for(i=1 ; i<sysenv.num_elements ; i++)
  {
/* only compare if lengths match */
  n = strlen(elements[i].symbol);
  if (n == m)
    {
    if (g_ascii_strncasecmp(elem1, elements[i].symbol, m) == 0)
      {
      j = i;
      break;
      }
    }
  }

#if DEBUG_ELEM_TYPE
if (j)
  printf("found.\n");
else
  printf("not found.\n");
#endif

/* done */
g_free(elem1);
return(j);
}

/*************************************************************/
/* does a given string represent an element symbol or number */
/*************************************************************/
gint elem_test(const gchar *input)
{
gint n;

n = elem_symbol_test(input);
if (!n)
  n = elem_number_test(input);

return(n);
}

/*********************************/
/* PDB hack for element matching */
/*********************************/
#define DEBUG_PDB_ELEM_TYPE 0
gint pdb_elem_type(const gchar *input)
{
gint i, j, len;
gchar *tmp;

#if DEBUG_PDB_ELEM_TYPE
printf("pdb1: [%s]\n", input);
#endif

/* checks */
if (input == NULL)
  return(0);
/* zero length string matches something, so force it to return 0 */
len = strlen(input);
if (!len)
  return(0);

/* if we encounter C, O, or H then ignore the rest of the string */
for (i=0 ; i<len ; i++)
  {
  switch(*(input+i))
    {
    case 'c':
    case 'C':
    case 'o':
    case 'O':
    case 'h':
    case 'H':
      tmp = g_strndup(input, i+1);
      j = elem_symbol_test(tmp);
      g_free(tmp);
      return(j);
    }
  }
/* default element matching */
return(elem_symbol_test(input));
}

/*************************/
/* retrieve element data */
/*************************/
/* will return non-default data (if it exists) otherwise the default */
#define DEBUG_GET_ELEM 0
gint get_elem_data(gint code, struct elem_pak *elem, struct model_pak *model)
{
GSList *list;
struct elem_pak *elem_data;

/* search local exceptions */
if (model)
  {
  for (list=model->elements ; list ; list=g_slist_next(list))
    {
    elem_data = list->data;
    if (code == elem_data->number)
      {
#if DEBUG_GET_ELEM
printf("Retrieving local exception.\n");
#endif
      memcpy(elem, elem_data, sizeof(struct elem_pak));
      return(0);
      }
    }
  }
/* search global exceptions */
for (list=sysenv.elements ; list ; list=g_slist_next(list))
  {
  elem_data = list->data;
  if (code == elem_data->number)
    {
#if DEBUG_GET_ELEM
printf("Retrieving global exception.\n");
#endif
    memcpy(elem, elem_data, sizeof(struct elem_pak));
    return(0);
    }
  }

/* return from the standard database */
if (code >= 0 && code < sysenv.num_elements)
  {
#if DEBUG_GET_ELEM
printf("Retrieving default.\n");
#endif
  memcpy(elem, &elements[code], sizeof(struct elem_pak));
  return(0);
  }

#if DEBUG_GET_ELEM
printf("ERROR: bad element code.\n");
#endif
return(1);
}

/**********************************/
/* store non-default element data */
/**********************************/
#define DEBUG_PUT_ELEM_DATA 0
void put_elem_data(struct elem_pak *elem, struct model_pak *model)
{
gint flag, create;
GSList *list, *elements;
struct elem_pak *elem_data;
struct core_pak *core;

/* checks */
g_return_if_fail(elem != NULL);

#if DEBUG_PUT_ELEM_DATA
printf(" *** Saving element ***\n");
printf("symbol: %s\n", elem->symbol);
printf("  name: %s\n", elem->name);
printf("number: %d\n", elem->number);
printf("weight: %f\n", elem->weight);
printf("  cova: %f\n", elem->cova);
printf("   vdw: %f\n", elem->vdw);
printf("charge: %f\n", elem->charge);
printf("colour: %f %f %f\n", elem->colour[0], elem->colour[1], elem->colour[2]);
#endif

/* global or local database */
if (model)
  {
  elements = model->elements;
  flag = TRUE;
  }
else
  {
  elements = sysenv.elements;
  flag = FALSE;
  }

/* replace element if exists */
create = TRUE;
for (list=elements ; list ; list=g_slist_next(list))
  {
  elem_data = list->data;
  if (elem_data->number == elem->number)
    {
#if DEBUG_PUT_ELEM_DATA
printf("replacing existing...\n");
#endif
    memcpy(elem_data, elem, sizeof(struct elem_pak));
    create = FALSE;
    break;
    }
  }

/* create new item if it didn't exist */
if (create)
  {
#if DEBUG_PUT_ELEM_DATA
printf("creating new exception...\n");
#endif
  elem_data = g_malloc(sizeof(struct elem_pak));
  memcpy(elem_data, elem, sizeof(struct elem_pak));

  if (flag)
    model->elements = g_slist_prepend(model->elements, elem_data);
  else
    sysenv.elements = g_slist_prepend(sysenv.elements, elem_data);
  }

/* update atom bond_cutoff */
if (flag)
  {
  for (list=model->cores ; list ; list=g_slist_next(list))
    {
    core = list->data;
    if (core->atom_code == elem->number)
      core->bond_cutoff = elem->cova;
    }
  }
}

/*************************************************************/
/* get the nth colour in a sequence (eg sample RGB spectrum) */
/*************************************************************/
void get_colour(gdouble *colour, gint n)
{
gint byte;
gdouble r, g, b, f;

byte = 1 + n % 7;

if (byte & 4)
  r = 1.0;
else
  r = 0.0;

if (byte & 2)
  g = 1.0;
else
  g = 0.0;

if (byte & 1)
  b = 1.0;
else
  b = 0.0;

f = n / 7;
f *= 0.5;
f += 1.0;
f = 1.0/f;
r *= f;
g *= f;
b *= f;

VEC3SET(colour, r, g, b);
}

/****************************************/
/* change the colour scheme for an atom */
/****************************************/
void atom_colour_scheme(gint type, struct core_pak *core, struct model_pak *model)
{
gdouble t;
struct elem_pak elem;

switch (type)
  {
  case ELEM:
    get_elem_data(core->atom_code, &elem, model);
    break;

  case MOL:
/* get the sequence number for the colour type */
/* NB: core->molecule is deprec. */
    get_colour(elem.colour, core->molecule);
    break;

  case OCCUPANCY:
    VEC3MUL(core->colour, core->sof);
    return;

  case REGION:
/* get the sequence number for the colour type */
    get_colour(elem.colour, core->region);
    break;

  case GROWTH_SLICE:
    if (core->growth)
      {
      VEC3SET(elem.colour, 0.5, 0.5, 0.5);
      }
    else
      {
      get_elem_data(core->atom_code, &elem, model);
      }
    break;

  case TRANSLATE:
    if (core->translate)
      {
      VEC3SET(elem.colour, 0.5, 0.5, 0.5);
      }
    else
      {
      get_elem_data(core->atom_code, &elem, model);
      }
    break;

  case VELOCITY:
#define BC 1.3806503
#define AN 6.0221420
#define TMIN 0.0
#define TRANGE 2000.0

    get_elem_data(core->atom_code, &elem, model);
    t = elem.weight * VEC3MAGSQ(core->v) * 10.0 / (3.0*AN*BC);

/*
printf("T = %f\n", t);
*/

/* clamp */
t -= TMIN;
if (t > TRANGE)
  {
/*
  printf("tmax = %f\n", t + TMIN);
*/
  t = TRANGE;
  }
if (t < 0.0)
  t = 0.0;

/* bright red when hot (TMAX) */
    elem.colour[0] = t/TRANGE;
    elem.colour[1] = 0.0;
/* dark blue when cold */
    elem.colour[2] = 0.5 * (1.0 - t/TRANGE);
    break;
  }

ARR3SET(core->colour, elem.colour);
VEC3MUL(core->colour, 65535.0);
}

/*************************************/
/* set the colour scheme for a model */
/*************************************/
void model_colour_scheme(gint type, struct model_pak *model)
{
GSList *list;
struct core_pak *core;

model->colour_scheme = type;
for (list=model->cores ; list ; list=g_slist_next(list))
  {
  core = (struct core_pak *) list->data;
  atom_colour_scheme(type, core, model);
  }
}

