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

#include "gdis.h"
#include "type.h"
#include "parse.h"
#include "coords.h"
#include "interface.h"
#include "ff.h"
#include "ff_pak.h"

extern struct sysenv_pak sysenv;

/*****************************************************/
/* create a new FF object from a GULP potential line */
/*****************************************************/
gpointer ff_gulp_new(const gchar *name)
{
gint type=-1, bond_units, data_units, bond_index=-1;
struct forcefield_pak *ff=NULL;

if (!name)
  return(NULL);

if (g_ascii_strncasecmp("harm", name, 4) == 0)
  {
  type = FF_HARMONIC;
  bond_units = FF_ANG;
  data_units = FF_EV;
  bond_index = 1;
  }
else if (g_ascii_strncasecmp("mors", name, 4) == 0)
  {
  type = FF_MORSE;
  bond_units = FF_ANG;
  data_units = FF_EV;
  bond_index = 2;
  }
else if (g_ascii_strncasecmp("buck", name, 4) == 0)
  {
  type = FF_BUCKINGHAM;
  bond_units = FF_ANG;
  data_units = FF_EV;
  }
else if (g_ascii_strncasecmp("lenn", name, 4) == 0)
  {
  type = FF_LENNARD;
  bond_units = FF_ANG;
  data_units = FF_EV;
  }
else if (g_ascii_strncasecmp("thre", name, 4) == 0)
  {
  type = FF_3B_HARMONIC;
  bond_units = FF_DEG;
  data_units = FF_EV;
  bond_index = 1;
  }
else if (g_ascii_strncasecmp("tors", name, 4) == 0)
  {
  type = FF_DIHEDRAL;  /* TODO - check connectivity for prop/improp? :-( */
  bond_units = FF_DEG;
  data_units = FF_EV;
  bond_index = 2;
  }

if (type != -1)
  {
  ff = ff_type_new(type);

  g_assert(ff != NULL);

  ff->bond_index = bond_index;
  ff->bond_units = bond_units;
  ff->data_units = data_units;

  if (g_strrstr(name, "k3"))
    {
    ff->data_expected++;
    ff->bond_index++;
    }
  if (g_strrstr(name, "k4"))
    {
    ff->data_expected++;
    ff->bond_index++;
    }
  if (g_strrstr(name, "kcal"))
    ff->data_units = FF_KCAL;
  }

return(ff);
}

/***********************************************/
/* process a block of text for GULP potentials */
/***********************************************/
#define DEBUG_FF_GULP_PARSE 0
GSList *ff_gulp_parse(const gchar *text)
{
gint i;
gchar **line;
gpointer test;
GSList *list=NULL;
struct forcefield_pak *ff=NULL;

#if DEBUG_FF_GULP_PARSE
printf("---------------\n");
printf("%s", text);
printf("---------------\n");
#endif

/* split into lines */
if (!text)
  return(NULL);
line = g_strsplit(text, "\n", -1);
if (!line)
  return(NULL);

/* loop over lines */
i=0;
while (*(line+i))
  {
/* recognizable FF type? */
  test = ff_gulp_new(*(line+i));
  if (test)
    {
/* yes - create new FF object */
    ff = test;
    list = g_slist_prepend(list, ff);
    }
  else
    {
/* no - enter as current FF data (if exists) */
    if (ff)
      ff_data_add(ff, *(line+i));
    }
  i++;
  }

g_strfreev(line);

list = g_slist_reverse(list);

#if DEBUG_FF_GULP_PARSE
ff_dump_all(list);
#endif

return(list);
}

/**************************************/
/* create a FF output string for GULP */
/**************************************/
#define DEBUG_FF_GULP 0
gchar *ff_gulp_string(gpointer data)
{
gint i;
gdouble x;
GString *text;
struct forcefield_pak *ff;

g_assert(data != NULL);

text = g_string_new(NULL);

/* duplicate forcefield so we can play with atom order etc for GULP output */
ff = ff_dup(data);

g_assert(ff != NULL);

/* line 1 - type */
switch (ff->type)
  {
  case FF_HARMONIC:
    switch (ff->data_expected)
      {
      case 3:
        text = g_string_assign(text, "harmonic k3 k4 bond ");
        break;
      case 2:
        text = g_string_assign(text, "harmonic k3 bond ");
        break;
      default:
        text = g_string_assign(text, "harmonic bond ");
      }
    break;
  case FF_MORSE:
    text = g_string_assign(text, "morse bond ");
    break;
  case FF_BUCKINGHAM:
    text = g_string_assign(text, "buckingham ");
    break;
  case FF_LENNARD:
    text = g_string_assign(text, "lennard ");
    break;
  case FF_3B_HARMONIC:
    switch (ff->data_expected)
      {
      case 3:
        text = g_string_assign(text, "three bond k3 k4 ");
        break;
      case 2:
        text = g_string_assign(text, "three bond k3 ");
        break;
      default:
        text = g_string_assign(text, "three bond ");
      }
/* GULP angles have the middle atom first */
/* GDIS uses connectivity trace order */
    ff_swap_atoms(ff, 0, 1);
    break;
  case FF_DIHEDRAL:
    text = g_string_assign(text, "torsion bond ");
    break;
  case FF_DIHEDRAL_RB:
    text = g_string_assign(text, "ryckaert ");
    break;
  case FF_IMPROPER:
/* HACK - implement as out of plane term, rather than muck around with cutoffs */
/* FIXME - k term need scaling? */
    text = g_string_assign(text, "outofplane intra ");
    ff_swap_atoms(ff, 0, 1);
    break;

  default:
    return(NULL);
  }

/* line 1 - units */
switch (ff->data_units)
  {
  case FF_KJ:
    text = g_string_append(text, "kjmol ");
    break;
  case FF_KCAL:
    text = g_string_append(text, "kcal ");
    break;
  case FF_EV:
  default:
    break;
  }

text = g_string_append(text, "\n");

/* line 2 - atoms */
for (i=0 ; i<ff->atoms_expected ; i++)
  g_string_append_printf(text, "%s ", ff->atom[i]);

/* line 2 - data */
for (i=0 ; i<ff->data_expected ; i++)
  g_string_append_printf(text, "%f ", ff->data[i]);
  
/* line 2 - bond */
if (ff->bond_expected)
  {
  x = ff->bond_value;
  switch (ff->bond_units)
    {
    case FF_AU:
      x *= AU2ANG;
      break;
    case FF_RAD:
      x *= R2D; 
      break;
    case FF_ANG:
    case FF_DEG:
    default: 
      break;
    }
  g_string_append_printf(text, "%f ", x);
  }

text = g_string_append(text, "\n");

#if DEBUG_FF_GULP
printf("-------------\n");
printf("%s", text->str);
printf("-------------\n");
#endif

g_free(ff);

return(g_string_free(text, FALSE));
}

/* NULL terminated variable argument char list */
/*
void dump_labels(gchar *a, ...)
{
va_list arg_list;
gchar *arg;
    
va_start(arg_list, a);
arg = a;
while (arg != NULL)
  {
  arg = va_arg(arg_list, gchar);

  printf("%s ", arg);
  }
printf("\n");
    
va_end(arg_list); 
}
*/

/***************************************************************/
/* build a list of unique two body terms from a list of labels */
/***************************************************************/
GSList *list_2b(GSList *labels)
{
gint i, j, m;
gchar **la;
GSList *list, *output;
struct forcefield_pak *ff;

m = g_slist_length(labels);

la = g_malloc(m * sizeof(gchar *));

i=0;
output = NULL;
for (list=labels ; list ; list=g_slist_next(list))
  la[i++] = g_strdup(list->data);

for (i=0 ; i<m ; i++)
  {
  for (j=i ; j<m ; j++)
    {
    ff = ff_type_new(-1);

    ff->atoms_expected = ff->atoms_current = 2;
    ff->bond_expected = FALSE;
    g_snprintf(ff->atom[0], FF_MAX_SYMBOL, "%s", la[i]);
    g_snprintf(ff->atom[1], FF_MAX_SYMBOL, "%s", la[j]);

    output = g_slist_prepend(output, ff);
    }
  }

return(output);
}

/*****************************************************************/
/* build a list of unique three body terms from a list of labels */
/*****************************************************************/
GSList *list_3b(GSList *labels)
{
gint i, j, k, m;
gchar **la;
GSList *list, *output;
struct forcefield_pak *ff;

m = g_slist_length(labels);

la = g_malloc(m * sizeof(gchar *));

i=0;
output = NULL;
for (list=labels ; list ; list=g_slist_next(list))
  la[i++] = g_strdup(list->data);

for (j=0 ; j<m ; j++)
  {
  for (i=0 ; i<m ; i++)
    {
    for (k=i ; k<m ; k++)
      {
      ff = ff_type_new(-1);

      ff->atoms_expected = ff->atoms_current = 3;
      ff->bond_expected = FALSE;
      g_snprintf(ff->atom[0], FF_MAX_SYMBOL, "%s", la[i]);
      g_snprintf(ff->atom[1], FF_MAX_SYMBOL, "%s", la[j]);
      g_snprintf(ff->atom[2], FF_MAX_SYMBOL, "%s", la[k]);

      output = g_slist_prepend(output, ff);
      }
    }
  }

return(output);
}

/**************************/
/* write non bonded terms */
/**************************/
gint ff_gulp_write_vdw(FILE *fp, GSList *types, struct model_pak *model)
{
gchar *label[FF_MAX_ATOMS], *text;
GSList *item1, *item2, *list2;
struct forcefield_pak *ff, *ff_match;

g_assert(model != NULL);

/* build unique two body interaction symbols */
list2 = list_2b(types);
for (item1=list2 ; item1 ; item1=g_slist_next(item1))
  {
  ff = item1->data;
  label[0] = ff->atom[0];
  label[1] = ff->atom[1];

/* search for valid FF term for current interaction */
  for (item2=model->ff_list ; item2 ; item2=g_slist_next(item2))
    {
    ff_match = item2->data;

    switch(ff_match->type)
      {
      case FF_BUCKINGHAM:
      case FF_LENNARD:
        if (ff_match_label(ff_match, label, 2))
          {
/*
printf("vdw: %s - %s\n", label[0], label[1]);
*/
          text = ff_gulp_string(ff_match);
          fprintf(fp, "%s\n", text);
          g_free(text);
          }
        break;
      }
    }
  }

return(0);
}

/**********************/
/* write bonded terms */
/**********************/
#define DEBUG_FF_GULP_WRITE_BONDS 1
gint ff_gulp_write_bonds(FILE *fp, GSList *types, struct model_pak *model)
{
gchar *label[FF_MAX_ATOMS], *text;
GSList *item1, *item2;
GSList *list2;
struct core_pak *c[4];
struct forcefield_pak *ff, *ff_match;

g_assert(model != NULL);


/*
ff_dump_type(FF_HARMONIC, model->ff_list);
*/


/* build unique two body interaction symbols */
list2 = list_2b(types);

/* enumerate bonds - look for instances */
for (item1=model->bonds ; item1 ; item1=g_slist_next(item1))
  {
  struct bond_pak *bond = item1->data;

  c[0] = bond->atom1;
  c[1] = bond->atom2;

  ff = ff_search(c, 2, list2);
  if (ff)
    ff->bond_expected = TRUE;
  }

/* enumerate instanced bonds - look for forcefield term */
for (item1=list2 ; item1 ; item1=g_slist_next(item1))
  {
  ff = item1->data;

#if DEBUG_FF_GULP_WRITE_BONDS
printf("[%s][%s] : %d : ", ff->atom[0], ff->atom[1], ff->bond_expected);
#endif

  if (ff->bond_expected)
    {
    label[0] = ff->atom[0];
    label[1] = ff->atom[1];

    for (item2=model->ff_list ; item2 ; item2=g_slist_next(item2))
      {
      ff_match = item2->data;

      switch(ff_match->type)
        {
        case FF_HARMONIC:
        case FF_MORSE:
          if (ff_match_label(ff_match, label, 2))
            {
            text = ff_gulp_string(ff_match);
            fprintf(fp, "%s\n", text);
            g_free(text);
#if DEBUG_FF_GULP_WRITE_BONDS
printf("*");
#endif
            }
          break;

        default:
          continue;
        }
      }
    }
#if DEBUG_FF_GULP_WRITE_BONDS
printf("\n");
#endif
  }
return(0);
}

/*********************************/
/* write bonded three body terms */
/*********************************/
#define DEBUG_FF_GULP_WRITE_ANGLES 1
gint ff_gulp_write_angles(FILE *fp, GSList *types, struct model_pak *model)
{
gint n;
gchar *label[FF_MAX_ATOMS], *text;
GSList *item1, *item2;
GSList *list3, *nlist;
struct core_pak *c[3];
struct forcefield_pak *ff, *ff_match;

g_assert(model != NULL);

/* build unique three body interaction symbols */
list3 = list_3b(types);

/* enumerate cores - look for instances */
for (item1=model->cores ; item1 ; item1=g_slist_next(item1))
  {
  nlist = connect_neighbours(item1->data);
  n = g_slist_length(nlist);

  if (n < 2)
    continue;

/* CURRENT - only checks for angle using first two bonded atoms */
/* FIXME - need to check all combinations if more than 2 bonded atoms */
  c[1] = item1->data;
  item2 = nlist;
  c[0] = item2->data;
  item2 = g_slist_next(item2);
  c[2] = item2->data;

  ff = ff_search(c, 3, list3);
  if (ff)
    ff->bond_expected = TRUE;
  }

/* enumerate instanced bonds - look for forcefield term */
for (item1=list3 ; item1 ; item1=g_slist_next(item1))
  {
  ff = item1->data;

#if DEBUG_FF_GULP_WRITE_ANGLES
printf("[%s][%s][%s] : %d : ", ff->atom[0], ff->atom[1], ff->atom[2], ff->bond_expected);
#endif

  if (ff->bond_expected)
    {
    label[0] = ff->atom[0];
    label[1] = ff->atom[1];
    label[2] = ff->atom[2];

    for (item2=model->ff_list ; item2 ; item2=g_slist_next(item2))
      {
      ff_match = item2->data;

      switch(ff_match->type)
        {
        case FF_3B_HARMONIC:
          if (ff_match_label(ff_match, label, 3))
            {
            text = ff_gulp_string(ff_match);
            fprintf(fp, "%s\n", text);
            g_free(text);
            }
          break;

        default:
          continue;
	}
      }
   }
#if DEBUG_FF_GULP_WRITE_ANGLES
printf("\n");
#endif
  }

return(0);
}

/******************************************************************/
/* write a GULP forcefield section using internal forcefield list */
/******************************************************************/
#define DEBUG_FF_GULP_WRITE 0
gint ff_gulp_write(FILE *fp, struct model_pak *model)
{
gint n;
GSList *list;

g_assert(model != NULL);

/* TODO - build equivalence list - atom_type <-> atom_label */
/* since GULP only does FF lookup via atom_label which must be a valid element */
/* whereas many FF's have to (stupidly) be referenced via atom_types (eg CA = carbon, not calcium) */
/* or for all unique types - forcibly assign unique atom_labels */
/* eg for each carbon FF type -> atom label: C1, C2, C3 ... */

n = type_check_list(model->cores);
if (n)
  {
  gui_text_show(WARNING, "Forcefield output may be incomplete; check your atom types.\n");
  return(1);
  }

list = find_unique(LABEL_FF, model);

if (list)
  {
  ff_gulp_write_vdw(fp, list, model);
  ff_gulp_write_bonds(fp, list, model);
  ff_gulp_write_angles(fp, list, model);

/* TODO - torsions , improper (oop) */
  }

return(0);
}
