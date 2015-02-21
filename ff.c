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
#include "parse.h"
#include "coords.h"
#include "interface.h"
#include "ff.h"
#include "ff_pak.h"

extern struct sysenv_pak sysenv;

/************/
/* debuging */
/************/
void ff_dump(gpointer data)
{
gint i;
struct forcefield_pak *ff = data;

printf("[%d] ", ff->type);

for (i=0 ; i<ff->atoms_current ; i++)
  printf("[%s] ", ff->atom[i]);

if (ff->bond_expected)
  printf("{%f} ", ff->bond_value);

for (i=0 ; i<ff->data_current ; i++)
  printf("[%f] ", ff->data[i]);

printf("\n");
}

void ff_dump_type(gint type, GSList *ff_list)
{
GSList *list;
struct forcefield_pak *ff;

for (list=ff_list ; list ; list=g_slist_next(list))
  {
  ff = list->data;
  if (ff->type == type)
    ff_dump(ff);
  }
}

void ff_dump_all(GSList *ff_list)
{
GSList *list;

for (list=ff_list ; list ; list=g_slist_next(list))
  ff_dump(list->data);
}

/***********************************/
/* swap atom labels in a FF object */
/***********************************/
void ff_swap_atoms(struct forcefield_pak *ff, gint i, gint j)
{
gint n;
gchar tmp;

/*
for (n=0 ; n<FF_MAX_ATOMS ; n++)
   printf("[%s]", ff->atom[n]);
*/

for (n=0 ; n<FF_MAX_SYMBOL ; n++)
  {
  tmp = ff->atom[i][n];
  ff->atom[i][n] = ff->atom[j][n];
  ff->atom[j][n] = tmp;
  }

/*
printf(" -> ");
for (n=0 ; n<FF_MAX_ATOMS ; n++)
   printf("[%s]", ff->atom[n]);
printf("\n");
*/
}

/************************************/
/* atomic number matching primitive */
/************************************/
gint ff_match_code(struct forcefield_pak *ff, gint *atoms, gint n)
{
gint i, match;

g_assert(ff != NULL);
g_assert(n <= FF_MAX_ATOMS);

/* number of atoms match forcefield? */
if (n != ff->atoms_current)
  return(FALSE);

/* forward test */
match = TRUE;
for (i=0 ; i<ff->atoms_current ; i++)
  {
  if (ff->atom[i][0] == 'x' || ff->atom[i][0] == 'X')
    continue;

  if (atoms[i] != elem_symbol_test(ff->atom[i]))
    {
    match = FALSE;
    break;
    }
  }

if (match)
  return(TRUE);

/* backward test */
match = TRUE;
for (i=0 ; i<ff->atoms_current ; i++)
  {
  if (ff->atom[i][0] == 'x' || ff->atom[i][0] == 'X')
    continue;

  if (atoms[n-i-1] != elem_symbol_test(ff->atom[i]))
    {
    match = FALSE;
    break;
    }
  }

return(match);
}

/***********************************/
/* exact label matching primitives */
/***********************************/
/* FIXME - can return false positives */
/* eg csc could be either [cs][c] or [c][sc] */
#define DEBUG_FF_MATCH_LABEL 0
gint ff_match_label(struct forcefield_pak *ff, gchar **atoms, gint n)
{
gint i;
GString *label, *test;

g_assert(ff != NULL);

/* number of atoms match forcefield? */
if (n != ff->atoms_current)
  return(FALSE);

/* create the atom label string */
label = g_string_new(NULL);
for (i=0 ; i<n ; i++)
  {
  g_assert(atoms[i] != NULL);
  g_string_append(label, atoms[i]);
  }

/* create forward match string */
test = g_string_new(NULL);
for (i=0 ; i<ff->atoms_current ; i++)
  {
  if (ff->atom[i][0] == 'x' || ff->atom[i][0] == 'X')
    g_string_append(test, atoms[i]);
  else
    g_string_append(test, ff->atom[i]);
  }

#if DEBUG_FF_MATCH_LABEL
printf("[%s] : [%s]\n", label->str, test->str);
#endif

/* forward match */
if (label->len == test->len)
  if (g_ascii_strncasecmp(label->str, test->str, label->len) == 0)
    {
    g_string_free(label, TRUE);
    g_string_free(test, TRUE);
    return(TRUE);
    }

/* truncate */
g_string_set_size(test, 0);

/* create backward match string */
test = g_string_new(NULL);
for (i=0 ; i<ff->atoms_current ; i++)
  {
  if (ff->atom[i][0] == 'x' || ff->atom[i][0] == 'X')
    g_string_prepend(test, atoms[i]);
  else
    g_string_prepend(test, ff->atom[i]);
  }

/* backward match */
if (label->len == test->len)
  if (g_ascii_strncasecmp(label->str, test->str, label->len) == 0)
    {
    g_string_free(label, TRUE);
    g_string_free(test, TRUE);
    return(TRUE);
    }

g_string_free(label, TRUE);
g_string_free(test, TRUE);
return(FALSE);
}

/*********************************/
/* duplicate a forcefield object */
/*********************************/
gpointer ff_dup(gpointer data)
{
struct forcefield_pak *ff_copy, *ff_orig=data;

g_assert(ff_orig != NULL);

ff_copy = g_malloc(sizeof(struct forcefield_pak));

g_assert(ff_copy != NULL);

memcpy(ff_copy, ff_orig, sizeof(struct forcefield_pak));

return(ff_copy);
}

/**********************************/
/* create a new forcefield object */
/**********************************/
gpointer ff_type_new(gint type)
{
gint i;
struct forcefield_pak *ff=NULL;

ff = g_malloc(sizeof(struct forcefield_pak));

/*
printf("Creating new ff of type: %d\n", type);
*/

ff->type = type;
ff->atoms_current = 0;
ff->data_current = 0;

ff->data_units = FF_UNKNOWN;
ff->bond_units = FF_UNKNOWN;
ff->bond_index = -1;
ff->bond_value = 0.0;

for (i=FF_MAX_ATOMS ; i-- ; )
  ff->atom[i][0] = '\0';

for (i=FF_MAX_DATA ; i-- ; )
  ff->data[i] = 0.0;

switch (type)
  {
  case FF_HARMONIC:
    ff->atoms_expected = 2;
    ff->data_expected = 1;
    ff->bond_expected = 1;
    break;
  case FF_MORSE:
    ff->atoms_expected = 2;
    ff->data_expected = 2;
    ff->bond_expected = 1;
    break;
  case FF_BUCKINGHAM:
    ff->atoms_expected = 2;
    ff->data_expected = 3;
    ff->bond_expected = 0;
    break;
  case FF_LENNARD:
    ff->atoms_expected = 2;
    ff->data_expected = 2;
    ff->bond_expected = 0;
    break;
  case FF_3B_HARMONIC:
    ff->atoms_expected = 3;
    ff->data_expected = 1;
    ff->bond_expected = 1;
    break;
  case FF_DIHEDRAL:
  case FF_IMPROPER:
    ff->atoms_expected = 4;
    ff->data_expected = 1;
    ff->bond_expected = 1;
    break;

  case FF_DIHEDRAL_RB:
    ff->atoms_expected = 4;
    ff->data_expected = 6;
    ff->bond_expected = 0;
    break;

  default:
/* print warning? */
    ff->atoms_expected = 0;
    ff->data_expected = 0;
    ff->bond_expected = 0;
  }

return(ff);
}

/****************************************/
/* add more data to a forcefield object */
/****************************************/
/* return true if complete? */
gint ff_data_add(gpointer forcefield, gchar *text)
{
gint n=0, len, num_tokens;
gchar **buff;
gdouble x;
struct forcefield_pak *ff = forcefield;

/* checks */
if (!ff)
  return(FALSE);
if (!text)
  return(FALSE);

buff = tokenize(text, &num_tokens);

/* while atoms still expected ... */
while (n<num_tokens && ff->atoms_current < ff->atoms_expected)
  {
  if (elem_symbol_test(*(buff+n)))
    {
    g_assert(ff->atoms_current < FF_MAX_ATOMS);

    len = strlen(*(buff+n));
    if (len > FF_MAX_SYMBOL-1)
      len = FF_MAX_SYMBOL-1;

    g_snprintf(ff->atom[ff->atoms_current++], FF_MAX_SYMBOL, "%s", *(buff+n));
    }
  else
    {
/* wildcard */
    if (g_ascii_strncasecmp(*(buff+n), "x", 1) == 0)
      {
      g_snprintf(ff->atom[ff->atoms_current++], FF_MAX_SYMBOL, "%s", *(buff+n));
      }
    }

  n++;
  }

/* while data still expected ... */
while (n<num_tokens)
  {
  if (str_is_float(*(buff+n)))
    {
    x = str_to_float(*(buff+n));

    if (ff->data_current == ff->bond_index)
      {
      ff->bond_value = x;
/* prevent further matches */
      ff->bond_index = -1;
      }
    else
      {
      if (ff->data_current < FF_MAX_DATA)
        ff->data[ff->data_current++] = x; 
      }
    }
  n++;
  }

g_strfreev(buff);

return(FALSE);
}

/*****************************/
/* general forcefield search */
/*****************************/
/* NB: it's left to the caller to work out the bonded/non-bonded issues */
#define DEBUG_FF_GET 0
gpointer ff_search(struct core_pak **c, gint n, GSList *ff_list)
{
gint i, code[FF_MAX_ATOMS];
gchar *label[FF_MAX_ATOMS];
GSList *list;
struct forcefield_pak *ff, *ff_match=NULL;

#if DEBUG_FF_GET
printf("searching for interaction term, order %d\n", n);
#endif

g_assert(c != NULL);
g_assert(n > 0 && n < FF_MAX_ATOMS+1);

/* init */
for (i=0 ; i<n ; i++)
  {
/* level 1 - search via atom_type if exists, otherwise atom_label */
  if (c[i]->atom_type)
    label[i] = c[i]->atom_type;
  else
    label[i] = c[i]->atom_label;

#if DEBUG_FF_GET
printf(" [%s]", label[i]);
#endif

/* level 2 - search via atomic number */
  code[i] = c[i]->atom_code;
  }

#if DEBUG_FF_GET
printf("\n");
#endif

for (list=ff_list ; list ; list=g_slist_next(list))
  {
  ff = list->data;

/* exact atom label match test */
  if (ff_match_label(ff, label, n))
    {
    ff_match = ff;
    break;
    }

/* element match */
  if (ff_match_code(ff, code, n))
    ff_match = ff;
  }

return(ff_match);
}

/*****************************************/
/* filter a FF list using the given type */
/*****************************************/
/* NB: FF objs are not duplicated */
GSList *ff_filter_list(gint type, GSList *ff_list)
{
GSList *flist=NULL, *list;
struct forcefield_pak *ff;

for (list=ff_list ; list ; list=g_slist_next(list))
  {
  ff = list->data;
  if (ff->type == type)
    flist = g_slist_prepend(flist, ff);
  }
return(flist);
}




