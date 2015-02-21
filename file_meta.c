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

#include "gdis.h"
#include "coords.h"
#include "model.h"
#include "file.h"
#include "parse.h"
#include "select.h"
#include "interface.h"

extern struct elem_pak elements[];

/* meta data for computational chemistry archiving */

/*****************************************************************/
/* split a string containing units into it's raw value and units */
/*****************************************************************/
void meta_parse_units(gchar **value, gchar **units, gchar *source)
{
gchar **buff;

buff = g_strsplit(source, " ", 2);

if (*buff+0)
  *value = g_strdup(*(buff+0));
else
  *value = NULL;

if (*buff+1)
  *units = g_strdup(*(buff+1));
else
  *units = NULL;

g_strfreev(buff);
}

/****************/
/* file writing */
/****************/
gint write_meta(gchar *filename, struct model_pak *model)
{
gchar *energy_package=NULL, *energy_value=NULL, *energy_gradient=NULL;
gchar *energy_method=NULL, *energy_basis=NULL, *energy_functional=NULL;
gchar *value=NULL, *units=NULL;
GSList *list;
FILE *fp;

/* checks */
g_return_val_if_fail(model != NULL, 1);
g_return_val_if_fail(filename != NULL, 2);

/* open the file */
fp = fopen(filename,"wt");
if (!fp)
  return(3);

/* fill in the blanks ... */
/* <key> separator <item1> separator <item 2> ... etc */
/* separator atm will be '|' */

/* structural information */
fprintf(fp, "formula|");
for (list=model->unique_atom_list ; list ; list=g_slist_next(list))
  {
  int code = GPOINTER_TO_INT(list->data);
  select_clear(model);
  select_all_elem(code, model);
  fprintf(fp, "%s %d ", elements[code].symbol, g_slist_length(model->selection));
  }
fprintf(fp, "\n");
fprintf(fp, "molecules|%d\n", g_slist_length(model->moles));
fprintf(fp, "dimensionality|%d\n", model->periodic);
switch(model->periodic)
  {
  case 3:
    fprintf(fp, "cell|%f %f %f %f %f %f\n",
                 model->pbc[0], model->pbc[1], model->pbc[2],
                 R2D*model->pbc[3], R2D*model->pbc[4], R2D*model->pbc[5]);
    fprintf(fp, "spacegroup|%s\n", model->sginfo.spacename);
    break;

  case 2:
    fprintf(fp, "cell|%f %f %f \n",
                 model->pbc[0], model->pbc[1], R2D*model->pbc[5]);
    break;

  case 1:
    fprintf(fp, "cell|%f\n", model->pbc[0]);
    break;
  }

/* energetics information */
/* only for recognized output types ... input files only have structual meta data stored */
/* check eligible types for energy etc output */
/* TODO - include a version number (if available) in package??? */
energy_value = property_lookup("Energy", model);

switch(model->id)
  {
  case ABINIT_OUT:
    energy_package = g_strdup("ABINIT");
    energy_gradient = property_lookup("RMS Gradient", model);
    break;
  case CASTEP_OUT:
    energy_package = g_strdup("CASTEP");
    energy_gradient = property_lookup("RMS Gradient", model);
    energy_functional = property_lookup("Functional", model);
    break;
  case GAMESS_OUT:
    energy_package = g_strdup("GAMESS");
    energy_gradient = property_lookup("RMS Gradient", model);
    energy_basis = property_lookup("Basis", model);
    break;
  case GAUSS_OUT:
    energy_package = g_strdup("GAUSSIAN");
    energy_method = property_lookup("Method", model);
    energy_gradient = property_lookup("RMS Gradient", model);
    energy_basis = property_lookup("Basis", model);
    break;
  case GULPOUT:
    energy_package = g_strdup("GULP");
    energy_method = g_strdup("Forcefield");
    energy_value = property_lookup("Energy", model);
    energy_gradient = property_lookup("Gnorm", model);
    break;
  case MVNOUT:
    energy_package = g_strdup("MARVIN");
    energy_method = g_strdup("Forcefield");
    energy_gradient = property_lookup("Gnorm", model);
    break;
  case NWCHEM_OUT:
    energy_package = g_strdup("NWCHEM");
    energy_gradient = property_lookup("RMS Gradient", model);
    break;
  case SIESTA_OUT:
    energy_package = g_strdup("SIESTA");
    energy_gradient = property_lookup("Maximum Gradient", model);
    energy_functional = property_lookup("Functional", model);
    break;
  }

if (energy_package)
  {
  fprintf(fp, "package|%s\n", energy_package);
  g_free(energy_package);
  }
if (energy_method)
  {
  fprintf(fp, "method|%s\n", energy_method);
  g_free(energy_method);
  }
if (energy_basis)
  {
  fprintf(fp, "basis|%s\n", energy_basis);
  g_free(energy_basis);
  }
if (energy_functional)
  {
  fprintf(fp, "functional|%s\n", energy_functional);
  g_free(energy_functional);
  }

if (energy_value)
  {
  meta_parse_units(&value, &units, energy_value);
  if (value)
    {
    fprintf(fp, "energy|%s", value);
    g_free(value);
    if (units)
      {
      fprintf(fp, "|%s", units);
      g_free(units);
      }
    }
  fprintf(fp, "\n");
  g_free(energy_value);
  }

if (energy_gradient)
  {
  meta_parse_units(&value, &units, energy_gradient);
  if (value)
    {
    fprintf(fp, "gradient|%s", value);
    g_free(value);
    if (units)
      {
      fprintf(fp, "|%s", units);
      g_free(units);
      }
    }
  fprintf(fp, "\n");
  g_free(energy_gradient);
  }

fclose(fp);
return(0);
}

