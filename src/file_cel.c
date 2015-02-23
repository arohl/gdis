/*
Copyright (C) 2005 by Marcin Wojdyr  

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

/* PowderCell structure data format - read-only 
 *
 * example:
 * 
 * cell 5.64009 5.64009 5.64009 90  90  90
 * Na      11      0       0      0
 * Cl      17      0.5     0      0
 * rgnr 225
 *
 * For detailed description of file format see:
 * http://users.omskreg.ru/~kolosov/bam/a_v/v_1/powder/details/strucdat.htm
 * or mirror:
 * http://www.ccp14.ac.uk/ccp/web-mirrors/powdcell/a_v/v_1/powder/details/powcell.htm
 *
 */

#include <stdio.h>
#include <string.h>

#include "gdis.h"
#include "coords.h"
#include "model.h"
#include "file.h"
#include "parse.h"


/****************/
/* file reading */
/****************/
gint read_cel(gchar *filename, struct model_pak *model)
{
gchar *line;
FILE *fp;
int i;
gint num_tokens, natom=0;
gchar **buff;

struct core_pak *core;

/* checks */
g_return_val_if_fail(model != NULL, 1);
g_return_val_if_fail(filename != NULL, 2);

fp = fopen(filename, "rt");
if (!fp)
  return 3;

/* 1st line  - cell parameters */
line = file_read_line(fp);
if (!line || strlen(line) < 5 || g_ascii_strncasecmp("cell", line, 4) != 0)
  {
  printf("The first line should start with the keyword CELL.\n");
  return 4;
  }
buff = tokenize(line+4, &num_tokens);
g_free(line);
if (num_tokens < 6)
  {
  g_strfreev(buff);
  printf("Keyword CELL should be followed by six numbers.\n");
  return 5;
  }
for (i=0; i<3; ++i)
  model->pbc[i] = str_to_float(buff[i]);
for (i=3; i<6; ++i)
  model->pbc[i] = str_to_float(buff[i]) * D2R;
g_strfreev(buff);

/* next lines - atomic positions */
for (;;)
  {
  line = file_read_line(fp);
  if (!line) /*the end of file*/
    {
    printf("No 'rgnr' symmetry line found.\n");
    return 6;
    }
  else if (g_ascii_strncasecmp("rgnr", line, 4) == 0) /*no more atomic pos.*/
    {
    break;
    }
  else if (g_ascii_strncasecmp("natom", line, 5) == 0)/*number of atoms */
    /*some old .cel files have second line with number of atoms eg. "natom 6"*/
    {
    buff = tokenize(line, &num_tokens);
    if (num_tokens > 1)
      natom = str_to_float(buff[1]);
    else
      printf("Warning: ignoring `natom' line:\n%s\n", line);
    g_free(line);
    g_strfreev(buff);
    }
  else if (strncmp("    ", line, 4) != 0) /* atomic position */
    {
    buff = tokenize(line, &num_tokens);
    g_free(line);
    if (num_tokens < 5)
      {
      g_strfreev(buff);
      continue;
      }
    core = new_core(*buff, model);
    core->atom_label = g_strdup(buff[0]);

    /* in second column there is either atomic number 
     * or something like "Mg2+" or "K+". The second form is for
     * " the use of different bonding states of one and the same 
     *   element (e.g. Fe2+ and Fe3+ in Fe3O4)"
     * FIXME how these bonding states can be interpreted in GDIS */
    if (g_ascii_isdigit(buff[1][0]))
      core->atom_code = str_to_float(buff[1]);
    else 
      {
      core->atom_code = elem_symbol_test(buff[1]);
      }
    for (i=0; i<3; ++i)
      core->x[i] = str_to_float(buff[2+i]);
    /* TODO interpret 2 next optional numbers:
     *   so-called multiplied substitution and replacement factor (SOF)
     *   and isotropic Debye-Waller factor
     * FIXME can they be interpreted by GDIS? -MW */
    model->cores = g_slist_prepend(model->cores, core);
    g_strfreev(buff);
    }
  else /* replacement atom */
    {
    /*FIXME replacement atoms are now silently ignored 
     *  how can I use this information in GDIS? - MW*/ 
    g_free(line);
    }
  }

/* last line - symmetry */
buff = tokenize(line, &num_tokens);
model->sginfo.spacenum = str_to_float(buff[1]);
/* FIXME/TODO how to interpret the second (optional) number?
 * From fileformat docs: 
 * "sometimes there exists more than one setting of a space-group type. 
 * Thus, a further number must be given if the structure hasn't been described 
 * using a conventional setting (standard setting)."
 * http://users.omskreg.ru/~kolosov/bam/a_v/v_1/powder/details/strucdat.htm
 * http://users.omskreg.ru/~kolosov/bam/a_v/v_1/powder/details/setting.htm
 * Unfortunatelly I'm ignorant about space-groups - MW
 */
g_free(line);
g_strfreev(buff);

if (natom>0 && natom != g_slist_length(model->cores))
  printf("Warning: expected %i atoms, have %i.", natom, 
		                                 g_slist_length(model->cores));

/* model setup */
model->fractional = TRUE;
model->periodic = 3;
strcpy(model->filename, filename);
g_free(model->basename);
model->basename = parse_strip(filename);

model_prep(model);

return 0;
}

