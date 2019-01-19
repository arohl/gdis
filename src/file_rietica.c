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
#include "file.h"
#include "parse.h"
#include "model.h"
#include "scan.h"
#include "matrix.h"
#include "interface.h"
#include "opengl.h"

#define DEBUG_MORE 0
#define MAX_KEYS 15

/* main structures */
extern struct sysenv_pak sysenv;
extern struct elem_pak elements[];

/****************/
/* file reading */
/****************/
#define DEBUG_READ_RIETICA 0
gint read_rietica(gchar *filename, struct model_pak *model)
{
gint i, phases=0, skip, num_tokens;
gchar **buff, *line;
float x, y, z;
gpointer scan;
GSList *list;
struct core_pak *core;

/* checks */
g_assert(model != NULL);
scan = scan_new(filename);
if (!scan)
  return(1);

/* FIXME - stop the previous file routines setting this */
model->id = -1;

while (!scan_complete(scan))
  {
  buff = scan_get_tokens(scan, &num_tokens);
  
/* search for structure start */
  if (num_tokens)
    {
    if (g_ascii_strncasecmp(*buff, "***", 3) == 0)
      {
      if (phases)
        model = model_new();
      phases++;

/* structure name - omit the 1st and last tokens (ie "***") */
      if (num_tokens > 1)
        {
        g_free(*(buff+num_tokens-1));
        *(buff+num_tokens-1) = NULL;
        g_free(model->basename);
        model->basename = g_strjoinv(" ", buff+1);
        }

/* parse spacegroup */
      if(scan_get_line(scan)){ /*FIX de93af*/
	line = scan_get_line(scan);
	model->sginfo.spacename = g_strstrip(g_strdup(line));
	model->sginfo.spacenum = -1;
      }

/* parse a structure */
      skip = 0;
      while (!scan_complete(scan))
        {
        g_strfreev(buff);
        buff = scan_get_tokens(scan, &num_tokens);

        if (num_tokens)
          {
          if (elem_symbol_test(*buff))
            {
/* new core */
/*
            if (num_tokens > 6)
*/
              {
              core = new_core(*buff, model);
              model->cores = g_slist_prepend(model->cores, core);

/* formatted output can result in -ve signs "joining" tokens */
              line = scan_cur_line(scan);
/* no doubt some fortran programmer thought this was a clever format */
              sscanf(line, "%*16c%8f%8f%8f", &x, &y, &z);
              VEC3SET(core->x, x, y, z);

/*
printf("> %s", line);
P3VEC(" - ", core->x);
              core->x[0] = str_to_float(*(buff+2));
              core->x[1] = str_to_float(*(buff+3));
              core->x[2] = str_to_float(*(buff+4));
              core->sof = str_to_float(*(buff+6));
*/

              skip = 0;
              }
            }
          else
            skip++;
          }

/* 4 lines after the last core - parse cell info and terminate structure */
        if (skip == 4)
          {
          if (num_tokens > 5)
            {
            for (i=6 ; i-- ; )
              model->pbc[i] = str_to_float(*(buff+i));
            model->pbc[3] *= D2R;
            model->pbc[4] *= D2R;
            model->pbc[5] *= D2R;
            }
          break;
          }
        }
      }
    }
  g_strfreev(buff);
  }

/* setup all new structures */
for (list=sysenv.mal ; list ; list=g_slist_next(list))
  {
  model = list->data;

  if (model->id == -1)
    {
    model->id = RIETICA;
    model->periodic = 3;
    model->fractional = TRUE;
    strcpy(model->filename, filename);
    model->cores = g_slist_reverse(model->cores);
    model_prep(model);
    }
  }

scan_free(scan);

return(0);
}

