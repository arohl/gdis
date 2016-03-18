/*
Copyright (C) 2015 by Andrew Lloyd Rohl

a.rohl@curtin.edu.au

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
#include "matrix.h"
#include "interface.h"
#include "parse.h"

/* main structures */
extern struct sysenv_pak sysenv;
extern struct elem_pak elements[];

/****************/
/* read routine */
/****************/
#define DEBUG_READ_MOL2 0
gint read_mol2(gchar *filename, struct model_pak *data)
{
  gint num_tokens;
  gchar **buff, line[LINELEN];
  struct core_pak *core;
  FILE *fp;
  
  /* checks */
  g_return_val_if_fail(data != NULL, 1);
  g_return_val_if_fail(filename != NULL, 1);
  
  fp = fopen(filename, "rt");
  if (!fp)
    return(1);
  
  
  while (!fgetline(fp, line))
    {
    if (g_strrstr(line, "@<TRIPOS>ATOM") != NULL)
      {
      while (!fgetline(fp, line))
        {
        if (line[0] == '@')
          {
          break;
          }
        buff = tokenize(line, &num_tokens);
        if (num_tokens < 6)
          {
          gui_text_show(ERROR, "Invalid ATOM line in mol2 file\n");
          return(2);
          }
        else
          {
          core = new_core(*(buff+1), data);
          data->cores = g_slist_prepend(data->cores, core);
          
          core->x[0] = str_to_float(*(buff+2));
          core->x[1] = str_to_float(*(buff+3));
          core->x[2] = str_to_float(*(buff+4));
          
          if (num_tokens >= 6)
            {
            core->atom_type = g_strdup(*(buff+5));
            }
          if (num_tokens >= 9)
            {
            core->charge = str_to_float(*(buff+8));
            core->lookup_charge = FALSE;
            }
#if DEBUG_READ_MOL2
            printf("Atom %d %f %f %f Charge %f\n",n, core->x[0],
                 core->x[1], core->x[2], core->charge);
#endif
          }
        g_strfreev(buff);
        }
      }
    }

/* model setup */
  strcpy(data->filename, filename);
  g_free(data->basename);
  data->basename = parse_strip(filename);
  
  model_prep(data);
  fclose(fp);
  
  return(0);
}
