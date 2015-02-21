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
#include "scan.h"
#include "interface.h"

/* main structures */
extern struct sysenv_pak sysenv;
extern struct elem_pak elements[];

/****************/
/* file writing */
/****************/
gint write_dlp(gchar *filename, struct model_pak *model)
{
return(0);
}

/****************/
/* file reading */
/****************/
#define DEBUG_READ_XYZ 0
gint read_dlp(gchar *filename, struct model_pak *model)
{
gint i, tokens;
gchar *line, **buff;
gpointer scan;
struct core_pak *core;

/* checks */
g_return_val_if_fail(model != NULL, 1);
g_return_val_if_fail(filename != NULL, 2);

scan = scan_new(filename);
if (!scan)
  return(3);

/* title line */
line = scan_get_line(scan);

/* config stuff */
line = scan_get_line(scan);

/* cell vectors */
model->construct_pbc = TRUE;
model->periodic = 3;
for (i=0 ; i<3 ; i++)
  {
  line = scan_get_line(scan);
  if (line)
    {
    buff = tokenize(line, &tokens);
    if (tokens > 2)
      {
      model->latmat[3*i] = str_to_float(*(buff));
      model->latmat[3*i+1] = str_to_float(*(buff+1));
      model->latmat[3*i+2] = str_to_float(*(buff+2));
      }
    g_strfreev(buff);
    }
  }

/* atoms & associated data */
model->fractional = FALSE;
while (!scan_complete(scan))
  {
  line = scan_get_line(scan);
  core = NULL;

/* element and id */
  buff = tokenize(line, &tokens);
  if (buff)
    core = core_new(*buff, NULL, model);
  g_strfreev(buff);

/* coordinates */
  line = scan_get_line(scan);
  buff = tokenize(line, &tokens);
  if (core && tokens > 2)
    {
    core->x[0] = str_to_float(*buff);
    core->x[1] = str_to_float(*(buff+1));
    core->x[2] = str_to_float(*(buff+2));
    }
  g_strfreev(buff);

/* velocity */
  line = scan_get_line(scan);
  buff = tokenize(line, &tokens);
  if (core && tokens > 2)
    {
    core->v[0] = str_to_float(*buff);
    core->v[1] = str_to_float(*(buff+1));
    core->v[2] = str_to_float(*(buff+2));
    }
  g_strfreev(buff);

/* TODO - acceleration??? */
  line = scan_get_line(scan);

/* finished reading atom data */
  if (core)
    model->cores = g_slist_prepend(model->cores, core);
  }


/*
strcpy(data->filename, filename);
g_free(data->basename);
data->basename = parse_strip(filename);
*/

model_prep(model);

scan_free(scan);

return(0);
}

