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
#include "matrix.h"
#include "model.h"
#include "parse.h"
#include "file.h"
#include "scan.h"

/* main structures */
extern struct sysenv_pak sysenv;
extern struct elem_pak elements[];

/****************/
/* file writing */
/****************/
gint write_dmol(gchar *filename, struct model_pak *model)
{
gint i;
gdouble x[3];
GSList *list;
struct core_pak *core;
FILE *fp;

/* checks */
g_return_val_if_fail(model != NULL, 1);
g_return_val_if_fail(filename != NULL, 2);

/* open the file */
fp = fopen(filename,"wt");
if (!fp)
  return(3);

if (model->periodic == 3)
  {
  fprintf(fp, "$cell vectors\n");

/* NB: DMOL matrices are transposed wrt gdis */
  for (i=0 ; i<3 ; i++)
    {
    fprintf(fp, "          %20.14f%20.14f%20.14f\n",
                          model->latmat[i]/AU2ANG,
                          model->latmat[i+3]/AU2ANG,
                          model->latmat[i+6]/AU2ANG);
    }
  }

fprintf(fp, "$coordinates\n");
for (list=model->cores ; list ; list=g_slist_next(list))
  {
  core = list->data;
  if (core->status & DELETED)
    continue;

/* everything is cartesian after latmat mult */
  ARR3SET(x, core->x);
  vecmat(model->latmat, x);
  fprintf(fp,"%-10s%20.14f%20.14f%20.14f%5d\n",
              elements[core->atom_code].symbol,
              x[0]/AU2ANG, x[1]/AU2ANG, x[2]/AU2ANG, core->atom_code);
  }
fprintf(fp, "$end\n");

fclose(fp);
return(0);
}

/****************/
/* file reading */
/****************/
gint read_dmol_frame(FILE *fp, struct model_pak *model)
{
gint i, num_tokens;
gchar *line, **buff;
struct core_pak *core;

g_assert(fp != NULL);

line = file_read_line(fp);

while (line)
  {
/* read cell vectors */
  if (g_ascii_strncasecmp(line, "$cell", 5) == 0 && model)
    {
    for (i=0 ; i<3 ; i++)
      {
      g_free(line);
      line = file_read_line(fp);
      buff = tokenize(line, &num_tokens);
      if (num_tokens > 2)
        {
        model->latmat[i] = AU2ANG*str_to_float(*(buff));
        model->latmat[i+3] = AU2ANG*str_to_float(*(buff+1));
        model->latmat[i+6] = AU2ANG*str_to_float(*(buff+2));
        }
      g_strfreev(buff);
      }
    model->periodic = 3;
    model->construct_pbc = TRUE;
    }

/* read coordinates */
  if (g_ascii_strncasecmp(line, "$coord", 5) == 0 && model)
    {
    g_free(line);
    line = file_read_line(fp);
    buff = tokenize(line, &num_tokens);
    while (num_tokens > 3)
      {
      if (elem_symbol_test(*buff))
        {
        core = new_core(*buff, model);
        model->cores = g_slist_prepend(model->cores, core);
        core->x[0] = AU2ANG*str_to_float(*(buff+1));
        core->x[1] = AU2ANG*str_to_float(*(buff+2));
        core->x[2] = AU2ANG*str_to_float(*(buff+3));
        }
      g_free(line);
      line = file_read_line(fp);
      g_strfreev(buff);
      buff = tokenize(line, &num_tokens);
      }
    g_strfreev(buff);
    model->fractional = FALSE;
    }

/* terminate frame read */
  if (g_ascii_strncasecmp(line, "$end", 4) == 0)
    return(0);

  g_free(line);
  line = file_read_line(fp);
  }

return(1);
}

/****************/
/* file reading */
/****************/
#define DEBUG_READ_XYZ 0
gint read_dmol(gchar *filename, struct model_pak *model)
{
gint flag;
FILE *fp;

/* checks */
g_return_val_if_fail(model != NULL, 1);
g_return_val_if_fail(filename != NULL, 2);

fp = fopen(filename,"rt");
if (!fp)
  return(3);

/* loop while there's data */
flag=0;
model->num_frames = 0;

read_dmol_frame(fp, model);

for (;;)
  {
  add_frame_offset(fp, model);

  if (read_dmol_frame(fp, NULL))
    break;

  model->num_frames++;
  }

/* get rid of frame list if only one frame */
if (model->num_frames == 1)
  {
  free_list(model->frame_list);
  model->frame_list = NULL;
  }

/* model setup */
strcpy(model->filename, filename);
g_free(model->basename);
model->basename = parse_strip(filename);
model_prep(model);

return(0);
}
