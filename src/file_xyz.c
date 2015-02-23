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
#include "matrix.h"
#include "interface.h"

#define DEBUG_MORE 0
#define MAX_KEYS 15

/* main structures */
extern struct sysenv_pak sysenv;
extern struct elem_pak elements[];

/****************/
/* file writing */
/****************/
gint write_xyz(gchar *filename, struct model_pak *data)
{
gdouble x[3];
GSList *list;
struct core_pak *core;
FILE *fp;

/* checks */
g_return_val_if_fail(data != NULL, 1);
g_return_val_if_fail(filename != NULL, 2);

/* open the file */
fp = fopen(filename,"wt");
if (!fp)
  return(3);

/* print header */
fprintf(fp,"%d\nXYZ\n", g_slist_length(data->cores));

for (list=data->cores ; list ; list=g_slist_next(list))
  {
  core = list->data;
  if (core->status & DELETED)
    continue;

/* everything is cartesian after latmat mult */
  ARR3SET(x, core->x);
  vecmat(data->latmat, x);
  fprintf(fp,"%-7s    %14.9f  %14.9f  %14.9f\n",
              elements[core->atom_code].symbol, x[0], x[1], x[2]);
  }

fclose(fp);
return(0);
}

/*************************************************/
/* read an xyz block into the model structure    */
/*************************************************/
#define READ_XYZ_BLOCK 0
gint read_xyz_block(FILE *fp, struct model_pak *data)
{
gint expect=0, num_tokens, n, orig=0, i=0, status=0;
gchar **buff, *line;
struct core_pak *core;
GSList *clist=NULL;

g_assert(fp != NULL);

/* init */
if (data)
  {
  clist = data->cores;
  orig = g_slist_length(clist);
  }

/* process the 2 line XYZ header */
line = file_read_line(fp);
buff = tokenize(line, &num_tokens);
g_free(line);
if (num_tokens)
  {
  if (str_is_float(*buff))
    expect = str_to_float(*buff);
  }
g_strfreev(buff);

line = file_read_line(fp);
buff = tokenize(line, &num_tokens);
g_free(line);
if (!expect && num_tokens)
  {
  if (str_is_float(*buff))
    expect = str_to_float(*buff);
  }
g_strfreev(buff);

#if DEBUG_READ_XYZ_BLOCK
printf("Expect = %d\n", expect);
#endif

/* loop while there's data */
for (;;)
  {
  buff = get_tokenized_line(fp, &num_tokens);

  if (!buff)
    {
    status = 1;
    g_strfreev(buff);
    break;
    }

/* add new atom if we have sufficient tokens */
/* TODO - check the last 3 tokesn are valild floats? */
  if (num_tokens > 3)
    {
    n = elem_test(*buff);

/* FIXME - type 'X' returns 0 ... dummy atom type */
/*
    if (n)
*/
      {
      core = NULL;
      if (clist)
        {
        core = clist->data;
        clist = g_slist_next(clist);

        g_free(core->atom_label);
        core->atom_label = g_strdup(*buff);

        g_free(core->atom_type);
        core->atom_type = g_strdup(*buff);
        core->atom_code = n;

        elem_init(core, data);
        }
      else
        {
        if (data)
          {
          core = new_core(*buff, data);
          data->cores = g_slist_append(data->cores, core);
          core->render_mode = data->default_render_mode;
          }
        }
      if (core)
        {
        core->x[0] = str_to_float(*(buff+1));
        core->x[1] = str_to_float(*(buff+2));
        core->x[2] = str_to_float(*(buff+3));
        }
      i++;

/* NEW - end this frame if we've got the number of atoms expected */
      if (expect && i == expect)
        {
        g_strfreev(buff);
        break;
        }
      }



    }

/* quit if we've already read in some atoms */
/* ie end of xyz frame */
  g_strfreev(buff);
  }

#if DEBUG_READ_XYZ_BLOCK
printf("Found cores: %d\n", i);
#endif

/* If there are less cores than previous - delete leftover cores */
if (data)
  {
  if (i < orig)
    {

#if DEBUG_READ_XYZ_BLOCK
printf("Truncating new core list to match new input length: %d\n", i);
#endif

    for (clist=g_slist_nth(data->cores, i) ; clist ; clist=g_slist_next(clist))
      {
      delete_core(clist->data);
      }
    delete_commit(data);
    }
  }
return(status);
}

/*********************/
/* xyz frame read    */
/*********************/
gint read_xyz_frame(FILE *fp, struct model_pak *model)
{
read_xyz_block(fp, model);
return(0);
}

/****************/
/* file reading */
/****************/
#define DEBUG_READ_XYZ 0
gint read_xyz(gchar *filename, struct model_pak *data)
{
gint flag;
FILE *fp;

/* checks */
g_return_val_if_fail(data != NULL, 1);
g_return_val_if_fail(filename != NULL, 2);

fp = fopen(filename,"rt");
if (!fp)
  return(3);

/* loop while there's data */
flag=0;
data->num_frames = 0;

read_xyz_block(fp, data);

for (;;)
  {
  add_frame_offset(fp, data);

  if (read_xyz_block(fp, NULL))
    break;

  data->num_frames++;
  }

/* get rid of frame list if only one frame */
if (data->num_frames == 1)
  {
  free_list(data->frame_list);
  data->frame_list = NULL;
  }

/* model setup */
strcpy(data->filename, filename);
g_free(data->basename);
data->basename = parse_strip(filename);
model_prep(data);

return(0);
}

