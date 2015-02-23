/*
Copyright (C) 2004 by Andrew Lloyd Rohl

andrew@ivec.org

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

#include "gdis.h"
#include "coords.h"
#include "file.h"
#include "parse.h"
#include "model.h"
#include "interface.h"

/* main structures */
extern struct sysenv_pak sysenv;
extern struct elem_pak elements[];

gint read_nw(gchar *filename, struct model_pak *data)
{
FILE *fp;
gchar line[LINELEN], *name;

fp = fopen(filename, "rt");
if (!fp)
  {
  sprintf(line, "Unable to open file %s\n", filename);
  gui_text_show(ERROR, line);
  return(1);
  }
else
  {
  name = g_path_get_basename(filename);
  sprintf(line, "Opening %s: \n", name);
  g_free(name);
  gui_text_show(STANDARD, line);
  }

strcpy(data->filename, filename);
g_free(data->basename);
data->basename = parse_strip(filename);
model_prep(data);
return(0);
}

/*******************************************/
/* read single NWChem output configuration */
/*******************************************/
#define DEBUG_READ_NWOUT 0
gint read_nwout_block(FILE *fp, struct model_pak *model)
{
gint num_tokens;
gchar **buff, line[LINELEN], *text;
GString *title, *gstring;
GSList *clist;
struct core_pak *core;

clist = model->cores;
  
/* skip until get a 1 in first column for first coordinate */
if (fgetline(fp, line))
  return(1);
while (g_ascii_strncasecmp(line, "    1", 5) != 0)
  {
  if (fgetline(fp, line))
    return(1);
  }

buff = tokenize(line, &num_tokens);
while (num_tokens > 0)
  {
  if (clist)
    {
    core = (struct core_pak *) clist->data;
    clist = g_slist_next(clist);
    }
  else
    {
    core = new_core(*(buff+1), model);
    model->cores = g_slist_append(model->cores, core);
    }

  core->x[0] = str_to_float(*(buff+3));
  core->x[1] = str_to_float(*(buff+4));
  core->x[2] = str_to_float(*(buff+5));
  #if DEBUG_READ_NWOUT
  printf("new coords %f %f %f\n", core->x[0],  core->x[1], core->x[2]);
  #endif

/* get next line */
  g_strfreev(buff);
  if (fgetline(fp, line))
    return(2);
  buff = tokenize(line, &num_tokens);
  }

g_strfreev(buff);

/* read until get the details of the current optimisation step */  
while (g_ascii_strncasecmp(line, "@", 1) != 0)
  {
  if (fgetline(fp, line))
    return(0);
  }

buff = tokenize(line, &num_tokens);
while (g_ascii_strncasecmp(line, "@", 1) == 0)
  {
  if (g_ascii_isdigit(buff[1][0]))
    {
    text = format_value_and_units(*(buff+2), 5);
    gstring = g_string_new(text);
    g_free(text);
    g_string_append(gstring, " a.u.");
    property_add_ranked(3, "Energy", gstring->str, model);
    g_string_free(gstring, TRUE);
    model->nwchem.energy = str_to_float(*(buff+2));
    
    model->nwchem.max_grad = str_to_float(*(buff+4));
    
    text = format_value_and_units(*(buff+5), 5);
    gstring = g_string_new(text);
    g_free(text);
    g_string_append(gstring, " a.u./A");
    property_add_ranked(4, "RMS Gradient", gstring->str, model);
    g_string_free(gstring, TRUE);
    
    model->nwchem.rms_grad = str_to_float(*(buff+5));
    model->nwchem.have_energy = model->nwchem.have_max_grad = model->nwchem.have_rms_grad = TRUE;
     }
  /* get next line */
  g_strfreev(buff);
  if (fgetline(fp, line))
    return(2);
  buff = tokenize(line, &num_tokens);
  }

/* has the optimiser finished successfully? Looking for 4 ok's */
if (num_tokens == 4)
  model->nwchem.min_ok = TRUE;

g_free(model->title);
title = g_string_new("");
g_string_append_printf(title, "E");
g_string_append_printf(title, "(DFT)");
g_string_append_printf(title, " = %.5f au", model->nwchem.energy);
g_string_append_printf(title, ", grad = %.5f", model->nwchem.rms_grad);
model->title = g_strdup(title->str);
g_string_free(title, TRUE); 

return(0);
}

/*******************************/
/* NWChem output frame reading */
/*******************************/
gint read_nwout_frame(FILE *fp, struct model_pak *model)
{
    return(read_nwout_block(fp, model));
}

/********************************/
/* Read in a NWChem output file */
/********************************/
gint read_nwout(gchar *filename, struct model_pak *model)
{
gint frame;
gchar line[LINELEN];
FILE *fp;

fp = fopen(filename, "rt");
if (!fp)
  return(1);
frame=0;
while (!fgetline(fp, line))
  {
	if (g_strrstr(line, "Geometry \"geometry\" -> \"geometry\"") != NULL)
		{
		/* go through all frames to count them */
		add_frame_offset(fp, model);
		read_nwout_block(fp, model);
		frame++;
    if (model->nwchem.min_ok)
      break;
    }
  }

/* done */
strcpy(model->filename, filename);
g_free(model->basename);
model->basename = parse_strip(filename);
model->num_frames = model->cur_frame = frame;
model->cur_frame--;

model_prep(model);

return(0);
}
