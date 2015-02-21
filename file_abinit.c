/*
Copyright (C) 2003 by Andrew Lloyd Rohl

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
#include "interface.h"

/* main structures */
extern struct sysenv_pak sysenv;
extern struct elem_pak elements[];

/* globals */
gint optcell;

/*******************************************/
/* read single ABINIT output configuration */
/*******************************************/
#define DEBUG_READ_ABOUT 0
gint read_about_block(FILE *fp, struct model_pak *model)
{
gint i, num_tokens;
gchar **buff, line[LINELEN];
gdouble acell[3];
GString *title;
GSList *clist;
struct core_pak *core;

clist = model->cores;
if (optcell > 0)
  {
  /* go to acell line and read it */
  if (fgetline(fp, line))
    {
    gui_text_show(ERROR, "unexpected end of file reading cell dimensions in animation\n");
    return(2);
    }
  buff = tokenize(line, &num_tokens);
  acell[0] = str_to_float(*(buff+1)) * AU2ANG;
  acell[1] = str_to_float(*(buff+2)) * AU2ANG;
  acell[2] = str_to_float(*(buff+3)) * AU2ANG;
  g_strfreev(buff);

  /* get the cell vectors i.e. rprim.acell */
  /* NB: gdis wants transposed ABINIT matrix */
  for (i=0; i<3; i++)
    {
    if (fgetline(fp, line))
      {
      gui_text_show(ERROR, "unexpected end of file reading cell dimensions\n");
      return(2);
      }
    buff = tokenize(line + 8, &num_tokens);
    model->latmat[0+i] = acell[i] * str_to_float(*(buff+0));
    model->latmat[3+i] = acell[i] * str_to_float(*(buff+1));
    model->latmat[6+i] = acell[i] * str_to_float(*(buff+2));
    g_strfreev(buff);
    }

  for (i=0; i<4; i++)
    {
    if (fgetline(fp, line))
      {
      gui_text_show(ERROR, "unexpected end of file skipping to coordinates in animation\n");
      return(2);
      }
    }
  }
  
/* get 1st line of coords */
if (fgetline(fp, line))
  return(1);
buff = tokenize(line, &num_tokens);

while (g_ascii_strncasecmp(line, " Cartesian forces", 17) != 0)
  {
  if (clist)
    {
    core = (struct core_pak *) clist->data;
    clist = g_slist_next(clist);
    }
  else
    {
    core = new_core(*(buff+4), model);
    model->cores = g_slist_append(model->cores, core);
    }

  core->x[0] = AU2ANG*str_to_float(*(buff+0));
  core->x[1] = AU2ANG*str_to_float(*(buff+1));
  core->x[2] = AU2ANG*str_to_float(*(buff+2));
  #if DEBUG_READ_ABOUT
  printf("new coords %f %f %f\n", core->x[0],  core->x[1], core->x[2]);
  #endif

/* get next line */
  g_strfreev(buff);
  if (fgetline(fp, line))
    return(2);
  buff = tokenize(line, &num_tokens);
  }
  
/* get gradients */
model->abinit.max_grad = str_to_float(*(buff+4));
model->abinit.rms_grad = str_to_float(*(buff+5));
g_strfreev(buff);

/* get energy */
if (fgetline(fp, line))
  return(2);
while (g_ascii_strncasecmp(line, " At the end of Broyden", 22) != 0)
  {
  if (fgetline(fp, line))
    return(2);
  }
buff = tokenize(line, &num_tokens);
model->abinit.energy = str_to_float(*(buff+9));

title = g_string_new("");
g_string_append_printf(title, "E");
g_string_append_printf(title, " = %.4f Ha, ", model->abinit.energy);
g_string_append_printf(title, "max grad = %.5f", model->abinit.max_grad);
model->title = g_strdup(title->str);
g_string_free(title, TRUE);
return(0);
}

/*******************************/
/* ABINIT output frame reading */
/*******************************/
gint read_about_frame(FILE *fp, struct model_pak *model)
{
    return(read_about_block(fp, model));
}

/********************************/
/* Read in a ABINIT output file */
/********************************/
gint read_about(gchar *filename, struct model_pak *model)
{
gint i, flag, frame, num_tokens, tot_tokens;
gint natom=0, ntype=0, *types;
GString *title;
gchar **buff, **curr_token, *ptr, line[LINELEN];
struct core_pak *core;
GSList *clist;
FILE *fp;

fp = fopen(filename, "rt");
if (!fp)
  return(1);
flag=frame=0;
while (!fgetline(fp, line))
  {
/* space group info */
  if (g_ascii_strncasecmp(" Symmetries :", line, 13) == 0)
    {
    ptr = g_strrstr(line, "(#");
    model->sginfo.spacenum = (gint) str_to_float(ptr+2);
    }

/* default cell dimensions */
/* NB: gdis wants transposed ABINIT matrix */
    if (g_ascii_strncasecmp(" Real(R)+Recip(G)", line, 17) == 0)
      {
      for (i=0; i<3; i++)
        {
        if (fgetline(fp, line))
          {
          gui_text_show(ERROR, "unexpected end of file reading cell dimensions\n");
          return(2);
          }
        buff = tokenize(line, &num_tokens);
        model->latmat[0+i] = str_to_float(*(buff+1))*AU2ANG;
        model->latmat[3+i] = str_to_float(*(buff+2))*AU2ANG;
        model->latmat[6+i] = str_to_float(*(buff+3))*AU2ANG;
        g_strfreev(buff);
        }
      model->construct_pbc = TRUE;
      model->periodic = 3;
      /* last part of data in output file before minimisation */
      break;
      }

/* optimisation type */  
    if (g_ascii_strncasecmp("   optcell", line, 10) == 0)
      {
      buff = tokenize(line, &num_tokens);
      optcell = (gint) str_to_float(*(buff+1));
      g_strfreev(buff);
      }
/* coordinates */
    if (g_ascii_strncasecmp("     natom", line, 10) == 0)
      {
      buff = tokenize(line, &num_tokens);
      natom = (gint) str_to_float(*(buff+1));
      g_strfreev(buff);
      }
    if (g_ascii_strncasecmp("     ntype", line, 10) == 0)
      {
      buff = tokenize(line, &num_tokens);
      ntype = (gint) str_to_float(*(buff+1));
      g_strfreev(buff);
      }
    if (g_ascii_strncasecmp("      type", line, 10) == 0)
      {
      tot_tokens = 0;
      buff = tokenize(line, &num_tokens);
      curr_token = buff + 1;
      while (tot_tokens < natom)
        {
        if (*curr_token == NULL)
          {
          g_strfreev(buff);
          buff = get_tokenized_line(fp, &num_tokens);
          curr_token = buff;
          }
        core = new_core("X", model);
        model->cores = g_slist_prepend(model->cores, core);
        core->atom_code = (gint) str_to_float(*curr_token) - 1;
        tot_tokens++;
        curr_token++;
        }
      model->cores = g_slist_reverse(model->cores);
      g_strfreev(buff);
      flag++;
      }
    if (g_ascii_strncasecmp("    xangst", line, 10) == 0)
      {
      clist = model->cores;
      if (clist == NULL)
        {
        gui_text_show(ERROR, "no atoms found\n");
        return(2);
        }
      buff = tokenize(line+10, &num_tokens);
      for (i=0; i<natom; i++)
        {
        core = (struct core_pak *) clist->data;
        core->x[0] = str_to_float(*(buff+0));
        core->x[1] = str_to_float(*(buff+1));
        core->x[2] = str_to_float(*(buff+2));
        g_strfreev(buff);
        buff = get_tokenized_line(fp, &num_tokens);
        clist = g_slist_next(clist);
        }
      }
    if (g_ascii_strncasecmp("     znucl", line, 10) == 0)
      {
      types = g_malloc(ntype * sizeof(gint));
      i = 0;
      tot_tokens = 0;
      buff = tokenize(line, &num_tokens);
      curr_token = buff + 1;
      while (tot_tokens < ntype)
        {
        if (*curr_token == NULL)
          {
          g_strfreev(buff);
          buff = get_tokenized_line(fp, &num_tokens);
          curr_token = buff;
          }
        types[i] = (gint) str_to_float(*curr_token);
        i++;
        tot_tokens++;
        curr_token++;
        }
      g_strfreev(buff);
      /* now update cores to get right element type */
      clist = model->cores;
      if (clist == NULL)
        {
        gui_text_show(ERROR, "no atoms found\n");
        return(2);
        }
      for (i=0; i<natom; i++)
        {
        core = clist->data;
        core_init(elements[types[core->atom_code]].symbol, core, model);
        clist = g_slist_next(clist);
        }
      }

  }

/* Additional info */
while (!fgetline(fp, line))
  {
  /* energy */
  if (g_ascii_strncasecmp(line, " Ewald  =", 9) == 0)
    {
    printf("reading energy\n");
    model->abinit.energy = str_to_float(&line[51]);
    printf("energy is %lf\n", model->abinit.energy);
    }
  
  /* gradient */
  else if ((g_ascii_strncasecmp(line, " frms,max,avg=", 14) == 0) && g_strrstr(line, "h/b"))
    {
    printf("reading gradient\n");
    buff = tokenize(line, &num_tokens);
    model->abinit.max_grad = str_to_float(*(buff+2));
    model->abinit.rms_grad = str_to_float(*(buff+1));
    g_strfreev(buff);
    printf("gradient is %lf\n", model->abinit.max_grad);
    }

  /* coordinates */
  else if ((g_ascii_strncasecmp(line, " Cartesian coordinates (bohr)", 29) == 0 && optcell == 0) || (g_ascii_strncasecmp(line, " Unit cell characteristics :", 28) == 0 && optcell > 0))
    {
    /* go through all frames to count them */
    add_frame_offset(fp, model);
    read_about_block(fp, model);
    frame++;
    }
  }

/* done */
if (flag)
  {
  strcpy(model->filename, filename);
  g_free(model->basename);
  model->basename = parse_strip(filename);
  model->num_frames = model->cur_frame = frame;
  model->cur_frame--;

  title = g_string_new("");
  g_string_append_printf(title, "E");
  g_string_append_printf(title, " = %.4f Ha, ", model->abinit.energy);
  g_string_append_printf(title, "max grad = %.5f", model->abinit.max_grad);
  model->title = g_strdup(title->str);
  g_string_free(title, TRUE);

  model_prep(model);
  }
else
  return(2);

return(0);
}

