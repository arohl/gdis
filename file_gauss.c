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

#include <stdio.h>
#include <string.h>

#include "gdis.h"
#include "coords.h"
#include "file.h"
#include "parse.h"
#include "matrix.h"
#include "model.h"
#include "interface.h"

/* main structures */
extern struct sysenv_pak sysenv;
extern struct elem_pak elements[];

/* global variables */
gint last_frame;

/************************/
/* GAUSSIAN file saving */
/************************/
gint write_gauss(gchar *filename, struct model_pak *model)
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

g_free(model->basename);
model->basename = parse_strip(filename);
fprintf(fp, "%%chk=%s.chk\n", model->basename);
fprintf(fp, "\n");
fprintf(fp, "# sp hf/sto-3g\n");
fprintf(fp, "\n");
gdis_blurb(fp);
if (model->title)
  fprintf(fp, "%s\n", model->title);
fprintf(fp, "\n");
fprintf(fp, "0 1\n");
for (list=model->cores ; list ; list=g_slist_next(list))
  {
  core = (struct core_pak *) list->data;
  if (core->status & DELETED)
    continue;
    
  /* everything is cartesian after latmat mult */
  ARR3SET(x, core->x);
  vecmat(model->latmat, x);
  fprintf(fp, "%-7s    %14.9f  %14.9f  %14.9f\n", elements[core->atom_code].symbol, x[0], x[1], x[2]);
  }
/* write out lattice vectors */
if (model->periodic > 0)
  for (i=0; i<model->periodic; i++)
    fprintf(fp, "%-7s    %14.9f  %14.9f  %14.9f\n", "Tv", model->latmat[0+i], model->latmat[3+i], model->latmat[6+i]);
fprintf(fp, "\n");  
fclose(fp);
return(0);
}

gchar **strip_empty_strings(gchar **str_array)
{
gchar **new_str_array = NULL;
gint i=0, j=0;

while (str_array[i] != NULL)
  {
  if (strlen(str_array[i]) != 0)
    {
    new_str_array = g_realloc(new_str_array, (j + 2)*sizeof(char *)); /* add one for NULL pointer */
    new_str_array[j] = str_array[i];
    j++;
    }
  else
    g_free(str_array[i]);
  i++;
  }
  new_str_array[j] = NULL;
  g_free(str_array);
  return(new_str_array);
}


gchar **read_gauss_keyword(gchar *string)
{
  gint i;
/* 
Options to keywords may be specified in any of the following forms:

            keyword = option  
            keyword(option) 
            keyword=(option1, option2, ...)  
            keyword(option1, option2, ...) 
*/

/* Spaces, tabs, commas, or forward slashes can be used in any combination to separate items within a line. Multiple spaces are treated as a single delimiter. */

for (i=0; i<strlen(string); i++)
  switch (string[i])
    {
    case '=':
    case '(':
    case ')':
    case '\t':
    case ',': string[i] = ' ';
    }
return(strip_empty_strings(g_strsplit(string, " ", 0)));
}


/*******************************************/
/* read single GAUSSIAN output configuration */
/*******************************************/
#define DEBUG_READ_GAUSSIAN_OUT 0
gint read_gauss_out_block(FILE *fp, struct model_pak *model)
{
gint i;
gint num_tokens;
gchar **buff, line[LINELEN], *text;
GSList *clist;
GString *energy_string, *grad_string;
struct core_pak *core;

clist = model->cores;

model->periodic = 0;
model->fractional = FALSE;
model->construct_pbc = FALSE;

/* read to end of iteration */
while (TRUE)
  {
  if (fgetline(fp, line))
    {
    gui_text_show(ERROR, "unexpected end of file reading to end of iteration\n");
    return(2);
    }
    
/*  if (g_ascii_strncasecmp(line, " GradGradGradGrad", 16) == 0)
      break;*/
 
/* read coordinates */
  if (g_strrstr(line, "Center     Atomic") != NULL)
    {
    for (i=0; i<3; i++)
      if (fgetline(fp, line))
        {
        gui_text_show(ERROR, "unexpected end of file skipping to coordinates\n");
        return(2);
        }

    buff = tokenize(line, &num_tokens);
    while (num_tokens == 6)
      {
      /* atomic number greater than 0 => real atom */
      if (str_to_float(*(buff+1)) > 0)
        {
        if (clist)
          {
          core = (struct core_pak *) clist->data;
          clist = g_slist_next(clist);
          }
        else
          {
          core = new_core(elements[(gint) str_to_float(*(buff+1))].symbol, model);
          model->cores = g_slist_append(model->cores, core);
          }
        core->x[0] = str_to_float(*(buff+3));
        core->x[1] = str_to_float(*(buff+4));
        core->x[2] = str_to_float(*(buff+5));
        #if DEBUG_READ_GAUSSIAN_OUT
        printf("new coords %f %f %f\n", core->x[0],  core->x[1], core->x[2]);
        #endif
        }
    /* atomic number = -2 => lattice vector */
      if (str_to_float(*(buff+1)) == -2)
        {
        model->latmat[0+model->periodic] = str_to_float(*(buff+3));
        model->latmat[3+model->periodic] = str_to_float(*(buff+4));
        model->latmat[6+model->periodic] = str_to_float(*(buff+5));
        model->periodic++;
        model->construct_pbc = TRUE;
        #if DEBUG_READ_GAUSSIAN_OUT
        printf(">> latmat: ");
        P3MAT(" ", model->latmat);
        #endif
        }
    /* get next line */
      g_strfreev(buff);
      if (fgetline(fp, line))
        {
        gui_text_show(ERROR, "unexpected end of file reading coordinates\n");
        return(2);
        }
      buff = tokenize(line, &num_tokens);
      }
    g_strfreev(buff);
    clist = model->cores;
    }

/* Energy */
  if (g_strrstr(line, "SCF Done:") != NULL)
    {
    buff = g_strsplit(line, "=", 2);
    text = format_value_and_units(g_ascii_strdown(*(buff+1), -1), 7);
    property_add_ranked(3, "Energy", text, model);
    g_free(text);
    g_strfreev(buff);
    }
  if (g_strrstr(line, "EUMP2") != NULL)
  	{
  	buff = g_strsplit(line, "=", 3);
    energy_string = g_string_new("");
    g_string_append_printf(energy_string, "%.7f %s", str_to_float(*(buff+2)), "a.u.");
    property_add_ranked(3, "Energy", energy_string->str, model);
    g_string_free(energy_string, TRUE); 
    g_strfreev(buff);
  	}
	/* TODO higher order MP */
  /* gradients */
  if (g_strrstr(line, "Internal  Forces:") != NULL)
    {
    buff = tokenize(line, &num_tokens);
    grad_string = g_string_new("");
    g_string_append_printf(grad_string, "%.7f %s", str_to_float(*(buff+5)), "a.u./Bohr");
    property_add_ranked(4, "RMS Gradient", grad_string->str, model);
    g_string_free(grad_string, TRUE); 
    g_strfreev(buff);
    }
    
  if (g_strrstr(line, "Predicted change in Energy") != NULL)
    {
    if (fgetline(fp, line))
      {
      gui_text_show(ERROR, "unexpected end of file checking for succesful optimization\n");
      return(2);
      }
      if (g_ascii_strcasecmp(property_lookup("Calculation", model), "Forces") == 0)
        last_frame = TRUE;
      if (g_strrstr(line, "Optimization completed") != NULL)
        last_frame = TRUE;
      /* fix me - better clue to stopping */
      break;
    }
    
  if (g_strrstr(line, "Job cpu time") != NULL)
    {
      last_frame = TRUE;
      /* fix me - better clue to stopping */
      break;
    }
  }    
return(0);
}

/*******************************/
/* GAUSSIAN output frame reading */
/*******************************/
gint read_gauss_out_frame(FILE *fp, struct model_pak *model)
{
    return(read_gauss_out_block(fp, model));
}

/********************************/
/* Read in a GAUSSIAN output file */
/********************************/
gint read_gauss_out(gchar *filename, struct model_pak *model)
{
gint frame, i, num_tokens;
gint have_scan=FALSE;
gchar **buff, **buff1, **tokens, line[LINELEN];
FILE *fp;

fp = fopen(filename, "rt");
if (!fp)
  return(1);

last_frame = FALSE;
frame=0;

/* defaults */
property_add_ranked(2, "Calculation", "Single Point", model);
property_add_ranked(6, "Method", "HF", model);
property_add_ranked(7, "Basis", "STO-3G", model);

while (!fgetline(fp, line))
  {
  /* find gaussian route section and store selected parameters in hash table */
  if (g_ascii_strncasecmp(line, " #", 2) == 0)
    {
    /* Spaces, tabs, commas, or forward slashes can be used in any combination to separate items within a line. Multiple spaces are treated as a single delimiter. */
    buff = tokenize((gchar *)(line) + 2, &num_tokens);
    
    for (i=0; i<num_tokens; i++)
      {
      /* handle iops */
      if (g_strrstr(buff[i], "iop") != NULL)
        continue;
      /* method and basis handled specially */
      /* hack to find method and basis but not general as one or other can be missing */
      if (g_strrstr(buff[i], "/") != NULL)
        {
        buff1 = g_strsplit(buff[i], "/", 2);
        property_add_ranked(6, "Method", buff1[0], model);
        property_add_ranked(7, "Basis", buff1[1], model);
        g_strfreev(buff1);
        continue;
        }
      tokens = read_gauss_keyword(buff[i]);
      if (g_ascii_strncasecmp(tokens[0], "Opt", 3) == 0)
        {
        property_add_ranked(2, "Calculation", "Optimisation", model);
        }
      else if (g_ascii_strncasecmp(tokens[0], "Force", 5) == 0)
        {
        property_add_ranked(2, "Calculation", "Forces", model);
        continue;
        }
      else if (g_ascii_strncasecmp(tokens[0], "SCRF", 4) == 0)
      	{
        property_add_ranked(6, "Solvation", "", model);
      	}
      g_strfreev(tokens);
      }
      g_strfreev(buff);
    }
   
  /* ModRedundant Section */
  if (g_ascii_strncasecmp(line, " The following ModRedundant input section has been read:", 56) == 0)
    {
    while (TRUE)
      {
      if (fgetline(fp, line))
        {
        gui_text_show(ERROR, "unexpected end of file reading ModRedundant input section\n");
        return(2);
        }
      if (g_ascii_strncasecmp(line, " Iteration", 10) == 0)
        break;
      buff = tokenize(line, &num_tokens);
      for (i=0; i<num_tokens; i++)
        {
        if (g_ascii_strncasecmp(buff[i], "S", 1) == 0)
          have_scan = TRUE;
        }
      g_strfreev(buff);  
      }
    }
    
  /* read coordinates */
	if ((g_strrstr(line, "Input orientation:") != NULL) || (g_strrstr(line, "Z-Matrix orientation:") != NULL) || (g_strrstr(line, "Standard orientation:") != NULL))
		{
		/* go through all frames to count them */
		add_frame_offset(fp, model);
		read_gauss_out_block(fp, model);
		frame++;
    if (!have_scan && last_frame)
      break;
    }
  }

/* read extra data from repeat of last coordinates here */    
  

/* done */
strcpy(model->filename, filename);
g_free(model->basename);
model->basename = parse_strip(filename);
model->num_frames = model->cur_frame = frame;
model->cur_frame--;

model_prep(model);

return(0);
}
