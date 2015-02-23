/*
Copyright (C) 2000 by Sean David Fleming

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
#include "matrix.h"
#include "interface.h"
#include "parse.h"

/* main structures */
extern struct sysenv_pak sysenv;
extern struct elem_pak elements[];

/**********************/
/* save in CSSR format */
/**********************/
gint write_cssr(gchar *filename, struct model_pak *data)
{
gint i;
gchar *label;
gdouble vec[3], charge;
GSList *list;
FILE *fp;
struct core_pak *core;

/* checks */
g_return_val_if_fail(data != NULL, 1);
g_return_val_if_fail(filename != NULL, 2);

/* open the file */
fp = fopen(filename,"w");
if (!fp)
  {
  gui_text_show(ERROR, "Bad filename!\n");
  return(1);
  }
if( data->periodic == 2 )
  {
  gui_text_show(ERROR, "Cannot save surface structure!\n");
  return(1);
  }
if( data->num_atoms > 9999)
  puts("Too many atoms for CSSR format: no save");
else
  gui_text_show(STANDARD, "Saving file in CSSR format!\n");

/* header */
if (data->periodic == 3)
  {
  /* Match sginfo list order with cssr options */
  if( data->sginfo.spacenum == 9 || data->sginfo.spacenum == 15)
    {
    if( data->sginfo.cellchoice > 3 && data->sginfo.cellchoice < 7 )
       data->sginfo.cellchoice += 6;
    else if( data->sginfo.cellchoice > 6 && data->sginfo.cellchoice < 10 )
       data->sginfo.cellchoice -= 3;
    else if( data->sginfo.cellchoice > 9 && data->sginfo.cellchoice < 13 )
       data->sginfo.cellchoice += 3;
    else if( data->sginfo.cellchoice > 12 && data->sginfo.cellchoice < 16 )
       data->sginfo.cellchoice -= 6;
    }
  if( data->sginfo.spacename)
    label = g_strdup(data->sginfo.spacename);
  else
    label = g_strdup("");
  fprintf(fp,"%38c %7.3f %7.3f %7.3f\n",' ',data->pbc[0],data->pbc[1],data->pbc[2]);
  fprintf(fp,"%21c %7.3f %7.3f %7.3f    SPGR =%3d %-11s", ' ',R2D*data->pbc[3],
         R2D*data->pbc[4], R2D*data->pbc[5],data->sginfo.spacenum, label);
  i=0;
  if( data->sginfo.cellchoice > 0 )
    fprintf(fp," OPT = %d",data->sginfo.cellchoice);
  fprintf(fp,"\n");
  }
else
  {
  fprintf(fp,"\n\n");
  i=1;
  }
fprintf(fp,"%4d   %d Created by GDIS\n", data->num_asym, i);
fprintf(fp,"\n");

/* coords */
i=0;
for (list=data->cores ; list ; list=g_slist_next(list))
  {
  core = (struct core_pak *) list->data;
                                                                                         
  if (core->status & DELETED)
    continue;
  if (core->primary)
    {
/* retrieve */
    ARR3SET(vec, core->x);
    i++;

  charge = atom_charge(core);

/* transformation needed? */
    fprintf(fp,"%4d %-4s  %9.5f %9.5f %9.5f",
              i, core->atom_label, vec[0], vec[1], vec[2]);
    fprintf(fp,"    0   0   0   0   0   0   0   0 %7.3f\n",
              charge);
    }
  }

/* done */
fclose(fp);
return(0);
}

/****************/
/* read routine */
/****************/
#define DEBUG_READ_CSSR 0
gint read_cssr(gchar *filename, struct model_pak *data)
{
gint n=0, num_tokens, sgopt;
gchar **buff, line[LINELEN];
gchar *spacename;
struct core_pak *core;
FILE *fp;

/* checks */
g_return_val_if_fail(data != NULL, 1);
g_return_val_if_fail(filename != NULL, 1);

fp = fopen(filename, "rt");
if (!fp)
  return(1);

/* periodicity search */
  fgetline(fp, line);
  g_strstrip(line);
  if (g_ascii_strcasecmp("", line) != 0)
    {
    buff = tokenize(line+38, &num_tokens);
    if (num_tokens == 3 )
      data->periodic = 3;
    data->pbc[0] = str_to_float(*(buff));
    data->pbc[1] = str_to_float(*(buff+1));
    data->pbc[2] = str_to_float(*(buff+2));
    g_strfreev(buff);
    }
  else
    data->periodic = 0;

  fgetline(fp, line);
  if (g_ascii_strcasecmp("", line) != 0)
    {
/* seek space group label */
    if( data->periodic )
      {
      spacename = g_strdup("");
      if( sscanf(line+23,"%lf%8lf%8lf    SPGR =%3d %11[a-zA-Z0-9-/ ] OPT =%2d",
        &data->pbc[3], &data->pbc[4], &data->pbc[5], &data->sginfo.spacenum,
          spacename, &data->sginfo.cellchoice) < 4)
             printf("Error reading cssr file %s\n", filename);
      data->pbc[3] *= D2R; data->pbc[4] *= D2R; data->pbc[5] *= D2R;
      g_strstrip(spacename);
      data->sginfo.spacename = g_strdup(spacename);
      }

/* If no space group given, or unrecognizable, use number */
    if( !data->sginfo.spacename || g_ascii_strcasecmp("",data->sginfo.spacename) == 0)
      {
      if( data->sginfo.spacenum > 1)
        {
        g_free(data->sginfo.spacename);
        data->sginfo.spacename = g_strdup_printf("%d",data->sginfo.spacenum);
        }
      else
        data->sginfo.spacename = g_strdup("P 1");
      }

    if( data->sginfo.cellchoice > 1)
      {
      /* Convert cssr options to axis choice for simple cases */
      if( data->sginfo.spacenum == 3 || data->sginfo.spacenum == 4 ||
          data->sginfo.spacenum == 6 || data->sginfo.spacenum == 10 ||
          data->sginfo.spacenum == 11 )
        {
        switch( data->sginfo.cellchoice )
          {
          case 1:
            data->sginfo.spacename = g_strconcat(data->sginfo.spacename, ":b", NULL);
            break;
          case 2:
            data->sginfo.spacename = g_strconcat(data->sginfo.spacename, ":c", NULL);
            break;
          case 3:
            data->sginfo.spacename = g_strconcat(data->sginfo.spacename, ":a", NULL);
            break;
          }
        }
      /* Convert cssr options to axis choice and origin choice for other cases */
      else if( data->sginfo.spacenum == 5 || data->sginfo.spacenum == 7 ||
          data->sginfo.spacenum == 8 || data->sginfo.spacenum == 12 ||
          data->sginfo.spacenum == 13 || data->sginfo.spacenum == 14 )
        {
          if( data->sginfo.cellchoice < 4)
            data->sginfo.spacename = g_strconcat(data->sginfo.spacename,":b", NULL);
          else if( data->sginfo.cellchoice > 3 && data->sginfo.cellchoice < 7)
            data->sginfo.spacename = g_strconcat(data->sginfo.spacename,":c", NULL);
          else if( data->sginfo.cellchoice > 6 && data->sginfo.cellchoice < 10)
            data->sginfo.spacename = g_strconcat(data->sginfo.spacename,":a", NULL);
          sgopt = data->sginfo.cellchoice % 3;
          if(!sgopt)
             sgopt = 3;
          sprintf(data->sginfo.spacename, "%s%d", data->sginfo.spacename, sgopt);
        }
      /* Match cssr options to sginfo name for two special cases */
      else if( data->sginfo.spacenum == 9 || data->sginfo.spacenum == 15)
       {
         if( data->sginfo.cellchoice < 4)
           data->sginfo.spacename = g_strconcat(data->sginfo.spacename,":b", NULL);
         else if( data->sginfo.cellchoice > 3 && data->sginfo.cellchoice < 7)
           data->sginfo.spacename = g_strconcat(data->sginfo.spacename,":c", NULL);
         else if( data->sginfo.cellchoice > 6 && data->sginfo.cellchoice < 10)
           data->sginfo.spacename = g_strconcat(data->sginfo.spacename,":a", NULL);
         else if( data->sginfo.cellchoice > 9 && data->sginfo.cellchoice < 13)
           data->sginfo.spacename = g_strconcat(data->sginfo.spacename,":-b", NULL);
         else if( data->sginfo.cellchoice > 12 && data->sginfo.cellchoice < 16)
           data->sginfo.spacename = g_strconcat(data->sginfo.spacename,":-c", NULL);
         else if( data->sginfo.cellchoice > 15 && data->sginfo.cellchoice < 19)
           data->sginfo.spacename = g_strconcat(data->sginfo.spacename,":-a", NULL);
          sgopt = data->sginfo.cellchoice % 3;
          if(!sgopt)
             sgopt = 3;
          sprintf(data->sginfo.spacename, "%s%d", data->sginfo.spacename, sgopt);
       }
      /* Adjust options to match name for SGs 50, 59, 68 */
      else if( data->sginfo.spacenum == 50 || data->sginfo.spacenum == 59 ||
               data->sginfo.spacenum == 68 )
        {
           if( data->sginfo.cellchoice % 2 == 0)
             data->sginfo.spacename = g_strconcat(data->sginfo.spacename, ":2", NULL);
        }
      /* Rhombohedral space groups */
      else if((data->sginfo.spacenum == 146 || data->sginfo.spacenum == 148 ||
               data->sginfo.spacenum == 155 ||
               data->sginfo.spacenum == 160 || data->sginfo.spacenum == 161 ||
               data->sginfo.spacenum == 166 || data->sginfo.spacenum == 167) &&
               data->sginfo.cellchoice > 1 )
                   data->sginfo.spacename = g_strconcat(data->sginfo.spacename, ":r", NULL);
      /* Groups with possible alternative origins */
      else if( data->sginfo.spacenum == 70 ||
              (data->sginfo.spacenum > 84 && data->sginfo.spacenum < 143) ||
               data->sginfo.spacenum > 201)
         sprintf(data->sginfo.spacename, "%s:%d", data->sginfo.spacename, data->sginfo.cellchoice);
    }

#if DEBUG_READ_CSSR
printf("Space group %d:%d %s\n",data->sginfo.spacenum, data->sginfo.cellchoice, data->sginfo.spacename);
#endif

    }
  else
    if( data->periodic == 3 )
      {
      gui_text_show(ERROR, "Incorrect format for CSSR file\n");
      fclose(fp);
      return(1);
      }

  buff = get_tokenized_line(fp, &num_tokens);
  if (g_ascii_strncasecmp("0", *(buff+1), 1) == 0 )
    data->fractional = TRUE;
  else
    data->fractional = FALSE;
  g_strfreev(buff);

/* coordinate search */
  for (;;)
    {
    buff = get_tokenized_line(fp, &num_tokens);
    if (!buff)
      break;
    if (num_tokens > 4)
      {
/* add new atom if first item is a valid atom type */
      if (elem_symbol_test(*(buff+1)))
        {
        core = new_core(*(buff+1), data);
        data->cores = g_slist_prepend(data->cores, core);
                                                                                         
        core->x[0] = str_to_float(*(buff+2));
        core->x[1] = str_to_float(*(buff+3));
        core->x[2] = str_to_float(*(buff+4));

        if (num_tokens > 12)
          {
          core->charge = str_to_float(*(buff+13));
          core->lookup_charge = FALSE;
          }
#if DEBUG_READ_CSSR
printf("Atom %d %f %f %f Charge %f\n",n, core->x[0],
        core->x[1], core->x[2], core->charge);
#endif
        n++;
        }
      }
    g_strfreev(buff);
    }

/* got everything */
data->num_asym = data->num_atoms = n;

#if DEBUG_READ_CSSR
printf("Found %d atoms.\n", n);
#endif

/* model setup */
strcpy(data->filename, filename);
g_free(data->basename);
data->basename = parse_strip(filename);

model_prep(data);
fclose(fp);

return(0);
}
