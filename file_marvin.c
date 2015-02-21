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
#include "matrix.h"
#include "model.h"
#include "space.h"
#include "scan.h"
#include "interface.h"

#define DEBUG_MORE 0
#define MAX_KEYS 15

/* main structures */
extern struct sysenv_pak sysenv;
extern struct elem_pak elements[];

/*******************************/
/* MARVIN specific file saving */
/*******************************/
gint write_marvin(gchar *filename, struct model_pak *model)
{
gint mol, region;
gdouble vec[3];
GSList *list;
struct core_pak *core;
struct shel_pak *shel;
FILE *fp;

/* open the file */
fp = fopen(filename,"w");
if (!fp)
  return(1);

gui_text_show(STANDARD, "Saving file in MARVIN restart format!\n");

/* save via the number of molecules in the */
/* model - since this allows mol numbering */
mol=0;

/* four region types printing */
for (region=REGION1A ; region<REGION2B ; region++)
  {
  if (model->region_empty[region])
    continue;
  switch(region)
    {
    case REGION1A:
      fprintf(fp,"coordinates 1 A\n");
      break;
    case REGION2A:
      fprintf(fp,"coordinates 2 A\n");
      break;
    case REGION1B:
      fprintf(fp,"coordinates 1 B\n");
      break;
    case REGION2B:
      fprintf(fp,"coordinates 2 B\n");
      break;
    }

/* do atoms */
  for (list=model->cores ; list ; list=g_slist_next(list))
    {
    core = (struct core_pak *) list->data;

/* don't print checks */
    if (core->region != region)
      continue;
    if (core->status & DELETED)
      continue;

/* NB: x,y,z have been fractionalized */
    ARR3SET(vec, core->x);
    vecmat(model->latmat, vec);

/* NEW - compute molecule number */
if (core->bonds)
  mol = g_slist_index(model->moles, core->mol);
else
  mol = 0;

    fprintf(fp,"%-7s   core   %14.9f  %14.9f  %14.9f  %4d  %14.9f\n",
               core->atom_label, vec[0], vec[1], vec[2], mol, core->charge);

/* shell */
    if (core->shell)
      {
      shel = (struct shel_pak *) core->shell;

      ARR3SET(vec, shel->x);
      vecmat(model->latmat, vec);
      fprintf(fp,"%-7s   shel   %14.9f  %14.9f  %14.9f  %4d  %14.9f\n",
             shel->shell_label, vec[0], vec[1], vec[2], mol, shel->charge);
      }
    }
/* done region */
  fprintf(fp,"end\n");
  }

/* done */
fclose(fp);
return(0);
}

/****************************/
/* marvin output file parse */
/****************************/
#define DEBUG_LOAD_MVNOUT 0
gint read_mvnout(gchar *filename, struct model_pak *model)
{
gboolean overwrite=FALSE;
gint i, j, region, offset, num_tokens;
gint count[REGION2B];
gchar *line, **buff, *text;
gdouble veca[2], vecb[2];
gdouble version;
gpointer scan;
GSList *core_list=NULL, *shell_list=NULL;
struct core_pak *core;
struct shel_pak *shell;

/* checks */
g_assert(model != NULL);
g_assert(filename != NULL);

scan = scan_new(filename);
if (!scan)
  return(1);

/* setup */
model->id = MVNOUT;
model->periodic = 2;
strcpy(model->filename, filename);
g_free(model->basename);
model->basename = parse_strip(filename);

for (i=REGION2B ; i-- ; )
  count[i] = 0;

/* NB: disallow image creation in this direction */
model->pbc[2] = 0.0;

/* scan the file */
offset=0;
line = scan_get_line(scan);
while (!scan_complete(scan))
  {
  buff = tokenize(line, &num_tokens);

/* YAMH - get version to determine a parsing parameter */
  if (g_ascii_strncasecmp(line, "                   Version", 27) == 0)
    {
    if (num_tokens > 1)
      {
      version = str_to_float(*(buff+1));
      offset = (version < 1.99) ? 0 : 1;
#if DEBUG_LOAD_MVNOUT
printf("Marvin version %lf, offset = %d\n", version, offset);
#endif 
      }
    }

/* get surface vectors */
  if (g_ascii_strncasecmp(line,"Surface Vectors:",16) == 0)
    {
/* read them in */
    if (num_tokens > 3)
      {
      veca[0] = str_to_float(*(buff+2));
      veca[1] = str_to_float(*(buff+3));
      g_strfreev(buff);
      buff = scan_get_tokens(scan, &num_tokens);
      if (num_tokens > 1)
        {
        vecb[0] = str_to_float(*buff);
        vecb[1] = str_to_float(*(buff+1));
/* NB: gdis stores vectors in reverse order */
        VEC3SET(&(model->latmat[0]), veca[0], vecb[0], 0.0);
        VEC3SET(&(model->latmat[3]), veca[1], vecb[1], 0.0);
        VEC3SET(&(model->latmat[6]), 0.0, 0.0, 1.0);
        model->construct_pbc = TRUE;
        model->periodic = 2;
        }
      }
    }

/* energy seaches */
  if (g_ascii_strncasecmp(line, "Total Energy", 12) == 0)
    {
    g_strfreev(buff);
    buff = g_strsplit(line, "=", 2);
    text = format_value_and_units(*(buff+1), 5);
    property_add_ranked(3, "Energy", text, model);
    g_free(text);
    }

  if (g_ascii_strncasecmp(line, "Attachment Energy", 17) == 0)
    {
    if (num_tokens > 4)
      {
      if (model->gulp.eatt[0] == 0.0)
        model->gulp.eatt[0] = str_to_float(*(buff+3));
      else
        model->gulp.eatt[1] = str_to_float(*(buff+3));
      g_free(model->gulp.eatt_units);
      model->gulp.eatt_units = g_strdup(*(buff+4));
      }
    }

  if (g_ascii_strncasecmp(line, "Surface Energy", 14) == 0)
    {
    if (num_tokens > 5)
      {
      if (model->gulp.esurf[0] == 0.0)
        model->gulp.esurf[0] = str_to_float(*(buff+3));
      else
        model->gulp.esurf[1] = str_to_float(*(buff+3));
      g_free(model->gulp.esurf_units);
      model->gulp.esurf_units = g_strdup(*(buff+5));
      }
    }

  if (g_ascii_strncasecmp(line, "Gnorm", 5) == 0)
    {
    g_strfreev(buff);
    buff = g_strsplit(line, "=", 2);
    text = format_value_and_units(*(buff+1), 5);
    property_add_ranked(6, "Gnorm", text, model);
    g_free(text);
    }

/* coordinate data search */
  if (num_tokens > 5 && g_ascii_strncasecmp(line, "BLOCK", 5) == 0)
    {
    if (g_ascii_strncasecmp(*(buff+2),"region",5) == 0)
      {

/* get region type */
    region = REGION1A;
    if (g_ascii_strncasecmp(*(buff+5),"1b",2) == 0)
      region = REGION1B;
    if (g_ascii_strncasecmp(*(buff+5),"2a",2) == 0)
      region = REGION2A;
    if (g_ascii_strncasecmp(*(buff+5),"2b",2) == 0)
      region = REGION2B;

/* check if it's coordinate data (don't want gradients) */
    line = scan_get_line(scan);

    g_strfreev(buff);
    buff = scan_get_tokens(scan, &num_tokens);

    if (num_tokens > 3)
    if (g_ascii_strncasecmp(*(buff+3),"coordinates",11) == 0)
      {
      model->region_empty[region] = FALSE;

/* only do this once as we don't want to overwrite from */
/* the start of the core/shell list at every new block */
      if (count[region] && !overwrite)
        {
#if DEBUG_LOAD_MVNOUT
printf("Repeat found - overwriting from now on.\n");
#endif
        core_list = model->cores = g_slist_reverse(model->cores);
        shell_list = model->shels = g_slist_reverse(model->shels);
        overwrite = TRUE;
        }

#if DEBUG_LOAD_MVNOUT
printf("Reading region: %d\n", region);
#endif

/* skip to data */
      line = scan_get_line(scan);
      line = scan_get_line(scan);

      while (!scan_complete(scan))
        {
        if (g_ascii_strncasecmp(line,"-------",7) == 0)
          break;

        g_strfreev(buff);
        buff = tokenize(line, &num_tokens);

        if (num_tokens < 7)
          break;

/* NEW - overwrite if 2nd occurance and existing, else create new */
/* NB: marvin 1.8 & 2.0 have this at different positions! */
      if (g_ascii_strncasecmp(*(buff+offset+2),"cor",3) == 0)
        {
        if (overwrite && core_list)
          {
          core = core_list->data;
          core_list = g_slist_next(core_list);
          }
        else
          {
          core = new_core(*(buff+offset+1), model);
          model->cores = g_slist_prepend(model->cores, core);
          core->charge = str_to_float(*(buff+offset+6));
          core->lookup_charge = FALSE;
          core->region = region;
          }
        core->x[0] = str_to_float(*(buff+offset+3));
        core->x[1] = str_to_float(*(buff+offset+4));
        core->x[2] = str_to_float(*(buff+offset+5));
        }

      if (g_ascii_strncasecmp(*(buff+offset+2),"she",3) == 0)
        {
        if (overwrite && shell_list)
          {
          shell = shell_list->data;
          shell_list = g_slist_next(shell_list);
          }
        else
          {
          shell = new_shell(*(buff+offset+1), model);
          model->shels = g_slist_prepend(model->shels, shell);
          shell->charge = str_to_float(*(buff+offset+6));
          shell->lookup_charge = FALSE;
          shell->region = region;
          }
        shell->x[0] = str_to_float(*(buff+offset+3));
        shell->x[1] = str_to_float(*(buff+offset+4));
        shell->x[2] = str_to_float(*(buff+offset+5));
        }
      line = scan_get_line(scan);
      }
      count[region]++;
      }
    }
    }

  line = scan_get_line(scan);
  g_strfreev(buff);
  }

/* init lighting */
j=0;
for (i=REGION1A ; i<REGION2B ; i++)
  {
  if (!model->region_empty[i] && !j)
    {
    model->region[i] = TRUE;
    j++;
    }
  else
    model->region[i] = FALSE;
  }
   
#if DEBUG_LOAD_MVNOUT
printf("Surface vectors: %lf x %lf\n",model->pbc[0],model->pbc[1]);
#endif

/* force non periodic in z (marvin is 2D) */
model->pbc[2] = 0;
model_prep(model);

scan_free(scan);

return(0);
}

/****************/
/* file reading */
/****************/
#define DEBUG_READ_MARVIN 0
gint read_marvin(gchar *filename, struct model_pak *model)
{
gint i, num_tokens, keyword, region;
gchar **buff, *tmp, line[LINELEN];
gdouble lat_mult = 1.0;
GSList *list;
struct core_pak *core;
struct shel_pak *shel;
FILE *fp;

/* checks */
g_return_val_if_fail(model != NULL, 1);
g_return_val_if_fail(filename != NULL, 1);

fp = fopen(filename, "rt");
if (!fp)
  return(1);

/* defaults */
/*
model->periodic = 2;
*/

model->region_empty[REGION1A] = TRUE;
model->region_empty[REGION2A] = TRUE;
model->region_empty[REGION1B] = TRUE;
model->region_empty[REGION2B] = TRUE;

keyword = COMPARE;
region = REGION1A;

while(!fgetline(fp, line))
  {
  list = get_keywords(line);
  if (list != NULL)
    {
    keyword = GPOINTER_TO_INT(list->data);
    switch(keyword)
      {
      case MARVIN_COORD:
        tmp = g_strstrip(g_strdup(get_token_pos(line, 1)));

        if (g_ascii_strncasecmp("1 a", tmp, 3) == 0)
          region = REGION1A;
        if (g_ascii_strncasecmp("2 a", tmp, 3) == 0)
          region = REGION2A;
        if (g_ascii_strncasecmp("1 b", tmp, 3) == 0)
          region = REGION1B;
        if (g_ascii_strncasecmp("2 b", tmp, 3) == 0)
          region = REGION2B;

        model->region_empty[region] = FALSE;

        g_free(tmp);
        break;
      }
    }
  else
    {
    buff = tokenize(line, &num_tokens);
/* have we got a valid element to read */
    if (buff)
      {
    if (elem_symbol_test(*buff))
      {
      switch(keyword)
        {
        case MARVIN_COORD:
        case MARVIN_BASIS:
          if (num_tokens > 4)
            {
            if (g_ascii_strncasecmp("cor", *(buff+1), 3) == 0)
              {
              core = new_core(*buff, model);
              model->cores = g_slist_prepend(model->cores, core);

              core->region = region;
              core->x[0] = str_to_float(*(buff+2));
              core->x[1] = str_to_float(*(buff+3));
              core->x[2] = str_to_float(*(buff+4));
              if (num_tokens > 6)
                {
                core->charge = str_to_float(*(buff+6));
                core->lookup_charge = FALSE;
                }
              }
            if (g_ascii_strncasecmp("shel", *(buff+1), 4) == 0)
              {
              shel = new_shell(*buff, model);
              model->shels = g_slist_prepend(model->shels, shel);

              shel->region = region;
              shel->x[0] = str_to_float(*(buff+2));
              shel->x[1] = str_to_float(*(buff+3));
              shel->x[2] = str_to_float(*(buff+4));
              if (num_tokens > 6)
                {
                shel->charge = str_to_float(*(buff+6));
                shel->lookup_charge = FALSE;
                }
              }
            }
          break;
        }
      }
    else
      {
/* Added by C. Fisher 2004 */
      if (g_ascii_strncasecmp(line, "lattice", 7) == 0)
        {
        buff = tokenize(line, &num_tokens);
        lat_mult = str_to_float(*(buff+1));
        model->fractional = TRUE;
        }
      if (g_ascii_strncasecmp(line, "latvec", 6) == 0)
        {
        model->periodic = 2; 
        if (!fgetline(fp, line))
          {
          g_strfreev(buff);
          buff = tokenize(line, &num_tokens);
          model->latmat[0] = str_to_float(*(buff+0));
          model->latmat[1] = str_to_float(*(buff+1));
          model->latmat[2] = str_to_float(*(buff+2));
          }
        if (!fgetline(fp, line))
          {
          g_strfreev(buff);
          buff = tokenize(line, &num_tokens);
          model->latmat[3] = str_to_float(*(buff+0));
          model->latmat[4] = str_to_float(*(buff+1));
          model->latmat[5] = str_to_float(*(buff+2));
          }
        if (!fgetline(fp, line))
          {
          g_strfreev(buff);
          buff = tokenize(line, &num_tokens);
          model->latmat[6] = str_to_float(*(buff+0));
          model->latmat[7] = str_to_float(*(buff+1));
          model->latmat[8] = str_to_float(*(buff+2));
          }
        model->construct_pbc = TRUE;
        }

      /* Added by C. Fisher 2004 */
      if (g_ascii_strncasecmp(line, "surfvec", 7) == 0)
        {
        if (!fgetline(fp, line))
          {
          g_strfreev(buff);
          buff = tokenize(line, &num_tokens);
          model->image_limit[1] = str_to_float(*(buff+0));
          }
        if (!fgetline(fp, line))
          {
          g_strfreev(buff);
          buff = tokenize(line, &num_tokens);
          model->image_limit[3] = str_to_float(*(buff+1));
          }
        }

      }
      }
    g_strfreev(buff);
    }
  }

/* Added by C. Fisher 2004 */
for(i=0;i<9;i++)
  model->latmat[i] *= lat_mult;

/* model setup */
strcpy(model->filename, filename);
g_free(model->basename);
model->basename = parse_strip(filename);

model_prep(model);

/* Added by C. Fisher 2004 */
if( model->image_limit[1] != 1.0 ||  model->image_limit[3] != 1.0)
  {
  space_make_images(CREATE, model);
  coords_init(CENT_COORDS, model);
  }

fclose(fp);
return(0);
}

