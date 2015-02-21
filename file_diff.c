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
#include "matrix.h"
#include "surface.h"
#include "interface.h"


/***********************/
/* DIFFaX read routine */
/***********************/
#define DEBUG_READ_DIFFAX 0
gint read_diffax(gchar *filename, struct model_pak *model)
{
gint num_tokens, num_layer, tot_layer;
gdouble offset;
gchar **buff;
GSList *list1, *list2;
struct core_pak *core;
struct layer_pak *layer;
FILE *fp;

/* checks */
g_return_val_if_fail(model != NULL, 1);
g_return_val_if_fail(filename != NULL, 2);
fp = fopen(filename, "rt");
if (!fp)
  return(3);

/* setup */
model->id = DIFFAX_INP;
model->fractional = TRUE;
model->periodic = 3;
model->colour_scheme = REGION;

strcpy(model->filename, filename);
g_free(model->basename);
model->basename = parse_strip(filename);

/* scan the file */
while ((buff = get_tokenized_line(fp, &num_tokens)))
  {
  diffax_keyword_search:;

/* restricted unit cell */
  if (g_ascii_strncasecmp("structural", *buff, 10) == 0)
    {
    g_strfreev(buff);
    buff = get_tokenized_line(fp, &num_tokens);
    if (num_tokens > 3)
      {
      model->pbc[0] = str_to_float(*(buff+0));
      model->pbc[1] = str_to_float(*(buff+1));
      model->pbc[2] = str_to_float(*(buff+2));
      model->pbc[3] = PI/2.0;
      model->pbc[4] = PI/2.0;
      model->pbc[5] = D2R*str_to_float(*(buff+3));
      }
    }

/* layer testing */
  if (g_ascii_strncasecmp("layer", *buff, 5) == 0)
    {
    layer = g_malloc(sizeof(struct model_pak));
    layer->width = 1.0;
    VEC3SET(layer->centroid, 0.5, 0.5, 0.5);
    layer->cores = NULL;
    model->layer_list = g_slist_prepend(model->layer_list, layer);

    g_strfreev(buff);
    buff = get_tokenized_line(fp, &num_tokens);
    if (buff)
      {
/* TODO - if centrosymmetric : add a -1 operation */
      }

/* get layer data */
    g_strfreev(buff);
    buff = get_tokenized_line(fp, &num_tokens);
    while (buff)
      {
      if (elem_symbol_test(*buff))
        {
        if (num_tokens > 6)
          {
/*
printf("[%s] [%s  %s  %s]\n", *buff, *(buff+2), *(buff+3), *(buff+4));
*/
          core = new_core(*buff, model);
          model->cores = g_slist_prepend(model->cores, core);
          layer->cores = g_slist_prepend(layer->cores, core);

          core->x[0] = str_to_float(*(buff+2));
          core->x[1] = str_to_float(*(buff+3));
          core->x[2] = str_to_float(*(buff+4));

          core->sof = str_to_float(*(buff+5));
          }
        }
      else
        goto diffax_keyword_search;

/* get next line of tokens */
      g_strfreev(buff);
      buff = get_tokenized_line(fp, &num_tokens);
      }
    }

  g_strfreev(buff);
  }

/* TODO - enumerate layers and scale, so they are stacked 1..n in a single cell */
/* also label the layers as different region types */
model->layer_list = g_slist_reverse(model->layer_list);
num_layer = 0;
tot_layer = g_slist_length(model->layer_list);

model->pbc[2] *= tot_layer;

#if DEBUG_READ_DIFFAX
printf("Read in %d layers.\n", tot_layer);
#endif

for (list1=model->layer_list ; list1 ; list1=g_slist_next(list1))
  {
  layer = (struct layer_pak *) list1->data;
  layer->width = 1.0 / (gdouble) tot_layer;

  offset = (gdouble) num_layer * layer->width;

  VEC3SET(layer->centroid, 0.0, 0.0, offset);

  for (list2=layer->cores ; list2 ; list2=g_slist_next(list2))
    {
    core = (struct core_pak *) list2->data;

/* scale to within the big cell (encloses all DIFFAX layers) */
    core->x[2] *= layer->width;

/* offset each particular layer */
    core->x[2] += offset;

    core->region = num_layer;
    }
  num_layer++;
  }

/* end of read */
fclose(fp);

/* post read setup */
model->cores = g_slist_reverse(model->cores);
model_prep(model);

return(0);
}

/***********************************/
/* DIFFaX input file write routine */
/***********************************/
gint write_diffax(gchar *filename, struct model_pak *data)
{
gint i, j, n, tot_layer;
gdouble pr, x[3];
GSList *list1, *list2;
struct layer_pak *layer;
struct core_pak *core;
struct elem_pak elem;
FILE *fp;

printf("write_diffax()\n");

/* checks */
g_return_val_if_fail(data != NULL, 1);
g_return_val_if_fail(filename != NULL, 2);

/* open the file */
fp = fopen(filename,"wt");
if (!fp)
  return(3);

/* setup the layers */
#ifdef WITH_GUI
diffract_layer_setup(data);
#endif
tot_layer = g_slist_length(data->layer_list);
if (!tot_layer)
  {
  printf("No layers to write.\n");
  return(0);
  }

/* output */
fprintf(fp, "\n");

fprintf(fp, "INSTRUMENTAL\nX-RAY\n1.5418\n");
fprintf(fp, "GAUSSIAN 0.4\n\n");

fprintf(fp, "STRUCTURAL\n%11.6f %11.6f %11.6f %11.6f\nUNKNOWN\n%d\ninfinite\n\n", 
            data->pbc[0], data->pbc[1], data->pbc[2]/tot_layer, R2D*data->pbc[5],
            tot_layer);

/* write out the layers */
i=1;
for (list1=data->layer_list ; list1 ; list1=g_slist_next(list1))
  {
  fprintf(fp, "\nLAYER %d\nnone\n", i);
/* write out the atoms within the layer */
  n=1;
  layer = (struct layer_pak *) list1->data;
  for (list2=layer->cores ; list2 ; list2=g_slist_next(list2))
    {
    core = (struct core_pak *) list2->data;
    if (core->status & DELETED)
      continue;

    ARR3SET(x, core->x);

/* center layer on origin */
    x[2] -= layer->centroid[2];
/* scale up to 0.0-1.0 range */
    x[2] *= tot_layer;
/* move centroid to middle of the layer */
    x[2] += 0.5;

    get_elem_data(core->atom_code, &elem, NULL);

    fprintf(fp, "%4s %4d  %11.6f %11.6f %11.6f  1.0 %8.6f\n",
                 elem.symbol, n, x[0], x[1], x[2], core->sof);

    n++;
    }
  i++;
  }

fprintf(fp, "\nSTACKING\nrecursive\ninfinite\n\n");

fprintf(fp, "TRANSITIONS\n");

/* omit region 0 - temporary impl for the later method where */
/* the desired layers will be in a linked list, rather than */
/* all 1..max regions included */
for (i=1 ; i<=tot_layer ; i++)
  {
  fprintf(fp, "\n");

/* uniform probability & no xlat */
  pr = 1.0 / (gdouble) tot_layer;
 
  for (j=1 ; j<=tot_layer ; j++)
    fprintf(fp, "%8.6f  0.0 0.0 1.0    {layer %d-%d}\n",pr,i,j);
  }

fclose(fp);
return(0);
}

/*******************************/
/* write a diffax control file */
/*******************************/
gint write_diffax_control(gchar *filename, struct model_pak *mdata)
{
FILE *fp;

/* open control */
fp = fopen("control.dat","wt");
if (!fp)
  return(3);

/* TODO - proper values from mdata */
fprintf(fp, "%s          {data file}\n", filename);
fprintf(fp, "%d          {dump file}\n", 0);
fprintf(fp, "%d          {symmetry file}\n", 0);
fprintf(fp, "%d          {pattern}\n", 3);
fprintf(fp, "%f %f %f    {theta min, max, delta}\n", 0.0, 90.0, 0.1);
fprintf(fp, "%d          {}\n", 1);
fprintf(fp, "%d          {}\n", 1);

fclose(fp);
return(0);
}

