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

#include <math.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#ifdef __APPLE__
#include <OpenGL/gl.h>
#else
#include <GL/gl.h>
#endif
#include "gdis.h"
#include "coords.h"
#include "model.h"
#include "file.h"
#include "parse.h"
#include "scan.h"
#include "matrix.h"
#include "spatial.h"
#include "interface.h"

/* main structures */
extern struct sysenv_pak sysenv;
extern struct elem_pak elements[];

/****************/
/* file reading */
/****************/
#define DEBUG_READ_OFF 0
gint read_off(gchar *filename, struct model_pak *model)
{
gint a, b, c, i, n, v, num_tokens, remaining, v_limit;
gint colour_flag, normal_flag, texture_flag;
gchar **buff, *text;
gdouble norm[3];
gpointer scan;
struct spatial_pak *spatial;
struct vec_pak **vec, *vec1;

/* checks */
g_return_val_if_fail(model != NULL, 1);

scan = scan_new(filename);
if (!scan)
  return(2);

/* setup flags */
colour_flag = normal_flag = texture_flag = FALSE;

/* parse header */
buff = scan_get_tokens(scan, &num_tokens);
if (num_tokens)
  {
  text = *buff;

n = strlen(text);
for (i=0 ; i<n ; i++)
  {
  switch (text[i])
    {
    case 'n':
      scan_get_line(scan); 
      break;

    case 'C':
      colour_flag = TRUE;
      break;

    case 'N':
      normal_flag = TRUE;
      break;

    case 'T':
      texture_flag = TRUE;
      break;
    }
  }
  }
else
  return(3);

g_strfreev(buff);
buff = scan_get_tokens(scan, &num_tokens);
v_limit = str_to_float(*buff);

/*
printf("expected vertices = %d [%d]\n", v_limit, normal_flag);
*/

g_assert(v_limit > 0);

g_strfreev(buff);
buff = scan_get_tokens(scan, &num_tokens);
vec = g_malloc(v_limit * sizeof(struct vec_pak *));

v=0;
while (!scan_complete(scan))
  {
  vec[v] = g_malloc(sizeof(struct vec_pak));

  i = 0; 
  remaining = num_tokens;

  if (remaining > 2)
    {
    vec[v]->x[0] = str_to_float(*(buff+i++));
    vec[v]->x[1] = str_to_float(*(buff+i++));
    vec[v]->x[2] = str_to_float(*(buff+i++));
    vecmat(model->ilatmat, vec[v]->x);
    remaining -= 3;
    }

/* TODO - if no normal - calculate? */
  if (normal_flag)
    {
    if (remaining > 2)
      {
      vec[v]->n[0] = str_to_float(*(buff+i++));
      vec[v]->n[1] = str_to_float(*(buff+i++));
      vec[v]->n[2] = str_to_float(*(buff+i++));
      }
    remaining -= 3;
    }

/* FIXME - can get quite tricky here with colour map indices etc. */
  if (colour_flag)
    {
/* rgb */
    if (remaining > 2)
      {
      vec[v]->colour[0] = str_to_float(*(buff+i++));
      vec[v]->colour[1] = str_to_float(*(buff+i++));
      vec[v]->colour[2] = str_to_float(*(buff+i++));

      if (VEC3MAGSQ(vec[v]->colour) < FRACTION_TOLERANCE)
        {
/* black -> GDIS's re-entrant surface colour */
        ARR3SET(vec[v]->colour, sysenv.render.rsurf_colour);
        }

      }
    remaining -= 3;
/* alpha */
/*
    if (remaining)
      vec[v]->colour[3] = str_to_float(*(buff+i++));
sysenv.render.transmit = str_to_float(*(buff+i++));
    remaining--;
*/
    }
  else
    {
/* default (no colour specified) */
    ARR3SET(vec[v]->colour, sysenv.render.rsurf_colour);
    }

  if (++v == v_limit)
    break;

  g_strfreev(buff);
  buff = scan_get_tokens(scan, &num_tokens);
  }

g_strfreev(buff);
buff = scan_get_tokens(scan, &num_tokens);

while (!scan_complete(scan))
  {
  if (num_tokens)
    {
    n = str_to_float(*buff);
    g_assert(num_tokens > n);

/* only create spatial for 3 or more vertices */
/* TODO - some files have lines with only 2 vertices - create line */
if (n > 2)
  {

/* TODO - only create one spatial - add vertices to this */
/* TODO - if no normals - use method that makes OpenGL calculate */
    spatial = spatial_new(NULL, SPATIAL_GENERIC, 0, TRUE, model);

/* check vertices have the correct ordering */
/* might not be necessary - geomview might already do this */
    a = str_to_float(*(buff+1));
    b = str_to_float(*(buff+2));
    c = str_to_float(*(buff+3));
/* TODO - convert to cartesian (general case) */
    calc_norm(norm, vec[a]->x, vec[b]->x, vec[c]->x);
    normalize(norm, 3);

/* fill out the vertices */
    for (i=0 ; i<n ; i++)
      {
      v = str_to_float(*(buff+i+1)); 

      g_assert(v >= 0);
      g_assert(v < v_limit);

/* copy vector, as normals may change depending on current polygon */
      vec1 = g_malloc(sizeof(struct vec_pak));
      ARR3SET(vec1->x, vec[v]->x);
      if (normal_flag)
        {
        ARR3SET(vec1->n, vec[v]->n);
        }
      else
        {
        ARR3SET(vec1->n, norm);
        }
      ARR3SET(vec1->colour, vec[v]->colour);

      spatial->list = g_slist_prepend(spatial->list, vec1);
      }

    spatial->list = g_slist_reverse(spatial->list);
  }


    }

  g_strfreev(buff);
  buff = scan_get_tokens(scan, &num_tokens);
  }

/* cleanup */
g_strfreev(buff);
for (i=v_limit ; i-- ; )
  g_free(*(vec+i));
g_free(vec);

/* model setup */
strcpy(model->filename, filename);
g_free(model->basename);
model->basename = parse_strip(filename);
model->id = GEOMVIEW_OFF;
model->fractional = FALSE;

model_prep(model);

return(0);
}

