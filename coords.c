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
#include <strings.h>
#include <stdlib.h>
#include <math.h>

#include "gdis.h"
#include "coords.h"
#include "edit.h"
#include "error.h"
#include "interface.h"
#include "matrix.h"
#include "measure.h"
#include "spatial.h"
#include "surface.h"
#include "numeric.h"
#include "morph.h"
#include "opengl.h"
#include "render.h"
#include "select.h"
#include "zone.h"

/* data structures */
extern struct sysenv_pak sysenv;
extern struct elem_pak elements[];

/**********************/
/* debugging routines */
/**********************/
void dump_links(GSList *list)
{
GSList *item;
struct core_pak *core;
struct shel_pak *shell;

for (item=list ; item ; item=g_slist_next(item))
  {
  core = item->data;

  printf("core: %p, shell: %p ", core, core->shell);
  P3VEC(" : ", core->x);
  }
for (item=list ; item ; item=g_slist_next(item))
  {
  shell = item->data;

  printf("shell: %p, core: %p\n", shell, shell->core);
  P3VEC(" : ", shell->x);
  }
}

void print_core_list(GSList *list)
{
gint i=0;
GSList *item;
struct core_pak *core;
struct shel_pak *shell;

for (item=list ; item ; item=g_slist_next(item))
  {
  core = item->data;

/*
if (i==4)
*/
  {
  printf("[%4s] %1d %1d c %10.4f %10.4f %10.4f : %6.2f\n",
          core->atom_label, core->primary, core->orig,
          core->x[0], core->x[1], core->x[2], core->charge);

  if (core->shell)
    {
    shell = core->shell;

    printf("[%4s] %1d %1d s %10.4f %10.4f %10.4f : %6.2f\n",
            shell->shell_label, shell->primary, shell->orig,
            shell->x[0], shell->x[1], shell->x[2], shell->charge);
    }
  }
  i++;
  }
}

void print_core_list_cart(GSList *list)
{
    gint i=0;
    GSList *item;
    struct core_pak *core;
    struct shel_pak *shell;
    
    for (item=list ; item ; item=g_slist_next(item))
    {
        core = item->data;
        
        /*
         if (i==4)
         */
        {
            printf("[%4s] %1d %1d c %10.4f %10.4f %10.4f : %6.2f\n",
                   core->atom_label, core->primary, core->orig,
                   core->rx[0], core->rx[1], core->rx[2], core->charge);
            
            if (core->shell)
            {
                shell = core->shell;
                
                printf("[%4s] %1d %1d s %10.4f %10.4f %10.4f : %6.2f\n",
                       shell->shell_label, shell->primary, shell->orig,
                       shell->rx[0], shell->rx[1], shell->rx[2], shell->charge);
            }
        }
        i++;
    }
}

/****************/
/* SPECIAL OBJS */
/****************/
#define DEBUG_INIT 0
void coords_init(gint type, struct model_pak *data)
{
g_return_if_fail(data != NULL);

switch(type)
  {
/* atomic coordinates processing */
  case INIT_COORDS:
/* find unique elements */
    g_slist_free(data->unique_atom_list);
    data->unique_atom_list = find_unique(ELEMENT, data);
/* init the charges */
    init_model_charges(data);
/* update electrostatic info */
    calc_emp(data);
/* repetition here, since there is interdependence */
/* between centroid (coords_center) and latmat (coords_compute) calc */
    coords_center(data);
    coords_compute(data);
  case CENT_COORDS:
    coords_center(data);
  case REDO_COORDS:
    coords_compute(data);
    break;

  default:
    printf("Unknown object type requested!\n");
  }
}

/******************************************/
/* convert fractional coords to cartesian */
/******************************************/
void coords_make_cartesian(struct model_pak *model)
{
GSList *list;
struct core_pak *core;
struct shel_pak *shel;

g_assert(model != NULL);

/* transform cores */
for (list=model->cores ; list ; list=g_slist_next(list))
  {
  core = list->data;
  vecmat(model->latmat, core->x);
  }
/* transform shells */
for (list=model->shels ; list ; list=g_slist_next(list))
  {
  shel = list->data;
  vecmat(model->latmat, shel->x);
  }

matrix_identity(model->latmat);
matrix_identity(model->ilatmat);
matrix_identity(model->rlatmat);
}

/******************************************/
/* convert cartesian coords to fractional */
/******************************************/
void coords_make_fractional(struct model_pak *model)
{
GSList *list;
struct core_pak *core;
struct shel_pak *shel;

g_assert(model != NULL);

/* transform cores */
for (list=model->cores ; list ; list=g_slist_next(list))
  {
  core = list->data;
  vecmat(model->ilatmat, core->x);
  }
/* transform shells */
for (list=model->shels ; list ; list=g_slist_next(list))
  {
  shel = list->data;
  vecmat(model->ilatmat, shel->x);
  }
}

/*************************************/
/* process coordinate scaling factor */
/*************************************/
void coords_init_units(struct model_pak *model)
{
gdouble scale=1.0;
GSList *list;
struct core_pak *core;
struct shel_pak *shel;

g_assert(model != NULL);

switch (model->coord_units)
  {
  case BOHR:
    scale = AU2ANG;
    break;
  default:
    return;
  }

/* transform cores */
for (list=model->cores ; list ; list=g_slist_next(list))
  {
  core = list->data;
  VEC3MUL(core->x, scale);
  }

/* transform shells */
for (list=model->shels ; list ; list=g_slist_next(list))
  {
  shel = list->data;
  VEC3MUL(shel->x, scale);
  }
}

/********/
/* AXES */
/********/
void make_axes(struct model_pak *data)
{
VEC3SET(data->axes[0].x, 1.0, 0.0, 0.0);
VEC3SET(data->axes[1].x, 0.0, 1.0, 0.0);
VEC3SET(data->axes[2].x, 0.0, 0.0, 1.0);
}

/********/
/* CELL */
/********/
void make_cell(struct model_pak *data)
{
gdouble b, c;

if (!data->periodic)
  return;

/* periodicity hacks */
b = (data->periodic < 2) ? 0.0 : 1.0;
c = (data->periodic < 3) ? 0.0 : 1.0;

/* end face 1 */
VEC3SET(data->cell[0].x, 0.0, 0.0, 0.0);
VEC3SET(data->cell[1].x, 0.0, b, 0.0);
VEC3SET(data->cell[2].x, 1.0, b, 0.0);
VEC3SET(data->cell[3].x, 1.0, 0.0, 0.0);

/* end face 2 */
VEC3SET(data->cell[4].x, 0.0, 0.0, c);
VEC3SET(data->cell[5].x, 0.0, b, c);
VEC3SET(data->cell[6].x, 1.0, b, c);
VEC3SET(data->cell[7].x, 1.0, 0.0, c);
}

/*********************************/
/* construct core-shell linkages */
/*********************************/
#define DEBUG_SHELLS 0
#define MAX_SHELL_DIST 0.6
void shell_make_links(struct model_pak *model)
{
gint match;
gdouble min, sep, x[3], r[3];
gpointer zone;
GSList *list1, *list2, *locality, *floating;
struct core_pak *core;
struct shel_pak *shel;

/* checks */
g_assert(model != NULL);
g_assert(model->zone_array != NULL);

/* enumerate all shells */
floating = NULL;
for (list1=model->shels ; list1 ; list1=g_slist_next(list1))
  {
  shel = list1->data;

/* get neighbourhood core list for the shell */
  zone = zone_get(shel->x, model->zone_array);
  locality = zone_area_cores(1, zone, model->zone_array);

/* upper bound for core-shell link */
  min = MAX_SHELL_DIST*MAX_SHELL_DIST;

/* enumerate cores in the shell's vicinity */
  match = 0;
  for (list2=locality ; list2 ; list2=g_slist_next(list2))
    {
    core = list2->data;

    if (shel->atom_code == core->atom_code)
      {
/* get the minimum fractional separation */
      ARR3SET(x, shel->x);
      ARR3SUB(x, core->x);
      fractional_min(x, model->periodic);
/* convert to cartesian & compare with cutoff */
      ARR3SET(r, x);
      vecmat(model->latmat, r);
      sep = VEC3MAGSQ(r);
      if (sep < min)
        {
/* relocate shell to be close to core */
        ARR3SET(shel->x, core->x);
        ARR3ADD(shel->x, x);
/*
        ARR3SET(shel->rx, core->rx);
        ARR3ADD(shel->rx, r);
*/

/* FIXME - what to do with previous matches? eg set (core->shell)->shell = NULL ? */
/* reference each other */
        core->shell = shel;
        shel->core = core;
        match++;
/* keep track of the minimum, so only the closest pair is recorded */
        min = sep;
        }
      }
    }
/* record any shells with no linkages */
  if (!match)
    {
    floating = g_slist_prepend(floating, shel);
    error_table_entry("Error: found shell with no core.\n");
    }

/* free the constructed locality */
  g_slist_free(locality);
  }

/* found some floating shells? */
if (floating)
  {
/*
  text = g_strdup_printf("Warning: found %d shell(s) with no core(s).\n",
                         g_slist_length(floating));
  gui_text_show(WARNING, text);
  g_free(text);
*/

/* flag shells (diff colour?) with no cores */
  for (list1=floating ; list1 ; list1=g_slist_next(list1))
    {
    shel = list1->data;

/*
printf("Floating shell=%p\n", shel);
*/

    VEC3SET(shel->colour, 1.0, 1.0, 1.0);
    }
  g_slist_free(floating);
  }
}

/**********/
/* COORDS */
/**********/
#define DEBUG_CENT 0
gint coords_center(struct model_pak *data)
{
gint i, n;
gdouble vec[3], r, r2, rmax;
GSList *list=NULL, *ilist;
struct image_pak *image;
struct spatial_pak *spatial;
struct core_pak *core;
struct vec_pak *p1;

/* initialize */
VEC3SET(data->offset, 0, 0, 0);
VEC3SET(data->centroid, 0.0, 0.0, 0.0);
r2 = 0.0;

/* atom centroid calc */
n=0;
for (list=data->cores ; list ; list=g_slist_next(list))
  {
  core = list->data;
  if (core->status & DELETED)
    continue;

  ARR3ADD(data->centroid, core->x);
  n++;
  }
if (n)
  {
  VEC3MUL(data->centroid, 1.0 / (gdouble) n);
  }

/* always make a pbc/2 the centroid in periodic cases */
for (i=data->periodic ; i-- ; )
  data->centroid[i] = 0.5; 

/* adjust for images */
n=1;
VEC3SET(vec, 0.0, 0.0, 0.0);
for (list=data->images ; list ; list=g_slist_next(list))
  {
  image = list->data;

  ARR3ADD(vec, image->pic);
  n++;
  }
VEC3MUL(vec, 1.0 / (gdouble) n);
ARR3ADD(data->centroid, vec);

/* get distance of furtherest displayed item from centroid */
for (list=data->cores ; list ; list=g_slist_next(list))
  {
  core = list->data;
  if (core->status & DELETED)
    continue;

  ilist=NULL;
  do
    {
    ARR3SET(vec, core->x);
    if (ilist)
      {
      image = ilist->data;

/* image */
      ARR3ADD(vec, image->pic);
      ilist = g_slist_next(ilist);
      }
    else
      {
      ilist = data->images;
      }

    ARR3SUB(vec, data->centroid);
    vecmat(data->latmat,vec);

/* dist squared */
    r = VEC3MAGSQ(vec);
    if (r > r2)
      r2 = r; 
    }
  while (ilist);
  }

/* check cell vertices */
if (data->periodic)
  {
  for (i=8 ; i-- ; )
    {
    ARR3SET(vec, data->cell[i].x);
/* centroid removal */
    ARR3SUB(vec, data->centroid);
/* transform */
    vecmat(data->latmat, vec);  
/* dist sq. */
    r = VEC3MAGSQ(vec);
    if (r > r2)
      r2 = r; 
    }
  }

/* check spatial objects (includes morphology) */
for (list=data->spatial ; list ; list=g_slist_next(list))
  {
  spatial = list->data;
/* FIXME - better scheme? this only counts 1st vertex */
  p1 = g_slist_nth_data(spatial->list, 0);
  if (p1)
    {
    ARR3SET(vec, p1->x);
    ARR3SUB(vec, data->centroid);
    vecmat(data->latmat, vec);  
    r = VEC3MAGSQ(vec);
    if (r > r2)
      r2 = r;
    }
  }

/* check camera waypoints */
for (list=data->waypoint_list ; list ; list=g_slist_next(list))
  {
  struct camera_pak *camera = list->data;

  r = VEC3MAGSQ(camera->x);
  if (r > r2)
    r2 = r;
  }

/* assign the calculated value, unless 0 */
rmax = sqrt(r2);
if (rmax)
  data->rmax = rmax;
else
  data->rmax = 5.0;   /* ~ correct axes lengths (completely empty model) */

data->zoom = data->rmax;

#if DEBUG_CENT
P3VEC("centroid = ", data->centroid);
printf("rmax = %f\n", data->rmax);
printf("zoom = %f\n", data->zoom);
#endif

return(0);
}

/******************************************/
/* sort coords by z - mainly for surfaces */
/******************************************/
gint sort_cores(struct core_pak *core1, struct core_pak *core2)
{
if (core1->x[2] < core2->x[2])
  return(1);
if (core1->x[2] > core2->x[2])
  return(-1);
return(0);
}
gint sort_shels(struct shel_pak *shel1, struct shel_pak *shel2)
{
if (shel1->x[2] < shel2->x[2])
  return(1);
if (shel1->x[2] > shel2->x[2])
  return(-1);
return(0);
}
void sort_coords(struct model_pak *data)
{
data->cores = g_slist_sort(data->cores, (gpointer) sort_cores);
data->shels = g_slist_sort(data->shels, (gpointer) sort_shels);
}

/*************************/
/* get coordinate limits */
/*************************/
void cor_calc_xlimits(gdouble *min, gdouble *max, GSList *list)
{
gint i;
GSList *item;
struct core_pak *core;

/* no cores? */
if (!list)
  return;

/* init */
core = list->data;
for (i=3 ; i-- ; )
  min[i] = max[i] = core->x[i];

/* NB: even tho' we check the "fractional" part - the values */
/* will still be cartesian for isolated (or the appropriate */
/* mixed frac/cart values for 1D/2D) models */
/* FIXME - not always true - newly generated surfaces have z values */
/* in fractional coords */
/* loop to get limits */
for (item=list ; item ; item=g_slist_next(item))
  {
  core = item->data;

  for (i=3 ; i-- ; )
    {
    if (core->x[i] < min[i])
      min[i] = core->x[i];
    if (core->x[i] > max[i])
      max[i] = core->x[i];
    }
  }
}

/**********************************/
/* compute cartesian world coords */
/**********************************/
/* TODO - make more fine grained ie cores/shells or spatials (or ALL) */
#define DEBUG_UPDATE_COORDS 0
void coords_compute(struct model_pak *data)
{
gint n, ghost, model;
gdouble vec[3];
gdouble mat4[16], vec4[4];
GSList *list=NULL, *list1, *list2;
struct core_pak *core;
struct shel_pak *shell;
struct model_pak *orig;
struct vec_pak *vector;
struct vertex_pak *v;
struct plane_pak *plane;
struct ribbon_pak *ribbon;
struct spatial_pak *spatial;
struct object_pak *odata;
struct image_pak *image;
GSList *plist, *glist=NULL, *olist;

#if DEBUG_UPDATE_COORDS
printf("coords_compute() start.\n");
#endif

g_return_if_fail(data != NULL);

/* update model & any overlayed (ghost) models */
ghost=0;
orig = data;
while(data)
  {
#if DEBUG_UPDATE_COORDS
P3MAT("instantaneous rot matrix", rot);
P3MAT("cummulative rot matrix", data->rotmat);
P3MAT("           irot matrix", data->irotmat);
#endif

/* precalc matrix products */
ARR3SET(&mat4[0], &data->latmat[0]);
ARR3SET(&mat4[4], &data->latmat[3]);
ARR3SET(&mat4[8], &data->latmat[6]);
ARR3SET(vec, data->centroid);
vecmat(data->latmat, vec);
mat4[3] = -vec[0];
mat4[7] = -vec[1];
mat4[11] = -vec[2];
VEC4SET(&mat4[12], 0.0, 0.0, 0.0, 1.0);

#if DEBUG_UPDATE_COORDS
P4MAT("mat4:", mat4);
#endif

memcpy(data->display_lattice, mat4, 16*sizeof(gdouble));

/* update image translation vectors */
for (list=data->images ; list ; list=g_slist_next(list))
  {
  image = list->data;
  ARR3SET(vec4, image->pic);
  vec4[3] = 0.0;
  vec4mat(mat4, vec4);
  ARR3SET(image->rx, vec4);
  }

/* calculate for atoms */
for (list=data->cores ; list ; list=g_slist_next(list))
  {
  core = list->data;
  if (core->status & DELETED)
    continue;

  ARR3SET(vec4, core->x);
  ARR3ADD(vec4, core->offset);
  vec4[3] = 1.0;
  vec4mat(mat4, vec4);
  ARR3SET(core->rx, vec4);
  }

/* calculate for shells */
for (list=data->shels ; list ; list=g_slist_next(list))
  {
  shell = list->data;
  if (shell->status & DELETED)
    continue;

  ARR3SET(vec4, shell->x);
  ARR3ADD(vec4, shell->offset);
  vec4[3] = 1.0;
  vec4mat(mat4, vec4);
  ARR3SET(shell->rx, vec4);
  }

/* cell */
for (n=8 ; n-- ; )
  {
  ARR3SET(vec4, data->cell[n].x);
  vec4[3] = 1.0;
  vec4mat(mat4, vec4);
  ARR3SET(data->cell[n].rx, vec4);
  }

/* axes */
for (n=6 ; n-- ; )
  {
  ARR3SET(vec4, data->axes[n].x);
  vec4[3] = 0.0;
  if (data->axes_type != CARTESIAN)
    {
    vec4mat(mat4, vec4);
    normalize(vec4, 3);
    }
  ARR3SET(data->axes[n].rx, vec4);
  VEC3MUL(data->axes[n].rx, 0.06*data->rmax);
  }

/* ribbon updates */
olist = data->ribbons;
while (olist)
  {
  odata = olist->data;

  switch (odata->type)
    {
    case RIBBON:
      plist = (GSList *) odata->data;
      while (plist != NULL)
        {
        ribbon = plist->data;

/* end point 1 */
        ARR3SET(vec4, ribbon->x1);
        vec4[3] = 1.0;
        vec4mat(mat4, vec4); 
        ARR3SET(ribbon->r1, vec4);

/* end point 2 */
        ARR3SET(vec4, ribbon->x2);
        vec4[3] = 1.0;
        vec4mat(mat4, vec4); 
        ARR3SET(ribbon->r2, vec4);

/* normal 1 */
        ARR3SET(vec4, ribbon->u1);
        vec4[3] = 0.0;
        vec4mat(mat4, vec4);
        ARR3SET(ribbon->n1, vec4);

/* normal 2 */
        ARR3SET(vec4, ribbon->u2);
        vec4[3] = 0.0;
        vec4mat(mat4, vec4);
        ARR3SET(ribbon->n2, vec4);

/* orientation vector 1 */
        ARR3SET(vec4, ribbon->v1);
        vec4[3] = 0.0;
        vec4mat(mat4, vec4);
        ARR3SET(ribbon->o1, vec4);

/* orientation vector 2 */
        ARR3SET(vec4, ribbon->v2);
        vec4[3] = 0.0;
        vec4mat(mat4, vec4);
        ARR3SET(ribbon->o2, vec4);

        plist = g_slist_next(plist);
        }
    break;
    }
  olist = g_slist_next(olist);
  }

/* spatials */
for (olist=data->spatial ; olist ; olist=g_slist_next(olist))
  {
  spatial = olist->data;
  plist = spatial->list;
  while (plist)
    {
    vector = plist->data;

/* rotate each point */
    ARR3SET(vec4, vector->x);
    vec4[3] = 1.0;
    vec4mat(mat4, vec4);
    ARR3SET(vector->rx, vec4);

/* rotate the normal */
    ARR3SET(vec4, vector->n);
    vec4[3] = 0.0;
    vec4mat(mat4, vec4);
    normalize(vec4, 3);
    ARR3SET(vector->rn, vec4);

    plist = g_slist_next(plist);
    }
  }

/* compute facet center (miller label placement) */
  for (list1=data->planes ; list1 ; list1=g_slist_next(list1))
    {
    plane = list1->data;
    if (plane->present)
      {
      VEC3SET(vec, 0.0, 0.0, 0.0);
      n = 0; 
      for (list2=plane->vertices ; list2 ; list2=g_slist_next(list2))
        {
        v = list2->data;
        ARR3ADD(vec, v->rx);
        n++;
        }
/* assign */
/* NB: n=0 can sometimes occur with a bad (ie non-closed) polyhedron */
      if (n)
        {
        VEC3MUL(vec, 1.0/(gdouble) n);
        }
      ARR3SET(plane->rx, vec);
      }
    }

/* ghost processsing */
  if (!ghost)
    {
    glist = data->ghosts;
    ghost++;
    }
  else
    glist = g_slist_next(glist);
/* get the ghost model's pointer */
  if (glist)
    {
    model = GPOINTER_TO_INT(glist->data);
    data = model_ptr(model, RECALL);
    }
  else
    data = NULL;
  }

#if DEBUG_UPDATE_COORDS
printf("coords_compute() done.\n");
#endif
}

/***********************************/
/* print atom coord data primitive */
/***********************************/
void print_core(struct core_pak *core)
{
gchar *txt;

g_assert(core != NULL);

txt = g_strdup_printf("%4s core at (%9.4f,%9.4f,%9.4f).\n",
                      core->atom_label, core->x[0], core->x[1], core->x[2]);
gui_text_show(STANDARD, txt);
g_free(txt);
}

/***********************************/
/* print atom coord data primitive */
/***********************************/
void print_core_cart(struct core_pak *core)
{
    gchar *txt;
    
    g_assert(core != NULL);
    
    txt = g_strdup_printf("%4s core at (%9.4f,%9.4f,%9.4f).\n",
                          core->atom_label, core->rx[0], core->rx[1], core->rx[2]);
    gui_text_show(STANDARD, txt);
    g_free(txt);
}

/************************************/
/* print shell coord data primitive */
/************************************/
void print_shell(struct shel_pak *shell)
{
gchar *txt;

g_assert(shell != NULL);

txt = g_strdup_printf("%4s shel at (%9.4f,%9.4f,%9.4f) [%9.4f, %9.4f, %9.4f]\n",
                      shell->shell_label, shell->x[0], shell->x[1], shell->x[2],
                      shell->rx[0], shell->rx[1], shell->rx[2]);
gui_text_show(STANDARD, txt);
g_free(txt);
}

void print_cores(struct model_pak *data)
{
GSList *list;
struct core_pak *core;

for (list=data->cores ; list ; list=g_slist_next(list))
  {
  core = list->data;

  printf("%4s core at (%9.4f,%9.4f,%9.4f).\n",
         core->atom_label, core->x[0], core->x[1], core->x[2]);
  }
}
void print_cores_cart(struct model_pak *data)
{
    GSList *list;
    struct core_pak *core;
    
    for (list=data->cores ; list ; list=g_slist_next(list))
    {
        core = list->data;
        
        printf("%4s core at (%9.4f,%9.4f,%9.4f).\n",
               core->atom_label, core->rx[0], core->rx[1], core->rx[2]);
    }
}

void print_shells(struct model_pak *data)
{
GSList *list;
struct shel_pak *shel;

for (list=data->shels ; list ; list=g_slist_next(list))
  {
  shel = list->data;

  printf("%4s shel at (%9.4f,%9.4f,%9.4f).\n",
         shel->shell_label, shel->x[0], shel->x[1], shel->x[2]);
  }
}

/****************************************/
/* print all connected cores and shells */
/****************************************/
void print_core_shell(struct model_pak *data)
{
gint n=0;
GSList *list;
struct core_pak *core;
struct shel_pak *shel;

for (list=data->cores ; list ; list=g_slist_next(list))
  {
  core = list->data;
  if (core->shell)
    {
    shel = core->shell;

if (n == 7 || n == 15)
  {
    printf("(%p) %4s core at (%9.4f,%9.4f,%9.4f).\n",
           core, core->atom_label, core->x[0], core->x[1], core->x[2]);
    printf("(%p) %4s shel at (%9.4f,%9.4f,%9.4f).\n",
           shel, shel->shell_label, shel->x[0], shel->x[1], shel->x[2]);
  }

    n++;
    }
  }
}

/***********************/
/* core free primitive */
/***********************/
void core_free(gpointer data)
{
struct core_pak *core = data;

g_free(core->atom_label);
g_free(core->atom_type);
g_free(core->res_name);
g_free(core->flags);
g_free(core->chain);
g_slist_free(core->bonds);

free_slist(core->vibx_list);
free_slist(core->viby_list);
free_slist(core->vibz_list);

g_free(core);
}

/****************************************************/
/* free all lists that reference cores or core data */
/****************************************************/
void free_core_list(struct model_pak *data)
{
GSList *list;

g_assert(data != NULL);

/* free all objects that contain core refs first */

/* free shells */
/* TODO - free shell data */
free_slist(data->shels);
data->shels = NULL;

/* free connectivity */
free_slist(data->bonds);
data->bonds = NULL;
free_slist(data->ubonds);
data->ubonds = NULL;
free_mol_list(data);
data->moles = NULL;

/* others */
g_slist_free(data->selection);
data->selection = NULL;
/* NB: don't free data */
g_slist_free(data->unique_atom_list);
data->unique_atom_list = NULL;

/* free the cores */
for (list=data->cores ; list ; list=g_slist_next(list))
  core_free(list->data);
g_slist_free(data->cores);
data->cores = NULL;
}

