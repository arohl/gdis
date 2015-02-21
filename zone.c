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
#ifdef __APPLE__
#include <OpenGL/gl.h>
#else
#include <GL/gl.h>
#endif

#include "gdis.h"
#include "coords.h"
#include "interface.h"
#include "gui_shorts.h"
#include "matrix.h"
#include "spatial.h"
#include "numeric.h"
#include "morph.h"
#include "opengl.h"
#include "select.h"
#include "zone.h"

/* data structures */
extern struct sysenv_pak sysenv;
extern struct elem_pak elements[];

/***********************************/
/* spatial partitioning structures */
/***********************************/

/* zone data */
struct zone_pak
{
gboolean visible;
gint grid[3];
GSList *cores;
GSList *shells;
};

/* zone array */
struct zone_array_pak
{
gdouble min[3];   /* min coord */
gdouble max[3];   /* max coord */
gdouble idx[3];   /* inverse dx - so can mult (not divide) to get index */
gint div[3];      /* divisions */
gint periodic;    /* follows model periodicity */
gint size;        /* num zones */
struct zone_pak **zone;
};

/********************************/
/* free space partitioning data */
/********************************/
void zone_free(gpointer data)
{
gint i;
struct zone_pak *zone;
struct zone_array_pak *za=data;

if (!za)
  return;

/* free zone array core lists */
for (i=za->size ; i-- ; )
  {
  zone = za->zone[i];
  g_slist_free(zone->cores);
  g_slist_free(zone->shells);
  g_free(zone);
  }

/* free the zone array */
g_free(za->zone);
g_free(za);
}

/**********************************/
/* construct a space partitioning */
/**********************************/
/* after this distance - connectivity is no longer considered (worst case) */
#define DEBUG_ZONE_MAKE 0
gpointer zone_make(gdouble zone_size, struct model_pak *model)
{
gint i, j, k, zi, num_zones;
gint div[3];
gdouble tmp, dx[3], min[3], max[3];
GSList *list;
struct core_pak *core;
struct shel_pak *shel;
struct zone_pak *zone;
struct zone_array_pak *za;

g_assert(model != NULL);

#if DEBUG_ZONE_MAKE
printf("zone_make(%s) [%d D]\n", model->filename, model->periodic);
#endif

za = g_malloc(sizeof(struct zone_array_pak));

/* FIXME - will there be a problem for periodic models */
/* if atoms are not constrained ie in the range [0-1) */
/* TODO - allow zones "outside" pbc, but when comparing -> clamp to within */

/* get the (possibly mixed) frac/cart coord limits */
VEC3SET(min, 0.0, 0.0, 0.0);
VEC3SET(max, 0.0, 0.0, 0.0);
cor_calc_xlimits(min, max, model->cores);

/* set number of divisions - CARTESIAN PART ONLY */
for (i=model->periodic ; i<3 ; i++)
  {
/* safety margin */
  max[i] += 0.1;
  min[i] -= 0.1;
  tmp = dx[i] = max[i] - min[i];
  tmp /= zone_size;
  div[i] = nearest_int(tmp);
  }

/* set number of divisions - FRACTIONAL PART ONLY */
for (i=model->periodic ; i-- ; )
  {
  min[i] = 0.0;
  max[i] = 1.0-G_MINDOUBLE;
  dx[i] = 1.0-G_MINDOUBLE;
  tmp = model->pbc[i] / zone_size;
  div[i] = nearest_int(tmp);
  }

/* division sizes */
for (i=3; i--; )
  {
  if (!div[i])
    div[i] = 1;
  dx[i] /= (gdouble) div[i];
  }

/* TODO - if periodic - extend zone by 1 in all directions (unfragment) */
/* ie div +=2 , min -= dx, max += dx */
/* TODO - OR (better) if we can just constrain to get the equivalent existing zone */
ARR3SET(za->min, min);
ARR3SET(za->max, max);
ARR3SET(za->div, div);
for (i=3 ; i-- ; )
  za->idx[i] = 1.0/dx[i];
za->periodic = model->periodic;
num_zones = div[0] * div[1] * div[2];
za->size = num_zones;

#if DEBUG_ZONE_MAKE
for (i=0 ; i<3 ; i++)
  printf("[%d] : %11.4f - %11.4f (div = %d)(dx=%f)(idx = %f)\n",
            i, za->min[i], za->max[i], za->div[i], dx[i], za->idx[i]);
printf("Allocating %dx%dx%d = %d zones...\n", div[0], div[1], div[2], num_zones);
#endif

za->zone = g_malloc(num_zones * sizeof(struct zone_pak *));

#if DEBUG_ZONE_MAKE
printf("Initializing zone array: %p, size : %d\n", za, num_zones);
#endif

for (k=0 ; k<div[2] ; k++)
  {
  for (j=0 ; j<div[1] ; j++)
    {
    for (i=0 ; i<div[0] ; i++)
      {
      zi = k*div[1]*div[0]; 
      zi += j*div[0]; 
      zi += i;

      za->zone[zi] = g_malloc(sizeof(struct zone_pak));

      zone = za->zone[zi];
      zone->cores = NULL;
      zone->shells = NULL;
      VEC3SET(zone->grid, i, j, k);
      zone->visible = TRUE;
      }
    }
  }

/* TODO - check if no zonal division is required */
/* NB: only for isolated molecules, as periodic will */
/* always have zones to make unfragmenting easier */

#if DEBUG_ZONE_MAKE
printf("Partioning %d cores, %d shells...\n", g_slist_length(model->cores),
                                              g_slist_length(model->shels));
#endif

/* core assignment loop */
for (list=model->cores ; list ; list=g_slist_next(list))
  {
  core = list->data;

/* zone index */
  zi = zone_index(core->x, za);
  zone = za->zone[zi];
  zone->cores = g_slist_prepend(zone->cores, core);

/* DEBUG - use region colouring to show zones */
/*
core->region = zi;
*/
  }

/* shell assignment loop */
for (list=model->shels ; list ; list=g_slist_next(list))
  {
  shel = list->data;

/* zone index */
  zi = zone_index(shel->x, za);
  if (zi < 0 || zi >= za->size)
    {
    zi = 0;
    }

/* add to the list */
  zone = za->zone[zi];
  zone->shells = g_slist_prepend(zone->shells, shel);

/* DEBUG - use region colouring to show zones */
/*
core->region = zi;
*/
  }

#if DEBUG_ZONE_MAKE
for (i=0 ; i<num_zones ; i++)
  {
  zone = za->zone[i];
  j = g_slist_length(zone->cores);
  k = g_slist_length(zone->shells);
  if (j || k)
    printf("zone %d (%p) [%d %d %d] : %d cores, %d shells.\n",
           i, zone, zone->grid[0], zone->grid[1], zone->grid[2], j, k);
  }
#endif

return(za);
}

/*********************************************/
/* initialize space partitioning for a model */
/*********************************************/
void zone_init(struct model_pak *model)
{
zone_free(model->zone_array);
model->zone_array = zone_make(4.0, model);
}

/******************/
/* zone debugging */
/******************/
void zone_info(gpointer data)
{
gint i, j, k;
struct zone_pak *zone;
struct zone_array_pak *za = data;

g_assert(za != NULL);

printf("zone_info() : %p\n", za);

for (i=0 ; i<3 ; i++)
  printf("[%d] : %11.4f - %11.4f  (%f)(%d)\n",
          i, za->min[i], za->max[i],
          za->idx[i], za->div[i]);

for (i=0 ; i<za->size ; i++)
  {
  zone = za->zone[i];
  j = g_slist_length(zone->cores);
  k = g_slist_length(zone->shells);
  if (j || k)
    printf("zone %d: %d cores, %d shells, visible = %d\n", i, j, k, zone->visible);
  }
}

/**************************************/
/* get the zone index of the position */
/**************************************/
/* FIXME - if out of bounds -> need a re-zone & then return the index as well */
#define DEBUG_ZONE_INDEX 0
gint zone_index(gdouble *position, gpointer data)
{
gint i, zi, grid[3];
gdouble x[3];
struct zone_array_pak *za=data;

g_assert(za != NULL);

/* pbc constrain */
ARR3SET(x, position);

#if DEBUG_ZONE_INDEX
printf("x1: %.20f %.20f %.20f\n", x[0], x[1], x[2]);
#endif

/* NB: grid is just used as a dummy variable here */
fractional_clamp(x, grid, za->periodic);

#if DEBUG_ZONE_INDEX
printf("x2: %.20f %.20f %.20f\n", x[0], x[1], x[2]);
#endif

/* get the grid index */
for (i=0 ; i<3 ; i++)
  grid[i] = (x[i]-za->min[i]) * za->idx[i];

/* clamp to enforce bounds */
for (i=3 ; i-- ; )
  grid[i] = CLAMP(grid[i], 0, za->div[i]-1); 

#if DEBUG_ZONE_INDEX
printf(" g: %d %d %d\n", grid[0], grid[1], grid[2]);
#endif

/* get the zone index */
zi = grid[2]*za->div[1]*za->div[0]; 
zi += grid[1]*za->div[0]; 
zi += grid[0];

return(zi);
}

/*****************************************/
/* get the zone at the supplied position */
/*****************************************/
gpointer zone_get(gdouble *x, gpointer data)
{
gint zi;
struct zone_array_pak *za=data;

g_assert(za != NULL);

zi = zone_index(x, za);

return(za->zone[zi]);
}

/*********************/
/* extract zone data */
/*********************/
GSList *zone_cores(gpointer data)
{
struct zone_pak *zone = data;

g_assert(zone != NULL);

return(zone->cores);
}

/*********************/
/* extract zone data */
/*********************/
GSList *zone_shells(gpointer data)
{
struct zone_pak *zone=data;

g_assert(zone != NULL);

return(zone->shells);
}

/**********************************************/
/* construct core list from surrounding zones */
/**********************************************/
#define DEBUG_ZONE_AREA_CORES 0
GSList *zone_area_cores(gint buffer, gpointer data1, gpointer data2)
{
guint i;
gint a, b, c, zi;
gint g[3], start[3], stop[3];
GSList *list;
struct zone_pak *zone=data1, *buffer_zone;
struct zone_array_pak *za=data2;

g_assert(buffer >= 0);
g_assert(zone != NULL);
g_assert(za != NULL);

#if DEBUG_ZONE_AREA_CORES
printf("zone_area_cores(%p)\n", zone);
printf("[%2d %2d %2d]\n", zone->grid[0], zone->grid[1], zone->grid[2]);
#endif

/* setup scan limits */
ARR3SET(start, zone->grid);
ARR3SET(stop, zone->grid);

/* add buffering zones around the edge */
VEC3SUB(start, buffer, buffer, buffer);
VEC3ADD(stop, buffer, buffer, buffer);

/* eliminate non-periodic parts that exceed the boundary */
for (i=za->periodic ; i<3 ; i++)
  {
  start[i] = CLAMP(start[i], 0, za->div[i]-1);
  stop[i] = CLAMP(stop[i], 0, za->div[i]-1);
  }

/* eliminate (redundant) periodic parts that exceed width */
for (i=0 ; i<za->periodic ; i++)
  {
  if ((stop[i] - start[i]) > za->div[i]-1)
    {
    start[i] = CLAMP(start[i], 0, za->div[i]-1);
    stop[i] = CLAMP(stop[i], 0, za->div[i]-1);
    }
  }

#if DEBUG_ZONE_AREA_CORES
printf("start: %d %d %d\n", start[0], start[1], start[2]);
printf(" stop: %d %d %d\n", stop[0], stop[1], stop[2]);
#endif

g_assert(stop[0] < 1000);
g_assert(stop[1] < 1000);
g_assert(stop[2] < 1000);

g_assert(start[0] > -1000);
g_assert(start[1] > -1000);
g_assert(start[2] > -1000);

/* zone sweep */
list=NULL;
for (c=start[2] ; c<=stop[2] ; c++)
  {
  for (b=start[1] ; b<=stop[1] ; b++)
    {
    for (a=start[0] ; a<=stop[0] ; a++)
      {
/* constrain zone indexed */
      VEC3SET(g, a, b, c);
      for (i=3 ; i-- ; )
        while (g[i] < 0)
          g[i] += za->div[i];
      for (i=3 ; i-- ; )
        while (g[i] >= za->div[i])
          g[i] -= za->div[i];

/* convert to a valid zone index */
      zi = g[2]*za->div[1]*za->div[0]; 
      zi += g[1]*za->div[0]; 
      zi += g[0];

#if DEBUG_ZONE_AREA_CORES
printf(" > [%2d %2d %2d] : zone = %d ", g[0], g[1], g[2], zi);
#endif

g_assert(zi >= 0);
g_assert(zi < za->size);

/* loop over cores in the corresponding zone */
      buffer_zone = za->zone[zi];
      g_assert(buffer_zone != NULL);

/* NB: only add if core clist is non-NULL */
      if (buffer_zone->cores)
        list = g_slist_concat(list, g_slist_copy(buffer_zone->cores));

#if DEBUG_ZONE_AREA_CORES
printf("(adding cores: %d)\n", g_slist_length(buffer_zone->cores));
#endif

      }
    }
  }
return(list);
}

/***************************************/
/* convert grid coords to world coords */
/***************************************/
void zone_coords_get(gdouble *x, gint a, gint b, gint c, struct zone_array_pak *za)
{
gint i;

x[0] = (gdouble) a;
x[1] = (gdouble) b;
x[2] = (gdouble) c;

for (i=3 ; i-- ; )
  {
  x[i] /= za->idx[i];
  x[i] += za->min[i];
  }
}

/********************************/
/* checks if a point is visible */
/********************************/
gint zone_vertex_visible(gdouble *x)
{
/* FIXME */
/*
gl_get_window_coords(x, p);

if (p[0] < 0 || p[0] > sysenv.width)
  return(FALSE);
if (p[1] < 0 || p[1] > sysenv.height)
  return(FALSE);
*/

return(TRUE);
}

/*******************************/
/* zone based visibility check */
/*******************************/
void zone_visible_init(struct model_pak *model)
{
gint i, n, flag, g[3];
gdouble x[4];
GSList *list;
struct core_pak *core;
struct zone_pak *zone;
struct zone_array_pak *za;

/* CURRENT - experimental */
/* probably far too slow to be useful, maybe a better approach would */
/* be to evaluate the viewing volume & loop over all atoms & hide/unhide */

g_assert(model != NULL);
za = model->zone_array;
g_assert(za != NULL);

n = 0;
x[3] = 1.0;

/* TODO - rewrite so zone vertex visibility is evaluated only once */
for (i=za->size ; i-- ; )
  {
  zone = za->zone[i];

zone->visible = TRUE;

  if (zone->cores)
    {
    ARR3SET(g, zone->grid);
    flag = 0;

/* get cube vertices */
    zone_coords_get(x, g[0], g[1], g[2], za);
vec4mat(model->display_lattice, x);
    if (zone_vertex_visible(x))
      continue;
    zone_coords_get(x, g[0]+1, g[1], g[2], za);
vec4mat(model->display_lattice, x);
    if (zone_vertex_visible(x))
      continue;
    zone_coords_get(x, g[0]+1, g[1]+1, g[2], za);
vec4mat(model->display_lattice, x);
    if (zone_vertex_visible(x))
      continue;
    zone_coords_get(x, g[0], g[1]+1, g[2], za);
vec4mat(model->display_lattice, x);
    if (zone_vertex_visible(x))
      continue;
    zone_coords_get(x, g[0], g[1], g[2]+1, za);
vec4mat(model->display_lattice, x);
    if (zone_vertex_visible(x))
      continue;
    zone_coords_get(x, g[0]+1, g[1], g[2]+1, za);
vec4mat(model->display_lattice, x);
    if (zone_vertex_visible(x))
      continue;
    zone_coords_get(x, g[0]+1, g[1]+1, g[2]+1, za);
vec4mat(model->display_lattice, x);
    if (zone_vertex_visible(x))
      continue;
    zone_coords_get(x, g[0], g[1]+1, g[2]+1, za);
vec4mat(model->display_lattice, x);
    if (zone_vertex_visible(x))
      continue;

zone->visible = FALSE;

n++;

/*
    printf("zone: %p [%d %d %d] is not visible.\n", zone, g[0], g[1], g[2]);
*/
    }
  }

/* flag off screen cores as hidden (ie don't bother sending to the renderer) */
for (i=za->size ; i-- ; )
  {
  zone = za->zone[i];

  if (zone->visible)
    {
    for (list=zone->cores ; list ; list=g_slist_next(list))
      {
      core = list->data;
      core->status &= ~HIDDEN;
      }
    }
  else
    {
    for (list=zone->cores ; list ; list=g_slist_next(list))
      {
      core = list->data;
      core->status |= HIDDEN;
      }
    }
  }

/*
if (n)
  printf("Off-screen zones found: %d\n", n);
zone_info(za);
*/
}

/***********************************/
/* create quads for occupied zones */
/***********************************/
/* TODO - only if visible (ie on surface) */
void zone_display_init(gpointer data, struct model_pak *model)
{
gint i, g[3];
gpointer spatial;
gdouble x1[3], x2[3], x3[3], x4[3], x5[3], x6[3], x7[3], x8[3];
gdouble n[3], c1[3], c2[3], c3[3];
struct zone_pak *zone;
struct zone_array_pak *za = data;

g_assert(model != NULL);
g_assert(za != NULL);

spatial = spatial_new("zones", SPATIAL_GENERIC, 4, TRUE, model);

VEC3SET(c1, 1.0, 0.0, 0.0);
VEC3SET(c2, 0.0, 1.0, 0.0);
VEC3SET(c3, 0.0, 0.0, 1.0);

for (i=za->size ; i-- ; )
  {
  zone = za->zone[i];

  if (zone->cores)
    {
    ARR3SET(g, zone->grid);

/* get cube vertices */
    zone_coords_get(x1, g[0], g[1], g[2], za);
    zone_coords_get(x2, g[0]+1, g[1], g[2], za);
    zone_coords_get(x3, g[0]+1, g[1]+1, g[2], za);
    zone_coords_get(x4, g[0], g[1]+1, g[2], za);
    zone_coords_get(x5, g[0], g[1], g[2]+1, za);
    zone_coords_get(x6, g[0]+1, g[1], g[2]+1, za);
    zone_coords_get(x7, g[0]+1, g[1]+1, g[2]+1, za);
    zone_coords_get(x8, g[0], g[1]+1, g[2]+1, za);

/* FIXME normals (vertex orders probably) seem screwed up */
/* face 1234 */
    VEC3SET(n, 0.0, 0.0, -1.0);
    spatial_vnorm_add(x1, n, c1, spatial);
    spatial_vnorm_add(x2, n, c1, spatial);
    spatial_vnorm_add(x3, n, c1, spatial);
    spatial_vnorm_add(x4, n, c1, spatial);

/* face 5678 */
    VEC3SET(n, 0.0, 0.0, 1.0);
    spatial_vnorm_add(x5, n, c1, spatial);
    spatial_vnorm_add(x6, n, c1, spatial);
    spatial_vnorm_add(x7, n, c1, spatial);
    spatial_vnorm_add(x8, n, c1, spatial);

/* face 4378 */
    VEC3SET(n, 0.0, 1.0, 0.0);
    spatial_vnorm_add(x4, n, c2, spatial);
    spatial_vnorm_add(x3, n, c2, spatial);
    spatial_vnorm_add(x7, n, c2, spatial);
    spatial_vnorm_add(x8, n, c2, spatial);

/* face 1265 */
    VEC3SET(n, 0.0, -1.0, 0.0);
    spatial_vnorm_add(x1, n, c2, spatial);
    spatial_vnorm_add(x2, n, c2, spatial);
    spatial_vnorm_add(x6, n, c2, spatial);
    spatial_vnorm_add(x5, n, c2, spatial);


/* face 4158 */
    VEC3SET(n, -1.0, 0.0, 0.0);
    spatial_vnorm_add(x4, n, c3, spatial);
    spatial_vnorm_add(x1, n, c3, spatial);
    spatial_vnorm_add(x5, n, c3, spatial);
    spatial_vnorm_add(x8, n, c3, spatial);

/* face 2376 */
    VEC3SET(n, 1.0, 0.0, 0.0);
    spatial_vnorm_add(x2, n, c3, spatial);
    spatial_vnorm_add(x3, n, c3, spatial);
    spatial_vnorm_add(x7, n, c3, spatial);
    spatial_vnorm_add(x6, n, c3, spatial);
    }

  }

}

