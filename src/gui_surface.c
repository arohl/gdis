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
#include <time.h>
#include <stdlib.h>
#include <unistd.h>

#include "gdis.h"
#include "coords.h"
#include "edit.h"
#include "file.h"
#include "parse.h"
#include "task.h"
#include "model.h"
#include "morph.h"
#include "numeric.h"
#include "sginfo.h"
#include "matrix.h"
#include "space.h"
#include "surface.h"
#include "gui_shorts.h"
#include "interface.h"
#include "dialog.h"
#include "opengl.h"

/* main pak structures */
extern struct sysenv_pak sysenv;
extern struct elem_pak elements[];

enum
{
SURF_TITLE,
SURF_REGIONS,
SURF_ESURF_UNRE,
SURF_EATT_UNRE,
SURF_ESURF_RE,
SURF_EATT_RE,
SURF_DIPOLE,
SURF_GNORM,
SURF_MODEL,
SURF_PLANE,
SURF_SHIFT,
SURF_NCOLS
};

/* globals */
GtkTreeStore *surf_tree_store;
GtkWidget *surf_tree_view;
GtkWidget *surf_morph_energy[4];
GtkWidget *surf_hkl_family;
gint all_planes = FALSE;
gdouble rank_value = 10.0;
struct model_pak surfdata;
gchar *titles[] = {"hkl/shift", "Dhkl/depth",
                   "Esurf (u)", "Eatt (u)", "Esurf (r)", "Eatt (r)",
                   "dipole  ", "gnorm  "};

/***************************/
/* attempt to match shifts */
/***************************/
gint compare_shifts(gconstpointer *s1, gconstpointer *s2)
{
struct shift_pak *s1data;
struct shift_pak *s2data;
s1data = (struct shift_pak *) s1;
s2data = (struct shift_pak *) s2;

/* seek mismatches & return failure */
if (s1data->shift != s2data->shift)
  return(1);
if (s1data->region[0] != s2data->region[0])
  return(1);
if (s1data->region[1] != s2data->region[1])
  return(1);

/* otherwise, return succesful match */
return(0);
}

/*************************************************/
/* set values for a shift in the treeview widget */
/*************************************************/
void fn_surf_set_shift(GtkTreeIter *iter, struct shift_pak *shift)
{
gint i;
gchar *text, *old_text;
GtkTreeModel *treemodel;

treemodel = gtk_tree_view_get_model(GTK_TREE_VIEW(surf_tree_view));
for (i=SURF_TITLE ; i<=SURF_GNORM ; i++)
  {
  switch(i)
    {
    case SURF_REGIONS:
      text = g_strdup_printf("%2d:%2d", shift->region[0], shift->region[1]);
      break;

    case SURF_EATT_UNRE:
      text = g_strdup_printf("%9.4f", shift->eatt[0]);
      break;

    case SURF_ESURF_UNRE:
      text = g_strdup_printf("%9.4f", shift->esurf[0]);
      break;

    case SURF_EATT_RE:
      text = g_strdup_printf("%9.4f", shift->eatt[1]);
      break;

    case SURF_ESURF_RE:
      text = g_strdup_printf("%9.4f", shift->esurf[1]);
      break;

    case SURF_DIPOLE:
      if (shift->dipole_computed)
        text = g_strdup_printf("%9.4f", shift->dipole);
      else
        text = g_strdup(" ? ");
      break;

    case SURF_GNORM:
      if (shift->gnorm < 0.0)
        text = g_strdup(" ? ");
      else
        text = g_strdup_printf("%9.4f", shift->gnorm);
      break;

    default:
      continue;
    }

  if (treemodel)
    {
    gtk_tree_model_get(treemodel, iter, i, &old_text, -1);
    g_free(old_text);
    }
  else
    printf("ERROR: can't free old plane/shift data.\n");

  gtk_tree_store_set(surf_tree_store, iter, i, text, -1);
  g_free(text);
  }
}

/**********************************************/
/* find the iterator corresponding to a plane */
/**********************************************/
gint surf_plane_iter(GtkTreeIter *iter, struct plane_pak *plane)
{
GtkTreeModel *treemodel;
struct plane_pak *comp;

/* checks */
if (!GTK_IS_TREE_VIEW(surf_tree_view))
  return(0);

treemodel = gtk_tree_view_get_model(GTK_TREE_VIEW(surf_tree_view));
if (!treemodel)
  return(0);

/* are there branches on the tree? */
if (gtk_tree_model_get_iter_first(treemodel, iter))
  {
  do
    {
    gtk_tree_model_get(treemodel, iter, SURF_PLANE, &comp, -1);
    if (comp == plane)
      return(1);
    }
  while (gtk_tree_model_iter_next(treemodel, iter)); 
  }
return(0);
}

/**********************************************/
/* find the iterator corresponding to a shift */
/**********************************************/
gint surf_shift_iter(GtkTreeIter *iter, struct shift_pak *shift)
{
gint n;
GtkTreeModel *treemodel;
GtkTreeIter parent;
struct shift_pak *data;

/* checks */
if (!GTK_IS_TREE_VIEW(surf_tree_view))
  return(0);

treemodel = gtk_tree_view_get_model(GTK_TREE_VIEW(surf_tree_view));
if (!treemodel)
  return(0);

/* are there branches on the tree? */
if (gtk_tree_model_get_iter_first(treemodel, &parent))
  {
/* loop over all parents */
  do
    {
/* loop over all children */
    n=0;
    while (gtk_tree_model_iter_nth_child(treemodel, iter, &parent, n))
      {
      gtk_tree_model_get(treemodel, iter, SURF_SHIFT, &data, -1);
      g_assert(data != NULL);
      if (data == shift)
        return(1);
      n++;
      }
    }
  while (gtk_tree_model_iter_next(treemodel, &parent)); 
  }
return(0);
}

/******************************************************/
/* remove all shifts listed in the tree under a plane */
/******************************************************/
void surf_prune_shifts(struct plane_pak *plane)
{
gint n;
GtkTreeIter child, parent;
GtkTreeModel *treemodel;

g_assert(plane != NULL);

/* find parent level iterator */
if (!surf_plane_iter(&parent, plane))
  return;

treemodel = gtk_tree_view_get_model(GTK_TREE_VIEW(surf_tree_view));
if (!treemodel)
  return;

n = gtk_tree_model_iter_n_children(treemodel, &parent) - 1;
while (n >= 0)
  {
  if (gtk_tree_model_iter_nth_child(treemodel, &child, &parent, n))
    {
    gtk_tree_store_remove(surf_tree_store, &child);
    }
  n--;
  }
}

/***********************************************/
/* search and update a shift's treeview values */
/***********************************************/
void surf_shift_refresh(struct shift_pak *shift)
{
GtkTreeIter iter;

if (surf_shift_iter(&iter, shift))
  fn_surf_set_shift(&iter, shift);
}

/***********************************************/
/* search and update a plane's treeview values */
/***********************************************/
void surf_plane_refresh(struct plane_pak *plane)
{
}

/*********************************/
/* add a shift to the tree store */
/*********************************/
void fn_graft_shift(GtkTreeIter *parent, struct shift_pak *shift, 
                                         struct plane_pak *plane,
                                         struct model_pak *model)
{
gchar *title, *dipole;
GtkTreeIter child;

g_assert(shift != NULL);

title = g_strdup_printf("%6.4f", shift->shift);
if (shift->dipole_computed)
  dipole = g_strdup_printf("%9.4f", shift->dipole);
else
  dipole = g_strdup(" ? ");

/* set child iterator data */
gtk_tree_store_append(surf_tree_store, &child, parent);
gtk_tree_store_set(surf_tree_store, &child,
                   SURF_TITLE, title,
                   SURF_DIPOLE, dipole,
                   SURF_MODEL, model,
                   SURF_PLANE, plane,
                   SURF_SHIFT, shift, -1);

g_free(title);
g_free(dipole);

/* set shift values */
fn_surf_set_shift(&child, shift);
}

/******************************************/
/* add a list of shifts to the tree store */
/******************************************/
void fn_graft_shift_list(struct plane_pak *plane, struct model_pak *model)
{
GSList *list;
GtkTreeIter parent;
struct shift_pak *shift;

g_assert(plane != NULL);

/* acquire parent level iterator */
if (!surf_plane_iter(&parent, plane))
  gtk_tree_store_append(surf_tree_store, &parent, NULL);

/* loop over all shifts */
for (list=plane->shifts ; list ; list=g_slist_next(list))
  {
  shift = list->data;

  fn_graft_shift(&parent, shift, plane, model);
  }
}

/*********************************/
/* add a plane to the tree store */
/*********************************/
void fn_graft_plane(struct plane_pak *plane, struct model_pak *model)
{
gchar *title, *depth;
GtkTreeIter root;
GtkTreeModel *treemodel;
GtkTreePath *treepath;

g_assert(model != NULL);
g_assert(plane != NULL);

/* acquire parent level iterator */
if (!surf_plane_iter(&root, plane))
  gtk_tree_store_append(surf_tree_store, &root, NULL);

/* set the parent iterator data */
title = g_strdup_printf("(%d%d%d)",plane->index[0],plane->index[1],plane->index[2]);
depth = g_strdup_printf("%f", plane->dhkl);
gtk_tree_store_set(surf_tree_store, &root,
                   SURF_TITLE, title,
                   SURF_REGIONS, depth,
                   SURF_MODEL, model,
                   SURF_PLANE, plane,
                   SURF_SHIFT, NULL, -1);
g_free(title);
g_free(depth);


/*
printf("Adding model: %p, plane: %p\n", model, plane);
*/

fn_graft_shift_list(plane, model);

/* get treepath of root */
treemodel = gtk_tree_view_get_model(GTK_TREE_VIEW(surf_tree_view));
if (treemodel)
  {
  treepath = gtk_tree_model_get_path(treemodel, &root);
  gtk_tree_view_expand_row(GTK_TREE_VIEW(surf_tree_view), treepath, TRUE);
  gtk_tree_path_free(treepath);
  }
}

/******************************************/
/* add a list of planes to the tree store */
/******************************************/
void fn_graft_plane_list(GSList *plist, struct model_pak *model)
{
GSList *list;
struct plane_pak *plane;

/* go through the plane list */
for (list=plist ; list ; list=g_slist_next(list))
  {
  plane = list->data;
  if (plane->primary)
    fn_graft_plane(plane, model);
  }
}

/*********************/
/* add shift to list */
/*********************/
gint fn_add_plane_with_shift(void)
{
GtkTreeIter iter;
struct model_pak *model, surf;
struct plane_pak *plane;
struct shift_pak *shift;

/* checks */
model = surfdata.surface.model;
if (!model)
  return(FALSE);
if (model->id == MORPH)
  return(FALSE);

/*
P3VEC("Adding: ", surfdata.surface.miller);
*/

/* create new shift with current region values */
shift = shift_new(surfdata.surface.shift);
shift->region[0] = surfdata.surface.region[0];
shift->region[1] = surfdata.surface.region[1];

/* init for dipole calculation */
model_init(&surf);
surf.surface.shift = shift->shift;
surf.surface.region[0] = 1;
surf.surface.region[1] = 0;

/* does the plane already exist? */
plane = plane_find(surfdata.surface.miller, model);
if (plane)
  {
  ARR3SET(surf.surface.miller, plane->index);
  generate_surface(model, &surf);

  shift->dipole = surf.gulp.sdipole;
  shift->dipole_computed = TRUE;

  plane->shifts = g_slist_append(plane->shifts, shift);

/* add to shift to treestore */
  if (surf_plane_iter(&iter, plane))
    fn_graft_shift(&iter, shift, plane, model);
  else
    printf("FIXME NOW - need to create plane iter here.\n");
  }
else
  {
/* create new plane with current hkl */
  plane = plane_new(surfdata.surface.miller, model);
  if (plane)
    { 
ARR3SET(surf.surface.miller, plane->index);
generate_surface(model, &surf);

shift->dipole = surf.gulp.sdipole;
shift->dipole_computed = TRUE;

    plane->shifts = g_slist_append(plane->shifts, shift);

/* append new plane to model */
    model->planes = g_slist_append(model->planes, plane);

/* add to whole plane to treestore */
    fn_graft_plane(plane, model);
    }
  }

model_free(&surf);

return(TRUE);
}

/**********************************************/
/* search for valid shifts - output to a file */
/**********************************************/
gint fn_add_valid_shifts(void)
{
struct model_pak *model;
struct surface_pak *surfdat = &surfdata.surface;
struct plane_pak *plane;

/* checks */
model = surfdat->model;
g_return_val_if_fail(model != NULL, 1);
if (model->periodic != 3)
  {
  gui_text_show(ERROR, "Base model is not 3D periodic.\n");
  return(1);
  }
if (!surfdata.surface.miller[0] &&
    !surfdata.surface.miller[1] &&
    !surfdata.surface.miller[2])
  {
  gui_text_show(ERROR, "Don't be silly.\n");
  return(2);
  }

/* does the plane already exist? */
plane = plane_find(surfdata.surface.miller, model);
if (plane)
  surf_prune_shifts(plane);
else
  {
  plane = plane_new(surfdata.surface.miller, model);
/* add to list */
  if (plane)
    model->planes = g_slist_append(model->planes, plane);
  else
    return(3);
  }

/* find valid shifts for the plane */
/* NB: this REPLACES existing shifts */
calc_valid_shifts(model, plane);

/* update treestore */
fn_graft_plane(plane, model);

return(FALSE);
}

/*************/
/* rank cuts */
/*************/
gint zsort(gdouble *z1, gdouble *z2)
{
if (*z1 > *z2)
  return(1);

if (*z2 > *z1)
  return(-1);

return(0);
}

/***************************************/
/* smallest common multiple for a pair */
/***************************************/
gint scf(gint a, gint b)
{
if (!b)
  return(abs(a));
return(scf(b, a%b));
}

/***************************************/
/* smallest common multiple for vector */
/***************************************/
gint simplify(gint *vec, gint size)
{
gint i, n;

n = scf(*vec, *(vec+1));  
for (i=2 ; i<size ; i++)
  n = scf(n, *(vec+i));  
return(n);
}

/***************************************************/
/* construct a morphology within the current model */
/***************************************************/
void surf_make_morph(GtkWidget *w, struct model_pak *model)
{
g_assert(model != NULL);

/* generate symmetry related faces */
/* FIXME - handle multiple calls (ie duplication) */
surf_symmetry_generate(model);

/* generate vertices */
morph_build(model);

/* init for display */
coords_init(CENT_COORDS, model);
redraw_canvas(SINGLE);
}

/*************************************/
/* construct a particular morphology */
/*************************************/
void make_morph(void)
{
gint status;
gchar *filename;
/*GSList *plist;*/
struct model_pak *data;
/*struct plane_pak *pdata;*/

/* get the current model */
data = surfdata.surface.model;
if (!data)
  return;

/* create a suitable name */
filename = g_strdup_printf("%s.gmf", data->basename);

/* store the planes */
write_gmf(filename, data);

/* read morphology file in & compute hull */
data = model_new();
/* default morphology computation type */
data->morph_type = DHKL;

/* only successful if properly loaded */
status = 2;
if (data)
  {
  status--;
  if (!read_gmf(filename, data))
    status--;
  }

if (!status)
  {
  tree_model_add(data);
  }

g_free(filename);
}

/* TODO - rewrite to generate specified number of faces (rank_spin) */
/* based on Dhkl ranking  - maybe! - calc valid shifts for each */
gint dhkl_compare(gpointer ptr1, gpointer ptr2)
{
struct plane_pak *plane1=ptr1, *plane2=ptr2;

if (plane1->dhkl > plane2->dhkl)
  return(-1);
if (plane1->dhkl < plane2->dhkl)
  return(1);

return(0);
}

/******************************************************/
/* adds ranked faces (from data) to the supplied list */
/******************************************************/
/* if num -> get num ranked faces, else get until Dhkl < min */
#define DEBUG_GET_RANKED_FACES 0
GSList *get_ranked_faces(gint num, gdouble min, struct model_pak *data)
{
gint n, h, k, l, limit;
gint f1[3], f2[3];
gdouble m[3];
GSList *list=NULL, *plist=NULL, *list1=NULL, *list2=NULL;
struct plane_pak *pdata;

/* loop over a large range */
/* FIXME - should loop until num/max is satisfied */
limit = 4;
for (h=limit ; h>=-limit ; h--)
  {
  for (k=limit ; k>=-limit ; k--)
    {
    for (l=limit ; l>=-limit ; l--)
      {
/* skip (000) */
      if (!h && !k && !l)
        continue;
      VEC3SET(f1, h, k, l);
/* skip things with a common multiple > 1, since */
/* fn_make_plane() takes care of systematic absences */
/* NB: this is ok for morph pred. but for XRD we may want such planes */
if (num)
  {
      n = simplify(f1, 3);
      if (n != 1)
        continue;
  }

/* add the plane */
      ARR3SET(m, f1);
      pdata = plane_new(m, data);
      if (pdata)
        plist = g_slist_prepend(plist, pdata);
      }
    }
  }
plist = g_slist_reverse(plist);

/* sort via Dhkl */
plist = g_slist_sort(plist, (gpointer) dhkl_compare);

/*
VEC3SET(f1, 1, 0, -4);
printf("f1 : [%d %d %d]\n", f1[0], f1[1], f1[2]);
n = simplify(f1, 3);
printf("n = %d\n", n);
*/

/* mark the planes that are symmetry related */
list1 = plist;
while (list1 != NULL)
  {
/* get a reference (primary) plane */
  pdata = list1->data;
  if (pdata->primary)
    {
    ARR3SET(f1, pdata->index);
    }
  else
    {
    list1 = g_slist_next(list1);
    continue;
    }
/* scan for symmetry related planes */
  list2 = g_slist_next(list1);
  while (list2 != NULL)
    {
    pdata = list2->data;
    ARR3SET(f2, pdata->index);
/* get the hkl's and test */
#if DEBUG_GET_RANKED_FACES
printf("testing: %d %d %d  &  %d %d %d\n", f1[0], f1[1], f1[2], f2[0], f2[1], f2[2]);
#endif
    if (facet_equiv(data,f1,f2))
      pdata->primary = FALSE;

    list2 = g_slist_next(list2);
    }
  list1 = g_slist_next(list1);
  }

/* write the planes list to the model */
/* write only up to the requested number */
list1 = plist;
n=0;
while (list1 != NULL)
  {
/* get a reference (primary) plane */
  pdata = list1->data;
/*
  if (pdata->primary)
*/
    {
/* automatically add valid shifts */
    calc_valid_shifts(data, pdata);

/* default single shift */
/*
    sdata = create_shift(0.0);
    pdata->shifts = g_slist_append(pdata->shifts, sdata);
*/

/* add plane */
/* TODO - don't prepend? (may be others already there) */
/*
    list = g_slist_prepend(list, pdata);
*/

/* number of unique planes */ 
if (pdata->primary)
  n++;

    if (num)
      {
/* got enough unique planes? */
/*
      if (n == num)
*/

      if (n > num)
        break;
      else
        list = g_slist_prepend(list, pdata);

      }
    else
      {
/*
      if (pdata->dhkl < min)
        break;
*/
      if (pdata->dhkl > min)
        list = g_slist_prepend(list, pdata);

      }
    }
  list1 = g_slist_next(list1);
  }
list = g_slist_reverse(list);
return(list);
}

/**********************/
/* GULP debugging aid */
/**********************/
void gulp_dump(struct model_pak *data)
{
printf("  num atoms = %d\n", data->num_atoms);
printf("        sbe = %f\n", data->gulp.sbulkenergy);
printf("     energy = %f\n", data->gulp.energy);
printf("E surf (un) = %f\n", data->gulp.esurf[0]);
printf("E surf (re) = %f\n", data->gulp.esurf[1]);
printf(" E att (un) = %f\n", data->gulp.eatt[0]);
printf(" E att (re) = %f\n", data->gulp.eatt[1]);
}

/**************************************************/
/* eliminate repeated z values in a (sorted) list */
/**************************************************/
GSList *trim_zlist(GSList *zlist)
{
gdouble dz, *z1, *z2;
GSList *list, *rlist;

/* get 1st item */
list = zlist;
if (list)
  {
  z1 = (gdouble *) list->data;
  list = g_slist_next(list);
  }
else
  return(NULL);

/* compare and eliminate repeats */
rlist = zlist;
while (list)
  {
  z2 = (gdouble *) list->data;
  list = g_slist_next(list);

  dz = fabs(*z1 - *z2);
  if (dz < FRACTION_TOLERANCE)
    {
    rlist = g_slist_remove(rlist, z2);
    g_free(z2);
    }
  else
    z1 = z2;
  }
return(rlist);
}

/*********************************************/
/* append a list of valid shifts for a plane */
/*********************************************/
#define DEBUG_CALC_SHIFTS 0
gint calc_valid_shifts(struct model_pak *data, struct plane_pak *pdata)
{
gint build, num_cuts, num_shifts, dummy[3];
gdouble gcd, zlim, *z1, *z2, vec[3];
GSList *slist, *list, *zlist;
struct model_pak *surf;
struct shift_pak *sdata;
struct core_pak *core;
struct mol_pak *mol;

/* allocate & template */
surf = g_malloc(sizeof(struct model_pak));
model_init(surf);

surf->mode = FREE;
surf->id = GULP;
surf->periodic = 2;

/* actual surface we want constructed */
ARR3SET(surf->surface.miller, pdata->index);
surf->surface.shift = 0.0;
surf->surface.region[0] = 1.0;
surf->surface.region[1] = 0.0;

gcd = GCD(pdata->index[0], GCD(pdata->index[1], pdata->index[2]));

#if DEBUG_CALC_SHIFTS
printf("*** calc shifts ***\n");
printf("miller: %f %f %f\n",surf->surface.miller[0]
                            ,surf->surface.miller[1]
                            ,surf->surface.miller[2]);
printf("region: %d %d\n", (gint) surf->surface.region[0],
                          (gint) surf->surface.region[1]);
printf("   gcd: %f\n", gcd);
printf(" Shift: %f\n", surf->surface.shift);
#endif

generate_surface(data, surf);
if (surf->surface.dspacing == 0.0)
  {
  printf("calc_valid_shifts() error, invalid dspacing.\n");
  goto calc_valid_shifts_cleanup;
  }

#if DEBUG_CALC_SHIFTS
printf("  Dhkl: %f\n", surf->surface.dspacing);
#endif

/* transfer appropriate setup data to new model */
gulp_data_copy(data, surf);
surf->gulp.run = E_SINGLE;
surf->gulp.method = CONV;

/* NEW */
pdata->area = surf->area;

/* only sample one dhkl's worth of cuts */
zlim = G_MINDOUBLE - 1.0 / gcd;

/* get all possible shifts from z coordinates */
zlist = NULL;
/*
z1 = g_malloc(sizeof(gdouble));
*z1 = 0.0;
zlist = g_slist_prepend(zlist, z1);
*/
if (data->surface.ignore_bonding)
  {
#if DEBUG_CALC_SHIFTS
printf("Ignoring bonding...\n");
#endif
/* core z coord */
  for (list=surf->cores ; list ; list=g_slist_next(list))
    {
    core = list->data;
    ARR3SET(vec, core->x);

/* helps get around some precision problems */
    vec[2] = decimal_round(vec[2], 6);

/* grab one slice [0, -1.0) -> [0, -Dhkl) */
    if (vec[2] < zlim)
      continue;
    if (vec[2] > 0.0)
      continue;

/* NB: cartesian z direction is opposite to shift value sign */
    z1 = g_malloc(sizeof(gdouble));
    *z1 = -vec[2];
    fractional_clamp(z1, dummy, 1);
    zlist = g_slist_prepend(zlist, z1);
    }
  }
else
  {
#if DEBUG_CALC_SHIFTS
printf("Using molecule centroids... [%d]\n", g_slist_length(surf->moles));
#endif
/* molecule centroid */
  for (list=surf->moles ; list ; list=g_slist_next(list))
    {
    mol = list->data;
    ARR3SET(vec, mol->centroid);

/* helps get around some precision problems */
    vec[2] = decimal_round(vec[2], 6);

/* grab one slice [0, -1.0) -> [0, -Dhkl) */
    if (vec[2] < zlim)
      continue;
    if (vec[2] > 0.0)
      continue;

/* NB: cartesian z direction is opposite to shift value sign */
    z1 = g_malloc(sizeof(gdouble));
    *z1 = -vec[2];
    fractional_clamp(z1, dummy, 1);
    zlist = g_slist_prepend(zlist, z1);
    }
  }
zlist = g_slist_sort(zlist, (gpointer) zsort);
num_cuts = g_slist_length(zlist);

/* NEW - some pathalogical cases have shift values just above 0.0 */
/* resulting in no cuts - so enforce at least one */
if (!num_cuts)
  {
  z1 = g_malloc(sizeof(gdouble));
  *z1 = 0.0;
  zlist = g_slist_prepend(zlist, z1);
  num_cuts = 1;
  }

#if DEBUG_CALC_SHIFTS
printf("Found %d shifts.\n", num_cuts);
for (list=zlist ; list ; list=g_slist_next(list))
  {
  z1 = (gdouble *) list->data;
  printf(" - %.20f\n", *z1);
  }
printf("\n");
#endif

zlist = trim_zlist(zlist);
num_cuts = g_slist_length(zlist);

#if DEBUG_CALC_SHIFTS
printf("Found %d unique shifts.\n", num_cuts);
for (list=zlist ; list ; list=g_slist_next(list))
  {
  z1 = (gdouble *) list->data;
  printf("%.4f ", *z1);
  }
printf("\n");
#endif

/* failsafe */
if (!zlist)
  goto calc_valid_shifts_cleanup;

/* use midpoints to avoid cutoff problems (also neater shift values) */
/* also, the unchanged first cut will be forced to 0.0 */
zlist = g_slist_reverse(zlist);
z1 = zlist->data;
list = g_slist_next(zlist);
while (list)
  {
  z2 = list->data;
  *z1 = 0.5*(*z1 + *z2);
  z1 = z2;
  list = g_slist_next(list);
  }
*z1 = 0.0;

#if DEBUG_CALC_SHIFTS
printf("Midpoint equivalent cuts:\n");
for (list=zlist ; list ; list=g_slist_next(list))
  {
  z1 = (gdouble *) list->data;
  printf(" %.20f", *z1);
  }
printf("\n");
#endif

/* re-rank */
zlist = g_slist_sort(zlist, (gpointer) zsort);

#if DEBUG_CALC_SHIFTS
printf("Ranked shifts: %d\n", num_cuts);
for (list=zlist ; list ; list=g_slist_next(list))
  {
  z1 = (gdouble *) list->data;
  printf(" %f", *z1);
  }
printf("\n");
#endif

zlist = trim_zlist(zlist);
num_shifts = g_slist_length(zlist);

#if DEBUG_CALC_SHIFTS
printf("Unique shifts: %d\n", num_shifts);
for (list=zlist ; list ; list=g_slist_next(list))
  {
  z1 = (gdouble *) list->data;
  printf(" %f", *z1);
  }
printf("\n");
#endif

/* turn bonding off */
build = data->build_molecules;
if (data->surface.ignore_bonding && build)
  {
  data->build_molecules = FALSE;
  connect_bonds(data);
  connect_molecules(data);
  }

/* generate & check each one */
/* NB: this'll overwrite our 1st model, saving having to delete it */
slist = NULL;
for (list=zlist ; list ; list=g_slist_next(list))
  {
  z1 = (gdouble *) list->data;

/* destroy surfs core/shell lists */
  free_core_list(surf);

/* make the surface (also computes dipole) */
  surf->surface.shift = *z1;
  generate_surface(data, surf);

#if DEBUG_CALC_SHIFTS
printf("Shift = %f, dipole = %f, bonds broken = %d / %f\n", *z1, surf->gulp.sdipole, surf->surface.bonds_cut, surf->area);
#endif

/* create new model if valid cut (unless user wants all cuts) */
  if (fabs(surf->gulp.sdipole) < data->surface.dipole_tolerance
      || data->surface.include_polar)
    {
    sdata = shift_new(*z1);
    sdata->dipole = surf->gulp.sdipole;
    sdata->dipole_computed = TRUE;

/* NEW - broken bonds per area calc */
    sdata->bbpa = (gdouble) surf->surface.bonds_cut;
    sdata->bbpa /= surf->area;

    slist = g_slist_prepend(slist, sdata);
    }
  }

/* new if shift list is non-empty, overwrite only if we found something */
if (slist)
  {
  slist = g_slist_reverse(slist);
  if (pdata->shifts)
    {
    shift_data_free(pdata->shifts);
    g_slist_free(pdata->shifts);
    }
  pdata->shifts = slist;
  }

/* restore source model's connectivity */
if (data->surface.ignore_bonding && build)
  {
  data->build_molecules = build;
  connect_bonds(data);
  connect_molecules(data);
  }

/* this has data only if the above was successful */
free_slist(zlist);

calc_valid_shifts_cleanup:

/* cleanup */
model_free(surf);
g_free(surf);

return(0);
}

/***********************************/
/* computes the dipole for a plane */
/***********************************/
#define DEBUG_TEST_VALID_SHIFTS 0
void test_valid_shifts(struct plane_pak *plane, struct model_pak *data)
{
GSList *list;
struct shift_pak *shift;
struct model_pak surf;

#if DEBUG_TEST_VALID_SHIFTS
printf("Testing for valid shifts [%d %d %d]\n",
        plane->index[0], plane->index[1], plane->index[2]); 
#endif

/* init the model to be used for generating shifts */
model_init(&surf);
gulp_data_copy(data, &surf);

/* test the shifts */
for (list=plane->shifts ; list ; list=g_slist_next(list))
  {
  shift = list->data;

/* init the surface we want */
  ARR3SET(surf.surface.miller, plane->index);
  surf.surface.region[0] = 1;
  surf.surface.region[1] = 0;
  surf.surface.shift = shift->shift;

/* destroy old cores & then create the surface */
  free_core_list(&surf);
  generate_surface(data, &surf);

#if DEBUG_TEST_VALID_SHIFTS
printf("[%f:%f]\n", shift->shift, surf.gulp.sdipole);
#endif

  shift->dipole = surf.gulp.sdipole;
  shift->dipole_computed = TRUE;
  }
}

/******************************/
/* do the energy calculations */
/******************************/
#define DEBUG_EXEC_ECALC_TASK 0
void exec_ecalc_task(struct model_pak *model, gpointer data)
{
gchar *inp;
struct task_pak *task = data;

g_assert(model != NULL);
g_assert(task != NULL);

/* build the full path filename */
inp = g_build_filename(sysenv.cwd, model->gulp.temp_file, NULL);
if (write_gulp(inp, model))
  printf("write failed.\n");
else
  {
  if (!model->gulp.no_exec)
    {
/* NEW */
    task->status_file = g_build_filename(sysenv.cwd, model->gulp.out_file, NULL);

    if (exec_gulp(model->gulp.temp_file, model->gulp.out_file))
      printf("exec failed.\n");
    }
  }
g_free(inp);
}

/*********************************************/
/* process result from an energy calculation */
/*********************************************/
#define DEBUG_PROC_ECALC_TASK 0
void proc_ecalc_task(struct model_pak *model)
{
gchar *out;
struct shift_pak *sdata;

g_assert(model != NULL);

/* don't process if gulp execution was turned off */
if (model->gulp.no_exec)
  return;

/* flag that gulp read routine shouldn't prep the model */
model->grafted = TRUE;

/* read gulp file */
out = g_build_filename(sysenv.cwd, model->gulp.out_file, NULL);
if (read_gulp_output(out, model))
  printf("proc_ecalc_task() error: read failed.\n");
g_free(out);

#if DEBUG_PROC_ECALC_TASK
gulp_dump(model);
#else
/*
unlink(model->gulp.temp_file);
unlink(model->gulp.out_file);
printf("inp: %s\n", model->gulp.temp_file);
printf("out: %s\n", model->gulp.out_file);
*/
#endif

/* feed energy into shift */
sdata = (model->planes)->data;
if (sdata)
  {
  sdata->esurf[0] = model->gulp.esurf[0];
  sdata->esurf[1] = model->gulp.esurf[1];
  sdata->eatt[0] = model->gulp.eatt[0];
  sdata->eatt[1] = model->gulp.eatt[1];
  sdata->gnorm = model->gulp.gnorm;
  }

/* gui update */
surf_shift_refresh(sdata);

/* free model */
/* NB: this was a borrowed pointer, so we don't want it freed */
model->planes = NULL;
model_free(model);
g_free(model);
}

/***********************************************/
/* submit a background energy calculation task */
/***********************************************/
void new_ecalc_task(struct model_pak *model,
                    struct plane_pak *pdata,
                    struct shift_pak *sdata)
{
gint r1size;
GSList *list;
struct model_pak *surf;
struct core_pak *core;

g_assert(model != NULL);
g_assert(pdata != NULL);
g_assert(sdata != NULL);

/* allocate & init new model for the surface */
surf = g_malloc(sizeof(struct model_pak));
model_init(surf);
gulp_data_copy(model, surf);
surf->id = GULP;
surf->gulp.method = CONV;
/* NEW - no dependance on this */
surf->surface.model = NULL;
surf->surface.shift = sdata->shift;
surf->surface.region[0] = sdata->region[0];
surf->surface.region[1] = sdata->region[1];
ARR3SET(surf->surface.miller, pdata->index);

/* add the cores */
generate_surface(model, surf);
gulp_files_init(surf);

/* compute surface bulk energy */
surf->gulp.sbulkenergy = model->gulp.energy;
surf->gulp.sbulkenergy /= (gdouble) model->num_atoms;
r1size=0;
for (list=surf->cores ; list ; list=g_slist_next(list))
  {
  core = list->data;
  if (core->region == REGION1A)
    r1size++;
  }
if (r1size)
  surf->gulp.sbulkenergy *= (gdouble) r1size;
else
  printf("Warning: empty region 1.\n");

/* FIXME - cheat, by tacking the shift to update on the *planes* slist */
g_assert(surf->planes == NULL);
surf->planes = g_slist_append(surf->planes, sdata);

task_new("Energy", &exec_ecalc_task, surf, &proc_ecalc_task, surf, NULL);
}

/***********************************/
/* single region convergence cycle */
/***********************************/
#define DEBUG_SURF_CONV 1
gint surf_conv(FILE *fp, struct model_pak *model, gint type)
{
gint i, n, flag, relax, region, r1size, status;
gdouble sbe, de, old_energy;
gchar *inp, *out, *full_inp, *full_out;
GSList *list;
struct model_pak *surf;
struct core_pak *core;
struct plane_pak *plane;
struct shift_pak *shift;

g_assert(model != NULL);

/* retrieve plane & shift */
plane = (model->planes)->data;
shift = (plane->shifts)->data;

switch(type)
  {
  case REGION1A:
    region = 0;
    break;
  case REGION2A:
    region = 1;
    break;
  default:
    printf("surf_conv() error: bad region type.\n");
    return(1);
  }

/* init surface */
surf = g_malloc(sizeof(struct model_pak));
model_init(surf);
gulp_data_copy(model, surf);
ARR3SET(surf->surface.miller, model->surface.miller);
surf->surface.region[0] = model->surface.region[0];
surf->surface.region[1] = model->surface.region[1];
surf->surface.shift = model->surface.shift;
surf->surface.converge_eatt = model->surface.converge_eatt;
surf->sginfo.lookup = FALSE;

/* which energy to converge */
if (model->gulp.run == E_OPTIMIZE)
  relax = 1;
else
  relax = 0;

/* surface bulk energy per atom */
sbe = model->gulp.energy;
sbe /= (gdouble) model->num_atoms;

/* converge region size */
n=flag=0;
old_energy = 0.0;
de = 99999999.9;
status = 0;
for (;;)
  {
/* FIXME - any other stuff to free? */
  free_core_list(surf);
  generate_surface(model, surf);

/* compute surface bulk energy */
  r1size=0;
  for (list=surf->cores ; list ; list=g_slist_next(list))
    {
    core = list->data;
    if (core->region == REGION1A)
      r1size++;
    }
  if (!r1size)
    printf("Warning: empty region 1.\n");

  surf->gulp.sbulkenergy = sbe * (gdouble) r1size;

/* create appropriate name */
  inp = g_strdup_printf("%s_%d_%d.gin", surf->basename,
                                        (gint) surf->surface.region[0],
                                        (gint) surf->surface.region[1]);

  out = g_strdup_printf("%s_%d_%d.got", surf->basename,
                                        (gint) surf->surface.region[0],
                                        (gint) surf->surface.region[1]);

  full_inp = g_build_filename(sysenv.cwd, inp, NULL);
  full_out = g_build_filename(sysenv.cwd, out, NULL);

/* command cycle: write input, execute, read output */
/*
if (0)
*/
for (i=0 ; i<3 ; i++)
  {
  if (!status)
    {
    switch(i)
      {
      case 0:
        status = write_gulp(full_inp, surf);
        break;

      case 1:
        status = exec_gulp(inp, out);
        if (!status)
          unlink(full_inp);
/*
printf("input: %s\n", full_inp);
*/
        break;

      case 2:
/* read_gulp_out() can call gui_text_show() which modifies */
/* the message widget, hence the thread lock is required */
/* FIXME - currently, this function surf_conv() is always */
/* called in a thread but what if it isn't? ie will the */
/* threads enter call deadlock if we're in the main GTK loop? */
        gdk_threads_enter();
        status =  read_gulp_output(full_out, surf);
        if (!status)
          unlink(full_out);
/*
printf("output: %s\n", full_out);
*/
        gdk_threads_leave();
        break;
      }
    }
  }

/* remove old files for next cycle */
  g_free(full_inp);
  g_free(full_out);
  g_free(inp);
  g_free(out);

/* get difference */
  if (surf->surface.converge_eatt)
    {
    de = surf->gulp.eatt[relax] - old_energy;
    old_energy = surf->gulp.eatt[relax];
    }
  else
    {
    de = surf->gulp.esurf[relax] - old_energy;
    old_energy = surf->gulp.esurf[relax];
    }

/* bad convergence checks */
  if (n > 10 && de > 0.1)
    {
    fprintf(fp, "ERROR: failed convergence cycle.\n");
    status = 1;
    break;
    }

if (surf->surface.converge_eatt)
  {
  fprintf(fp, "[%d:%d]  Eatt = %f (%f)\n", (gint) surf->surface.region[0],
                                      (gint) surf->surface.region[1],
                                      surf->gulp.eatt[relax], de);
  }
else
  {
  fprintf(fp, "[%d:%d] Esurf = %f (%f)\n", (gint) surf->surface.region[0],
                                      (gint) surf->surface.region[1],
                                      surf->gulp.esurf[relax], de);
  }
fflush(fp);

/* convergence check */
  if (surf->surface.converge_eatt)
    {
    if (fabs(de) < MAX_DEEATT)
      break;
    }
  else
    {
    if (fabs(de) < MAX_DESURF)
      break;
    }

/* next size */
  surf->surface.region[region]++;
  n++;
  }
if (!status)
  fprintf(fp, "Region %d converged.\n", region+1);
fflush(fp);

/* we can go back one size, due to the way convergence is checked */
surf->surface.region[region]--;

/* transfer data to source model */
model->surface.region[region] = surf->surface.region[region];
model->gulp.eatt[region] = surf->gulp.eatt[region];
model->gulp.esurf[region] = surf->gulp.esurf[region];
model->gulp.gnorm = surf->gulp.gnorm;

/* cleanup */
model_free(surf);
g_free(surf);

return(status);
}

/*****************************************/
/* region convergence calculation (task) */
/*****************************************/
#define DEBUG_EXEC_REGCON_TASK 0
void exec_regcon_task(struct model_pak *model, struct task_pak *task)
{
gint m, r, run;
gchar *name, *tmp;
struct plane_pak *plane;
FILE *fp;

/* checks */
g_assert(model != NULL);
g_assert(model->planes != NULL);
g_assert(model->periodic == 3);

/* NEW - status file */
tmp = gun("txt");
name = g_build_filename(sysenv.cwd, tmp, NULL);
g_free(tmp);
fp = fopen(name, "wt");
if (fp)
  task->status_file = name;
else
  {
  fp = stdout;
  g_free(name);
  }

/* retrieve plane & shift */
plane = (model->planes)->data;

/* estimate starting region sizes from the dspacing for the plane */
g_return_if_fail(plane != NULL);

/* mult dhkl by the gcd */
m = GCD(plane->index[0], GCD(plane->index[1], plane->index[2]));

if (plane->dhkl < MIN_THICKNESS)
  {
/* FIXME - what to do if dhkl is very small or even zero */
  r = 1 + (gint) (MIN_THICKNESS / (m*plane->dhkl));

  if (r > model->surface.region[0])
    model->surface.region[0] = r;
  if (r > model->surface.region[1])
    model->surface.region[1] = r;
  }

fprintf(fp, "----------------------------------------------\n");
fprintf(fp, "         Miller: %d %d %d \n", plane->index[0] , plane->index[1] , plane->index[2]);
fprintf(fp, "           Dhkl: %f\n", plane->dhkl);
fprintf(fp, "          Shift: %f\n", model->surface.shift);
fprintf(fp, "Initial regions: %d , %d\n", (gint) model->surface.region[0],
                                          (gint) model->surface.region[1]);
fprintf(fp, "    Convergence: ");
if (model->surface.converge_eatt)
  fprintf(fp, "Attachment Energy based\n");
else
  fprintf(fp, "Surface Energy based\n");
fprintf(fp, "----------------------------------------------\n");

/* save run type */
run = model->gulp.run;

fprintf(fp, "Beginning unrelaxed convergence...\n");

/* converge region 2 size */
model->gulp.run = E_SINGLE;
if (model->surface.converge_r2)
  {
  if (surf_conv(fp, model, REGION2A))
    {
    task->message = g_strdup("Failed unrelaxed convergence of region 2.\n");
    return;
    }
  }
else
  fprintf(fp, "Skipping region 2 convergence...\n");
  
/* converge region 1 size */
if (model->surface.converge_r1)
  {
  if (surf_conv(fp, model, REGION1A))
    {
    task->message = g_strdup("Failed unrelaxed convergence of region 1.\n");
    return;
    }
  }
else
  fprintf(fp, "Skipping region 1 convergence...\n");

/* if opti is specified - repeat with unrelaxed regions as starting point */
if (run == E_OPTIMIZE)
  {
  model->gulp.run = E_OPTIMIZE;

  fprintf(fp, "Beginning relaxed convergence...\n");
  fflush(fp);

/* converge region 2 size */
  if (model->surface.converge_r2)
    {
    if (surf_conv(fp, model, REGION2A))
      {
      if (task)
        task->message = g_strdup("Failed relaxed convergence of region 2.\n");
      return;
      }
    fflush(fp);
    }
  else
    fprintf(fp, "Skipping region 2 convergence...\n");

/* converge region 1 size */
  if (model->surface.converge_r1)
    {
    if (surf_conv(fp, model, REGION1A))
      {
      if (task)
        task->message = g_strdup("Failed relaxed convergence of region 1.\n");
      return;
      }
    fflush(fp);
    }
  else
    printf("Skipping region 1 convergence...\n");
  }
fclose(fp);
}

/**************************************/
/* region convergence task processing */
/**************************************/
#define DEBUG_PROC_REGON_TASK 0
void proc_regcon_task(struct model_pak *model)
{
struct plane_pak *plane;
struct shift_pak *shift;

/* checks */
g_assert(model  != NULL);
plane = (model->planes)->data;
g_assert(plane != NULL);
shift = (plane->shifts)->data;
g_assert(shift != NULL);

/* the shift pointer belongs to the main model, so */
/* we don't want it free'd when model is destroyed */
plane->shifts = NULL;

#if DEBUG_PROC_REGON_TASK
printf("Final regions: %d , %d\n", (gint) model->surface.region[0],
                                   (gint) model->surface.region[1]);
printf("Final Esurf: %f , %f\n", model->gulp.esurf[0], model->gulp.esurf[1]);
printf("Final Eatt: %f , %f\n", model->gulp.eatt[0], model->gulp.eatt[1]);
printf("Final gnorm: %f\n", model->gulp.gnorm);
#endif

shift->region[0] = model->surface.region[0];
shift->region[1] = model->surface.region[1];
shift->esurf[0] = model->gulp.esurf[0];
shift->esurf[1] = model->gulp.esurf[1];
shift->eatt[0] = model->gulp.eatt[0];
shift->eatt[1] = model->gulp.eatt[1];
shift->gnorm = model->gulp.gnorm;

/* update gui */
surf_shift_refresh(shift);
shift->locked = FALSE;

/* cleanup */
model_free(model);
g_free(model);
}

/*********************************/
/* region convergence task setup */
/*********************************/
void new_regcon_task(struct model_pak *model,
                     struct plane_pak *plane,
                     struct shift_pak *shift)
{
GSList *list;
struct core_pak *core;
struct plane_pak *plane2;
struct model_pak *temp;

/* checks */
g_assert(model != NULL);
g_assert(plane != NULL);
g_assert(shift != NULL);

/* NEW - prevent deletion */
shift->locked = TRUE;

/* duplicate the data passed, since it may be changed */
/* or destroyed by the user while the task is still queued */
temp = g_malloc(sizeof(struct model_pak));
model_init(temp);

/* duplicate source model data */
temp->gulp.energy = model->gulp.energy;
gulp_data_copy(model, temp);
temp->cores = dup_core_list(model->cores);
/* FIXME - this is ugly, but the only way (curently) to do it, */
/* as dup_core_list dups shells - but doesn't add them to a list */
for (list=temp->cores ; list ; list=g_slist_next(list))
  {
  core = list->data;
  if (core->shell)
    temp->shels = g_slist_prepend(temp->shels, core->shell);
  }

/* TODO - implement a copy_lattice_info() primitive */
temp->periodic = model->periodic;
temp->fractional = model->fractional;
memcpy(temp->pbc, model->pbc, 6 * sizeof(gdouble));
memcpy(temp->latmat, model->latmat, 9 * sizeof(gdouble));
memcpy(temp->ilatmat, model->ilatmat, 9 * sizeof(gdouble));
temp->sginfo.spacename = g_strdup(model->sginfo.spacename);
temp->sginfo.cellchoice = model->sginfo.cellchoice;

/* always use name for lookup */
temp->sginfo.spacenum = -1;

if (model->surface.ignore_bonding)
  temp->build_molecules = FALSE;

/* TODO - enforce unfragment() ??? */
model_prep(temp);

/*
temp->fractional = TRUE;
zone_init(temp);
connect_bonds(temp);
connect_molecules(temp);
*/

/* NEW - no dependence once we've spawned the task */
temp->surface.converge_eatt = model->surface.converge_eatt;
temp->surface.model = NULL;
temp->surface.shift = shift->shift;
temp->surface.region[0] = shift->region[0];
temp->surface.region[1] = shift->region[1];
ARR3SET(temp->surface.miller, plane->index);

/* append the required plane and shift to the temporary model */
plane2 = plane_new(temp->surface.miller, model);
plane2->shifts = g_slist_append(plane2->shifts, shift);
g_assert(temp->planes == NULL);
temp->planes = g_slist_append(temp->planes, plane2);

task_new("Regcon", &exec_regcon_task, temp, &proc_regcon_task, temp, NULL);
}

/*************************/
/* surface creation task */
/*************************/
void cb_surf_create(GtkWidget *w, struct model_pak *model)
{
struct plane_pak *plane;
struct shift_pak *shift;

/* setup plane */
plane = plane_new(surfdata.surface.miller, model);
g_return_if_fail(plane != NULL);

/* setup shift */
shift = shift_new(surfdata.surface.shift);
shift->region[0] = surfdata.surface.region[0];
shift->region[1] = surfdata.surface.region[1];

/* create */
make_surface(model, plane, shift);  

/* avoid BSOD */
coords_init(INIT_COORDS, model);

redraw_canvas(SINGLE);

g_free(plane);
g_free(shift);
}

/*****************************************/
/* act on all selected tasks in the list */
/*****************************************/
void surf_task_selected(GtkWidget *w, gint type)
{
gint refresh=FALSE;
GList *list, *row;
GtkTreeIter iter;
GtkTreeModel *treemodel;
GtkTreeSelection *selection;
struct model_pak *model;
struct plane_pak *plane;
struct shift_pak *shift;

/* checks */
treemodel = GTK_TREE_MODEL(surf_tree_store);
selection = gtk_tree_view_get_selection(GTK_TREE_VIEW(surf_tree_view));
if (!selection || !treemodel)
  return;

/* acquire a list of selected rows */
list = gtk_tree_selection_get_selected_rows(selection, &treemodel);
for (row=list ; row ; row=g_list_next(row))
  {
  if (gtk_tree_model_get_iter(treemodel, &iter, row->data))
    {
    gtk_tree_model_get(treemodel, &iter, SURF_MODEL, &model,
                                         SURF_PLANE, &plane,
                                         SURF_SHIFT, &shift, -1);

/* perform the required operation */
    switch (type)
      {
      case CALC_SHIFTS:
        if (plane && !shift)
          {
          surf_prune_shifts(plane);
          calc_valid_shifts(model, plane);
          fn_graft_plane(plane, model);
          }
        break;

      case CALC_ENERGY:
        if (shift)
          new_ecalc_task(model, plane, shift);
        break;

      case CONV_REGIONS:
        if (shift)
          new_regcon_task(model, plane, shift);  
        break;

      case MAKE_FACES:
        if (shift)
          make_surface(model, plane, shift);  
        refresh = TRUE;
        break;
      }
    }
  }
g_list_foreach(list, (gpointer) gtk_tree_path_free, NULL);
g_list_free(list);

if (refresh)
  redraw_canvas(SINGLE);
}

/****************************************/
/* save the current surface calculation */
/****************************************/
void export_planes(gchar *name)
{
struct model_pak *data;
gchar *filename, *a, *b;

data = surfdata.surface.model;
if (!data)
  return;

/* process filename to force .gmf extension */
a = parse_strip(name); 
b = g_strconcat(a, ".gmf", NULL);

filename = g_build_filename(sysenv.cwd, b, NULL);

write_gmf(filename, data);

dialog_destroy_type(FILE_SELECT);

g_free(filename);
g_free(b);
g_free(a);
}

/******************************************/
/* apply a new morphology type to a model */
/******************************************/
void change_morph_type(struct model_pak *data, gint type)
{
gpointer camera;

g_return_if_fail(data != NULL);

data->morph_type = type;

/* save camera */
camera = camera_dup(data->camera);

morph_build(data);
model_prep(data);
coords_init(REDO_COORDS, data);

/* rescale & restore camera */
camera_rescale(data->rmax, camera);
camera_copy(data->camera, camera);
g_free(camera);

redraw_canvas(SINGLE);
}

/*****************************/
/* delete a shift in a plane */
/*****************************/
void surf_shift_delete(struct shift_pak *shift, struct plane_pak *plane)
{
g_assert(shift != NULL);
g_assert(plane != NULL);

plane->shifts = g_slist_remove(plane->shifts, shift);
shift_free(shift);
}

/*******************************************/
/* delete a plane and symmetry equivalents */
/*******************************************/
void surf_plane_delete(struct plane_pak *plane, struct model_pak *data)
{
gint n=1;
GSList *list;
struct plane_pak *pcomp;

g_assert(plane != NULL);
g_assert(data != NULL);

/* remove equivalents */
list = data->planes;
while (list)
  {
  pcomp = list->data;
  list = g_slist_next(list);

/* NB: don't delete the reference plane until the end */
  if (pcomp == plane)
    continue;

  if (facet_equiv(data, plane->index, pcomp->index))
    {
    data->planes = g_slist_remove(data->planes, pcomp);
    g_free(pcomp);
    n++;
    }
  }
/* remove main plane */
data->planes = g_slist_remove(data->planes, plane);
plane_free(plane);
}

/***************************/
/* shift deletion callback */
/***************************/
void surf_prune_selected(void)
{
GList *list, *row;
GtkTreeIter iter;
GtkTreeModel *treemodel;
GtkTreeSelection *selection;
struct model_pak *model=NULL;
struct plane_pak *plane=NULL;
struct shift_pak *shift=NULL;

/* checks */
treemodel = GTK_TREE_MODEL(surf_tree_store);
selection = gtk_tree_view_get_selection(GTK_TREE_VIEW(surf_tree_view));
if (!selection || !treemodel)
  return;

/* delete all selected shifts */
/* reverse enumeration so the tree store doesn't change when an iter is removed */
list = gtk_tree_selection_get_selected_rows(selection, &treemodel);
for (row=g_list_last(list) ; row ; row=g_list_previous(row))
  {
  if (gtk_tree_model_get_iter(treemodel, &iter, row->data))
    {
    gtk_tree_model_get(treemodel, &iter, SURF_MODEL, &model,
                                         SURF_SHIFT, &shift,
                                         SURF_PLANE, &plane, -1);
    if (plane)
      {
      if (shift)
        {
        gtk_tree_store_remove(surf_tree_store, &iter);
        surf_shift_delete(shift, plane);
        }
      else
        {
        gtk_tree_store_remove(surf_tree_store, &iter);
        surf_plane_delete(plane, model);
        }
      }
    }
  }

if (model)
  if (model->id == MORPH)
    change_morph_type(model, model->morph_type);

g_list_foreach(list, (gpointer) gtk_tree_path_free, NULL);
g_list_free(list);
}

/************************/
/* free the entire list */
/************************/
void surf_prune_all(void)
{
struct model_pak *model=surfdata.surface.model;

g_assert(model != NULL);

/* clear the tree view widget */ 
gtk_tree_store_clear(surf_tree_store);

/* free the underlying data */
plane_data_free(model->planes);
g_slist_free(model->planes);
model->planes = NULL;
model->num_planes = 0;
}

/*****************************/
/* remove all invalid shifts */
/*****************************/
void surf_prune_invalid(void)
{
gint m, n;
GtkTreeIter iter, parent;
GtkTreeModel *treemodel;
struct model_pak *model=surfdata.surface.model;
struct plane_pak *plane=NULL;
struct shift_pak *shift=NULL;

g_assert(model != NULL);
treemodel = GTK_TREE_MODEL(surf_tree_store);
if (!treemodel)
  return;

/* loop backwards over planes */
m = gtk_tree_model_iter_n_children(treemodel, NULL) - 1;
while (m >= 0)
  {
  if (!gtk_tree_model_iter_nth_child(treemodel, &parent, NULL, m))
    break;
  m--;
  gtk_tree_model_get(treemodel, &parent, SURF_PLANE, &plane, -1);

/* loop backwards over shifts in current plane */
  n = gtk_tree_model_iter_n_children(treemodel, &parent) - 1;
  while (n >= 0)
    {
    if (!gtk_tree_model_iter_nth_child(treemodel, &iter, &parent, n))
      break;
    n--;
    gtk_tree_model_get(treemodel, &iter, SURF_SHIFT, &shift, -1);

/* test dipole */
    if (fabs(shift->dipole) >= model->surface.dipole_tolerance)
      {
      gtk_tree_store_remove(surf_tree_store, &iter);
      surf_shift_delete(shift, plane);
      }
    }
  }
}

/*****************************************************/
/* callbacks to change the morphology type displayed */
/*****************************************************/
void bfdh_morph(struct model_pak *data)
{
change_morph_type(data, DHKL);
}
void equn_morph(struct model_pak *data)
{
change_morph_type(data, EQUIL_UN);
}
void grun_morph(struct model_pak *data)
{
change_morph_type(data, GROWTH_UN);
}
void eqre_morph(struct model_pak *data)
{
change_morph_type(data, EQUIL_RE);
}
void grre_morph(struct model_pak *data)
{
change_morph_type(data, GROWTH_RE);
}

/***********************************/
/* shift energy value modification */
/***********************************/
void shift_commit(GtkWidget *w, gpointer dummy)
{
gint i, empty[4];
const gchar *text[4];
gdouble value[4];
GList *list, *row;
GtkTreeIter iter;
GtkTreeModel *treemodel;
GtkTreeSelection *selection;
struct model_pak *model;
struct plane_pak *plane;
struct shift_pak *shift;

/* checks */
treemodel = GTK_TREE_MODEL(surf_tree_store);
selection = gtk_tree_view_get_selection(GTK_TREE_VIEW(surf_tree_view));
if (!selection || !treemodel)
  return;

/* acquire the energy values */
for (i=4 ; i-- ; )
  {
  text[i] = gtk_entry_get_text(GTK_ENTRY(surf_morph_energy[i]));
  if (strlen(text[i]))
    empty[i] = FALSE;
  else
    empty[i] = TRUE;
  value[i] = str_to_float(text[i]);
  }

/* modify all selected rows */
list = gtk_tree_selection_get_selected_rows(selection, &treemodel);
for (row=list ; row ; row=g_list_next(row))
  {
  if (gtk_tree_model_get_iter(treemodel, &iter, row->data))
    {
    gtk_tree_model_get(treemodel, &iter, SURF_MODEL, &model,
                                         SURF_PLANE, &plane,
                                         SURF_SHIFT, &shift, -1);
    if (model && plane && shift)
      {
      if (!empty[0])
        shift->esurf[0] = value[0];
      if (!empty[1])
        shift->eatt[0] = value[1];
      if (!empty[2])
        shift->esurf[1] = value[2];
      if (!empty[3])
        shift->eatt[1] = value[3];
      fn_surf_set_shift(&iter, shift);
      update_plane_energy(plane, model);
      }
    }
  }

if (model)
  if (model->id == MORPH)
    change_morph_type(model, model->morph_type);

g_list_foreach(list, (gpointer) gtk_tree_path_free, NULL);
g_list_free(list);
}

/************************************************/
/* added ranked faces based on dhkl to the list */
/************************************************/
#define DEBUG_ADD_RANKED_FACES 0
void add_ranked_faces(GtkWidget *w, gpointer dummy)
{
struct model_pak *data;
GSList *list;

/* checks */
data = surfdata.surface.model;
if (!data)
  return;
if (data->periodic != 3)
  return;

list = get_ranked_faces((gint) rank_value, 0.0, data);

#if DEBUG_ADD_RANKED_FACES
printf("Maximum faces: %f\n", rank_value);
printf("  Added faces: %d\n", g_slist_length(list));
#endif

data->planes = g_slist_concat(data->planes, list);

fn_graft_plane_list(list, data);
}

/*****************************************/
/* import a previous surface calculation */
/*****************************************/
void import_planes(gchar *name)
{
gchar *filename;
GSList *plist;
struct model_pak *data;
struct plane_pak *pdata;

data = surfdata.surface.model;
if (!data)
  return;
if (data->periodic != 3)
  {
  dialog_destroy_type(FILE_SELECT);
  gui_text_show(ERROR, "Source structure is not 3D periodic.");
  return;
  }

filename = g_build_filename(sysenv.cwd, name, NULL);

/* get the planes */
if (load_planes(filename, data))
  {
  gui_text_show(ERROR, "Bad planes file.");
  return;
  }

plist = data->planes;
while (plist != NULL)
  {
  pdata = plist->data;
 
/* if shift list is empty (primary plane only!) */
/* then add the result of a calc valid shifts */
/* otherwise verify the supplied shift */
  if (pdata->primary)
    {
    if (pdata->shifts)
      test_valid_shifts(pdata, data);
    else
      calc_valid_shifts(data, pdata);
    }

  plist = g_slist_next(plist);
  }

/* clean up */
dialog_destroy_type(FILE_SELECT);

fn_graft_plane_list(data->planes, data);

g_free(filename);
}

/***********************************************/
/* handle import/export file selection request */
/***********************************************/
void surf_load_planes(void)
{
file_dialog("Load planes", NULL, FILE_LOAD, (gpointer) import_planes, MORPH);
}
void surf_save_planes(void)
{
file_dialog("Save planes", NULL, FILE_SAVE, (gpointer) export_planes, MORPH);
}

/**********************************/
/* tree selection change callback */
/**********************************/
void surf_selection_changed(GtkTreeSelection *selection, gpointer data)
{
gint i, j;
gchar *text;
GList *list, *row;
GtkTreeIter iter, child;
GtkTreePath *last;
GtkTreeModel *treemodel;
struct model_pak *model=NULL;
struct plane_pak *plane=NULL;
struct shift_pak *shift=NULL;

/* checks */
treemodel = GTK_TREE_MODEL(surf_tree_store);
if (!treemodel)
  return;

/* blank the hkl family entry */
gtk_entry_set_text(GTK_ENTRY(surf_hkl_family), "");

/* only update input widget if one row selected */
list = gtk_tree_selection_get_selected_rows(selection, &treemodel);
if (g_list_length(list) == 1)
  {
  if (gtk_tree_model_get_iter(treemodel, &iter, list->data))
    {
    gtk_tree_model_get(treemodel, &iter, SURF_MODEL, &model,
                                         SURF_PLANE, &plane, 
                                         SURF_SHIFT, &shift, -1);
    g_assert(model != NULL);


/* CURRENT - print hkl family */
{
gint *m;
GSList *l1, *l2;
GString *family;

family = g_string_new(NULL);

l2 = get_facet_equiv(model, plane->index);
for (l1=l2 ; l1 ; l1=g_slist_next(l1))
  {
  m = l1->data;
  g_string_sprintfa(family, "(%d %d %d) ", m[0], m[1], m[2]);
  }

gtk_entry_set_text(GTK_ENTRY(surf_hkl_family), family->str);

g_string_free(family, TRUE);

free_slist(l2);
}


/* disallow parent selection? */
  if (model->id == MORPH)
    {
    if (plane && shift)
      {
      j=0;
      for (i=SURF_ESURF_UNRE ; i<=SURF_EATT_RE ; i++)
        {
        gtk_tree_model_get(treemodel, &iter, i, &text, -1);
        g_strstrip(text);

        if (text)
          gtk_entry_set_text(GTK_ENTRY(surf_morph_energy[j]), text);
         g_free(text);
        j++;
        }
      }
    }
  else
    {
    if (plane)
      {
      ARR3SET(surfdata.surface.miller, plane->index);
      }
    if (shift)
      {
      surfdata.surface.shift = shift->shift;
      surfdata.surface.region[0] = shift->region[0];
      surfdata.surface.region[1] = shift->region[1];
      }
/* update widgets */
    gui_relation_update(model);
    }
  }
  }

/* for each selected plane - expand row - select all shifts */
for (row=list ; row ; row=g_list_next(row))
  {
  if (gtk_tree_model_get_iter(treemodel, &iter, row->data))
    {
    gtk_tree_view_expand_row(GTK_TREE_VIEW(surf_tree_view), row->data, TRUE);
    gtk_tree_model_get(treemodel, &iter, SURF_MODEL, &model,
                                         SURF_PLANE, &plane,
                                         SURF_SHIFT, &shift, -1);
/* if (and only if) row is a plane - select self and all children */
    if (plane && !shift)
      {
      j = gtk_tree_model_iter_n_children(treemodel, &iter);
      if (j)
        {
        if (gtk_tree_model_iter_nth_child(treemodel, &child, &iter, j-1))
          {
          last = gtk_tree_model_get_path(treemodel, &child);
          gtk_tree_selection_select_range(selection, row->data, last);
          gtk_tree_path_free(last);
          }
        }
      }
    }
  }

g_list_foreach(list, (gpointer) gtk_tree_path_free, NULL);
g_list_free(list);
}

/*********************/
/* collapse all rows */
/*********************/
void surf_collapse_all(void)
{
gtk_tree_view_collapse_all(GTK_TREE_VIEW(surf_tree_view));
}

/*******************/
/* expand all rows */
/*******************/
void surf_expand_all(void)
{
gtk_tree_view_expand_all(GTK_TREE_VIEW(surf_tree_view));
}

/*******************/
/* select all rows */
/*******************/
void surf_select_all(void)
{
GtkTreeSelection *selection;

selection = gtk_tree_view_get_selection(GTK_TREE_VIEW(surf_tree_view));
if (!selection)
  return;

gtk_tree_selection_select_all(selection);
}

/**************************/
/* morphology type change */
/**************************/
void surf_morph_type(GtkWidget *w, gpointer data)
{
const gchar *text;
struct model_pak *model=data;

g_assert(model != NULL);

text = gtk_entry_get_text(GTK_ENTRY(w));

model->morph_type = DHKL;
if (g_ascii_strncasecmp(text, "Esurf (u", 8) == 0)
  model->morph_type = EQUIL_UN;
if (g_ascii_strncasecmp(text, "Esurf (r", 8) == 0)
  model->morph_type = EQUIL_RE;
if (g_ascii_strncasecmp(text, "Eatt (u", 7) == 0)
  model->morph_type = GROWTH_UN;
if (g_ascii_strncasecmp(text, "Eatt (r", 7) == 0)
  model->morph_type = GROWTH_RE;
if (g_ascii_strncasecmp(text, "broken", 6) == 0)
  model->morph_type = MORPH_BBPA;
}

/*******************************/
/* shift perturbation callback */
/*******************************/
void gui_shift_perturb(GtkWidget *w, struct model_pak *model)
{
surf_shift_explore(model, &surfdata.surface);
}

/******************/
/* MENU structure */
/******************/
static GtkItemFactoryEntry surface_menu[] = 
{
  { "/File",                     NULL, NULL, 0, "<Branch>" },
/*
  { "/File/Create morphology",   NULL, make_morph, 1, NULL },
*/
  { "/File/sep1",                NULL, NULL, 0, "<Separator>" },
  { "/File/Load planes...",      NULL, surf_load_planes, 1, NULL },
  { "/File/Save planes...",      NULL, surf_save_planes, 1, NULL },
  { "/Edit",                     NULL, NULL, 0, "<Branch>" },
  { "/Edit/Collapse all",        NULL, surf_collapse_all, 1, NULL },
  { "/Edit/Expand all",          NULL, surf_expand_all, 1, NULL },
  { "/Edit/sep1",                NULL, NULL, 0, "<Separator>" },
  { "/Edit/Delete invalid",      NULL, surf_prune_invalid, 1, NULL },
  { "/Edit/Delete selected",     NULL, surf_prune_selected, 1, NULL },
  { "/Edit/Delete all",          NULL, surf_prune_all, 1, NULL },
  { "/Edit/sep1",                NULL, NULL, 0, "<Separator>" },
  { "/Edit/Select all",          NULL, surf_select_all, 1, NULL }
};

/************************************/
/* surface creation/cut elimination */
/************************************/
void surface_dialog(void)
{
gint i, j, surface_entries;
gchar *title;
gpointer dialog, entry;
GList *list;
GtkWidget *window, *scr_win, *frame, *menu_bar;
GtkWidget *hbox1, *hbox, *vbox1, *vbox2, *vbox3, *vbox;
GtkWidget *label, *button, *spin;
GtkCellRenderer *renderer;
GtkTreeViewColumn *column;
GtkTreeSelection *select;
GtkItemFactory *item;
struct model_pak *data;

data = sysenv.active_model;

/* checks */
if (!data)
  return;
if (data->periodic != 3)
  {
  if (data->id != MORPH)
    {
    gui_text_show(ERROR, "Your model is not 3D periodic.\n");
    return;
    }
  }

/* new dialog */
dialog = dialog_request(GENSURF, "Surfaces", NULL, NULL, data);
if (!dialog)
  return;
window = dialog_window(dialog);

if (fabs(data->gulp.energy) < FRACTION_TOLERANCE && data->id != MORPH)
  {
  gui_text_show(WARNING, "Has the total energy been calculated?\n");
  }

/* FIXME - what if bonding is switched off? */
data->surface.bonds_full = g_slist_length(data->bonds);

/* init the current surface */
all_planes = FALSE;
model_init(&surfdata);
surfdata.surface.model = data;
surfdata.surface.shift = data->surface.shift;
surfdata.surface.region[0] = data->surface.region[0];
surfdata.surface.region[1] = data->surface.region[1];
ARR3SET(surfdata.surface.miller, data->surface.miller);

if (data->id == MORPH)
  title = g_strdup_printf("%s morphology", data->basename);
else
  {
  title = g_strdup_printf("%s surfaces", data->basename);

/* enforce whole molecules */
  model_colour_scheme(data->colour_scheme, data); 
  coords_init(CENT_COORDS, data);
  }

gtk_window_set_default_size(GTK_WINDOW(window), 800, 600);
gtk_window_set_title(GTK_WINDOW(window), title);
g_free(title);

/* main vbox */
vbox = gtk_vbox_new(FALSE,0);
gtk_container_add(GTK_CONTAINER(GTK_DIALOG(window)->vbox), vbox);


/* NEW - menu */
surface_entries = sizeof (surface_menu) / sizeof (surface_menu[0]);
item = gtk_item_factory_new(GTK_TYPE_MENU_BAR, "<surface>", NULL);
gtk_item_factory_create_items(item, surface_entries, surface_menu, NULL);
menu_bar = gtk_item_factory_get_widget(item, "<surface>");
gtk_box_pack_start(GTK_BOX(vbox), menu_bar, FALSE, FALSE, 0);


/* hbox for split pane */
hbox1 = gtk_hbox_new(FALSE,0);
gtk_box_pack_start(GTK_BOX(vbox), hbox1, TRUE, TRUE, 10);

if (data->id == MORPH)
  {
/* left vbox */
  vbox1 = gtk_vbox_new(FALSE, PANEL_SPACING);
  gtk_box_pack_start(GTK_BOX(hbox1), vbox1, FALSE, FALSE, PANEL_SPACING);

  frame = gtk_frame_new("Morphology type");
  gtk_box_pack_start(GTK_BOX(vbox1),frame,FALSE,FALSE,0);
 
  vbox = gtk_vbox_new(FALSE,0);
  gtk_container_add(GTK_CONTAINER(frame), vbox);
  gtk_container_set_border_width(GTK_CONTAINER(GTK_BOX(vbox)), PANEL_SPACING);
 
/* display type */
/* TODO - for loop ??? */
  new_radio_group(0, vbox, TT);
  button = add_radio_button("BFDH", (gpointer) bfdh_morph, data);
  if (data->morph_type == DHKL)
    gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(button), TRUE);
 
  button = add_radio_button("Equilibrium unrelaxed", (gpointer) equn_morph, data);
  if (data->morph_type == EQUIL_UN)
    gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(button), TRUE);
 
  button = add_radio_button("Growth unrelaxed", (gpointer) grun_morph, data);
  if (data->morph_type == GROWTH_UN)
    gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(button), TRUE);
 
  button = add_radio_button("Equilibrium relaxed", (gpointer) eqre_morph, data);
  if (data->morph_type == EQUIL_RE)
    gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(button), TRUE);
 
  button = add_radio_button("Growth relaxed", (gpointer) grre_morph, data);
  if (data->morph_type == GROWTH_RE)
    gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(button), TRUE);
 
/* editable shift values */
  frame = gtk_frame_new("Shift values");
  gtk_box_pack_start(GTK_BOX(vbox1),frame,FALSE,FALSE,0);

  vbox = gtk_vbox_new(FALSE, PANEL_SPACING);
  gtk_container_add(GTK_CONTAINER(frame), vbox);
  gtk_container_set_border_width(GTK_CONTAINER(vbox), PANEL_SPACING);
 
  j=0;
  for (i=SURF_ESURF_UNRE ; i<=SURF_EATT_RE ; i++)
    {
    hbox = gtk_hbox_new(FALSE, PANEL_SPACING);
    gtk_box_pack_start(GTK_BOX(vbox), hbox, FALSE, FALSE, 0);
    switch (i)
      {
      case SURF_ESURF_UNRE:
        label = gtk_label_new("Esurf (unrelaxed)");
        break;
      case SURF_EATT_UNRE:
        label = gtk_label_new("Eatt (unrelaxed)");
        break;
      case SURF_ESURF_RE:
        label = gtk_label_new("Esurf (relaxed)");
        break;
      case SURF_EATT_RE:
      default:
        label = gtk_label_new("Eatt (relaxed)");
        break;
      }
/* row */
    gtk_box_pack_start(GTK_BOX(hbox), label, FALSE, FALSE, 0);
    surf_morph_energy[j] = gtk_entry_new();
    gtk_box_pack_end(GTK_BOX(hbox), surf_morph_energy[j], FALSE, FALSE, 0);
    gtk_widget_set_size_request(surf_morph_energy[j], 9*sysenv.gtk_fontsize, -1);
    j++;
    }
  }
else
  {
/* left vbox */
vbox = gtk_vbox_new(FALSE,0);
gtk_box_pack_start(GTK_BOX(hbox1), vbox, FALSE, FALSE, 10);

frame = gtk_frame_new(NULL);
gtk_box_pack_start(GTK_BOX(vbox),frame,FALSE,FALSE,0);
vbox3 = gtk_vbox_new(FALSE, PANEL_SPACING);
gtk_container_add(GTK_CONTAINER(frame), vbox3);
gtk_container_set_border_width(GTK_CONTAINER(vbox3), PANEL_SPACING);

/* two vboxes; one for labels, the other for the spinners */
hbox = gtk_hbox_new(FALSE, PANEL_SPACING);
gtk_box_pack_start(GTK_BOX(vbox3), hbox, TRUE, FALSE, 0);
vbox1 = gtk_vbox_new(TRUE, PANEL_SPACING);
gtk_box_pack_start(GTK_BOX(hbox), vbox1, TRUE, FALSE, 0);
vbox2 = gtk_vbox_new(TRUE, PANEL_SPACING);
gtk_box_pack_start(GTK_BOX(hbox), vbox2, TRUE, FALSE, 0);

/* labels */
label = gtk_label_new("Miller");
gtk_box_pack_start(GTK_BOX(vbox1), label, TRUE, FALSE, 0);
label = gtk_label_new("Shift");
gtk_box_pack_start(GTK_BOX(vbox1), label, TRUE, FALSE, 0);
label = gtk_label_new("Depths");
gtk_box_pack_start(GTK_BOX(vbox1), label, TRUE, FALSE, 0);

/* miller */
hbox = gtk_hbox_new(TRUE, 0);
gtk_box_pack_start(GTK_BOX(vbox2), hbox, TRUE, FALSE, 0);

spin = gui_direct_spin(NULL, 
                 &surfdata.surface.miller[0], -99, 99, 1,
                  NULL, NULL, hbox);
gtk_widget_set_size_request(spin, 6*sysenv.gtk_fontsize, -1);

spin = gui_direct_spin(NULL, 
                 &surfdata.surface.miller[1], -99, 99, 1,
                  NULL, NULL, hbox);
gtk_widget_set_size_request(spin, 6*sysenv.gtk_fontsize, -1);

spin = gui_direct_spin(NULL, 
                 &surfdata.surface.miller[2], -99, 99, 1,
                  NULL, NULL, hbox);
gtk_widget_set_size_request(spin, 6*sysenv.gtk_fontsize, -1);

/* shift */
hbox = gtk_hbox_new(FALSE, 0);
gtk_box_pack_start(GTK_BOX(vbox2), hbox, TRUE, FALSE, 0);
spin = gui_direct_spin(NULL, 
                        &surfdata.surface.shift, 0.0, 1.0, 0.05,
                         NULL, NULL, hbox);
gtk_spin_button_set_digits(GTK_SPIN_BUTTON(spin), 4);
gtk_widget_set_size_request(spin, 6*sysenv.gtk_fontsize, -1);


/* regions */
hbox = gtk_hbox_new(FALSE,0);
gtk_box_pack_start(GTK_BOX(vbox2), hbox, TRUE, FALSE, 0);

gui_direct_spin(NULL, 
                 &surfdata.surface.region[0], 0, 99, 1,
                  NULL, NULL, hbox);

gui_direct_spin(NULL, 
                 &surfdata.surface.region[1], 0, 99, 1,
                  NULL, NULL, hbox);

/* create this surface */
gui_button(" Create ", cb_surf_create, (gpointer) data, hbox, TF);


/* frame for adding a faces */
frame = gtk_frame_new(NULL);
gtk_box_pack_start(GTK_BOX(vbox), frame, FALSE, FALSE, PANEL_SPACING);
vbox1 = gtk_vbox_new(TRUE, PANEL_SPACING);
gtk_container_add(GTK_CONTAINER(frame), vbox1);
gtk_container_set_border_width(GTK_CONTAINER(vbox1), PANEL_SPACING);

gui_button_x("Add the current surface", fn_add_plane_with_shift, NULL, vbox1);
gui_button_x("Add all valid shifts", fn_add_valid_shifts, NULL, vbox1);

hbox = gtk_hbox_new(FALSE, PANEL_SPACING);
gtk_box_pack_start(GTK_BOX(vbox1), hbox, TRUE, TRUE, 0);
gui_direct_spin("Add ", &rank_value, 1, 50, 1, NULL, NULL, hbox);
label = gtk_label_new(" Dhkl ranked faces ");
gtk_box_pack_start(GTK_BOX(hbox), label, FALSE, FALSE, 0);
gui_button_x(NULL, add_ranked_faces, NULL, hbox);

/* CURRENT */
/*
gui_button_x("Shift perturbation (exp) ", gui_shift_perturb, data, vbox1);
*/


/* frame for morphology construction */
frame = gtk_frame_new(NULL);
gtk_box_pack_start(GTK_BOX(vbox), frame, FALSE, FALSE, PANEL_SPACING);
vbox1 = gtk_vbox_new(TRUE, PANEL_SPACING);
gtk_container_add(GTK_CONTAINER(frame), vbox1);
gtk_container_set_border_width(GTK_CONTAINER(vbox1), PANEL_SPACING);

/* FIXME - this seems to be broken  */
/*
gui_button_x("Overlay morphology", surf_make_morph, data, vbox1);
*/
gui_button_x("Create morphology", make_morph, NULL, vbox1);

hbox = gtk_hbox_new(FALSE, PANEL_SPACING);
gtk_box_pack_start(GTK_BOX(vbox1), hbox, TRUE, TRUE, 0);
spin = new_spinner("Create nuclei with size ", 1.0, 1000.0, 1.0, NULL, NULL, hbox);

g_object_set_data(G_OBJECT(spin), "model", data);
gui_button_x(NULL, morph_sculpt, spin, hbox);


/* CURRENT */
gui_direct_check("Cleave to match shift (EXP)", &data->sculpt_shift_use, NULL, NULL, vbox1);




hbox = gtk_hbox_new(FALSE, PANEL_SPACING);
gtk_box_pack_start(GTK_BOX(vbox1), hbox, TRUE, TRUE, 0);

list = NULL;
list = g_list_append(list, "Dhkl");
list = g_list_append(list, "Esurf (u)");
list = g_list_append(list, "Esurf (r)");
list = g_list_append(list, "Eatt (u)");
list = g_list_append(list, "Eatt (r)");
list = g_list_append(list, "Broken bonds");

entry = gui_pulldown_new("Morphology type  ", list, FALSE, hbox);
g_signal_connect(GTK_OBJECT(entry), "changed", GTK_SIGNAL_FUNC(surf_morph_type), data);


/* frame for convergence check buttons */
  frame = gtk_frame_new("Convergence");
  gtk_box_pack_start(GTK_BOX(vbox),frame,FALSE,FALSE,PANEL_SPACING);
  vbox1 = gtk_vbox_new(FALSE,0);
  gtk_container_add(GTK_CONTAINER(frame), vbox1);
  gtk_container_set_border_width (GTK_CONTAINER(GTK_BOX(vbox1)), PANEL_SPACING);
  gui_direct_check("Attachment energy based", &data->surface.converge_eatt, NULL, NULL, vbox1);
  gui_direct_check("Converge region 1", &data->surface.converge_r1, NULL, NULL, vbox1);
  gui_direct_check("Converge region 2", &data->surface.converge_r2, NULL, NULL, vbox1);

/* frame for option check buttons */
  frame = gtk_frame_new("Construction");
  gtk_box_pack_start(GTK_BOX(vbox),frame,FALSE,FALSE,PANEL_SPACING);
  vbox1 = gtk_vbox_new(FALSE,0);
  gtk_container_add(GTK_CONTAINER(frame), vbox1);
  gtk_container_set_border_width (GTK_CONTAINER(GTK_BOX(vbox1)), PANEL_SPACING);
  gui_direct_check("Allow polar surfaces", &data->surface.include_polar,
                     NULL, NULL, vbox1);
  gui_direct_check("Allow bond cleaving", &data->surface.ignore_bonding,
                     NULL, NULL, vbox1);
  gui_direct_check("Ignore model symmetry", &data->surface.ignore_symmetry,
                     NULL, NULL, vbox1);
  gui_direct_check("Preserve depth periodicity", &data->surface.true_cell,
                     NULL, NULL, vbox1);
  gui_direct_check("Preserve atom ordering", &data->surface.keep_atom_order,
                     NULL, NULL, vbox1);

/* misc options */
  frame = gtk_frame_new(NULL);
  gtk_box_pack_start(GTK_BOX(vbox),frame,FALSE,FALSE,PANEL_SPACING);
  vbox1 = gtk_vbox_new(FALSE,0);
  gtk_container_add(GTK_CONTAINER(frame), vbox1);
  gtk_container_set_border_width(GTK_CONTAINER(GTK_BOX(vbox1)), PANEL_SPACING);
  hbox = gtk_hbox_new(FALSE, PANEL_SPACING);
  gtk_box_pack_start(GTK_BOX(vbox1), hbox, TRUE, TRUE, 0);
  gui_direct_spin("Surface dipole cutoff ", &data->surface.dipole_tolerance,
                     0.001, 99.999, 0.001, NULL, NULL, hbox);
  }

/* right vbox */
vbox1 = gtk_vbox_new(FALSE, 0);
gtk_box_pack_start(GTK_BOX(hbox1), vbox1, TRUE, TRUE, 0);

/* valid cut listing */
frame = gtk_frame_new(NULL);
gtk_box_pack_start(GTK_BOX(vbox1),frame,TRUE,TRUE,0);

vbox = gtk_vbox_new(FALSE, PANEL_SPACING);
gtk_container_add(GTK_CONTAINER(frame), vbox);
gtk_container_set_border_width(GTK_CONTAINER(GTK_BOX(vbox)), PANEL_SPACING);

/* scrolled model pane */
scr_win = gtk_scrolled_window_new (NULL, NULL);
gtk_scrolled_window_set_policy (GTK_SCROLLED_WINDOW (scr_win),
                                GTK_POLICY_AUTOMATIC, GTK_POLICY_AUTOMATIC);
gtk_box_pack_start(GTK_BOX(vbox), scr_win, TRUE, TRUE, 0);

/* planes storage */
surf_tree_store = gtk_tree_store_new(SURF_NCOLS, G_TYPE_STRING, G_TYPE_STRING,
                                     G_TYPE_STRING, G_TYPE_STRING, G_TYPE_STRING,
                                     G_TYPE_STRING, G_TYPE_STRING, G_TYPE_STRING,
                                     G_TYPE_POINTER, G_TYPE_POINTER, G_TYPE_POINTER);

/* planes viewing widget */
surf_tree_view = gtk_tree_view_new_with_model(GTK_TREE_MODEL(surf_tree_store));
gtk_scrolled_window_add_with_viewport(GTK_SCROLLED_WINDOW(scr_win), surf_tree_view);

/* set up column renderers */
for (i=0 ; i<=SURF_GNORM ; i++)
  {
  renderer = gtk_cell_renderer_text_new();
  column = gtk_tree_view_column_new_with_attributes(titles[i], renderer, "text", i, NULL);
  gtk_tree_view_append_column(GTK_TREE_VIEW(surf_tree_view), column);
  }

/* setup the selection handler */
select = gtk_tree_view_get_selection(GTK_TREE_VIEW(surf_tree_view));
gtk_tree_selection_set_mode(select, GTK_SELECTION_MULTIPLE);
g_signal_connect(G_OBJECT(select), "changed",
                 G_CALLBACK(surf_selection_changed),
                 NULL);

fn_graft_plane_list(data->planes, data);


/* CURRENT - list of hkl family planes */

surf_hkl_family = gtk_entry_new();
gtk_box_pack_start(GTK_BOX(vbox), surf_hkl_family, FALSE, FALSE, 0);


if (data->id == MORPH)
  {
/* shift operation button row */
  hbox = gtk_hbox_new(TRUE, PANEL_SPACING);
  gtk_box_pack_start(GTK_BOX(vbox), hbox, FALSE, FALSE, 0);
  gui_button("Commit ", shift_commit, NULL, hbox, TT);
  gui_button("Delete ", surf_prune_selected, data, hbox, TT);
  }
else
  {
  hbox = gtk_hbox_new(TRUE, PANEL_SPACING);
  gtk_box_pack_start(GTK_BOX(vbox), hbox, FALSE, FALSE, 0);
  gui_button("Create surface", surf_task_selected, (gpointer) MAKE_FACES, hbox, TT);
  gui_button("Converge regions", surf_task_selected, (gpointer) CONV_REGIONS, hbox, TT);
  gui_button("Calculate energy", surf_task_selected, (gpointer) CALC_ENERGY, hbox, TT);
  gui_button("Calculation setup", gulp_dialog, data, hbox, TT);
/*
  hbox = gtk_hbox_new(TRUE, PANEL_SPACING);
  gtk_box_pack_start(GTK_BOX(vbox), hbox, FALSE, FALSE, 0);
  gui_button("Find valid shifts", surf_task_selected, (gpointer) CALC_SHIFTS, hbox, TT);
  gui_button("Remove selected", surf_prune_selected, data, hbox, TT);
  gui_button("Remove invalid", surf_prune_invalid, data, hbox, TT);
  gui_button("Remove all", surf_prune_all, data, hbox, TT);
*/
  }

/* terminating button */
gui_stock_button(GTK_STOCK_CLOSE, dialog_destroy, dialog,
                   GTK_DIALOG(window)->action_area);

/* done */
gtk_widget_show_all(window);
}

