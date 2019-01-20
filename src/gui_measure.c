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
#include <stdlib.h>
#include <string.h>
#include <strings.h>

#ifndef __WIN32
#include <sys/times.h>
#endif

#include "gdis.h"
#include "coords.h"
#include "edit.h"
#include "matrix.h"
#include "opengl.h"
#include "parse.h"
#include "render.h"
#include "spatial.h"
#include "select.h"
#include "gui_shorts.h"
#include "interface.h"
#include "dialog.h"
#include "measure.h"

extern struct sysenv_pak sysenv;
extern struct elem_pak elements[];

gdouble measure_min_12 = 0.1, measure_max_12 = 2.0;
gdouble measure_min_23 = 0.1, measure_max_23 = 2.0;
gdouble measure_min_a = 0.1, measure_max_a = 180.0;
GtkWidget *meas_tv;
GtkWidget *match1, *match2, *match3, *match4, *search, *entry;
GtkWidget *spin1, *spin2, *spin3, *spin4, *spin5, *spin6;
GtkTreeStore *meas_ts=NULL;

enum {MEAS_NAME, MEAS_ATOMS, MEAS_VALUE, MEAS_MODEL, MEAS_POINTER, MEAS_NCOLS};

/************************/
/* measurement printout */
/************************/
void meas_dump_all(void)
{
struct model_pak *model;

model = sysenv.active_model;
if (model)
  measure_dump_all(model);
}

/*****************************/
/* safely acquire tree model */
/*****************************/
GtkTreeModel *meas_treemodel(void)
{
if (!meas_tv)
  return(NULL);
if (!GTK_IS_TREE_VIEW(meas_tv))
  return(NULL);
return(gtk_tree_view_get_model(GTK_TREE_VIEW(meas_tv)));
}

/*********************************************/
/* free data pointers assoc with an iterator */
/*********************************************/
void meas_free_iter_data(GtkTreeIter *iter)
{
gchar *text;
GtkTreeModel *treemodel;

treemodel = meas_treemodel();
if (!treemodel)
  return;

gtk_tree_model_get(treemodel, iter, MEAS_NAME, &text, -1);
g_free(text);
gtk_tree_model_get(treemodel, iter, MEAS_ATOMS, &text, -1);
g_free(text);
gtk_tree_model_get(treemodel, iter, MEAS_VALUE, &text, -1);
g_free(text);
}

/************************/
/* graft model iterator */
/************************/
gint meas_model_iter(GtkTreeIter *iter, struct model_pak *model)
{
GtkTreeModel *treemodel;
struct model_pak *target=NULL;

/* checks */
g_assert(model != NULL);
if (!meas_ts)
  return(0);
treemodel = meas_treemodel();
if (!treemodel)
  return(0);

/* get the first parent iterator */
if (gtk_tree_model_get_iter_first(treemodel, iter))
  {
/* search for target model, until we run out of iterators */
  gtk_tree_model_get(treemodel, iter, MEAS_MODEL, &target, -1);

  while (target != model)
    {
    if (gtk_tree_model_iter_next(treemodel, iter))
      gtk_tree_model_get(treemodel, iter, MEAS_MODEL, &target, -1);
    else
      break;
    }
  }

/* append new if model is not on the tree */
if (target != model)
  {
  gtk_tree_store_append(meas_ts, iter, NULL);
  gtk_tree_store_set(meas_ts, iter, MEAS_NAME, model->basename,
                                    MEAS_MODEL, model,
                                    MEAS_POINTER, NULL,
                                    -1);
  }
return(1);
}

/**************************************/
/* locate a single measurement's iter */
/**************************************/
void measure_iter_find(GtkTreeIter *iter, gpointer measure, struct model_pak *model)
{
/* TODO */
}

/*******************************/
/* update a single measurement */
/*******************************/
void measure_update_selected(gpointer measure, struct model_pak *model)
{
gchar *value;
GtkTreeIter iter;
GtkTreeModel *treemodel;
GtkTreeSelection *selection;

/* checks */
treemodel = meas_treemodel();
if (!treemodel)
  return;

/* get selected iterator */
if (!(selection = gtk_tree_view_get_selection(GTK_TREE_VIEW(meas_tv))))
  return;
if (!gtk_tree_selection_get_selected(selection, &treemodel, &iter))
  return;


value = g_strdup(measure_value_get(measure));

gtk_tree_store_set(meas_ts, &iter, MEAS_VALUE, value, -1);

g_free(value);
}

/******************************/
/* graft a single measurement */
/******************************/
void meas_graft_single(gpointer measure, struct model_pak *model)
{
gchar *info[3];
GtkTreeIter parent, iter;

/* checks */
g_assert(model != NULL);
g_assert(measure != NULL);

/* update dialog (if it exists) */
if (!meas_ts)
  return;
if (!meas_model_iter(&parent, model))
  return;

/* set child iterator data */
gtk_tree_store_append(meas_ts, &iter, &parent);

info[0] = measure_type_label_create(measure);
info[1] = measure_constituents_create(measure);
info[2] = g_strdup(measure_value_get(measure));

/* add the info */
gtk_tree_store_set(meas_ts, &iter, MEAS_NAME, info[0],
                                   MEAS_ATOMS, info[1],
                                   MEAS_VALUE, info[2],
                                   MEAS_MODEL, model,
                                   MEAS_POINTER, measure,
                                   -1);
g_free(info[0]);
g_free(info[1]);
g_free(info[2]);

gtk_tree_view_expand_all(GTK_TREE_VIEW(meas_tv));
}

/***********************************************************/
/* remove a model's measurements, leaving the model's iter */
/***********************************************************/
void meas_prune_all(struct model_pak *model)
{
GtkTreeIter parent, child;
GtkTreeModel *treemodel;

/* checks */
g_assert(model != NULL);
if (!meas_ts)
  return;
treemodel = meas_treemodel();
if (!treemodel)
  return;

if (meas_model_iter(&parent, model))
  {
/* remove all child iterators and free their data */ 
  while (gtk_tree_model_iter_children(treemodel, &child, &parent))
    {
    meas_free_iter_data(&child);
    gtk_tree_store_remove(meas_ts, &child); 
    }
  }
}

/**********************************/
/* graft measurements to the tree */
/**********************************/
void meas_graft_all(struct model_pak *model)
{
GSList *list;
GtkTreeModel *treemodel;
GtkTreeIter parent;
GtkTreePath *path;

/* return */
g_assert(model != NULL);
if (!meas_ts)
  return;
if (!meas_tv)
  return;

/* get model's iterator (if any) */
if (!meas_model_iter(&parent, model))
  return;

if (!g_slist_length(model->measure_list))
  {
/* no measurements - remove parent iterator */
  gtk_tree_store_remove(meas_ts, &parent); 
  }
else
  {
/* add all measurements */
  for (list=model->measure_list ; list ; list=g_slist_next(list))
    meas_graft_single(list->data, model);

/* expand model's iterator row */
  treemodel = meas_treemodel();
  if (treemodel)
    {
    path = gtk_tree_model_get_path(treemodel, &parent);
    gtk_tree_view_expand_row(GTK_TREE_VIEW(meas_tv), path, FALSE);
    gtk_tree_path_free(path);
    }
  }
}

/***********************************/
/* update measurements for a model */
/***********************************/
void meas_graft_model(struct model_pak *model)
{
g_assert(model != NULL);

/* update underlying measurements */
measure_update_global(model);
meas_prune_all(model);
meas_graft_all(model);
}

/*************************************/
/* remove model and its measurements */
/*************************************/
void meas_prune_model(struct model_pak *model)
{
GtkTreeIter parent, child;
GtkTreeModel *treemodel;

/* checks */
g_assert(model != NULL);
if (!meas_ts)
  return;
treemodel = meas_treemodel();
if (!treemodel)
  return;

if (meas_model_iter(&parent, model))
  {
/* free all child iterator data (iterators themselves are auto freed after) */
  if (gtk_tree_model_iter_children(treemodel, &child, &parent))
    {
    do
      {
      meas_free_iter_data(&child);
      }
    while (gtk_tree_model_iter_next(treemodel, &child));
    }
/* remove parent (auto removes all child iterators) */
  gtk_tree_store_remove(meas_ts, &parent); 
  }
}

/***************************/
/* geometry label deletion */
/***************************/
void meas_prune_selected(void)
{
GtkTreeIter iter;
GtkTreeModel *treemodel;
GtkTreeSelection *selection;
gpointer measurement=NULL;
struct model_pak *model=NULL;

/* checks */
treemodel = meas_treemodel();
if (!treemodel)
  return;

/* get selected iterator */
if (!(selection = gtk_tree_view_get_selection(GTK_TREE_VIEW(meas_tv))))
  return;
if (!gtk_tree_selection_get_selected(selection, &treemodel, &iter))
  return;

gtk_tree_model_get(treemodel, &iter, MEAS_MODEL, &model, -1);
gtk_tree_model_get(treemodel, &iter, MEAS_POINTER, &measurement, -1);

g_assert(model != NULL);

/* if measurement selected - delete single, else delete all */
if (measurement)
  {
/* free data */
  meas_free_iter_data(&iter);
  gtk_tree_store_remove(meas_ts, &iter); 
  measure_free(measurement, model);
  }
else
  {
  meas_prune_model(model);
  measure_free_all(model);
  }

redraw_canvas(SINGLE);
}

/*****************/
/* info on bonds */
/*****************/
void info_bond(GtkWidget *w, gint x, gint y, struct model_pak *model)
{
gpointer measure;
struct core_pak *core[2];
struct bond_pak *bond;

/* seek */
bond = gl_seek_bond(x, y, model);
if (!bond)
  return;

/* FIXME - have to ignore hydrogen bonds to stop dodgy bond measurements being drawn */
if (bond->type == BOND_HBOND)
  return;

core[0] = bond->atom1;
core[1] = bond->atom2;

measure = measure_bond_test(core, 0.0, 0.0, model);
if (measure)
  {
  measure_update_single(measure, model);
  meas_graft_single(measure, model);
  }
redraw_canvas(SINGLE);
}

/*****************/
/* info on dists */
/*****************/
void info_dist(gint x, gint y, struct model_pak *model)
{
gpointer measure;
static struct core_pak *core[2];

/* attempt to find core */
if (!(core[model->state] = gl_seek_core(x, y, model)))
  return;

/* are we looking at the 1st or the 2nd call? */
switch (model->state)
  {
  case 0:
/* select */
    select_add_core(core[0], model);
    model->state++;
    break;

  case 1:
/* remove from selection */
    select_del_core(core[0], model);
    select_del_core(core[1], model);

/* create the measurement (if it doesn't exist) */
    measure = measure_distance_test(MEASURE_DISTANCE, core, 0.0, 0.0, model);
    if (measure)
      {
      measure_update_single(measure, model);
      meas_graft_single(measure, model);
      }

    model->state--;
    break;

  default:
    g_assert_not_reached();
    model->state=0;
    return;
  }
redraw_canvas(SINGLE);
}

/******************/
/* info on angles */
/******************/
void info_angle(gint x, gint y, struct model_pak *model)
{
gpointer measure;
static struct core_pak *core[3];

/* attempt to find core */
if (!(core[model->state] = gl_seek_core(x, y, model)))
  return;

/* are we looking at the 1st or the 2nd call? */
switch (model->state)
  {
  case 0:
  case 1:
    select_add_core(core[model->state], model);
    model->state++;
    break;

  case 2:
/* remove highlighting */
    select_del_core(core[0], model);
    select_del_core(core[1], model);
    select_del_core(core[2], model);

/* create the measurement (if it doesn't exist) */
    measure = measure_angle_test(core, 0.0, 180.0, model);
    if (measure)
      {
      measure_update_single(measure, model);
      meas_graft_single(measure, model);
      }

    model->state=0;
    break;

  default:
    g_assert_not_reached();
    model->state=0;
    return;
  }
redraw_canvas(SINGLE);
}

/****************************/
/* info on torsional angles */
/****************************/
void info_torsion(gint x, gint y, struct model_pak *model)
{
gpointer measure;
static struct core_pak *core[4];

/* attempt to find core */
if (!(core[model->state] = gl_seek_core(x, y, model)))
  return;

/* what stage are we at? */
switch (model->state)
  {
  case 0:
  case 1:
  case 2:
    select_add_core(core[model->state], model);
    model->state++;
    break;

  case 3:
/* remove highlighting */
    select_del_core(core[0], model);
    select_del_core(core[1], model);
    select_del_core(core[2], model);
    select_del_core(core[3], model);

/* create the measurement (if it doesn't exist) */
    measure = measure_torsion_test(core, 0.0, 180.0, model);
    if (measure)
      {
      measure_update_single(measure, model);
      meas_graft_single(measure, model);
      }

    model->state=0;
    break;

  default:
    g_assert_not_reached();
  }
redraw_canvas(SINGLE);
}

/***************************/
/* geometry search - bonds */
/***************************/
#define DEBUG_BOND_SEARCH 0
void measure_bond_setup(struct model_pak *model)
{
const gchar *label[2];

/* checks */
g_assert(model != NULL);

label[0] = gtk_entry_get_text(GTK_ENTRY(GTK_COMBO(match1)->entry));
label[1] = gtk_entry_get_text(GTK_ENTRY(GTK_COMBO(match2)->entry));
measure_bond_search(label, measure_min_12, measure_max_12, model);
meas_graft_model(model);
}

/*******************************/
/* geometry search - distances */
/*******************************/
#define DEBUG_DIST_SEARCH 0
void measure_dist_setup(gint type, struct model_pak *model)
{
const gchar *label[2];

/* checks */
g_assert(model != NULL);

/* get pattern matching entries */
label[0] = gtk_entry_get_text(GTK_ENTRY(GTK_COMBO(match1)->entry));
label[1] = gtk_entry_get_text(GTK_ENTRY(GTK_COMBO(match2)->entry));
measure_distance_search(label, type, measure_min_12, measure_max_12, model);

meas_graft_model(model);
}

/*********************************/
/* geometry search - bond angles */
/*********************************/
#define DEBUG_ANGLE_SEARCH 0
void measure_bangle_setup(struct model_pak *model)
{
const gchar *label[3];

g_assert(model != NULL);

/* get pattern matching entries */
label[0] = gtk_entry_get_text(GTK_ENTRY(GTK_COMBO(match1)->entry));
label[1] = gtk_entry_get_text(GTK_ENTRY(GTK_COMBO(match2)->entry));
label[2] = gtk_entry_get_text(GTK_ENTRY(GTK_COMBO(match3)->entry));
measure_bangle_search(label, measure_min_a, measure_max_a, model);
meas_graft_model(model);
}

/********************************/
/* geometry search - any angles */
/********************************/
void measure_angle_setup(struct model_pak *model)
{
gdouble measure_range[6];
const gchar *label[3];

/* checks */
g_assert(model != NULL);

/* get pattern matching entries */
label[0] = gtk_entry_get_text(GTK_ENTRY(GTK_COMBO(match1)->entry));
label[1] = gtk_entry_get_text(GTK_ENTRY(GTK_COMBO(match2)->entry));
label[2] = gtk_entry_get_text(GTK_ENTRY(GTK_COMBO(match3)->entry));

VEC2SET(&measure_range[0], measure_min_12, measure_max_12);
VEC2SET(&measure_range[2], measure_min_23, measure_max_23);
VEC2SET(&measure_range[4], measure_min_a, measure_max_a);
measure_angle_search(label, measure_range, model);
meas_graft_model(model);
}

/******************************/
/* geometry search - torsions */
/******************************/
void measure_dihedral_setup(struct model_pak *data)
{
printf("Not implemented yet!\n");
}

/********************************************/
/* geometry search callback - general setup */
/********************************************/
void measure_search_setup(void)
{
const gchar *text;
struct model_pak *model;

/* checks */
model = sysenv.active_model;
if (!model)
  return;

/* get the measurement type (ensure match mask ignores irrelevant bits) */
text = gtk_entry_get_text(GTK_ENTRY(GTK_COMBO(search)->entry));
if (g_ascii_strncasecmp(text, "Bonds", 5) == 0)
  measure_bond_setup(model);
if (g_ascii_strncasecmp(text, "Distance", 8) == 0)
  measure_dist_setup(MEASURE_DISTANCE, model);
if (g_ascii_strncasecmp(text, "Intermolecular", 14) == 0)
  measure_dist_setup(MEASURE_INTER, model);
if (g_ascii_strncasecmp(text, "Bond Angles", 11) == 0)
  measure_bangle_setup(model);
if (g_ascii_strncasecmp(text, "Angles", 6) == 0)
  measure_angle_setup(model);

redraw_canvas(SINGLE);
}

/********************/
/* clean up on exit */
/********************/
void meas_cleanup(void)
{
/* NB: since the tree store is independant of the model's geom list */
/* it must be completely removed (and then restored) with the dialog */
gtk_tree_store_clear(meas_ts);
meas_tv = NULL;
meas_ts = NULL;
}

/**************************************/
/* measurement commit change callback */
/**************************************/
#define DEBUG_MEASURE_VALUE_CHANGED 0
void measure_value_changed(GtkWidget *w, gpointer dummy)
{
gdouble old, new, delta;
gdouble v[3], o[3], mat[9];
gpointer measure = NULL;
GSList *list, *list1, *list2;
GtkTreeSelection *selection;
GtkTreeModel *treemodel;
GtkTreeIter iter;
struct core_pak *core, *core1, *core2, *core3, *core4;
struct shel_pak *shel;
struct model_pak *model;

/* CURRENT */

selection = gtk_tree_view_get_selection(GTK_TREE_VIEW(meas_tv));

if (gtk_tree_selection_get_selected(selection, &treemodel, &iter))
  {
  gtk_tree_model_get(treemodel, &iter, MEAS_MODEL, &model, MEAS_POINTER, &measure, -1);
  if (measure && model)
    {
#if DEBUG_MEASURE_VALUE_CHANGED
printf("old value: [%s]\n", measure_value_get(measure));
printf("new value: [%s]\n", gtk_entry_get_text(GTK_ENTRY(entry)));
#endif

old = str_to_float(measure_value_get(measure));
new = str_to_float(gtk_entry_get_text(GTK_ENTRY(entry)));
delta = new - old;

    switch (measure_type_get(measure))
      {
      case MEASURE_BOND:

/* TODO - check that atoms are bonded */

list1 = measure_cores_get(measure);
core1 = g_slist_nth_data(list1, 0);
core2 = g_slist_nth_data(list1, 1);

ARR3SET(v, core2->x);
ARR3SUB(v, core1->x);

/* TODO - get the PBC adjusted offset vector */
/*
ARR3SUB(v, bond->offset);
*/

vecmat(model->latmat, v);
normalize(v, 3);

VEC3MUL(v, delta);

vecmat(model->ilatmat, v);


#if DEBUG_MEASURE_VALUE_CHANGED
P3VEC("offset: ", v);
printf("[%s] - [%s]\n", core1->atom_label, core2->atom_label);
#endif

connect_fragment_init(model);

list2 = connect_fragment_get(core1, core2, model);

#if DEBUG_MEASURE_VALUE_CHANGED
printf("Adjusting position of %d atoms.\n", g_slist_length(list2));
#endif

for (list=list2 ; list ; list=g_slist_next(list))
  {
  core = list->data;

/* adjust core position via delta scaled vector */
/* NB: the vector needs to take PBCs into account */

#if DEBUG_MEASURE_VALUE_CHANGED
printf(" + [%s]\n", core->atom_label);
#endif

  ARR3ADD(core->x, v);
  }

        break;

      case MEASURE_ANGLE:

list1 = measure_cores_get(measure);

core = g_slist_nth_data(list1, 0);
core1 = g_slist_nth_data(list1, 1);
core2 = g_slist_nth_data(list1, 2);

compute_normal(v, core->x, core1->x, core2->x);

/* cartesian offset to move middle core to origin */
ARR3SET(o, core1->x);
vecmat(model->latmat, o);

/* convert to radians and build an appropriate rotation matrix */
delta *= D2R;

matrix_v_rotation(mat, v, delta);

connect_fragment_init(model);

list2 = connect_fragment_get(core1, core2, model);

#if DEBUG_MEASURE_VALUE_CHANGED
printf("Adjusting position of %d atoms.\n", g_slist_length(list2));
#endif

for (list=list2 ; list ; list=g_slist_next(list))
  {
  core = list->data;

/* transform */
  vecmat(model->latmat, core->x);
  ARR3SUB(core->x, o);
  vecmat(mat, core->x);
  ARR3ADD(core->x, o);
  vecmat(model->ilatmat, core->x);

/* assoc. shell? */
  if (core->shell)
    {
    shel = core->shell;

/* transform */
    vecmat(model->latmat, shel->x);
    ARR3SUB(core->x, o);
    vecmat(mat, shel->x);
    ARR3ADD(core->x, o);
    vecmat(model->ilatmat, shel->x);
    }
  }

        break;


      case MEASURE_TORSION:

list1 = measure_cores_get(measure);
core2 = g_slist_nth_data(list1, 1);
core3 = g_slist_nth_data(list1, 2);
core4 = g_slist_nth_data(list1, 3);

ARR3SET(v, core3->x);
ARR3SUB(v, core2->x);
normalize(v, 3);

/* cartesian offset to move middle core to origin */
ARR3SET(o, core3->x);
vecmat(model->latmat, o);

/* convert to radians and build an appropriate rotation matrix */
delta *= -D2R;

matrix_v_rotation(mat, v, delta);

connect_fragment_init(model);

list2 = connect_fragment_get(core3, core4, model);

for (list=list2 ; list ; list=g_slist_next(list))
  {
  core = list->data;

/* transform */
  vecmat(model->latmat, core->x);
  ARR3SUB(core->x, o);
  vecmat(mat, core->x);
  ARR3ADD(core->x, o);
  vecmat(model->ilatmat, core->x);

/* assoc. shell? */
  if (core->shell)
    {
    shel = core->shell;

/* transform */
    vecmat(model->latmat, shel->x);
    ARR3SUB(core->x, o);
    vecmat(mat, shel->x);
    ARR3ADD(core->x, o);
    vecmat(model->ilatmat, shel->x);
    }
  }

        break;



      default:
        gui_text_show(WARNING, "Sorry, can't adjust this measurement.\n");
        return;

      }
    }
  }
if (measure && model){/*FIX c315f7*/
/* update */
coords_compute(model);
connect_refresh(model);

measure_update_single(measure, model);
measure_update_selected(measure, model);
}/*FIX c315f7*/
gui_refresh(GUI_CANVAS);
}

/*****************************/
/* manually add measurements */
/*****************************/
void meas_manual_page(GtkWidget *box)
{
GtkWidget *frame, *vbox, *hbox, *label;

/* Frame - type */
frame = gtk_frame_new (NULL);
gtk_box_pack_start(GTK_BOX(box),frame,FALSE,FALSE,0); 
gtk_container_set_border_width(GTK_CONTAINER(frame), PANEL_SPACING);

vbox = gtk_vbox_new(TRUE, PANEL_SPACING);
gtk_container_add(GTK_CONTAINER(frame), vbox);
gtk_container_set_border_width(GTK_CONTAINER(vbox), PANEL_SPACING);

gui_button_x("Measure distances",
               gtk_mode_switch, GINT_TO_POINTER(DIST_INFO), vbox);
gui_button_x("Measure bonds",
               gtk_mode_switch, GINT_TO_POINTER(BOND_INFO), vbox);
gui_button_x("Measure angles",
               gtk_mode_switch, GINT_TO_POINTER(ANGLE_INFO), vbox);
gui_button_x("Measure torsions",
               gtk_mode_switch, GINT_TO_POINTER(DIHEDRAL_INFO), vbox);

frame = gtk_frame_new (NULL);
gtk_box_pack_start(GTK_BOX(box),frame,FALSE,FALSE,0); 
gtk_container_set_border_width(GTK_CONTAINER(frame), PANEL_SPACING);

vbox = gtk_vbox_new(FALSE, PANEL_SPACING);
gtk_container_add(GTK_CONTAINER(frame), vbox);
gtk_container_set_border_width(GTK_CONTAINER(vbox), PANEL_SPACING);

hbox = gtk_hbox_new(FALSE, PANEL_SPACING);
gtk_box_pack_start(GTK_BOX(vbox),hbox,FALSE,FALSE,0); 

label = gtk_label_new("Measurement value");
gtk_box_pack_start(GTK_BOX(hbox),label,FALSE,FALSE,0); 

entry = gtk_entry_new();
gtk_box_pack_end(GTK_BOX(hbox),entry,FALSE,FALSE,0); 

g_signal_connect(GTK_OBJECT(entry), "activate", GTK_SIGNAL_FUNC(measure_value_changed), NULL);
}

/*******************************************/
/* search mode changed - update dist range */
/*******************************************/
void search_type_changed(GtkWidget *w, gpointer *data)
{
const gchar *tmp;

tmp = gtk_entry_get_text(GTK_ENTRY(GTK_COMBO(search)->entry));

/* default atom pair matching */
gtk_widget_set_sensitive(GTK_WIDGET(match3), FALSE);
gtk_widget_set_sensitive(GTK_WIDGET(spin3), FALSE);
gtk_widget_set_sensitive(GTK_WIDGET(spin4), FALSE);
gtk_widget_set_sensitive(GTK_WIDGET(spin5), FALSE);
gtk_widget_set_sensitive(GTK_WIDGET(spin6), FALSE);

if (g_ascii_strncasecmp(tmp, "Angles", 6) == 0)
  {
  gtk_widget_set_sensitive(GTK_WIDGET(match3), TRUE);
  gtk_widget_set_sensitive(GTK_WIDGET(spin3), TRUE);
  gtk_widget_set_sensitive(GTK_WIDGET(spin4), TRUE);
  gtk_widget_set_sensitive(GTK_WIDGET(spin5), TRUE);
  gtk_widget_set_sensitive(GTK_WIDGET(spin6), TRUE);
  }
if (g_ascii_strncasecmp(tmp, "Bond Angles", 11) == 0)
  {
  gtk_widget_set_sensitive(GTK_WIDGET(match3), TRUE);
  gtk_widget_set_sensitive(GTK_WIDGET(spin3), TRUE);
  gtk_widget_set_sensitive(GTK_WIDGET(spin4), TRUE);
  gtk_widget_set_sensitive(GTK_WIDGET(spin5), TRUE);
  gtk_widget_set_sensitive(GTK_WIDGET(spin6), TRUE);
  }
}

/***************/
/* search page */
/***************/
void meas_search_page(GtkWidget *box)
{
GtkWidget *frame, *hbox, *table, *button, *label;
GList *match_list=NULL, *search_list=NULL;

/* main frame */
frame = gtk_frame_new (NULL);
gtk_box_pack_start(GTK_BOX(box),frame,FALSE,FALSE,0); 
gtk_container_set_border_width(GTK_CONTAINER(frame), PANEL_SPACING/2);
hbox = gtk_hbox_new(FALSE, 0);
gtk_container_add(GTK_CONTAINER(frame), hbox);

search_list = g_list_append(search_list, "Bonds");
search_list = g_list_append(search_list, "Distances");
search_list = g_list_append(search_list, "Intermolecular");
search_list = g_list_append(search_list, "Bond Angles");
search_list = g_list_append(search_list, "Angles");

search = gtk_combo_new();
gtk_box_pack_start(GTK_BOX(hbox),search,FALSE,FALSE,0); 
gtk_entry_set_text(GTK_ENTRY(GTK_COMBO(search)->entry), "Bonds");
gtk_entry_set_editable(GTK_ENTRY(GTK_COMBO(search)->entry), FALSE);
gtk_combo_set_popdown_strings(GTK_COMBO(search), search_list);
g_signal_connect(GTK_OBJECT(GTK_COMBO(search)->entry), "changed", 
                 GTK_SIGNAL_FUNC(search_type_changed), (gpointer *) search);

button = gui_icon_button(GTK_STOCK_FIND, "Search ", measure_search_setup, NULL, NULL);
gtk_box_pack_end(GTK_BOX(hbox),button,FALSE,FALSE,0); 

/* setup frame */
frame = gtk_frame_new(NULL);
gtk_box_pack_start(GTK_BOX(box),frame,FALSE,FALSE,0); 
gtk_container_set_border_width(GTK_CONTAINER(frame), PANEL_SPACING/2);
table = gtk_table_new(3, 4, FALSE);
gtk_container_add(GTK_CONTAINER(frame), table);

/* labels for top row */
label = gtk_label_new("Atom 1");
gtk_table_attach_defaults(GTK_TABLE(table),label,0,1,0,1);
label = gtk_label_new("Atom 2");
gtk_table_attach_defaults(GTK_TABLE(table),label,1,2,0,1);
label = gtk_label_new("Atom 3");
gtk_table_attach_defaults(GTK_TABLE(table),label,2,3,0,1);

/* construct combo options */
match_list = g_list_append(match_list, "Any");
match_list = g_list_append(match_list, "Selected");

/* atom match combo boxes */
match1 = gtk_combo_new();
gtk_entry_set_text(GTK_ENTRY(GTK_COMBO(match1)->entry), "Any");
gtk_combo_set_popdown_strings(GTK_COMBO(match1), match_list);
gtk_table_attach_defaults(GTK_TABLE(table),match1,0,1,1,2);

match2 = gtk_combo_new();
gtk_entry_set_text(GTK_ENTRY(GTK_COMBO(match2)->entry), "Any");
gtk_combo_set_popdown_strings(GTK_COMBO(match2), match_list);
gtk_table_attach_defaults(GTK_TABLE(table),match2,1,2,1,2);

match3 = gtk_combo_new();
gtk_entry_set_text(GTK_ENTRY(GTK_COMBO(match3)->entry), "Any");
gtk_combo_set_popdown_strings(GTK_COMBO(match3), match_list);
gtk_table_attach_defaults(GTK_TABLE(table),match3,2,3,1,2);

/* measurement ranges */
label = gtk_label_new(" 1 - 2 cutoffs ");
gtk_table_attach_defaults(GTK_TABLE(table), label, 0, 1, 2, 3);
hbox = gtk_hbox_new(TRUE, 0);
gtk_table_attach_defaults(GTK_TABLE(table),hbox,0,1,3,4);
spin1 = gui_direct_spin(NULL, &measure_min_12, 0.1, 100.0, 0.1, NULL, NULL, NULL);
gtk_box_pack_start(GTK_BOX(hbox), spin1, TRUE, TRUE, 0);
spin2 = gui_direct_spin(NULL, &measure_max_12, 0.1, 100.0, 0.1, NULL, NULL, NULL);
gtk_box_pack_start(GTK_BOX(hbox), spin2, TRUE, TRUE, 0);

label = gtk_label_new(" 2 - 3 cutoffs ");
gtk_table_attach_defaults(GTK_TABLE(table), label, 1, 2, 2, 3);
hbox = gtk_hbox_new(TRUE, 0);
gtk_table_attach_defaults(GTK_TABLE(table),hbox,1,2,3,4);
spin3 = gui_direct_spin(NULL, &measure_min_23, 0.1, 100.0, 0.1, NULL, NULL, NULL);
gtk_box_pack_start(GTK_BOX(hbox), spin3, TRUE, TRUE, 0);
spin4 = gui_direct_spin(NULL, &measure_max_23, 0.1, 100.0, 0.1, NULL, NULL, NULL);
gtk_box_pack_start(GTK_BOX(hbox), spin4, TRUE, TRUE, 0);

label = gtk_label_new(" Angle cutoffs ");
gtk_table_attach_defaults(GTK_TABLE(table), label, 2, 3, 2, 3);
hbox = gtk_hbox_new(TRUE, 0);
gtk_table_attach_defaults(GTK_TABLE(table),hbox,2,3,3,4);
spin5 = gui_direct_spin(NULL, &measure_min_a, 0.0, 180.0, 5.0, NULL, NULL, NULL);
gtk_box_pack_start(GTK_BOX(hbox), spin5, TRUE, TRUE, 0);
spin6 = gui_direct_spin(NULL, &measure_max_a, 0.0, 180.0, 5.0, NULL, NULL, NULL);
gtk_box_pack_start(GTK_BOX(hbox), spin6, TRUE, TRUE, 0);

/* I don't like doing this, but the default size is */
/* much larger than it needs to be & wastes too much space */
gtk_widget_set_size_request(match1, 13*sysenv.gtk_fontsize, -1);
gtk_widget_set_size_request(match2, 13*sysenv.gtk_fontsize, -1);
gtk_widget_set_size_request(match3, 13*sysenv.gtk_fontsize, -1);
gtk_widget_set_size_request(search, 13*sysenv.gtk_fontsize, -1);

gtk_table_set_row_spacings(GTK_TABLE(table), PANEL_SPACING);
gtk_table_set_col_spacings(GTK_TABLE(table), PANEL_SPACING);

/* default measurement is bonds */
gtk_widget_set_sensitive(GTK_WIDGET(match3), FALSE);
gtk_widget_set_sensitive(GTK_WIDGET(spin3), FALSE);
gtk_widget_set_sensitive(GTK_WIDGET(spin4), FALSE);
gtk_widget_set_sensitive(GTK_WIDGET(spin5), FALSE);
gtk_widget_set_sensitive(GTK_WIDGET(spin6), FALSE);
}

/***********************************/
/* select a particular measurement */
/***********************************/
void measure_select(gpointer measurement, struct model_pak *model)
{
GSList *list;

for (list=model->measure_list ; list ; list=g_slist_next(list))
  {
  if (measurement == list->data)
    {
/* selected measurement */
    measure_colour_set(1.0, 1.0, 0.0, list->data);
    }
  else
    {
/* reset other measurements */
    measure_colour_set(0.8, 0.8, 0.8, list->data);
    }
  }
redraw_canvas(SINGLE);
}

/**********************************/
/* measurement selection callback */
/**********************************/
void measure_select_changed(GtkTreeSelection *selection, gpointer data)
{
GtkTreeModel *treemodel;
GtkTreeIter iter;
gpointer measure;
struct model_pak *model;

gtk_entry_set_text(GTK_ENTRY(entry), "");

if (gtk_tree_selection_get_selected(selection, &treemodel, &iter))
  {
  gtk_tree_model_get(treemodel, &iter, MEAS_MODEL, &model, -1);
  gtk_tree_model_get(treemodel, &iter, MEAS_POINTER, &measure, -1);
  if (measure && model)
    {
    measure_select(measure, model);

    gtk_entry_set_text(GTK_ENTRY(entry), measure_value_get(measure));

/* CURRENT - only allow adjustments to supported editables */
    switch (measure_type_get(measure))
      {
      case MEASURE_BOND:
      case MEASURE_ANGLE:
      case MEASURE_TORSION:
        gtk_widget_set_sensitive(GTK_WIDGET(entry), TRUE);
        break;

      default:
        gtk_widget_set_sensitive(GTK_WIDGET(entry), FALSE);
        break;
      }
    }
  }
}

/**************************/
/* geometric measurements */
/**************************/
void gui_measure_dialog(void)
{
gint i;
gpointer dialog;
gchar *titles[] = {"  Name  ", " Constituent atoms ", " Value "};
GSList *list;
GtkWidget *window, *swin, *notebook;
GtkWidget *frame, *vbox, *hbox, *label;
GtkTreeSelection *select;
GtkCellRenderer *renderer;
GtkTreeViewColumn *column;
struct model_pak *model;

/* request dialog */
dialog = dialog_request(GEOMETRY, "Measurements", NULL, meas_cleanup, NULL);
if (!dialog)
  return;
window = dialog_window(dialog);
gtk_window_set_default_size(GTK_WINDOW(window), 250, 550);

/* notebook */
notebook = gtk_notebook_new();
gtk_notebook_set_tab_pos(GTK_NOTEBOOK(notebook), GTK_POS_TOP);
gtk_box_pack_start(GTK_BOX(GTK_DIALOG(window)->vbox),notebook,FALSE,FALSE,0);
gtk_notebook_set_show_border(GTK_NOTEBOOK(notebook), TRUE);

/* manual measurement */
vbox = gtk_vbox_new(FALSE, PANEL_SPACING);
label = gtk_label_new(" Manual ");
gtk_notebook_append_page(GTK_NOTEBOOK(notebook), vbox, label);
meas_manual_page(vbox);

/* searching */
vbox = gtk_vbox_new(FALSE, PANEL_SPACING);
label = gtk_label_new(" Search ");
gtk_notebook_append_page(GTK_NOTEBOOK(notebook), vbox, label);
meas_search_page(vbox);

/* measurements viewing pane */
frame = gtk_frame_new ("Label list");
gtk_box_pack_start(GTK_BOX(GTK_DIALOG(window)->vbox),frame,TRUE,TRUE,0); 
gtk_container_set_border_width(GTK_CONTAINER(frame), PANEL_SPACING/2);
vbox = gtk_vbox_new(FALSE, 0);
gtk_container_add(GTK_CONTAINER(frame), vbox);
gtk_container_set_border_width(GTK_CONTAINER(GTK_BOX(vbox)), PANEL_SPACING);
swin = gtk_scrolled_window_new (NULL, NULL);
gtk_scrolled_window_set_policy(GTK_SCROLLED_WINDOW(swin),
                               GTK_POLICY_AUTOMATIC, GTK_POLICY_AUTOMATIC);
gtk_box_pack_start(GTK_BOX(vbox), swin, TRUE, TRUE, 0);

/* NEW -  tree view */
/* underlying data storage */
/* 3 strings - data, 1 pointer - the measurement itself */
/* NB: model assoc. is passed with the selection callback */
if (!meas_ts)
  meas_ts = gtk_tree_store_new(MEAS_NCOLS, 
                               G_TYPE_STRING,
                               G_TYPE_STRING,
                               G_TYPE_STRING,
                               G_TYPE_POINTER,
                               G_TYPE_POINTER);

meas_tv = gtk_tree_view_new_with_model(GTK_TREE_MODEL(meas_ts));
gtk_scrolled_window_add_with_viewport(GTK_SCROLLED_WINDOW(swin), meas_tv);

select = gtk_tree_view_get_selection(GTK_TREE_VIEW(meas_tv));
g_signal_connect(G_OBJECT(select), "changed",
                 G_CALLBACK(measure_select_changed), NULL);

/* set up column renderers */
for (i=0 ; i<=MEAS_VALUE ; i++)
  {
  renderer = gtk_cell_renderer_text_new();
  column = gtk_tree_view_column_new_with_attributes
              (titles[i], renderer, "text", i, NULL);
  gtk_tree_view_append_column(GTK_TREE_VIEW(meas_tv), column);
  }

/* list modification buttons */
hbox = gtk_hbox_new(TRUE, 10);
gtk_box_pack_start(GTK_BOX(vbox), hbox, FALSE, FALSE, PANEL_SPACING);

gui_button(" Select ", measure_select_all, NULL, hbox, TT);
gui_button(" Delete ", meas_prune_selected, NULL, hbox, TT);
gui_button(" Dump ", meas_dump_all, NULL, hbox, TT);

/* terminating button */
gui_stock_button(GTK_STOCK_CLOSE, dialog_destroy, dialog,
                   GTK_DIALOG(window)->action_area);

/* done */
gtk_widget_show_all(window);

/* refresh all measurements */
for (list=sysenv.mal ; list ; list=g_slist_next(list))
  {
  model = list->data;
  meas_graft_model(model);
  }
gtk_tree_view_expand_all(GTK_TREE_VIEW(meas_tv));
}

