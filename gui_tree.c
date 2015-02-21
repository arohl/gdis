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

#include "gdis.h"
#include "coords.h"
#include "model.h"
#include "space.h"
#include "graph.h"
#include "select.h"
#include "matrix.h"
#include "project.h"
#include "gui_shorts.h"
#include "interface.h"
#include "dialog.h"

#include "methane.xpm"
#include "box.xpm"
#include "surface.xpm"
#include "polymer.xpm"
#include "diamond2.xpm"
#include "graph.xpm"

/* top level data structure */
extern struct sysenv_pak sysenv;

/**********************/
/* model tree globals */
/**********************/
enum
{
TREE_PIXMAP,
TREE_NAME,
TREE_POINTER,
TREE_DATA,
TREE_NCOLS
};

/*****************************************************/
/* select new iter, assuming current will be deleted */
/*****************************************************/
void tree_select_next(GtkTreeIter *iter)
{
GtkTreeIter next;
GtkTreePath *treepath;
GtkTreeModel *treemodel;
GtkTreeSelection *selection;

selection = gtk_tree_view_get_selection(GTK_TREE_VIEW(sysenv.tree)); 
if (!selection)
  return;

/* attempt to select previous iterator */
treemodel = GTK_TREE_MODEL(sysenv.tree_store);
treepath = gtk_tree_model_get_path(treemodel, iter);
if (gtk_tree_path_prev(treepath))
  {
  if (gtk_tree_model_get_iter(treemodel, &next, treepath))
    gtk_tree_selection_select_iter(selection, &next);
  }
else
  {
/* attempt to select next iterator, else - select parent */
  next = *iter;
  if (gtk_tree_model_iter_next(treemodel, &next))
    gtk_tree_selection_select_iter(selection, &next);
  else
    {
    if (gtk_tree_model_iter_parent(treemodel, &next, iter))
      gtk_tree_selection_select_iter(selection, &next);
    else
      sysenv.active_model = NULL;
    }
  }
gtk_tree_path_free(treepath);
}

/********************************/
/* replacement deletion routine */
/********************************/
void tree_select_delete(void)
{
gint depth;
gpointer data;
GtkTreeIter iter, next;
GtkTreeSelection *selection;
GtkTreeModel *treemodel;
struct model_pak *model;

selection = gtk_tree_view_get_selection(GTK_TREE_VIEW(sysenv.tree)); 
if (!selection)
  return;
if (!gtk_tree_selection_get_selected(selection, &treemodel, &iter))
  return;

depth = gtk_tree_store_iter_depth(sysenv.tree_store, &iter);

gtk_tree_model_get(GTK_TREE_MODEL(sysenv.tree_store), &iter,
                   TREE_POINTER, &model,
                   TREE_DATA, &data,
                   -1);

/* alter the selection iter and then remove old iter */
next = iter;
tree_select_next(&next);
gtk_tree_store_remove(sysenv.tree_store, &iter);

/* data pointer cleanup */
switch (depth)
  {
  case 0:
    if (model)
      model_delete(model);
    break;

  case 1:
    if (model && data)
      graph_free(data, model);
    break;

  default:
    g_assert_not_reached();
  }
/*
tree_select_model(sysenv.active_model);
*/
}

/**************************/
/* tree traverse function */
/**************************/
gint func(GtkTreeModel *model, GtkTreePath *path, GtkTreeIter *iter, gpointer data)
{
GtkTreeSelection *selection;
struct model_pak *mdata;

/* get the model pointer from the tree store */
gtk_tree_model_get(model, iter, TREE_POINTER, &mdata, -1);

/* if it matches the model supplied, select it */
if (data == mdata)
  {
  selection = gtk_tree_view_get_selection(GTK_TREE_VIEW(sysenv.tree)); 
  gtk_tree_selection_select_iter(selection, iter);
/* return TRUE to indicate we're done  */
  return(TRUE);
  }
/* no match, keep traversing */
return(FALSE);
}

/******************************/
/* select the specified model */
/******************************/
void tree_select_model(struct model_pak *model)
{
if (!model)
  return;
sysenv.active_model = model;

gtk_tree_model_foreach(GTK_TREE_MODEL(sysenv.tree_store), &func, model);
}

/***************************/
/* select the active model */
/***************************/
void tree_select_active(void)
{
/* TODO - clear if no models */
tree_select_model(sysenv.active_model);
}

/**************************/
/* get model's iterator */
/**************************/
gboolean tree_model_iter(GtkTreeIter *iter, gpointer model)
{
gpointer data;

if (!gtk_tree_model_get_iter_first(GTK_TREE_MODEL(sysenv.tree_store), iter))
  return(FALSE);
do
  {
  gtk_tree_model_get(GTK_TREE_MODEL(sysenv.tree_store), iter, TREE_POINTER, &data, -1);
  if (data == model)
    return(TRUE);
  }
while (gtk_tree_model_iter_next(GTK_TREE_MODEL(sysenv.tree_store), iter));

return(FALSE);
}

/*******************/
/* refresh a model */
/*******************/
void tree_model_refresh(struct model_pak *model)
{
gboolean flag=FALSE;
GSList *list;
GtkTreeIter root, branch, active;
GdkPixbuf *pixbuf;
GtkTreeSelection *selection;

/* checks */
if (!model)
  return;
if (!tree_model_iter(&root, model))
  return;

selection = gtk_tree_view_get_selection(GTK_TREE_VIEW(sysenv.tree)); 

/* update name */
gtk_tree_store_set(sysenv.tree_store, &root, TREE_NAME, model->basename, -1);

/* update any graphs */
pixbuf = gdk_pixbuf_new_from_xpm_data((const gchar **) graph_xpm);
for (list=model->graph_list ; list ; list=g_slist_next(list))
  {
  if (graph_grafted(list->data))
    continue;

  gtk_tree_store_append(sysenv.tree_store, &branch, &root);
  gtk_tree_store_set(sysenv.tree_store, &branch,
                     TREE_PIXMAP, pixbuf, 
                     TREE_NAME, graph_treename(list->data),
                     TREE_POINTER, model,
                     TREE_DATA, list->data,
                     -1);
  graph_set_grafted(TRUE, list->data);

  if (selection)
    {
    active = branch;
    flag = TRUE;
    }
  }
/*
gdk_pixbuf_unref(pixbuf);
*/

/* NB: have to expand first and THEN select - as by default */
/* appending to the tree store is done in unexpanded fashion */
gtk_tree_view_expand_all(GTK_TREE_VIEW(sysenv.tree));
if (flag)
  gtk_tree_selection_select_iter(selection, &active);
}

/**************************************/
/* add (or refresh if exists) a model */
/**************************************/
void tree_model_add(struct model_pak *model)
{
GSList *list;
GtkTreeIter root, branch;
GdkPixbuf *pixbuf;

/* checks */
g_assert(model != NULL);
if (model->grafted)
  {
  tree_model_refresh(model);
  return;
  }

/* graft the model */
gtk_tree_store_append(sysenv.tree_store, &root, NULL);
model->grafted = TRUE;

/* setup pixmap */
if (model->id == MORPH)
  pixbuf = gdk_pixbuf_new_from_xpm_data((const gchar **) diamond2_xpm);
else
  {
  switch(model->periodic)
    {
    case 3:
      pixbuf = gdk_pixbuf_new_from_xpm_data((const gchar **) box_xpm);
      break;
    case 2:
      pixbuf = gdk_pixbuf_new_from_xpm_data((const gchar **) surface_xpm);
      break;
    case 1:
      pixbuf = gdk_pixbuf_new_from_xpm_data((const gchar **) polymer_xpm);
      break;
    default:
      pixbuf = gdk_pixbuf_new_from_xpm_data((const gchar **) methane_xpm);
      break;
    }
  }

/* set the parent iterator data */
gtk_tree_store_set(sysenv.tree_store, &root,
                   TREE_PIXMAP, pixbuf, 
                   TREE_NAME, model->basename,
                   TREE_POINTER, model,
                   TREE_DATA, NULL,
                   -1);
/*
gdk_pixbuf_unref(pixbuf);
*/

/* add any attached graphs */
pixbuf = gdk_pixbuf_new_from_xpm_data((const gchar **) graph_xpm);
for (list=model->graph_list ; list ; list=g_slist_next(list))
  {
  gtk_tree_store_append(sysenv.tree_store, &branch, &root);

/* TODO - better naming eg graph_1, graph_2 etc */
/* eg a way to convert graph pointer to a number? (that way stays the same) */
  gtk_tree_store_set(sysenv.tree_store, &branch,
                     TREE_PIXMAP, pixbuf, 
                     TREE_NAME, graph_treename(list->data),
                     TREE_POINTER, model,
                     TREE_DATA, list->data,
                     -1);
  graph_set_grafted(TRUE, list->data);

/* if we have a current graph - make it the selection */
/* NB: graph_flag should not be set for a REPLACE with a name change */
/*
  if (model->graph_active == list->data)
    {
    iter = branch;
    active_flag = TRUE;
    }
*/
  }
/*
gdk_pixbuf_unref(pixbuf);
*/

/* add any attached pictures */
/*
for (list=model->picture_list ; list ; list=g_slist_next(list))
  {
  gtk_tree_store_append(sysenv.tree_store, &branch, &root);
  gtk_tree_store_set(sysenv.tree_store, &branch,
                     TREE_NAME, g_strdup(list->data),
                     TREE_POINTER, list->data,
                     -1);

  if (model->picture_active == list->data)
    {
    iter = branch;
    active_flag = TRUE;
    }
  }
*/
}

/********************************/
/* handle tree selection events */
/********************************/
#define DEBUG_TREE_SELECTION_CHANGED 0
void tree_selection_changed(GtkTreeSelection *selection, gpointer data)
{
gint depth;
gpointer graph;
GtkTreeIter iter;
GtkTreeModel *treemodel;
struct model_pak *model=NULL;

if (gtk_tree_selection_get_selected(selection, &treemodel, &iter))
  {
  depth = gtk_tree_store_iter_depth(sysenv.tree_store, &iter);
  switch (depth)
    {
    case 0:
      gtk_tree_model_get(treemodel, &iter, TREE_POINTER, &model, -1);
      if (model)
        gui_model_select(model);
      break;

    case 1:
      gtk_tree_model_get(treemodel, &iter, TREE_POINTER, &model, TREE_DATA, &graph, -1);
      if (graph)
        {
        gui_model_select(model);
        model->graph_active = graph;
        redraw_canvas(SINGLE);
        }
      break;

    default:
      g_assert_not_reached();
    }
  }
}

/*********************************/
/* create the model viewing pane */
/*********************************/
void tree_init(GtkWidget *box)
{
GtkWidget *swin;
GtkCellRenderer *renderer;
GtkTreeViewColumn *column;
GtkTreeSelection *select;

/* scrolled window for the model pane */
swin = gtk_scrolled_window_new(NULL, NULL);
gtk_scrolled_window_set_policy(GTK_SCROLLED_WINDOW(swin),
                               GTK_POLICY_AUTOMATIC,
                               GTK_POLICY_AUTOMATIC);
gtk_container_add(GTK_CONTAINER(box), swin);

/* underlying data storage */
sysenv.tree_store = gtk_tree_store_new(TREE_NCOLS,
                                       GDK_TYPE_PIXBUF,
                                       G_TYPE_STRING,
                                       G_TYPE_POINTER,
                                       G_TYPE_POINTER);
/* actual tree widget */
sysenv.tree = gtk_tree_view_new_with_model(GTK_TREE_MODEL(sysenv.tree_store));
gtk_scrolled_window_add_with_viewport(GTK_SCROLLED_WINDOW(swin), sysenv.tree);

/* setup the model pixmap rendering colum */
renderer = gtk_cell_renderer_pixbuf_new();
column = gtk_tree_view_column_new_with_attributes("image", renderer,
                                                  "pixbuf", TREE_PIXMAP, NULL);
gtk_tree_view_append_column(GTK_TREE_VIEW(sysenv.tree), column);

/* setup the model name rendering colum */
renderer = gtk_cell_renderer_text_new();
column = gtk_tree_view_column_new_with_attributes("name", renderer,
                                                  "text", TREE_NAME, NULL);
gtk_tree_view_append_column(GTK_TREE_VIEW(sysenv.tree), column);

/* setup the selection handler */
select = gtk_tree_view_get_selection(GTK_TREE_VIEW(sysenv.tree));
gtk_tree_selection_set_mode(select, GTK_SELECTION_BROWSE);
g_signal_connect(G_OBJECT(select), "changed",
                 G_CALLBACK(tree_selection_changed),
                 NULL);

gtk_tree_view_set_headers_visible(GTK_TREE_VIEW(sysenv.tree), FALSE);

gtk_widget_show(swin);
}

