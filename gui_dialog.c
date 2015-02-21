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
#include <unistd.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <gdk/gdk.h>
#include <gtk/gtk.h>

#ifndef __WIN32
#include <sys/times.h>
#endif

#include "gdis.h"
#include "file.h"
#include "matrix.h"
#include "parse.h"
#include "gui_shorts.h"
#include "interface.h"
#include "dialog.h"

/* top level data structure */
extern struct sysenv_pak sysenv;

extern GtkWidget *window;

GSList *dir_list=NULL;
GtkListStore *file_path_ts=NULL, *file_name_ts=NULL;
GtkWidget *file_path_tv=NULL, *file_name_tv=NULL;
GtkWidget *curr_path, *file_name;

enum {FILE_PATH, FILE_PATH_NCOLS};
enum {FILE_NAME, FILE_NAME_NCOLS};

/* replacement for dialog_pak */
struct dialog_pak
{
gint type;
gpointer model;
gpointer data;
GHashTable *children;
void (*refresh)(gpointer);
void (*cleanup)(gpointer);
GtkWidget *window;
};

/************************************/
/* callback for destroying a dialog */
/************************************/
void destroy_widget(GtkWidget *w, gpointer target)
{
if (GTK_IS_WIDGET(target))
  gtk_widget_destroy(target);
}

/*******************************/
/* dialog extraction primitive */
/*******************************/
gpointer dialog_window(gpointer data)
{
struct dialog_pak *dialog = data;

return(dialog->window);
}

gpointer dialog_model(gpointer data)
{
struct dialog_pak *dialog = data;

return(dialog->model);
}

/************************************/
/* dialog data extraction primitive */
/************************************/
gpointer dialog_child_get(gpointer data, const gchar *key)
{
struct dialog_pak *dialog = data;
gpointer child = NULL;

if (dialog->children)
  child = g_hash_table_lookup(dialog->children, key);

return(child);
}

/**********************************/
/* dialog data addition primitive */
/**********************************/
/* used for associating a label with a child widget pointer within a dialog */
void dialog_child_set(gpointer data, const gchar *key, gpointer value)
{
struct dialog_pak *dialog = data;

if (!dialog->children)
  dialog->children = g_hash_table_new_full(g_str_hash, g_str_equal, g_free, NULL);

/* use replace, so old g_strdup() key gets freed */
g_hash_table_replace(dialog->children, g_strdup(key), value);
}

/******************/
/* debugging dump */
/******************/
void dialog_dump(void)
{
gint i=0;
struct dialog_pak *dialog;
GSList *list;

for (list=sysenv.dialog_list ; list ; list=g_slist_next(list))
  {
  dialog = list->data;
  printf("[%d] model: %p, type: %d\n", i++, dialog->model, dialog->type);
  }
}

/********************************************/
/* cleanup when a dialog has been destroyed */
/********************************************/
#define DEBUG_DIALOG_FREE 0
void dialog_free(GtkWidget *w, gpointer data)
{
struct dialog_pak *dialog = data;

#if DEBUG_DIALOG_FREE
printf("freeing dialog : %p...\n", data);
#endif

if (dialog->cleanup)
  dialog->cleanup(dialog->model);

if (dialog->children)
  g_hash_table_destroy(dialog->children); 

sysenv.dialog_list = g_slist_remove(sysenv.dialog_list, dialog);
g_free(dialog);
}

/***************************/
/* dialog destroy callback */
/***************************/
void dialog_destroy(GtkWidget *w, gpointer data)
{
struct dialog_pak *dialog = data;

gtk_widget_destroy(dialog->window);
}

/***************************/
/* destroy dialogs by type */
/***************************/
void dialog_destroy_type(gint type)
{
GSList *list;
struct dialog_pak *dialog;

list = sysenv.dialog_list;
while (list)
  {
  dialog = list->data;
  list = g_slist_next(list);

  if (dialog->type == type)
    dialog_destroy(NULL, dialog);
  }
}

/***********************************************/
/* destroy all dialogs associated with a model */
/***********************************************/
void dialog_destroy_model(struct model_pak *model)
{
GSList *list;
struct dialog_pak *dialog;

list = sysenv.dialog_list;
while (list)
  {
  dialog = list->data;
  list = g_slist_next(list);

  if (dialog->model == model)
    dialog_destroy(NULL, dialog);
  }
}

/*******************************************************/
/* destroy a given dialog type associated with a model */
/*******************************************************/
void dialog_destroy_single(gint type, struct model_pak *model)
{
GSList *list;
struct dialog_pak *dialog;

list = sysenv.dialog_list;
while (list)
  {
  dialog = list->data;
  list = g_slist_next(list);

  if (dialog->type == type && dialog->model == model)
    dialog_destroy(NULL, dialog);
  }
}

/*****************************/
/* sanity check for a dialog */
/*****************************/
gint dialog_valid(gpointer dialog)
{
if (g_slist_find(sysenv.dialog_list, dialog))
  return(TRUE);
return(FALSE);
}

/********************************************/
/* test if a specific type of dialog exists */
/********************************************/
gint dialog_exists(gint type, struct model_pak *model)
{
GSList *list;
struct dialog_pak *dialog;

for (list=sysenv.dialog_list ; list ; list=g_slist_next(list))
  {
  dialog = list->data;

  if (model)
    {
    if (dialog->model == model && dialog->type == type)
      return(TRUE);
    }
  else
    {
    if (dialog->type == type)
      return(TRUE);
    }
  }
return(FALSE);
}

/******************************/
/* update all current dialogs */ 
/******************************/
void dialog_refresh_all(void)
{
GSList *list;
struct dialog_pak *dialog;

for (list=sysenv.dialog_list ; list ; list=g_slist_next(list))
  {
  dialog = list->data;
  if (dialog->refresh)
    dialog->refresh(dialog);
  }
}

/***********************************/
/* update dialog of specified type */
/***********************************/
void dialog_refresh_type(gint type)
{
GSList *list;
struct dialog_pak *dialog;

for (list=sysenv.dialog_list ; list ; list=g_slist_next(list))
  {
  dialog = list->data;

  if (dialog->type == type)
    if (dialog->refresh)
      dialog->refresh(dialog);
  }
}

/*************************************/
/* update dialogs of specified model */
/*************************************/
void dialog_refresh_model(struct model_pak *model)
{
GSList *list;
struct dialog_pak *dialog;

for (list=sysenv.dialog_list ; list ; list=g_slist_next(list))
  {
  dialog = list->data;

  if (dialog->model == model)
    if (dialog->refresh)
      dialog->refresh(dialog);
  }
}

/*****************************/
/* dialog creation primitive */
/*****************************/
#define DEBUG_DIALOG_REQUEST 0
gpointer dialog_request(gint type,
                        gchar *title,
                        gpointer refresh,
                        gpointer cleanup,
                        struct model_pak *model)
{
GSList *list;
struct dialog_pak *dialog;

#if DEBUG_DIALOG_REQUEST
printf("request to open dialog type : %d, model : %p\n", type, model);
#endif

/* is another of the same type active? */
for (list=sysenv.dialog_list ; list ; list=g_slist_next(list))
  {
  dialog = list->data;

/* allow only one instance of these dialogs (eg because of globals) */
  switch (type)
    {
    case CREATOR:
    case DOCKING:
    case GENSURF:
    case SURF:
      if (type == dialog->type)
        {
        gdk_window_show((dialog->window)->window);
        return(NULL);
        }
    }

/* type dependent testing */
  if (dialog->type == type && dialog->model == model)
    {
#if DEBUG_DIALOG_REQUEST
printf(" - Exists...\n");
#endif
/* dialog already open - raise it to the fore & select assoc. model */
    gdk_window_show((dialog->window)->window);
    if (model)
      gui_model_select(model);
    return(NULL);
    }
  }
#if DEBUG_DIALOG_REQUEST
printf(" - Creating...\n");
#endif

/* dialog setup */
dialog = g_malloc(sizeof(struct dialog_pak));
sysenv.dialog_list = g_slist_prepend(sysenv.dialog_list, dialog);
dialog->type = type;
dialog->model = model;
dialog->data = NULL;
dialog->children = NULL;
dialog->refresh = refresh;
dialog->cleanup = cleanup;

/* widget setup */
dialog->window = gtk_dialog_new();
gtk_window_set_title(GTK_WINDOW(dialog->window), title);
g_signal_connect(GTK_OBJECT(dialog->window), "destroy",
                 GTK_SIGNAL_FUNC(dialog_free), (gpointer) dialog);

return(dialog);
}

/******************************/
/* safely acquire tree models */
/******************************/
GtkTreeModel *file_path_treemodel(void)
{
if (!file_path_tv)
  return(NULL);
if (!GTK_IS_TREE_VIEW(file_path_tv))
  return(NULL);
return(gtk_tree_view_get_model(GTK_TREE_VIEW(file_path_tv)));
}

GtkTreeModel *file_name_treemodel(void)
{
if (!file_name_tv)
  return(NULL);
if (!GTK_IS_TREE_VIEW(file_name_tv))
  return(NULL);
return(gtk_tree_view_get_model(GTK_TREE_VIEW(file_name_tv)));
}

/***************************************/
/* update contents of file load widget */
/***************************************/
#define DEBUG_UPDATE_FILE_PANE 0
gint update_file_pane(void)
{
gint filter;
gchar *name, *full;
GSList *list;
GtkTreeIter iter;
struct stat buff;

/* NEW */
filter = sysenv.file_type;

/* getting from this directory */
gtk_label_set_text(GTK_LABEL(curr_path), sysenv.cwd);

/* clear old data */
gtk_list_store_clear(file_path_ts);
gtk_list_store_clear(file_name_ts);
free_slist(dir_list);

/* get directory listing */
/* NB: success will always return something, even if only ".." */
dir_list = file_dir_list(sysenv.cwd, TRUE);
if (!dir_list)
  {
  gui_text_show(ERROR, "No directory listing; check your permissions.\n");
  gtk_list_store_append(file_path_ts, &iter);
  gtk_list_store_set(file_path_ts, &iter, FILE_PATH, "..", -1);
  }

list = dir_list;
while (list)
  {
/* stat the file (MUST have the full path) */
  name = (gchar *) list->data;
  full = g_build_filename(sysenv.cwd, list->data, NULL);
  stat(full, &buff);

#if DEBUG_UPDATE_FILE_PANE
printf("[%s] : %d\n",full,buff.st_mode);
#endif

/* convert and check if directory */
  if ((buff.st_mode & S_IFMT) ==  S_IFDIR)
    {
    gtk_list_store_append(file_path_ts, &iter);
    gtk_list_store_set(file_path_ts, &iter, FILE_PATH, name, -1);
    }
  else
    {
/* is it a recognized type? */
    if (file_extension_valid(list->data))
      {
      gtk_list_store_append(file_name_ts, &iter);
      gtk_list_store_set(file_name_ts, &iter, FILE_NAME, name, -1);
      }

#if OLD_GOK /* FIXME: can it be deleted? I haven't changed this code 
               below to handle empty extensions */
    file_data = get_file_info((gpointer *) list->data, BY_EXTENSION);
    if (file_data)
      {
/* display name if matches current file filter */
      if (filter == DATA) 
        {
/* don't add filetypes with no read/write routines - pictures */
        if (file_data->read_file || file_data->write_file)
          {
          gtk_list_store_append(file_name_ts, &iter);
          gtk_list_store_set(file_name_ts, &iter, FILE_NAME, name, -1);
          }
        }
      else
        {
/* add if specifically requested (ie even if not in menu) */
        if (filter == file_data->group)
          {
          gtk_list_store_append(file_name_ts, &iter);
          gtk_list_store_set(file_name_ts, &iter, FILE_NAME, name, -1);
          }
        }
      }
#endif

    }
  list = g_slist_next(list);
  g_free(full);
  }
return(TRUE);
}

/*****************************/
/* file type/filter handlers */
/*****************************/
void type_change(GtkWidget *w)
{
G_CONST_RETURN gchar *tmp;
struct file_pak *file_data;

tmp = gtk_entry_get_text(GTK_ENTRY(w));

/* update if valid type */
file_data = get_file_info((gpointer *) tmp, BY_LABEL);
if (file_data)
  {
  sysenv.file_type = file_data->id;
  update_file_pane();
  }
}

/********************************************/
/* primary event handler for file selection */
/********************************************/
void file_event_handler(GtkWidget *w, 
                        gpointer secondary_handler(gchar *, struct model_pak *))
{
gchar *name, *fullname;

/* get the current name */
/* NB: we strdup, so we don't have to make name a gconst pointer */
name = g_strdup(gtk_entry_get_text(GTK_ENTRY(file_name)));

/* attempt to use name to change the path */
fullname = g_build_filename(sysenv.cwd, name, NULL);
if (set_path(fullname))
  secondary_handler(name, NULL);

g_free(fullname);
g_free(name);
}

/********************************/
/* cartesian/fractional toggles */
/********************************/
void toggle_save_type(GtkWidget *w, struct model_pak *data)
{
data->fractional ^= 1;
}

/**********************/
/* selection handlers */
/**********************/
static void file_path_activate(GtkTreeView *treeview, GtkTreePath *treepath)
{
gchar *text, *path;
GtkTreeIter iter;
GtkTreeModel *treemodel;

treemodel = gtk_tree_view_get_model(treeview);
gtk_tree_model_get_iter(treemodel, &iter, treepath);

gtk_tree_model_get(treemodel, &iter, FILE_PATH, &text, -1);

path = g_build_path(DIR_SEP, sysenv.cwd, text, NULL);

set_path(path);

g_free(text);
g_free(path);

update_file_pane();
}

static void cb_file_name_changed(GtkTreeSelection *selection, gpointer data)
{
gchar *text;
GtkTreeIter iter;
GtkTreeModel *treemodel;

treemodel = file_name_treemodel();
if (!treemodel)
  return;

if (gtk_tree_selection_get_selected(selection, &treemodel, &iter))
  {
  gtk_tree_model_get(treemodel, &iter, FILE_NAME, &text, -1);
  gtk_entry_set_text(GTK_ENTRY(file_name), text);
  g_free(text);
  }
}

static void file_name_activate(GtkTreeView *treeview, GtkTreePath *treepath)
{
gchar *text;
GtkTreeIter iter;
GtkTreeModel *treemodel;

treemodel = gtk_tree_view_get_model(treeview);
gtk_tree_model_get_iter(treemodel, &iter, treepath);
gtk_tree_model_get(treemodel, &iter, FILE_PATH, &text, -1);
file_load(text, NULL);
g_free(text);
}

/*******************/
/* cleanup on exit */
/*******************/
void file_cleanup(void)
{
/* NB: since the tree store is independant of the model's geom list */
/* it must be completely removed (and then restored) with the dialog */
gtk_list_store_clear(file_path_ts);
gtk_list_store_clear(file_name_ts);
free_slist(dir_list);

/* TODO - free row data (ie FILE_PATH/FILE_NAME strings) */
file_path_tv = NULL;
file_path_ts = NULL;
file_name_tv = NULL;
file_name_ts = NULL;
}

/*************************************************/
/* customized load widget (with filetype filter) */
/*************************************************/
void 
file_dialog(gchar *title, 
            gchar *name, 
            gint type,
            gpointer secondary_handler(gchar *, struct model_pak *),
            gint filter)
{
gpointer dialog;
GSList *flist;
GList *elist;
GtkWidget *window, *swin, *hbox, *combo, *label;
GtkCellRenderer *renderer;
GtkTreeViewColumn *column;
GtkTreeSelection *select;
struct model_pak *data=NULL;
struct file_pak *fdata;

/* if save type, check we have a loaded model */
if (secondary_handler == (gpointer) file_save_as)
  {
  data = sysenv.active_model;
  if (!data)
    return;
  }

/* get a dialog if possible */
dialog = dialog_request(FILE_SELECT, "File dialog", NULL, NULL, NULL);
if (!dialog)
  return;
window = dialog_window(dialog);

/* make and set up the dialog window */
gtk_window_set_title(GTK_WINDOW(window), title);
gtk_window_set_default_size(GTK_WINDOW(window), 600, 400);
gtk_window_set_position(GTK_WINDOW(window), GTK_WIN_POS_CENTER);
gtk_container_set_border_width(GTK_CONTAINER(GTK_BOX(GTK_DIALOG(window)->vbox)),10);

/* TOP ROW - cwd printed */
hbox = gtk_hbox_new(FALSE,0);
gtk_box_pack_start(GTK_BOX(GTK_DIALOG(window)->vbox),hbox,FALSE,FALSE,0);
label = gtk_label_new("Current path: ");
gtk_box_pack_start(GTK_BOX(hbox), label, FALSE, FALSE, 0);

curr_path = gtk_label_new("dummy");
gtk_box_pack_start(GTK_BOX(hbox), curr_path, FALSE, FALSE, 0);

/* option to save as cartesian or fractional coords */
if (secondary_handler == (gpointer) file_save_as)
  {
  g_assert(data != NULL);
  if (data->periodic)
    {
    hbox = gtk_hbox_new(FALSE,0);
    gtk_box_pack_start(GTK_BOX(GTK_DIALOG(window)->vbox),
                                         hbox,FALSE,FALSE,0);

    new_check_button("Save as cartesian coordinates", toggle_save_type, data,
                                                    !data->fractional, hbox);
    }
  }

/* SECOND ROW - sub directory & file listings */
hbox = gtk_hbox_new(FALSE, PANEL_SPACING);
gtk_box_pack_start(GTK_BOX(GTK_DIALOG(window)->vbox), hbox, TRUE, TRUE, 0);
gtk_container_set_border_width(GTK_CONTAINER(hbox), PANEL_SPACING);

/* scrolled model pane */
swin = gtk_scrolled_window_new(NULL, NULL);
gtk_scrolled_window_set_policy(GTK_SCROLLED_WINDOW (swin),
                               GTK_POLICY_AUTOMATIC, GTK_POLICY_AUTOMATIC);
gtk_box_pack_start(GTK_BOX(hbox), swin, TRUE, TRUE, 0);

/* path treestore/treeview */
file_path_ts = gtk_list_store_new(FILE_PATH_NCOLS, G_TYPE_STRING);
file_path_tv = gtk_tree_view_new_with_model(GTK_TREE_MODEL(file_path_ts));
gtk_scrolled_window_add_with_viewport(GTK_SCROLLED_WINDOW(swin), file_path_tv);
renderer = gtk_cell_renderer_text_new();
column = gtk_tree_view_column_new_with_attributes("Tree", renderer, "text", FILE_PATH, NULL);
gtk_tree_view_append_column(GTK_TREE_VIEW(file_path_tv), column);

/* NB: use this method instead of selections to handle events */
/* as selections will core dump if the handler clears the store */
g_signal_connect(file_path_tv, "row_activated", G_CALLBACK(file_path_activate), NULL);

/* scrolled model pane */
swin = gtk_scrolled_window_new(NULL, NULL);
gtk_scrolled_window_set_policy(GTK_SCROLLED_WINDOW (swin),
                               GTK_POLICY_AUTOMATIC, GTK_POLICY_AUTOMATIC);
gtk_box_pack_start(GTK_BOX(hbox), swin, TRUE, TRUE, 0);

/* filename treestore/treeview */
file_name_ts = gtk_list_store_new(FILE_NAME_NCOLS, G_TYPE_STRING);
file_name_tv = gtk_tree_view_new_with_model(GTK_TREE_MODEL(file_name_ts));
gtk_scrolled_window_add_with_viewport(GTK_SCROLLED_WINDOW(swin), file_name_tv);
renderer = gtk_cell_renderer_text_new();
column = gtk_tree_view_column_new_with_attributes("Files", renderer, "text", FILE_NAME, NULL);
gtk_tree_view_append_column(GTK_TREE_VIEW(file_name_tv), column);
/* setup the selection handler */
select = gtk_tree_view_get_selection(GTK_TREE_VIEW(file_name_tv));
gtk_tree_selection_set_mode(select, GTK_SELECTION_BROWSE);
/* single click - update entry */
g_signal_connect(G_OBJECT(select), "changed",
                 G_CALLBACK(cb_file_name_changed),
                 NULL);

/* double click - automatic load */
if (secondary_handler == (gpointer) file_load)
  g_signal_connect(file_name_tv, "row_activated", G_CALLBACK(file_name_activate), NULL);


/* THIRD ROW - filename currently selected & file extension filter */
hbox = gtk_hbox_new(FALSE, PANEL_SPACING);
gtk_box_pack_start(GTK_BOX(GTK_DIALOG(window)->vbox), hbox, FALSE, FALSE, 0);

/* filename */
file_name = gtk_entry_new();
gtk_box_pack_start(GTK_BOX(hbox), file_name, TRUE, TRUE, 0);
if (name)
  gtk_entry_set_text(GTK_ENTRY(file_name), name);
gtk_entry_set_editable(GTK_ENTRY(file_name), TRUE);

/* hook a <CR> event to the load action */
g_signal_connect(GTK_OBJECT(file_name), "activate", 
                 GTK_SIGNAL_FUNC(file_event_handler), secondary_handler);

/* build the recognized extension list */
elist = NULL;
flist = sysenv.file_list;
while(flist != NULL)
  {
  fdata = flist->data;
/* include in menu? */
  if (fdata->menu)
    elist = g_list_append(elist, fdata->label);
  flist = g_slist_next(flist);
  }

/* combo box for file filter */
combo = gtk_combo_new();
gtk_entry_set_editable(GTK_ENTRY(GTK_COMBO(combo)->entry), FALSE);
gtk_combo_set_popdown_strings(GTK_COMBO(combo), elist);

/* set the currently selected type (BEFORE the changed event is connected) */
fdata = get_file_info(GINT_TO_POINTER(filter), BY_FILE_ID);
if (fdata)
  gtk_entry_set_text(GTK_ENTRY(GTK_COMBO(combo)->entry), fdata->label);

gtk_box_pack_start(GTK_BOX (hbox), combo, TRUE, TRUE, 0);
g_signal_connect(GTK_OBJECT(GTK_COMBO(combo)->entry), "changed", 
                 GTK_SIGNAL_FUNC(type_change), NULL);
gtk_widget_set_size_request(combo, 6*sysenv.gtk_fontsize, -1);

/* terminating buttons */
switch (type)
  {
  case FILE_LOAD:
    gui_stock_button(GTK_STOCK_OPEN, file_event_handler, secondary_handler, 
                       GTK_DIALOG(window)->action_area);
    break;

  case FILE_SAVE:
    gui_stock_button(GTK_STOCK_SAVE, file_event_handler, secondary_handler, 
                       GTK_DIALOG(window)->action_area);
    break;
  }

gui_stock_button(GTK_STOCK_CANCEL, dialog_destroy, dialog,
                   GTK_DIALOG(window)->action_area);

/* all done */
gtk_widget_show_all(window);
sysenv.file_type = filter;
update_file_pane();
}

/*************************************************/
/* callback to update a widget background colour */
/*************************************************/
void set_colour_widget(GtkWidget *cs, GtkWidget *w)
{
GdkColor colour;

if (GTK_IS_WIDGET(cs))
  gtk_color_selection_get_current_color(GTK_COLOR_SELECTION(cs), &colour);
if (GTK_IS_WIDGET(w))
  gtk_widget_modify_bg(w, GTK_STATE_NORMAL, &colour);
}

/***************************************/
/* callback to update the colour value */
/***************************************/
void set_colour_value(GtkWidget *cs, gdouble *x)
{
GdkColor colour;

gtk_color_selection_get_current_color(GTK_COLOR_SELECTION(cs), &colour);

x[0] = colour.red;
x[1] = colour.green;
x[2] = colour.blue;
VEC3MUL(x, 1.0/65535.0);
}

/**********************************/
/* create a colour editing dialog */
/**********************************/
void dialog_colour_new(GtkWidget *w, gdouble *colour)
{
GtkWidget *csd;

/* create colour selection dialog */
csd = gtk_color_selection_dialog_new("Edit colour");

/* set the initial colour */
gtk_color_selection_set_color(GTK_COLOR_SELECTION(GTK_COLOR_SELECTION_DIALOG(csd)->colorsel), colour);

/* value colour update */
g_signal_connect(GTK_OBJECT(GTK_COLOR_SELECTION_DIALOG(csd)->colorsel), "color_changed",
                 GTK_SIGNAL_FUNC(set_colour_value), (gpointer) colour);
/* widget colour update */
if (w)
  g_signal_connect(GTK_OBJECT(GTK_COLOR_SELECTION_DIALOG(csd)->colorsel), "color_changed",
                   GTK_SIGNAL_FUNC(set_colour_widget), (gpointer) w);

g_signal_connect(GTK_OBJECT(GTK_COLOR_SELECTION_DIALOG(csd)->ok_button), "clicked",
                 GTK_SIGNAL_FUNC(destroy_widget), (gpointer) csd);

gtk_widget_hide(GTK_COLOR_SELECTION_DIALOG(csd)->cancel_button);

gtk_widget_show(csd);

/*
g_signal_connect(GTK_OBJECT(w), "destroy", GTK_SIGNAL_FUNC(destroy_widget), csd);
*/
}

