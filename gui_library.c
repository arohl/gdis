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
#include <time.h>

#include "gdis.h"
#include "model.h"
#include "library.h"
#include "interface.h"
#include "gui_shorts.h"

extern struct sysenv_pak sysenv;

gpointer active_folder = NULL;
gpointer active_entry = NULL;

/****************************/
/* fill out the folder list */
/****************************/
void gui_folder_populate(gpointer key, gpointer value, gpointer data)
{
gchar *name = key;
GtkTreeIter iter;
GtkListStore *list_store = data;

gtk_list_store_append(list_store, &iter);
gtk_list_store_set(list_store, &iter, 0, name, -1);
}

/****************************************************/
/* fill out the entry list (based on active folder) */
/****************************************************/
void gui_library_entry_populate(GtkListStore *list_store)
{
GSList *list;
GtkTreeIter iter;
struct folder_pak *folder;
struct entry_pak *entry;

gtk_list_store_clear(list_store);

if (active_folder)
  folder = g_hash_table_lookup(sysenv.library, active_folder);
else
  folder = g_hash_table_lookup(sysenv.library, "default");

if (folder)
  {
  for (list=folder->list ; list ; list=g_slist_next(list))
    {
    entry = list->data;

    gtk_list_store_append(list_store, &iter);
    gtk_list_store_set(list_store, &iter, 0, entry->name, 1, entry, -1);
    }
  }
}

/****************************/
/* folder selection handler */
/****************************/
void gui_library_folder_selected(GtkTreeSelection *selection, gpointer list)
{
GtkTreeIter iter;
GtkTreeModel *treemodel;

/* record selection as the active folder */
if (gtk_tree_selection_get_selected(selection, &treemodel, &iter))
  {
  gtk_tree_model_get(treemodel, &iter, 0, &active_folder, -1);
/*
printf("select: [%s]\n", (gchar *) active_folder);
*/
  }

/* redraw the entry list */
gui_library_entry_populate(list);
}

/***************************/
/* entry selection handler */
/***************************/
void gui_library_entry_selected(GtkTreeSelection *selection, gpointer data)
{
gint n;
gchar *text;
GtkTreeIter iter;
GtkTreeModel *treemodel;
GtkTextBuffer *buffer;
struct entry_pak *entry;

/* record selection as the active entry */
if (gtk_tree_selection_get_selected(selection, &treemodel, &iter))
  {
  gtk_tree_model_get(treemodel, &iter, 1, &entry, -1);
  active_entry = entry;
  buffer = gtk_text_view_get_buffer(GTK_TEXT_VIEW(data));
  text = (entry->info)->str;
  if (text)
    {
    n = strlen(text);
    gtk_text_buffer_set_text(buffer, text, n);
    }
  }
}

/*************************************/
/* import the selected library entry */
/*************************************/
void gui_library_entry_import(void)
{
struct entry_pak *entry;
struct model_pak *model;

if (active_entry)
  {
  entry = active_entry;

/*
printf("importing: %s [%p]\n", entry->name, entry->offset);
*/

  model = model_new();
  if (library_entry_get(entry->offset, model))
    {
    gui_text_show(ERROR, "Failed to load library file.\n");
    model_delete(model);
    return;
    }

/* set the name */
  g_free(model->basename);
  model->basename = g_strdup(entry->name);

/* update/redraw */
  tree_model_add(model);
  tree_select_model(model);
  redraw_canvas(SINGLE);
  }
}

/*****************************/
/* display the model library */
/*****************************/
void gui_library_window(GtkWidget *box)
{
GtkCellRenderer *r;
GtkTreeViewColumn *c;
GtkListStore *list1, *list2;
GtkTreeSelection *select;
GtkWidget *frame, *hbox, *vbox, *swin, *tree1, *tree2, *view;

g_assert(box != NULL);

/* frame with split pane (folders : entries) */
frame = gtk_frame_new(NULL);
gtk_box_pack_start(GTK_BOX(box),frame,TRUE,TRUE,0); 
hbox = gtk_hbox_new(FALSE, PANEL_SPACING);
gtk_container_add(GTK_CONTAINER(frame), hbox);
gtk_container_set_border_width(GTK_CONTAINER(GTK_BOX(hbox)), PANEL_SPACING);

/* list of library folders */
swin = gtk_scrolled_window_new(NULL, NULL);
gtk_scrolled_window_set_policy(GTK_SCROLLED_WINDOW(swin),
                               GTK_POLICY_AUTOMATIC, GTK_POLICY_AUTOMATIC);
gtk_box_pack_start(GTK_BOX(hbox), swin, TRUE, TRUE, 0);

/* list */
list1 = gtk_list_store_new(1, G_TYPE_STRING);
tree1 = gtk_tree_view_new_with_model(GTK_TREE_MODEL(list1));
gtk_scrolled_window_add_with_viewport(GTK_SCROLLED_WINDOW(swin), tree1);
r = gtk_cell_renderer_text_new();
c = gtk_tree_view_column_new_with_attributes(" ", r, "text", 0, NULL);
gtk_tree_view_append_column(GTK_TREE_VIEW(tree1), c);
gtk_tree_view_set_headers_visible(GTK_TREE_VIEW(tree1), FALSE);

g_hash_table_foreach(sysenv.library, gui_folder_populate, list1);

/* list of model entries */
swin = gtk_scrolled_window_new(NULL, NULL);
gtk_scrolled_window_set_policy(GTK_SCROLLED_WINDOW(swin),
                               GTK_POLICY_AUTOMATIC, GTK_POLICY_AUTOMATIC);
gtk_box_pack_start(GTK_BOX(hbox), swin, TRUE, TRUE, 0);

/* list */
list2 = gtk_list_store_new(2, G_TYPE_STRING, G_TYPE_POINTER);
tree2 = gtk_tree_view_new_with_model(GTK_TREE_MODEL(list2));
gtk_scrolled_window_add_with_viewport(GTK_SCROLLED_WINDOW(swin), tree2);
r = gtk_cell_renderer_text_new();
c = gtk_tree_view_column_new_with_attributes(" ", r, "text", 0, NULL);
gtk_tree_view_append_column(GTK_TREE_VIEW(tree2), c);
gtk_tree_view_set_headers_visible(GTK_TREE_VIEW(tree2), FALSE);

gui_library_entry_populate(list2);

/* library entry description text */
vbox = gtk_vbox_new(FALSE, PANEL_SPACING);
gtk_box_pack_start(GTK_BOX(box),vbox,FALSE,FALSE,PANEL_SPACING); 
view = gtk_text_view_new();
gtk_text_view_set_editable(GTK_TEXT_VIEW(view), FALSE);
gtk_box_pack_start(GTK_BOX(vbox), view, TRUE, TRUE, 0);

/* actions */
/*
gui_button_x(" Import model", gui_library_entry_import, NULL, vbox);
*/

/* folder selection handler */
select = gtk_tree_view_get_selection(GTK_TREE_VIEW(tree1));
gtk_tree_selection_set_mode(select, GTK_SELECTION_SINGLE);
g_signal_connect(G_OBJECT(select), "changed",
                 G_CALLBACK(gui_library_folder_selected),
                 list2);

/* entry selection handler */
select = gtk_tree_view_get_selection(GTK_TREE_VIEW(tree2));
gtk_tree_selection_set_mode(select, GTK_SELECTION_SINGLE);
g_signal_connect(G_OBJECT(select), "changed",
                 G_CALLBACK(gui_library_entry_selected),
                 view);

g_signal_connect(G_OBJECT(tree2), "row-activated",
                 G_CALLBACK(gui_library_entry_import),
                 NULL);
}

