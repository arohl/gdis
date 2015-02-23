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

#ifndef _WIN32

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "gdis.h"
#include "host.h"
#include "dialog.h"
#include "interface.h"
#include "gui_shorts.h"

extern struct sysenv_pak sysenv;

GtkListStore *host_ls=NULL;
GtkWidget *host_tv=NULL;

/* TODO - list of queued/running/etc jobs */
/* TODO - intended to replace task manager ... since those are just jobs running on localhost */
void gui_job_page(GtkWidget *box)
{
}

/**********************************************/
/* host program list GUI population primitive */
/**********************************************/
void gui_host_service_add(gchar *name, gpointer service, gpointer data)
{
gint flags;
gchar *available, *invoke;
GtkTreeIter iter;

flags = host_service_flags(service);

invoke = "";
switch (flags)
  {
  case SERVICE_BACKGROUND:
    invoke = "background job"; 
    break;
  case SERVICE_MPI:
    invoke = "mpi job"; 
    break;
  case SERVICE_QSUB:
    invoke = "queued job"; 
    break;
  case SERVICE_QSUB_MPI:
    invoke = "queued mpi job"; 
    break;
  }

gtk_list_store_append(host_ls, &iter);

/* test if the service configuration is available */
if (host_service_available(service))
  available = "yes";
else
  available = "no";

gtk_list_store_set(host_ls, &iter, 0, name, 1, invoke, 2, available, -1);
}

/*****************************************/
/* write out host info to host info list */
/*****************************************/
void gui_host_update(gpointer host, gpointer dialog)
{
if (!host_ls)
  return;

gtk_list_store_clear(host_ls);

host_service_foreach(host, gui_host_service_add, dialog);
}

/****************************/
/* connect to host callback */
/****************************/
void gui_host_connect(GtkWidget *w, gpointer dialog)
{
gchar *text;
gpointer host;
GtkEntry *user_entry, *host_entry;

user_entry = dialog_child_get(dialog, "user");
host_entry = dialog_child_get(dialog, "host");

text = g_strdup_printf("%s@%s", gtk_entry_get_text(user_entry), gtk_entry_get_text(host_entry));
host = host_new(text);
g_free(text);

if (host_connect(host))
  printf("host_connect(): success!\n");
else
  printf("host_connect(): failed!\n");

gui_host_update(host, dialog);
}

/*********************************/
/* service method cycle callback */
/*********************************/
void gui_job_method_change(GtkWidget *w, gpointer dialog)
{
gint flags;
gchar *name, *invoke;
gpointer host, service;
GtkTreeModel *treemodel;
GtkTreeSelection *selection;
GtkTreeIter iter;

/* TODO - active host, rather than 1st on list */
if (sysenv.host_list)
  host = sysenv.host_list->data;
else
  {
  printf("No active host connections.\n");
  return;
  }

treemodel = gtk_tree_view_get_model(GTK_TREE_VIEW(host_tv));
selection = gtk_tree_view_get_selection(GTK_TREE_VIEW(host_tv));

if (gtk_tree_selection_get_selected(selection, &treemodel, &iter))
  {
  gtk_tree_model_get (treemodel, &iter, 0, &name, -1);

  host_service_cycle(host, name);

  service = host_service_get(host, name);

flags = host_service_flags(service);
invoke = "";
switch (flags)
  {
  case SERVICE_BACKGROUND:
    invoke = "background job"; 
    break;
  case SERVICE_MPI:
    invoke = "mpi job"; 
    break;
  case SERVICE_QSUB:
    invoke = "queued job"; 
    break;
  case SERVICE_QSUB_MPI:
    invoke = "queued mpi job"; 
    break;
  }

  if (host_service_available(service))
    gtk_list_store_set(host_ls, &iter, 1, invoke, 2, "yes", -1);
  else
    gtk_list_store_set(host_ls, &iter, 1, invoke, 2, "no", -1);
  }
}

/************************/
/* host setup/info page */
/************************/
void gui_host_page(GtkWidget *box, gpointer dialog)
{
gchar *titles[] = {" Service ", " Method ", " Available "};
GtkWidget *swin, *table, *label, *entry, *button;
GtkCellRenderer *renderer;
GtkTreeViewColumn *column;

table = gtk_table_new(3, 3, FALSE);
gtk_box_pack_start(GTK_BOX(box), table, FALSE, FALSE, 0);

label = gtk_label_new("Hostname");
gtk_table_attach_defaults(GTK_TABLE(table),label,0,1,0,1);
label = gtk_label_new("Username");
gtk_table_attach_defaults(GTK_TABLE(table),label,0,1,1,2);
/*
label = gtk_label_new("Password");
gtk_table_attach_defaults(GTK_TABLE(table),label,0,1,2,3);
*/

entry = gtk_entry_new();
gtk_table_attach_defaults(GTK_TABLE(table),entry,1,2,0,1);
dialog_child_set(dialog, "host", entry);

entry = gtk_entry_new();
gtk_table_attach_defaults(GTK_TABLE(table),entry,1,2,1,2);
dialog_child_set(dialog, "user", entry);

/* TODO - leave this out until the security/privacy implications have been dealt with */
/*
entry = gtk_entry_new();
gtk_table_attach_defaults(GTK_TABLE(table),entry,1,2,2,3);
*/

button = gtk_button_new_with_label(" Connect ");
gtk_table_attach_defaults(GTK_TABLE(table),button,2,3,0,1);
g_signal_connect(GTK_OBJECT(button), "clicked", GTK_SIGNAL_FUNC(gui_host_connect), dialog);

/* current connected host properties */
swin = gtk_scrolled_window_new(NULL, NULL);
gtk_scrolled_window_set_policy(GTK_SCROLLED_WINDOW (swin),
                               GTK_POLICY_AUTOMATIC, GTK_POLICY_AUTOMATIC);
gtk_box_pack_start(GTK_BOX(box), swin, TRUE, TRUE, 0);

host_ls = gtk_list_store_new(3, G_TYPE_STRING, G_TYPE_STRING, G_TYPE_STRING);
host_tv = gtk_tree_view_new_with_model(GTK_TREE_MODEL(host_ls));
gtk_scrolled_window_add_with_viewport(GTK_SCROLLED_WINDOW(swin), host_tv);

renderer = gtk_cell_renderer_text_new();
column = gtk_tree_view_column_new_with_attributes(titles[0], renderer, "text", 0, NULL);
gtk_tree_view_append_column(GTK_TREE_VIEW(host_tv), column);
gtk_tree_view_column_set_expand(column, FALSE);

renderer = gtk_cell_renderer_text_new();
column = gtk_tree_view_column_new_with_attributes(titles[1], renderer, "text", 1, NULL);
gtk_tree_view_append_column(GTK_TREE_VIEW(host_tv), column);
gtk_tree_view_column_set_expand(column, TRUE);
gtk_tree_view_column_set_clickable(column, TRUE);
g_signal_connect(GTK_OBJECT(column), "clicked", GTK_SIGNAL_FUNC(gui_job_method_change), dialog);

renderer = gtk_cell_renderer_text_new();
column = gtk_tree_view_column_new_with_attributes(titles[2], renderer, "text", 2, NULL);
gtk_tree_view_append_column(GTK_TREE_VIEW(host_tv), column);
gtk_tree_view_column_set_expand(column, FALSE);
}

/*****************************/
/* the gdis job setup dialog */
/*****************************/
void gui_job_dialog(void)
{
gpointer dialog;
GtkWidget *window, *notebook, *page, *label;

dialog = dialog_request(777, "Remote job manager", NULL, NULL, NULL);
if (!dialog)
  return;
window = dialog_window(dialog);

gtk_widget_set_size_request(window, 400, 400);

notebook = gtk_notebook_new();
gtk_notebook_set_tab_pos(GTK_NOTEBOOK(notebook), GTK_POS_TOP);
gtk_container_add(GTK_CONTAINER(GTK_DIALOG(window)->vbox), notebook);
gtk_notebook_set_show_border(GTK_NOTEBOOK(notebook), FALSE);

/* run type page */
page = gtk_vbox_new(FALSE, 0);
label = gtk_label_new(" Hosts ");
gtk_notebook_append_page(GTK_NOTEBOOK(notebook), page, label);
gui_host_page(page, dialog);

/* run type page */
page = gtk_vbox_new(FALSE, 0);
label = gtk_label_new(" Jobs ");
gtk_notebook_append_page(GTK_NOTEBOOK(notebook), page, label);
gui_job_page(page);

gtk_widget_show_all(window);

/* FIXME - active host */
if (sysenv.host_list)
  gui_host_update(sysenv.host_list->data, dialog);
}
#endif
