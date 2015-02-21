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
#include <stdlib.h>
#include <string.h>

#include "gdis.h"
#include "host.h"
#include "job.h"
#include "task.h"
#include "file.h"
#include "parse.h"
#include "dialog.h"
#include "interface.h"
#include "gui_shorts.h"
#include "grid.h"
#include "grisu_client.h"

extern struct sysenv_pak sysenv;

GtkListStore *grisu_ls=NULL;
GtkWidget *grisu_tv=NULL;
GtkListStore *grid_ls=NULL;
GtkWidget *grid_tv=NULL;
GtkListStore *job_ls=NULL;
GtkWidget *job_tv=NULL;
GtkWidget *grid_id_label, *grid_id_entry, *grid_password_entry, *grid_status_entry;
GtkWidget *grid_vo_combo;
GtkWidget *grid_configure_box;
GtkWidget *grid_service_entry;
GtkWidget *grid_myproxy_server_entry;
GtkWidget *grid_myproxy_port_entry;
GSList *grid_vo_list=NULL;
GtkListStore *grid_vo_store=NULL;
gint grid_method = GRID_MYPROXY;

/* CURRENT - used to enforce 1 background gui task at a time */
GMutex *gui_grid_mutex=NULL;
GtkWidget *grid_status_bar;


/* CURRENT - data struct for background tasks */

struct grid_data_pak
{
gchar *application;
GSList *queues;
GSList *jobname;
GSList *status;
GSList *files;
};

gpointer grid_data_new(void)
{
struct grid_data_pak *data;

data = g_malloc(sizeof(struct grid_data_pak));

data->application=NULL;
data->queues=NULL;
data->jobname = NULL;
data->status = NULL;
data->files = NULL;

return(data);
}

void grid_data_free(struct grid_data_pak *data)
{
g_assert(data != NULL);

g_free(data->application);

free_slist(data->queues);
free_slist(data->jobname);
free_slist(data->status);
free_slist(data->files);

g_free(data);
}

/***************************************************/
/* primtives for activity with unknown work length */
/***************************************************/
gint grid_task_start(void)
{
if (g_mutex_trylock(gui_grid_mutex))
  {
  if (GTK_IS_ENTRY(grid_status_bar))
    gtk_entry_set_text(GTK_ENTRY(grid_status_bar), "Operation in progress, please wait ...");
  return(TRUE);
  }
else
  printf("I'm busy, please wait.\n");

return(FALSE);
}

void grid_task_stop(void)
{
if (GTK_IS_ENTRY(grid_status_bar))
  gtk_entry_set_text(GTK_ENTRY(grid_status_bar), "");
g_mutex_unlock(gui_grid_mutex);
}

void grid_task_feedback(const gchar *text)
{
if (GTK_IS_ENTRY(grid_status_bar))
  gtk_entry_set_text(GTK_ENTRY(grid_status_bar), text);
}

/*********************************/
/* GUI vo entry changed callback */
/*********************************/
void gui_vo_changed(GtkWidget *w, gpointer data)
{
#ifdef WITH_GRISU
const gchar *vo;

vo = gtk_entry_get_text(GTK_ENTRY(GTK_BIN(grid_vo_combo)->child));

if (g_strncasecmp(vo, "no", 2) == 0)
  grisu_vo_set("");
else
  grisu_vo_set(vo);
#endif
}

/*****************************/
/* refresh VO combo pulldown */
/*****************************/
void gui_vo_refresh(void)
{
#if GTK_MAJOR_VERSION >= 2 && GTK_MINOR_VERSION >= 4
GSList *list;
GtkTreeIter iter;

g_assert(grid_vo_store != NULL);

gtk_list_store_clear(grid_vo_store);

for (list=grid_vo_list ; list ; list=g_slist_next(list))
  {
  gtk_list_store_append(grid_vo_store, &iter);
  gtk_list_store_set(grid_vo_store, &iter, 0, list->data, -1);
  }

/* TODO - print this text as the default on startup */
gtk_list_store_prepend(grid_vo_store, &iter);
gtk_list_store_set(grid_vo_store, &iter, 0, "No VO", -1);

#else
/* use deprec combo */
#endif
}

/*************************************/
/* change grid authentication method */
/*************************************/
void gui_grid_method(GtkWidget *w, gpointer dummy)
{
#ifdef WITH_GRISU
gchar *subject;
const gchar *text;

if (w)
  {
/* change event only */
  text = gtk_label_get_text(GTK_LABEL(grid_id_label));
  if (g_strncasecmp(text, "grid", 4) == 0)
    grid_method = GRID_MYPROXY;
  else
    grid_method = GRID_LOCALPROXY;
  }
else
  {
/* initial setup */
  switch (grid_method)
    {
    case GRID_MYPROXY:
/* want the user to enter this, so scrub the default */
      grisu_keypair_set("","");
      break;

    case GRID_LOCALPROXY:
/* auto generate these */
      grisu_keypair_set(grid_random_alpha(10),grid_random_alphanum(12));
      break;
    }
  }

/* common GUI update */
switch (grid_method)
  {
  case GRID_MYPROXY:
    gtk_label_set_text(GTK_LABEL(grid_id_label), "Username");
    text = grisu_username_get();
    if (text)
      gtk_entry_set_text(GTK_ENTRY(grid_id_entry), text);
    else
      gtk_entry_set_text(GTK_ENTRY(grid_id_entry), "");
    break;

  case GRID_LOCALPROXY:
    gtk_label_set_text(GTK_LABEL(grid_id_label), "Grid Identity");
    subject = grid_get_DN("/home/sean/.globus/usercert.pem");
    gtk_entry_set_text(GTK_ENTRY(grid_id_entry), subject);
    g_free(subject);
    break;
  }
#endif
}

/*************************************************/
/* toggle the display of the configuration panel */
/*************************************************/
void gui_grid_configure(gpointer dummy)
{
static gint visible=TRUE;

if (visible)
  {
  gtk_widget_hide(grid_configure_box);
  visible = FALSE;
  }
else
  {
  gtk_widget_show(grid_configure_box);
  visible = TRUE;

#ifdef WITH_GRISU
/* update boxes with data */
  gtk_entry_set_text(GTK_ENTRY(grid_service_entry), grisu_ws_get());
  gtk_entry_set_text(GTK_ENTRY(grid_myproxy_server_entry), grisu_myproxy_get());
  gtk_entry_set_text(GTK_ENTRY(grid_myproxy_port_entry), grisu_myproxyport_get());
#endif
  }
}

/********************************/
/* grid connect background task */
/********************************/
void grid_connect_start(gpointer dummy)
{
#ifdef WITH_GRISU
grid_vo_list = grisu_fqans_get();

/* FIXME - is a NULL vo list proof of failed auth? */
if (grid_vo_list)
  grid_auth_set(TRUE);
else
  grid_auth_set(FALSE);
#endif
}

/********************************/
/* grid connect GUI update task */
/********************************/
void grid_connect_stop(gpointer dummy)
{
gui_vo_refresh();

gtk_entry_set_text(GTK_ENTRY(grid_status_entry), grid_get_status());

grid_task_stop();
}

/***************************************/
/* apply current authentication method */
/***************************************/
void gui_grid_connect(GtkWidget *w, gpointer dummy)
{
#ifdef WITH_GRISU
/* CURRENT - mutex lock - 1 grid op at a time */
if (!grid_task_start())
  return;

switch (grid_method)
  {
  case GRID_MYPROXY:
    grisu_keypair_set(gtk_entry_get_text(GTK_ENTRY(grid_id_entry)),
                      gtk_entry_get_text(GTK_ENTRY(grid_password_entry)));
    break; 

  case GRID_LOCALPROXY:
    grisu_keypair_set(NULL, gtk_entry_get_text(GTK_ENTRY(grid_password_entry)));
/* upload proxy */
    grid_connect(grid_method);
    break;
  }

/* get VO's for display (also serves as auth check) */
grisu_auth_set(TRUE);

if (grid_vo_list)
  free_slist(grid_vo_list);

task_new("Authenticate", grid_connect_start, NULL, grid_connect_stop, NULL, NULL);
#endif
}

/**********************/
/* grid configuration */
/**********************/
void gui_grid_setup_page(GtkWidget *box, gpointer dialog)
{
GtkWidget *table, *label, *button;

table = gtk_table_new(3, 4, FALSE);
gtk_box_pack_start(GTK_BOX(box), table, FALSE, FALSE, 0);

/* reflect the authentication mechanism (grid proxy, myproxy, shib?) */
grid_id_label = gtk_label_new(" xxxxxxxxxxxxxxxx ");
gtk_table_attach(GTK_TABLE(table),grid_id_label,0,1,0,1, GTK_SHRINK, GTK_SHRINK, 0, 0);

label = gtk_label_new(" Password ");
gtk_table_attach(GTK_TABLE(table),label,0,1,1,2, GTK_SHRINK, GTK_SHRINK, 0, 0);

label = gtk_label_new(" Status ");
gtk_table_attach(GTK_TABLE(table),label,0,1,2,3, GTK_SHRINK, GTK_SHRINK, 0, 0);

grid_id_entry = gtk_entry_new();
gtk_table_attach_defaults(GTK_TABLE(table),grid_id_entry,1,2,0,1);

grid_password_entry = gtk_entry_new();
gtk_entry_set_visibility(GTK_ENTRY(grid_password_entry), FALSE);
gtk_table_attach_defaults(GTK_TABLE(table),grid_password_entry,1,2,1,2);
g_signal_connect(GTK_OBJECT(grid_password_entry), "activate",
                 GTK_SIGNAL_FUNC(gui_grid_connect), NULL);

grid_status_entry = gtk_entry_new();
gtk_table_attach_defaults(GTK_TABLE(table),grid_status_entry,1,2,2,3);
gtk_entry_set_text(GTK_ENTRY(grid_status_entry), grid_get_status());

//button = gui_button("Authenticate", gpointer cb, gpointer arg, NULL, 0);

/* TODO - this is inteded for specifying certificate location (assuming not in default location) */
/* TODO - also for specifying what myproxy server / port / etc to use */
button = gui_button(" Configure ", gui_grid_configure, NULL, NULL, 0);
gtk_table_attach(GTK_TABLE(table),button,2,3,0,1, GTK_FILL, GTK_FILL, 0, 0);

/* this is for uploading a credential to the myproxy server */
button = gui_button(" Connect ", gui_grid_connect, NULL, NULL, 0);
gtk_table_attach(GTK_TABLE(table),button,2,3,1,2, GTK_FILL, GTK_FILL, 0, 0);

#if GTK_MAJOR_VERSION >= 2 && GTK_MINOR_VERSION >= 4

/* create store and fill it */
grid_vo_store = gtk_list_store_new(1, G_TYPE_STRING);
gui_vo_refresh();

/* create combo box and set start value */
grid_vo_combo = gtk_combo_box_entry_new_with_model(GTK_TREE_MODEL(grid_vo_store), 0);
gtk_entry_set_text(GTK_ENTRY(GTK_BIN(grid_vo_combo)->child), "No VO");

gtk_table_attach(GTK_TABLE(table),grid_vo_combo,2,3,2,3, GTK_FILL, GTK_FILL, 0, 0);

g_signal_connect(GTK_OBJECT(GTK_BIN(grid_vo_combo)->child), "changed",
                 GTK_SIGNAL_FUNC(gui_vo_changed), NULL);

#else
// < GTK 2.2
// what to do here?
#endif

/* init for default method */
gui_grid_method(NULL, NULL);
}

/********************/
/* back ground task */
/********************/
void grid_refresh_all_start(gpointer data)
{
#ifdef WITH_GRISU
const gchar *status;
GSList *list;
struct grid_data_pak *grid_data = data;

g_assert(grid_data != NULL);

/* query grisu for all job names */
grid_data->jobname = grisu_job_names();
grid_data->status = NULL;

for (list=grid_data->jobname ; list ; list=g_slist_next(list))
  {
  status = grisu_job_status(list->data);
  grid_data->status = g_slist_prepend(grid_data->status, g_strdup(status));
  }

grid_data->status = g_slist_reverse(grid_data->status);
#endif
}

/*******************/
/* GUI update task */
/*******************/
void grid_refresh_all_stop(gpointer data)
{
#ifdef WITH_GRISU
const gchar *name, *status;
GSList *list1, *list2;
GtkTreeIter iter;
struct grid_data_pak *grid_data = data;

g_assert(grid_data != NULL);

gtk_list_store_clear(job_ls);

list2 = grid_data->status;
for (list1=grid_data->jobname ; list1 ; list1=g_slist_next(list1))
  {
  name = list1->data;
  if (list2)
    status = list2->data;
  else
    status = "undefined";  // shouldn't happen

  gtk_list_store_append(job_ls, &iter);
  gtk_list_store_set(job_ls, &iter, 0, name, 1, status, -1);

  list2 = g_slist_next(list2);
  }

grid_data_free(grid_data);

grid_task_stop();
#endif
}

/**********************/
/* refresh everything */
/**********************/
void gui_job_refresh_all(void)
{
#ifdef WITH_GRISU
struct grid_data_pak *data;

if (g_strncasecmp(grid_get_status(), "not", 3) == 0)
  {
  grid_task_feedback("Not authenticated.");
  return;
  }

/* CURRENT - mutex lock - 1 grid op at a time */
if (!grid_task_start())
  return;

data = grid_data_new();

task_new("Job list", grid_refresh_all_start, data, grid_refresh_all_stop, data, NULL);
#endif
}

/**************************************/
/* Refresh the status of selected job */
/**************************************/
void gui_job_refresh_selected(void)
{
#ifdef WITH_GRISU
const gchar *name, *status;
GList *list, *row;
GtkTreeModel *treemodel;
GtkTreeSelection *selection;
GtkTreeIter iter;

/* CURRENT - mutex lock - 1 grid op at a time */
if (!grid_task_start())
  return;

treemodel = gtk_tree_view_get_model(GTK_TREE_VIEW(job_tv));
selection = gtk_tree_view_get_selection(GTK_TREE_VIEW(job_tv));

list = gtk_tree_selection_get_selected_rows(selection, &treemodel);

for (row=list ; row ; row=g_list_next(row))
  {
  if (gtk_tree_model_get_iter(treemodel, &iter, row->data))
    {
    gtk_tree_model_get(treemodel, &iter, 0, &name, -1);

/* TODO - get job details? eg vo etc */
    status = grisu_job_status(name);
    if (status)
      gtk_list_store_set(job_ls, &iter, 1, status, -1);
    }
  }
grid_task_stop();
#endif
}

/****************************/
/* Download background task */
/****************************/
void grid_job_download_start(gpointer data)
{
GSList *list;
struct grid_data_pak *grid_data = data;

g_assert(grid_data != NULL);
g_assert(grid_data->jobname != NULL);

list = grid_data->jobname; 

grid_data->files = grid_download_job(list->data);
}

/*******************/
/* GUI update task */
/*******************/
void grid_job_download_stop(gpointer data)
{
GSList *list;
struct grid_data_pak *grid_data = data;

g_assert(grid_data != NULL);
g_assert(grid_data->jobname != NULL);

list = grid_data->jobname;

grid_process_job(list->data, grid_data->files);

grid_data_free(grid_data);
}

/*****************************************************/
/* Retrieve output from a selected job (if complete) */
/*****************************************************/
void gui_job_download(GtkWidget *w, gpointer dummy)
{
const gchar *name, *status;
GList *list, *row;
GtkTreeModel *treemodel;
GtkTreeSelection *selection;
GtkTreeIter iter;
struct grid_data_pak *data;

/* FIXME - need to make grisu_client (especially soap stuff) thread safe */
if (!grid_task_start())
  return;

/* get selected job names for downloading */
treemodel = gtk_tree_view_get_model(GTK_TREE_VIEW(job_tv));
selection = gtk_tree_view_get_selection(GTK_TREE_VIEW(job_tv));

list = gtk_tree_selection_get_selected_rows(selection, &treemodel);

for (row=list ; row ; row=g_list_next(row))
  {
  if (gtk_tree_model_get_iter(treemodel, &iter, row->data))
    {
    gtk_tree_model_get(treemodel, &iter, 0, &name, 1, &status, -1);
    if (g_ascii_strncasecmp(status, "done", 4) == 0)
      {
/* FIXME - won't behave correctly for multiple downloads */
/* CURRENT - mutex lock - 1 grid op at a time */
      data = grid_data_new();
      data->jobname = g_slist_prepend(data->jobname, g_strdup(name));

      task_new("Download", grid_job_download_start, data, grid_job_download_stop, data, NULL);
      }
    }
  }

/* release the mutex lock because downloads are allowed to be queued as long */
/* as they are not done simulataneously (until soap stuff is thread safe) */
grid_task_stop();
}

/******************************/
/* background job remove task */
/******************************/
void grid_job_remove_start(gpointer data)
{
#ifdef WITH_GRISU
GSList *list;
struct grid_data_pak *grid_data = data;

g_assert(grid_data != NULL);

for (list=grid_data->jobname ; list ; list=g_slist_next(list))
  grisu_job_remove(list->data);
#endif
}

/******************/
/* completed task */
/******************/
void grid_job_remove_stop(gpointer data)
{
#ifdef WITH_GRISU
struct grid_data_pak *grid_data = data;

g_assert(grid_data != NULL);

grid_data_free(grid_data);
grid_task_stop();
#endif
}

/*************************************/
/* stop and cleanup the selected job */
/*************************************/
void gui_job_remove(GtkWidget *w, gpointer dummy)
{
#ifdef WITH_GRISU
const gchar *name;
GList *list, *row;
GtkTreeModel *treemodel;
GtkTreeSelection *selection;
GtkTreeIter iter;
struct grid_data_pak *grid_data;

/* CURRENT - mutex lock - 1 grid op at a time */
if (!grid_task_start())
  return;

grid_data = grid_data_new();

treemodel = gtk_tree_view_get_model(GTK_TREE_VIEW(job_tv));
selection = gtk_tree_view_get_selection(GTK_TREE_VIEW(job_tv));

/* NB: remove in reverse order so we dont get the iterators confused */
list = gtk_tree_selection_get_selected_rows(selection, &treemodel);
list = g_list_reverse(list);

for (row=list ; row ; row=g_list_next(row))
  {
  if (gtk_tree_model_get_iter(treemodel, &iter, row->data))
    {
    gtk_tree_model_get(treemodel, &iter, 0, &name, -1);
//    grisu_job_remove(name);
    grid_data->jobname = g_slist_prepend(grid_data->jobname, g_strdup(name));

    gtk_list_store_remove(job_ls, &iter);
    }
  }

task_new("Job remove", grid_job_remove_start, grid_data, grid_job_remove_stop, grid_data, NULL);
#endif
}

/***********************/
/* job monitoring page */
/***********************/
void gui_grid_monitoring_page(GtkWidget *box, gpointer dialog)
{
gchar *titles[] = {" Name ", " Status " };
GtkWidget *swin, *hbox, *vbox;
GtkCellRenderer *renderer;
GtkTreeViewColumn *column;
GtkTreeSelection *selection;

g_assert(box != NULL);

vbox = gui_frame_vbox(NULL, TRUE, TRUE, box);

swin = gtk_scrolled_window_new(NULL, NULL);
gtk_scrolled_window_set_policy(GTK_SCROLLED_WINDOW(swin),
                               GTK_POLICY_AUTOMATIC, GTK_POLICY_AUTOMATIC);
gtk_box_pack_start(GTK_BOX(vbox), swin, TRUE, TRUE, 0);

job_ls = gtk_list_store_new(2, G_TYPE_STRING, G_TYPE_STRING);
job_tv = gtk_tree_view_new_with_model(GTK_TREE_MODEL(job_ls));
gtk_scrolled_window_add_with_viewport(GTK_SCROLLED_WINDOW(swin), job_tv);

renderer = gtk_cell_renderer_text_new();
column = gtk_tree_view_column_new_with_attributes(titles[0], renderer, "text", 0, NULL);
gtk_tree_view_append_column(GTK_TREE_VIEW(job_tv), column);
gtk_tree_view_column_set_expand(column, FALSE);

renderer = gtk_cell_renderer_text_new();
column = gtk_tree_view_column_new_with_attributes(titles[1], renderer, "text", 1, NULL);
gtk_tree_view_append_column(GTK_TREE_VIEW(job_tv), column);
gtk_tree_view_column_set_expand(column, TRUE);

selection = gtk_tree_view_get_selection(GTK_TREE_VIEW(job_tv));
gtk_tree_selection_set_mode(selection, GTK_SELECTION_MULTIPLE);

hbox = gtk_hbox_new(TRUE, PANEL_SPACING);
gtk_box_pack_start(GTK_BOX(vbox), hbox, FALSE, FALSE, 0);

/* TODO ... then automate */
gui_button(" Refresh all", gui_job_refresh_all, NULL, hbox, 0);
/* TODO - refresh selected */
gui_button(" Refresh ", gui_job_refresh_selected, NULL, hbox, 0);
/* upload is part of submit job */
//gui_button(" Upload ", gui_job_upload, NULL, hbox, 0);
gui_button(" Download ", gui_job_download, NULL, hbox, 0);
//gui_button(" Stop ", NULL, NULL, hbox, 0);
gui_button(" Remove ", gui_job_remove, NULL, hbox, 0);

}

/******************************************/
/* background task for application search */
/******************************************/
void grid_find_start(gpointer data)
{
struct grid_data_pak *grid_data = data;

g_assert(grid_data != NULL);

grid_data->queues = grid_search_by_name(grid_data->application);
}

/*******************/
/* GUI update task */
/*******************/
void grid_find_stop(gpointer data)
{
GSList *list;
GtkTreeIter iter;
struct grid_data_pak *grid_data = data;

g_assert(grid_data != NULL);

gtk_list_store_clear(grisu_ls);

for (list=grid_data->queues ; list ; list=g_slist_next(list))
  {
  gtk_list_store_append(grisu_ls, &iter);
  gtk_list_store_set(grisu_ls, &iter, 0, grid_data->application, 1, list->data, -1);
  }

grid_data_free(grid_data);

grid_task_stop();
}

/****************************/
/* connect to host callback */
/****************************/
void gui_grid_find(GtkWidget *w, gpointer dialog)
{
GtkEntry *entry;
struct grid_data_pak *grid_data;

/* CURRENT - mutex lock - 1 grid op at a time */
if (!grid_task_start())
  return;

entry = dialog_child_get(dialog, "application");

grid_data = grid_data_new();
grid_data->application = g_strdup(gtk_entry_get_text(entry));

/* clear old list */
task_new("Discovery", grid_find_start, grid_data, grid_find_stop, grid_data, NULL);
}

/****************************/
/* populate registered list */
/****************************/
void gui_grid_register_refresh(void)
{
gchar *key, *value;
GList *item, *list;
GtkTreeIter iter;

list = grid_application_all();

/* clear old GUI list */
gtk_list_store_clear(grid_ls);

for (item=list ; item ; item=g_list_next(item))
  {
  key = item->data;
  value = grid_application_get(key);

  gtk_list_store_append(grid_ls, &iter);
  gtk_list_store_set(grid_ls, &iter, 0, key, 1, value, -1);
  }
}

/*********************************/
/* set a new submission location */
/*********************************/
void gui_grid_register(GtkWidget *w, gpointer dummy)
{
const gchar *text1, *text2;
GtkTreeModel *treemodel;
GtkTreeSelection *selection;
GtkTreeIter iter;

treemodel = gtk_tree_view_get_model(GTK_TREE_VIEW(grisu_tv));
selection = gtk_tree_view_get_selection(GTK_TREE_VIEW(grisu_tv));

if (gtk_tree_selection_get_selected(selection, &treemodel, &iter))
  {
/* get application and location from selected item */
  gtk_tree_model_get(treemodel, &iter, 0, &text1, 1, &text2, -1);

//printf("Registering [%s][%s]\n", text1, text2);
//printf("Site = %s\n", grisu_site_name(text2));

/* register in hash table */
  grid_application_set(text1, text2);

/* refresh registered GUI list */
  gui_grid_register_refresh();
  }
else
  printf("Please select an entry from the Discovery list!\n");
}

/***********************************/
/* remove a registered application */
/***********************************/
void gui_grid_remove(GtkWidget *w, gpointer dummy)
{
const gchar *key;
GtkTreeModel *treemodel;
GtkTreeSelection *selection;
GtkTreeIter iter;

treemodel = gtk_tree_view_get_model(GTK_TREE_VIEW(grid_tv));
selection = gtk_tree_view_get_selection(GTK_TREE_VIEW(grid_tv));

if (gtk_tree_selection_get_selected(selection, &treemodel, &iter))
  {
/* get application from selected item */
  gtk_tree_model_get(treemodel, &iter, 0, &key, -1);

//printf("Removing: %s\n", key);

  grid_application_remove(key);

  gui_grid_register_refresh();
  }
 else
  printf("Please select an entry from the Registered list!\n");
}

/*****************************************/
/* setup/submit a registered application */
/*****************************************/
void gui_grid_job(GtkWidget *w, gpointer dummy)
{
#ifdef WITH_GRISU
gint type;
const gchar *name;
GtkTreeModel *treemodel;
GtkTreeSelection *selection;
GtkTreeIter iter;

if (!sysenv.active_model)
  {
printf("No active model.\n");
  return;
  }

treemodel = gtk_tree_view_get_model(GTK_TREE_VIEW(grid_tv));
selection = gtk_tree_view_get_selection(GTK_TREE_VIEW(grid_tv));

if (gtk_tree_selection_get_selected(selection, &treemodel, &iter))
  {
/* get application name from selected item */
  gtk_tree_model_get(treemodel, &iter, 0, &name, -1);

  type = grisu_job_type(name);

  switch (type)
    {
    case JOB_GULP:
      gulp_dialog();
      break;

    case JOB_GAMESS:
      gamess_dialog();
      break;
    }
  }
else
  {
  printf("Nothing selected.\n");
  }
#endif
}

/************************/
/* host setup/info page */
/************************/
void gui_grid_discover_page(GtkWidget *box, gpointer dialog)
{
gchar *titles[] = {" Application ", " Submission location " };
GtkWidget *swin, *hbox2, *vbox, *entry, *button;
GtkCellRenderer *renderer;
GtkTreeViewColumn *column;

/* left pane */
vbox = gui_frame_vbox(NULL, TRUE, TRUE, box);

/* available applications/queues list */
swin = gtk_scrolled_window_new(NULL, NULL);
gtk_scrolled_window_set_policy(GTK_SCROLLED_WINDOW(swin),
                               GTK_POLICY_AUTOMATIC, GTK_POLICY_AUTOMATIC);
gtk_box_pack_start(GTK_BOX(vbox), swin, TRUE, TRUE, 0);

grisu_ls = gtk_list_store_new(2, G_TYPE_STRING, G_TYPE_STRING);
grisu_tv = gtk_tree_view_new_with_model(GTK_TREE_MODEL(grisu_ls));
gtk_scrolled_window_add_with_viewport(GTK_SCROLLED_WINDOW(swin), grisu_tv);

renderer = gtk_cell_renderer_text_new();
column = gtk_tree_view_column_new_with_attributes(titles[0], renderer, "text", 0, NULL);
gtk_tree_view_append_column(GTK_TREE_VIEW(grisu_tv), column);
gtk_tree_view_column_set_expand(column, FALSE);

renderer = gtk_cell_renderer_text_new();
column = gtk_tree_view_column_new_with_attributes(titles[1], renderer, "text", 1, NULL);
gtk_tree_view_append_column(GTK_TREE_VIEW(grisu_tv), column);
gtk_tree_view_column_set_expand(column, TRUE);

/* discovery functionality */
hbox2 = gtk_hbox_new(FALSE, 0);
gtk_box_pack_end(GTK_BOX(vbox), hbox2, FALSE, FALSE, 0);

/* search pattern */
entry = gtk_entry_new();
dialog_child_set(dialog, "application", entry);
gtk_box_pack_start(GTK_BOX(hbox2), entry, TRUE, TRUE, 0);
/* hit enter to trigger search */
g_signal_connect(GTK_OBJECT(entry), "activate", GTK_SIGNAL_FUNC(gui_grid_find), dialog);

/* search action */
button = gtk_button_new_with_label(" Search ");
g_signal_connect(GTK_OBJECT(button), "clicked", GTK_SIGNAL_FUNC(gui_grid_find), dialog);
gtk_box_pack_start(GTK_BOX(hbox2), button, FALSE, FALSE, 0);

/* register action */
button = gtk_button_new_with_label(" Register ");
gtk_box_pack_start(GTK_BOX(hbox2), button, FALSE, FALSE, 0);
g_signal_connect(GTK_OBJECT(button), "clicked", GTK_SIGNAL_FUNC(gui_grid_register), NULL);
}

/*********************************/
/* job setup and submission page */
/*********************************/
void gui_grid_submit_page(GtkWidget *box, gpointer dialog)
{
gchar *titles[] = {" Application ", " Submission location " };
GtkWidget *hbox2, *vbox, *button;
GtkCellRenderer *renderer;
GtkTreeViewColumn *column;

vbox = gui_frame_vbox(NULL, TRUE, TRUE, box);

/* registered applications */
grid_ls = gtk_list_store_new(2, G_TYPE_STRING, G_TYPE_STRING);
grid_tv = gtk_tree_view_new_with_model(GTK_TREE_MODEL(grid_ls));
gtk_box_pack_start(GTK_BOX(vbox), grid_tv, TRUE, TRUE, 0);

renderer = gtk_cell_renderer_text_new();
column = gtk_tree_view_column_new_with_attributes(titles[0], renderer, "text", 0, NULL);
gtk_tree_view_append_column(GTK_TREE_VIEW(grid_tv), column);
gtk_tree_view_column_set_expand(column, TRUE);

renderer = gtk_cell_renderer_text_new();
column = gtk_tree_view_column_new_with_attributes(titles[1], renderer, "text", 1, NULL);
gtk_tree_view_append_column(GTK_TREE_VIEW(grid_tv), column);
gtk_tree_view_column_set_expand(column, TRUE);

/* registered functionality */
hbox2 = gtk_hbox_new(TRUE, 0);
gtk_box_pack_end(GTK_BOX(vbox), hbox2, FALSE, FALSE, 0);

button = gtk_button_new_with_label(" Unregister ");
gtk_box_pack_start(GTK_BOX(hbox2), button, FALSE, FALSE, 0);
g_signal_connect(GTK_OBJECT(button), "clicked", GTK_SIGNAL_FUNC(gui_grid_remove), NULL);

/* TODO - configure job button (or auto pop up on submit) */
/* CURRENT - submit a job that tests the selected object */
button = gtk_button_new_with_label(" Setup job ");
gtk_box_pack_start(GTK_BOX(hbox2), button, FALSE, FALSE, 0);
g_signal_connect(GTK_OBJECT(button), "clicked", GTK_SIGNAL_FUNC(gui_grid_job), NULL);

/* init GUI */
gui_grid_register_refresh();
}

/***************************/
/* configuration callbacks */
/***************************/
void gui_grid_ws_change(GtkWidget *w, gpointer dummy)
{
#ifdef WITH_GRISU
const gchar *ws;
ws = gtk_entry_get_text(GTK_ENTRY(grid_service_entry));
grisu_ws_set(ws);
#endif
}

void gui_grid_server_change(GtkWidget *w, gpointer dummy)
{
#ifdef WITH_GRISU
const gchar *server;
server = gtk_entry_get_text(GTK_ENTRY(grid_myproxy_server_entry));
grisu_myproxy_set(server, NULL);
#endif
}

void gui_grid_port_change(GtkWidget *w, gpointer dummy)
{
#ifdef WITH_GRISU
const gchar *port;
port = gtk_entry_get_text(GTK_ENTRY(grid_myproxy_port_entry));
grisu_myproxy_set(NULL, port);
#endif
}

/***********************/
/* configuration panel */
/***********************/
void gui_grid_configure_page(GtkWidget *box)
{
GtkWidget *table, *label;

g_assert(box != NULL);

table = gtk_table_new(2, 3, FALSE);
gtk_container_add(GTK_CONTAINER(box), table);

label = gtk_label_new("Service interface");
gtk_table_attach(GTK_TABLE(table),label,0,1,0,1, GTK_SHRINK, GTK_SHRINK, 0, 0);

label = gtk_label_new("Myproxy server ");
gtk_table_attach(GTK_TABLE(table),label,0,1,1,2, GTK_SHRINK, GTK_SHRINK, 0, 0);

label = gtk_label_new("Mpyroxy port ");
gtk_table_attach(GTK_TABLE(table),label,0,1,2,3, GTK_SHRINK, GTK_SHRINK, 0, 0);

grid_service_entry = gtk_entry_new();
gtk_table_attach_defaults(GTK_TABLE(table),grid_service_entry,1,2,0,1);
g_signal_connect(GTK_OBJECT(grid_service_entry), "changed",
                 GTK_SIGNAL_FUNC(gui_grid_ws_change), NULL);

grid_myproxy_server_entry = gtk_entry_new();
gtk_table_attach_defaults(GTK_TABLE(table),grid_myproxy_server_entry,1,2,1,2);
g_signal_connect(GTK_OBJECT(grid_myproxy_server_entry), "changed",
                 GTK_SIGNAL_FUNC(gui_grid_server_change), NULL);

grid_myproxy_port_entry = gtk_entry_new();
gtk_table_attach_defaults(GTK_TABLE(table),grid_myproxy_port_entry,1,2,2,3);
g_signal_connect(GTK_OBJECT(grid_myproxy_port_entry), "changed",
                 GTK_SIGNAL_FUNC(gui_grid_port_change), NULL);
}

/********************************************/
/* close dialog (if nothing is in progress) */
/********************************************/
void gui_grid_close(GtkWidget *w, gpointer data)
{
/* check for background task mutex */
if (!grid_task_start())
 return;

/* destroy */
dialog_destroy(NULL, data);

/* clear mutex */
grid_task_stop();
}

/*****************************/
/* the gdis job setup dialog */
/*****************************/
void gui_grid_dialog(void)
{
gpointer dialog;
GtkWidget *window, *vbox, *notebook, *page, *label;

if (!gui_grid_mutex)
  gui_grid_mutex = g_mutex_new();

/* CURRENT */
grid_setup();

dialog = dialog_request(666, "Grid Manager", NULL, NULL, NULL);
if (!dialog)
  return;
window = dialog_window(dialog);

vbox = gtk_vbox_new(FALSE, PANEL_SPACING);
gtk_container_add(GTK_CONTAINER(GTK_DIALOG(window)->vbox), vbox);

gtk_widget_set_size_request(window, -1, 440);

/* authentication display */
page = gui_frame_vbox("Authentication", FALSE, FALSE, vbox);
gui_grid_setup_page(page, dialog);

/* configuration display */
grid_configure_box = gtk_frame_new("Configuration");
gtk_box_pack_start(GTK_BOX(vbox), grid_configure_box, FALSE, FALSE, 0);
gui_grid_configure_page(grid_configure_box);

/* notebook functionality splitting */
notebook = gtk_notebook_new();
gtk_notebook_set_tab_pos(GTK_NOTEBOOK(notebook), GTK_POS_TOP);
gtk_box_pack_start(GTK_BOX(vbox), notebook, TRUE, TRUE, 0);
gtk_notebook_set_show_border(GTK_NOTEBOOK(notebook), FALSE);

/* discovery page */
page = gtk_vbox_new(FALSE, 0);
label = gtk_label_new(" Discovery ");
gtk_notebook_append_page(GTK_NOTEBOOK(notebook), page, label);
gui_grid_discover_page(page, dialog);

/* submission/job configuration page */
page = gtk_vbox_new(FALSE, 0);
label = gtk_label_new(" Submission ");
gtk_notebook_append_page(GTK_NOTEBOOK(notebook), page, label);
gui_grid_submit_page(page, dialog);

/* job monitoring page */
page = gtk_vbox_new(FALSE, 0);
label = gtk_label_new(" Monitoring ");
gtk_notebook_append_page(GTK_NOTEBOOK(notebook), page, label);
gui_grid_monitoring_page(page, dialog);

/* CURRENT - indicate action in progress */
/* TODO - not really intended for download progress (as there may be more than 1) */
/* but rather as user feedback during (slow) grid responses to even simple commands */
grid_status_bar = gtk_entry_new();
gtk_entry_set_editable(GTK_ENTRY(grid_status_bar), FALSE);
gtk_box_pack_end(GTK_BOX(vbox), grid_status_bar, FALSE, FALSE, 0);

/* TODO - disallow even X close button if in the middle of a grid operation */
gui_stock_button(GTK_STOCK_CLOSE, gui_grid_close, dialog, GTK_DIALOG(window)->action_area);

gtk_widget_show_all(window);

/* default to hiden */
gtk_widget_hide(grid_configure_box);
}
