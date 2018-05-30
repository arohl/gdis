/*
Copyright (C) 2000 by Sean David Fleming

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

/* irix */
#define _BSD_SIGNALS 1

#include <fcntl.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#include <signal.h>
#include <time.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>


#include "gdis.h"
#include "task.h"
#include "file.h"
#include "parse.h"
#include "gui_shorts.h"
#include "interface.h"
#include "dialog.h"

/* top level data structure */
extern struct sysenv_pak sysenv;

#ifndef __WIN32
#include <sys/wait.h>
#endif

/* task manager globals */
gdouble max_threads=1;
gint task_show_running = 1;
gint task_show_queued = 0;
gint task_show_completed = 1;

GtkWidget *task_list=NULL;
gint selected_task = -1;
gint task_info_pid = 0;
gint task_thread_pid = -1;
GArray *psarray = NULL;
GtkListStore *task_list_ls=NULL;
GtkWidget *task_list_tv=NULL, *task_text_view;
char *strptime(const char *, const char *, struct tm *);

/**************************/
/* a simple sleep routine */
/**************************/
void delay(gint usecs)
{
#ifndef __WIN32
struct timespec req;

req.tv_sec = usecs / 1000000;
req.tv_nsec = 1000*(usecs - 1000000*req.tv_sec);

/* TODO - if terminated early - call again with remainder */
nanosleep(&req, NULL);
#endif
}

#define EXIST(idx) ((idx) != -1)

gint find_pid(struct task_pak *curr_process, gint *pid)
{
if (curr_process->ppid == GPOINTER_TO_INT(pid))
  return (0);
else
  return(-1);
}

gint find_pid_index(gint pid)
{
gint i;

  for (i = psarray->len - 1; i >= 0 && g_array_index(psarray, struct task_pak, i).pid != pid; i--);
  return(i);
}

#define DEBUG_ADD_FROM_TREE 0
void add_from_tree(struct task_pak *tdata, gint idx)
{
struct task_pak *process_ptr;
gint child;
div_t h_sec, sec, min;
gchar *h_sec_pos, *sec_pos, *min_pos, *hour_pos;

/* add in current process */  
process_ptr = &g_array_index(psarray, struct task_pak, idx);
tdata->pcpu += process_ptr->pcpu;
tdata->pmem += process_ptr->pmem;
/* parse the cpu time */
#if defined(__APPLE__) && defined(__MACH__) 
h_sec_pos = g_strrstr(process_ptr->time, ".") + 1;
sec_pos = h_sec_pos - 3;
min_pos = process_ptr->time;
hour_pos = NULL;
#else
hour_pos = process_ptr->time;
min_pos = hour_pos + 3;
sec_pos = min_pos + 3;
h_sec_pos = NULL;
#endif
tdata->h_sec += (gint) str_to_float(h_sec_pos);
h_sec = div(tdata->h_sec, 100);
tdata->h_sec = h_sec.rem;
tdata->sec =tdata->sec + h_sec.quot + (gint) str_to_float(sec_pos);
sec = div(tdata->sec, 60);
tdata->sec = sec.rem;
tdata->min = tdata->min + sec.quot + (gint) str_to_float(min_pos);
min = div(tdata->min, 60);
tdata->min = min.rem;
tdata->hour = tdata->hour + min.quot + (gint) str_to_float(hour_pos);
#if DEBUG_ADD_FROM_TREE
printf("%s hour=%d min=%d, sec=%d, hsec=%d\n", process_ptr->time, tdata->hour, tdata->min, tdata->sec, tdata->h_sec);
#endif

/* process children */
for (child = process_ptr->child; EXIST(child); child = g_array_index(psarray, struct task_pak, child).sister)
  add_from_tree(tdata, child);
}

/*****************************************************/
/* get the results of a ps on the requested process */
/*****************************************************/
#define DEBUG_CALC_TASK_INFO 0
void calc_task_info(struct task_pak *tdata)
{
gint num_tokens;
gint me, parent, sister;
gint found;
gchar line[LINELEN], *line_no_time, **buff, *cmd;
struct task_pak curr_process, *process_ptr;
struct tm start_tm;
time_t oldest_time;
FILE *fp;

#ifdef _WIN32
/* it's just not worth it ... */
return;
#endif

/* dispose of task list if it exists */
if (psarray != NULL)
  g_array_free(psarray, TRUE);

/* initialise process record */
curr_process.parent = curr_process.child = curr_process.sister = -1;

/* setup ps command */

/*
 * SG's - 'ps -Ao "pid ppid pcpu vsz etime comm"
 * 			(unfortunately, no %mem - and vsz is absolute :-(
 */

#ifdef __sgi
cmd = g_strdup_printf("ps -Ao \"pid ppid pcpu vsz etime comm\"");
#else
cmd = g_strdup_printf("ps axo \"lstart pid ppid %%cpu %%mem time ucomm\"");
#endif

/* run ps command */
fp = popen(cmd, "r");
g_free(cmd);
if (!fp)
  {
  gui_text_show(ERROR, "unable to launch ps command\n");
  return;
  }

/* skip title line */
if (fgetline(fp, line))
  {
  gui_text_show(ERROR, "unable to read first line of output from ps command\n");
  return;
  }
    
/* load data into array */
psarray = g_array_new(FALSE, FALSE, sizeof(struct task_pak));
while (!fgetline(fp, line))
  {
#if DEBUG_CALC_TASK_INFO
printf("%s", line);
#endif
/* extract what we want */
#ifdef __WIN32
/* TODO - win32 replacement? */
  line_no_time = NULL;
#else
/* first get start time */
  line_no_time = strptime(line, "%c", &start_tm);
#endif
  if (line_no_time == NULL)
    {
    gui_text_show(ERROR, "file produced by ps not understood\n");
    return;
    }
  curr_process.start_time = mktime(&start_tm);
  buff = tokenize(line_no_time, &num_tokens);
  if (num_tokens >= 5)
    {
    curr_process.pid = (int) str_to_float(*(buff+0));
    curr_process.ppid = (int) str_to_float(*(buff+1));
    curr_process.pcpu = str_to_float(*(buff+2));
    curr_process.pmem = str_to_float(*(buff+3));
    curr_process.time = g_strdup(*(buff+4));
    g_array_append_val(psarray, curr_process);
    }
  g_strfreev(buff);
  }
pclose(fp);

/* Build the process hierarchy. Every process marks itself as first child */
/* of it's parent or as sister of first child of its parent */
/* algorithm taken from pstree.c by Fred Hucht */

for (me = 0; me < psarray->len; me++)
  {
  parent = find_pid_index(g_array_index(psarray, struct task_pak, me).ppid);
  if (parent != me && parent != -1)
    { /* valid process, not me */
    g_array_index(psarray, struct task_pak, me).parent = parent;
    if (g_array_index(psarray, struct task_pak, parent).child == -1) /* first child */
      g_array_index(psarray, struct task_pak, parent).child = me;
    else
      {
      for (sister = g_array_index(psarray, struct task_pak, parent).child; EXIST(g_array_index(psarray, struct task_pak, sister).sister); sister = g_array_index(psarray, struct task_pak, sister).sister);
      g_array_index(psarray, struct task_pak, sister).sister = me;
      }
    }
  }

#if DEBUG_CALC_TASK_INFO
printf("process list\n");
for (me = 0; me < psarray->len; me++)
  {
  process_ptr = &g_array_index(psarray, struct task_pak, me);
  printf("pid: %d ppid: %d pcpu: %.1f pmem: %.1f date: %s parent: %d child: %d sister: %d time: %s", process_ptr->pid,  process_ptr->ppid, process_ptr->pcpu, process_ptr->pmem, process_ptr->time, process_ptr->parent, process_ptr->child, process_ptr->sister, ctime(&process_ptr->start_time));
  }
#endif

/* set tdata's pid to 0 so if this routine fails, it won't kill gdis! */
tdata->pcpu =  tdata->pmem = 0;
tdata->h_sec = tdata->sec = tdata->min = tdata->hour = 0;

/* scan through the running processes and find the task thread */
oldest_time = time(NULL);
found = FALSE;
for (me = 0; me < psarray->len; me++)
  {
  process_ptr = &g_array_index(psarray, struct task_pak, me);
  if (process_ptr->pid == tdata->pid)
    {
    found = TRUE;

      #if DEBUG_CALC_TASK_INFO
      printf("pid: %d time:%s\n", child_pid, ctime(&oldest_time));
      #endif
    break;
    }
  }
if (found)
  {
    /* find child */
/*    process_ptr = &g_array_index(psarray, struct task_pak, child_idx);
    child_idx = process_ptr->child;
    process_ptr = &g_array_index(psarray, struct task_pak, child_idx);
    child_pid = tdata->pid = process_ptr->pid;
    #if DEBUG_CALC_TASK_INFO
    printf("pid: %d time:%s\n", child_pid, ctime(&oldest_time));
    #endif*/

  add_from_tree(tdata, me);
  }

/* FIXME - core dumps here in analysis (esp. if >1 jobs queued) */
if (tdata->time)
  g_free(tdata->time);

/*
sprintf(line, "%02d:%02d:%02d", tdata->hour, tdata->min, tdata->sec);
tdata->time = g_strdup(line);
*/
tdata->time = g_strdup_printf("%02d:%02d:%02d", tdata->hour, tdata->min, tdata->sec);
 

/* search through any children of the the oldest child */
/*
pid = child_pid;
while ((list = g_slist_find_custom(palist, GINT_TO_POINTER(pid), (gpointer) find_pid)) != NULL)
  {
  curr_process = (struct task_pak *) list->data;
  pid = tdata->pid = curr_process->pid;
  tdata->pcpu += curr_process->pcpu;
  tdata->pmem += curr_process->pmem;
  if (tdata->time)
      g_free(tdata->time);
  tdata->time = g_strdup(curr_process->time);
  }
#if DEBUG_CALC_TASK_INFO
printf("FINAL TOTALS\n");
printf("pid: %d pcpu: %.1f pmem: %.1f time: %s\n", tdata->pid, tdata->pcpu, tdata->pmem, tdata->time);
#endif

free_slist(palist);
*/

/* prevent rounding errors giving >100% */
if (tdata->pcpu > 100.0)
  tdata->pcpu = 100.0;
if (tdata->pmem > 100.0)
  tdata->pmem = 100.0;
}

/**********************************************/
/* get pointer to the currently selected task */
/**********************************************/
gpointer task_selected(void)
{
gpointer task;
GtkTreeModel *treemodel;
GtkTreeSelection *selection;
GtkTreeIter iter;

treemodel = gtk_tree_view_get_model(GTK_TREE_VIEW(task_list_tv));
selection = gtk_tree_view_get_selection(GTK_TREE_VIEW(task_list_tv));

if (gtk_tree_model_get_iter_first(treemodel, &iter))
  {
  do
    {
    if (gtk_tree_selection_iter_is_selected(selection, &iter))
      {
      gtk_tree_model_get(treemodel, &iter, 6, &task, -1);
      return(task);
      }
    }
  while (gtk_tree_model_iter_next(treemodel, &iter));
  }
return(NULL);
}

/**********************************/
/* get the currently selected row */
/**********************************/
gint task_selected_row(void)
{
gint row=0;
GtkTreeModel *treemodel;
GtkTreeSelection *selection;
GtkTreeIter iter;

treemodel = gtk_tree_view_get_model(GTK_TREE_VIEW(task_list_tv));
selection = gtk_tree_view_get_selection(GTK_TREE_VIEW(task_list_tv));

if (gtk_tree_model_get_iter_first(treemodel, &iter))
  {
  do
    {
    if (gtk_tree_selection_iter_is_selected(selection, &iter))
      return(row);
    row++;
    }
  while (gtk_tree_model_iter_next(treemodel, &iter));
  }
return(-1);
}

/*********************************/
/* rewrite the current task list */
/*********************************/
gint update_task_info(void)
{
gint i, type, update, row, num_rows=0;
gchar *txt;
GSList *tlist=NULL;
GtkTreeIter iter;
GtkTreeModel *treemodel;
GtkTreeSelection *selection;
GtkTextMark *mark;
GtkTextIter text_iter;
GtkTextBuffer *buffer;
struct task_pak *task=NULL;
static gpointer previous_task=NULL;

/* checks */
if (!dialog_exists(TASKMAN, NULL))
  return(TRUE);
if (!task_list_ls)
  return(TRUE);

/* NEW - update using control file */
task = task_selected();
buffer = gtk_text_view_get_buffer(GTK_TEXT_VIEW(task_text_view));
if (task)
  {
  if (previous_task != task)
    {
/* complete refresh */
    task->status_index = 0;
    gtk_text_buffer_set_text(buffer, "", 0);
    update = TRUE;
    }
  else
    update = FALSE;

  switch (task->status)
    {
    case RUNNING:
/* running - always update */
      update = TRUE;
      break;

    case KILLED:
    case COMPLETED:
      if (task->status_fp)
        {
/* completed - only update if status file is still open */
        task->status_index = 0;
        update = TRUE;
        }
      break;
    }

  if (update)
    {
/* refresh the status message and print */
    task_status_update(task);
    txt = (task->status_text)->str;
    if (task->status_index >= 0)
      {
      i = strlen(txt+task->status_index);
      if (i)
        {
/* display only most recent lines */
/*
        gtk_text_buffer_set_text(buffer, txt+task->status_index, i);
*/

/* append to the whole thing */
        gtk_text_buffer_insert_at_cursor(buffer, txt+task->status_index, i);
        task->status_index += i;

/* force scroll to end of buffer by default */
        mark = gtk_text_buffer_get_insert(buffer);
        gtk_text_buffer_get_end_iter(buffer, &text_iter);
        gtk_text_buffer_move_mark(buffer, mark, &text_iter);
        gtk_text_view_scroll_to_mark(GTK_TEXT_VIEW(task_text_view),
                                     mark, 0.0, FALSE, 0, 0);
        }
      }
    }
  previous_task = task;
  }
else
  {
  gtk_text_buffer_set_text(buffer, "", 0);
  previous_task = NULL;
  }

/* remove deleted tasks */
/* CURRENT - NEVER allow this? so user can also look at status file of completed jobs */
tlist = sysenv.task_list;
while (tlist)
  {
  task = tlist->data;
  tlist = g_slist_next(tlist);
  if (task->status == REMOVED)
    {
    sysenv.task_list = g_slist_remove(sysenv.task_list, task);
    task_free(task);
    }
  }

/* clear list and re-populate according to status */
row = task_selected_row();
gtk_list_store_clear(task_list_ls);
for (type=RUNNING ; type<=COMPLETED ; type++)
  {
/* fill out the list with current status type */
  for (tlist=sysenv.task_list ; tlist ; tlist=g_slist_next(tlist))
    {
    task = tlist->data;
    if (task->status != type)
      continue;

    gtk_list_store_append(task_list_ls, &iter);

    if (task->label)
      gtk_list_store_set(task_list_ls, &iter, 1, task->label, -1);
    else
      gtk_list_store_set(task_list_ls, &iter, 1, "Unknown", -1);

    switch(task->status)
      {
      case QUEUED:
        gtk_list_store_set(task_list_ls, &iter, 2, "Queued", -1);
        break;

      case RUNNING:
        calc_task_info(task);
        txt = g_strdup_printf("%d", task->pid);
        gtk_list_store_set(task_list_ls, &iter, 0, txt, -1);
        g_free(txt);
        if (task->progress > 0.0)
          txt = g_strdup_printf("%5.1f%%", task->progress);
        else
          txt = g_strdup("Running");
        gtk_list_store_set(task_list_ls, &iter, 2, txt, -1);
        g_free(txt);
        txt = g_strdup_printf("%5.1f", task->pcpu);
        gtk_list_store_set(task_list_ls, &iter, 3, txt, -1);
        g_free(txt);
        txt = g_strdup_printf("%5.1f", task->pmem);
        gtk_list_store_set(task_list_ls, &iter, 4, txt, -1);
        g_free(txt);
        txt = g_strdup_printf(" %s", task->time);
        gtk_list_store_set(task_list_ls, &iter, 5, txt, -1);
        g_free(txt);
        break;

      case KILLED:
        gtk_list_store_set(task_list_ls, &iter, 2, "Killed", -1);
        break;

      case COMPLETED:
        gtk_list_store_set(task_list_ls, &iter, 2, "Completed", -1);
        break;
      }
  
/* NEW - add the task pointer itself */
    gtk_list_store_set(task_list_ls, &iter, 6, task, -1);
    num_rows++;
    }
  }

/* restore selected row, or the previous (if possible) if deleted  */
if (row >= num_rows && row)
  row--;
if (row >= 0)
  {
  treemodel = gtk_tree_view_get_model(GTK_TREE_VIEW(task_list_tv));
  if (gtk_tree_model_iter_nth_child(treemodel, &iter, NULL, row))
    {
    selection = gtk_tree_view_get_selection(GTK_TREE_VIEW(task_list_tv)); 
    if (selection)
      gtk_tree_selection_select_iter(selection, &iter);
    }
  }

/* NEW */
sysenv.max_threads = max_threads;
g_thread_pool_set_max_threads(sysenv.thread_pool, sysenv.max_threads, NULL);
g_thread_pool_stop_unused_threads();

return(TRUE);
}

/***********************************/
/* recursively kill a process tree */
/***********************************/
void task_kill_tree(gint child)
{
gchar *cmd;
struct task_pak *process;

while (EXIST(child))
  {
  process = &g_array_index(psarray, struct task_pak, child);
/* recurse down the tree */
  task_kill_tree(process->child);
/* terminate process */
  cmd = g_strdup_printf("kill TERM %d >& /dev/null", process->pid);
  system(cmd);
  g_free(cmd);
/* get next process at this level */
  child = process->sister;
  }
}

/***********************/
/* stop a running task */
/***********************/
void task_kill_running(struct task_pak *task)
{
gint i;
struct task_pak *process;

/* get task and bottom level process pid */
if (task && psarray)
  {
  for (i=psarray->len ; i-- ; )
    {
    process = &g_array_index(psarray, struct task_pak, i);
    if (process->pid == task->pid)
      {
      task_kill_tree(process->child);
      break;
      }
    }
  task->status = KILLED;
  }
}

/***************************************/
/* remove and free all completed tasks */
/***************************************/
void task_remove_completed(GtkWidget *w, gpointer data)
{
GSList *list;
struct task_pak *task;

list = sysenv.task_list;
while (list)
  {
  task = list->data;
  list = g_slist_next(list);
  if (task->status == COMPLETED || task->status == KILLED)
    {
    sysenv.task_list = g_slist_remove(sysenv.task_list, task);
    task_free(task);
    }
  }
}

/**********************************************/
/* terminate selected running or queued tasks */
/**********************************************/
void task_kill_selected(void)
{
struct task_pak *task;

task = task_selected();
if (task)
  {
  switch (task->status)
    {
    case RUNNING:
      task_kill_running(task);
      break;

    case QUEUED:
      task->status = KILLED;
      break;
    }
  }
}

/*****************************************/
/* terminate all running or queued tasks */
/*****************************************/
void task_kill_all(void)
{
GSList *list;
struct task_pak *task;

/* remove all queued first - to ensure they aren't executed */
for (list=sysenv.task_list ; list ; list=g_slist_next(list))
  {
  task = list->data;
  if (task->status == QUEUED)
    task->status = KILLED;
  }

for (list=sysenv.task_list ; list ; list=g_slist_next(list))
  {
  task = list->data;
  if (task->status == RUNNING)
    task_kill_running(task);
  }
}

/***********************/
/* task dialog cleanup */
/***********************/
void task_cleanup(void)
{
gtk_list_store_clear(task_list_ls);
task_list_ls = NULL;
task_list_tv = NULL;
}

/****************************************************/
/* a task dialog for viewing/killing existing tasks */
/****************************************************/
void task_dialog(void)
{
gint i;
gchar *titles[6] = {"  PID  ", "       Job       ", "   Status   ", "  \% CPU  ", "  \% Mem  ", "  Time  "};
gpointer dialog;
GtkWidget *window, *swin, *frame, *vbox;
GtkCellRenderer *renderer;
GtkTreeViewColumn *column;

/* request a new dialog */
dialog = dialog_request(TASKMAN, "Task Manager", NULL, task_cleanup, NULL);
if (!dialog)
  return;
window = dialog_window(dialog);
gtk_window_set_default_size(GTK_WINDOW(window), 600, 600);

/* CURRENT - locking this to 1 so grid stuff is limited to 1 task */
/* this is because Im using global vars (esp soap struct) in the */
/* grisu_client code - so the whole thing isn't thread safe */
/* task thread setup */
//GtkWidget *hbox;
//
//hbox = gtk_hbox_new(FALSE, 0);
//gtk_box_pack_start(GTK_BOX(GTK_DIALOG(window)->vbox), hbox, FALSE, FALSE, 0); 
//gtk_container_set_border_width(GTK_CONTAINER(hbox), PANEL_SPACING);
//gui_direct_spin("Maximum number of running tasks", &max_threads, 1, 32, 1, NULL, NULL, hbox);

#if _WIN32
gtk_widget_set_sensitive(GTK_WIDGET(hbox), FALSE);
#endif

/* job list frame */
frame = gtk_frame_new(NULL);
gtk_box_pack_start(GTK_BOX(GTK_DIALOG(window)->vbox), frame, TRUE, TRUE, 0); 
gtk_container_set_border_width(GTK_CONTAINER(frame), PANEL_SPACING);
/* create a vbox in the frame */
vbox = gtk_vbox_new(FALSE, PANEL_SPACING/2);
gtk_container_add(GTK_CONTAINER(frame), vbox);
gtk_container_set_border_width(GTK_CONTAINER(GTK_BOX(vbox)), PANEL_SPACING);

/* scrolled model pane */
swin = gtk_scrolled_window_new(NULL, NULL);
gtk_scrolled_window_set_policy(GTK_SCROLLED_WINDOW (swin),
                               GTK_POLICY_AUTOMATIC, GTK_POLICY_AUTOMATIC);
gtk_box_pack_start(GTK_BOX(vbox), swin, TRUE, TRUE, 0);

task_list_ls = gtk_list_store_new(7, G_TYPE_STRING, G_TYPE_STRING, G_TYPE_STRING,
                                     G_TYPE_STRING, G_TYPE_STRING, G_TYPE_STRING,
                                     G_TYPE_POINTER);
task_list_tv = gtk_tree_view_new_with_model(GTK_TREE_MODEL(task_list_ls));
gtk_scrolled_window_add_with_viewport(GTK_SCROLLED_WINDOW(swin), task_list_tv);
for (i=0 ; i<6 ; i++)
  {
  renderer = gtk_cell_renderer_text_new();
  column = gtk_tree_view_column_new_with_attributes(titles[i], renderer, "text", i, NULL);
  gtk_tree_view_append_column(GTK_TREE_VIEW(task_list_tv), column);
  }

/* setup the selection handler */
/*
{
GtkTreeSelection *select;
select = gtk_tree_view_get_selection(GTK_TREE_VIEW(task_list_tv));
g_signal_connect(G_OBJECT(select), "changed",
                 G_CALLBACK(task_selection_changed),
                 NULL);
}
*/


/* selected task status file frame */
gtk_box_pack_start(GTK_BOX(vbox), gtk_hseparator_new(), FALSE, FALSE, 0);
swin = gtk_scrolled_window_new(NULL, NULL);
gtk_box_pack_start(GTK_BOX(vbox), swin, TRUE, TRUE, 0);
gtk_scrolled_window_set_policy(GTK_SCROLLED_WINDOW(swin),
                               GTK_POLICY_AUTOMATIC, GTK_POLICY_AUTOMATIC);
task_text_view = gtk_text_view_new();
gtk_text_view_set_editable(GTK_TEXT_VIEW(task_text_view), FALSE);
gtk_text_view_set_cursor_visible(GTK_TEXT_VIEW(task_text_view), FALSE);
gtk_container_add(GTK_CONTAINER(swin), task_text_view);


/* control buttons */
gui_icon_button(GTK_STOCK_REMOVE, "Kill selected",
                  task_kill_selected, NULL,
                  GTK_DIALOG(window)->action_area);

gui_icon_button(GTK_STOCK_DELETE, "Kill all",
                  task_kill_all, NULL,
                  GTK_DIALOG(window)->action_area);

gui_icon_button(GTK_STOCK_CLEAR, "Remove completed",
                  task_remove_completed, NULL,
                  GTK_DIALOG(window)->action_area);

gui_stock_button(GTK_STOCK_CLOSE,
                   dialog_destroy, dialog,
                   GTK_DIALOG(window)->action_area);

/* done */
gtk_widget_show_all(window);

/* refresh labels */
update_task_info();
}

