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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#include "gdis.h"
#include "file.h"
#include "task.h"
#include "job.h"
#include "grid.h"
#include "interface.h"

/* top level data structure */
extern struct sysenv_pak sysenv;

#ifndef __WIN32
#include <sys/wait.h>
#endif


/**********************************************/
/* execute task in thread created by the pool */
/**********************************************/
#define DEBUG_TASK_PROCESS 0
void task_process(struct task_pak *task, gpointer data)
{
/* checks */
if (!task)
  return;
if (task->status != QUEUED)
  return;

/* setup for current task */
if(!task->is_async) task->pid = getpid();
/* TODO - should be mutex locking this */
task->status = RUNNING;

/* NB: the primary task needs to not do anything that could */
/* cause problems eg update GUI elements since this can cause conflicts */

/* execute the primary task */
task->primary(task->ptr1, task);
if(task->is_async){
  pid_t w;
  gint status;
   /*here we NEED to wait for async task to die*/
  do{
    w = waitpid(task->pid, &status,WUNTRACED|WCONTINUED);
    if(w==-1) fprintf(stdout,"TASK %i incorrect termination!\n",task->pid);
    if(WIFSIGNALED(status)) task->status=KILLED;/*<-necessary?*/
#if DEBUG_TASK_PROCESS
  /*do we have checkpoint/restart?*/
    if(WIFSTOPPED(status))
  	fprintf(stdout,"TASK %i stopped!\n",task->pid);
    else if(WIFCONTINUED(status))
  	fprintf(stdout,"TASK %i continued!\n",task->pid);
#endif
  } while((!WIFEXITED(status))&&(!WIFSIGNALED(status)));
}

/* NB: GUI updates should ALL be done in the cleanup task, since we */
/* do the threads enter to avoid problems - this locks the GUI */
/* until we call the threads leave function at the end */

gdk_threads_enter();

if (task->status != KILLED)
  {
/* execute the cleanup task */
  if (task->cleanup)
    task->cleanup(task->ptr2);
  task->status = COMPLETED;
  }

gdk_flush();
gdk_threads_leave();

/* job completion notification */
gdk_beep();
}

/***************************************************/
/* set up the thread pool to process task requests */
/***************************************************/
void task_queue_init(void)
{
#ifdef G_THREADS_ENABLED
gdk_threads_init();
if (!g_thread_supported())
  {
/* TODO - disallow queueing of background tasks if this happens */
  gui_text_show(ERROR, "Task queue initialization failed.\n");
  }
else
  sysenv.thread_pool = g_thread_pool_new((GFunc) task_process, NULL,
                                         sysenv.max_threads,
                                         FALSE, NULL);
#endif
}

/*****************************/
/* terminate the thread pool */
/*****************************/
void task_queue_free(void)
{
g_thread_pool_free(sysenv.thread_pool, TRUE, FALSE);
}

/*************************/
/* free a task structure */
/*************************/
void task_free(gpointer data)
{
struct task_pak *task = data;

g_assert(task != NULL);

g_free(task->label);
g_free(task->time);
g_free(task->message);
g_free(task->status_file);

if (task->status_fp)
  fclose(task->status_fp);

g_string_free(task->status_text, TRUE);

g_free(task);
}

/****************************/
/* submit a background task */
/****************************/
/* TODO - only show certain tasks in the manager, since this */
/* could be used to do any tasks in the background - some of */
/* which may be slow GUI tasks we dont want to be cancellable */
void task_new(const gchar *label,
              gpointer func1, gpointer arg1,
              gpointer func2, gpointer arg2,
              gpointer model)
{
struct task_pak *task;

/* duplicate the task data */
task = g_malloc(sizeof(struct task_pak));
sysenv.task_list = g_slist_prepend(sysenv.task_list, task);
task->is_async = FALSE;
task->pid = -1;
task->status = QUEUED;
task->time = NULL;
task->message = NULL;
task->pcpu = 0.0;
task->pmem = 0.0;
task->progress = 0.0;
task->locked_model = model;

task->status_file = NULL;
task->status_fp = NULL;
task->status_index = -1;
task->status_text = g_string_new(NULL);

task->label = g_strdup(label);
task->primary = func1;
task->cleanup = func2;
task->ptr1 = arg1;
task->ptr2 = arg2;
/*
if (model)
  ((struct model_pak *) model)->locked = TRUE;
*/

/* queue the task */
g_thread_pool_push(sysenv.thread_pool, task, NULL);
}

/**************************************/
/* submit a system task:              */
/* tentative to get some tasks hidden */
/* and executed immediately.. (OVHPA) */
/**************************************/
void task_system_new(const gchar *label,
              gpointer func1, gpointer arg1,
              gpointer func2, gpointer arg2,
              gpointer model)
{
struct task_pak *task;

/* duplicate the task data */ 
task = g_malloc(sizeof(struct task_pak));
sysenv.task_list = g_slist_prepend(sysenv.task_list, task);
task->is_async = FALSE;
task->pid = -1;
task->status = QUEUED;
task->time = NULL;
task->message = NULL;
task->pcpu = 0.0;
task->pmem = 0.0;
task->progress = 0.0;
task->locked_model = model;

task->status_file = NULL;
task->status_fp = NULL;
task->status_index = -1;
task->status_text = g_string_new(NULL);

task->label = g_strdup_printf("SYS: %s",label);
task->primary = func1;
task->cleanup = func2;
task->ptr1 = arg1;
task->ptr2 = arg2;

/* queue the task */
g_thread_pool_push(sysenv.thread_pool, task, NULL);
	/*push it to front*/
	if(g_thread_pool_move_to_front(sysenv.thread_pool,task)){
		/*task is selected and front*/
		if(task->status==QUEUED){
			/*temporarily increase the pool size to immediately get the system task to start*/
			gint max=g_thread_pool_get_max_threads(sysenv.thread_pool);
			g_thread_pool_set_max_threads(sysenv.thread_pool,max+2,NULL);
			while(task->status==QUEUED) usleep(50*1000);/*sleep 50ms until sys task is started*/
			g_thread_pool_set_max_threads(sysenv.thread_pool,max,NULL);/*go back to former thread pool size*/
		}
	}
}
/*********************************************/
/* NEW - async call gives correct thread PID */
/*********************************************/
#define DEBUG_TASK_ASYNC 0
gint task_async(const gchar *command,pid_t *pid)
{
gint status;
gchar **argv;
GError *error=NULL;

/* checks */
if (!command)
  return(1);

#if _WIN32
chdir(sysenv.cwd);
system(command);
#else
/* setup the command vector */
argv = g_malloc(4 * sizeof(gchar *));
*(argv) = g_strdup("/bin/sh");
*(argv+1) = g_strdup("-c");
*(argv+2) = g_strdup(command);
*(argv+3) = NULL;

status=g_spawn_async(sysenv.cwd,argv,NULL,G_SPAWN_DO_NOT_REAP_CHILD,NULL,NULL,pid,&error);

g_strfreev(argv);
#endif

if (!status)
  printf("task_async() error: %s\n", error->message);

return(status);
}


/**************************************/
/* platform independant task spawning */
/**************************************/
#define DEBUG_TASK_SYNC 0
gint task_sync(const gchar *command) 
{
gint status;
gchar **argv;
GError *error=NULL;

/* checks */
if (!command)
  return(1);

#if _WIN32
chdir(sysenv.cwd);
system(command);
#else
/* setup the command vector */
argv = g_malloc(4 * sizeof(gchar *));
*(argv) = g_strdup("/bin/sh");
*(argv+1) = g_strdup("-c");
*(argv+2) = g_strdup(command);
*(argv+3) = NULL;
status = g_spawn_sync(sysenv.cwd, argv, NULL, 0, NULL, NULL, NULL, NULL, NULL, &error);
g_strfreev(argv);
#endif

if (!status)
  printf("task_sync() error: %s\n", error->message);

return(status);
}

/********************************************/
/* filter out unwanted lines in status file */
/********************************************/
gint task_status_keep(gint type, const gchar *line)
{
switch (type)
  {
  case GULP:
    if (strstr(line, "CPU"))
      return(1);
    if (strstr(line, " **"))
      return(1);
/*
    if (strstr(line, "="))
      if (strstr(line, "energy"))
        return(1);
*/
    break;

  default:
    return(1);
  }
return(0);
}

/**************************************************/
/* create descriptive string from the status file */
/**************************************************/
void task_status_update(struct task_pak *task)
{
/*gint filter;*/
gchar *line;

g_assert(task != NULL);

/* setup any status file filtering */
/*
if (g_ascii_strncasecmp("gulp", task->label, 4) == 0)
  filter = GULP;
else
  filter = 0;
*/

/* read in the status file */
if (task->status_file)
  {
  if (!task->status_fp)
    {
    task->status_index = 0;
/* exit if we've read in the file and closed it (due to completion) */
    if (strlen((task->status_text)->str))
      return;
    task->status_fp = fopen(task->status_file, "rt");
    }

  line = file_read_line(task->status_fp);
  while (line)
    {
/*
    if (task_status_keep(filter, line))
*/
      g_string_append(task->status_text, line);

    g_free(line);
    line = file_read_line(task->status_fp);
    }

  if (task->status == COMPLETED || task->status == KILLED)
    {
    fclose(task->status_fp);
    task->status_fp = NULL;
    }
  }
}

