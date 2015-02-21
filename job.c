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
#include <unistd.h>

#include "gdis.h"
#include "file.h"
#include "host.h"
#include "job.h"

/* top level data structure */
extern struct sysenv_pak sysenv;

/* the essential aim of the job structure is simply to establish a link between 
input/ouput file pairs (local <-> remote) so that jobs can be tracked/updated.
NB: Place each remotely executed job in its own directory ( job1 job2 etc ) and give
the user an option to copy ALL files in that directory back to the local machine 
(by default copy only the output file) This way we get addition data
such as e density and anything else copied back if desired ... possible refinements
include a list of all files with size and type that user can check/uncheck as 
transaction transfers ... binary files may be a problem */

struct job_pak 
{
gint id;
gint type;
gint count;
gint status;

gchar *local_dir;
gchar *local_input;
gchar *local_output;

gchar *remote_dir;
gchar *remote_input;
gchar *remote_output;

gpointer host;
};

/**********************/
/* free job structure */
/**********************/
void job_free(struct job_pak *job)
{
g_free(job->local_dir);
g_free(job->local_input);
g_free(job->local_output);
g_free(job->remote_dir);
g_free(job->remote_input);
g_free(job->remote_output);
g_free(job);
}

/******************************/
/* initialize a job structure */
/******************************/
gpointer job_setup(const gchar *name, struct model_pak *model)
{
static gint count=0;
struct job_pak *job;

job = g_malloc(sizeof(struct job_pak));

/* use this to give unique job names (within cwd of an ssh session) */
/* since the cwd is date/time based we wont get overwrites */
job->count = count++;
job->host = NULL;
job->local_dir = NULL;
job->local_input = NULL;
job->local_output = NULL;
job->remote_dir = NULL;
job->remote_input = NULL;
job->remote_output = NULL;

if (g_ascii_strncasecmp("gulp", name, 4) == 0)
  job->type = JOB_GULP;

switch (job->type)
  {
  case JOB_GULP:
    job->local_dir = g_strdup(sysenv.cwd);
    job->local_input = g_strdup(model->gulp.temp_file);
    job->local_output = g_strdup(model->gulp.out_file);

/* FIXME - filename is not the full path name ... which strictly it should be */
/* however - job->local_input is expected to contain only the filename (without */
/* the path) by the subsequent job_start() call */
write_gulp(job->local_input, model);

    break;

  default:
    printf("Error: unknown job type: %s\n", name);
    g_free(job);
    job = NULL;
  }

return(job);
}

/********************************/
/* begin the execution of a job */
/********************************/
#define DEBUG_JOB_START 1
void job_start(gpointer host, gpointer data)
{
gchar *tmp;
struct job_pak *job = data;

/* assign host to job, and initialize */
job->host = host;
job->remote_dir = g_strdup(host_cwd(host));

/* build full path remote names */
job->remote_input = g_build_filename(job->remote_dir, job->local_input, NULL);
job->remote_output = g_build_filename(job->remote_dir, job->local_output, NULL);

/* convert local names to full path names */
tmp = g_build_filename(job->local_dir, job->local_input, NULL);
g_free(job->local_input);
job->local_input = tmp;

tmp = g_build_filename(job->local_dir, job->local_output, NULL);
g_free(job->local_output);
job->local_output = tmp;

#if DEBUG_JOB_START
printf("local dir: %s\n", job->local_dir);
printf("local inp: %s\n", job->local_input);
printf("local out: %s\n", job->local_output);
printf("remote dir: %s\n", job->remote_dir);
printf("remote inp: %s\n", job->remote_input);
printf("remote out: %s\n", job->remote_output);
#endif

host_file_write(host, job->local_input, job->remote_input);


/* TODO - execute ... bg or queued ... on remote host */


}

/***************************/
/* create a new job object */
/***************************/
gint job_new(const gchar *name, struct model_pak *model)
{
gpointer host, service, job;

if (!sysenv.host_list)
  {
  printf("No active connections.\n");
  return(FALSE);
  }
host = sysenv.host_list->data;

service = host_service_get(host, name);
if (!service)
  {
  printf("Unknown service: %s.\n", name);
  return(FALSE);
  }

if (!host_service_available(service))
  {
  printf("Service: %s not available on active host.\n", name);
  return(FALSE);
  }

printf("Requesting service: %s on host %s\n", name, host_name(host));

job = job_setup(name, model);

/* TODO - put elsewhere (some sort of scheduler decides the host part?) */
job_start(host, job);


/* TODO - this'll get placed in job cleanup when jobs are actually run */
job_free(job);


return(TRUE);
}

/********************************/
/* check on the status of a job */
/********************************/
void job_status_get(gpointer data)
{
struct job_pak *job = data;

g_assert(job != NULL);

/* TODO - more than one way of doing this ... eg could have qstat for queued jobs */
/* etc ... scanning output is more portable ... but wont determine if a job has crashed */

/* could record the last time the output file had new stuff added ... */

/* get output file */
host_file_read(job->host, job->local_output, job->remote_output);

/* scan output */


}

/************************************/
/* eg use to kill/stop/restart jobs */
/************************************/
void job_status_set()
{
}

#endif

