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
#include <stdlib.h>
#include <unistd.h>
#include <sys/types.h>

#ifndef __WIN32
#include <sys/wait.h>
#endif

#include "gdis.h"
#include "file.h"
#include "grid.h"
#include "job.h"
#include "task.h"
#include "parse.h"
#include "model.h"
#include "grisu_client.h"
#include "interface.h"

extern struct sysenv_pak sysenv;

/* CURRENT - register apps for a GRID destination */
/* NULL -> not on grid */
/* TODO - consider storing structs containing more info */
/*        ie site, app details etc - for faster job building */

GHashTable *grid_table=NULL;
/* CURRENT - assume credential already uploaded for the time being */
gint grid_authenticated=FALSE;
gchar *myproxy_init=NULL;

void grid_application_set(const gchar *name, const gchar *value)
{
#ifdef WITH_GRISU
if (grid_table)
  g_hash_table_replace(grid_table, g_strdup(name), g_strdup(value));
#else
return;
#endif
}

gchar *grid_application_get(const gchar *name)
{
#ifdef WITH_GRISU
if (grid_table)
  return(g_hash_table_lookup(grid_table, name));
else
  return(NULL);
#else
return(NULL);
#endif
}

void grid_application_remove(const gchar *name)
{
#ifdef WITH_GRISU
if (grid_table)
  g_hash_table_remove(grid_table, name);
#else
return;
#endif
}

GList *grid_application_all(void)
{
#ifdef WITH_GRISU
return(g_hash_table_get_keys(grid_table));
#else
return(NULL);
#endif
}

GSList *grid_search_by_name(const gchar *name)
{
#ifdef WITH_GRISU
return(grisu_submit_find(name));
#else
return(NULL);
#endif
}

/*********************************************/
/* allocate and initialize a grid job object */
/*********************************************/
gpointer grid_new(void)
{
struct grid_pak *grid;

grid = g_malloc(sizeof(struct grid_pak));

grid->user_vo = NULL;
grid->jobname = NULL;
grid->exename = NULL;
grid->exe_version = NULL;
grid->jobcode = JOB_UNKNOWN;

grid->remote_q = NULL;
grid->remote_root = NULL;
grid->remote_exe = NULL;
grid->remote_exe_module = NULL;
grid->remote_site = NULL;

grid->remote_exe_type=-1;
grid->remote_exe_np=1;

grid->local_cwd = NULL;
grid->local_input = NULL;
grid->local_output = NULL;

return(grid);
}

/**************************/
/* free a grid job object */
/**************************/
void grid_free(gpointer data)
{
struct grid_pak *grid=data;

g_assert(grid != NULL);

g_free(grid->user_vo);
g_free(grid->jobname);
g_free(grid->exename);
g_free(grid->exe_version);

g_free(grid->remote_q);
g_free(grid->remote_root);
g_free(grid->remote_exe);
g_free(grid->remote_exe_module);
g_free(grid->remote_site);

g_free(grid->local_cwd);
g_free(grid->local_input);
g_free(grid->local_output);

g_free(grid);
}

/***************************************/
/* get the DN from an x509 certificate */
/***************************************/
// openssl x509 -noout -in usercert.pem -subject
gchar *grid_get_DN(gchar *certificate_fullpath)
{
gint i, j, status;
gchar *cmd, *out, *err, *subject=NULL;
GError *error=NULL;

cmd = g_strdup_printf("openssl x509 -noout -in %s -subject", certificate_fullpath);

/* TODO - get the output */
if (g_spawn_command_line_sync(cmd, &out, &err, &status, &error))
  {
  for (i=0 ; i<8000 ; i++)
    {
    if (out[i] == '/')
      {
/* not sure why but plain g_strdup adds a crappy char (null?) */
/* to the string, so scan for the end and use strndup instead */
      for (j=i ; j<9000 ; j++)
        {
        if (out[j] == '\0')
          break;
        }
      subject = g_strndup(&out[i], j-i-1);
      break;
      }
    }
  }
g_free(cmd);

return(subject);
}

gboolean grid_credential_get()
{
return(FALSE);
}

/*******************************************************/
/* set valid auth in case of external means (eg) grisu */
/*******************************************************/
/* TODO - could chuck into sysenv eventually */
void grid_auth_set(gint value)
{
grid_authenticated=value;
}

/*********************************************/
/* verify current authentication information */
/*********************************************/
gint grid_auth_check(void)
{
#ifdef WITH_GRISU
GSList *list;

/* perform a grisu operation that requires authentication */
/* NB: use anything fast (as long as auth info is in the header) */
//list = grisu_fqans_get();
list = grisu_site_list();
if (list)
  {
  g_slist_free(list);
  grid_authenticated = TRUE;
  return(TRUE);
  }
#endif

grid_authenticated = FALSE;

return(FALSE);
}

/****************************************************/
/* upload a temporary credential to a remote server */
/****************************************************/
/* TODO - if grid_passwd is NULL -> use stdin */
#define DEBUG_GRID_CREDENTIAL_INIT 0
void grid_credential_init(const gchar *grid_password)
{
#ifdef WITH_GRISU
gint inp, out, err;
gchar **argv, **envp;
const gchar *username, *password, *server, *port;
GError *error=NULL;
GPid pid;
FILE *fpi, *fpo, *fpe;

/* TODO - informative user error messages */
if (!grid_password)
  return;
if (grid_authenticated)
  return;
if (!myproxy_init)
  return;

username = grisu_username_get();
password = grisu_password_get();
server = grisu_myproxy_get();
port = grisu_myproxyport_get();

#if DEBUG_GRID_CREDENTIAL_INIT
printf("Credential for: "); 
printf("%s : %s -> %s : %s\n", username, password, server, port);
#endif

/* FIXME - authentication issues with myproxy/myproxy2 */
/* solution is to point at myproxy2 and set */
/* MYPROXY_SERVER_DN /C=AU/O=APACGrid/OU=VPAC/CN=myproxy2.arcs.org.au */

/* FIXME - need to ensure enviroment has GT_PROXY_MODE set to "old" */
/* unless markus fixes grisu's credential retrieval */
/* build argument vector */

/* FIXME - implement this */
envp = g_malloc(3 * sizeof(gchar *));
*(envp) = g_strdup("GT_PROXYMODE=old");
*(envp+1) = g_strdup("MYPROXY_SERVER_DN=/C=AU/O=APACGrid/OU=VPAC/CN=myproxy2.arcs.org.au");
*(envp+2) = NULL;


argv = g_malloc(10 * sizeof(gchar *));
/* NB: need full path to executable */
*(argv) = g_strdup(myproxy_init);

*(argv+1) = g_strdup("-l");
*(argv+2) = g_strdup(username);

/* small lifetime while we debug */
*(argv+3) = g_strdup("-c");
*(argv+4) = g_strdup("1");

*(argv+5) = g_strdup("-a");

*(argv+6) = g_strdup("-S");

*(argv+7) = g_strdup("-s");
*(argv+8) = g_strdup(server);
*(argv+9) = NULL;


/* spawn process */
if (g_spawn_async_with_pipes(NULL, argv, NULL, G_SPAWN_DO_NOT_REAP_CHILD,
                             NULL, NULL, &pid, &inp, &out, &err, &error))
  {

  if (pid)
    {
/* parent - wait until child is done */
//printf("parent: waiting for child...\n");

fpi = fdopen(inp, "a");
fpo = fdopen(out, "r");
fpe = fdopen(err, "r");

/* NB: -S will make myproxy prompt for grid passphrase */
/* and then the credential password once only */
/* CURRENT - not sure why this isnt working - it seems */
/* to work if run manually from the command line */

fprintf(fpi, "%s\n", grid_password);
fflush(fpi);

/* CURRENT - this seems to fix the problem! */
/* obviously a delay is needed while myproxy-init processes the */
/* grid_password to generate a local proxy */
/* the problem is - how to we poll stdout to know when to send the */
/* passphrase ... which is probably better than an arbitrary sleep */

#ifndef WIN32
sleep(5);
#endif

/* test if process has terminated */
/* note that if it has - child is reaped and we cant call waitpid() again */

#ifndef __WIN32
if (waitpid(pid, NULL, WNOHANG) > 0)
  {
printf("Bad password.\n");

    g_spawn_close_pid(pid);
    grid_authenticated = FALSE;

  return;
  }
#endif

//printf("waitpid 4: %d\n", waitpid(pid, NULL, WEXITED | WNOHANG));


fprintf(fpi, "%s\n", password);
fflush(fpi);
fclose(fpi);

#if DEBUG_GRID_CREDENTIAL_INIT
gchar *line;
do
  {
  line = file_read_line(fpo);
  printf("FORK, stdout: %s\n", line);
  g_free(line);
  }
while(line != NULL);
do
  {
  line = file_read_line(fpe);
  printf("FORK, stderr : %s\n", line);
  g_free(line);
  }
while(line != NULL);
#endif

/* cleanup */
fclose(fpo);
fclose(fpe);

close(inp);
close(out);
close(err);

#ifndef __WIN32
    waitpid(pid, NULL, WEXITED);
#endif

    g_spawn_close_pid(pid);

    grid_authenticated = TRUE;


/* sleep - give the uploaded proxy some time to "stick" to the server */
/* FIXME - better to poll (myproxy-logon?) for the credential */
printf("Uploading Credential, please wait ...\n");

#ifndef WIN32
sleep(5);
#endif

printf("Done.\n");

    }
  }
else
  {
  printf("Failed to upload credential, myproxy-init not found or bad connectivity?\n");
  }

g_strfreev(argv);
g_strfreev(envp);
#endif

return;
}

/****************************************************/
/* get short text message describing current status */
/****************************************************/
const gchar *grid_get_status()
{
if (grid_authenticated)
  return("Authenticated");
else
  return("Not authenticated");
}

/***********************************************************/
/* generate a random alphabetical string of specified size */
/***********************************************************/
gchar *grid_random_alpha(gint length)
{
gint i;
GRand *r;
GString *text;

r = g_rand_new();
text = g_string_new(NULL);

for (i=length ; i-- ; )
  {
  g_string_sprintfa(text, "%c", g_rand_int_range(r, 'a', 'z'));
  }

g_rand_free(r);

return(g_string_free(text, FALSE));
}

/************************************************************/
/* generate a random alpha-numeric string of specified size */
/************************************************************/
gchar *grid_random_alphanum(gint length)
{
gint i, decide;
GRand *r;
GString *text;

r = g_rand_new();
text = g_string_new(NULL);

/* NB: some systems require an alphabetical character first */
g_string_sprintfa(text, "%c", g_rand_int_range(r, 'a', 'z'));
i=1;
while (i<length)
  {
  decide = g_rand_int_range(r, 0, 99);
  if (decide < 50)
    g_string_sprintfa(text, "%c", g_rand_int_range(r, 'a', 'z'));
  else
    g_string_sprintfa(text, "%c", g_rand_int_range(r, '0', '9'));
  i++;
  }

g_rand_free(r);

return(g_string_free(text, FALSE));
}

/**********************************/
/* attempt to upload a credential */
/**********************************/
gint grid_connect(gint method)
{
#ifdef WITH_GRISU
switch (method)
  {
  case GRID_MYPROXY:
    grid_auth_check();
    break;

  case GRID_LOCALPROXY:
/* if successful, send info to grisu client for soap header */
/* CURRENT - myproxy-init fed with grid_passwd */
    grid_credential_init(grisu_password_get());
    break;
  }

grisu_auth_set(TRUE);

return(TRUE);

#else
return(FALSE);
#endif
}

/**************/
/* initialize */
/**************/
void grid_setup(void)
{
if (!grid_table)
  {
#ifdef WITH_GRISU
  grisu_init();

  if (!myproxy_init)
    myproxy_init = g_find_program_in_path("myproxy-init");
  if (!myproxy_init)
    printf("grid_setup() Warning. Failed to locate: myproxy-init\n");

#endif

  grid_table = g_hash_table_new(g_str_hash, g_str_equal);
  }
}

/***********/
/* cleanup */
/***********/
void grid_cleanup(void)
{
g_hash_table_destroy(grid_table);
}

/*****************************************/
/* process completed jobs after download */
/*****************************************/
/* this could be generic, since we assume output is on local filestore */
void grid_process_job(const gchar *job_name, GSList *file_list)
{
#ifdef WITH_GRISU
gint type;
gchar *ext, *name;
gchar *input=NULL, *output=NULL;
GSList *item;
struct model_pak *model;

type = grisu_job_type(job_name);

switch (type)
  {
  case JOB_GULP:
    printf("processing GULP job ... \n");

    for (item=file_list ; item ; item=g_slist_next(item))
      {
      name = item->data;

/* skip if no extension */
      ext = file_extension_get(name);
      if (!ext)
        continue;

/* load results file as preference for coords */
      if (g_strncasecmp(ext, "res", 3) == 0)
        input = name;
/* only load original input file if no results file found */
      if (g_strncasecmp(ext, "gin", 3) == 0)
        {
        if (!input)
          input = name;
        }
/* output file containing energies etc */
      if (g_strncasecmp(ext, "go", 2) == 0)
        output = name;
      }

if (!output)
  {
  printf("No GULP output file found. Failed job?\n");
  break;
  }
if (!input)
  {
  printf("No GULP input file found. Failed authentication?\n");
  break;
  }

/* create new model for the results */
    model = model_new();
/* read main data from the res file (correct charges etc.) */
    read_gulp(input, model);
/* graft to the model tree, so subsequent GULP read doesn't replace coords */

/* CURRENT - give the loaded model the same name as the job ... better method? */
    if (model->basename)
      {
      g_free(model->basename);
//      model->basename = g_strdup(job_name);
//      model->basename = g_strdup(input);

      model->basename = parse_strip(input);
      }

sysenv.active_model = model;

    tree_model_add(model);
/* get the output energy/phonons etc. */
    read_gulp_output(output, model);

/* FIXME - if the GULP job fails - model_prep() isnt called, and the */
/* model camera isnt initialized => error trap when gdis tries to visualize */
    if (!model->camera)
      {
printf("WARNING: GULP calculation has possibly failed.\n");
      model_prep(model);
      }
/* TODO - return model - so gui part decides whether to force a select/display or not */
//gui_model_select(model);

    break;

  case JOB_GAMESS:
    printf("processsing GAMESS job ... \n");

    for (item=file_list ; item ; item=g_slist_next(item))
      {
      name = item->data;
      ext = file_extension_get(name);
/* skip if no extension */
      if (!ext)
        continue;

      if (g_strncasecmp(ext, "gmot", 4) == 0)
        output = name;
      }

    if (output)
      {
/* create new model for the results */
      model = model_new();
/* read main data from the res file (correct charges etc.) */
      read_gms_out(output, model);

      if (model->basename)
        {
        g_free(model->basename);
        model->basename = g_strdup(job_name);
        }

sysenv.active_model = model;

      tree_model_add(model);

/* FIXME - if the job fails - model_prep() isnt called, and the */
/* model camera isnt initialized => error trap when gdis tries to visualize */
    if (!model->camera)
      {
printf("WARNING: GAMESS calculation has possibly failed.\n");
      model_prep(model);
      }

/* TODO - return model - so gui part decides whether to force a select/display or not */
//gui_model_select(model);
//sysenv.active_model = model;

      }
    else
      {
printf("No GAMESS output file found. Failed job?\n");
      }

    break;

  default:
    printf("Error, unknown job\n");
    return;
  }
#endif
}

/***************************/
/* download completed jobs */
/***************************/
GSList *grid_download_job(const gchar *name)
{
#ifdef WITH_GRISU
gchar *remote_cwd, *local, *remote;
GSList *item, *list, *results=NULL;

printf("Getting result files for: %s\n", name);

remote_cwd = grisu_job_cwd(name);

/* TODO - get job type (GULP etc) from name (eg grisu_job_type()) and filter */
//printf("Working directory = %s\n", cwd);

if (remote_cwd)
  {
  list = grisu_file_list(remote_cwd);

  for (item=list ; item ; item=g_slist_next(item))
    {
    gchar *file = item->data;
    printf("[%s]\n", file);
// TODO - strictly - should be the job's local_output destination here */
    local =  g_build_filename(sysenv.cwd, file, NULL);
    remote =  g_build_filename(remote_cwd, file, NULL);

//printf("Transferring [%s] -> [%s]\n", remote, local);

/* if download is successful add to results list for post-processing */
    if (grisu_file_download(remote, local))
      results = g_slist_prepend(results, local);

    g_free(remote);
    }
  g_free(remote_cwd);
  }

return(results);
#else
return(NULL);
#endif
}

/**********************************/
/* background task job submission */
/**********************************/
void grid_job_start(gpointer data)
{
#ifdef WITH_GRISU
gint alen, rlen;
gchar *jobxml, *value, *tmp;
gchar *fs_absolute, *fs_relative;
gchar *local, *remote;
GHashTable *details;
GSList *list;
struct grid_pak *grid=data;

g_assert(grid != NULL);
g_assert(grid->exename != NULL);
g_assert(grid->remote_q != NULL);

grid->jobcode = grisu_job_type(grid->exename);
grid->user_vo = g_strdup(grisu_vo_get());

grid->remote_site = grisu_site_name(grid->remote_q);
if (!grid->remote_site)
  {
  printf("Could not determine submission site.\n");
  return;
  }

list = grisu_application_versions(grid->exename, grid->remote_site);
if (list)
  {
/* CURRENT - just use the first entry */
  if (grid->exe_version)
    g_free(grid->exe_version);
  grid->exe_version = g_strdup(list->data); 
  }
else
  {
/* strictly, not all programs would require a version */
/* eg gamess does, but gulp doesnt */
  printf("Warning: could not determine executable version of [%s] at [%s]\n", grid->exename, grid->remote_site);
  }

details = grisu_application_details(grid->exename, grid->remote_site);
if (details)
  {
  tmp = g_hash_table_lookup(details, "Executables");
  if (!tmp)
    {
    printf("Failed to locate executable on remote site.\n");
    return;
    }
  grid->remote_exe = g_strdup(tmp);

  tmp = g_hash_table_lookup(details, "Module");
  if (tmp)
    grid->remote_exe_module = g_strdup(tmp);

//printf("Requiring module [%s]\n", grid->remote_exe_module);

/* CURRENT - prefer serial over parallel */
  value = g_hash_table_lookup(details, "ParallelAvail");
  if (g_strncasecmp(value, "true", 4) == 0)
    {
    grid->remote_exe_type = GRID_MPI;
/* CURRENT - hack to force // ie mpi */
    grid->remote_exe_np = 2;
    }

/* CURRENT - prefer serial over parallel */
  value = g_hash_table_lookup(details, "SerialAvail");
  if (g_strncasecmp(value, "true", 4) == 0)
    {
    grid->remote_exe_type = GRID_SERIAL;
    grid->remote_exe_np = 1;
    }

  if (grid->remote_exe_type == -1)
    {
    printf("Failed to determine serial/parallel invocation method on remote site.\n");
    return;
    }
  }
else
  {
  printf("Failed to get application details from MDS.\n");
  return;
  }

// -----------------

grid->jobname = grisu_jobname_request(grid->jobcode);
if (!grid->jobname)
  {
  printf("Failed to create job, not authenticated?\n");
  return;
  }

/* calculate mount point */
fs_absolute = grisu_absolute_job_dir(grid->jobname, grid->remote_q, grid->user_vo);
if (!fs_absolute)
  {
  printf("Failed to locate a valid remote filesystem for VO [%s] at [%s]\n", grid->user_vo, grid->remote_site);
  return;
  }

fs_relative = grisu_relative_job_dir(grid->jobname);
g_assert(fs_relative != NULL);

alen = strlen(fs_absolute);
rlen = strlen(fs_relative);
if (alen > rlen)
  grid->remote_root = g_strndup(fs_absolute, alen-rlen);
else
  {
  printf("Bad filesystem returned by grisu.\n");
// cope as best we can
  grid->remote_root = g_strdup(fs_absolute);
  }

//printf("Remote root: %s\n", grid->remote_root);

/* upload stage */

/* TODO - stat the file to ensure it was written */
if (!grid->local_input || !grid->local_cwd)
  {
  printf("Error, job input file not correctly written.\n");
  return;
  }

/* TODO - could do this earlier as well ... */
switch (grid->jobcode)
  {
  case JOB_GAMESS:
    grid->local_output = parse_extension_set(grid->local_input, "gmot");
    break;
  case JOB_GULP:
    grid->local_output = parse_extension_set(grid->local_input, "got");
    break;
  }

local = g_build_filename(grid->local_cwd, grid->local_input, NULL);
remote = g_build_filename(fs_absolute, grid->local_input, NULL);

//printf("upload [%s] -> [%s]\n", local, remote);

grisu_file_upload(local, remote);

g_free(local);
g_free(remote);


jobxml = grisu_xml_create(grid);

//printf(" *** XML:\n%s\n", jobxml);

grisu_job_configure(grid->jobname, jobxml);
grisu_job_submit(grid->jobname, grid->user_vo);

g_free(fs_absolute);
g_free(fs_relative);
g_free(jobxml);

#endif
}

/*******************/
/* completion task */
/*******************/
void grid_job_stop(gpointer data)
{
struct grid_pak *grid=data;

/* TODO - free the grid structre */
grid_free(grid);

#ifdef WITH_GUI
grid_task_stop();
#endif

/* refresh all jobs */
/* FIXME - too expensive to keep doing all the time */
/* FIXME - do this at start (eg when connecting ... and only add/del subsequently) */
//gui_job_refresh_all();
}

/*******************************************************/
/* submit a grid job (potentially from setup a dialog) */
/*******************************************************/
/* name = application name (eg "gulp") */
gint grid_job_submit(const gchar *name, gpointer data)
{
#ifdef WITH_GRISU
gchar *remote_q;
struct grid_pak *grid=data;

/* checks */
g_assert(grid != NULL);

/* test we have a valid grid destination */
remote_q = grid_application_get(name);
if (!remote_q)
  {
  printf("No valid destination queue registered for: %s\n", name);
  return(FALSE);
  }

#ifdef WITH_GUI
if (!grid_task_start())
  return(FALSE);
#endif

grid->exename = g_strdup(name);
grid->remote_q = g_strdup(remote_q);

/* submit and lock model */
task_new("grid", grid_job_start, grid, grid_job_stop, grid, NULL);

return(TRUE);
#else
printf("No grid support installed.\n");

return(FALSE);
#endif
}

