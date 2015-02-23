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
#include "scan.h"
#include "host.h"

/* top level data structure */
extern struct sysenv_pak sysenv;

/* NB: this host stuff could equally well be the localhost, */
/* instead of some remote host ... in which case the system() */
/* call could be used instead of piping commands via the input */
/* and output file descriptors */

/* NB: remote stuff, where we issue commands - assumes unix */
/* style host ie echo, date, cd ... etc are all assumed to work */

struct host_pak 
{
gchar *name;
gchar *cwd;
gint type;         /* unix only at the moment */
gint connected;
gint input;
gint output;

/* NEW */
/* if qsub & mpirun services available - permits alternate invocation */
gpointer services;
};

struct service_pak 
{
gpointer host;
gint flags;
gchar *fullpath;
};

/**************************/
/* line reading primitive */
/**************************/
int readstr(int fd, char *a, int n)
{
int i;
	
for (i=0; i<n; i++)
  {
  if (read(fd, a+i, 1) != 1)
    return -1;

  if (a[i] == '\n')
    {
    a[i] = 0;
    return i;
    }
  }
return n;
}

/**********************/
/* free a service pak */
/**********************/
void host_service_free(gpointer data)
{
struct service_pak *service = data;

g_free(service->fullpath);
g_free(service);
}

/**********************************/
/* allocate a new host connection */
/**********************************/
gpointer host_new(gchar *name)
{
struct host_pak *host;

host = g_malloc(sizeof(struct host_pak));

host->name = g_strdup(name);
host->cwd = NULL;
host->connected = FALSE;
host->input = -1;
host->output = -1;
host->services = g_hash_table_new_full(g_str_hash, g_str_equal, g_free, host_service_free);

return(host);
}

/****************************************/
/* disconnect and free a host structure */
/****************************************/
void host_free(gpointer data)
{
struct host_pak *host = data;

if (host)
  {
  if (host->connected)
    host_disconnect(host);

  g_free(host->name);
  g_free(host->cwd);
  g_hash_table_destroy(host->services);
  g_free(host);
  }
}

/*************************/
/* free all host entries */
/*************************/
void host_free_all(void)
{
GSList *list;

list = sysenv.host_list;
while (list)
  {
  host_free(list->data);
  list = g_slist_next(list);
  }
g_slist_free(sysenv.host_list);
}

/*************************************************/
/* perform a foreach on the host's program table */
/*************************************************/
void host_service_foreach(gpointer data, gpointer func, gpointer arg)
{
struct host_pak *host = data;

g_assert(host != NULL);

g_hash_table_foreach(host->services, func, arg);
}

/*******************************/
/* structure access primitives */
/*******************************/
const gchar *host_name(gpointer data)
{
struct host_pak *host = data;

g_assert(host != NULL);

return(host->name);
}

const gchar *host_cwd(gpointer data)
{
struct host_pak *host = data;

g_assert(host != NULL);

return(host->cwd);
}

gpointer host_service_get(gpointer data, const gchar *name)
{
struct host_pak *host = data;

g_assert(host != NULL);

return(g_hash_table_lookup(host->services, name));
}

gchar *host_service_fullpath(gpointer data)
{
struct service_pak *service = data;

return(service->fullpath);
}

gint host_service_flags(gpointer data)
{
struct service_pak *service = data;

return(service->flags);
}

/*****************************************/
/* cycle through service execution types */
/*****************************************/
void host_service_cycle(gpointer data, const gchar *name)
{
struct host_pak *host = data;
struct service_pak *service;

service = g_hash_table_lookup(host->services, name);
if (!service)
  return;

/* don't cycle support services (eg qsub/mpi) */
if (service->flags == SERVICE_SECONDARY)
  return;

/* cycle the service type */
service->flags++;
if (service->flags >= SERVICE_PRIMARY)
  service->flags = SERVICE_BACKGROUND;
}

/**********************************************************************************/
/* check if a service name (and current flags) is available for the supplied host */
/**********************************************************************************/
gint host_service_available(gpointer data)
{
struct service_pak *service = data, *mpi, *qsub;
struct host_pak *host;

g_assert(service != NULL);
if (!service->fullpath)
  return(FALSE);

host = service->host;
mpi = g_hash_table_lookup(host->services, "mpirun");
qsub = g_hash_table_lookup(host->services, "qsub");

switch (service->flags)
  {
  case SERVICE_MPI:
    if (!mpi->fullpath)
      return(FALSE);
    break;
  case SERVICE_QSUB:
    if (!qsub->fullpath)
      return(FALSE);
    break;
  case SERVICE_QSUB_MPI:
    if (!mpi->fullpath || !qsub->fullpath)
      return(FALSE);
    break;
  }

return(TRUE);
}

/***************************************************/
/* primitive for populating the host program table */
/***************************************************/
/* NB: like all, assumes unix host */
void host_service_find(struct host_pak *host, const gchar *name)
{
gint flags;
gchar *text, *reply;
struct service_pak *service;

g_assert(host != NULL);
g_assert(name != NULL);

/* sepcial case services */
if (g_strrstr(name, "mpi") || g_strrstr(name, "qsub"))
  flags = SERVICE_SECONDARY;
else
  flags = SERVICE_BACKGROUND;

/* attempt to locate executable */
text = g_strdup_printf("which %s\n", name);
reply = host_talk(host, text);
g_free(text);

if (g_strrstr(reply, "not found"))
  {
  g_free(reply);
  reply = NULL;
  }

/* add the service */
service = g_malloc(sizeof(struct service_pak));
service->host = host;
service->flags = flags;
service->fullpath = reply;

g_hash_table_insert(host->services, g_strdup(name), service);
}

/**********************************************************/
/* setup session working directory for job file transfers */
/**********************************************************/
/* NB: assumes unix style remote host */
/* TODO - check for host type etc etc - which can be reported via the job/host gui */
void host_init(struct host_pak *host)
{
gchar *date, *text, *reply;

date = host_talk(host, "date\n");

/* force existence of gdis ssh session top level directory */
reply = host_talk(host, "md gdis_ssh\n");
g_free(reply);
reply = host_talk(host, "cd gdis_ssh\n");
g_free(reply);

/* create a working directory using current date */
text = g_strdup_printf("md \"%s\"\n", date);
reply = host_talk(host, text);
g_free(reply);
g_free(text);
text = g_strdup_printf("cd \"%s\"\n", date);
reply = host_talk(host, text);
g_free(reply);
g_free(text);


/* NEW */
host->cwd = host_talk(host, "pwd\n");


/* create a (dummy) control file */
/* TODO (maybe) ... could be used to save/restore info about submitted jobs */
reply = host_talk(host, "echo \"1\" > control.txt\n");
g_free(reply);

g_free(date);

/* initialize available host services */
host_service_find(host, "gulp");

/* secondary (enabling) services */
host_service_find(host, "qsub");
host_service_find(host, "mpirun");

/* register the initialized host */
sysenv.host_list = g_slist_prepend(sysenv.host_list, host);
}

/**********************************/
/* activate a host ssh connection */
/**********************************/
/* NB: assumes unix style remote host */
gint host_connect(gpointer data)
{
int n;
int fdto[2];
int fdfrom[2];
char cmd[1000];
struct host_pak *host = data;

g_assert(host != NULL);
if (!host->name)
  {
  perror("Fatal: no host name.");
  return(FALSE);
  }
if (pipe(fdto) == -1 || pipe(fdfrom) == -1) 
  {
  perror("Fatal: failed to create pipes.");
  return(FALSE);
  }

switch (fork()) 
  {
  case -1:
    perror("Fatal: failed to fork.");
    return(FALSE);

/* child */
  case 0:
/* setup child io streams */
    dup2(fdto[0], fileno(stdin));
    dup2(fdfrom[1], fileno(stdout));
/* close unwanted streams */
    close(fdto[1]); 
    close(fdfrom[0]);

/* FIXME - can we get rid of the stdin not a terminal complaint? */
/*
    execlp("ssh", "ssh", "-e", "none", "-T", "-q", host->name, (char *) 0);
*/

    execlp("ssh", "ssh", "-T", "-q", host->name, (char *) 0);

    exit(3);

/* parent */
  default:
/* close unwanted streams */
    close(fdto[0]);
    close(fdfrom[1]);

/* send a command */
    write(fdto[1], "echo \"ssh:connect\"\n", 19);

/* read/write to child (ssh) via file desc. */
/* read -> fdfrom[0] */
/* write -> fdto[1] */
    for (;;)
      {
      n = readstr(fdfrom[0], cmd, sizeof(cmd));

      if (n > 1)
        {
/* if we get here ... we've established a connection */
/*
    printf("%s\n", cmd);
*/
        if (strstr(cmd, "ssh:connect"))
          {
          printf("Connection established.\n");

          host->input = fdto[1];
          host->output = fdfrom[0];
          host->connected = TRUE;

/* set up session working directory */
          host_init(host);

          break;
          }
        }
      else
        {
        host->connected = FALSE;
        return(FALSE);
        }
      }
  }

return(TRUE);
}

/************************************/
/* deactivate a host ssh connection */
/************************************/
void host_disconnect(gpointer data)
{
struct host_pak *host = data;

printf("Closing connection to: %s\n", host->name);

host->connected = FALSE;
write(host->input, "exit\n", 5);
}

/************************************************************/
/* send a message to remote host shell, and return response */
/************************************************************/
/* NB: assumes unix style remote host */
gchar *host_talk(gpointer data, const gchar *message)
{
gint n;
char cmd[1000];
GString *response;
struct host_pak *host = data;

g_assert(host != NULL);

if (!host->connected)
  return(NULL); 

write(host->input, message, strlen(message));
write(host->input, "echo \"host:done\"\n", 17);

response = g_string_new(NULL);

/* TODO - implement some form of timeout response check? */
/* if error, either we didnt connect, or host is not unix */
/* and didnt understand the echo command */
for (;;)
  {
  n = readstr(host->output, cmd, sizeof(cmd));

  if (n > 1)
    {
    if (strstr(cmd, "host:done"))
      break;
    else
      g_string_sprintfa(response, "%s", cmd);
    }
  }

return(g_string_free(response, FALSE));
}

/*********************************************/
/* ensure input string is safe to print/echo */
/*********************************************/
/* NB: unix shell specific */
gchar *host_safe_text(const gchar *input)
{
gint i;
GString *output;

if (!input)
  return(NULL);

output = g_string_new(NULL);

for (i=0 ; i<strlen(input) ; i++)
  {
  switch (input[i])
    {
/* escape special chars */
    case '$':
    case '!':
      g_string_sprintfa(output, "\\%c", input[i]);
      break;

    default:
/* disallow control chars */
      if (g_ascii_iscntrl(input[i]))
        break;
      g_string_sprintfa(output, "%c", input[i]);
    }
  }

return(g_string_free(output, FALSE));
}

/*****************************************/
/* write a locally stored file to a host */
/*****************************************/
/* FIXME - this works ... but is ugly and slow ... the experimental version is better */
/* (ie using ascii-xfr) but doesnt terminate properly */
gint host_file_write(gpointer data, gchar *local, gchar *remote)
{
gchar *text, *tmp, *reply;
gpointer scan;
struct host_pak *host = data;

g_assert(host != NULL);

scan = scan_new(local);
if (!scan)
  {
  printf("Error: missing local file: %s\n", local);
  return(1);
  }

/* create the remote file */
text = g_strdup_printf("echo -n \"\" > \"%s\"\n", remote);

reply = host_talk(host, text);

g_free(reply);
g_free(text);

/* concat line by line */
/* FIXME - is there a better way to do this? */
for (;;)
  {
/* escape special characters ... remove control codes */
  tmp = host_safe_text(scan_get_line(scan));
  if (scan_complete(scan))
    break;

  text = g_strdup_printf("echo \"%s\" >> \"%s\"\n", tmp, remote);
  g_free(tmp);

  reply = host_talk(host, text);
  g_free(reply);
  g_free(text);
  }

return(0);
}

/*****************************************/
/* write a locally stored file to a host */
/*****************************************/
/* PROBLEM is - doesnt seem to terminate properly */
gint host_file_write_experimental(gpointer data, gchar *local, gchar *remote)
{
gchar *text;
struct host_pak *host = data;
gpointer scan;

if (!host->connected)
  return(FALSE); 

scan = scan_new(local);
if (!scan)
  return(1);

text = g_strdup_printf("ascii-xfr -v -r %s\n", remote);
write(host->input, text, strlen(text));
g_free(text);

for (;;)
  {
  text = scan_get_line(scan);
  if (scan_complete(scan))
    {
/* FIXME - does not seem to be terminating the transfer ...  */
    write(host->input, "", 1);
    write(host->input, "", 1);
    break;
    }
  write(host->input, text, strlen(text));
  }

scan_free(scan);

return(FALSE);
}

/************************************************/
/* read a file on a remote host to a local file */
/************************************************/
gint host_file_read(gpointer data, gchar *local, gchar *remote)
{
gchar *text, cmd[10];
struct host_pak *host = data;
FILE *fp;

if (!host->connected)
  return(FALSE); 

text = g_strdup_printf("ascii-xfr -v -s -e %s\n", remote);

write(host->input, text, strlen(text));

fp = fopen(local, "w");

for (;;)
  {
  read(host->output, cmd, 1);

fprintf(fp, "%c", cmd[0]);

  if (cmd[0] == 26 || cmd[0] == '')
    {
    break;
    }
  }

return(FALSE);
}

/****************************************************/
/* simple test of the host communication primitives */
/****************************************************/
/* NB: assumes unix style remote host */
void host_test(void)
{
gpointer host;
gchar *reply;

host = host_new("sean@onyx.ivec.org");

if (host_connect(host))
  printf("host_connect(): success!\n");
else
  {
  printf("host_connect(): failed!\n");
  return;
  }


{
/*
struct host_pak *h = host;

write(h->input, "cat > gok.txt\n", 14);
write(h->input, "abcdef\n", 7);
*/

/* CURRENT - how do we terminate the redirect ... ie simulate a control-D or control-C */
/* appears there is no portable way for doing it ... can send stuff to /dev/tty in linux */
/* none of these work ... */
/*
write(h->input, (char *) '~', 1);
write(h->input, (char *) 'd', 1);
write(h->input, (char *) "^C", 2);
write(h->input, (char *) "^D", 2);
write(h->input, (char *) "", 1);
write(h->input, (char *) "", 1);
write(h->input, (char *) 27, 1);
write(h->input, (char *) "", 1);
write(h->input, (char *) 16, 1);
write(h->input, (char *) "", 1);
*/

}

/*
printf("sending: [date]\n");
reply = host_talk(host, "date\n");
printf("response: [%s]\n", reply);
g_free(reply);
*/


printf("sending: [which gulp]\n");
reply = host_talk(host, "which gulp\n");
printf("response: [%s]\n", reply);
g_free(reply);

/* CURRENT - using ascii-xfr to accomplish these */
/* FIXME - read then write works ... but reversing order causes a hang ... why??? */
/* seems that after a write has been done ... things are messed up */

/*
host_file_read(host, "dummy.txt", "/home/sean/.cshrc");
host_file_write(host, "aloh4-.car", "aloh4-copy.car");
host_file_write(host, "aloh4-.car", "aloh4-copy2.car");
*/



host_disconnect(host);

host_free(host);
}
#endif
