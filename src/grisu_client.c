/*
Copyright (C) 2008 by Sean David Fleming

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

#include <glib.h>
#include <time.h>
#include <sys/stat.h>

#include "gdis.h"
#include "grid.h"
#include "job.h"
#include "parse.h"
#include "logging.h"

#include "soapH.h"
#include "grisuHttpBinding.nsmap"

#define DEBUG_GRISU_CLIENT 1

extern struct sysenv_pak sysenv;

//const char *ws = "";  // send SOAP xml to stdout (deprec use logging instead)
const char ivec_ws[] = "https://grisu.ivec.org:8443/grisu-ws/services/grisu";
const char vpac_ws[] = "https://ngportal.vpac.org/grisu-ws/services/grisu";
const char vdev_ws[] = "https://ngportaldev.vpac.org/grisu-ws/services/grisu";
const char markus[] = "http://localhost:8080/grisu-ws/services/grisu";
const char *ws = vdev_ws;

struct soap soap;
struct SOAP_ENV__Header *grisu_soap_header=NULL;
gboolean grisu_auth = FALSE;
gchar *grisu_username=NULL;
gchar *grisu_password=NULL;
gchar *grisu_vo=NULL;
gchar *grisu_server=NULL;
gchar *grisu_port=NULL;
gchar *grisu_cert_path=NULL;

/******************************************************/
/* get time/date in suitable format for job timestamp */
/******************************************************/
/* TODO - re-locate? */
const gchar *timestamp(void)
{
time_t t1;
struct tm *t2;
static gchar text[30];  // FIXME - a little dodgy

t1 = time(NULL);
t2 = localtime(&t1);

if (strftime(text, 30, "%d-%m-%Y_%H:%M:%S", t2))
  return(text);

return("xx:yy");
}

/*************************************/
/* retrieve current credential setup */
/*************************************/
const gchar *grisu_username_get(void)
{
return(grisu_username);
}
const gchar *grisu_password_get(void)
{
return(grisu_password);
}
const gchar *grisu_vo_get(void)
{
return(grisu_vo);
}
const gchar *grisu_myproxy_get(void)
{
return(grisu_server);
}
const gchar *grisu_myproxyport_get(void)
{
return(grisu_port);
}
const gchar *grisu_ws_get(void)
{
return(ws);
}

/****************************************/
/* set the temporary credential keypair */
/****************************************/
void grisu_keypair_set(const gchar *username, const gchar *password)
{
if (username)
  {
  g_free(grisu_username);
  grisu_username = g_strdup(username);
  }
if (password)
  {
  g_free(grisu_password);
  grisu_password = g_strdup(password);
  }
}

/*******************************/
/* user's virtual organization */
/*******************************/
void grisu_vo_set(const gchar *vo)
{
if (vo)
  {
  g_free(grisu_vo);
  grisu_vo = g_strdup(vo);
  }
}

/********************************/
/* the service interface to use */
/********************************/
void grisu_ws_set(const gchar *grisu_ws)
{
ws = grisu_ws;
}

/********************************************************/
/* set the myproxy server:port to upload credentials to */
/********************************************************/
void grisu_myproxy_set(const gchar *server, const gchar *port)
{
if (server)
  {
  g_free(grisu_server);
  grisu_server = g_strdup(server);
  }
if (port)
  {
/* TODO - check it's a number? */
  g_free(grisu_port);
  grisu_port = g_strdup(port);
  }
}

/*************************************************************/
/* control if authentication tokens are placed in the header */
/*************************************************************/
void grisu_auth_set(gboolean value)
{
grisu_auth = value;
}

/*****************************************************/
/* build SOAP header with authentication information */
/*****************************************************/
/* NB: this should be called before EVERY soap call since, for reasons I don't understand, */
/* WS related soap calls seem to reset the header to NULL after they are sent off */
#define DEBUG_GRISU_AUTH_HEADER 0
void grisu_auth_header(void)
{
/* CURRENT - fixing crappy windows soap failure */
/* WINDOWS seems to require a reinit EVERY time - or it forgets stuff and soap calls fail */
/* TODO - test if this bollocks up linux or not ... */
soap_init1(&soap, SOAP_ENC_MTOM);

if (soap_ssl_client_context(&soap, SOAP_SSL_SKIP_HOST_CHECK, NULL, NULL, NULL, grisu_cert_path, NULL))
  printf("grisu_auth_header(): SSL init failed.\n");

/* allocate if necessary */
if (!grisu_soap_header)
  grisu_soap_header = g_malloc(sizeof(struct SOAP_ENV__Header));

/* fill out the header with current values */
grisu_soap_header->username = grisu_username;
grisu_soap_header->password = grisu_password;
grisu_soap_header->myproxyserver = grisu_server;
grisu_soap_header->myproxyport = grisu_port;

#if DEBUG_GRISU_AUTH_HEADER
printf("SOAP Header:      ptr = %p\n", grisu_soap_header);
printf("SOAP Header: username = %s\n", grisu_soap_header->username);
printf("SOAP Header: password = %s\n", grisu_soap_header->password);
printf("SOAP Header:   server = %s\n", grisu_soap_header->myproxyserver);
printf("SOAP Header:     port = %s\n", grisu_soap_header->myproxyport);
#endif

/* use the header if auth is turned on */
if (grisu_auth)
  soap.header = grisu_soap_header;
else
  soap.header = NULL;
}

/***********************************/
/* turn soap logging plugin on/off */
/***********************************/
#define GSOAP_LOGGING_SENT "SENT.log"
#define GSOAP_LOGGING_RECV "RECV.log"
void grisu_soap_logging(gboolean state)
{
struct logging_data *logdata;

/* retrieve data struct */
logdata = (struct logging_data *) soap_lookup_plugin(&soap, logging_id);
if (!logdata)
  {
  if (soap_register_plugin(&soap, logging))
     soap_print_fault(&soap, stderr);
  else
    {
/* register plugin if not loaded */
    printf("Loaded gsoap logging plugin...\n");
    logdata = (struct logging_data *) soap_lookup_plugin(&soap, logging_id);
    }
  }

/* set logging state */
if (logdata)
  {
  if (state)
    {
/* TODO - open and write something to indicate logging on/off switches? */
/* then append rather than overwrite (for clarity) */
//    logdata->inbound = stdout; 
//    logdata->outbound = stdout;
//    logdata->inbound = fopen("/home/sean/prog/gdis/RECV.log", "at");
//    logdata->outbound = fopen("/home/sean/prog/gdis/SENT.log", "at");
    logdata->outbound = fopen(GSOAP_LOGGING_SENT, "wt");
    logdata->inbound = fopen(GSOAP_LOGGING_RECV, "wt");
    }
  else
    {
    logdata->inbound = NULL; 
    logdata->outbound = NULL;
    }
// process messages 
//  size_t bytes_in = logdata->stat_recv;
//  size_t bytes_out = logdata->stat_sent;
  }
else
  printf("Gsoap logging plugin not found.\n");
}

/******************/
/* xml primitives */
/******************/
void grisu_xml_start(GMarkupParseContext *context,
                     const gchar *element_name,
                     const gchar **attribute_names,
                     const gchar **attribute_values,
                     gpointer data,
                     GError **error)
{
gint i;
GSList **list = data;

/* process the file element */
if (g_ascii_strncasecmp(element_name, "file\0", 5) == 0)
  {
  i=0;
  while (*(attribute_names+i))
    {
    if (g_ascii_strncasecmp(*(attribute_names+i), "name", 4) == 0)
      {
      *list = g_slist_prepend(*list, g_strdup(*(attribute_values+i)));
      }
    if (g_ascii_strncasecmp(*(attribute_names+i), "size", 4) == 0)
      {
/* TODO - do somethin with the file size (could be useful) */
      }
    i++;
    }
  }
}

void grisu_xml_end(GMarkupParseContext *context,
                   const gchar *element_name,
                   gpointer data,
                   GError **error)
{
}

void grisu_xml_text(GMarkupParseContext *context,
                    const gchar *text,
                    gsize text_len,  
                    gpointer data,
                    GError **error)
{
}

/******************/
/* xml primitives */
/******************/
/* CURRENT - grisu xml parsing (1) file lists & 2) job description) */
gpointer grisu_xml_parse(const gchar *xml_text)
{
GMarkupParser xml_parser;
GMarkupParseContext *xml_context;
GError *error;
GSList *file_list=NULL;

xml_parser.start_element = &grisu_xml_start;
xml_parser.end_element = &grisu_xml_end;
xml_parser.text = &grisu_xml_text;
xml_parser.passthrough = NULL;
xml_parser.error = NULL;

xml_context = g_markup_parse_context_new(&xml_parser, 0, &file_list, NULL);

if (!g_markup_parse_context_parse(xml_context, xml_text, strlen(xml_text), &error))
    printf("grisu_xml_parse() : parsing error.\n");

/* cleanup */
if (!g_markup_parse_context_end_parse(xml_context, &error))
  printf("grisu_xml_parse() : errors occurred processing xml.\n");

g_markup_parse_context_free(xml_context);

return(file_list);
}

/************************************************/
/* remote directory listing (url or virtual fs) */
/************************************************/
GSList *grisu_file_list(const gchar *dir)
{
GSList *list=NULL;
struct _ns1__ls_USCOREstring put;
struct _ns1__ls_USCOREstringResponse get;

g_assert(dir != NULL);

put.in0 = g_strdup(dir);
put.in1 = 1;
put.in2 = 1;

/* prereqs */
grisu_auth_set(TRUE);
grisu_auth_header();

if (soap_call___ns1__ls_USCOREstring(&soap, ws, "", &put, &get) == SOAP_OK)
  {
  list = grisu_xml_parse(get.out);
  }
else
  {
  printf("grisu_file_list(): SOAP ERROR\n");
  soap_print_fault(&soap, stderr);
  }

g_free(put.in0);

return(list);
}

/*******************************/
/* obtain current mount points */
/*******************************/
// FIXME - often seem to get compile complaints about type missmatches on this fn
#define GRISU_DF 0
#if GRISU_DF
void grisu_df(void)
{
gint i;
struct _ns1__df put;
struct _ns1__dfResponse get;

printf("*** GRISU WS - scanning mountpoints ***\n");

/* prereqs */
grisu_auth_set(TRUE);
grisu_auth_header();

if (soap_call___ns1__df(&soap, ws, "", &put, &get) == SOAP_OK)
  {
  struct ns3__ArrayOfMountPoint *amp = get.out;	

printf("returned entries: %d\n", amp->__sizeMountPoint);

for (i=0 ; i<amp->__sizeMountPoint ; i++)
  {
/* FIXME - the type in the stub seems to change when generated */
/* on different versions of ubunutu */
  struct ns3__MountPoint *mp = *(amp->MountPoint+i);

//  printf("[%s] [%s] [%s] [%s]\n", mp->dn, mp->fqan, mp->mountpoint, mp->rootUrl);
  printf("[%s] [%s]\n", mp->fqan, mp->rootUrl);
  }

  }
else
  {
  printf(" *** soap error ***\n");
  soap_print_fault(&soap, stderr);
  }

}
#endif

/*****************************************/
/* mount a remote filesystem for staging */
/*****************************************/
/* currently unused */
gboolean grisu_mount(const gchar *mount_point, const gchar *url)
{
gint value;
struct _ns1__mount put;
struct _ns1__mountResponse get;

printf("*** GRISU WS - creating mountpoint ***\n");

g_assert(mount_point != NULL);
g_assert(url != NULL);

put.in0 = g_strdup(url);
put.in1 = g_strdup(mount_point);
put.in2 = 1;

/* prereqs */
grisu_auth_set(TRUE);
grisu_auth_header();

if (soap_call___ns1__mount(&soap, ws, "", &put, &get) == SOAP_OK)
  {
  value = TRUE;
  }
else
  {
  printf("grisu_mount(): SOAP ERROR\n");
  soap_print_fault(&soap, stderr);
  value = FALSE;
  }

g_free(put.in0);
g_free(put.in1);

return(value);
}

/*************************************************************/
/* mount a remote filesystem (under supplied VO) for staging */
/*************************************************************/
/* currently unused */
gboolean grisu_mount_vo(const gchar *mount_point, const gchar *url, const gchar *vo)
{
gint value;
struct _ns1__mount1 put;
struct _ns1__mount1Response get;

printf("*** GRISU WS - creating mountpoint with VO ***\n");

g_assert(mount_point != NULL);
g_assert(url != NULL);
g_assert(vo != NULL);

put.in0 = g_strdup(url);
put.in1 = g_strdup(mount_point);
put.in2 = g_strdup(vo);
put.in3 = 1;

/* prereqs */
grisu_auth_set(TRUE);
grisu_auth_header();

if (soap_call___ns1__mount1(&soap, ws, "", &put, &get) == SOAP_OK)
  {
  value = TRUE;
  }
else
  {
  printf("grisu_mount_vo(): SOAP ERROR\n");
  soap_print_fault(&soap, stderr);
  value = FALSE;
  }

g_free(put.in0);
g_free(put.in1);
g_free(put.in2);

return(value);
}

/******************************/
/* streaming upload callbacks */
/******************************/
void *mime_read_open(struct soap *soap, void *handle, const char *id, const char *type, const char *description)
{
//printf("mime: opening outbound (%s)\n", type);
return handle;
}
void mime_read_close(struct soap *soap, void *handle)
{
//printf("mime: closing outbound\n");
fclose((FILE*)handle);
}
size_t mime_read(struct soap *soap, void *handle, char *buf, size_t len)
{ 
return fread(buf, 1, len, (FILE*)handle);
}

/***********************************************/
/* upload a file to the jobs working directory */
/***********************************************/
gint grisu_file_upload(const gchar *fullpath, const gchar *remotedir)
{
gint value;
struct _ns1__upload put;
struct _ns1__uploadResponse get;
struct xsd__base64Binary image;
struct stat sb;
FILE *fd;

g_assert(fullpath != NULL);
g_assert(remotedir != NULL);

fd = fopen(fullpath, "rb");
if (!fd)
  {
  printf("grisu_file_upload() ERROR, could not find: %s\n", fullpath);
  return(FALSE);
  }

/* pre-req */
grisu_auth_set(TRUE);
grisu_auth_header();

/* try to get file length */
if (!fstat(fileno(fd), &sb) && sb.st_size > 0)
   {
   soap.fmimereadopen = mime_read_open;
   soap.fmimereadclose = mime_read_close;
   soap.fmimeread = mime_read;
   image.__ptr = (unsigned char *) fd;
/* NB: size must be set */
   image.__size = sb.st_size;
   }
 else
   {
   printf("grisu_file_upload() ERROR, couldn't get size of file: %s\n", fullpath);
   return(FALSE);

/* TODO - don't know/can't get the size - buffer it? - deal with this later */
/* see gsoap docs on Streaming Chunked MTOM/MIME maybe? */
/*
printf("unknown size ...\n");
   size_t i;
   int c;
   image.__ptr = (unsigned char*)soap_malloc(&soap, MAX_FILE_SIZE);
   for (i = 0; i < MAX_FILE_SIZE; i++)
      {
         if ((c = fgetc(fd)) == EOF)
            break;
         image.__ptr[i] = c;
      }
   fclose(fd);
   image.__size = i;
*/
   }

/* NB - must set id & options to NULL */
/* or else crap characters can be inserted that screw up transmission */
image.id = NULL;
image.options = NULL;
/* MIME type */
/* FIXME - do we need to use image/xxx for binary uploads? */
image.type = "text/html";

/* source is the MIME/MTOM attachments */
put.in0 = &image;
put.in1 = g_strdup(remotedir);
put.in2 = TRUE;

if (soap_call___ns1__upload(&soap, ws, "", &put, &get) == SOAP_OK)
  {
  value = TRUE;
  }
else
  {
  printf("grisu_file_upload(): SOAP ERROR\n");
  soap_print_fault(&soap, stderr);
  value = FALSE;
  }

g_free(put.in1);

return(value);
}

/********************************/
/* streaming download callbacks */
/********************************/
void *mime_write_open(struct soap *soap,
                      void *unused_handle,
                      const char *id, const char *type, const char *desc,
                      enum soap_mime_encoding encoding)
{
g_assert(soap != NULL);

/* TODO - cope with untested MIME encoding types */
/*
printf("mime: opening inbound (%s)\n", type);
printf("mime: encoding type ");
*/

switch (encoding)
  {
/* known to work */
  case SOAP_MIME_BINARY:
//    printf("ok\n");
    break;
/* untested */
  case SOAP_MIME_NONE:
    printf("untested 1\n");
    break;
  case SOAP_MIME_7BIT:
    printf("untested 2\n");
    break;
  case SOAP_MIME_8BIT:
    printf("untested 3\n");
    break;
  case SOAP_MIME_QUOTED_PRINTABLE:
    printf("untested 4\n");
    break;
  case SOAP_MIME_BASE64:
    printf("untested 5\n");
    break;
  case SOAP_MIME_IETF_TOKEN:
    printf("untested 6\n");
    break;
  case SOAP_MIME_X_TOKEN:
    printf("untested 7\n");
    break;
/* unknown */
  default:
    printf("unknown\n");
  }

/* CURRENT - using a pointer in the soap struct for setting the destination file */

return(soap->os);
}

/********************************/
/* streaming download callbacks */
/********************************/
void mime_write_close(struct soap *soap, void *handle)
{
//printf("mime: closing inbound\n");
fclose((FILE*)handle);
}

/********************************/
/* streaming download callbacks */
/********************************/
int mime_write(struct soap *soap, void *handle, const char *buf, size_t len)
{
size_t nwritten;

while (len)
  {
  nwritten = fwrite(buf, 1, len, (FILE*)handle);
  if (!nwritten)
    {
    soap->errnum = errno;
    return SOAP_EOF;
    }
  len -= nwritten;
  buf += nwritten;
  }
return SOAP_OK;
}

/***************************************************/
/* download a file from the jobs working directory */
/***************************************************/
gint grisu_file_download(const gchar *remote_src, const gchar *local_dest)
{
gint value;
struct _ns1__download put;
struct _ns1__downloadResponse get;
FILE *fd;

g_assert(remote_src != NULL);
g_assert(local_dest != NULL);

/* prereqs */
grisu_auth_set(TRUE);
grisu_auth_header();

/* CURRENT - using a pointer in the current soap struct to pass destination file for writing */
fd = fopen(local_dest, "wb");
if (fd)
  {
/* FIXME - a bit naughty as I don't know what soap.os (output stream?) is for */
  if (soap.os)
    printf("grisu_file_download() warning, soap.os pointer was not NULL\n");

  soap.os = (void *) fd;
  }
else
  {
  printf("grisu_file_download() ERROR, couldn't open [%s] for writing\n", local_dest);
  return(FALSE);
  }

soap.fmimewriteopen = mime_write_open;
soap.fmimewriteclose = mime_write_close;
soap.fmimewrite = mime_write; 
/*
soap.connect_timeout = 10;
soap.send_timeout = 30;
soap.recv_timeout = 30;
*/

/* remote source URL */
put.in0 = g_strdup(remote_src);

if (soap_call___ns1__download(&soap, ws, "", &put, &get) == SOAP_OK)
  {
  value = TRUE;
  }
else
  {
  printf("grisu_file_download(): SOAP ERROR\n");
  soap_print_fault(&soap, stderr);
  value = FALSE;
  }

g_free(put.in0);

return(value);
}

/*********************/
/* get all job names */
/*********************/
GSList *grisu_job_names(void)
{
gint i;
struct _ns1__getAllJobnames put;
struct _ns1__getAllJobnamesResponse get;
struct ns1__ArrayOfString *out; 
GSList *list=NULL;

/* prereqs */
grisu_auth_set(TRUE);
grisu_auth_header();

if (soap_call___ns1__getAllJobnames(&soap, ws, "", &put, &get) == SOAP_OK)
  {
  out = get.out;
  for (i=out->__sizestring ; i-- ; )
    list = g_slist_prepend(list, g_strdup(out->string[i]));
  }
else
  {
  printf("grisu_job_names(): SOAP ERROR\n");
  soap_print_fault(&soap, stderr);
  }

return(list);
}

/**************************/
/* kill and cleanup a job */
/**************************/
void grisu_job_remove(const gchar *name)
{
struct _ns1__kill put;
struct _ns1__killResponse get;

/* prereq */
g_assert(name != NULL);

grisu_auth_set(TRUE);
grisu_auth_header();

put.in0 = g_strdup(name);
/* TODO - true = delete files? if so implement a stop */
/* and a stop + cleanup */
put.in1 = 1;

if (soap_call___ns1__kill(&soap, ws, "", &put, &get) == SOAP_OK)
  {
//printf("Job deleted.\n");
  }
else
  {
  printf("grisu_job_remove(): SOAP ERROR\n");
  soap_print_fault(&soap, stderr);
  }

g_free(put.in0);

}

/************************************/
/* submit a job to the grisu server */
/************************************/
gint grisu_job_submit(const gchar *name, const gchar *vo)
{
gint value;
struct _ns1__submitJob put;
struct _ns1__submitJobResponse get;

g_assert(name != NULL);

/* prereqs */
grisu_auth_set(TRUE);
grisu_auth_header();

/* input */
put.in0 = g_strdup(name);
if (vo)
  put.in1 = g_strdup(vo);
else
  put.in1 = g_strdup("/ARCS/Startup");

/* submit the job */
if (soap_call___ns1__submitJob(&soap, ws, "", &put, &get) == SOAP_OK)
  {
  value = TRUE;
  }
else
  {
  printf("grisu_job_submit(): SOAP ERROR\n");
  soap_print_fault(&soap, stderr);
  value = FALSE;
  }

g_free(put.in0);
g_free(put.in1);

return(value);
}

/****************************************/
/* get job type from the jobname string */
/****************************************/
/* NB: may be deprecated if a better way of recording the job type is developed */
gint grisu_job_type(const gchar *name)
{
g_assert(name != NULL);

if (g_strncasecmp(name, "gulp", 4) == 0)
  return(JOB_GULP);

if (g_strncasecmp(name, "gamess", 6) == 0)
  return(JOB_GAMESS);

printf("grisu_job_type() - unknown type [%s]\n", name);

return(JOB_UNKNOWN);
}

/*************************************************/
/* generate executable specific job descriptions */
/*************************************************/
gchar *grisu_xml_exe(gpointer data)
{
gint type;
gchar *input;
GString *xml;
struct grid_pak *grid = data;

g_assert(grid != NULL);

xml = g_string_new(NULL);

type = grisu_job_type(grid->jobname);

switch (type)
  {
  case JOB_GULP:
    g_string_append_printf(xml, "<Executable>%s</Executable>\n", grid->remote_exe);
    g_string_append_printf(xml, "<Input>%s</Input>\n", grid->local_input);
    g_string_append_printf(xml, "<Output>%s</Output>\n", grid->local_output);
    g_string_append_printf(xml, "<Error>stderr.txt</Error>\n");
    break;

  case JOB_GAMESS:
    g_string_append_printf(xml, "<Executable>%s</Executable>\n", grid->remote_exe);
/* GAMESS likes to be different - strip off input file's extension */
    input = parse_strip_extension(grid->local_input);
    g_string_append_printf(xml, "<Argument>%s</Argument>\n", input);
    g_free(input);
    g_string_append_printf(xml, "<Argument>%s</Argument>\n", grid->exe_version);
    g_string_append_printf(xml, "<Argument>%d</Argument>\n", grid->remote_exe_np); 
    g_string_append_printf(xml, "<Output>%s</Output>\n", grid->local_output);
    g_string_append_printf(xml, "<Error>stderr.txt</Error>\n");
    break;

  default:
    printf("Error - unknown job type.\n");
  }

if (grid->remote_exe_module)
 g_string_append_printf(xml, "<Module>%s</Module>\n", grid->remote_exe_module);

return(g_string_free(xml, FALSE));
}

/***********************************************/
/* generate a job description string for grisu */
/***********************************************/
/* TODO - more sophisticated XML generator, eg for CPU time etc values */
gchar *grisu_xml_create(gpointer data)
{
GString *xml;
struct grid_pak *grid = data;

g_assert(grid != NULL);
g_assert(grid->jobname != NULL);
g_assert(grid->remote_q != NULL);
g_assert(grid->remote_exe != NULL);
g_assert(grid->remote_root != NULL);

if (grid->remote_exe_type == GRID_MPI)
  {
//  printf("TODO - mpi flag in job XML\n");
  }

xml = g_string_new(NULL);

g_string_append_printf(xml, "<JobDefinition xmlns=\"http://schemas.ggf.org/jsdl/2005/11/jsdl\">\n");
g_string_append_printf(xml, "<JobDescription>\n");
g_string_append_printf(xml, "<JobIdentification><JobName>%s</JobName></JobIdentification>\n", grid->jobname);
g_string_append_printf(xml, "<Application>\n");
g_string_append_printf(xml, "<ApplicationName>uname</ApplicationName>\n");
g_string_append_printf(xml, "<POSIXApplication xmlns=\"http://schemas.ggf.org/jsdl/2005/11/jsdl-posix\">\n");
g_string_append_printf(xml, "<WorkingDirectory filesystemName=\"userExecutionHostFs\">");
g_string_append_printf(xml, "grisu-jobs-dir/%s</WorkingDirectory>\n", grid->jobname);

/* CURRENT - exectuable specific xml */
g_string_append_printf(xml, grisu_xml_exe(grid));

g_string_append_printf(xml, "</POSIXApplication>\n");

/* TODO - change this to suit job type later */
g_string_append_printf(xml, "<TotalCPUTime>60</TotalCPUTime>\n");
g_string_append_printf(xml, "<TotalCPUCount>%d</TotalCPUCount>\n", grid->remote_exe_np);
g_string_append_printf(xml, "</Application>\n");

g_string_append_printf(xml, "<Resources>\n");
g_string_append_printf(xml, "<CandidateHosts><HostName>%s</HostName></CandidateHosts>\n", grid->remote_q);
g_string_append_printf(xml, "<FileSystem name=\"userExecutionHostFs\">\n");
g_string_append_printf(xml, "<MountSource>%s</MountSource>\n", grid->remote_root);
g_string_append_printf(xml, "<FileSystemType>normal</FileSystemType>\n");
g_string_append_printf(xml, "</FileSystem></Resources></JobDescription></JobDefinition>\n");

return(g_string_free(xml, FALSE));
}

/*************************************/
/* get a valid grisu job name to use */
/*************************************/
gchar *grisu_jobname_request(gint type)
{
gint i;
GString *name;
struct _ns1__createJob put;
struct _ns1__createJobResponse get;

/* prereqs */
grisu_auth_set(TRUE);
grisu_auth_header();

name = g_string_new(NULL);

/* repeat for max tries */
for (i=3 ; i-- ; )
  {
  switch (type)
    {
    case JOB_GULP:
      g_string_printf(name, "gulp_%s", timestamp());
      break;

    case JOB_GAMESS:
      g_string_printf(name, "gamess_%s", timestamp());
      break;

    default:
printf("Error, unknown job type.\n");
      return(NULL);
    }

  put.in0 = name->str;
  put.in1 = 0;

  if (soap_call___ns1__createJob(&soap, ws, "", &put, &get) == SOAP_OK)
    return(g_string_free(name, FALSE));
  else
    {
    printf("grisu_jobname_request() FAILED, trying a new name\n");
    soap_print_fault(&soap, stderr);
    }
  }

g_string_free(name, TRUE);

return(NULL);
}

/*******************************************/
/* attach a job description to a grisu job */
/*******************************************/
gint grisu_job_configure(gchar *job_name, gchar *job_xml)
{
struct _ns1__setJobDescription_USCOREstring put;
struct _ns1__setJobDescription_USCOREstringResponse get;

g_assert(job_name != NULL);
g_assert(job_xml != NULL);

/* prereqs */
grisu_auth_set(TRUE);
grisu_auth_header();

/* args */
put.in0 = job_name;
put.in1 = job_xml;

if (soap_call___ns1__setJobDescription_USCOREstring(&soap, ws, "", &put, &get) == SOAP_OK)
  return(TRUE);
else
  {
  printf("grisu_job_configure(): SOAP ERROR\n");
  soap_print_fault(&soap, stderr);
  }

return(FALSE);
}

/***********************************/
/* get the job's working directory */
/***********************************/
gchar *grisu_job_cwd(const gchar *name)
{
gchar *cwd=NULL;
struct _ns1__getJobDirectory put;
struct _ns1__getJobDirectoryResponse get;

/* prereqs */
g_assert(name != NULL);

grisu_auth_set(TRUE);
grisu_auth_header();

put.in0 = g_strdup(name);

if (soap_call___ns1__getJobDirectory(&soap, ws, "", &put, &get) == SOAP_OK)
  {
/* FIXME - understand return values - should this be strdup'ed? */
  cwd = get.out;
  }
else
  {
  printf("grisu_job_cwd(): SOAP ERROR\n");
  soap_print_fault(&soap, stderr);
  }

g_free(put.in0);

return(cwd);
}

/***************************/
/* get the status of a job */
/***************************/
const gchar *grisu_job_status(const gchar *name)
{
struct _ns1__getJobStatus put;
struct _ns1__getJobStatusResponse get;

/* prereqs */
g_assert(name != NULL);

grisu_auth_set(TRUE);
grisu_auth_header();

put.in0 = g_strdup(name);

if (soap_call___ns1__getJobStatus(&soap, ws, "", &put, &get) == SOAP_OK)
  {

g_free(put.in0);

  if (get.out > 1000)
    {
/* TODO - exit code = (out - 1000) */
    return("Done");
    }

  switch (get.out)
    {
    case -1002:
      return("Loading...");
    case -1001:
      return("N/A");
    case -1000:
      return("Undefined");
    case -999:
      return("No such job");
    case -100:
      return("Job created");
    case -4:
      return("Ready to submit");
    case -3:
      return("External handle ready");
    case -2:
      return("Staging in");
    case -1:
      return("Unsubmitted");
    case 0:
      return("Pending");
    case 1:
      return("Active");
    case 101:
      return("Cleaning up"); 
    case 900:
      return("Unknown");
    case 998:
      return("Killed");
    case 999:
      return("Failed");
    case 1000:
      return("Done");
    }
  }
else
  {
  printf("grisu_job_status(): SOAP ERROR\n");
  soap_print_fault(&soap, stderr);

g_free(put.in0);
  }

return(NULL);
}

/***************************************/
/* get job details - grisu XML jobfile */
/***************************************/
gchar *grisu_job_details(const gchar *name)
{
gchar *xml=NULL;
struct _ns1__getJobDetails_USCOREstring put;
struct _ns1__getJobDetails_USCOREstringResponse get;

/* prereqs */
g_assert(name != NULL);

grisu_auth_set(TRUE);
grisu_auth_header();

put.in0 = g_strdup(name);

if (soap_call___ns1__getJobDetails_USCOREstring(&soap, ws, "", &put, &get) == SOAP_OK)
  {
  xml = get.out;
  }
else
  {
  printf("grisu_job_details(): SOAP ERROR\n");
  soap_print_fault(&soap, stderr);
  }

g_free(put.in0);

return(xml);
}

/**********************************/
/* get current list of grisu jobs */
/**********************************/
/* Currently unused */
gchar *grisu_ps(void)
{
struct _ns1__ps_USCOREstring put;
struct _ns1__ps_USCOREstringResponse get;

printf(" *** GRISU WS - get job information ***\n");

/* prereqs */
grisu_auth_set(TRUE);
grisu_auth_header();

if (soap_call___ns1__ps_USCOREstring(&soap, ws, "", &put, &get) == SOAP_OK)
  {
printf("Response: %s\n", get.out);
  }
else
  {
  printf(" *** soap error ***\n");
  soap_print_fault(&soap, stderr);
  }

return(NULL);
}

/*****************************************/
/* retreive a list of available user VOs */
/*****************************************/
GSList *grisu_fqans_get(void)
{
struct _ns1__getFqans put;
struct _ns1__getFqansResponse get;
GSList *list=NULL;

/* prereqs */
grisu_auth_set(TRUE);
grisu_auth_header();

if (soap_call___ns1__getFqans(&soap, ws, "", &put, &get) == SOAP_OK)
  {
  int i;
  struct ns1__ArrayOfString *out = get.out;

  for (i=0 ; i<out->__sizestring ; i++)
    {
    list = g_slist_prepend(list, g_strdup(out->string[i]));
//    printf("[%d] %s\n", i, out->string[i]);
    }
  }
else
  {
  printf("grisu_fqans_get(): SOAP ERROR\n");
  soap_print_fault(&soap, stderr);
  }

return(list);
}

/*******************************************************/
/* request the remote job's absolute working directory */
/*******************************************************/
gchar *grisu_absolute_job_dir(const gchar *jobname, const gchar *submit, const gchar *vo)
{
gchar *dir=NULL;
struct _ns1__calculateAbsoluteJobDirectory put;
struct _ns1__calculateAbsoluteJobDirectoryResponse get;

g_assert(jobname != NULL);
g_assert(submit != NULL);

/* prereqs */
grisu_auth_set(TRUE);
grisu_auth_header();

/* input */
put.in0 = g_strdup(jobname);
put.in1 = g_strdup(submit);

if (vo)
  put.in2 = g_strdup(vo);
else
  put.in2 = g_strdup("");

if (soap_call___ns1__calculateAbsoluteJobDirectory(&soap, ws, "", &put, &get) == SOAP_OK)
  {
  dir = get.out;
  }
else
  {
  printf("grisu_absolute_job_dir(): SOAP ERROR\n");
  soap_print_fault(&soap, stderr);
  }

g_free(put.in0);
g_free(put.in1);
g_free(put.in2);

return(dir);
}

/**********************************************/
/* request the remote job's working directory */
/**********************************************/
gchar *grisu_relative_job_dir(const gchar *jobname)
{
gchar *dir=NULL;
struct _ns1__calculateRelativeJobDirectory put;
struct _ns1__calculateRelativeJobDirectoryResponse get;

g_assert(jobname != NULL);

/* prereqs */
grisu_auth_set(TRUE);
grisu_auth_header();

put.in0 = g_strdup(jobname);

if (soap_call___ns1__calculateRelativeJobDirectory(&soap, ws, "", &put, &get) == SOAP_OK)
  {
  dir = get.out;
  }
else
  {
  printf("grisu_relative_job_dir(): SOAP ERROR\n");
  soap_print_fault(&soap, stderr);
  }

g_free(put.in0);

return(dir);
}

/*************************/
/* query available sites */
/*************************/
GSList *grisu_site_list(void)
{
gint i;
GSList *list=NULL;
struct _ns1__getAllSites put;
struct _ns1__getAllSitesResponse get;
struct ns1__ArrayOfString *out;

/* prereqs */
grisu_auth_set(TRUE);
grisu_auth_header();

/* query WS for available sites */
if (soap_call___ns1__getAllSites(&soap, ws, "", &put, &get) == SOAP_OK)
  {
  out = get.out;
  for (i=out->__sizestring ; i-- ; )
    list = g_slist_prepend(list, g_strdup(out->string[i]));
  }
else
  {
  printf(" *** soap error ***\n");
  soap_print_fault(&soap, stderr); 
  }

return(list);
}

/******************************************/
/* query available applications at a site */
/******************************************/
GSList *grisu_application_list(const gchar *site)
{
gint i;
GSList *list=NULL;
struct _ns1__getAllAvailableApplications put;
struct _ns1__getAllAvailableApplicationsResponse get;
struct ns1__ArrayOfString tmp, *out;

g_assert(site != NULL);

/* prereqs */
grisu_auth_set(TRUE);
grisu_auth_header();

/* args */
tmp.__sizestring = 1;
tmp.string = malloc(sizeof(char *));
tmp.string[0] = g_strdup(site); 
put.in0 = &tmp;

if (soap_call___ns1__getAllAvailableApplications(&soap, ws, "", &put, &get) == SOAP_OK)
  {
  out = get.out;
  for (i=out->__sizestring ; i-- ; )
    list = g_slist_prepend(list, g_strdup(out->string[i]));
  }
else
  {
  printf("grisu_application_list(): SOAP ERROR\n");
  soap_print_fault(&soap, stderr); 
  }

g_free(tmp.string[0]);
g_free(tmp.string);

return(list);
}

/*****************************************************/
/* retrieve the site name from a submission location */
/*****************************************************/
/* TODO - cache - getAllHosts() value */
/* TODO - during grid connect phase - grab & cache a must as we can */
gchar *grisu_site_name(const gchar *submit)
{
gchar **host, *site=NULL;
struct _ns1__getSite put;
struct _ns1__getSiteResponse get;

g_assert(submit != NULL);

/* CURRENT - bit of a hack to get the site name of a submission location */
/* can only query the site of a host ... so we process the submit location */
/* to get this for the query */

host = g_strsplit(submit, ":", 3);

/* TODO - some sites have #Loadlevel crap after the hostname - get rid of this? */
/* NEW - grisu expected to be updated to cope with submission location as input */
//printf(" *** GRISU WS - requesting site name for [%s] ***\n", host[1]);

/* prereqs */
grisu_auth_set(TRUE);
grisu_auth_header();

put.in0 = host[1];

if (soap_call___ns1__getSite(&soap, ws, "", &put, &get) == SOAP_OK)
  {
  if (get.out)
    site = g_strdup(get.out);
  }
else
  {
  printf("grisu_site_name(): SOAP ERROR\n");
  soap_print_fault(&soap, stderr); 
  }

g_strfreev(host);

return(site);
}

/*******************************************************************/
/* find all queue/host pairs where an application can be submitted */
/*******************************************************************/
GSList *grisu_submit_find(const gchar *name)
{
gint i;
GSList *list=NULL;
struct _ns1__getSubmissionLocationsForApplication put;
struct _ns1__getSubmissionLocationsForApplicationResponse get;
struct ns1__ArrayOfString *out;

g_assert(name != NULL);

/* prereqs */
grisu_auth_header();

put.in0 = g_strdup(name);

if (soap_call___ns1__getSubmissionLocationsForApplication(&soap, ws, "", &put, &get) == SOAP_OK)
  {
  out = get.out;
  for (i=out->__sizestring ; i-- ; )
    list = g_slist_prepend(list, g_strdup(out->string[i]));
  }
else
  {
  printf("grisu_submit_find(): soap error\n");
  soap_print_fault(&soap, stderr); 
  }

g_free(put.in0);

return(list);
}

/****************************************/
/* get version names of site executable */
/****************************************/
GSList *grisu_application_versions(const gchar *name, const gchar *site)
{
gint i;
GSList *list=NULL;
struct _ns1__getVersionsOfApplicationOnSite put;
struct _ns1__getVersionsOfApplicationOnSiteResponse get;
struct ns1__ArrayOfString *out;

g_assert(name != NULL);
g_assert(site != NULL);

put.in0 = g_strdup(name);
put.in1 = g_strdup(site);

/* prereqs */
grisu_auth_set(TRUE);
grisu_auth_header();

if (soap_call___ns1__getVersionsOfApplicationOnSite(&soap, ws, "", &put, &get) == SOAP_OK)
  {
  out = get.out;
  for (i=out->__sizestring ; i-- ; )
    list = g_slist_prepend(list, g_strdup(out->string[i]));
  }
else
  {
  printf(" *** soap error ***\n");
  soap_print_fault(&soap, stderr); 
  }

g_free(put.in0);
g_free(put.in1);

return(list);
}

/*******************************************************/
/* get information about how the exectuable is invoked */
/*******************************************************/
#define DEBUG_GRISU_APPLICATION_DETAILS 0
GHashTable *grisu_application_details(const gchar *name, const gchar *site)
{
gint i;
GHashTable *details=NULL;
struct _ns1__getApplicationDetails put;
struct _ns1__getApplicationDetailsResponse get;
struct ns1__anyType2anyTypeMap *map;	
struct _ns1__anyType2anyTypeMap_entry *entry;

g_assert(name != NULL);
g_assert(site != NULL);

details = g_hash_table_new_full(g_str_hash, g_str_equal, g_free, g_free);

put.in0 = g_strdup(name);
put.in1 = g_strdup(site);

/* prereqs */
grisu_auth_set(TRUE);
grisu_auth_header();

if (soap_call___ns1__getApplicationDetails(&soap, ws, "", &put, &get) == SOAP_OK)
  {
  map = get.out;

  for (i=0 ; i<map->__sizeentry ; i++)
    {
    entry = map->entry+i;

    g_hash_table_insert(details, g_strdup(entry->key), g_strdup(entry->value));

#if DEBUG_GRISU_APPLICATION_DETAILS
printf("[%s] = [%s]\n", entry->key, entry->value); 
#endif
    }
  }
else
  {
  printf(" *** soap error ***\n");
  soap_print_fault(&soap, stderr); 
  }

g_free(put.in0);
g_free(put.in1);

return(details);
}

/******************************************/
/* initialize subsequent webservice calls */
/******************************************/
int grisu_init(void)
{
gint value=0;
const gchar *home;

/* soap setup */
//soap_init(&soap);
/* enable MTOM for file xfer -> attachment/streaming */
soap_init1(&soap, SOAP_ENC_MTOM);

/* grisu init */
if (!grisu_username)
  grisu_username = grid_random_alpha(10);
if (!grisu_password)
  grisu_password = grid_random_alphanum(12);
if (!grisu_server)
  grisu_server=g_strdup("myproxy.arcs.org.au");

// FIXME - sometimes needed when myproxy plays up also might need to setenv
// MYPROXY_SERVER_DN = /C=AU/O=APACGrid/OU=VPAC/CN=myproxy2.arcs.org.au
//  grisu_server=g_strdup("myproxy2.arcs.org.au");

if (!grisu_port)
  grisu_port=g_strdup("7512");


/* certificate path init */
home = g_get_home_dir();
if (home)
  grisu_cert_path = g_build_filename(home, ".globus", "certificates", NULL);
if (!g_file_test(grisu_cert_path, G_FILE_TEST_IS_DIR))
  {
  printf("ERROR in grisu_init(), bad certificate location: %s\n", grisu_cert_path);
  g_free(grisu_cert_path);

/* don't return if no certs - we can soldier on */
/* as it's only myproxy-init that won't work */
  grisu_cert_path=NULL;
  }

#if DEBUG_GRISU_CLIENT
printf("Using certificate location: %s\n", grisu_cert_path);
#endif

/* use SOAP_SSL_DEFAULT in production code */
/* keyfile: required only when client must authenticate to server */
/* password to read the keyfile */
/* optional cacert file to store trusted certificates */
/* optional capath to directory with trusted certificates */
/* if randfile!=NULL: use a file with random data to seed randomness */ 
/*
if (soap_ssl_client_context(&soap, SOAP_SSL_NO_AUTHENTICATION, NULL, NULL, NULL, NULL, NULL))
*/

/* CURRENT - to verify server we use the standard (CA signed) certificates */

/* CURRENT - if no certs is this secure? */
/*
if (soap_ssl_client_context(&soap, SOAP_SSL_DEFAULT, NULL, NULL, NULL, path, NULL))
*/


/* CURRENT - this has to be done before all calls anyway ... */
/* create header for soap struct with grisu authentication info */
/*
grisu_auth_header();
*/

/*
if (soap_ssl_client_context(&soap, SOAP_SSL_SKIP_HOST_CHECK, NULL, NULL, NULL, path, NULL))
  { 
  soap_print_fault(&soap, stderr);
  value = 1;
  }
else
  {
#if DEBUG_GRISU_CLIENT
  printf("Requiring secure context.\n");
#endif
  }
*/

return(value);
}

/**********************/
/* webservice cleanup */
/**********************/
void grisu_stop(void)
{
/* grisu */
g_free(grisu_username);
g_free(grisu_password);
g_free(grisu_server);
g_free(grisu_port);
g_free(grisu_soap_header);
g_free(grisu_cert_path);

/* soap */
soap_destroy(&soap); 
soap_end(&soap); 
soap_done(&soap); 
}

