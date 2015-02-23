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

void grid_application_set(const gchar *, const gchar *);
gchar *grid_application_get(const gchar *);
void grid_application_remove(const gchar *);
GList *grid_application_all(void);
GSList *grid_search_by_name(const gchar *);
gchar *grid_get_DN(gchar *);
const gchar *grid_get_status(void);
gchar *grid_random_alpha(gint);
gchar *grid_random_alphanum(gint);
gint grid_task_start(void);
void grid_task_stop(void);
void grid_job_start(gpointer);
void grid_job_stop(gpointer);
gint grid_job_new(const gchar *, gpointer);
gint grid_connect(gint);
gint grid_auth_check(void);
void grid_auth_set(gint);
void grid_credential_init(const gchar *);
gint grid_job_submit(const gchar *, gpointer);
GSList *grid_download_job(const gchar *);
void grid_process_job(const gchar *, GSList *);
// TODO grid_test -> grid_submit
void grid_test(const gchar *, const gchar *);
void grid_setup(void);
void grid_cleanup(void);
void grid_free(gpointer);
gpointer grid_new(void);

enum {GRID_MYPROXY, GRID_LOCALPROXY, GRID_SHIBBOLETH};
enum {GRID_SERIAL, GRID_MPI}; // other parallel types?

/* stuff in a job xml file - can be filled out by querying the server */
struct grid_pak
{
/* CURRENT - bit of a hack for MDS query/setup */
gint jobcode;
gchar *exename;
gchar *exe_version;

/* auth */
gchar *user_vo;

/* submission */
gchar *jobname;
gchar *remote_q;
gchar *remote_root;
gchar *remote_exe;
gchar *remote_exe_module;
gchar *remote_site;

/* processor specific */
gint remote_exe_type;
gint remote_exe_np;

/* source/output destination for the job */
gchar *local_cwd; // if NULL on download - ask user
gchar *local_input;
gchar *local_output;

/* TODO - local_upload_list if multiple data files need to be sent eg basis */
};

