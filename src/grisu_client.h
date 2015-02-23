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

int grisu_init(void);
void grisu_stop(void);
void grisu_auth_set(gboolean);
void grisu_keypair_set(const gchar *, const gchar *);
void grisu_vo_set(const gchar *);
void grisu_ws_set(const gchar *);
void grisu_myproxy_set(const gchar *, const gchar *);
const gchar *grisu_username_get(void);
const gchar *grisu_password_get(void);
const gchar *grisu_vo_get(void);
const gchar *grisu_ws_get(void);
const gchar *grisu_myproxy_get(void);
const gchar *grisu_myproxyport_get(void);
GHashTable *grisu_application_details(const gchar *, const gchar *);
GSList *grisu_application_versions(const gchar *, const gchar *);
GSList *grisu_site_list(void);
GSList *grisu_submit_find(const gchar *);
GSList *grisu_fqans_get(void);
GSList *grisu_job_names(void);
gchar *grisu_site_name(const gchar *);
gchar *grisu_jobname_request(gint);
gint grisu_job_type(const gchar *);
gint grisu_job_configure(gchar *, gchar *);
gchar *grisu_job_create(const gchar *);
gchar *grisu_job_submit(const gchar *, const gchar *);
gchar *grisu_job_cwd(const gchar *);
gchar *grisu_xml_create(gpointer);
const gchar *grisu_job_status(const gchar *);
gchar *grisu_job_details(const gchar *);
void grisu_job_remove(const gchar *);
gchar *grisu_ps(void);
GSList *grisu_file_list(const gchar *);
gchar *grisu_absolute_job_dir(const gchar *, const gchar *, const gchar *);
gchar *grisu_relative_job_dir(const gchar *);
gint grisu_file_upload(const gchar *, const gchar *);
gint grisu_file_download(const gchar *, const gchar *);

