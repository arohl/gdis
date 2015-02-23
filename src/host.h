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

gpointer host_new(gchar *);
void host_free(gpointer);
void host_free_all(void);
gint host_connect(gpointer);
void host_disconnect(gpointer);
gchar *host_talk(gpointer, const gchar *);
const gchar *host_name(gpointer);
const gchar *host_cwd(gpointer);

gint host_file_write(gpointer, gchar *, gchar *);
gint host_file_read(gpointer, gchar *, gchar *);

gpointer host_service_get(gpointer, const gchar *);
void host_service_foreach(gpointer, gpointer, gpointer);
void host_service_cycle(gpointer, const gchar *);
gint host_service_available(gpointer);
gint host_service_flags(gpointer);
gchar *host_service_fullpath(gpointer);

enum {SERVICE_BACKGROUND, SERVICE_MPI, SERVICE_QSUB,
      SERVICE_QSUB_MPI, SERVICE_PRIMARY, SERVICE_SECONDARY};

