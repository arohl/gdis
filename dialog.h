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

/* dialogs */
gpointer dialog_request(gint, gchar *, gpointer, gpointer, struct model_pak *);
gpointer dialog_window(gpointer);
gpointer dialog_model(gpointer);
void dialog_dump(void);
gint dialog_valid(gpointer);
gint dialog_exists(gint, struct model_pak *);
void dialog_refresh_all(void);
void dialog_refresh_type(gint);
void dialog_refresh_model(struct model_pak *);
void dialog_close_model(struct model_pak *);
void dialog_close(gint, struct model_pak *);
void dialog_destroy(GtkWidget *, gpointer);
void dialog_destroy_type(gint);
void dialog_destroy_model(struct model_pak *);
void dialog_destroy_single(gint, struct model_pak *);

gpointer dialog_child_get(gpointer, const gchar *);
void dialog_child_set(gpointer, const gchar *, gpointer);

void dialog_colour_new(GtkWidget *, gdouble *);

/* deprec */
/* cleanup control */
void file_cleanup(void);
void meas_cleanup(void);
void render_cleanup(void);
void symmetry_cleanup(void);
void task_cleanup(void);

/* individual dialogs */
void periodicity_dialog(void);
void surface_dialog(void);
void gulp_dialog(void);
void gamess_dialog(void);
void moldy_dialog(void);
void monty_dialog(void);
gint calculate_crystal_graph(struct model_pak * model);


