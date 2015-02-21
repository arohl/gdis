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

#include <math.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>
#include <unistd.h>
#ifdef __APPLE__
#include <OpenGL/gl.h>
#else
#include <GL/gl.h>
#endif

#include "gdis.h"
#include "file.h"
#include "parse.h"
#include "coords.h"
#include "matrix.h"
#include "molsurf.h"
#include "spatial.h"
#include "surface.h"
#include "sginfo.h"
#include "task.h"
#include "gui_shorts.h"
#include "interface.h"
#include "dialog.h"
#include "opengl.h"

/* main pak structures */
extern struct sysenv_pak sysenv;
extern struct elem_pak elements[];

/* molsurf globals */
gdouble ms_prad=1.0;
gdouble ms_blur=0.3;
gdouble ms_eden=0.1;
gint ms_method=MS_MOLECULAR, ms_colour=MS_TOUCH;

gpointer pulldown_colour;

GtkWidget *epot_vbox, *surf_epot_min, *surf_epot_max, *surf_epot_div;
GtkWidget *epot_pts;

/***************************************/
/* setup and run a surface calculation */
/***************************************/
#define DEBUG_CALC_MOLSURF 0
void ms_calculate(void)
{
gchar *text;
gdouble value;
const gchar *tmp;
struct model_pak *model;

model = sysenv.active_model;
g_assert(model != NULL);

/*
spatial_destroy_by_label("molsurf", model);
*/

/* get colouring type */
tmp = gtk_entry_get_text(GTK_ENTRY(pulldown_colour));

/* FIXME - can we rename De Josh? It's too close to Default and */
/* will cause the "default" to be "de" if the strcmp order is reversed */
if (g_ascii_strncasecmp(tmp,"De", 2) == 0)
  ms_colour = MS_DE;
if (g_ascii_strncasecmp(tmp,"Default", 7) == 0)
  ms_colour = MS_TOUCH;
if (g_ascii_strncasecmp(tmp,"AFM", 3) == 0)
  ms_colour = MS_AFM;
if (g_ascii_strncasecmp(tmp,"Electrostatic", 13) == 0)
  ms_colour = MS_EPOT;
if (g_ascii_strncasecmp(tmp,"Curvedness", 10) == 0)
  ms_colour = MS_CURVEDNESS;
if (g_ascii_strncasecmp(tmp,"Shape Index", 11) == 0)
  ms_colour = MS_SHAPE_INDEX;

/* get the scale */
if (ms_colour == MS_EPOT && !model->epot_autoscale)
  {
  model->epot_min = str_to_float(gtk_entry_get_text(GTK_ENTRY(surf_epot_min)));
  model->epot_max = str_to_float(gtk_entry_get_text(GTK_ENTRY(surf_epot_max)));
  model->epot_div = str_to_float(gtk_entry_get_text(GTK_ENTRY(surf_epot_div)));
  }
else
  {
  model->epot_min = G_MAXDOUBLE;
  model->epot_max = -G_MAXDOUBLE;
  }

/* main call */
value = ms_blur;
switch (ms_method)
  {
  case MS_EDEN:
  case MS_SSATOMS:
    value = ms_eden;
    break;
  }
ms_cube(value, ms_method, ms_colour, model);

/* update widget */
if (ms_colour == MS_EPOT && model->epot_autoscale)
  {
  text = g_strdup_printf("%f", model->epot_min);
  gtk_entry_set_text(GTK_ENTRY(surf_epot_min), text);
  g_free(text);

  text = g_strdup_printf("%f", model->epot_max);
  gtk_entry_set_text(GTK_ENTRY(surf_epot_max), text);
  g_free(text);

  text = g_strdup_printf("%d", model->epot_div);
  gtk_entry_set_text(GTK_ENTRY(surf_epot_div), text);
  g_free(text);
  }

coords_init(CENT_COORDS, model);

sysenv.refresh_dialog = TRUE;

redraw_canvas(SINGLE);
}

/************************/
/* simple deletion hook */
/************************/
void ms_delete(void)
{
struct model_pak *model;

model = sysenv.active_model;
g_assert(model != NULL);

/* remove any previous surfaces */
spatial_destroy_by_label("molsurf", model);
model->ms_colour_scale = FALSE;
coords_init(CENT_COORDS, model);
redraw_canvas(SINGLE);
}

/*******************************************/
/* Molecular surface colour mode selection */
/*******************************************/
void ms_iso_method(GtkWidget *entry)
{
const gchar *tmp;

g_assert(entry != NULL);

tmp = gtk_entry_get_text(GTK_ENTRY(entry));

ms_method = MS_MOLECULAR;
if (g_ascii_strncasecmp(tmp, "Electron dens", 13) == 0)
  ms_method = MS_EDEN;
if (g_ascii_strncasecmp(tmp, "Hirshfeld", 9) == 0)
  ms_method = MS_HIRSHFELD;
if (g_ascii_strncasecmp(tmp,"Promolecule", 11) == 0)
	ms_method = MS_SSATOMS;

redraw_canvas(ALL);
}

/*******************************************/
/* Molecular surface colour mode selection */
/*******************************************/
/*
void ms_colour_mode(GtkWidget *entry)
{
const gchar *tmp;

g_assert(entry != NULL);

tmp = gtk_entry_get_text(GTK_ENTRY(entry));

gtk_widget_set_sensitive(epot_vbox, FALSE);

if (g_ascii_strncasecmp(tmp,"Nearest atom", 12) == 0)
  ms_colour = MS_TOUCH;
if (g_ascii_strncasecmp(tmp,"AFM", 3) == 0)
  ms_colour = MS_AFM;
if (g_ascii_strncasecmp(tmp,"Electrostatic", 13) == 0)
  {
  ms_colour = MS_EPOT;
  gtk_widget_set_sensitive(epot_vbox, TRUE);
  }
if (g_ascii_strncasecmp(tmp,"Hirshfeld", 9) == 0)
  ms_colour = MS_HIRSHFELD;

redraw_canvas(ALL);
}
*/

/************************************************************/
/* callback to update electrostatic autoscaling sensitivity */
/************************************************************/
void gui_epot_scale_sensitive(GtkWidget *w, gpointer data)
{
if (gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(w)))
  gtk_widget_set_sensitive(epot_vbox, FALSE);
else
  gtk_widget_set_sensitive(epot_vbox, TRUE);
}

/***************************/
/* molecular surface setup */
/***************************/
void gui_isosurf_dialog()
{
gchar *text;
gpointer dialog, ptr;
GList *list;
GtkWidget *window, *frame, *vbox, *vbox2, *hbox, *label;
struct model_pak *data;

/* checks */
data = sysenv.active_model;
if (!data)
  return;
if (data->id == MORPH)
  return;

/* create dialog */
dialog = dialog_request(SURF, "Iso-Surfaces", NULL, NULL, data);
if (!dialog)
  return;
window = dialog_window(dialog);

/* isosurface type */
frame = gtk_frame_new(NULL);
gtk_box_pack_start(GTK_BOX(GTK_DIALOG(window)->vbox),frame,FALSE,FALSE,0); 
gtk_container_set_border_width(GTK_CONTAINER(frame), PANEL_SPACING);
vbox = gtk_vbox_new(FALSE, PANEL_SPACING);
gtk_container_add(GTK_CONTAINER(frame), vbox);
hbox = gtk_hbox_new(FALSE, 0);
gtk_box_pack_start(GTK_BOX(vbox), hbox, FALSE, TRUE, 0);

/* combo box */
list = NULL;
list = g_list_prepend(list, "Molecular surface");
list = g_list_prepend(list, "Hirshfeld surface");
list = g_list_prepend(list, "Electron density");
list = g_list_prepend(list, "Promolecule isosurface");
list = g_list_reverse(list);

ptr = gui_pulldown_new("Iso-surface type", list, FALSE, hbox);

g_signal_connect(GTK_OBJECT(ptr), "changed", GTK_SIGNAL_FUNC(ms_iso_method), NULL);

/* colour mode */
hbox = gtk_hbox_new(FALSE, 0);
gtk_box_pack_start(GTK_BOX(vbox), hbox, FALSE, TRUE, 0);

/* NB: Default != Nearest atom - eg electron density */
list = NULL;
list = g_list_prepend(list, "Default");
list = g_list_prepend(list, "AFM");
list = g_list_prepend(list, "Electrostatic");
list = g_list_prepend(list, "Curvedness");
list = g_list_prepend(list, "Shape Index");
list = g_list_prepend(list, "De");
list = g_list_reverse(list);
pulldown_colour = gui_pulldown_new("Colour method", list, FALSE, hbox);
/* redo when colour mode changes can be done without recalculating */
/*
g_signal_connect(GTK_OBJECT(GTK_COMBO(ms_colour_combo)->entry), "changed", 
                 GTK_SIGNAL_FUNC(ms_colour_mode), NULL);
*/

/* frame for spinner setup */
frame = gtk_frame_new(NULL);
gtk_box_pack_start(GTK_BOX(GTK_DIALOG(window)->vbox),frame,FALSE,FALSE,0); 
vbox = gtk_vbox_new(FALSE, PANEL_SPACING);
gtk_container_add(GTK_CONTAINER(frame), vbox);
gtk_container_set_border_width(GTK_CONTAINER(frame), PANEL_SPACING);

gui_direct_spin("Triangulation grid size",
                  &sysenv.render.ms_grid_size, 0.05, 10.0, 0.05, NULL, NULL, vbox);

gui_direct_spin("Molecular surface blurring", &ms_blur, 0.05, 1.0, 0.05, NULL, NULL, vbox);
gui_direct_spin("Electron density value", &ms_eden, 0.001, 1.0, 0.001, NULL, NULL, vbox);

/* electrostatic potential scale setup */
frame = gtk_frame_new(NULL);
gtk_box_pack_start(GTK_BOX(vbox), frame, FALSE, TRUE,0);

vbox2 = gtk_vbox_new(FALSE, 0);
gtk_container_add(GTK_CONTAINER(frame), vbox2);

gui_auto_check("Electrostatic autoscaling", gui_epot_scale_sensitive, NULL, &data->epot_autoscale, vbox2);

epot_vbox = gtk_vbox_new(FALSE, 0);
gtk_box_pack_start(GTK_BOX(vbox2), epot_vbox, FALSE, TRUE,0);

hbox = gtk_hbox_new(FALSE, 0);
gtk_box_pack_start(GTK_BOX(epot_vbox), hbox, FALSE, TRUE,0);
label = gtk_label_new("maximum ");
gtk_box_pack_start(GTK_BOX(hbox), label, FALSE, FALSE, 0);
surf_epot_max = gtk_entry_new();
gtk_box_pack_end(GTK_BOX(hbox), surf_epot_max, FALSE, FALSE, 0);
text = g_strdup_printf("%f", data->epot_max);
gtk_entry_set_text(GTK_ENTRY(surf_epot_max), text);
g_free(text);

hbox = gtk_hbox_new(FALSE, 0);
gtk_box_pack_start(GTK_BOX(epot_vbox), hbox, FALSE, TRUE,0);
label = gtk_label_new("minimum ");
gtk_box_pack_start(GTK_BOX(hbox), label, FALSE, FALSE, 0);
surf_epot_min = gtk_entry_new();
gtk_box_pack_end(GTK_BOX(hbox), surf_epot_min, FALSE, FALSE, 0);
text = g_strdup_printf("%f", data->epot_min);
gtk_entry_set_text(GTK_ENTRY(surf_epot_min), text);
g_free(text);

hbox = gtk_hbox_new(FALSE, 0);
gtk_box_pack_start(GTK_BOX(epot_vbox), hbox, FALSE, TRUE,0);
label = gtk_label_new("divisions ");
gtk_box_pack_start(GTK_BOX(hbox), label, FALSE, FALSE, 0);
surf_epot_div = gtk_entry_new();
gtk_box_pack_end(GTK_BOX(hbox), surf_epot_div, FALSE, FALSE, 0);
text = g_strdup_printf("%d", data->epot_div);
gtk_entry_set_text(GTK_ENTRY(surf_epot_div), text);
g_free(text);

/* make, hide, close - terminating buttons */
gui_stock_button(GTK_STOCK_EXECUTE, ms_calculate, NULL,
                   GTK_DIALOG(window)->action_area);

gui_stock_button(GTK_STOCK_REMOVE, ms_delete, NULL,
                   GTK_DIALOG(window)->action_area);

gui_stock_button(GTK_STOCK_CLOSE, dialog_destroy, dialog,
                   GTK_DIALOG(window)->action_area);

/* done */
gtk_widget_show_all(window);
gtk_widget_set_sensitive(epot_vbox, FALSE);
}

