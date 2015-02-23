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
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "gdis.h"
#include "defect.h"
#include "dialog.h"
#include "matrix.h"
#include "gui_shorts.h"
#include "interface.h"

/* globals */
extern struct sysenv_pak sysenv;
struct defect_pak defect;

/**********************/
/* setup gui defaults */
/**********************/
void gui_defect_default(void)
{
defect.cleave = FALSE;
defect.neutral = FALSE;
defect.cluster = FALSE;
VEC2SET(defect.region, 10, 10);
VEC3SET(defect.orient, 1, 0, 0);
VEC3SET(defect.burgers, 1, 0, 0);
VEC2SET(defect.origin, 0, 0);
VEC2SET(defect.center, 0, 0);
}

/*****************************************/
/* callback to start the defect building */
/*****************************************/
void gui_defect_build(void)
{
struct model_pak *model;

/* checks and setup */
model = sysenv.active_model;
if (!model)
  return;
if (model->periodic != 3)
  return;

/* TODO - if the building is long - run in background */
defect_new(&defect, model);
}

/**************************/
/* update on model switch */
/**************************/
void gui_defect_refresh(gpointer data)
{
}

/**********************************/
/* neutralization scheme callback */
/**********************************/
void gui_defect_neutralization(GtkWidget *entry)
{
const gchar *text;

g_assert(entry != NULL);
g_assert(GTK_IS_ENTRY(entry));

text = gtk_entry_get_text(GTK_ENTRY(entry));
}

/************************************/
/* setup the defect creation dialog */
/************************************/
void gui_defect_dialog(void)
{
gpointer dialog;
GList *list;
GtkWidget *window, *vbox, *hbox, *table, *label, *spin, *combo;

/* create new dialog */
/* TODO - get rid of the need for a unique code in this call */
/* replace with lookups based on the name & single/multiple instance allowed */
dialog = dialog_request(100, "Dislocation builder", gui_defect_refresh, NULL, NULL);
if (!dialog)
  return;
window = dialog_window(dialog);

/* init */
gui_defect_default();

/* defect build setup frame */
vbox = gui_frame_vbox("Geometry", FALSE, FALSE, GTK_DIALOG(window)->vbox);
table = gtk_table_new(4, 5, FALSE);
gtk_box_pack_start(GTK_BOX(vbox), table, TRUE, TRUE, 0); 

label = gtk_label_new("Orientation vector ");
gtk_misc_set_alignment(GTK_MISC(label), 0, 0.5);
gtk_table_attach_defaults(GTK_TABLE(table),label,0,1,0,1);
spin = gui_direct_spin(NULL, &defect.orient[0], -20, 20, 1, NULL, NULL, NULL);
gtk_table_attach_defaults(GTK_TABLE(table),spin,1,2,0,1);
spin = gui_direct_spin(NULL, &defect.orient[1], -20, 20, 1, NULL, NULL, NULL);
gtk_table_attach_defaults(GTK_TABLE(table),spin,2,3,0,1);
spin = gui_direct_spin(NULL, &defect.orient[2], -20, 20, 1, NULL, NULL, NULL);
gtk_table_attach_defaults(GTK_TABLE(table),spin,3,4,0,1);

label = gtk_label_new("Burgers vector ");
gtk_misc_set_alignment(GTK_MISC(label), 0, 0.5);
gtk_table_attach_defaults(GTK_TABLE(table),label,0,1,1,2);
spin = gui_direct_spin(NULL, &defect.burgers[0], -20, 20, 0.05, NULL, NULL, NULL);
gtk_table_attach_defaults(GTK_TABLE(table),spin,1,2,1,2);
spin = gui_direct_spin(NULL, &defect.burgers[1], -20, 20, 0.05, NULL, NULL, NULL);
gtk_table_attach_defaults(GTK_TABLE(table),spin,2,3,1,2);
spin = gui_direct_spin(NULL, &defect.burgers[2], -20, 20, 0.05, NULL, NULL, NULL);
gtk_table_attach_defaults(GTK_TABLE(table),spin,3,4,1,2);

label = gtk_label_new("Defect origin ");
gtk_misc_set_alignment(GTK_MISC(label), 0, 0.5);
gtk_table_attach_defaults(GTK_TABLE(table),label,0,1,2,3);
spin = gui_direct_spin(NULL, &defect.origin[0], 0.0, 1.0, 0.05, NULL, NULL, NULL);
gtk_table_attach_defaults(GTK_TABLE(table),spin,2,3,2,3);
spin = gui_direct_spin(NULL, &defect.origin[1], 0.0, 1.0, 0.05, NULL, NULL, NULL);
gtk_table_attach_defaults(GTK_TABLE(table),spin,3,4,2,3);

label = gtk_label_new("Defect center ");
gtk_misc_set_alignment(GTK_MISC(label), 0, 0.5);
gtk_table_attach_defaults(GTK_TABLE(table),label,0,1,3,4);
spin = gui_direct_spin(NULL, &defect.center[0], 0.0, 1.0, 0.05, NULL, NULL, NULL);
gtk_table_attach_defaults(GTK_TABLE(table),spin,2,3,3,4);
spin = gui_direct_spin(NULL, &defect.center[1], 0.0, 1.0, 0.05, NULL, NULL, NULL);
gtk_table_attach_defaults(GTK_TABLE(table),spin,3,4,3,4);

label = gtk_label_new("Region sizes ");
gtk_misc_set_alignment(GTK_MISC(label), 0, 0.5);
gtk_table_attach_defaults(GTK_TABLE(table),label,0,1,4,5);
spin = gui_direct_spin(NULL, &defect.region[0], 0, 999, 1, NULL, NULL, NULL);
gtk_table_attach_defaults(GTK_TABLE(table),spin,2,3,4,5);
spin = gui_direct_spin(NULL, &defect.region[1], 0, 999, 1, NULL, NULL, NULL);
gtk_table_attach_defaults(GTK_TABLE(table),spin,3,4,4,5);

/* construct setup frame */
vbox = gui_frame_vbox("Options", FALSE, FALSE, GTK_DIALOG(window)->vbox);

hbox = gtk_hbox_new(FALSE, 0);
gtk_box_pack_start(GTK_BOX(vbox), hbox, TRUE, TRUE, 0); 
label = gtk_label_new("Charge neutralization ");
gtk_box_pack_start(GTK_BOX(hbox), label, FALSE, FALSE, 0); 
list = NULL;
list = g_list_prepend(list, "none");
combo = gtk_combo_new();
gtk_entry_set_editable(GTK_ENTRY(GTK_COMBO(combo)->entry), FALSE);
gtk_combo_set_popdown_strings(GTK_COMBO(combo), list);
g_signal_connect(GTK_OBJECT(GTK_COMBO(combo)->entry), "changed", 
                 GTK_SIGNAL_FUNC(gui_defect_neutralization), NULL);
gtk_box_pack_end(GTK_BOX(hbox), combo, FALSE, FALSE, 0); 

gui_direct_check("Allow bond cleaving ", &defect.cleave, NULL, NULL, vbox);
gui_direct_check("Discard periodicity ", &defect.cluster, NULL, NULL, vbox);

/* terminating buttons */
gui_stock_button(GTK_STOCK_EXECUTE, gui_defect_build, NULL,
                   GTK_DIALOG(window)->action_area);

gui_stock_button(GTK_STOCK_CLOSE, dialog_destroy, dialog,
                   GTK_DIALOG(window)->action_area);

gtk_widget_show_all(window);
}


