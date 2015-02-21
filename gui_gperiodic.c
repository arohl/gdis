/*
Copyright (C) 2000 by Sean David Fleming

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
#include <math.h>
#include <string.h>
#include <unistd.h>
#include <assert.h>
#include <errno.h>
#include <gtk/gtk.h>

#include "gdis.h"
#include "coords.h"
#include "matrix.h"
#include "render.h"
#include "gperiodic.h"
#include "gui_shorts.h"
#include "interface.h"
#include "dialog.h"
#include "opengl.h"

extern struct sysenv_pak sysenv;
extern struct elem_pak elements[];
extern GtkWidget *window;

/* currently modifiable elem  */
/* FIXME - since it's a global, only one dialog can be allowed at a time */
struct elem_pak elem_edit;
GtkWidget *elem_edit_colour;

/***********************************/
/* update element database locally */
/***********************************/
void elem_change_current(GtkWidget *w, struct elem_pak *element)
{
struct model_pak *model;

model = sysenv.active_model;
if (model)
  {
  put_elem_data(element, model);
  model_colour_scheme(model->colour_scheme, model); 
  connect_refresh(model);
  redraw_canvas(SINGLE);
  }
}

/************************************/
/* update element database globally */
/************************************/
void elem_change_global(GtkWidget *w, struct elem_pak *element)
{
struct model_pak *model;

put_elem_data(element, NULL);
refresh_table();

model = sysenv.active_model;
if (model)
  {
  put_elem_data(element, model);
  connect_refresh(model);
  model_colour_scheme(model->colour_scheme, model); 
  redraw_canvas(SINGLE);
  }
}

/******************************/
/* reset the element database */
/******************************/
void elem_reset(GtkWidget *w, struct elem_pak *element)
{
GSList *list1, *list2;
GdkColor colour;
struct elem_pak *elem;
struct model_pak *model;

/* update global database */
list1 = sysenv.elements;
while (list1 != NULL)
  {
  elem = list1->data;
  list1 = g_slist_next(list1);
  if (elem->number == element->number)
    sysenv.elements = g_slist_remove(sysenv.elements, elem);
  }
refresh_table();

/* update model database */
for (list1=sysenv.mal ; list1 ; list1=g_slist_next(list1))
  {
  model = list1->data;
  list2 = model->elements;
  while (list2)
    {
    elem = list2->data;
    list2 = g_slist_next(list2);
    if (elem->number == element->number)
      model->elements = g_slist_remove(model->elements, elem);
    }
  }

/* update model display */
model = sysenv.active_model;
if (model)
  {
  connect_refresh(model);
  model_colour_scheme(model->colour_scheme, model); 
  redraw_canvas(SINGLE);
  }

/* widget updates */
memcpy(&elem_edit, &elements[elem_edit.number], sizeof(struct elem_pak));
gui_relation_update_widget(&elem_edit.cova);
gui_relation_update_widget(&elem_edit.vdw);
gui_relation_update_widget(&elem_edit.charge);

colour.red   = elem_edit.colour[0]*65535.0;
colour.green = elem_edit.colour[1]*65535.0;
colour.blue  = elem_edit.colour[2]*65535.0;
gtk_widget_modify_bg(elem_edit_colour, GTK_STATE_NORMAL, &colour);
}

/***************************/
/* Single element dialogue */
/***************************/
void display_element_dialog(GtkWidget *w, gint i)
{
gchar *text;
gpointer dialog;
GtkWidget *window, *frame, *vbox, *hbox, *label;
struct table_entry *entry;

entry = (struct table_entry *) &table[i];

/* get global elem data */
get_elem_data(i+1, &elem_edit, NULL);

/* NEW - close any other element dialog */
dialog_destroy_type(ELEM_EDIT);

/* retrieve a suitable dialog */
dialog = dialog_request(ELEM_EDIT, "Element edit", NULL, NULL, NULL);
if (!dialog)
  return;
window = dialog_window(dialog);

/* top frame */
frame = gtk_frame_new(NULL);
gtk_box_pack_start(GTK_BOX(GTK_DIALOG(window)->vbox), frame, TRUE, FALSE, 0);
gtk_container_set_border_width(GTK_CONTAINER(frame), PANEL_SPACING);
vbox = gtk_vbox_new(TRUE,1);
gtk_container_add(GTK_CONTAINER(frame), vbox);
gtk_container_set_border_width(GTK_CONTAINER(vbox), PANEL_SPACING);

hbox = gtk_hbox_new(FALSE,0);
gtk_box_pack_start(GTK_BOX(vbox), hbox, FALSE, FALSE, 0);
label = gtk_label_new("Name:");
gtk_box_pack_start(GTK_BOX(hbox), label, FALSE, FALSE, 0);
label = gtk_label_new(elem_edit.name);
gtk_box_pack_end(GTK_BOX(hbox), label, FALSE, FALSE, 0);

hbox = gtk_hbox_new(FALSE,0);
gtk_box_pack_start(GTK_BOX(vbox), hbox, FALSE, FALSE, 0);
label = gtk_label_new("Symbol:");
gtk_box_pack_start(GTK_BOX(hbox), label, FALSE, FALSE, 0);
label = gtk_label_new(elem_edit.symbol);
gtk_box_pack_end(GTK_BOX(hbox), label, FALSE, FALSE, 0);

hbox = gtk_hbox_new(FALSE,0);
gtk_box_pack_start(GTK_BOX(vbox), hbox, FALSE, FALSE, 0);
label = gtk_label_new("Number:");
gtk_box_pack_start(GTK_BOX(hbox), label, FALSE, FALSE, 0);
text = g_strdup_printf("%d", elem_edit.number);
label = gtk_label_new(text);
g_free(text);
gtk_box_pack_end(GTK_BOX(hbox), label, FALSE, FALSE, 0);

hbox = gtk_hbox_new(FALSE,0);
gtk_box_pack_start(GTK_BOX(vbox), hbox, FALSE, FALSE, 0);
label = gtk_label_new("Weight:");
gtk_box_pack_start(GTK_BOX(hbox), label, FALSE, FALSE, 0);
text = g_strdup_printf("%8.4f", elem_edit.weight);
label = gtk_label_new(text);
g_free(text);
gtk_box_pack_end(GTK_BOX(hbox), label, FALSE, FALSE, 0);

hbox = gtk_hbox_new(FALSE,0);
gtk_box_pack_start(GTK_BOX(vbox), hbox, FALSE, FALSE, 0);
label = gtk_label_new("Ionic/Covalent radius:    ");
gtk_box_pack_start(GTK_BOX(hbox), label, FALSE, FALSE, 0);
gui_direct_spin(NULL, &elem_edit.cova, 0.1, 3.0, 0.01, NULL, NULL, hbox);

hbox = gtk_hbox_new(FALSE,0);
gtk_box_pack_start(GTK_BOX(vbox), hbox, FALSE, FALSE, 0);
label = gtk_label_new("VdW radius:");
gtk_box_pack_start(GTK_BOX(hbox), label, FALSE, FALSE, 0);
gui_direct_spin(NULL, &elem_edit.vdw, 0.1, 4.0, 0.01, NULL, NULL, hbox);

gui_colour_box("Colour: ", elem_edit.colour, vbox);

/* application buttons/check boxes */
frame = gtk_frame_new(NULL);
gtk_box_pack_start(GTK_BOX(GTK_DIALOG(window)->vbox), frame, TRUE, FALSE, 0);
gtk_container_set_border_width(GTK_CONTAINER(frame), PANEL_SPACING);

/* two column element data display */
vbox = gtk_vbox_new(FALSE,PANEL_SPACING);
gtk_container_add(GTK_CONTAINER(frame), vbox);
gtk_container_set_border_width(GTK_CONTAINER(vbox), PANEL_SPACING);
gui_button_x("Apply to current model", elem_change_current, &elem_edit, vbox);
gui_button_x("Apply globally", elem_change_global, &elem_edit, vbox);
gui_button_x("Reset values", elem_reset, &elem_edit, vbox);

gui_stock_button(GTK_STOCK_CLOSE, dialog_destroy, dialog,
                   GTK_DIALOG(window)->action_area);

gtk_widget_show_all(window);
}

/***********************************/
/* on the fly table colour refresh */
/***********************************/
void refresh_table(void)
{
gint i;
GdkColor colour;
struct elem_pak element;

/* checks */
if (!dialog_exists(GPERIODIC, NULL))
  return;

for(i=0 ; i<sizeof(table) ; i++) 
  {
/* stop if no more element data (NB: gperiodic/gdis index mismatch) */
  if(i >= sysenv.num_elements-1) 
    break;
/* stop if no more element data (NB: gperiodic/gdis index mismatch) */
  if (get_elem_data(i+1, &element, NULL))
    break;
/* jump the gaps */
  if (GTK_IS_WIDGET(table[i].button)) 
    {
    colour.red   = element.colour[0]*65535.0;
    colour.green = element.colour[1]*65535.0;
    colour.blue  = element.colour[2]*65535.0;
    gtk_widget_modify_bg(table[i].button, GTK_STATE_NORMAL, &colour);
    }
  }
}

/*****************************************************************/
/* update the current model with values from the global database */
/*****************************************************************/
void refresh_model_from_table(GtkWidget *w, gpointer dummy)
{
struct model_pak *model;

model = sysenv.active_model;
if (model)
  {
/* TODO - refresh bonding/connectivity & other data? */
/* refresh colour */
  model_colour_scheme(model->colour_scheme, model); 
  redraw_canvas(SINGLE);
  }
}

/******************************/
/* Construct a periodic table */
/******************************/
void gui_gperiodic_dialog(void)
{
gint i;
gchar *text=NULL;
gpointer dialog;
GtkWidget *w, *label, *vbox, *periodic_table;
struct elem_pak elem;
GdkColor colour;

/* retrieve a suitable dialog */
dialog = dialog_request(GPERIODIC, "GPeriodic", NULL, NULL, NULL);
if (!dialog)
  return;
w = dialog_window(dialog);

/* background colour for periodic table */
colour.red = 65535;
colour.green = 65535;
colour.blue = 65535;
gtk_widget_modify_bg(w, GTK_STATE_NORMAL, &colour);

/* use a vbox for the menubar and the table of elements... */
vbox = gtk_vbox_new(FALSE, PANEL_SPACING);
gtk_container_border_width(GTK_CONTAINER(vbox), PANEL_SPACING);
gtk_container_add(GTK_CONTAINER(GTK_DIALOG(w)->vbox), vbox);

/* nice heading */
label = gtk_label_new("Periodic Table of the Elements");
gtk_box_pack_start(GTK_BOX(vbox), label, FALSE, FALSE, 0);

/* create the table widget to hold the periodic table */
periodic_table = gtk_table_new(11, 18, TRUE);
gtk_box_pack_end(GTK_BOX(vbox), periodic_table, TRUE, TRUE, 0);

/* now for each element in the table of elements, create a display */
/* item for it, and add it to the table... */
for (i=0 ; i<sizeof(table) ; i++ ) 
  {
/* stop if no more element data (NB: gperiodic/gdis index mismatch) */
  if (get_elem_data(i+1, &elem, NULL))
    break;

/* create the button */
  table[i].button = gtk_button_new_with_label(elem.symbol);
  colour.red   = elem.colour[0]*65535.0;
  colour.green = elem.colour[1]*65535.0;
  colour.blue  = elem.colour[2]*65535.0;
  gtk_widget_modify_bg(table[i].button, GTK_STATE_NORMAL, &colour);

/* set up a string for the tooltips */
  text = g_strdup_printf("%s  n:%d", elem.symbol, elem.number);

/* create a new tooltips object... */
  table[i].tooltip = gtk_tooltips_new();
  gtk_tooltips_set_delay(table[i].tooltip,100);
  gtk_tooltips_set_tip(table[i].tooltip,table[i].button,text,NULL);
  g_free(text);

/* attach the button to the table */
  gtk_table_attach_defaults(GTK_TABLE(periodic_table), table[i].button,
               table[i].x - 1, table[i].x, table[i].y - 1, table[i].y);

/* connect the destroy method to it */
  g_signal_connect(GTK_OBJECT(table[i].button), "clicked",
                   GTK_SIGNAL_FUNC(display_element_dialog), GINT_TO_POINTER(i));
  }

gtk_table_set_row_spacings(GTK_TABLE(periodic_table), PANEL_SPACING);
gtk_table_set_col_spacings(GTK_TABLE(periodic_table), PANEL_SPACING);

/* terminating buttons */
gui_stock_button(GTK_STOCK_REFRESH, refresh_model_from_table, NULL, 
                   GTK_DIALOG(w)->action_area);

gui_stock_button(GTK_STOCK_CLOSE, dialog_destroy, dialog,
                   GTK_DIALOG(w)->action_area);

gtk_widget_show_all(w);
}
