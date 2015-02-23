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
#include "coords.h"
#include "file.h"
#include "graph.h"
#include "analysis.h"
#include "dialog.h"
#include "interface.h"
#include "gui_shorts.h"

extern struct elem_pak elements[];
extern struct sysenv_pak sysenv;

/******************************/
/* export graphs as text data */
/******************************/
void analysis_export(gchar *name)
{
struct model_pak *model;

model = sysenv.active_model;
g_assert(model != NULL);

dialog_destroy_type(FILE_SELECT);

if (model->graph_active)
  {
  graph_write(name, model->graph_active);
  gui_text_show(STANDARD, "Successfully exported graph.\n");
  }
else
  gui_text_show(WARNING, "Current model is not a graph.\n");
}

/******************************/
/* export graph file selector */
/******************************/
void analysis_export_dialog(void)
{
/* FIXME - dialog + model are linked */
if (sysenv.active_model)
  file_dialog("Export graph as text", NULL, FILE_SAVE, 
              (gpointer) analysis_export, MD_ANALYSIS);
}

/**********************/
/* submit an rdf task */
/**********************/
void exec_analysis_task(GtkWidget *w, gpointer dialog)
{
const gchar *type;
gpointer gui_menu_calc;
struct analysis_pak *analysis;
struct model_pak *model;

/* duplicate the analysis structure at time of request */
g_assert(dialog != NULL);
model = dialog_model(dialog);
g_assert(model != NULL);
analysis = analysis_dup(model->analysis);

/* get the task type requested from the dialog */
gui_menu_calc = dialog_child_get(dialog, "gui_menu_calc");
type = gui_pulldown_text(gui_menu_calc);

if (g_ascii_strncasecmp(type, "RDF", 3) == 0)
  {
  analysis->rdf_normalize = TRUE;
  task_new("RDF", &analysis_plot_rdf, analysis, &analysis_show, model, model);
  return;
  }
if (g_ascii_strncasecmp(type, "Pair", 4) == 0)
  {
  analysis->rdf_normalize = FALSE;
  task_new("Pair", &analysis_plot_rdf, analysis, &analysis_show, model, model);
  return;
  }
if (g_ascii_strncasecmp(type, "VACF", 4) == 0)
  {
  task_new("VACF", &analysis_plot_vacf, analysis, &analysis_show, model, model);
  return;
  }
if (g_ascii_strncasecmp(type, "Temp", 4) == 0)
  {
  task_new("Temp", &analysis_plot_temp, analysis, &analysis_show, model, model);
  return;
  }
if (g_ascii_strncasecmp(type, "Pot", 3) == 0)
  {
  task_new("Energy", &analysis_plot_pe, analysis, &analysis_show, model, model);
  return;
  }
if (g_ascii_strncasecmp(type, "Kin", 3) == 0)
  {
  task_new("Energy", &analysis_plot_ke, analysis, &analysis_show, model, model);
  return;
  }
if (g_ascii_strncasecmp(type, "Meas", 4) == 0)
  {
  task_new("Measure", &analysis_plot_meas, analysis, &analysis_show, model, model);
  return;
  }
}

/**********************/
/* atom entry changed */
/**********************/
void md_atom_changed(GtkWidget *w, gchar **atom)
{
g_assert(w != NULL);

/* NB: we pass a pointer to the character pointer */
if (*atom)
  g_free(*atom);

*atom = g_strdup(gtk_entry_get_text(GTK_ENTRY(w)));
}

/*******************************/
/* multi-frame analysis dialog */
/*******************************/
void gui_analysis_dialog(void)
{
gpointer dialog;
GtkWidget *window, *combo, *label;
GtkWidget *frame, *main_hbox, *main_vbox, *hbox, *vbox;
GtkWidget *gui_menu_calc, *gui_vbox_rdf;
GSList *list, *item;
GList *match_list=NULL, *calc_list=NULL;
struct model_pak *model;
struct analysis_pak *analysis;

/* checks and setup */
model = sysenv.active_model;
if (!model)
  return;
if (!model->animation)
  return;
if (analysis_init(model))
  return;

analysis = model->analysis;

/* create new dialog */
dialog = dialog_request(MD_ANALYSIS, "MD analysis", NULL, NULL, model);
if (!dialog)
  return;
window = dialog_window(dialog);

/* --- main box */
main_hbox = gtk_hbox_new(TRUE, 0);
gtk_box_pack_start(GTK_BOX(GTK_DIALOG(window)->vbox), main_hbox, TRUE, TRUE, 0);

/* --- left pane */
main_vbox = gtk_vbox_new(FALSE, 0);
gtk_box_pack_start(GTK_BOX(main_hbox), main_vbox, TRUE, TRUE, 0);

/* --- analysis dialogs are specific to the active model when initiated */
frame = gtk_frame_new("Model");
gtk_box_pack_start(GTK_BOX(main_vbox), frame, FALSE, FALSE, 0);
gtk_container_set_border_width(GTK_CONTAINER(frame), PANEL_SPACING);
hbox = gtk_hbox_new(TRUE, PANEL_SPACING);
gtk_container_add(GTK_CONTAINER(frame), hbox);
label = gtk_label_new(model->basename);
gtk_box_pack_start(GTK_BOX(hbox), label, FALSE, FALSE, 0);

/* --- calculation menu */
frame = gtk_frame_new("Calculate");
gtk_box_pack_start(GTK_BOX(main_vbox), frame, FALSE, FALSE, 0);
gtk_container_set_border_width(GTK_CONTAINER(frame), PANEL_SPACING);
vbox = gtk_vbox_new(TRUE, PANEL_SPACING);
gtk_container_add(GTK_CONTAINER(frame), vbox);

/* --- menu items */
if (g_file_test(model->gulp.trj_file, G_FILE_TEST_EXISTS))
  {
  calc_list = g_list_prepend(calc_list, "Potential E");
  calc_list = g_list_prepend(calc_list, "Kinetic E");
  calc_list = g_list_prepend(calc_list, "Temperature");
  calc_list = g_list_prepend(calc_list, "VACF");
  }
else
  {
  calc_list = g_list_prepend(calc_list, "Measurements");
  }
calc_list = g_list_prepend(calc_list, "RDF");
calc_list = g_list_prepend(calc_list, "Pair count");
gui_menu_calc = gui_pulldown_new("Perform ", calc_list, FALSE, vbox);
dialog_child_set(dialog, "gui_menu_calc", gui_menu_calc);

/* --- right pane */
gui_vbox_rdf = gtk_vbox_new(FALSE, 0);
dialog_child_set(dialog, "gui_vbox_rdf", gui_vbox_rdf);
gtk_box_pack_start(GTK_BOX(main_hbox), gui_vbox_rdf, TRUE, TRUE, 0);

/* --- interval setup */
frame = gtk_frame_new("Analysis Interval");
gtk_box_pack_start(GTK_BOX(gui_vbox_rdf), frame, FALSE, FALSE, 0);
gtk_container_set_border_width(GTK_CONTAINER(frame), PANEL_SPACING);
vbox = gtk_vbox_new(TRUE, PANEL_SPACING);
gtk_container_add(GTK_CONTAINER(frame), vbox);

/* TODO - some warning about the RDF being valid only for < L/2? */
/* see Allan & Tildesley pg 199 */
gui_direct_spin("Start", &analysis->start, 0.0, 10.0*model->rmax, 0.1, NULL, NULL, vbox);
gui_direct_spin("Stop", &analysis->stop, 0.1, 10.0*model->rmax, 0.1, NULL, NULL, vbox);
gui_direct_spin("Step", &analysis->step, 0.1, model->rmax, 0.1, NULL, NULL, vbox);

/* --- RDF atom setup */
frame = gtk_frame_new("Analysis Atoms");
gtk_box_pack_start(GTK_BOX(gui_vbox_rdf), frame, FALSE, FALSE, 0);
gtk_container_set_border_width(GTK_CONTAINER(frame), PANEL_SPACING);
vbox = gtk_vbox_new(TRUE, PANEL_SPACING);
gtk_container_add(GTK_CONTAINER(frame), vbox);

/* setup unique element label list */
list = find_unique(LABEL, model);
match_list = g_list_append(match_list, "Any");
for (item=list ; item  ; item=g_slist_next(item))
  {
  match_list = g_list_append(match_list, item->data);
  }
g_slist_free(list);

/* match 1 */
combo = gtk_combo_new();
gtk_combo_set_popdown_strings(GTK_COMBO(combo), match_list);
gtk_entry_set_editable(GTK_ENTRY(GTK_COMBO(combo)->entry), TRUE);
gtk_box_pack_start(GTK_BOX(vbox), combo, TRUE, TRUE, 0);
g_signal_connect(GTK_OBJECT(GTK_COMBO(combo)->entry), "changed",
                 GTK_SIGNAL_FUNC(md_atom_changed), (gpointer) &analysis->atom1);

/* match 2 */
combo = gtk_combo_new();
gtk_combo_set_popdown_strings(GTK_COMBO(combo), match_list);
gtk_entry_set_editable(GTK_ENTRY(GTK_COMBO(combo)->entry), TRUE);
gtk_box_pack_start(GTK_BOX(vbox), combo, TRUE, TRUE, 0);
g_signal_connect(GTK_OBJECT(GTK_COMBO(combo)->entry), "changed",
                 GTK_SIGNAL_FUNC(md_atom_changed), (gpointer) &analysis->atom2);

/* terminating buttons */
gui_stock_button(GTK_STOCK_EXECUTE, exec_analysis_task, dialog,
                   GTK_DIALOG(window)->action_area);

gui_stock_button(GTK_STOCK_CLOSE, dialog_destroy, dialog,
                   GTK_DIALOG(window)->action_area);

gtk_widget_show_all(window);
}

