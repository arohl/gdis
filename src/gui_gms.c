/*
Copyright (C) 2003 by Andrew Lloyd Rohl and Sean David Fleming

andrew@ivec.org
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
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <unistd.h>

#include "gdis.h"
#include "gamess.h"
#include "coords.h"
#include "file.h"
#include "task.h"
#include "grid.h"
#include "matrix.h"
#include "parse.h"
#include "spatial.h"
#include "gui_shorts.h"
#include "dialog.h"
#include "interface.h"
#include "opengl.h"

#define CONTROL_TAB_LABEL " Control "
#define BASIS_TAB_LABEL " Basis Set "
#define OPTIMISATION_TAB_LABEL " Optimisation "

extern struct basis_pak basis_sets[];

enum {
  TIMLIM, MWORDS, WIDE_OUTPUT
};

struct callback_pak run_control[] = {
{"Time limit (minutes) ", TIMLIM},
{"Memory limit (megawords) ", MWORDS},
{NULL, 0}
};

enum {
  TOTAL_Q, MULTIPLICITY
};

struct callback_pak econfig[] = {
{"Total charge ", TOTAL_Q},
{"Multiplicity ", MULTIPLICITY},
{NULL, 0}
};

enum {
  NUM_P, NUM_D, NUM_F
};

struct callback_pak pol_funcs[] = {
{"p ", NUM_P},
{"d ", NUM_D},
{"f ", NUM_F},
{NULL, 0}
};

enum {
  HEAVY_DIFFUSE, HYDROGEN_DIFFUSE, NSTEP, GAMESS_TITLE, GAMESS_ENERGY, INPUTFILE, MODELNAME
};

extern struct sysenv_pak sysenv;

/*********************/
/* run a gamess file */
/*********************/
#define DEBUG_EXEC_GAMESS 0
gint exec_gamess(gchar *input,struct task_pak *task)
{
gchar *cmd, hostname[256];

/* checks */
if (!sysenv.gamess_path)
  return(-1);

/* NEW - acquire hostname for GAMESS job submission */
/*
if (gethostname(hostname, sizeof(hostname)))
*/
  sprintf(hostname, "localhost");

#if __WIN32
{
gchar *basename, *bat, *fubar, *fname, *temp, *tmp, *dir, *inp, *out;
FILE *fp;

/* setup variables */
basename = parse_strip(input);
fname = g_strdup_printf("%s.F05", basename);
fubar = g_build_filename(sysenv.gamess_path, "scratch", fname, NULL);
temp = g_build_filename(sysenv.gamess_path, "temp", basename, NULL);
dir = g_build_filename(sysenv.gamess_path, NULL);
inp = g_build_filename(sysenv.cwd, input, NULL);
out = g_build_filename(sysenv.cwd, basename, NULL);

/* create a batch file to run GAMESS */
bat = gun("bat");
fp = fopen(bat, "wt");
fprintf(fp, "@echo off\n");
fprintf(fp, "echo Running GAMESS using 1 CPU, job: %s\n", input);

/* remove old crap in the temp directory */
fprintf(fp, "del %s.*\n", temp);

/* put the input file in the scratch directory */
fprintf(fp, "copy \"%s\" \"%s\"\n", inp, fubar);

/* run the execution script */
fprintf(fp, "cd %s\n", sysenv.gamess_path);
fprintf(fp, "csh -f runscript.csh %s 04 1 %s %s > \"%s.gmot\"\n",
             basename, dir, hostname, out);
fprintf(fp, "echo GAMESS job completed\n");
fclose(fp);

tmp = g_build_filename(sysenv.cwd, bat, NULL);
cmd = g_strdup_printf("\"%s\"", tmp);

g_free(basename);
g_free(fname);
g_free(fubar);
g_free(temp);
g_free(tmp);
g_free(inp);
g_free(out);
g_free(bat);
}
#else
cmd = g_strdup_printf("%s/%s %s", sysenv.gamess_path, sysenv.gamess_exe, input); 
#endif

#if DEBUG_EXEC_GAMESS
printf("executing: [%s]\n",cmd);
#else

task->is_async = TRUE;
//task->status_file = g_build_filename(sysenv.cwd,???,NULL);/*what should be the status file?*/
task_async(cmd,&(task->pid));
#endif

g_free(cmd);

return(0);
}

/*****************************/
/* execute a gamess run (task) */
/*****************************/
#define DEBUG_EXEC_GAMESS_TASK 0
void exec_gamess_task(struct model_pak *model,struct task_pak *task)
{
gchar *out;

/* construct output filename */
g_free(model->gamess.out_file);

out = parse_extension_set(model->gamess.temp_file, "gmot");
model->gamess.out_file = g_build_filename(sysenv.cwd, out, NULL);
g_free(out);

#if __WIN32
/*
out = g_strdup_printf("\"%s\"", model->gamess.out_file);
*/
out = g_shell_quote(model->gamess.out_file);
g_free(model->gamess.out_file);
model->gamess.out_file = out;
#endif

/* save input file and execute */
file_save_as(model->gamess.temp_file, model);
exec_gamess(model->gamess.temp_file,task);
}

/*****************************/
/* process a gamess run (task) */
/*****************************/
#define DEBUG_PROC_GAMESS 0
void proc_gamess_task(gpointer ptr)
{
gchar *filename;
GString *line;
struct model_pak *data;

/* TODO - model locking (moldel_ptr RO/RW etc) to prevent screw ups */
g_return_if_fail(ptr != NULL);
data = ptr;

/* win32 fix - which can't make up its mind if it wants to quote or not */
filename = g_shell_unquote(data->gamess.out_file, NULL);

if (!g_file_test(filename, G_FILE_TEST_EXISTS))
  {
  printf("Missing output file [%s]\n", filename);
  return;
  }

if (data->gamess.run_type < GMS_OPTIMIZE && !data->animation)
  {
/* same model (ie with current energetics dialog) so update */
    file_load(filename, data);
/* update energy (TODO - only if successful) */
    line = g_string_new(NULL);
    if (data->gamess.have_energy)
      g_string_sprintf(line,"%f",data->gamess.energy);
    else
      g_string_sprintf(line, "not yet calculated");
/* is there a dialog entry to be updated? */
    if (GTK_IS_ENTRY(data->gamess.energy_entry))
      gtk_entry_set_text(GTK_ENTRY(data->gamess.energy_entry), line->str);
    g_string_free(line, TRUE);
  }
else
  {
/* TODO - make it possile to get dialog data by request */
/* so that we can check if a dialog exsits to be updated */
/* get new coords */
/* create new model for the minimized result */
  file_load(filename, NULL);
  }

g_free(filename);

redraw_canvas(ALL);
}

/*************************************************/
/* controls the use of GAMESS to optimise coords */
/*************************************************/
#define DEBUG_RUN_GAMESS 0
void gui_gamess_task(GtkWidget *w, struct model_pak *data)
{
/* checks */
g_return_if_fail(data != NULL);
if (!sysenv.gamess_path)
  {
  gui_text_show(ERROR, "GAMESS executable was not found.\n");
  return;
  }

#if DEBUG_RUN_GAMESS
printf("output file: %s\n", data->gamess.out_file);
#endif

task_new("Gamess", &exec_gamess_task, data, &proc_gamess_task, data, data);
}

/*********************************/
/* GAMESS execution type toggles */
/*********************************/
void gamess_exe_run(struct model_pak *mdata)
{
mdata->gamess.exe_type = GMS_RUN;
}
void gamess_exe_check(struct model_pak *mdata)
{
mdata->gamess.exe_type = GMS_CHECK;
}
void gamess_exe_debug(struct model_pak *mdata)
{
mdata->gamess.exe_type = GMS_DEBUG;
}

/***************************/
/* GAMESS run type toggles */
/***************************/
void gamess_run_single(struct model_pak *data)
{
data->gamess.run_type = GMS_ENERGY;
}
void gamess_run_gradient(struct model_pak *data)
{
data->gamess.run_type = GMS_GRADIENT;
}
void gamess_run_hessian(struct model_pak *data)
{
data->gamess.run_type = GMS_HESSIAN;
}
void gamess_run_optimize(struct model_pak *data)
{
data->gamess.run_type = GMS_OPTIMIZE;
}

/***************************/
/* GAMESS SCF type toggles */
/***************************/
void gamess_scf_rhf(struct model_pak *data)
{
data->gamess.scf_type = GMS_RHF;
}
void gamess_scf_uhf(struct model_pak *data)
{
data->gamess.scf_type = GMS_UHF;
}
void gamess_scf_rohf(struct model_pak *data)
{
data->gamess.scf_type = GMS_ROHF;
}
void gamess_units_angs(struct model_pak *data)
{
data->gamess.units = GMS_ANGS;
}
void gamess_units_bohr(struct model_pak *data)
{
data->gamess.units = GMS_BOHR;
}
void gamess_optimise_qa(struct model_pak *data)
{
data->gamess.opt_type = GMS_QA;
}
void gamess_optimise_nr(struct model_pak *data)
{
data->gamess.opt_type = GMS_NR;
}
void gamess_optimise_rfo(struct model_pak *data)
{
data->gamess.opt_type = GMS_RFO;
}
void gamess_optimise_schlegel(struct model_pak *data)
{
data->gamess.opt_type = GMS_SCHLEGEL;
}
void event_destroy_user_data_from_object(GtkObject *object, gpointer user_data)
{
g_free(user_data);
}

/*****************************/
/* basis set selection event */
/*****************************/
void basis_selected(GtkWidget *w, struct model_pak *model)
{
gint i=0;
const gchar *label;

label = gtk_entry_get_text(GTK_ENTRY(w));

while (basis_sets[i].label)
  {
  if (g_ascii_strcasecmp(label, basis_sets[i].label) == 0)
    {
    model->gamess.basis = basis_sets[i].basis;
    model->gamess.ngauss = basis_sets[i].ngauss;
    }
  i++;
  }
}

/******************************/
/* functional selection event */
/******************************/
void gms_functional_select(GtkWidget *w, struct model_pak *model)
{
const gchar *label;

label = gtk_entry_get_text(GTK_ENTRY(w));

if (g_ascii_strncasecmp(label, "svwn", 4) == 0)
  {
  model->gamess.dft_functional = SVWN;
  return;
  }
if (g_ascii_strncasecmp(label, "blyp", 4) == 0)
  {
  model->gamess.dft_functional = BLYP;
  return;
  }
if (g_ascii_strncasecmp(label, "b3ly", 4) == 0)
  {
  model->gamess.dft_functional = B3LYP;
  return;
  }
}

/***************************/
/* text handling callbacks */
/***************************/
gint event_pol_text_fields(GtkWidget *w, gpointer *obj)
{
gint id;
struct model_pak *data;

/* checks */
g_return_val_if_fail(w != NULL, FALSE);

/* ascertain type of modification required */
id = GPOINTER_TO_INT(g_object_get_data(G_OBJECT(w), "id"));

data = sysenv.active_model;
sysenv.moving = FALSE;

switch(id)
  {
  case GAMESS_TITLE:
    g_free(data->gamess.title);
    data->gamess.title = g_strdup(gtk_entry_get_text(GTK_ENTRY(w)));
    break;
  case INPUTFILE:
    g_free(data->gamess.temp_file);
    data->gamess.temp_file = g_strdup(gtk_entry_get_text(GTK_ENTRY(w)));
    break; 
  }

return(FALSE);
}

/*********************************/
/* GAMESS job setup/results page */
/*********************************/
void gui_gamess_widget(GtkWidget *box, struct model_pak *data)
{
gint i, j;
GString *line;
GList *list;
GtkWidget *hbox, *vbox, *vbox1, *vbox2, *vbox3, *vbox4, *page;
GtkWidget *frame, *button, *label, *entry, *notebook;
GtkWidget *combo;

/* checks */
g_assert(box != NULL);
g_assert(data != NULL);

/* string manipulation scratchpad */
line = g_string_new(NULL);

/* create notebook */
notebook = gtk_notebook_new();
gtk_notebook_set_tab_pos(GTK_NOTEBOOK(notebook), GTK_POS_TOP);
gtk_container_add(GTK_CONTAINER(box), notebook);
gtk_notebook_set_show_border(GTK_NOTEBOOK(notebook), FALSE);

/* Control page */
page = gtk_vbox_new(FALSE, 0);
label = gtk_label_new (CONTROL_TAB_LABEL);
gtk_notebook_append_page(GTK_NOTEBOOK(notebook), page, label);
hbox = gtk_hbox_new(FALSE, 0);
gtk_container_add(GTK_CONTAINER(page), hbox);

/* split panel */
vbox1 = gtk_vbox_new(FALSE, 0);
gtk_box_pack_start(GTK_BOX(hbox), vbox1, FALSE, FALSE, 0);
vbox2 = gtk_vbox_new(FALSE, 0);
gtk_box_pack_start(GTK_BOX(hbox), vbox2, FALSE, FALSE, 0);
vbox3 = gtk_vbox_new(FALSE, 0);
gtk_box_pack_start(GTK_BOX(hbox), vbox3, FALSE, FALSE, 0);

/* frame for execution type */
frame = gtk_frame_new("Execution type");
gtk_box_pack_start(GTK_BOX(vbox1), frame, FALSE, FALSE, 0);
gtk_container_set_border_width(GTK_CONTAINER(frame), PANEL_SPACING);
vbox = gtk_vbox_new(TRUE, 0);
gtk_container_add(GTK_CONTAINER(frame), vbox);

/* radio buttons for execution type */
new_radio_group(0, vbox, TT);
button = add_radio_button("Run", (gpointer) gamess_exe_run, data);
if (data->gamess.exe_type == GMS_RUN)
  gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(button), TRUE);

button = add_radio_button("Check", (gpointer) gamess_exe_check, data);
if (data->gamess.exe_type == GMS_CHECK)
  gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(button), TRUE);

button = add_radio_button("Debug", (gpointer) gamess_exe_debug, data);
if (data->gamess.exe_type == GMS_DEBUG)
  gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(button), TRUE);

/* frame for run type */
frame = gtk_frame_new("Run type");
gtk_box_pack_start(GTK_BOX(vbox1), frame, FALSE, FALSE, 0);
gtk_container_set_border_width(GTK_CONTAINER(frame), PANEL_SPACING);
vbox = gtk_vbox_new(TRUE, 0);
gtk_container_add(GTK_CONTAINER(frame), vbox);

/* radio buttons for run type */
new_radio_group(0, vbox, TT);
button = add_radio_button("Single point", (gpointer) gamess_run_single, data);
if (data->gamess.run_type == GMS_ENERGY)
  gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(button), TRUE);

button = add_radio_button("Gradient", (gpointer) gamess_run_gradient, data);
if (data->gamess.run_type == GMS_GRADIENT)
  gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(button), TRUE);

button = add_radio_button("Hessian", (gpointer) gamess_run_hessian, data);
if (data->gamess.run_type == GMS_HESSIAN)
  gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(button), TRUE);

button = add_radio_button("Optimize", (gpointer) gamess_run_optimize, data);
if (data->gamess.run_type == GMS_OPTIMIZE)
  gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(button), TRUE);

/* frame for SCF type */
frame = gtk_frame_new("SCF type");
gtk_box_pack_start(GTK_BOX(vbox2), frame, FALSE, FALSE, 0);
gtk_container_set_border_width(GTK_CONTAINER(frame), PANEL_SPACING);
vbox = gtk_vbox_new(TRUE, 0);
gtk_container_add(GTK_CONTAINER(frame), vbox);

/* radio buttons for SCF type */
new_radio_group(0, vbox, TT);
button = add_radio_button("RHF", (gpointer) gamess_scf_rhf, data);
if (data->gamess.scf_type == GMS_RHF)
  gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(button), TRUE);

button = add_radio_button("UHF", (gpointer) gamess_scf_uhf, data);
if (data->gamess.scf_type == GMS_UHF)
  gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(button), TRUE);

button = add_radio_button("ROHF", (gpointer) gamess_scf_rohf, data);
if (data->gamess.scf_type == GMS_ROHF)
  gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(button), TRUE);

/* frame for input units */
frame = gtk_frame_new("Input units");
gtk_box_pack_start(GTK_BOX(vbox2), frame, FALSE, FALSE, 0);
gtk_container_set_border_width(GTK_CONTAINER(frame), PANEL_SPACING);
vbox = gtk_vbox_new(TRUE, 0);
gtk_container_add(GTK_CONTAINER(frame), vbox);

/* radio buttons for SCF type */
new_radio_group(0, vbox, TT);
button = add_radio_button("Angstrom", (gpointer) gamess_units_angs, data);
if (data->gamess.units == GMS_ANGS)
  gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(button), TRUE);

button = add_radio_button("Bohr", (gpointer) gamess_units_bohr, data);
if (data->gamess.units == GMS_BOHR)
  gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(button), TRUE);

/* frame for SCF options */
frame = gtk_frame_new("SCF options");
gtk_box_pack_start(GTK_BOX(vbox2), frame, FALSE, FALSE, 0);
gtk_container_set_border_width(GTK_CONTAINER(frame), PANEL_SPACING);
vbox = gtk_vbox_new(TRUE, 0);
gtk_container_add(GTK_CONTAINER(frame), vbox);
gui_direct_spin("Maximum iterations", &data->gamess.maxit, 1, 150, 1,
                  NULL, NULL, vbox);

/* frame for DFT */
frame = gtk_frame_new(NULL);
gtk_box_pack_start(GTK_BOX(vbox3), frame, FALSE, FALSE, 0);
gtk_container_set_border_width(GTK_CONTAINER(frame), PANEL_SPACING);
vbox = gtk_vbox_new(TRUE, 0);
gtk_container_add(GTK_CONTAINER(frame), vbox);

vbox4 = gtk_vbox_new(TRUE, 0);
/* TODO - set sensitive stuff -> candidate for a short cut */
gui_direct_check("DFT calculation", &data->gamess.dft, gui_checkbox_refresh, vbox4, vbox);
gtk_box_pack_start(GTK_BOX(vbox), vbox4, FALSE, FALSE, 0);


/* TODO - the pulldown list is another candidate for a short cut */
hbox = gtk_hbox_new(FALSE, PANEL_SPACING);
gtk_box_pack_start(GTK_BOX(vbox4), hbox, FALSE, FALSE, 0);

label = gtk_label_new("Functional");
gtk_box_pack_start(GTK_BOX(hbox), label, FALSE, FALSE, 0);

list = NULL;
list = g_list_append(list, "SVWN (LDA)");
list = g_list_append(list, "BLYP");
list = g_list_append(list, "B3LYP");

combo = gtk_combo_new();
gtk_entry_set_editable(GTK_ENTRY(GTK_COMBO(combo)->entry), FALSE);
gtk_combo_set_popdown_strings(GTK_COMBO(combo), list);
gtk_box_pack_start(GTK_BOX(hbox), combo, FALSE, FALSE, 0);

g_signal_connect(GTK_OBJECT(GTK_COMBO(combo)->entry), "changed", 
                 GTK_SIGNAL_FUNC(gms_functional_select), data);




/* frame for time, memory and output */
frame = gtk_frame_new("Run control");
gtk_box_pack_start(GTK_BOX(vbox3), frame, FALSE, FALSE, 0);
gtk_container_set_border_width(GTK_CONTAINER(frame), PANEL_SPACING);
vbox = gtk_vbox_new(TRUE, 0);
gtk_container_add(GTK_CONTAINER(frame), vbox);

/* Add spinners */
i=0;
while (run_control[i].label)
  {
  switch (run_control[i].id)
    {
    case TIMLIM:
      gui_direct_spin(run_control[i].label,
                        &data->gamess.time_limit, 0, 30000, 10,
                        NULL, NULL, vbox);
      break;

    case MWORDS:
      gui_direct_spin(run_control[i].label,
                        &data->gamess.mwords, 1, 150, 1,
                        NULL, NULL, vbox);
      break;

    default:
      g_assert_not_reached();
    }
  i++;
  }

/* wide output button */
gui_direct_check("Wide output", &data->gamess.wide_output, NULL, NULL, vbox);

/* frame for charge and multiplicity */
frame = gtk_frame_new("Electronic configuration");
gtk_box_pack_start(GTK_BOX(vbox3), frame, FALSE, FALSE, 0);
gtk_container_set_border_width(GTK_CONTAINER(frame), PANEL_SPACING);
vbox = gtk_vbox_new(TRUE, 0);
gtk_container_add(GTK_CONTAINER(frame), vbox);

/* Add spinners */
i=0;
while (econfig[i].label)
  {
  switch (econfig[i].id)
    {
    case TOTAL_Q:
      gui_direct_spin(econfig[i].label,
                        &data->gamess.total_charge, -5, 5, 1,
                        NULL, NULL, vbox);
      break;

    case MULTIPLICITY:
      gui_direct_spin(econfig[i].label,
                        &data->gamess.multiplicity, 1, 5, 1,
                        NULL, NULL, vbox);
      break;

    default:
      g_assert_not_reached();
    }
  i++;
  }
  
/* basis set page */
page = gtk_vbox_new(FALSE,0);
label = gtk_label_new (BASIS_TAB_LABEL);
gtk_notebook_append_page(GTK_NOTEBOOK(notebook), page, label);
hbox = gtk_hbox_new(FALSE, 0);
gtk_box_pack_start(GTK_BOX(page), hbox, FALSE, FALSE, 0);

/* split panel */
vbox1 = gtk_vbox_new(FALSE, 0);
gtk_box_pack_start(GTK_BOX(hbox), vbox1, TRUE, TRUE, 0);
vbox2 = gtk_vbox_new(FALSE, 0);
gtk_box_pack_start(GTK_BOX(hbox), vbox2, FALSE, FALSE, 0);

/* frame for basis set */
frame = gtk_frame_new("Basis set");
gtk_box_pack_start(GTK_BOX(vbox1), frame, FALSE, FALSE, 0);
gtk_container_set_border_width(GTK_CONTAINER(frame), PANEL_SPACING);
vbox = gtk_vbox_new(FALSE, 0);
gtk_container_add(GTK_CONTAINER(frame), vbox);
j = -1;
i = 0;
list = NULL;
while (basis_sets[i].label)
  {
  list = g_list_append(list, basis_sets[i].label);
  if ((basis_sets[i].basis == data->gamess.basis) && (basis_sets[i].ngauss == data->gamess.ngauss))
    j = i;
  i++;
  }

combo = gtk_combo_new();
gtk_entry_set_editable(GTK_ENTRY(GTK_COMBO(combo)->entry), FALSE);
gtk_combo_set_popdown_strings(GTK_COMBO(combo), list);
gtk_box_pack_start(GTK_BOX(vbox), combo, FALSE, FALSE, 0);
/* NB: set text BEFORE the changed event is connected */
if (j >= 0)
  gtk_entry_set_text(GTK_ENTRY(GTK_COMBO(combo)->entry), basis_sets[j].label);
g_signal_connect(GTK_OBJECT(GTK_COMBO(combo)->entry), "changed", 
                 GTK_SIGNAL_FUNC(basis_selected), data);

/* frame for polarization functions set */
frame = gtk_frame_new("Polarization functions");
gtk_box_pack_start(GTK_BOX(vbox2), frame, FALSE, FALSE, 0);
gtk_container_set_border_width(GTK_CONTAINER(frame), PANEL_SPACING);
vbox = gtk_vbox_new(TRUE, 0);
gtk_container_add(GTK_CONTAINER(frame), vbox);

/* Add spinners */
i=0;
while (pol_funcs[i].label)
  {
  switch (pol_funcs[i].id)
    {
    case NUM_P:
      gui_direct_spin(pol_funcs[i].label, &data->gamess.num_p, 0, 3, 1,
                        NULL, NULL, vbox);
      break;
    case NUM_D:
      gui_direct_spin(pol_funcs[i].label, &data->gamess.num_d, 0, 3, 1,
                        NULL, NULL, vbox);
      break;
    case NUM_F:
      gui_direct_spin(pol_funcs[i].label, &data->gamess.num_f, 0, 1, 1,
                        NULL, NULL, vbox);
      break;
    default:
      g_assert_not_reached();
    }
  i++;
  }

/* frame for polarization functions set */
frame = gtk_frame_new("Diffuse functions");
gtk_box_pack_start(GTK_BOX(vbox2), frame, FALSE, FALSE, 0);
gtk_container_set_border_width(GTK_CONTAINER(frame), PANEL_SPACING);
vbox = gtk_vbox_new(TRUE, 0);
gtk_container_add(GTK_CONTAINER(frame), vbox);
gui_direct_check("Heavy atoms (s & p)", &data->gamess.have_heavy_diffuse,
                    NULL, NULL, vbox);
gui_direct_check("Hydrogen (s only)", &data->gamess.have_hydrogen_diffuse,
                    NULL, NULL, vbox);

/* frame for optimisation settings */
page = gtk_vbox_new(FALSE,0);
label = gtk_label_new (OPTIMISATION_TAB_LABEL);
gtk_notebook_append_page(GTK_NOTEBOOK(notebook), page, label);
hbox = gtk_hbox_new(FALSE, 0);
gtk_box_pack_start(GTK_BOX(page), hbox, FALSE, FALSE, 0);

/* split panel */
vbox1 = gtk_vbox_new(FALSE, 0);
gtk_box_pack_start(GTK_BOX(hbox), vbox1, FALSE, FALSE, 0);
vbox2 = gtk_vbox_new(FALSE, 0);
gtk_box_pack_start(GTK_BOX(hbox), vbox2, FALSE, FALSE, 0);

/* optimiser to use */
frame = gtk_frame_new("Optimiser");
gtk_box_pack_start(GTK_BOX(vbox1), frame, FALSE, FALSE, 0);
gtk_container_set_border_width(GTK_CONTAINER(frame), PANEL_SPACING);
vbox = gtk_vbox_new(TRUE, 0);
gtk_container_add(GTK_CONTAINER(frame), vbox);

/* radio buttons for run type */
new_radio_group(0, vbox, TT);
button = add_radio_button("Quadratic approximation", (gpointer) gamess_optimise_qa, data);
if (data->gamess.opt_type == GMS_QA)
  gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(button), TRUE);
button = add_radio_button("Newton-Raphson", (gpointer) gamess_optimise_nr, data);
if (data->gamess.opt_type == GMS_NR)
  gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(button), TRUE);
button = add_radio_button("RFO", (gpointer) gamess_optimise_rfo, data);
if (data->gamess.opt_type == GMS_RFO)
  gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(button), TRUE);
button = add_radio_button("Quasi-NR", (gpointer) gamess_optimise_schlegel, data);
if (data->gamess.opt_type == GMS_SCHLEGEL)
  gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(button), TRUE);

/* frame for optimization cycles */
frame = gtk_frame_new(NULL);
gtk_box_pack_start(GTK_BOX(vbox1), frame, FALSE, FALSE, 0);
gtk_container_set_border_width(GTK_CONTAINER(frame), PANEL_SPACING);
vbox = gtk_vbox_new(TRUE, 0);
gtk_container_add(GTK_CONTAINER(frame), vbox);

/* optimization cycles */
hbox = gtk_hbox_new(FALSE, 0);
gtk_box_pack_start(GTK_BOX(vbox), hbox, FALSE, FALSE, 0);
gui_direct_spin("Maximum cycles", &data->gamess.nstep,
                  0, 500, 10, NULL, NULL, hbox);


/* setup some gulp file name defaults */
if (g_ascii_strncasecmp(data->gamess.temp_file, "none", 4) == 0)
  {
  g_free(data->gamess.temp_file);
  data->gamess.temp_file = g_strdup_printf("%s.inp", data->basename); 
  }
if (g_ascii_strncasecmp(data->gamess.out_file, "none", 4) == 0)
  {
  g_free(data->gamess.out_file);
  data->gamess.out_file = g_strdup_printf("%s.gmot", data->basename); 
  }


frame = gtk_frame_new("Files");
gtk_box_pack_start(GTK_BOX(box),frame,FALSE,FALSE,0);
gtk_container_set_border_width(GTK_CONTAINER(frame), PANEL_SPACING);
hbox = gtk_hbox_new(FALSE, 0);
gtk_container_add (GTK_CONTAINER(frame),hbox);
gtk_container_set_border_width(GTK_CONTAINER(hbox), PANEL_SPACING);

/* left vbox */
vbox = gtk_vbox_new(TRUE, 2);
gtk_box_pack_start (GTK_BOX (hbox), vbox, TRUE, TRUE, 0);

/* GAMESS input filename */
label = gtk_label_new(" Job input file ");
gtk_box_pack_start (GTK_BOX (vbox), label, TRUE, TRUE, 0);

/* right vbox */
vbox = gtk_vbox_new(TRUE, 2);
gtk_box_pack_start (GTK_BOX (hbox), vbox, TRUE, TRUE, 0);

/* temporary input */
entry = gtk_entry_new_with_max_length(LINELEN);
gtk_entry_set_text(GTK_ENTRY(entry), data->gamess.temp_file);
gtk_box_pack_start(GTK_BOX(vbox), entry, TRUE, TRUE, 0);
g_signal_connect(GTK_OBJECT (entry), "changed",
                          GTK_SIGNAL_FUNC (event_pol_text_fields),
                         NULL);
g_object_set_data(G_OBJECT(entry), "id", (gpointer) INPUTFILE);

/* details */
frame = gtk_frame_new("Details");
gtk_box_pack_start(GTK_BOX(box),frame,FALSE,FALSE,0);
gtk_container_set_border_width(GTK_CONTAINER(frame), PANEL_SPACING);
hbox = gtk_hbox_new(FALSE, 0);
gtk_container_add (GTK_CONTAINER(frame),hbox);

/* left vbox */
vbox = gtk_vbox_new(TRUE, 2);
gtk_box_pack_start (GTK_BOX (hbox), vbox, TRUE, TRUE, 0);

if (g_ascii_strncasecmp(data->gamess.title, "none", 4) == 0)
  {
  g_free(data->gamess.title);
  data->gamess.title = g_strdup(data->basename); 
  }

/* temporary GAMESS input filename */
label = gtk_label_new("Title");
gtk_box_pack_start (GTK_BOX (vbox), label, TRUE, TRUE, 0);
label = gtk_label_new("Total energy (Hartrees)");
gtk_box_pack_start (GTK_BOX (vbox), label, TRUE, TRUE, 0);

/* right vbox */
vbox = gtk_vbox_new(TRUE, 2);
gtk_box_pack_start (GTK_BOX (hbox), vbox, TRUE, TRUE, 0);

entry = gtk_entry_new_with_max_length(LINELEN);
gtk_entry_set_text(GTK_ENTRY(entry), data->gamess.title);
gtk_box_pack_start (GTK_BOX (vbox), entry, TRUE, TRUE, 0);
g_signal_connect(GTK_OBJECT (entry), "changed",
                          GTK_SIGNAL_FUNC (event_pol_text_fields),
                         NULL);
g_object_set_data(G_OBJECT(entry), "id", (gpointer) GAMESS_TITLE);

data->gamess.energy_entry = gtk_entry_new_with_max_length(LINELEN);
if (data->gamess.have_energy)
  g_string_sprintf(line,"%f",data->gamess.energy);
else
  g_string_sprintf(line,"not yet calculated");

gtk_entry_set_text(GTK_ENTRY(data->gamess.energy_entry),line->str);
gtk_box_pack_start (GTK_BOX (vbox), data->gamess.energy_entry, TRUE, TRUE, 0);
gtk_entry_set_editable(GTK_ENTRY(data->gamess.energy_entry), FALSE);

/* done */
gtk_widget_show_all(box);

if (data->gamess.dft)
  gtk_widget_set_sensitive(vbox4, TRUE);
else
  gtk_widget_set_sensitive(vbox4, FALSE);

g_string_free(line, TRUE);
}

/***********************/
/* grid job submission */
/***********************/
void gui_gms_grid_submit(GtkWidget *w, gpointer data)
{
#ifdef WITH_GRISU
gchar *tmp;
struct model_pak *model = data;
struct grid_pak *grid;

g_assert(model != NULL);

grid = grid_new();
grid->local_cwd = g_strdup(sysenv.cwd);
grid->local_input = g_strdup(model->gamess.temp_file);

/* immediately write to disk so job is preserved even if user deletes model */
tmp = g_build_filename(grid->local_cwd, grid->local_input, NULL);
write_gms(tmp, model);
g_free(tmp);

/* TODO - deal with gamess vs gamess-us problems that could arise */
if (!grid_job_submit("gamess-us", grid))
  {
/* TODO - failed message */
  grid_free(grid);
  }
#endif
}

/****************************************/
/* run the approriate energetics method */
/****************************************/
void gamess_execute(GtkWidget *w, gpointer data)
{
struct model_pak *model = data;

g_assert(model != NULL);

gui_gamess_task(NULL, model);
}

/****************************/
/* GAMESS energetics dialog */
/****************************/
void gamess_dialog(void)
{
gpointer dialog;
GtkWidget *window, *frame, *vbox;
struct model_pak *model;

model = sysenv.active_model;
if (!model)
  return;

/* request an energetics dialog */
dialog = dialog_request(GAMESS, "GAMESS configuration", NULL, NULL, NULL);
if (!dialog)
  return;
window = dialog_window(dialog);

frame = gtk_frame_new(NULL);
gtk_box_pack_start(GTK_BOX(GTK_DIALOG(window)->vbox), frame, FALSE, FALSE, 0);
gtk_container_set_border_width(GTK_CONTAINER(frame), PANEL_SPACING);

vbox = gtk_vbox_new(FALSE,0);
gtk_container_add(GTK_CONTAINER(frame), vbox);
gui_gamess_widget(vbox, model);

/* terminating button */
gui_stock_button(GTK_STOCK_EXECUTE, gamess_execute, model,
                   GTK_DIALOG(window)->action_area);

#ifdef WITH_GRISU
gui_icon_button(GTK_STOCK_APPLY, "Grid Submit", gui_gms_grid_submit, model, GTK_DIALOG(window)->action_area);
#endif

gui_stock_button(GTK_STOCK_CLOSE, dialog_destroy, dialog,
                   GTK_DIALOG(window)->action_area);

gtk_widget_show_all(window);
}

