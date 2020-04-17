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
#include <unistd.h>
#include <math.h>
#include <sys/stat.h>
#include <sys/types.h>

#include "gdis.h"
#include "coords.h"
#include "model.h"
#include "file.h"
#include "matrix.h"
#include "module.h"
#include "parse.h"
#include "project.h"
#include "quaternion.h"
#include "task.h"
#include "interface.h"
#include "dialog.h"
#include "gui_shorts.h"

/* globals */
extern struct sysenv_pak sysenv;
gint dock_no_execute=TRUE, dock_grid_on=TRUE, dock_rotate_on=TRUE;
gint dock_rigid_on=TRUE, dock_rigid_x=FALSE, dock_rigid_y=FALSE, dock_rigid_z=TRUE;
gdouble dock_grid[3] = {2,2,2}, dock_rotate[3] = {1, 1, 4};
gdouble dock_cell[2] = {1.0, 1.0};

struct dock_pak
{
gchar *path;
};

/*****************************/
/* module configure function */
/*****************************/
/*
gint module_libdock_config(gpointer module)
{
module_label_set("docking", module);
module_symbol_add("docking_setup", module);
module_resident_set(TRUE, module);
return(0);
}
*/

/*************************/
/* docking task analysis */
/*************************/
void docking_cleanup(struct dock_pak *dock, struct task_pak *task)
{
gchar *path, *ext;
GSList *list, *dir_list;
struct file_pak *file;
struct model_pak model;

/* temporary model data storage */
model_init(&model);

/* load gulp file energies into project */
/* TODO - scan for .res files */
/* TODO - if none - scan .got file for possible error reports */
path = dock->path;

/* scan project directory for .res files */
dir_list = file_dir_list(path, 0);
for (list=dir_list ; list ; list=g_slist_next(list))
  {
  file = get_file_info(list->data, BY_EXTENSION);
  if (file)
    {
    if (file->id == GULP)
      {
      ext  = find_char(list->data, '.', LAST);
      if (g_ascii_strncasecmp(ext, ".res", 3) == 0)
        {
/* process dump file */
        printf("%s\n", (gchar *) list->data);
        }
      }
    if (file->id == GULPOUT)
      {
/* process output file */

read_gulp_output(list->data, &model);

      printf("%s : %f\n", (gchar *) list->data, model.gulp.energy);
      }
    }
  }
/* cleanup */
g_slist_free(dir_list);
model_free(&model);
}

/***************************/
/* docking task execution */
/**************************/
void docking_execute(struct dock_pak *dock, struct task_pak *task)
{
gchar *path, *cmd, *base, *input, *output;
GSList *list, *dir_list;
struct file_pak *file;

path = dock->path;

/*
printf("expected files: %d\n", dock->n);
*/

/* NB: have to change to the project directory in order to run GULP jobs in it */
/* TODO - save sysenv.cwd and restore after? */
if (chdir(path))
  {
  printf("Failed to change to project directory.\n");
  return;
  }

/* submit all GULP jobs in the project directory */
dir_list = file_dir_list(path, 0);
task->progress = 0.0;
for (list=dir_list ; list ; list=g_slist_next(list))
  {
  file = get_file_info(list->data, BY_EXTENSION);
  if (file)
    {
    if (file->id == GULP)
      {
      input = list->data;
      base = parse_strip(input);
      output = g_strconcat(base, ".got", NULL);

/* this doesn't work - why? */
/*
if (exec_gulp(input, output))
  printf(" >>> exec failed!!!\n");
*/

/* run the GULP jobs (NOT in bg as we're already a task) */
      cmd = g_strdup_printf("%s < %s > %s", sysenv.gulp_path, input, output);
      IGNORE_RETURN(system(cmd));

      g_free(cmd);
      g_free(base);
      g_free(output);
      }
    }
  }
g_slist_free(dir_list);
}

/**************************/
/* create docking project */
/**************************/
void docking_project_create(GtkWidget *w, struct model_pak *model)
{
gint a, b, i, m, n, rx, ry, rz, size, rigid_save;
gint a_max, b_max, rx_max, ry_max, rz_max;
gchar *file, *dump, *dump_save, *rigid_move_save;
gdouble dx, dy, dz, x[3], scale[3], mat[9], dock_centroid[3], q[4];
GString *name, *rigid;
GSList *list, *core_list, *shell_list;
struct dock_pak *dock;
struct core_pak *core, *core2;
struct shel_pak *shell, *shell2;
FILE *fp;

/* checks */
g_assert(model != NULL);
size = g_slist_length(model->selection);
if (!size)
  {
  gui_text_show(WARNING, "Please select the subset you wish to dock.\n");
  return;
  }

/* create new docking project */
dock = g_malloc(sizeof(struct dock_pak));

/* NEW - setup project path */
/*
g_path_get_dirname(model->fullpath);
g_get_current_dir();
*/

/* seek a file name that doesn't exist (avoid background overwriting) */
name = g_string_new(NULL);
i=0;
do
  {
  g_string_sprintf(name, "project_%06d", i);
  i++;
  }
while (g_file_test(name->str, G_FILE_TEST_EXISTS));

dock->path = g_build_path(sysenv.cwd, name->str, NULL);

printf("creating new project: [%s]\n", dock->path);

#if WIN32
if (mkdir(dock->path))
#else
if (mkdir(dock->path, 0700))
#endif
  {
  gui_text_show(ERROR, "Failed to create project directory.\n");
  g_free(dock->path);
  g_free(dock);
  return;
  }

/* project control file */
g_string_sprintf(name, "%s%sproject.pcf", dock->path, DIR_SEP);
fp = fopen(name->str, "wt");

/* save original variables */
dump_save = model->gulp.dump_file;
model->gulp.dump_file = NULL;
rigid_save = model->gulp.rigid;
model->gulp.rigid = dock_rigid_on;
rigid_move_save = model->gulp.rigid_move;
model->gulp.rigid_move = NULL;
if (model->gulp.rigid)
  {
  rigid = g_string_new(NULL);
  if (dock_rigid_x)
    g_string_sprintf(rigid, "x");
  if (dock_rigid_y)
    g_string_sprintfa(rigid, "y");
  if (dock_rigid_z)
    g_string_sprintfa(rigid, "z");
  model->gulp.rigid_move = g_string_free(rigid, FALSE);
  }

/* duplicate selection for docking */
core_list = NULL;
shell_list = NULL;
VEC3SET(dock_centroid, 0.0, 0.0, 0.0);
for (list=model->selection ; list ; list=g_slist_next(list))
  {
  core2 = dup_core(list->data);
  core_list = g_slist_prepend(core_list, core2);
  if (core2->shell)
    shell_list = g_slist_prepend(shell_list, core2->shell);
/* compute centroid */
  ARR3ADD(dock_centroid, core2->x);
  }

/* NB: lists must have the same order as original selection */
core_list = g_slist_reverse(core_list);
shell_list = g_slist_reverse(shell_list);
VEC3MUL(dock_centroid, 1.0/(gdouble) size);

/* fractional translation grid units */
scale[0] = dock_cell[0] / dock_grid[0];
scale[1] = dock_cell[1] / dock_grid[1];

/* rotational increments */
dx = PI/dock_rotate[0];
dy = PI/dock_rotate[1];
dz = PI/dock_rotate[2];

/* translational sampling */
if (dock_grid_on)
  {
  a_max = dock_grid[0];
  b_max = dock_grid[1];
  }
else
  {
  a_max = 1;
  b_max = 1;
  }

/* rotational sampling */
if (dock_rotate_on)
  {
  rx_max = dock_rotate[0];
  ry_max = dock_rotate[1];
  rz_max = dock_rotate[2];
  }
else
  {
  rx_max = 1;
  ry_max = 1;
  rz_max = 1;
  }

/* project header */
fprintf(fp, "%%title solvent mapping project\n");
fprintf(fp, "%%set %d %d %f %f\n", a_max, b_max, dock_cell[0], dock_cell[1]);

/* loop over all grid translations */
m = n = 0;
for (a=0 ; a<a_max ; a++)
  {
  for (b=0 ; b<b_max ; b++)
    {
    VEC3SET(x, a, b, 0.0);
    x[0] *= scale[0];
    x[1] *= scale[1];

/* loop over rotations */
    VEC4SET(q, 1.0, 0.0, 0.0, 0.0);
    for (rx=0 ; rx<rx_max ; rx++)
      {
      if (rx)
        quat_concat_euler(q, PITCH, dx);

    for (ry=0 ; ry<ry_max ; ry++)
      {
      if (ry)
        quat_concat_euler(q, ROLL, dy);

    for (rz=0 ; rz<rz_max ; rz++)
      {
      if (rz)
        quat_concat_euler(q, YAW, dz);

/* build total rotation matrix */
      quat_matrix(mat, q);

/* transform the cores and shells */
      i = 0;
      for (list=model->selection ; list ; list=g_slist_next(list))
        {
        core = list->data;

/* FIXME - should we restore this after? how? */
core->region = 2;

/* get original selection core coordinates */
        core2 = g_slist_nth_data(core_list, i);
        ARR3SET(core->x, core2->x);
/* perform the rotation (NB: must be done about origin in cartesian space) */
        ARR3SUB(core->x, dock_centroid);
        vecmat(model->latmat, core->x);
        vecmat(mat, core->x);
        vecmat(model->ilatmat, core->x);
        ARR3ADD(core->x, dock_centroid);
/* add the current translation offset */
        ARR3ADD(core->x, x);

/* as above, for the associated shell */
        if (core->shell)
          {
          shell = core->shell;
shell->region = 2;
          shell2 = core2->shell;
g_assert(shell2 != NULL);
          ARR3SET(shell->x, shell2->x);
          ARR3SUB(shell->x, dock_centroid);
          vecmat(model->latmat, shell->x);
          vecmat(mat, shell->x);
          vecmat(model->ilatmat, shell->x);
          ARR3ADD(shell->x, dock_centroid);
          ARR3ADD(shell->x, x);
          }
        i++;
        }
/* write docking configuration */
/*
      file = g_strdup_printf("%s_%06d.gin", model->basename, n);
*/

/* m identifies grid points (for later minimum check) */
fprintf(fp, "%s_%06d.gin %f %f %d\n", model->basename, n, x[0], x[1], m);

      file = g_strdup_printf("%s%s%s_%06d.gin", dock->path, DIR_SEP, model->basename, n);
      dump = g_strdup_printf("%s_%06d.res", model->basename, n);

      model->gulp.dump_file = dump;

      write_gulp(file, model);

      g_free(file);
      g_free(dump);
      n++;

      }
      }
      }
    m++;
    }
  }

/* restore original variables */
model->gulp.dump_file = dump_save;
model->gulp.rigid = rigid_save;
g_free(model->gulp.rigid_move);
model->gulp.rigid_move = rigid_move_save;

/* restore original selection (delete, then concat saved list) */
i = 0;
for (list=model->selection ; list ; list=g_slist_next(list))
  {
  core = list->data;
  core2 = g_slist_nth_data(core_list, i);
  ARR3SET(core->x, core2->x);
  if (core->shell)
    {
    shell = core->shell;
    shell2 = core2->shell;
g_assert(shell2 != NULL);
    ARR3SET(shell->x, shell2->x);
    }
  i++;
  }

/* free docking core/shell lists */
free_slist(core_list);
free_slist(shell_list);
g_free(dock);/*FIX f04cde*/
g_string_free(name, TRUE);
fclose(fp);
/* run docking in background unless told to stop after setup */
/*
if (!dock_no_execute)
  submit_task("Docking", &docking_execute, dock, &docking_cleanup, dock, model);
*/
}

/*****************************/
/* docking widget refreshing */
/*****************************/
void dock_toggle_refresh(GtkWidget *w, GtkWidget *box)
{
if (gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(w)))
  gtk_widget_set_sensitive(box, TRUE);
else
  gtk_widget_set_sensitive(box, FALSE);
}

/************************/
/* docking setup dialog */
/************************/
void gui_dock_dialog(void)
{
gchar *text;
gpointer dialog;
GtkWidget *window, *label, *spin;
GtkWidget *frame, *vbox, *vbox1, *hbox, *hbox2, *vbox_left, *vbox_right;
struct model_pak *model;

/* checks */
model = sysenv.active_model;
if (!model)
  {
  gui_text_show(ERROR, "Please load a surface first.\n");
  return;
  }
if (model->periodic != 2)
  {
  gui_text_show(ERROR, "Your model is not a surface.\n");
  return;
  }

/* request a dialog */
dialog = dialog_request(DOCKING, "Docking setup", NULL, NULL, model);
if (!dialog)
  return;
window = dialog_window(dialog);

/* title display */
frame = gtk_frame_new(NULL);
gtk_box_pack_start(GTK_BOX(GTK_DIALOG(window)->vbox), frame, FALSE, FALSE, 0);
vbox = gtk_vbox_new(TRUE, 0);
gtk_container_add(GTK_CONTAINER(frame), vbox);

gtk_container_set_border_width(GTK_CONTAINER(vbox), PANEL_SPACING);

hbox = gtk_hbox_new(FALSE, PANEL_SPACING);
gtk_box_pack_start(GTK_BOX(vbox), hbox, FALSE, FALSE, 0);
text = g_strdup_printf("Docking setup: %s", model->basename);
label = gtk_label_new(text);
gtk_box_pack_start(GTK_BOX(hbox), label, FALSE, FALSE, 0);
g_free(text);

/* NEW - split pane display */
hbox2 = gtk_hbox_new(FALSE, 0);
gtk_box_pack_start(GTK_BOX(GTK_DIALOG(window)->vbox), hbox2, TRUE, TRUE, 0);
vbox_left = gtk_vbox_new(FALSE, 0);
gtk_box_pack_start(GTK_BOX(hbox2), vbox_left, FALSE, FALSE, PANEL_SPACING);
vbox_right = gtk_vbox_new(FALSE, 0);
gtk_box_pack_start(GTK_BOX(hbox2), vbox_right, FALSE, FALSE, PANEL_SPACING);

/* translational sampling */
frame = gtk_frame_new(NULL);
gtk_box_pack_start(GTK_BOX(vbox_left), frame, FALSE, FALSE, 0);
vbox = gtk_vbox_new(FALSE, 0);
gtk_container_add(GTK_CONTAINER(frame), vbox);

vbox1 = gtk_vbox_new(FALSE, 0);
gui_direct_check("Translational sampling", &dock_grid_on,
                   dock_toggle_refresh, vbox1, vbox);
gtk_box_pack_start(GTK_BOX(vbox), vbox1, FALSE, FALSE, PANEL_SPACING);

/* axes sampling fraction */
hbox = gtk_hbox_new(FALSE, 0);
gtk_box_pack_start(GTK_BOX(vbox1), hbox, FALSE, FALSE, PANEL_SPACING);

label = gtk_label_new("axes fractions");
gtk_box_pack_start(GTK_BOX(hbox), label, FALSE, FALSE, PANEL_SPACING);

spin = gui_direct_spin(NULL, &dock_cell[0], 0.0, 1.0, 0.05, NULL, NULL, NULL);
gtk_box_pack_end(GTK_BOX(hbox), spin, FALSE, FALSE, PANEL_SPACING);

spin = gui_direct_spin(NULL, &dock_cell[1], 0.0, 1.0, 0.05, NULL, NULL, NULL);
gtk_box_pack_end(GTK_BOX(hbox), spin, FALSE, FALSE, PANEL_SPACING);

/* grid points */
hbox = gtk_hbox_new(FALSE, 0);
gtk_box_pack_start(GTK_BOX(vbox1), hbox, FALSE, FALSE, PANEL_SPACING);

label = gtk_label_new("grid size");
gtk_box_pack_start(GTK_BOX(hbox), label, FALSE, FALSE, PANEL_SPACING);

spin = gui_direct_spin(NULL, &dock_grid[0], 1, 10, 1, NULL, NULL, NULL);
gtk_box_pack_end(GTK_BOX(hbox), spin, FALSE, FALSE, PANEL_SPACING);

spin = gui_direct_spin(NULL, &dock_grid[1], 1, 10, 1, NULL, NULL, NULL);
gtk_box_pack_end(GTK_BOX(hbox), spin, FALSE, FALSE, PANEL_SPACING);

/* rotational sampling */
frame = gtk_frame_new(NULL);
gtk_box_pack_start(GTK_BOX(vbox_left), frame, TRUE, TRUE, 0);
vbox = gtk_vbox_new(FALSE, 0);
gtk_container_add(GTK_CONTAINER(frame), vbox);

vbox1 = gtk_vbox_new(TRUE, PANEL_SPACING);
gui_direct_check("Rotational sampling", &dock_rotate_on, dock_toggle_refresh, vbox1, vbox);
gtk_box_pack_start(GTK_BOX(vbox), vbox1, FALSE, FALSE, PANEL_SPACING);
gui_direct_spin("  x axis ", &dock_rotate[0], 1, 60, 1, NULL, NULL, vbox1);
gui_direct_spin("  y axis ", &dock_rotate[1], 1, 60, 1, NULL, NULL, vbox1);
gui_direct_spin("  z axis ", &dock_rotate[2], 1, 60, 1, NULL, NULL, vbox1);

/* rigid body docking */
frame = gtk_frame_new(NULL);
gtk_box_pack_start(GTK_BOX(vbox_right), frame, FALSE, FALSE, 0);
vbox = gtk_vbox_new(FALSE, 0);
gtk_container_add(GTK_CONTAINER(frame), vbox);

vbox1 = gtk_vbox_new(TRUE, PANEL_SPACING);
gui_direct_check("Treat as a rigid body", &dock_rigid_on, dock_toggle_refresh, vbox1, vbox);
gtk_box_pack_start(GTK_BOX(vbox), vbox1, FALSE, FALSE, PANEL_SPACING);

gui_direct_check("Allow translation in x", &dock_rigid_x, NULL, NULL, vbox1);
gui_direct_check("Allow translation in y", &dock_rigid_y, NULL, NULL, vbox1);
gui_direct_check("Allow translation in z", &dock_rigid_z, NULL, NULL, vbox1);

/* control options */
frame = gtk_frame_new(NULL);
gtk_box_pack_start(GTK_BOX(vbox_right), frame, FALSE, FALSE, 0);
vbox = gtk_vbox_new(TRUE, 0);
gtk_container_add(GTK_CONTAINER(frame), vbox);
/*hbox =*/ gtk_hbox_new(FALSE, PANEL_SPACING);/*FIX aef23d*/
gui_direct_check("Create project files then stop", &dock_no_execute, NULL, NULL, vbox);

/* CURRRENT - executing the project has not been implemented yet */
gtk_widget_set_sensitive(GTK_WIDGET(vbox), FALSE);

/* control buttons */
gui_stock_button(GTK_STOCK_EXECUTE, docking_project_create, model,
                   GTK_DIALOG(window)->action_area);
gui_stock_button(GTK_STOCK_CLOSE, dialog_destroy, dialog,
                   GTK_DIALOG(window)->action_area);

gtk_widget_show_all(window);
}

