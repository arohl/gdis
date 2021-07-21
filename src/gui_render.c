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
#include <time.h>
#include <string.h>
#include <unistd.h>
#include <math.h>
#ifdef __APPLE__
#include <OpenGL/gl.h>
#else
#include <GL/gl.h>
#endif

#ifndef __WIN32
#include <sys/wait.h>
#endif

#include "gdis.h"
#include "coords.h"
#include "edit.h"
#include "file.h"
#include "parse.h"
#include "task.h"
#include "morph.h"
#include "matrix.h"
#include "opengl.h"
#include "render.h"
#include "select.h"
#include "spatial.h"
#include "zone.h"
#include "gui_shorts.h"
#include "interface.h"
#include "dialog.h"

#include <glib/gstdio.h>

#define DEBUG 0
#define DELETE_POV_FILE 0

extern struct sysenv_pak sysenv;
extern struct elem_pak elements[];
extern GtkWidget *window;

/* global render parameters */
struct light_pak current_light;
GtkWidget *scale_spin;
gdouble render_animate_frames = 100;
gint render_animate_overwrite = TRUE;

#define SCALE_MAG 10
#define PIX2ANG 0.03

/* this is an OpenGL limitation (ie shininess) */
#define MAX_HL 128

/****************/
/* render setup */
/****************/
void render_setup(GtkWidget *w, gpointer dummy)
{
/* checks */
if (!sysenv.povray_path)
  {
  gui_text_show(ERROR, "POVRay executable was not found.\n");
  return;
  }

/* single model background job */
povray_task();
}

/******************************/
/* font cancel button handler */
/******************************/
void font_dialog_close(GtkWidget *w, gpointer fsd)
{
gtk_widget_destroy(GTK_WIDGET(fsd));
}

/**************************/
/* font ok button handler */
/**************************/
void gl_font_selected(GtkWidget *w, gpointer fsd)
{
gchar *fn;

fn = gtk_font_selection_dialog_get_font_name ((GtkFontSelectionDialog *) fsd);
strcpy(sysenv.gl_fontname, fn);
redraw_canvas(ALL);
}

/**********************/
/* drawing font selection */
/**********************/
void gl_font_dialog(void)
{
GtkWidget *fsd;

fsd =  gtk_font_selection_dialog_new("Font selection");

/* setup events for the file selection widget */
g_signal_connect(GTK_OBJECT(GTK_FONT_SELECTION_DIALOG(fsd)->ok_button),
                 "clicked", (GtkSignalFunc) gl_font_selected, fsd);
g_signal_connect(GTK_OBJECT(GTK_FONT_SELECTION_DIALOG(fsd)->cancel_button),
                 "clicked", (GtkSignalFunc) font_dialog_close, fsd);

gtk_widget_show(fsd);
}

/************************************/
/* event handler for display dialog */
/************************************/
void render_refresh(void)
{
redraw_canvas(ALL);
}

/********************/
/* Property toggles */
/********************/
void toggle_axes_type(void)
{
struct model_pak *model;

model = sysenv.active_model;
if (!model)
 return;

if (model->periodic)
  {
  if (model->axes_type == CARTESIAN)
    model->axes_type = OTHER;
  else
    model->axes_type = CARTESIAN;
  }
coords_compute(model);
redraw_canvas(SINGLE);
}

void morph_toggle(GtkWidget *w, gpointer data)
{
redraw_canvas(SINGLE);
}

/****************/
/* unified hook */
/****************/
gint event_render_modify(GtkWidget *w, gpointer *obj)
{
gint id, refresh=0;
const gchar *entry;
#ifdef UNUSED_BUT_SET
struct model_pak *data;
#endif

/* checks */
g_return_val_if_fail(obj != NULL, FALSE);

/* ascertain type of modification required */
id = GPOINTER_TO_INT(g_object_get_data(G_OBJECT(obj), "id"));

#ifdef UNUSED_BUT_SET
data = sysenv.active_model;
#endif

sysenv.moving = FALSE;

switch(id)
  {
  case ANIM_GIF:
    sysenv.render.animate_type = ANIM_GIF;
    break;
  case ANIM_MPEG:
    sysenv.render.animate_type = ANIM_MPEG;
    break;
  case ANIM_NAME:
    entry = gtk_entry_get_text(GTK_ENTRY(obj));
    g_free(sysenv.render.animate_file);
    sysenv.render.animate_file = g_strdup(entry);
    break;
  case LIGHT_TYPE:
    entry = gtk_entry_get_text(GTK_ENTRY(GTK_COMBO(obj)->entry));
    if (g_ascii_strncasecmp(entry, "Positional", 10) == 0)
      current_light.type = POSITIONAL;
    else
      current_light.type = DIRECTIONAL;
    break;
  case MORPH_FINISH:
    entry = (gchar *) gtk_entry_get_text(GTK_ENTRY(obj));
    g_free(sysenv.render.morph_finish);
    sysenv.render.morph_finish = g_strdup(entry);
    refresh++;
    break;
  }

if (refresh)
  redraw_canvas(ALL);

return(FALSE);
}

/****************************/
/* rendering mode for model */
/****************************/
void render_mode_set(GtkWidget *w, gpointer data)
{
gint mode;
struct model_pak *model;

model = sysenv.active_model;
if (!model)
  return;

mode = GPOINTER_TO_INT(data);

spatial_destroy_by_label("polyhedra", model);
spatial_destroy_by_label("zones", model);

/* new mode becomes default for the model */
model->default_render_mode = mode;

/* selection based rendering */
if (model->selection)
  core_render_mode_set(mode, model->selection);
else
  core_render_mode_set(mode, model->cores);

redraw_canvas(SINGLE);
}

/*************************************/
/* display polyhedral representation */
/*************************************/
void render_mode_polyhedral(void)
{
struct model_pak *model;

model = sysenv.active_model;
if (!model)
  return;

create_polyhedra(model);
coords_init(REDO_COORDS, model);
redraw_canvas(SINGLE);
}

/*************************************/
/* display zone based representation */
/*************************************/
void render_mode_zone(void)
{
gpointer za;
struct model_pak *model;

model = sysenv.active_model;
if (!model)
  return;

core_render_mode_set(ZONE, model->cores);

/* CURRENT */
za = zone_make(sysenv.render.zone_size, model);
zone_display_init(za, model);
zone_free(za);

coords_init(REDO_COORDS, model);
redraw_canvas(SINGLE);
}

/*************************************/
/* wire-frame / solid atom rendering */
/*************************************/
void render_wire_atoms(void)
{
struct model_pak *model;

model = sysenv.active_model;
if (!model)
  return;

if (model->selection)
  core_render_wire_set(TRUE, model->selection);
else
  core_render_wire_set(TRUE, model->cores);

redraw_canvas(SINGLE);
}

void render_solid_atoms(void)
{
struct model_pak *model;

model = sysenv.active_model;
if (!model)
  return;

if (model->selection)
  core_render_wire_set(FALSE, model->selection);
else
  core_render_wire_set(FALSE, model->cores);

redraw_canvas(SINGLE);
}

/************************************/
/* toggle the render to file option */
/************************************/
/*
void toggle_animate(GtkWidget *w, GtkWidget *a_frame)
{
if (sysenv.render.animate)
  gtk_widget_set_sensitive(GTK_WIDGET(a_frame), TRUE);
else
  gtk_widget_set_sensitive(GTK_WIDGET(a_frame), FALSE);
}
*/

/*****************************/
/* task orientated rendering */
/*****************************/
void povray_exec_task(gpointer ptr,struct task_pak *task)
{
GString *cmd;
gchar *file;

g_return_if_fail(ptr != NULL);
cmd = g_string_new(NULL);

file = g_shell_quote(ptr);

/* build the command line */
g_string_printf(cmd, "%s +I%s -Ga -P +W%d +H%d +FT",//g_string_sprintf deprecated
                       sysenv.povray_path, file, 
                       (gint) sysenv.render.width, (gint) sysenv.render.height); 
if (sysenv.render.antialias)
  g_string_append_printf(cmd, " +A +AM2");//g_string_sprintfa deprecated

g_free(file);

/* don't do continuous display updates */
g_string_append_printf(cmd, " -D");//g_string_sprintfa deprecated

task->is_async = TRUE;
//task->status_file = g_build_filename(sysenv.cwd,???,NULL);/*what should be the status file?*/
task_async(cmd->str,&(task->pid));

g_string_free(cmd, TRUE);
}

/***************************/
/* task orientated viewing */
/***************************/
void exec_img_task(gpointer ptr)
{
gchar *tmp, *basename, *fullpath, *file_img;

/* post povray command */
g_return_if_fail(ptr != NULL);

/* construct basename (no extension) */
basename = parse_strip((gchar *) ptr);
fullpath = g_build_filename(sysenv.cwd, basename, NULL);
g_free(basename);

/* remove .pov input file */
if (sysenv.render.no_keep_tempfiles)
  {
  tmp = g_strdup_printf("%s.pov", fullpath);
  g_remove(tmp);
  g_free(tmp);
  }

/* construct image filename */
tmp = g_strdup_printf("%s.tga", fullpath);
file_img = g_shell_quote(tmp);
g_free(tmp);

/* construct viewing task command */
tmp = g_strdup_printf("%s %s", sysenv.viewer_path, file_img);
g_spawn_command_line_async(tmp, NULL);
g_free(tmp);

/* cleanup */
g_free(file_img);
g_free(fullpath);

/* this was the temporary filename which now must be cleaned up */
g_free(ptr);
}

/*************************/
/* foreground rendering  */
/*************************/
void povray_exec(gchar *name)
{
gchar *filename;
GString *cmd;

g_return_if_fail(name != NULL);
  return;

cmd = g_string_new(NULL);

filename = g_shell_quote(name);

/* build the command line */
g_string_printf(cmd, "%s +I%s -GA -P +W%d +H%d +FT ",//g_string_sprintf deprecated
                      sysenv.povray_path, filename,
                      (gint) sysenv.render.width , (gint) sysenv.render.height);
if (sysenv.render.antialias)
  g_string_append_printf(cmd,"+A +AM2 ");//g_string_sprintfa deprecated

/* after rendering delete input file, */
if (sysenv.render.no_keep_tempfiles)
  g_string_append_printf(cmd, "; rm -f %s ", filename);//g_string_sprintfa deprecated

/* don't do continuous display updates */
g_string_append_printf(cmd, " -D");//g_string_sprintfa deprecated

/* execute */
printf("system: %s\n", cmd->str);


IGNORE_RETURN(system(cmd->str));
printf("\n");

/* cleanup */
g_string_free(cmd, TRUE);
g_free(filename);
}

/************************/
/* background rendering */
/************************/
void povray_task(void)
{
gchar *basename, *fullname;
struct model_pak *model;

model = sysenv.active_model;
if (!model)
  return;

/* make an input file */
basename = gun("pov");
if (!basename)
  {
  printf("povray_task() error: failed to get an unused filename.\n");
  return;
  }

/* create the input file */
fullname = g_build_filename(sysenv.cwd, basename, NULL);
write_povray(fullname, model);
g_free(fullname);

if (sysenv.render.no_povray_exec)
  {
/* no further tasks, can free basename */
  g_free(basename);
  }
else
  {
/* basename should be freed by second (cleanup) routine */
  task_new("POVRay", &povray_exec_task, basename, &exec_img_task, basename, NULL);
  }
}

/*************/
/* callbacks */
/*************/
void geom_label_toggle(GtkWidget *w, gpointer dummy)
{
struct model_pak *model = sysenv.active_model;

if (model)
  {
  model->show_geom_labels ^= 1;
  redraw_canvas(SINGLE);
  }
}

void update_geom_line_width(GtkWidget *w, gpointer *ptr)
{
sysenv.render.geom_line_width = GTK_ADJUSTMENT(w)->value;
redraw_canvas(SINGLE);
}

/****************************************/
/* active model display property widget */
/****************************************/
void gui_display_widget(GtkWidget *box)
{
GtkWidget *vbox;

g_assert(box != NULL);

vbox = gui_frame_vbox(NULL, FALSE, FALSE, box);

gui_button_x("Ball & Stick", render_mode_set, GINT_TO_POINTER(BALL_STICK), vbox);
gui_button_x("CPK", render_mode_set, GINT_TO_POINTER(CPK), vbox);
gui_button_x("Liquorice", render_mode_set, GINT_TO_POINTER(LIQUORICE), vbox);
gui_button_x("Polyhedral", render_mode_polyhedral, NULL,  vbox);
gui_button_x("Stick", render_mode_set, GINT_TO_POINTER(STICK), vbox);
gui_button_x("Zone based", render_mode_zone, NULL, vbox);
gui_button_x("Wire frame", render_wire_atoms, NULL, vbox);
gui_button_x("Solid", render_solid_atoms, NULL, vbox);
}

/***********************/
/* main rendering page */
/***********************/
void render_main_page(GtkWidget *box)
{
GtkWidget *vbox1, *vbox2, *vbox, *hbox;

/* left & right pane split */
hbox = gtk_hbox_new(FALSE, PANEL_SPACING);
gtk_container_add(GTK_CONTAINER(box), hbox);
gtk_container_set_border_width(GTK_CONTAINER(hbox), PANEL_SPACING);
vbox1 = gtk_vbox_new(FALSE, PANEL_SPACING);
gtk_box_pack_start(GTK_BOX(hbox), vbox1, TRUE, TRUE, 0);
vbox2 = gtk_vbox_new(FALSE, PANEL_SPACING);
gtk_box_pack_start(GTK_BOX(hbox), vbox2, TRUE, TRUE, 0);

/* left pane */
vbox = gui_frame_vbox(NULL, FALSE, FALSE, vbox1);

gui_button_x("Ball & Stick", render_mode_set, GINT_TO_POINTER(BALL_STICK), vbox);
gui_button_x("CPK", render_mode_set, GINT_TO_POINTER(CPK), vbox);
gui_button_x("Liquorice", render_mode_set, GINT_TO_POINTER(LIQUORICE), vbox);
gui_button_x("Polyhedral", render_mode_polyhedral, NULL,  vbox);
gui_button_x("Stick", render_mode_set, GINT_TO_POINTER(STICK), vbox);
gui_button_x("Zone based (experimental)", render_mode_zone, NULL, vbox);
gui_button_x("Wire frame molecules", render_wire_atoms, NULL, vbox);
gui_button_x("Solid molecules", render_solid_atoms, NULL, vbox);
gui_button_x("Switch axes type", toggle_axes_type, NULL, vbox);

/* extra toggles */
vbox = gui_frame_vbox(NULL, FALSE, FALSE, vbox1);

gui_direct_check("Antialias", &sysenv.render.antialias,
                   render_refresh, NULL, vbox);
gui_direct_check("CPK scaling for B&S atoms", &sysenv.render.scale_ball_size,
                   render_refresh, NULL, vbox);
gui_direct_check("Wire frame surfaces", &sysenv.render.wire_surface,
                   render_refresh, NULL, vbox);
gui_direct_check("Show hidden surfaces", &sysenv.render.wire_show_hidden,
                   render_refresh, NULL, vbox);

/* surface transmission frame */
vbox = gui_frame_vbox(NULL, FALSE, FALSE, vbox1);

/* ghost atoms and shells */
gui_direct_spin("Ghost atom opacity",
                  &sysenv.render.ghost_opacity, 0.0, 1.0, 0.1,
                  render_refresh, NULL, vbox);
/* molecular surfaces */
gui_direct_spin("Surface opacity",
                  &sysenv.render.transmit, 0.0, 1.0, 0.1,
                  render_refresh, NULL, vbox);

/* EXP - zone based display */
vbox = gui_frame_vbox(NULL, FALSE, FALSE, vbox1);

gui_direct_spin("Zone grid size ",
                  &sysenv.render.zone_size, 1.0, 1000.0, 1.0,
                  render_refresh, NULL, vbox);

/* radii frame */
vbox = gui_frame_vbox(NULL, FALSE, FALSE, vbox2);

gui_direct_spin("Ball radius",
                  &sysenv.render.ball_radius, 0.1, 0.5, 0.02,
                  render_refresh, NULL, vbox);

gui_direct_spin("Cylinder radius",
                  &sysenv.render.stick_radius, 0.02, 0.5, 0.01,
                  render_refresh, NULL, vbox);

gui_direct_spin("Stick thickness",
                  &sysenv.render.stick_thickness, 0.1, 9.0, 0.05,
                  render_refresh, NULL, vbox);

gui_direct_spin("Line drawing width ",
                  &sysenv.render.frame_thickness, 0.1, 9.0, 0.05,
                  render_refresh, NULL, vbox);

gui_direct_spin("CPK scaling",
                  &sysenv.render.cpk_scale, 0.1, 3.0, 0.02,
                  render_refresh, NULL, vbox);

/* highlighting frame */
vbox = gui_frame_vbox(NULL, FALSE, FALSE, vbox2);

gui_direct_spin("Atom highlight power",
                   &sysenv.render.ahl_strength, 0.0, 1.0, 0.05,
                   render_refresh, NULL, vbox);

gui_direct_spin("Atom highlight focus",
                   &sysenv.render.ahl_size, 0.0, MAX_HL, 5.0,
                   render_refresh, NULL, vbox);

gui_direct_spin("Surface highlight power",
                   &sysenv.render.shl_strength, 0.0, 1.0, 0.05,
                   render_refresh, NULL, vbox);

gui_direct_spin("Surface highlight focus",
                   &sysenv.render.shl_size, 0.0, MAX_HL, 5.0,
                   render_refresh, NULL, vbox);

/* ribbon control frame */
vbox = gui_frame_vbox(NULL, FALSE, FALSE, vbox2);

gui_direct_spin("Ribbon curvature control",
                  &sysenv.render.ribbon_curvature, 0.0, 1.0, 0.1,
                  render_refresh, NULL, vbox);

gui_direct_spin("Ribbon thickness",
                  &sysenv.render.ribbon_thickness, 0.4, 6.0, 0.2,
                  render_refresh, NULL, vbox);

gui_direct_spin("Ribbon quality",
                  &sysenv.render.ribbon_quality, 1.0, 30.0, 1.0,
                  render_refresh, NULL, vbox);
}

/************************************/
/* colouring method change callback */
/************************************/
void set_colour_scheme(GtkWidget *w, gpointer *data)
{
const gchar *tmp;
struct model_pak *model;

tmp = gtk_entry_get_text(GTK_ENTRY(data));
model = sysenv.active_model;
if (!model)
  return;

/* default to atom selection */
if (g_ascii_strncasecmp(tmp, "element", 7) == 0)
  model_colour_scheme(ELEM, model);

if (g_ascii_strncasecmp(tmp, "growth", 6) == 0)
  model_colour_scheme(GROWTH_SLICE, model);

if (g_ascii_strncasecmp(tmp, "translation", 11) == 0)
  model_colour_scheme(TRANSLATE, model);

if (g_ascii_strncasecmp(tmp, "molecule", 8) == 0)
  model_colour_scheme(MOL, model);

if (g_ascii_strncasecmp(tmp, "occupancy", 8) == 0)
  model_colour_scheme(OCCUPANCY, model);

if (g_ascii_strncasecmp(tmp, "region", 6) == 0)
  model_colour_scheme(REGION, model);

if (g_ascii_strncasecmp(tmp, "temp", 4) == 0)
  model_colour_scheme(VELOCITY, model);

redraw_canvas(SINGLE);
}

/**************************/
/* colour control options */
/**************************/
void render_colours_section(GtkWidget *box)
{
gpointer entry;
GList *list=NULL;
GtkWidget *vbox1, *vbox2, *vbox, *hbox;

/* left & right pane split */
hbox = gtk_hbox_new(FALSE, PANEL_SPACING);
gtk_container_add(GTK_CONTAINER(box), hbox);
vbox1 = gtk_vbox_new(FALSE, PANEL_SPACING);
gtk_box_pack_start(GTK_BOX(hbox), vbox1, FALSE, FALSE, 0);
vbox2 = gtk_vbox_new(FALSE, PANEL_SPACING);
gtk_box_pack_start(GTK_BOX(hbox), vbox2, TRUE, TRUE, 0);

/* next frame */
vbox = gui_frame_vbox(NULL, FALSE, FALSE, vbox1);

gui_colour_box("Background colour ", sysenv.render.bg_colour, vbox);
/*
gui_colour_box("Crystal colour ", sysenv.render.morph_colour, vbox);
*/
gui_colour_box("Re-entrant colour ", sysenv.render.rsurf_colour, vbox);
gui_colour_box("Ribbon colour ", sysenv.render.ribbon_colour, vbox);

/* next frame */
vbox = gui_frame_vbox(NULL, FALSE, FALSE, vbox1);

/* colouring scheme */
hbox = gtk_hbox_new(FALSE, PANEL_SPACING);
gtk_box_pack_start(GTK_BOX(vbox), hbox, FALSE, FALSE, 0);

list = NULL;
list = g_list_append(list, "element");
list = g_list_append(list, "growth slice");
list = g_list_append(list, "molecule");
list = g_list_append(list, "occupancy");
list = g_list_append(list, "region");
list = g_list_append(list, "temperature");
list = g_list_append(list, "translation");
entry = gui_pulldown_new("Colouring scheme ", list, FALSE, hbox);
gui_button_x(NULL, set_colour_scheme, entry, hbox);
}

/******************/
/* camera globals */
/******************/
GtkListStore *camera_list=NULL;
GtkWidget *camera_tree=NULL;
GtkWidget *camera_mode_entry=NULL, *camera_zoom_entry=NULL, *camera_proj_entry=NULL;
GtkWidget *camera_fov_entry=NULL, *camera_axis_entry=NULL;
struct camera_pak *camera_current=NULL;

gdouble render_animate_angle[3] = {0.0, 360.0, 5.0};

/*****************************/
/* update camera info values */
/*****************************/
void gui_camera_update(gpointer data)
{
gchar *text;
struct model_pak *model;
struct camera_pak *camera = data;

/* checks */
model = sysenv.active_model;
if (!model)
  return;
if (!camera)
  camera = model->camera;
g_assert(camera != NULL);

camera_current = camera;

if (GTK_IS_ENTRY(camera_proj_entry))
  {
  if (camera->perspective)
    gtk_entry_set_text(GTK_ENTRY(camera_proj_entry), "perspective");
  else
    gtk_entry_set_text(GTK_ENTRY(camera_proj_entry), "orthographic");
  }

if (GTK_IS_ENTRY(camera_fov_entry))
  {
  if (camera->perspective)
    {
    text = g_strdup_printf("%.0f", camera->fov);
    gtk_entry_set_text(GTK_ENTRY(camera_fov_entry), text);
    g_free(text);
    gtk_widget_set_sensitive(GTK_WIDGET(camera_fov_entry), TRUE);
    }
  else
    {
    gtk_entry_set_text(GTK_ENTRY(camera_fov_entry), "180");
    gtk_widget_set_sensitive(GTK_WIDGET(camera_fov_entry), FALSE);
    }
  }

if (GTK_IS_ENTRY(camera_mode_entry))
  {
  if (camera->mode == LOCKED)
    gtk_entry_set_text(GTK_ENTRY(camera_mode_entry), "origin");
  else
    gtk_entry_set_text(GTK_ENTRY(camera_mode_entry), "variable");
  }

if (GTK_IS_ENTRY(camera_zoom_entry))
  {
  if (camera->perspective)
    {
    gtk_entry_set_text(GTK_ENTRY(camera_zoom_entry), "not applicable");
    gtk_widget_set_sensitive(GTK_WIDGET(camera_zoom_entry), FALSE);
    }
  else
    {
    text = g_strdup_printf("%.1f", 1000.0 / (model->rmax * camera->zoom));
    gtk_entry_set_text(GTK_ENTRY(camera_zoom_entry), text);
    g_free(text);
    gtk_widget_set_sensitive(GTK_WIDGET(camera_zoom_entry), TRUE);
    }
  }
}

/*******************************/
/* manual zoom change callback */
/*******************************/
void gui_camera_zoom_changed(GtkWidget *w, gpointer data)
{
gdouble z;
const gchar *text;
struct model_pak *model;
struct camera_pak *camera;

/* checks */
model = sysenv.active_model;
if (!model)
  return;
camera = model->camera;

if (GTK_IS_ENTRY(w) && camera)
  {
  text = gtk_entry_get_text(GTK_ENTRY(w));

  z = str_to_float(text);

  camera->zoom = 1000.0 / (model->rmax * z);

  gui_refresh(GUI_CANVAS);
  }
}

/*****************************/
/* waypoint animate callback */
/*****************************/
void render_waypoint_animate(GtkWidget *w, gpointer data)
{
gint frames = render_animate_frames;
struct model_pak *model;

model = sysenv.active_model;
if (model)
  camera_waypoint_animate(frames, render_animate_overwrite, model);
}

/*****************************/
/* rotation animate callback */
/*****************************/
void render_rotate_animate(GtkWidget *w, gpointer data)
{
gint axis=2;
const gchar *text;
struct model_pak *model;

model = sysenv.active_model;

if (model)
  {
  text = gtk_entry_get_text(GTK_ENTRY(camera_axis_entry));

  if (g_ascii_strncasecmp(text, "x ", 2) == 0)
    axis = 0;
  if (g_ascii_strncasecmp(text, "y ", 2) == 0)
    axis = 1;
  if (g_ascii_strncasecmp(text, "z ", 2) == 0)
    axis = 2;

  camera_rotate_animate(axis, render_animate_angle, render_animate_overwrite, model);
  }
}

/***********************************/
/* camera waypoint delete callback */
/***********************************/
void gui_delete_waypoint(GtkWidget *w, gpointer data)
{
GtkTreeIter iter;
GtkTreeModel *treemodel;
GtkTreeSelection *selection;
struct camera_pak *camera;
struct model_pak *model;

/* checks */
if (GTK_IS_TREE_VIEW(data))
  selection = gtk_tree_view_get_selection(GTK_TREE_VIEW(data));
else
  return;
model = sysenv.active_model;
if (!model)
  return;

/* record selection as the active folder */
if (gtk_tree_selection_get_selected(selection, &treemodel, &iter))
  {
  gtk_tree_model_get(treemodel, &iter, 1, &camera, -1);
  if (camera == model->camera_default)
    {
    gui_text_show(ERROR, "You cannot delete the default camera.\n");
    return;
    }
  else
    {
    if (camera == model->camera)
      model->camera = model->camera_default;
    model->waypoint_list = g_slist_remove(model->waypoint_list, camera);
    gui_camera_refresh();
    redraw_canvas(SINGLE);
    }
  }
}

/**************************************/
/* fill out the list of model cameras */
/**************************************/
void gui_camera_populate(void)
{
gint n;
gchar *text;
gpointer camera;
GSList *list;
GtkTreeIter iter;
GtkTreeSelection *selection;
struct model_pak *model;

gtk_list_store_clear(camera_list);

model = sysenv.active_model;
if (!model)
  return;

/* default camera */
n=0;
gtk_list_store_append(camera_list, &iter);
gtk_list_store_set(camera_list, &iter, 0, "default", 1, model->camera_default, -1);

/* select camera if it's currently active */
selection = gtk_tree_view_get_selection(GTK_TREE_VIEW(camera_tree)); 
if (selection)
  {
  if (model->camera == model->camera_default)
    {
    gtk_tree_selection_select_iter(selection, &iter);
/* don't bother checking waypoint list (below) */
    selection = NULL;
    }
  }

/* waypoints */
for (list=model->waypoint_list ; list ; list=g_slist_next(list))
  {
  camera = list->data;
  text = g_strdup_printf("waypoint %d", ++n);
  gtk_list_store_append(camera_list, &iter);
  gtk_list_store_set(camera_list, &iter, 0, text, 1, camera, -1);
  g_free(text);

/* select camera if it's currently active */
  if (selection)
    {
    if (model->camera == camera)
      gtk_tree_selection_select_iter(selection, &iter);
    }
  }
}

/******************************/
/* select a particular camera */
/******************************/
void gui_camera_select(GtkTreeSelection *selection, gpointer data)
{
GtkTreeIter iter;
GtkTreeModel *treemodel;
struct camera_pak *camera;
struct model_pak *model;

/* checks */
model = sysenv.active_model;
if (!model)
  return;

/* record selection as the active folder */
if (gtk_tree_selection_get_selected(selection, &treemodel, &iter))
  {
  gtk_tree_model_get(treemodel, &iter, 1, &camera, -1);

  gui_camera_update(camera);
  model->camera = camera;
  redraw_canvas(SINGLE);
  }
}

/******************************/
/* camera display page update */
/******************************/
void gui_camera_refresh(void)
{
if (dialog_exists(POVRAY, NULL))
  {
  gui_camera_populate();
  gui_camera_update(NULL);
  }
}

/*************************************/
/* camera projection modify callback */
/*************************************/
void gui_camera_projection_changed(GtkWidget *w, gpointer data)
{
const gchar *text;

/* checks */
g_assert(GTK_IS_ENTRY(w));
if (!camera_current)
  return;

text = gtk_entry_get_text(GTK_ENTRY(w));

/* NB: sometimes an empty string is passed - don't update unless we have to */
if (g_ascii_strncasecmp(text, "ortho", 5) == 0)
  {
  camera_current->perspective = FALSE;
  redraw_canvas(SINGLE);
  gui_camera_update(camera_current);
  }
if (g_ascii_strncasecmp(text, "persp", 5) == 0)
  {
  camera_current->perspective = TRUE;
  redraw_canvas(SINGLE);
  gui_camera_update(camera_current);
  }
}

/*******************************/
/* camera mode modify callback */
/*******************************/
void gui_camera_mode_changed(GtkWidget *w, gpointer data)
{
const gchar *text;

/* checks */
g_assert(GTK_IS_ENTRY(w));
if (!camera_current)
  return;

text = gtk_entry_get_text(GTK_ENTRY(w));

/* NB: sometimes an empty string is passed - don't update unless we have to */
if (g_ascii_strncasecmp(text, "orig", 4) == 0)
  {
  camera_current->mode = LOCKED;
  redraw_canvas(SINGLE);
  }
if (g_ascii_strncasecmp(text, "vari", 4) == 0)
  {
  camera_current->mode = FREE;
  redraw_canvas(SINGLE);
  }
}

/***********************************/
/* camera waypoint create callback */
/***********************************/
void gui_create_waypoint(void)
{
struct model_pak *model;
gpointer camera;

model = sysenv.active_model;
if (model)
  {
/* allocate */
  camera = camera_dup(model->camera);
  model->waypoint_list = g_slist_append(model->waypoint_list, camera);

/* move new camera slightly along the viewing vector */
/*
  ARR3SET(v, camera->v);
  VEC3MUL(v, 0.1);
  ARR3ADD(camera->x, v);
*/

  model->camera = camera;
  
  gui_camera_refresh();
  }
}

/****************************************************/
/* camera manipulation (waypoints, projection etc.) */
/****************************************************/
void render_camera_section(GtkWidget *box)
{
GList *list;
GtkCellRenderer *r;
GtkTreeViewColumn *c;
GtkTreeSelection *select;
GtkWidget *vbox1, *vbox2, *vbox, *hbox;
GtkWidget *swin, *label, *combo;

/* left & right pane split */
hbox = gtk_hbox_new(FALSE, PANEL_SPACING);
gtk_container_add(GTK_CONTAINER(box), hbox);
vbox1 = gtk_vbox_new(FALSE, PANEL_SPACING);
gtk_box_pack_start(GTK_BOX(hbox), vbox1, TRUE, TRUE, 0);
vbox2 = gtk_vbox_new(FALSE, PANEL_SPACING);
gtk_box_pack_start(GTK_BOX(hbox), vbox2, TRUE, TRUE, 0);

/* CURRENT - camera info */
vbox = gui_frame_vbox(NULL, FALSE, FALSE, vbox1);

/* projection */
list = NULL;
list = g_list_prepend(list, "perspective");
list = g_list_prepend(list, "orthographic");
camera_proj_entry = gui_pulldown_new("projection ", list, FALSE, vbox);

/* focal point */
list = NULL;
list = g_list_prepend(list, "variable");
list = g_list_prepend(list, "origin");
camera_mode_entry = gui_pulldown_new("focal point", list, FALSE, vbox);

/* field of view */
camera_fov_entry =  gui_text_entry("field of view", NULL, FALSE, FALSE, vbox);

/* zoom */
camera_zoom_entry =  gui_text_entry("zoom factor", NULL, TRUE, FALSE, vbox);
g_signal_connect(GTK_OBJECT(camera_zoom_entry), "activate",
                 GTK_SIGNAL_FUNC(gui_camera_zoom_changed),  NULL);

/* frame - animation from waypoints */
vbox = gui_frame_vbox(NULL, FALSE, FALSE, vbox1);

gui_button_x("Add waypoint ", gui_create_waypoint, NULL, vbox);
gui_button_x("Delete waypoint ", gui_delete_waypoint, camera_tree, vbox);

gui_button_x("Create animation from waypoints ", render_waypoint_animate, NULL, vbox);

gui_direct_spin("Frames for each waypoint traversal ",
                  &render_animate_frames, 1.0, 1000.0, 1,
                  NULL, NULL, vbox);


/* frame - animation from axis rotations */
vbox = gui_frame_vbox(NULL, FALSE, FALSE, vbox1);

hbox = gtk_hbox_new(FALSE, PANEL_SPACING);
gtk_box_pack_start(GTK_BOX(vbox), hbox, FALSE, FALSE, 0);

label = gtk_label_new("Rotation vector ");
gtk_box_pack_start(GTK_BOX(hbox), label, FALSE, FALSE, 0);

list = NULL;
list = g_list_prepend(list, "z axis");
list = g_list_prepend(list, "y axis");
list = g_list_prepend(list, "x axis");
combo = gtk_combo_new();
camera_axis_entry = GTK_COMBO(combo)->entry;
gtk_entry_set_editable(GTK_ENTRY(camera_axis_entry), FALSE);
gtk_combo_set_popdown_strings(GTK_COMBO(combo), list);
gtk_box_pack_end(GTK_BOX(hbox), combo, FALSE, FALSE, 0);

gui_direct_spin("Start ", &render_animate_angle[0], 0, 359, 5, NULL, NULL, vbox);
gui_direct_spin("Stop", &render_animate_angle[1], 1, 360, 5, NULL, NULL, vbox);
gui_direct_spin("Increment", &render_animate_angle[2], 1, 10, 1, NULL, NULL, vbox);

gui_button_x("Create animated rotation ", render_rotate_animate, NULL, vbox);


/* frame - animation from axis rotations */
vbox = gui_frame_vbox(NULL, FALSE, FALSE, vbox1);

gui_direct_check("Overwrite old frames ", &render_animate_overwrite, NULL, NULL, vbox);

/* right frame */
vbox = gui_frame_vbox(NULL, FALSE, FALSE, vbox2);

/* shortcut buttons for camera positioning */
/*
hbox = gtk_hbox_new(TRUE, PANEL_SPACING);
gtk_box_pack_start(GTK_BOX(vbox), hbox, FALSE, FALSE, 0);
gui_button(" x ", gui_view_x, NULL, hbox, TT);
gui_button(" y ", gui_view_y, NULL, hbox, TT);
gui_button(" z ", gui_view_z, NULL, hbox, TT);
gui_button(" a ", gui_view_a, NULL, hbox, TT);
gui_button(" b ", gui_view_b, NULL, hbox, TT);
gui_button(" c ", gui_view_c, NULL, hbox, TT);
*/

/* list of cameras & waypoints */
swin = gtk_scrolled_window_new(NULL, NULL);
gtk_scrolled_window_set_policy(GTK_SCROLLED_WINDOW(swin),
                               GTK_POLICY_AUTOMATIC, GTK_POLICY_AUTOMATIC);
gtk_box_pack_start(GTK_BOX(vbox), swin, TRUE, TRUE, 0);

/* list */
camera_list = gtk_list_store_new(2, G_TYPE_STRING, G_TYPE_POINTER);
camera_tree = gtk_tree_view_new_with_model(GTK_TREE_MODEL(camera_list));
gtk_scrolled_window_add_with_viewport(GTK_SCROLLED_WINDOW(swin), camera_tree);
r = gtk_cell_renderer_text_new();
c = gtk_tree_view_column_new_with_attributes(" ", r, "text", 0, NULL);
gtk_tree_view_append_column(GTK_TREE_VIEW(camera_tree), c);
gtk_tree_view_set_headers_visible(GTK_TREE_VIEW(camera_tree), FALSE);

/* camera selection handler */
select = gtk_tree_view_get_selection(GTK_TREE_VIEW(camera_tree));
gtk_tree_selection_set_mode(select, GTK_SELECTION_SINGLE);
g_signal_connect(G_OBJECT(select), "changed",
                 G_CALLBACK(gui_camera_select),
                 NULL);

/* add/delete waypoint buttons */
/*
hbox = gtk_hbox_new(FALSE, PANEL_SPACING);
gtk_box_pack_start(GTK_BOX(vbox), hbox, FALSE, FALSE, 0);
gui_button("Add", gui_create_waypoint, NULL, hbox, TT);
gui_button("Delete", gui_delete_waypoint, camera_tree, hbox, TT);
*/

/* TODO - a mechanism so that whenever a switch_model happens a */
/* function can be automatically called to update stuff */
/* AND is destroyed whenver the widget it's attached to is destroyed */

/* fill out the dialog */
gui_camera_refresh();

/* now attach modify callbacks */
g_signal_connect(GTK_OBJECT(camera_proj_entry), "changed", 
                 GTK_SIGNAL_FUNC(gui_camera_projection_changed), NULL);
g_signal_connect(GTK_OBJECT(camera_mode_entry), "changed", 
                 GTK_SIGNAL_FUNC(gui_camera_mode_changed), NULL);
}

/**********************/
/* light list globals */
/**********************/
enum 
{
RENDER_LIGHT_COLOUR,
RENDER_LIGHT_VECTOR,
RENDER_LIGHT_TYPE,
RENDER_LIGHT_ATTRIBUTE,
RENDER_LIGHT_NCOLS
};

GtkWidget *render_light_tv;
GtkListStore *render_light_ls;


/****************************************************/
/* convert 3 doubles [0:1] into a hex colour string */
/****************************************************/
gchar *get_hex_colour(gdouble *colour)
{
guint r, g, b;
gdouble rgb[3];

ARR3SET(rgb, colour);
VEC3MUL(rgb, 255.0);

r = rgb[0];
g = rgb[1];
b = rgb[2];

return(g_strdup_printf("#%2X%2X%2X", r, g, b));
}

/**********************************/
/* get the currently selected row */
/**********************************/
gint render_light_selected(void)
{
gint row=0;
GtkTreeModel *treemodel;
GtkTreeSelection *selection;
GtkTreeIter iter;

treemodel = gtk_tree_view_get_model(GTK_TREE_VIEW(render_light_tv));
selection = gtk_tree_view_get_selection(GTK_TREE_VIEW(render_light_tv));

if (gtk_tree_model_get_iter_first(treemodel, &iter))
  {
  do
    {
    if (gtk_tree_selection_iter_is_selected(selection, &iter))
      return(row);
    row++;
    }
  while (gtk_tree_model_iter_next(treemodel, &iter));
  }
return(-1);
}

/***********************************/
/* construct the light source list */
/***********************************/
void update_light_list(void)
{
gint row, num_rows=0;
gchar *text;
GSList *list;
GtkTreeIter iter;
GtkTreeModel *treemodel;
GtkTreeSelection *selection;
struct light_pak *light;

/* checks */
g_assert(render_light_ls != NULL);

/* store */
row = render_light_selected();

/* re-populate */
gtk_list_store_clear(render_light_ls);
for (list=sysenv.render.light_list ; list ; list=g_slist_next(list))
  {
  light = list->data;
  gtk_list_store_append(render_light_ls, &iter);

/* position or direction vector */
  text = g_strdup_printf(" (%5.1f, %5.1f, %5.1f) ", light->x[0],
                                                    light->x[1],
                                                    light->x[2]);
  gtk_list_store_set(render_light_ls, &iter, RENDER_LIGHT_VECTOR, text, -1);
  g_free(text);

/* type */
  if (light->type == DIRECTIONAL)
    text = g_strdup(" Directional");
  else
    text = g_strdup(" Positional");
  gtk_list_store_set(render_light_ls, &iter, RENDER_LIGHT_TYPE, text, -1);
  g_free(text);

/* FIXME - easier way to specify the RGB colour? */
  text = get_hex_colour(light->colour);
  gtk_list_store_set(render_light_ls, &iter, RENDER_LIGHT_ATTRIBUTE, text, -1);
  g_free(text);

  num_rows++;
  }

/* restore selected row, or the previous (if possible) if deleted  */
if (row >= num_rows && row)
  row--;
if (row >= 0)
  {
  treemodel = gtk_tree_view_get_model(GTK_TREE_VIEW(render_light_tv));
  if (gtk_tree_model_iter_nth_child(treemodel, &iter, NULL, row))
    {
    selection = gtk_tree_view_get_selection(GTK_TREE_VIEW(render_light_tv)); 
    if (selection)
      gtk_tree_selection_select_iter(selection, &iter);
    }
  }
}

/*******************************************/
/* light list manipulation, add new source */
/*******************************************/
void add_current_light(GtkWidget *w, gpointer dummy)
{
struct light_pak *light;

/* duplicate data for the list */
light = g_malloc(sizeof(struct light_pak));
memcpy(light, &current_light, sizeof(struct light_pak));

/* append & update */
sysenv.render.light_list = g_slist_append(sysenv.render.light_list, light);
update_light_list();

redraw_canvas(SINGLE);
}

/******************************************/
/* light list manipulation, modify source */
/******************************************/
void mod_selected_light(GtkWidget *w, gpointer dummy)
{
gint row;
struct light_pak *light;

row = render_light_selected();
if (row < 0)
  return;

/* get the light's data from the list */
light = g_slist_nth_data(sysenv.render.light_list, row);

/* overwrite with the current data */
memcpy(light, &current_light, sizeof(struct light_pak));

update_light_list();

redraw_canvas(SINGLE);
}

/******************************************/
/* light list manipulation, delete source */
/******************************************/
void del_selected_light(GtkWidget *w, gpointer dummy)
{
gint row;
struct light_pak *light;

row = render_light_selected();
if (row < 0)
  return;

light = g_slist_nth_data(sysenv.render.light_list, row);

sysenv.render.light_list = g_slist_remove(sysenv.render.light_list, light);

update_light_list();

redraw_canvas(SINGLE);
}

/******************/
/* render cleanup */
/******************/
void render_cleanup(void)
{
gtk_list_store_clear(render_light_ls);
render_light_ls = NULL;
render_light_tv = NULL;
}

/*********************************************/
/* callback to fill dialog with light values */
/*********************************************/
void render_light_activate(GtkTreeView *treeview, GtkTreePath *treepath)
{
gint n;
struct light_pak *light;

/* get selected light */
n = render_light_selected();
if (n < 0)
  return;

/* get the light's data from the list */
light = g_slist_nth_data(sysenv.render.light_list, n);

/* transfer to current displayed values */
if (light)
  {
/* overwrite */
  memcpy(&current_light, light, sizeof(struct light_pak));

/* update */
  gui_relation_update(NULL);
  }
}

/****************************/
/* lighting control options */
/****************************/
void render_lighting_section(GtkWidget *box)
{
gint i;
gchar *titles[] = {" Colour ", "  Absolute vector  ", " Type "};
GList *list=NULL;
GtkWidget *swin, *vbox1, *vbox2, *vbox, *hbox;
GtkWidget *label, *combo;
GtkCellRenderer *renderer;
GtkTreeViewColumn *column;

/* left & right pane split */
hbox = gtk_hbox_new(FALSE, PANEL_SPACING);
gtk_container_add(GTK_CONTAINER(box), hbox);
vbox1 = gtk_vbox_new(FALSE, PANEL_SPACING);
gtk_box_pack_start(GTK_BOX(hbox), vbox1, FALSE, FALSE, 0);
vbox2 = gtk_vbox_new(FALSE, PANEL_SPACING);
gtk_box_pack_start(GTK_BOX(hbox), vbox2, TRUE, TRUE, 0);

/* box 1 - colour (and edit button) */
vbox = gui_frame_vbox(NULL, FALSE, FALSE, vbox1);

gui_colour_box("Current light source colour ", current_light.colour, vbox);

/* box 2 - light components */
vbox = gui_frame_vbox(NULL, FALSE, FALSE, vbox1);

hbox = gtk_hbox_new(FALSE, 0);
gtk_box_pack_start(GTK_BOX(vbox), hbox, FALSE, TRUE, 0);

gui_direct_spin("Ambient component", &current_light.ambient,
                  0.0, 1.0, 0.1, NULL, NULL, vbox);

gui_direct_spin("Diffuse component", &current_light.diffuse,
                  0.0, 1.0, 0.1, NULL, NULL, vbox);

gui_direct_spin("Specular component", &current_light.specular,
                  0.0, 1.0, 0.1, NULL, NULL, vbox);


/* left box - light type */
vbox = gui_frame_vbox(NULL, FALSE, FALSE, vbox1);

/* type label */
hbox = gtk_hbox_new(FALSE, 0);
gtk_box_pack_start(GTK_BOX(vbox), hbox, FALSE, TRUE, 0);
label = gtk_label_new("  Type  ");
gtk_box_pack_start(GTK_BOX(hbox), label, FALSE, FALSE, 0);
/* TODO - pulldown dir/pos */
/* NEW - combo box for the selection mode */
list = g_list_append(list, "Directional");
list = g_list_append(list, "Positional");
combo = gtk_combo_new();
gtk_entry_set_editable(GTK_ENTRY(GTK_COMBO(combo)->entry), FALSE);
gtk_combo_set_popdown_strings(GTK_COMBO(combo), list);
gtk_box_pack_start(GTK_BOX(hbox), combo, FALSE, FALSE, 0);
g_signal_connect(GTK_OBJECT(GTK_COMBO(combo)->entry), "changed", 
                 GTK_SIGNAL_FUNC(event_render_modify), (gpointer) combo);
g_object_set_data(G_OBJECT(combo), "id", (gpointer) LIGHT_TYPE);


/* position/direction vector input */
hbox = gtk_hbox_new(FALSE, 0);
gtk_box_pack_start(GTK_BOX(vbox), hbox, TRUE, TRUE, 0);

/* x component */
vbox = gtk_vbox_new(TRUE, 0);
gtk_box_pack_start(GTK_BOX(hbox), vbox, TRUE, TRUE, 0);
label = gtk_label_new("X");
gtk_box_pack_start(GTK_BOX(vbox), label, TRUE, TRUE, 0);

gui_direct_spin(NULL, &current_light.x[0],
                  -1000.0, 1000.0, 0.1, NULL, NULL, vbox);

/* y component */
vbox = gtk_vbox_new(TRUE, 0);
gtk_box_pack_start(GTK_BOX(hbox), vbox, TRUE, TRUE, 0);
label = gtk_label_new("Y");
gtk_box_pack_start(GTK_BOX(vbox), label, TRUE, TRUE, 0);

gui_direct_spin(NULL, &current_light.x[1],
                  -1000.0, 1000.0, 0.1, NULL, NULL, vbox);

/* z component */
vbox = gtk_vbox_new(TRUE, 0);
gtk_box_pack_start(GTK_BOX(hbox), vbox, TRUE, TRUE, 0);
label = gtk_label_new("Z");
gtk_box_pack_start(GTK_BOX(vbox), label, TRUE, TRUE, 0);

gui_direct_spin(NULL, &current_light.x[2],
                  -1000.0, 1000.0, 0.1, NULL, NULL, vbox);


/* next frame - actions */
vbox = gui_frame_vbox(NULL, FALSE, FALSE, vbox1);

hbox = gtk_hbox_new(TRUE, PANEL_SPACING);
gtk_box_pack_start(GTK_BOX(vbox), hbox, TRUE, TRUE, 0);
gui_button(" Add ", add_current_light, NULL, hbox, TT);
gui_button(" Modify ", mod_selected_light, NULL, hbox, TT);
gui_button(" Delete ", del_selected_light, NULL, hbox, TT);

/* right box 1 - light source listing */
vbox = gui_frame_vbox(NULL, TRUE, TRUE, vbox2);

/* scrolled win for the list */
swin = gtk_scrolled_window_new (NULL, NULL);
gtk_scrolled_window_set_policy(GTK_SCROLLED_WINDOW(swin),
                               GTK_POLICY_AUTOMATIC, GTK_POLICY_AUTOMATIC);
gtk_box_pack_start(GTK_BOX(vbox), swin, TRUE, TRUE, 0);
/*
gtk_widget_set_size_request(swin, 24*sysenv.gfontsize, -1);
*/

/* NEW - light list in a list store */
render_light_ls = gtk_list_store_new(RENDER_LIGHT_NCOLS, 
                                     G_TYPE_STRING, G_TYPE_STRING, G_TYPE_STRING, G_TYPE_STRING);
render_light_tv = gtk_tree_view_new_with_model(GTK_TREE_MODEL(render_light_ls));
gtk_scrolled_window_add_with_viewport(GTK_SCROLLED_WINDOW(swin), render_light_tv);
for (i=0 ; i<=RENDER_LIGHT_TYPE ; i++)
  {
  renderer = gtk_cell_renderer_text_new();
  column = gtk_tree_view_column_new_with_attributes(titles[i], renderer, "text", i, NULL);

if (!i) gtk_tree_view_column_add_attribute(column, renderer, "background", RENDER_LIGHT_ATTRIBUTE);

  gtk_tree_view_append_column(GTK_TREE_VIEW(render_light_tv), column);
  }

g_signal_connect(render_light_tv, "row_activated",
                 G_CALLBACK(render_light_activate), NULL);
}

/***************************/
/* OpenGL specific options */
/***************************/
void render_opengl_section(GtkWidget *box)
{
GtkWidget *vbox1, *vbox2, *vbox, *hbox;

/* left & right pane split */
hbox = gtk_hbox_new(FALSE, PANEL_SPACING);
gtk_container_add(GTK_CONTAINER(box), hbox);
gtk_container_set_border_width(GTK_CONTAINER(hbox), PANEL_SPACING);
vbox1 = gtk_vbox_new(FALSE, PANEL_SPACING);
gtk_box_pack_start(GTK_BOX(hbox), vbox1, TRUE, TRUE, 0);
vbox2 = gtk_vbox_new(FALSE, PANEL_SPACING);
gtk_box_pack_start(GTK_BOX(hbox), vbox2, TRUE, TRUE, 0);

/* depth queuing */
vbox = gui_frame_vbox(NULL, FALSE, FALSE, vbox1);

gui_direct_check("Depth queueing", &sysenv.render.fog,
                   render_refresh, NULL, vbox);

gui_direct_spin("Depth queueing strength",
                  &sysenv.render.fog_density, 0.0, 1.0, 0.1,
                  render_refresh, NULL, vbox);

/* new frame */
vbox = gui_frame_vbox(NULL, FALSE, FALSE, vbox1);

/* sphere quality */
gui_direct_spin("Sphere quality",
                  &sysenv.render.sphere_quality, 0.0, 8.0, 1.0,
                  render_refresh, NULL, vbox);

/* cylinder quality */
gui_direct_spin("Cylinder quality",
                  &sysenv.render.cylinder_quality, 4.0, 20.0, 1.0,
                  render_refresh, NULL, vbox);

/* quality */
vbox = gui_frame_vbox(NULL, FALSE, FALSE, vbox1);

gui_direct_check("Automatic quality adjustment",
                   &sysenv.render.auto_quality,
                   NULL, NULL, vbox);
gui_direct_check("Fast rotation", &sysenv.render.fast_rotation,
                   NULL, NULL, vbox);
gui_direct_check("Selection halos", &sysenv.render.halos,
                   render_refresh, NULL, vbox);

/* new frame */
vbox = gui_frame_vbox(NULL, FALSE, FALSE, vbox2);
hbox = gtk_hbox_new(FALSE, 0);
gtk_box_pack_start(GTK_BOX(vbox), hbox, FALSE, FALSE, PANEL_SPACING);

gui_button_x("Label geometric measurements", geom_label_toggle, NULL, hbox);

gui_direct_spin("Line width",
                  &sysenv.render.geom_line_width, 0.5, 9.0, 0.5,
                  render_refresh, NULL, vbox);

/* next frame */
vbox = gui_frame_vbox(NULL, FALSE, FALSE, vbox2);

/* change font options */
gui_button_x("Change graphics font", gl_font_dialog, NULL, vbox);
}

/***************************/
/* POVRay specific options */
/***************************/
void render_povray_section(GtkWidget *box)
{
GtkWidget *vbox1, *vbox2, *vbox, *hbox;

/* left & right pane split */
hbox = gtk_hbox_new(FALSE, PANEL_SPACING);
gtk_container_add(GTK_CONTAINER(box), hbox);
gtk_container_set_border_width(GTK_CONTAINER(hbox), PANEL_SPACING);
vbox1 = gtk_vbox_new(FALSE, PANEL_SPACING);
gtk_box_pack_start(GTK_BOX(hbox), vbox1, TRUE, TRUE, 0);
vbox2 = gtk_vbox_new(FALSE, PANEL_SPACING);
gtk_box_pack_start(GTK_BOX(hbox), vbox2, TRUE, TRUE, 0);

/* new frame */
vbox = gui_frame_vbox("Image size", FALSE, FALSE, vbox1);

gui_direct_spin("Width",
                  &sysenv.render.width, 100, 2000, 100,
                  NULL, NULL, vbox);

gui_direct_spin("Height",
                  &sysenv.render.height, 100, 2000, 100,
                  NULL, NULL, vbox);

/* new frame */
vbox = gui_frame_vbox(NULL, FALSE, FALSE, vbox1);

gui_direct_check("Shadowless", &sysenv.render.shadowless, NULL, NULL, vbox);

gui_direct_check("Create povray input files then stop",
                   &sysenv.render.no_povray_exec,
                   NULL, NULL, vbox);
gui_direct_check("Delete intermediate input/image files",
                   &sysenv.render.no_keep_tempfiles,
                   NULL, NULL, vbox);


/* FIXME - this is broken with the new spatial morphology display */
/*
frame = gtk_frame_new("Glass morphology");
gtk_box_pack_start(GTK_BOX(vbox2), frame, FALSE, FALSE, 0);
vbox = gtk_vbox_new(TRUE, 5);
gtk_container_add(GTK_CONTAINER(frame), vbox);
gtk_container_set_border_width(GTK_CONTAINER(GTK_BOX(vbox)), PANEL_SPACING);

gui_direct_spin("Refractive index",
                  &sysenv.render.ref_index, 1.0, 3.0, 0.1,
                  NULL, NULL, vbox);

hbox = gtk_hbox_new (FALSE, 0);
gtk_box_pack_start(GTK_BOX(vbox),hbox,FALSE,TRUE,0);
label = gtk_label_new("Surface finish ");
gtk_box_pack_start(GTK_BOX(hbox), label, FALSE, TRUE, 0);
entry = gtk_entry_new_with_max_length(LINELEN-1);
gtk_box_pack_end(GTK_BOX(hbox), entry, FALSE, TRUE, 0);

g_signal_connect(GTK_OBJECT(entry), "activate",
                 GTK_SIGNAL_FUNC(event_render_modify), (gpointer) entry);
g_object_set_data(G_OBJECT(entry), "id", (gpointer) MORPH_FINISH);
*/

/* run POVRay job in the background */
hbox = gtk_hbox_new(FALSE, 0);
gtk_box_pack_end(GTK_BOX(box), hbox, TRUE, FALSE, PANEL_SPACING);
gui_button("    Render    ", render_setup, NULL, hbox, TF);
}

/***************************/
/* turn stereo display off */
/***************************/
void gui_stereo_off(GtkWidget *w, gpointer dummy)
{
sysenv.stereo = FALSE;
sysenv.render.perspective = FALSE;
stereo_close_window();
redraw_canvas(SINGLE);
}

/********************************/
/* activate desired stereo mode */
/********************************/
void gui_stereo_on(GtkWidget *w, gpointer mode)
{
switch(GPOINTER_TO_INT(mode))
  {
  case 0:
/* windowed */
    if (sysenv.stereo_fullscreen)
      stereo_close_window();
    sysenv.stereo = TRUE;
    sysenv.stereo_fullscreen = FALSE;
    sysenv.render.perspective = TRUE;
    break;

  case 1:
/* fullscreen */
    sysenv.stereo = TRUE;
    sysenv.stereo_fullscreen = TRUE;
    sysenv.render.perspective = TRUE;
    stereo_open_window();
    break;
  }
redraw_canvas(SINGLE);
}

/************************************/
/* update eye separation percentage */
/************************************/
void gui_stereo_refresh(GtkWidget *w, gpointer dummy)
{
gui_refresh(GUI_CANVAS);
}

/*************************/
/* stereo config options */
/*************************/
void render_stereo_section(GtkWidget *box)
{
GtkWidget *vbox, *label;

g_assert(box != NULL);

vbox = gui_frame_vbox(NULL, FALSE, FALSE, box);

if (sysenv.render.stereo_quadbuffer)
  label = gtk_label_new("Quad buffered OpenGL visual");
else
  label = gtk_label_new("Standard OpenGL visual");

gtk_box_pack_start(GTK_BOX(vbox), label, FALSE, FALSE, 0);

gui_button_x("Windowed stereo ", gui_stereo_on, GINT_TO_POINTER(0), vbox);
gui_button_x("Fullscreen stereo ", gui_stereo_on, GINT_TO_POINTER(1), vbox);
gui_button_x("Stereo off", gui_stereo_off, NULL, vbox);

gui_direct_spin("Eye separation ", &sysenv.render.stereo_eye_offset, 0.01, 4.0, 0.01,
                gui_stereo_refresh, NULL, vbox);

gui_direct_spin("Frustum asymmetry ", &sysenv.render.stereo_parallax, 0.01, 4.0, 0.01,
                gui_stereo_refresh, NULL, vbox);

gui_direct_check("Left eye active ", &sysenv.render.stereo_left,
                 gui_stereo_refresh, NULL, vbox);
gui_direct_check("Right eye active ", &sysenv.render.stereo_right,
                 gui_stereo_refresh, NULL, vbox);
}

/*****************/
/* misc. options */
/*****************/
void render_misc_section(GtkWidget *box)
{
GtkWidget *vbox1, *vbox2, *vbox, *hbox;
GtkWidget *frame;
//more sensitive boxes
GtkWidget *region_render_sensitive_box;
struct model_pak *model;

/* TODO - allow render dialog even when no models, but don't draw */
/* or make insensitive model specific buttons */
model = sysenv.active_model;
/* FIXME - dont draw anything if no models are loaded (needs a better/dynamic fix) */
if (!model)
  return;

/* left & right pane split */
hbox = gtk_hbox_new(FALSE, PANEL_SPACING);
gtk_container_add(GTK_CONTAINER(box), hbox);
gtk_container_set_border_width(GTK_CONTAINER(hbox), PANEL_SPACING);
vbox1 = gtk_vbox_new(FALSE, PANEL_SPACING);
gtk_box_pack_start(GTK_BOX(hbox), vbox1, TRUE, TRUE, 0);
vbox2 = gtk_vbox_new(FALSE, PANEL_SPACING);
gtk_box_pack_start(GTK_BOX(hbox), vbox2, TRUE, TRUE, 0);

/* left pane */
vbox = gui_frame_vbox(NULL, FALSE, FALSE, vbox1);

/* relation check buttons */
gui_auto_check("Show cell", render_refresh, NULL, &model->show_cell, vbox);
gui_auto_check("Show cell images", render_refresh, NULL, &model->show_cell_images, vbox);
gui_auto_check("Show cell lengths", render_refresh, NULL, &model->show_cell_lengths, vbox);
gui_auto_check("Show cores", render_refresh, NULL, &model->show_cores, vbox);
gui_auto_check("Show shells", render_refresh, NULL, &model->show_shells, vbox);
gui_auto_check("Show core-shell links", render_refresh, NULL, &model->show_links, vbox);
gui_auto_check("Show core indices", render_refresh, NULL, &model->show_atom_index, vbox);
gui_auto_check("Show core labels", render_refresh, NULL, &model->show_atom_labels, vbox);
gui_auto_check("Show core types", render_refresh, NULL, &model->show_atom_types, vbox);
gui_auto_check("Show core charges", render_refresh, NULL, &model->show_core_charges, vbox);
gui_auto_check("Show shell charges", render_refresh, NULL, &model->show_shell_charges, vbox);
gui_auto_check("Show atom charges", render_refresh, NULL, &model->show_atom_charges, vbox);

/*VZ*/
gui_auto_check("Show NMR shielding", render_refresh, NULL, &model->show_nmr_shifts, vbox);
gui_auto_check("Show NMR CSA", render_refresh, NULL, &model->show_nmr_csa, vbox);
gui_auto_check("Show NMR EFG", render_refresh, NULL, &model->show_nmr_efg, vbox);

/*
gui_auto_check("Show normal bonds", connect_refresh_global, NULL, &model->build_molecules, vbox);
*/
gui_auto_check("Show hydrogen bonds", connect_refresh_global, NULL, &model->build_hydrogen, vbox);
gui_auto_check("Show zeolite bonds", connect_refresh_global, NULL, &model->build_zeolite, vbox);

gui_auto_check("Show labels on selection ", render_refresh, NULL, &model->show_selection_labels, vbox);

/* next frame */
vbox = gui_frame_vbox(NULL, FALSE, FALSE, vbox2);

gui_auto_check("Show axes", render_refresh, NULL, &model->show_axes, vbox);

gui_direct_check("Show energy", &sysenv.render.show_energy, render_refresh, NULL, vbox);

gui_auto_check("Show spatial text ", morph_toggle, NULL, &model->morph_label, vbox);
gui_auto_check("Show camera waypoints ", morph_toggle, NULL, &model->show_waypoints, vbox);

/* Region options for surfaces */
frame = gtk_frame_new("Regions");
gtk_box_pack_start(GTK_BOX(vbox2), frame, FALSE, FALSE, 0);
region_render_sensitive_box = gtk_vbox_new(TRUE, PANEL_SPACING);
model->custom_regions_frame = region_render_sensitive_box;

gtk_container_add(GTK_CONTAINER(frame), region_render_sensitive_box);

gtk_container_set_border_width(GTK_CONTAINER(region_render_sensitive_box), PANEL_SPACING);

/* gui_auto_check("Show growth slice only", render_refresh, NULL, &model->show_growth_slice, region_render_sensitive_box); */

gui_auto_check("Show region 1", render_refresh, NULL, &model->show_region1A, region_render_sensitive_box);

gui_auto_check("Show region 2", render_refresh, NULL, &model->show_region2A, region_render_sensitive_box);

gui_auto_check("Show region 3", render_refresh, NULL, &model->show_region1B, region_render_sensitive_box);

gui_auto_check("Show region 4", render_refresh, NULL, &model->show_region2B, region_render_sensitive_box);

gtk_widget_set_sensitive(GTK_WIDGET(model->custom_regions_frame), TRUE);
}

/***************************/
/* dialog refresh function */
/***************************/
void render_dialog_refresh(struct model_pak *model)
{
gui_camera_refresh();
}

/***********************/
/* render setup widget */
/***********************/
void gui_render_dialog()
{
gpointer dialog;
GtkWidget *window, *frame;
GtkWidget *label, *notebook, *page;
struct light_pak *ldata;
/*
struct model_pak *model;
*/

/* NB: some check boxes depend on having a valid model */
/*
model = sysenv.active_model;
if (!model)
  return;
*/

/* request dialog slot */
dialog = dialog_request(POVRAY, "Display properties", render_dialog_refresh, render_cleanup, NULL);
if (!dialog)
  return;
window = dialog_window(dialog);

/* init current light */
ldata = ((GSList *) sysenv.render.light_list)->data;

if (ldata)
  memcpy(&current_light, ldata, sizeof(struct light_pak));
else
  {
  gui_text_show(ERROR, "Empty light list.\n");
  return;
  }

/* notebook frame */
frame = gtk_frame_new(NULL);
gtk_box_pack_start(GTK_BOX(GTK_DIALOG(window)->vbox), frame, FALSE, FALSE, 0);
gtk_container_set_border_width(GTK_CONTAINER(frame), PANEL_SPACING);

/* create notebook */
notebook = gtk_notebook_new();
gtk_notebook_set_tab_pos(GTK_NOTEBOOK(notebook), GTK_POS_TOP);
gtk_container_add(GTK_CONTAINER(frame), notebook);
gtk_notebook_set_show_border(GTK_NOTEBOOK(notebook), TRUE);

/* main page */
page = gtk_vbox_new(FALSE, PANEL_SPACING);
label = gtk_label_new("Main");
gtk_notebook_append_page(GTK_NOTEBOOK(notebook), page, label);
render_main_page(page);

/* colour page */
page = gtk_vbox_new(FALSE, PANEL_SPACING);
label = gtk_label_new("Colours");
gtk_notebook_append_page(GTK_NOTEBOOK(notebook), page, label);
render_colours_section(page);

/* NEW */
/* camera page */
page = gtk_vbox_new(FALSE, PANEL_SPACING);
label = gtk_label_new("Camera");
gtk_notebook_append_page(GTK_NOTEBOOK(notebook), page, label);
render_camera_section(page);

/* lighting page */
page = gtk_vbox_new(FALSE, PANEL_SPACING);
label = gtk_label_new("Lights");
gtk_notebook_append_page(GTK_NOTEBOOK(notebook), page, label);
render_lighting_section(page);

/* OpenGL page */
page = gtk_vbox_new(FALSE, PANEL_SPACING);
label = gtk_label_new("OpenGL");
gtk_notebook_append_page(GTK_NOTEBOOK(notebook), page, label);
render_opengl_section(page);

/* POVRay page */
page = gtk_vbox_new(FALSE, PANEL_SPACING);
label = gtk_label_new("POVRay");
gtk_notebook_append_page(GTK_NOTEBOOK(notebook), page, label);
render_povray_section(page);

/* Stereo config page */
page = gtk_vbox_new(FALSE, PANEL_SPACING);
label = gtk_label_new("Stereo");
gtk_notebook_append_page(GTK_NOTEBOOK(notebook), page, label);
render_stereo_section(page);

/* misc page */
page = gtk_vbox_new(FALSE, PANEL_SPACING);
label = gtk_label_new("Toggles");
gtk_notebook_append_page(GTK_NOTEBOOK(notebook), page, label);
render_misc_section(page);

/* terminating button */
gui_stock_button(GTK_STOCK_CLOSE, dialog_destroy, dialog,
                   GTK_DIALOG(window)->action_area);

/* display the dialog */
gtk_widget_show_all(window);

update_light_list();
}
