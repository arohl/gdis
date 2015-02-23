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

#include "gdis.h"
#include "coords.h"
#include "model.h"
#include "edit.h"
#include "file.h"
#include "render.h"
#include "matrix.h"
#include "quaternion.h"
#include "measure.h"
#include "numeric.h"
#include "opengl.h"
#include "gui_shorts.h"
#include "interface.h"
#include "dialog.h"

#define DEBUG 0
#define DEBUG2 0

/* global pak structures */
extern struct sysenv_pak sysenv;
extern struct elem_pak elements[];

GtkWidget *cf_spin;

struct frame_pak
{
gdouble latmat[9];
GSList *core_list;
GSList *shell_list;
};

/******************************/
/* free all pre-loaded frames */
/******************************/
void anim_free_all_frames(struct model_pak *model)
{
GSList *list;
struct frame_pak *frame;

g_assert(model != NULL);

for (list = model->frame_data_list ; list ; list=g_slist_next(list))
  {
  frame = list->data;

  free_slist(frame->core_list);
  free_slist(frame->shell_list);
  g_free(frame);
  }
g_slist_free(model->frame_data_list);
model->frame_data_list = NULL;
}

/**********************************/
/* attempt to pre-load all frames */
/**********************************/
gint anim_read_all_frames(struct model_pak *model)
{
gint i;
struct frame_pak *frame;
FILE *fp;

/* checks */
g_assert(model != NULL);
g_assert(model->animation == TRUE);

/* TODO - test is we have enough memory by using g_try_malloc() */

/* read in all frames and store relevant data */
  fp = fopen(model->filename, "r");
  if (!fp)
    {
    printf("Could not open source file.");
    return(1);
    }
  for (i=0 ; i<model->num_frames ; i++)
    {
/* NB: frames start at 1 */
    read_raw_frame(fp, i+1, model);
    matrix_lattice_init(model);

    frame = g_malloc(sizeof(struct frame_pak *));
    memcpy(frame->latmat, model->latmat, 9*sizeof(gdouble));
    frame->core_list = dup_core_list(model->cores);
    frame->shell_list = dup_shell_list(model->shels);

    model->frame_data_list = g_slist_prepend(model->frame_data_list, frame);
    }
  model->frame_data_list = g_slist_reverse(model->frame_data_list);

return(1);
}

/***********************************/
/* retrieve the nth transformation */
/***********************************/
gint read_transform(gint n, struct model_pak *model)
{
gpointer camera;

g_assert(model != NULL);
g_assert(model->camera != NULL);

camera = g_slist_nth_data(model->transform_list, n);

g_assert(camera != NULL);

model->camera = camera;

return(0);
}

/**************************/
/* retrieve the nth frame */
/**************************/
#define DEBUG_READ_FRAME 0
gint read_frame(FILE *fp, gint n, struct model_pak *model)
{
gint status;
gdouble rmax, v1[3], v2[3];
gpointer camera=NULL;
GString *err_text;

g_assert(model != NULL);

rmax = model->rmax;

/* conventional or transformation style animation */
if (model->transform_list)
  {
/* NEW - process transformation list as an animation */
  status = read_transform(n, model);
  }
else
  {
g_assert(fp != NULL);

ARR3SET(v1, model->centroid);
ARR3SET(v2, model->offset);

status = read_raw_frame(fp, n, model);
if (status)
  {
  err_text = g_string_new("");
  g_string_printf(err_text, "Error reading frame: %d\n", n);
  gui_text_show(ERROR, err_text->str);
  g_string_free(err_text, TRUE);
  return(1);
  }

/* setup (have to save camera) */
camera = camera_dup(model->camera);
model_prep(model);
model->camera = camera;
g_free(model->camera_default);
model->camera_default = camera;

ARR3SET(model->centroid, v1);
ARR3SET(model->offset, v2);
  }

/* apply desired constraint */
if (model->periodic && !model->anim_fix)
  {
  switch (model->anim_confine)
    {
    case PBC_CONFINE_ATOMS:
      coords_confine_cores(model->cores, model);

    case PBC_CONFINE_MOLS:
      coords_compute(model);
      connect_bonds(model);
      connect_molecules(model);
      break;
    }
  }

if (model->anim_noscale)
  model->rmax = rmax;

return(0);
}

/***************************************/
/* current frame modification callback */
/***************************************/
#define DEBUG_SELECT_FRAME 0
void select_frame(GtkWidget *w, struct model_pak *model)
{
FILE *fp=NULL;

g_assert(w != NULL);
g_assert(model != NULL);

#if DEBUG_SELECT_FRAME
printf("request frame: %d, model: %p\n", model->cur_frame, model);
#endif

/* FIXME - better way to cope with transform animations */
if (!model->transform_list)
  fp = fopen(model->filename, "r");

read_frame(fp, model->cur_frame, model);

meas_graft_model(model);

gui_active_refresh();

redraw_canvas(SINGLE);

if (!model->transform_list)
  fclose(fp);
}

/*********************************************************/
/* read and draw the next frame in the current animation */
/*********************************************************/
#define DEBUG_DISPLAY_NEXT_FRAME 0
gint anim_next_frame(struct model_pak *model)
{
gulong time;
gchar *text, *name;

g_assert(model != NULL);

/* increment and test if should we return to the start */
model->cur_frame += model->anim_step;
if (model->cur_frame >= model->num_frames)
  if (model->anim_loop)
    model->cur_frame = 1;

/* continue until we run out of frames (or a stop is flagged) */
if (model->cur_frame < model->num_frames && model->animating)
  {
#if DEBUG_DISPLAY_NEXT_FRAME
printf("displaying [%d]\n", model->cur_frame);
#endif

time = mytimer();

/* if a dialog exists - update via the current frame spinner */
  if (dialog_exists(ANIM, model))
    gui_relation_update(model);
  else
    {
/* otherwise, update manually */
    read_frame(model->afp, model->cur_frame, model);
    meas_graft_model(model);
    gui_active_refresh();
    redraw_canvas(SINGLE);
    }

/* animation adjusted redraw time */
time = mytimer() - time;
model->redraw_cumulative += time;

/* NEW - render to file */
  if (sysenv.render.animate)
    {
    text = g_strdup_printf("%s_%06d.pov", sysenv.render.animate_file, model->cur_frame);
    name = g_build_filename(sysenv.cwd, text, NULL);

    write_povray(name, model);

/* NB: added this as jago keeps locking up on multi-frame renders */
    if (!sysenv.render.no_povray_exec)
      povray_exec(name);

    g_free(text);
    g_free(name);
    }

  return(TRUE);
  }

/* FIXME - find a better way to do this... */
if (!model->transform_list)
  fclose(model->afp);

/* done animation */
model->animating = FALSE;
model->cur_frame--;

/* create movie? */
if (sysenv.render.animate && !sysenv.render.no_povray_exec)
  {
  text = NULL; 
  switch (sysenv.render.animate_type)
    {
    case ANIM_GIF:
      text = g_strdup_printf("%s -delay %d %s_*.tga %s.gif",
                             sysenv.convert_path,
                             (gint) sysenv.render.delay,
                             sysenv.render.animate_file, sysenv.render.animate_file);
      break;

    case ANIM_MPEG:
      text = g_strdup_printf("%s -quality %d -delay %d %s_*.tga %s.mpg",
                             sysenv.convert_path, (gint) sysenv.render.mpeg_quality,
                             (gint) sysenv.render.delay,
                             sysenv.render.animate_file, sysenv.render.animate_file);
      break;
    }

  if (text)
    {
    system(text);
    g_free(text);
    gui_text_show(DEFAULT, "Completed movie creation.\n");
    }
  }

/* cleanup */
if (sysenv.render.no_keep_tempfiles)
  {
#ifndef __WIN32
  text = g_strdup_printf("rm -rf %s_*.pov", sysenv.render.animate_file);
  system(text);
  g_free(text);
  text = g_strdup_printf("rm -rf %s_*.tga", sysenv.render.animate_file);
  system(text);
  g_free(text);
#endif
/* TODO - windows equivalents */
  }

/* done - return FALSE to terminate the timer */
return(FALSE);
}

/**********************************/
/* set current frame to beginning */
/**********************************/
void anim_rewind(GtkWidget *w, gpointer data)
{
struct model_pak *model;

model = dialog_model(data);
g_assert(model != NULL);

model->cur_frame = 0;

gui_relation_update(model);
}

/****************************/
/* set current frame to end */
/****************************/
void anim_fastforward(GtkWidget *w, gpointer data)
{
struct model_pak *model;

model = dialog_model(data);
g_assert(model != NULL);

model->cur_frame = model->num_frames-1;

gui_relation_update(model);
}

/**********************/
/* one frame backward */
/**********************/
void anim_step_backward(GtkWidget *w, gpointer data)
{
struct model_pak *model;

model = dialog_model(data);
g_assert(model != NULL);

if (model->cur_frame > 0)
  {
  model->cur_frame--;
  gui_relation_update(model);
  }
}

/*********************/
/* one frame forward */
/*********************/
void anim_step_forward(GtkWidget *w, gpointer data)
{
struct model_pak *model;

model = dialog_model(data);
g_assert(model != NULL);

if (model->cur_frame < model->num_frames-1)
  {
  model->cur_frame++;
  gui_relation_update(model);
  }
}

/****************************/
/* start an animation timer */
/****************************/
void anim_start(GtkWidget *w, gpointer data)
{
gint freq;
struct model_pak *model;

/* don't allow more than one */
model = dialog_model(data);
g_assert(model != NULL);
if (model->animating)
  return;

/* timeout frequency */
/* NB: we can't refresh any faster than units of 25 (ie the canvas redraw frequency) */
freq = 25 * model->anim_speed;

if (!model->transform_list)
  {
  model->afp = fopen(model->filename, "r");
  if (!model->afp)
    {
    gui_text_show(ERROR, "Failed to open animation stream.\n");
    return;
    }
  }

g_timeout_add(freq, (GSourceFunc) &anim_next_frame, model);
model->animating = TRUE;
}

/****************************/
/* stop the animation timer */
/****************************/
void anim_stop(GtkWidget *w, gpointer data)
{
struct model_pak *model;

model = dialog_model(data);
g_assert(model != NULL);
model->animating = FALSE;
}

/************************************/
/* toggle the render to file option */
/************************************/
void anim_render(GtkWidget *w, GtkWidget *a_frame)
{
if (sysenv.render.animate)
  gtk_widget_set_sensitive(GTK_WIDGET(a_frame), TRUE);
else
  gtk_widget_set_sensitive(GTK_WIDGET(a_frame), FALSE);
}

/***********************************/
/* pbc confinement option handlers */
/***********************************/
void atom_pbc(struct model_pak *model)
{
model->anim_confine = PBC_CONFINE_ATOMS;
}
void mol_pbc(struct model_pak *model)
{
model->anim_confine = PBC_CONFINE_MOLS;
}
void no_pbc(struct model_pak *model)
{
model->anim_confine = PBC_CONFINE_NONE;
}

/********************/
/* cleanup callback */
/********************/
void animate_cleanup(struct model_pak *model)
{
g_assert(model != NULL);
model->animating = FALSE;
}

/***********************/
/* configure animation */
/***********************/
void gui_animate_dialog(void)
{
gchar *tmp, *title;
gpointer dialog;
GtkWidget *window, *frame, *vbox, *hbox, *hbox2, *anim_box;
GtkWidget *notebook, *page, *button, *label, *entry, *scale;
GSList *group=NULL;
struct model_pak *data;

/* checks */
data = sysenv.active_model;
if (!data)
  return;
if (!data->animation)
  return;

/* CURRENT */
/*
anim_read_all_frames(data);
*/

/* prevent recording whilst in animation */
gui_mode_switch(FREE);

/* dialog setup */
title = g_strdup_printf("Animation: %s", data->basename);
dialog = dialog_request(ANIM, title, NULL, animate_cleanup, data);
g_free(title);
if (!dialog)
  return;
window = dialog_window(dialog);

/* notebook frame */
frame = gtk_frame_new(NULL);
gtk_box_pack_start(GTK_BOX(GTK_DIALOG(window)->vbox), frame, FALSE, FALSE, 0);
gtk_container_set_border_width(GTK_CONTAINER(frame), PANEL_SPACING);

/* create notebook */
notebook = gtk_notebook_new();
gtk_notebook_set_tab_pos(GTK_NOTEBOOK(notebook), GTK_POS_TOP);
gtk_container_add(GTK_CONTAINER(frame), notebook);
gtk_notebook_set_show_border(GTK_NOTEBOOK(notebook), TRUE);

/* page 1 */
page = gtk_vbox_new(FALSE, PANEL_SPACING);
label = gtk_label_new("Control");
gtk_notebook_append_page(GTK_NOTEBOOK(notebook), page, label);

/* general animation stuff */
frame = gtk_frame_new(NULL);
gtk_box_pack_start(GTK_BOX(page),frame,FALSE,FALSE,0);
/* create a vbox in the frame */
vbox = gtk_vbox_new(TRUE, PANEL_SPACING);
gtk_container_add(GTK_CONTAINER(frame), vbox);

/* num frames */
hbox = gtk_hbox_new(FALSE, 0);
gtk_container_set_border_width(GTK_CONTAINER(hbox), PANEL_SPACING);
gtk_box_pack_start(GTK_BOX(vbox), hbox, FALSE, FALSE, 0);
label = gtk_label_new("Number of frames:");
gtk_box_pack_start(GTK_BOX(hbox), label, FALSE, FALSE, 0);
tmp = g_strdup_printf("%9d", data->num_frames);
label = gtk_label_new(tmp);
g_free(tmp);
gtk_box_pack_end(GTK_BOX(hbox), label, FALSE, FALSE, 0);

/* desired delay */
gui_direct_spin("Animation delay", &data->anim_speed,
                  1.0, 40.0, 1.0, NULL, NULL, vbox);

gui_direct_spin("Animation step size ", &data->anim_step,
                  1.0, data->num_frames, 1.0, NULL, NULL, vbox);

gui_direct_check("Don't recalculate connectivity", &data->anim_fix, NULL, NULL, vbox);

gui_direct_check("Don't recalculate scale", &data->anim_noscale, NULL, NULL, vbox);

gui_direct_check("Loop", &data->anim_loop, NULL, NULL, vbox);


/* page 2 */
page = gtk_vbox_new(FALSE, PANEL_SPACING);
label = gtk_label_new("Processing");
gtk_notebook_append_page(GTK_NOTEBOOK(notebook), page, label);

/* cycle options */
frame = gtk_frame_new(NULL);
gtk_box_pack_start(GTK_BOX(page),frame,FALSE,FALSE,0);

vbox = gtk_vbox_new(FALSE, 0);
gtk_container_add(GTK_CONTAINER(frame), vbox);
gtk_container_set_border_width(GTK_CONTAINER(vbox), PANEL_SPACING);

/* actions at start of each cycle */
new_radio_group(0, vbox, FF);

button = add_radio_button("Confine atoms to PBC", (gpointer) atom_pbc, data);
if (data->anim_confine == PBC_CONFINE_ATOMS)
  gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(button), TRUE);

button = add_radio_button("Confine mols to PBC", (gpointer) mol_pbc, data);
if (data->anim_confine == PBC_CONFINE_MOLS)
  gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(button), TRUE);

button = add_radio_button("Cell confinement off", (gpointer) no_pbc, data);
if (data->anim_confine == PBC_CONFINE_NONE)
  gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(button), TRUE);


/* page 3 */
page = gtk_vbox_new(FALSE, PANEL_SPACING);
label = gtk_label_new("Rendering");
gtk_notebook_append_page(GTK_NOTEBOOK(notebook), page, label);

/* cycle options */
frame = gtk_frame_new(NULL);
gtk_box_pack_start(GTK_BOX(page), frame, FALSE, FALSE, 0);

vbox = gtk_vbox_new(FALSE, 0);
gtk_container_add(GTK_CONTAINER(frame), vbox);
gtk_container_set_border_width(GTK_CONTAINER(vbox), PANEL_SPACING);


anim_box = gtk_vbox_new(TRUE, 0);
gui_direct_check("Create movie", &sysenv.render.animate,
                   anim_render, anim_box, vbox);
gtk_box_pack_end(GTK_BOX(vbox), anim_box, TRUE, TRUE, 0);


/* sensitive box */
vbox = gtk_vbox_new(FALSE,0);
gtk_box_pack_start(GTK_BOX(anim_box), vbox, TRUE, TRUE, 0);
gtk_container_set_border_width(GTK_CONTAINER(vbox), PANEL_SPACING);


/* start off with a button */
button = gtk_radio_button_new_with_label (NULL, "Animated GIF");
g_signal_connect(GTK_OBJECT(button), "clicked",
                   GTK_SIGNAL_FUNC(event_render_modify), (gpointer) button);
gtk_box_pack_start(GTK_BOX(vbox), button, FALSE, TRUE, 0);
g_object_set_data(G_OBJECT(button), "id", (gpointer) ANIM_GIF);

/* make a radio group */
group = gtk_radio_button_group(GTK_RADIO_BUTTON(button));
/* do the rest of the buttons */
button = gtk_radio_button_new_with_label(group, "MPEG");
g_signal_connect(GTK_OBJECT(button), "clicked",
                 GTK_SIGNAL_FUNC(event_render_modify), (gpointer) button);
gtk_box_pack_start (GTK_BOX (vbox), button, FALSE, TRUE, 0);
g_object_set_data(G_OBJECT(button), "id", (gpointer) ANIM_MPEG);

/* movie creation parameters */
gui_direct_spin("MPEG quality",
                  &sysenv.render.mpeg_quality, 1, 100, 1,
                  NULL, NULL, vbox);
gui_direct_spin("Delay (ms)",
                  &sysenv.render.delay, 0, 100, 5,
                  NULL, NULL, vbox);

/* file entry */
hbox = gtk_hbox_new (FALSE, 0);
gtk_box_pack_start(GTK_BOX(vbox),hbox,FALSE,TRUE,0);
label = gtk_label_new("Filename ");
gtk_box_pack_start (GTK_BOX (hbox), label, FALSE, TRUE, 0);
entry = gtk_entry_new();
gtk_box_pack_end(GTK_BOX (hbox), entry, FALSE, TRUE, 0);
gtk_entry_set_text(GTK_ENTRY(entry), sysenv.render.animate_file);

/* update hook */
g_signal_connect(GTK_OBJECT(entry), "changed",
                 GTK_SIGNAL_FUNC(event_render_modify), (gpointer) entry);
g_object_set_data(G_OBJECT(entry), "id", (gpointer) ANIM_NAME);

/* misc options */
/*
gui_direct_check("Create povray input files then stop",
                   &sysenv.render.no_povray_exec,
                   NULL, NULL, vbox);
gui_direct_check("Delete intermediate input/image files",
                   &sysenv.render.no_keep_tempfiles,
                   NULL, NULL, vbox);
*/


/* set initial state */
if (!sysenv.render.animate)
  gtk_widget_set_sensitive(anim_box, FALSE);


/* NEW - slider for current frame */
frame = gtk_frame_new(NULL);
gtk_box_pack_start(GTK_BOX(GTK_DIALOG(window)->vbox), frame, FALSE, FALSE, 0);

vbox = gtk_vbox_new(FALSE, PANEL_SPACING);
gtk_container_add(GTK_CONTAINER(frame), vbox);

scale = gui_direct_hscale(0, data->num_frames-1, 1, &data->cur_frame, select_frame, data, vbox);
gtk_range_set_update_policy(GTK_RANGE(scale), GTK_UPDATE_DISCONTINUOUS);

/* control buttons */
hbox2 = gtk_hbox_new(FALSE, 0);
gtk_box_pack_start(GTK_BOX(vbox), hbox2, FALSE, FALSE, PANEL_SPACING);
hbox = gtk_hbox_new(TRUE, PANEL_SPACING);
gtk_box_pack_start(GTK_BOX(hbox2), hbox, TRUE, FALSE, 0);

gui_icon_button("GDIS_REWIND", NULL,
                  anim_rewind, dialog,
                  hbox);

gui_icon_button("GDIS_STEP_BACKWARD", NULL,
                  anim_step_backward, dialog,
                  hbox);

gui_icon_button("GDIS_PLAY", NULL,
                  anim_start, dialog,
                  hbox);

gui_icon_button("GDIS_PAUSE", NULL,
                   anim_stop, dialog,
                  hbox);

gui_icon_button("GDIS_STEP_FORWARD", NULL,
                  anim_step_forward, dialog,
                  hbox);

gui_icon_button("GDIS_FASTFORWARD", NULL,
                  anim_fastforward, dialog,
                  hbox);

gui_stock_button(GTK_STOCK_CLOSE, dialog_destroy, dialog,
                   GTK_DIALOG(window)->action_area);

/* display the dialog */
gtk_widget_show_all(window);

gui_relation_update(data);
}

