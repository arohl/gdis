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

#include <fcntl.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <sys/stat.h>

#include "gdis.h"
#include "coords.h"
#include "edit.h"
#include "file.h"
#include "graph.h"
#include "model.h"
#include "parse.h"
#include "scan.h"
#include "task.h"
#include "grid.h"
#include "matrix.h"
#include "surface.h"
#include "spatial.h"
#include "gui_shorts.h"
#include "interface.h"
#include "dialog.h"
#include "opengl.h"
#include "render.h"
#include "gui_image.h"

extern struct sysenv_pak sysenv;
extern struct elem_pak elements[];

/* TODO - avoid using these globals */
GtkWidget *phonon_freq, *phonon_ir, *phonon_raman;

/*****************************/
/* phonon selection handlers */
/*****************************/
void update_phonon_range(struct model_pak *model)
{
/* checks */
g_assert(model != NULL);
if (!model->phonon_slider)
  return;

/* update the slider */
if (model->num_phonons > 0)
  gtk_range_set_range(GTK_RANGE(model->phonon_slider), 1, model->num_phonons);
}

/*****************************/
/* execute a gulp run (task) */
/*****************************/
#define DEBUG_EXEC_GULP_TASK 0
void exec_gulp_task(gpointer ptr, gpointer data)
{
gchar *inpfile;
struct model_pak *model = ptr;
struct task_pak *task = data;

/* checks */
g_assert(model != NULL);
g_assert(task != NULL);

/* construct fullpath input filename - required for writing */
inpfile = g_build_filename(sysenv.cwd, model->gulp.temp_file, NULL);
task->status_file = g_build_filename(sysenv.cwd, model->gulp.out_file, NULL);

#if DEBUG_EXEC_GULP_TASK
printf(" input file: %s\n", inpfile);
printf("output file: %s\n", model->gulp.out_file);
#endif

write_gulp(inpfile, model);
g_free(inpfile);

/* are we supposed to execute GULP? */
if (model->gulp.no_exec)
  {
#if DEBUG_EXEC_GULP_TASK
printf("Skipping GULP execution on user request.\n");
#endif
  return;
  }

exec_gulp(model->gulp.temp_file, model->gulp.out_file);
}

/*****************************/
/* process a gulp run (task) */
/*****************************/
#define DEBUG_PROC_GULP 0
void proc_gulp_task(gpointer ptr)
{
gchar *inp, *res, *out;
GString *line;
struct model_pak *dest, *data;

/* TODO - model locking (moldel_ptr RO/RW etc) to prevent screw ups */
g_return_if_fail(ptr != NULL);
data = ptr;

/* don't attempt to process if gulp execution was turned off */
if (data->gulp.no_exec)
  return;

/* FIXME - if the user has moved directories since submitting */
/* the gulp job then cwd will have changed and this will fail */
inp = g_build_filename(sysenv.cwd, data->gulp.temp_file, NULL);
res = g_build_filename(sysenv.cwd, data->gulp.dump_file, NULL);

out = g_build_filename(sysenv.cwd, data->gulp.out_file, NULL);

switch (data->gulp.run)
  {
  case E_SINGLE:
/* same model (ie with current energetics dialog) so update */
    read_gulp_output(out, data);
/* update energy (TODO - only if successful) */
    line = g_string_new(NULL);
    if (data->gulp.free)
      g_string_sprintf(line, "%f (free energy)", data->gulp.energy);
    else
      g_string_sprintf(line, "%f", data->gulp.energy);

property_add_ranked(2, "Energy", line->str, data);

/* is there a dialog entry to be updated? */
    if (GTK_IS_ENTRY(data->gulp.energy_entry))
      gtk_entry_set_text(GTK_ENTRY(data->gulp.energy_entry), line->str);
    if (data->periodic == 2)
      {
/* update surface dipole dialog entry */
      g_string_sprintf(line,"%f e.Angs",data->gulp.sdipole);
      if (GTK_IS_ENTRY(data->gulp.sdipole_entry))
        gtk_entry_set_text(GTK_ENTRY(data->gulp.sdipole_entry), line->str);

/* update surface energy dialog entry */
      g_string_sprintf(line,"%f    %s",data->gulp.esurf[0], data->gulp.esurf_units);
      if (GTK_IS_ENTRY(data->gulp.esurf_entry))
        gtk_entry_set_text(GTK_ENTRY(data->gulp.esurf_entry), line->str);

/* update attachment energy dialog entry */
      g_string_sprintf(line,"%f    %s",data->gulp.eatt[0], data->gulp.eatt_units);
      if (GTK_IS_ENTRY(data->gulp.eatt_entry))
        gtk_entry_set_text(GTK_ENTRY(data->gulp.eatt_entry), line->str);
      }

    update_phonon_range(data);

    gui_active_refresh();

    g_string_free(line, TRUE);
    break;

  case E_OPTIMIZE:
/* TODO - make it possile to get dialog data by request */
/* so that we can check if a dialog exsits to be updated */
/* get new coords */
/* create new model for the minimized result */
    dest = model_new();
    g_return_if_fail(dest != NULL);

/* read main data from the res file (correct charges etc.) */
    read_gulp(res, dest);
/* graft to the model tree, so subsequent GULP read doesn't replace coords */
    tree_model_add(dest);
/* get the output energy/phonons etc. */
    read_gulp_output(out, dest);


/* FIXME - if the GULP job fails - model_prep() isnt called, and the */
/* model camera isnt initialized => error trap when gdis tries to visualize */
if (!dest->camera)
  {
printf("WARNING: GULP calculation has possibly failed.\n");
  model_prep(dest);
  }


    break;

/* MD */
  default:
    break;
  }

g_free(inp);
g_free(res);
g_free(out);

redraw_canvas(ALL);
return;
}

/***********************************************/
/* controls the use of GULP to optimise coords */
/***********************************************/
#define DEBUG_RUN_GULP 0
void gui_gulp_task(GtkWidget *w, struct model_pak *data)
{
/* checks */
g_return_if_fail(data != NULL);
if (!sysenv.gulp_path)
  {
  gui_text_show(ERROR, "GULP executable was not found.\n");
  return;
  }

#if DEBUG_RUN_GULP
printf("output file: %s\n", data->gulp.out_file);
#endif

task_new("gulp", &exec_gulp_task, data, &proc_gulp_task, data, data);
}

/*************************/
/* GULP run type toggles */
/*************************/
void gulp_run_single(struct model_pak *data)
{
data->gulp.run = E_SINGLE;
}
void gulp_run_optimize(struct model_pak *data)
{
data->gulp.run = E_OPTIMIZE;
}
void gulp_run_dynamics(struct model_pak *data)
{
data->gulp.run = MD;
}


/************************/
/* GULP keyword toggles */
/************************/
gint cb_gulp_keyword(GtkWidget *w, gpointer *obj)
{
gint keyword;
struct model_pak *data;

keyword = GPOINTER_TO_INT(g_object_get_data(G_OBJECT(obj), "key"));

data = (struct model_pak *) g_object_get_data(G_OBJECT(obj), "ptr");
g_return_val_if_fail(data != NULL, FALSE);

switch(keyword)
  {
  case NVE:
  case NVT:
  case NPT:
    if (data->gulp.run != MD)
      {
      gui_text_show(WARNING, "Ignoring ensemble as it's not a dynamics run.\n");
      return(TRUE);
      }
    data->gulp.ensemble = keyword;
    break;
  case CONP:
  case CONV:
    data->gulp.method = keyword;
    break;
  case BFGS_OPT:
  case CONJ_OPT:
  case RFO_OPT:
    data->gulp.optimiser = keyword;
    break;
  case MOLE:
  case MOLMEC:
  case MOLQ:
  case NOBUILD:
    data->gulp.coulomb = keyword;
    break;
  case TEMP_FILE:
    g_free(data->gulp.temp_file);
    data->gulp.temp_file = g_strdup(gtk_entry_get_text(GTK_ENTRY(obj)));
    parse_space_replace(data->gulp.temp_file,'_');

    break;
  case DUMP_FILE:
    g_free(data->gulp.dump_file);
    data->gulp.dump_file = g_strdup(gtk_entry_get_text(GTK_ENTRY(obj)));
    parse_space_replace(data->gulp.dump_file,'_');
    break;
  case SBULK_ENERGY:
    data->gulp.sbulkenergy = str_to_float(gtk_entry_get_text(GTK_ENTRY(obj)));
    break;
  }

return(FALSE);
}

/*******************************************************/
/* alternate event parser (overlap in optimiser stuff) */
/*******************************************************/
gint switch_keyword(GtkWidget *w, gpointer *obj)
{
gint keyword;
struct model_pak *data;

keyword = GPOINTER_TO_INT(g_object_get_data(G_OBJECT(obj), "key"));

data = (struct model_pak *) g_object_get_data(G_OBJECT(obj), "ptr");
g_return_val_if_fail(data != NULL, FALSE);

switch(keyword)
  {
  case SWITCH_OFF:
  case BFGS_OPT:
  case CONJ_OPT:
  case RFO_OPT:
    data->gulp.optimiser2 = keyword;
    break;
  case CYCLE:
  case GNORM:
    data->gulp.switch_type = keyword;
    break;
  }
return(FALSE);
}

/*********************************/
/* change optimiser switch value */
/*********************************/
gint switch_value_changed(gpointer *entry)
{
struct model_pak *data;

/* retrieve the appropriate model */
g_return_val_if_fail(entry != NULL, FALSE);
data = (struct model_pak *) g_object_get_data(G_OBJECT(entry), "ptr");
g_return_val_if_fail(data != NULL, FALSE);

data->gulp.switch_value = str_to_float(gtk_entry_get_text(GTK_ENTRY(entry)));
return(FALSE);
}

/***********************************/
/* register structure name changes */
/***********************************/
void gulp_jobname_changed(GtkWidget *w, struct model_pak *model)
{
const gchar *text;

g_assert(w != NULL);
g_assert(model != NULL);

g_free(model->gulp.dump_file);
g_free(model->gulp.temp_file);
g_free(model->gulp.trj_file);
g_free(model->gulp.out_file);

text = gtk_entry_get_text(GTK_ENTRY(w));
model->gulp.temp_file = g_strdup_printf("%s.gin", text);
model->gulp.dump_file = g_strdup_printf("%s.res", text);
model->gulp.trj_file = g_strdup_printf("%s.trg", text);
model->gulp.out_file = g_strdup_printf("%s.got", text);

/* disallow spaces (bad filename) */
parse_space_replace(model->gulp.temp_file, '_');
parse_space_replace(model->gulp.dump_file, '_');
parse_space_replace(model->gulp.trj_file, '_');
parse_space_replace(model->gulp.out_file, '_');

gui_relation_update_widget(&model->gulp.dump_file);
gui_relation_update_widget(&model->gulp.temp_file);
gui_relation_update_widget(&model->gulp.trj_file);
}

/**************************/
/* phonon scaling handler */
/**************************/
void phonon_scaling_changed(GtkWidget *w, gpointer *ptr)
{
sysenv.render.phonon_scaling = GTK_ADJUSTMENT(w)->value;
}

#define DEBUG_PHONON_DISPLAY 0

/************************************/
/* vibrational mode display handler */
/************************************/
void phonon_mode_hide(GtkWidget *w, struct model_pak *model)
{
spatial_destroy_by_label("phonons", model);
redraw_canvas(SINGLE);
}

/************************************/
/* vibrational mode display handler */
/************************************/
void phonon_mode_show(GtkWidget *w, struct model_pak *data)
{
gdouble x1[3], x2[3], colour[3];
GSList *list, *xlist, *ylist, *zlist;
gpointer *freq;
struct spatial_pak *spatial;
struct core_pak *core, *prim;

/* get rid of previous phonon objects */
spatial_destroy_by_label("phonons", data);
VEC3SET(colour, 0.0, 1.0, 0.0);

/* find and display current */
freq = g_slist_nth_data(data->phonons, data->current_phonon-1);
if (freq)
  {
/* create & init the spatial object */
  spatial = spatial_new("phonons", SPATIAL_VECTOR, 2, TRUE, data);

/* get eigenvectors from all atoms */
  for (list=data->cores ; list; list=g_slist_next(list))
    {
    core = list->data;
    ARR3SET(x1, core->x);

/* get eigenvector list */
    if (core->primary)
      {
      xlist = core->vibx_list; 
      ylist = core->viby_list; 
      zlist = core->vibz_list; 
      }
    else
      {
      prim = core->primary_core;
      xlist = prim->vibx_list; 
      ylist = prim->viby_list; 
      zlist = prim->vibz_list; 
      }
    g_assert(xlist != NULL);
    g_assert(ylist != NULL);
    g_assert(zlist != NULL);

/* get current eigenvector */
    x2[0] = *((gdouble *) g_slist_nth_data(xlist, data->current_phonon-1));
    x2[1] = *((gdouble *) g_slist_nth_data(ylist, data->current_phonon-1));
    x2[2] = *((gdouble *) g_slist_nth_data(zlist, data->current_phonon-1));
 
#if DEBUG_PHONON_DISPLAY
P3VEC("vec: ", x2);
#endif

/* compute coords */
    VEC3MUL(x2, sysenv.render.phonon_scaling);
    vecmat(data->ilatmat, x2);
    ARR3ADD(x2, x1);
/* add to spatial */
    spatial_vertex_add(x2, colour, spatial);
    spatial_vertex_add(x1, colour, spatial);
    }

/* drawing update */
  coords_compute(data);
  redraw_canvas(SINGLE);
  }
}

/********************************/
/* eigenvector show/hide toggle */
/********************************/
void phonon_mode_toggle(GtkWidget *w, struct model_pak *model)
{
if (model->show_eigenvectors)
  phonon_mode_show(w, model);
else
  phonon_mode_hide(w, model);
}

/***********************************/
/* cleanup a vibrational animation */
/***********************************/
void phonon_timeout_cleanup(struct model_pak *model)
{
GSList *list;
struct core_pak *core;

g_assert(model != NULL);

for (list=model->cores ; list; list=g_slist_next(list))
  {
  core = list->data;

  VEC3SET(core->offset, 0.0, 0.0, 0.0);
  }
}

/************************************/
/* timeout to control the animation */
/************************************/
/* NB: lots of sanity checks that return with FALSE (ie stop timeout) */
/* so if the model is deleted during an animation we don't segfault */
#define MAX_PULSE_COUNT 10.0
gint phonon_mode_timeout(struct model_pak *model)
{
static gint count=0;
gdouble f;
gchar *name, *text;
gpointer ptr;
GSList *list, *xlist, *ylist, *zlist;
struct core_pak *core, *prim;

/* checks */
g_assert(model != NULL);
ptr = g_slist_nth_data(model->phonons, model->current_phonon-1);
if (!ptr)
  {
  phonon_timeout_cleanup(model);
  return(FALSE);
  }
if (!model->pulse_direction)
  return(FALSE);

/* setup scaling for this step */
model->pulse_count += model->pulse_direction;
if (model->pulse_count <= -model->pulse_max)
  {
  model->pulse_count = -model->pulse_max;
  model->pulse_direction = 1;
  }
if (model->pulse_count >= model->pulse_max)
  {
  model->pulse_count = model->pulse_max;
  model->pulse_direction = -1;
  }

f = sysenv.render.phonon_scaling;
f *= (gdouble) model->pulse_count;
f /= model->pulse_max;

/* get eigenvectors from all atoms */
for (list=model->cores ; list; list=g_slist_next(list))
  {
  core = list->data;

/* get eigenvector list */
  if (core->primary)
    {
    xlist = core->vibx_list; 
    ylist = core->viby_list; 
    zlist = core->vibz_list; 
    }
  else
    {
    prim = core->primary_core;
    xlist = prim->vibx_list; 
    ylist = prim->viby_list; 
    zlist = prim->vibz_list; 
    }
  g_assert(xlist != NULL);
  g_assert(ylist != NULL);
  g_assert(zlist != NULL);

/* vibrational eigenvector */
  ptr = g_slist_nth_data(xlist, model->current_phonon-1);
  if (!ptr)
    {
    phonon_timeout_cleanup(model);
    return(FALSE);
    }
  core->offset[0] = *((gdouble *) ptr);

  ptr = g_slist_nth_data(ylist, model->current_phonon-1);
  if (!ptr)
    {
    phonon_timeout_cleanup(model);
    return(FALSE);
    }
  core->offset[1] = *((gdouble *) ptr);

  ptr = g_slist_nth_data(zlist, model->current_phonon-1);
  if (!ptr)
    {
    phonon_timeout_cleanup(model);
    return(FALSE);
    }
  core->offset[2] = *((gdouble *) ptr);

/* pulse offset scaling */
  VEC3MUL(core->offset, f);
  vecmat(model->ilatmat, core->offset);

/* TODO - shell also? */
  }

/* recalc coords */
coords_compute(model);

/* CURRENT - output to povray for movie rendering */
if (model->phonon_movie)
  {

  if (!model->pulse_count && model->pulse_direction==1)
    {
    model->phonon_movie = FALSE;
    count=0;

    text = g_strdup_printf("%s -delay 20 %s_*.tga %s.gif", sysenv.convert_path,
                           model->phonon_movie_name, model->phonon_movie_name);

    system(text);
    g_free(text);

    return(FALSE);
    }
  else
    {
    text = g_strdup_printf("%s_%06d.pov", model->phonon_movie_name, count++);
    name = g_build_filename(sysenv.cwd, text, NULL);
    write_povray(name, model);

    povray_exec(name);

    g_free(text);
    g_free(name);
    }

  }
else
  redraw_canvas(SINGLE);

return(TRUE);
}

/****************************/
/* animate the current mode */
/****************************/
void phonon_mode_animate(GtkWidget *w, struct model_pak *model)
{
g_assert(model != NULL);

model->pulse_count = 0;
model->pulse_direction = 1;

g_timeout_add(100, (GSourceFunc) phonon_mode_timeout, model);
}

/****************************/
/* create a movie of the current mode */
/****************************/
void phonon_mode_movie(GtkWidget *w, struct model_pak *model)
{
g_assert(model != NULL);

model->phonon_movie = TRUE;
phonon_mode_animate(NULL, model);
}

/******************************/
/* stop phonon mode animation */
/*****************************/
void phonon_mode_stop(GtkWidget *w, struct model_pak *model)
{
/* reset */
model->pulse_direction = 0;
model->pulse_count = 0;
phonon_timeout_cleanup(model);

/* redraw */
coords_compute(model);
redraw_canvas(SINGLE);
}

/*********************************/
/* phonon slider changed updates */
/*********************************/
void update_phonon_frequency(GtkWidget *w, struct model_pak *model)
{
gpointer freq, ir, raman;

g_assert(model != NULL);

freq = g_slist_nth_data(model->phonons, model->current_phonon-1);
if (freq)
  {
  ir = g_slist_nth_data(model->ir_list, model->current_phonon-1);
  raman = g_slist_nth_data(model->raman_list, model->current_phonon-1);

/* update displayed frequency */
  gtk_entry_set_text(GTK_ENTRY(phonon_freq), (gchar *) freq);
  gtk_entry_set_text(GTK_ENTRY(phonon_ir), (gchar *) ir);
  gtk_entry_set_text(GTK_ENTRY(phonon_raman), (gchar *) raman);


/* update displayed vectors */
  if (model->show_eigenvectors)
    {
    phonon_mode_hide(NULL, model);
    phonon_mode_show(NULL, model);
    }
  }
}

/*****************************/
/* phonon intensity plotting */
/*****************************/
void phonon_graph(GSList *mode, struct model_pak *model)
{
gint n;
gdouble *spectrum;
gpointer graph;
GSList *list;

g_assert(model != NULL);

/* checks */
if (model->num_phonons < 1)
  return;

spectrum = g_malloc(model->num_phonons * sizeof(gdouble));

/* extract data into array */
n=0;
for (list=mode ; list ; list=g_slist_next(list))
  {
g_assert(n < model->num_phonons);

  spectrum[n] = str_to_float(list->data);
  n++;
  }

/* create a new graph */
if (mode == model->ir_list)
  graph = graph_new("IR spectrum", model);
else
  graph = graph_new("Raman spectrum", model);

graph_add_data(model->num_phonons, spectrum, 1, model->num_phonons, graph);
graph_set_yticks(FALSE, 2, graph);
tree_model_add(model);

g_free(spectrum);
}

/***************************************/
/* phonon intensity plotting callbacks */
/***************************************/
void cb_graph_ir(GtkWidget *w, struct model_pak *model)
{
g_assert(model != NULL);
phonon_graph(model->ir_list, model);
}

void cb_graph_raman(GtkWidget *w, struct model_pak *model)
{
g_assert(model != NULL);
phonon_graph(model->raman_list, model);
}


/* CURRENT - invoke test module */
/*
void phonon_module_invoke(void)
{
struct model_pak *model;

model = sysenv.active_model;

module_function_invoke("test", model);
}
*/

/*******************************/
/* decrement current frequency */
/*******************************/
void phonon_mode_prev(GtkWidget *w, struct model_pak *model)
{
if (model->current_phonon > 1)
  model->current_phonon--;

gtk_range_set_value(GTK_RANGE(model->phonon_slider),
                   (gdouble) model->current_phonon);
}

/*******************************/
/* increment current frequency */
/*******************************/
void phonon_mode_next(GtkWidget *w, struct model_pak *model)
{
if (model->current_phonon < model->num_phonons)
  model->current_phonon++;

gtk_range_set_value(GTK_RANGE(model->phonon_slider),
                   (gdouble) model->current_phonon);
}

/****************/
/* phonon page  */
/****************/
void gulp_phonon_box(GtkWidget *box, struct model_pak *data)
{
GtkWidget *swin, *table, *frame, *hbox, *hbox2, *vbox, *label;

hbox = gtk_hbox_new(FALSE, 0);
gtk_box_pack_start(GTK_BOX(box), hbox, TRUE, TRUE, 0);
gtk_container_set_border_width(GTK_CONTAINER(hbox), PANEL_SPACING);

frame = gtk_frame_new(NULL);
gtk_box_pack_start(GTK_BOX(hbox), frame, FALSE, FALSE, 0);
vbox = gtk_vbox_new(FALSE, 0);
gtk_container_add(GTK_CONTAINER(frame), vbox);
hbox2 = gtk_hbox_new(FALSE, 0);
gtk_box_pack_start(GTK_BOX(vbox), hbox2, FALSE, FALSE, 0);
gui_direct_check("Compute vibrational modes", &data->gulp.phonon, NULL, NULL, hbox2);
hbox2 = gtk_hbox_new(FALSE, 0);
gtk_box_pack_start(GTK_BOX(vbox), hbox2, FALSE, FALSE, 0);
gui_direct_check("Compute eigenvectors", &data->gulp.eigen, NULL, NULL, hbox2);

if (data->periodic)
  {
  frame = gtk_frame_new("kpoints");
  gtk_box_pack_end(GTK_BOX(hbox), frame, TRUE, TRUE, PANEL_SPACING);
  swin = gui_text_window(&data->gulp.kpoints, TRUE);
  gtk_container_add(GTK_CONTAINER(frame), swin);
  }

/* phonon slide selector */
frame = gtk_frame_new("Eigenvectors");
gtk_box_pack_start(GTK_BOX(box), frame, FALSE, FALSE, 0);
gtk_container_set_border_width(GTK_CONTAINER(frame), PANEL_SPACING);
vbox = gtk_vbox_new(FALSE, PANEL_SPACING);
gtk_container_add(GTK_CONTAINER(frame), vbox);
hbox = gtk_hbox_new(FALSE, 0);
gtk_box_pack_start(GTK_BOX(vbox), hbox, FALSE, FALSE, 0);
label = gtk_label_new("Number  ");
gtk_box_pack_start(GTK_BOX(hbox), label, FALSE, FALSE, 0);
data->phonon_slider = gui_direct_hscale(0, 1, 1, &data->current_phonon, 
                                    update_phonon_frequency, data, hbox);

/* CURRENT */
gui_button(" < ", phonon_mode_prev, data, hbox, FF);
gui_button(" > ", phonon_mode_next, data, hbox, FF);


/* phonon info */
table = gtk_table_new(2, 3, TRUE);
gtk_box_pack_start(GTK_BOX(vbox), table, TRUE, TRUE, 0);

label = gtk_label_new("Frequency");
gtk_table_attach_defaults(GTK_TABLE(table), label, 0, 1, 0, 1);

label = gui_button("IR intensity", cb_graph_ir, data, NULL, 0);
gtk_table_attach_defaults(GTK_TABLE(table), label, 1, 2, 0, 1);

label = gui_button("Raman intensity", cb_graph_raman, data, NULL, 0);
gtk_table_attach_defaults(GTK_TABLE(table), label, 2, 3, 0, 1);

phonon_freq = gtk_entry_new_with_max_length(LINELEN);
gtk_entry_set_text(GTK_ENTRY(phonon_freq), " ");
gtk_entry_set_editable(GTK_ENTRY(phonon_freq), FALSE);
gtk_table_attach_defaults(GTK_TABLE(table), phonon_freq, 0, 1, 1, 2);

phonon_ir = gtk_entry_new_with_max_length(LINELEN);
gtk_entry_set_text(GTK_ENTRY(phonon_ir), " ");
gtk_entry_set_editable(GTK_ENTRY(phonon_ir), FALSE);
gtk_table_attach_defaults(GTK_TABLE(table), phonon_ir, 1, 2, 1, 2);

phonon_raman = gtk_entry_new_with_max_length(LINELEN);
gtk_entry_set_text(GTK_ENTRY(phonon_raman), " ");
gtk_entry_set_editable(GTK_ENTRY(phonon_raman), FALSE);
gtk_table_attach_defaults(GTK_TABLE(table), phonon_raman, 2, 3, 1, 2);


gui_direct_check("Display Eigenvectors", &data->show_eigenvectors, 
                 phonon_mode_toggle, data, vbox);

gui_direct_spin("Eigenvector scaling", &sysenv.render.phonon_scaling,
                0.1, 9.9, 0.1, update_phonon_frequency, data, vbox);

gui_direct_spin("Animation resolution", &data->pulse_max,
                10.0, 100.0, 1.0, NULL, NULL, vbox);

/* animation */
hbox2 = gtk_hbox_new(FALSE, 0);
gtk_box_pack_start(GTK_BOX(vbox), hbox2, TRUE, FALSE, 0);

gui_text_entry(NULL, &data->phonon_movie_name, TRUE, FALSE, hbox2);


hbox = gtk_hbox_new(TRUE, PANEL_SPACING);
gtk_box_pack_end(GTK_BOX(hbox2), hbox, FALSE, FALSE, 0);

gui_icon_button("image_camera", NULL, phonon_mode_movie, (gpointer) data, hbox);

gui_icon_button("GDIS_PLAY", NULL,
                  phonon_mode_animate, (gpointer) data, hbox);

gui_icon_button("GDIS_STOP", NULL, phonon_mode_stop, (gpointer) data, hbox);


gtk_widget_show_all(box);

update_phonon_range(data);
}

/*********************************************/
/* read in solvent accessible surface points */
/*********************************************/
/* TODO - put elsewhere? */
void gulp_cosmic_read(gchar *filename, struct model_pak *model)
{
gint i, n=0, dim, expect=-1, num_tokens;
gchar **buff;
gdouble x[3], colour[3];
gpointer scan, spatial;

g_assert(filename != NULL);
g_assert(model != NULL);

VEC3SET(colour, 0.0, 1.0, 0.0);

n=0;
scan = scan_new(filename);
if (!scan)
  return;

/* header */
buff = scan_get_tokens(scan, &num_tokens);
if (num_tokens == 1)
  expect = str_to_float(*buff);

/* new GULP version  - dimension + cell vectors */
buff = scan_get_tokens(scan, &num_tokens);
if (num_tokens == 1)
  {
  dim = str_to_float(*buff);
  g_strfreev(buff);

/* TODO - process cell vectors */
  for (i=dim ; i-- ; )
    {
    buff = scan_get_tokens(scan, &num_tokens);
    g_strfreev(buff);
    }
  }
else
  scan_put_line(scan);

/* body */
while (!scan_complete(scan))
  {
  buff = scan_get_tokens(scan, &num_tokens);

  if (num_tokens == 6)
    {
/* coordinates */
    x[0] = str_to_float(*(buff+2));
    x[1] = str_to_float(*(buff+3));
    x[2] = str_to_float(*(buff+4));
    vecmat(model->ilatmat, x);
/* add a new point */
    spatial = spatial_new("SAS", SPATIAL_GENERIC, 1, TRUE, model);
    spatial_vertex_add(x, colour, spatial);
    n++;
    }
  else
    break;

  g_strfreev(buff);
  }

/*
printf("points: %d / %d\n", n, expect);
*/

scan_free(scan);
}

void gulp_cosmic_display(GtkWidget *w, struct model_pak *model)
{
gchar *name;

name = parse_extension_set(model->filename, "sas");

/*
printf("cosmic: %s\n", name);
*/

gulp_cosmic_read(name, model);

/* updates */
coords_init(CENT_COORDS, model);

redraw_canvas(SINGLE);

g_free(name);
}

/*************************************/
/* set COSMO widget sensitive status */
/*************************************/
void gulp_cosmo_panel_refresh(GtkWidget *w, gpointer dialog)
{
gchar *text;
GtkWidget *box;
GtkEntry *solvation;
struct model_pak *model;

box = dialog_child_get(dialog, "cosmo_box");
g_assert(box != NULL);

model = dialog_model(dialog);
g_assert(model != NULL);

solvation = dialog_child_get(dialog, "solvation_model");
text = gui_pd_text(solvation);
if (!text)
  return;

if (g_ascii_strncasecmp(text, "none", 4) == 0)
  {
  gtk_widget_set_sensitive(box, FALSE); 
  model->gulp.solvation_model = GULP_SOLVATION_NONE;
  }

if (g_ascii_strncasecmp(text, "cosmo", 5) == 0)
  {
  gtk_widget_set_sensitive(box, TRUE); 
  model->gulp.solvation_model = GULP_SOLVATION_COSMO;
  }

if (g_ascii_strncasecmp(text, "cosmic", 6) == 0)
  {
  gtk_widget_set_sensitive(box, TRUE); 
  model->gulp.solvation_model = GULP_SOLVATION_COSMIC;
  }

g_free(text);
}

/********************************************/
/* calculate points for a COSMO calculation */
/********************************************/
gint gulp_cosmo_points(gint shape, gint k, gint l)
{
gint points;

if (shape)
  points = 2 + 10 * pow(3, k) * pow(4, l);
else
  points = 2 + 4 * pow(3, k) * pow(4, l);

return(points);
}

/********************************************/
/* update points for a COSMO calculation */
/********************************************/
void gulp_cosmo_points_refresh(GtkWidget *w, gpointer dialog)
{
gint k, l;
gchar *text;
GtkLabel *label;
GtkSpinButton *spin;
GtkEntry *shape;
struct model_pak *model;

model = dialog_model(dialog);
g_assert(model != NULL);

shape = dialog_child_get(dialog, "cosmo_shape");
if (g_ascii_strncasecmp(gtk_entry_get_text(shape), "dodec", 5) == 0)
  model->gulp.cosmo_shape = 1;
else
  model->gulp.cosmo_shape = 0;

spin = dialog_child_get(dialog, "cosmo_index_k");
k = SPIN_IVAL(spin);

spin = dialog_child_get(dialog, "cosmo_index_l");
l = SPIN_IVAL(spin);

model->gulp.cosmo_points = gulp_cosmo_points(model->gulp.cosmo_shape, k, l);

label = dialog_child_get(dialog, "cosmo_points");
g_assert(label != NULL);
text = g_strdup_printf("%6d", model->gulp.cosmo_points);
gtk_label_set_text(label, text);
g_free(text);

/* CURRENT - enforce segments <= points */
spin = dialog_child_get(dialog, "cosmo_segments");
g_assert(spin != NULL);

gtk_spin_button_set_range(spin, 1.0, (gdouble) model->gulp.cosmo_points);
}

/*******************************************/
/* update segments for a COSMO calculation */
/*******************************************/
void gulp_cosmo_segments_refresh(GtkSpinButton *w, gpointer dialog)
{
struct model_pak *model;

model = dialog_model(dialog);
g_assert(model != NULL);

model->gulp.cosmo_segments = SPIN_IVAL(w);
}

/**************************************************/
/* toggle solvent accessible surface (SAS) output */
/**************************************************/
void gulp_cosmo_sas_toggle(GtkWidget *w, gpointer data)
{
struct model_pak *model = data;

g_assert(model != NULL);

if (gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(w)))
  model->gulp.cosmo_sas = TRUE;
else
  model->gulp.cosmo_sas = FALSE;
}

/*********************************************************/
/* options for displaying the solvent accessible surface */
/*********************************************************/
void gulp_solvation_box(GtkWidget *box, gpointer dialog)
{
gint active;
GSList *slist;
GList *list=NULL;
gpointer cosmo, ptr;
GtkWidget *hbox, *vbox, *vbox2, *label, *shape_k, *shape_l, *segments;
struct model_pak *model;

model = dialog_model(dialog);
g_assert(model != NULL);

hbox = gtk_hbox_new(FALSE, 0);
gtk_box_pack_start(GTK_BOX(box), hbox, FALSE, FALSE, 0);
gtk_container_set_border_width(GTK_CONTAINER(hbox), PANEL_SPACING);

vbox = gtk_vbox_new(FALSE, PANEL_SPACING);
gtk_box_pack_start(GTK_BOX(hbox), vbox, FALSE, FALSE, 0);
gtk_container_set_border_width(GTK_CONTAINER(vbox), PANEL_SPACING);

vbox2 = gtk_vbox_new(FALSE, PANEL_SPACING);
gtk_box_pack_end(GTK_BOX(vbox), vbox2, TRUE, TRUE, 0);
dialog_child_set(dialog, "cosmo_box", vbox2);

/*
cosmo = new_check_button("COSMO calculation ", gulp_cosmo_panel_refresh, dialog,
                         model->gulp.cosmo, vbox);
*/

switch(model->gulp.solvation_model)
  {
  case GULP_SOLVATION_COSMIC:
    active = 1;
    break;
  case GULP_SOLVATION_COSMO:
    active = 2;
    break;
  default:
    active = 0;
  }

/* FIXME - pulldown needs a way to initialize */
slist = NULL;
slist = g_slist_append(slist, "None");
slist = g_slist_append(slist, "COSMIC");
slist = g_slist_append(slist, "COSMO");
cosmo = gui_pd_new(slist, active, gulp_cosmo_panel_refresh, dialog);

gui_hbox_pack(vbox, "Solvation model ", cosmo, 0);

dialog_child_set(dialog, "solvation_model", cosmo);

g_slist_free(slist);


/* solvent parameters */

gui_direct_spin("Solvent epsilon ", &model->gulp.cosmo_solvent_epsilon,
                                   1.0, 1000.0, 0.1, NULL, NULL, vbox2);
gui_direct_spin("Solvent radius ", &model->gulp.cosmo_solvent_radius,
                                   0.1, 9.9, 0.1, NULL, NULL, vbox2);
gui_direct_spin("Solvent delta ", &model->gulp.cosmo_solvent_delta,
                                   0.1, 9.9, 0.1, NULL, NULL, vbox2);
gui_direct_spin("Solvent rmax ", &model->gulp.cosmo_solvent_rmax,
                                 1.0, 99.0, 1.0, NULL, NULL, vbox2);
gui_direct_spin("Smoothing ", &model->gulp.cosmo_smoothing,
                              0.0, 2.0, 0.1, NULL, NULL, vbox2);

/* solvent surface construction geometry */
list = NULL;
if (model->gulp.cosmo_shape)
  {
  list = g_list_append(list, "Dodecahedron");
  list = g_list_append(list, "Octahedron");
  }
else
  {
  list = g_list_append(list, "Octahedron");
  list = g_list_append(list, "Dodecahedron");
  }
ptr = gui_pulldown_new("Shape approximation  ", list, FALSE, vbox2);
g_signal_connect(GTK_OBJECT(ptr), "changed",
                 GTK_SIGNAL_FUNC(gulp_cosmo_points_refresh), dialog);

dialog_child_set(dialog, "cosmo_shape", ptr);
g_list_free(list);

hbox = gtk_hbox_new(FALSE, 0);
gtk_box_pack_start(GTK_BOX(vbox2), hbox, FALSE, FALSE, 0);

shape_k = new_spinner("Shape indices ", 0, 99, 1, gulp_cosmo_points_refresh, dialog, hbox);
shape_l = new_spinner(NULL, 0, 99, 1, gulp_cosmo_points_refresh, dialog, hbox);
dialog_child_set(dialog, "cosmo_index_k", shape_k);
dialog_child_set(dialog, "cosmo_index_l", shape_l);

hbox = gtk_hbox_new(FALSE, 0);
gtk_box_pack_start(GTK_BOX(vbox2), hbox, FALSE, FALSE, 0);
label = gtk_label_new("Points per atom ");
gtk_box_pack_start(GTK_BOX(hbox), label, FALSE, FALSE, 0);
label = gtk_label_new("?");
gtk_box_pack_end(GTK_BOX(hbox), label, FALSE, FALSE, 0);
dialog_child_set(dialog, "cosmo_points", label);

/* CURRENT */
hbox = gtk_hbox_new(FALSE, 0);
gtk_box_pack_start(GTK_BOX(vbox2), hbox, FALSE, FALSE, 0);
segments = new_spinner("Segments per atom ", 1, model->gulp.cosmo_points, 1,
                       gulp_cosmo_segments_refresh, dialog, hbox);
dialog_child_set(dialog, "cosmo_segments", segments);


/* solvent surface */
new_check_button("Output SAS points ", gulp_cosmo_sas_toggle, model,
                         model->gulp.cosmo_sas, vbox2);

gui_button_x("Visualize SAS points ", gulp_cosmic_display, model, vbox2);



/* initialize values */
gtk_spin_button_set_value(GTK_SPIN_BUTTON(shape_k), model->gulp.cosmo_shape_index[0]);
gtk_spin_button_set_value(GTK_SPIN_BUTTON(shape_l), model->gulp.cosmo_shape_index[1]);
gtk_spin_button_set_value(GTK_SPIN_BUTTON(segments), model->gulp.cosmo_segments);

gulp_cosmo_points_refresh(NULL, dialog);
gulp_cosmo_panel_refresh(cosmo, dialog);

gtk_widget_show_all(box);
}

/*******************************/
/* GULP job setup/results page */
/*******************************/
void gui_gulp_widget(GtkWidget *box, gpointer dialog)
{
gchar *text;
GSList *keyword[5];
GString *line;
GtkWidget *hbox, *vbox, *vbox1, *vbox2, *page, *swin;
GtkWidget *frame, *button, *label, *entry, *notebook;
struct model_pak *model;

model = dialog_model(dialog);
g_assert(model != NULL);

/* string manipulation scratchpad */
line = g_string_new(NULL);

/* create notebook */
notebook = gtk_notebook_new();
gtk_notebook_set_tab_pos(GTK_NOTEBOOK(notebook), GTK_POS_TOP);
gtk_container_add(GTK_CONTAINER(box), notebook);
gtk_notebook_set_show_border(GTK_NOTEBOOK(notebook), FALSE);

/* run type page */
page = gtk_vbox_new(FALSE, 0);
label = gtk_label_new (" Control ");
gtk_notebook_append_page(GTK_NOTEBOOK(notebook), page, label);
hbox = gtk_hbox_new(FALSE, 0);
gtk_container_add(GTK_CONTAINER(page), hbox);

/* split panel */
vbox1 = gtk_vbox_new(FALSE, 0);
gtk_box_pack_start(GTK_BOX(hbox), vbox1, TRUE, TRUE, 0);
vbox2 = gtk_vbox_new(FALSE, 0);
gtk_box_pack_start(GTK_BOX(hbox), vbox2, TRUE, TRUE, 0);

/* frame */
frame = gtk_frame_new("Run type");
gtk_box_pack_start(GTK_BOX(vbox1), frame, FALSE, FALSE, 0);
gtk_container_set_border_width(GTK_CONTAINER(frame), PANEL_SPACING/2);
vbox = gtk_vbox_new(TRUE, 0);
gtk_container_add(GTK_CONTAINER(frame), vbox);

/* do the radio buttons */
new_radio_group(0, vbox, TT);
button = add_radio_button("Single point", (gpointer) gulp_run_single, model);
if (model->gulp.run == E_SINGLE)
  gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(button), TRUE);

button = add_radio_button("Optimise", (gpointer) gulp_run_optimize, model);
if (model->gulp.run == E_OPTIMIZE)
  gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(button), TRUE);

button = add_radio_button("Dynamics", (gpointer) gulp_run_dynamics, model);
if (model->gulp.run == MD)
  gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(button), TRUE);

/* method constraint */
frame = gtk_frame_new(" Constraint ");
gtk_box_pack_start(GTK_BOX(vbox1), frame, FALSE, FALSE, 0);
gtk_container_set_border_width(GTK_CONTAINER(frame), PANEL_SPACING/2);
vbox = gtk_vbox_new(TRUE, 0);
gtk_container_add(GTK_CONTAINER(frame),vbox);

/* do the first (DEFAULT MODE) radio button */
button = gtk_radio_button_new_with_label(NULL, "Constant pressure");
gtk_box_pack_start (GTK_BOX (vbox), button, TRUE, TRUE, 0);
g_signal_connect(GTK_OBJECT(button), "clicked",
                 GTK_SIGNAL_FUNC (cb_gulp_keyword), (gpointer) button);
g_object_set_data(G_OBJECT(button), "key", (gpointer) CONP);
g_object_set_data(G_OBJECT(button), "ptr", (gpointer) model);
if (model->gulp.method == CONP)
  gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(button), TRUE);

/* make a radio group */
keyword[1] = gtk_radio_button_group(GTK_RADIO_BUTTON (button));
/* do the next button */
button = gtk_radio_button_new_with_label (keyword[1], "Constant volume");
gtk_box_pack_start(GTK_BOX (vbox), button, TRUE, TRUE, 0);
g_signal_connect(GTK_OBJECT(button), "clicked",
                 GTK_SIGNAL_FUNC (cb_gulp_keyword), (gpointer) button);
g_object_set_data(G_OBJECT(button), "key", (gpointer) CONV);
g_object_set_data(G_OBJECT(button), "ptr", (gpointer) model);
if (model->gulp.method == CONV)
  gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(button), TRUE);

/* frame */
frame = gtk_frame_new(" Molecule options ");
gtk_box_pack_start(GTK_BOX(vbox1), frame, FALSE, FALSE, 0);
gtk_container_set_border_width(GTK_CONTAINER(frame), PANEL_SPACING/2);
vbox = gtk_vbox_new(TRUE, 0);
gtk_container_add(GTK_CONTAINER(frame),vbox);
/* do the first radio button */
button = gtk_radio_button_new_with_label(NULL,
       "Coulomb subtract all intramolecular");
gtk_box_pack_start (GTK_BOX (vbox), button, TRUE, TRUE, 0);
g_signal_connect(GTK_OBJECT (button), "clicked",
                 GTK_SIGNAL_FUNC (cb_gulp_keyword), (gpointer) button);
g_object_set_data(G_OBJECT(button), "key", (gpointer) MOLE);
g_object_set_data(G_OBJECT(button), "ptr", (gpointer) model);
if (model->gulp.coulomb == MOLE)
  gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(button), TRUE);

/* make a radio group */
keyword[2] = gtk_radio_button_group(GTK_RADIO_BUTTON (button));
/* do the next button */
button = gtk_radio_button_new_with_label(keyword[2], 
       "Coulomb subtract 1-2 and 1-3 intramolecular");
gtk_box_pack_start (GTK_BOX (vbox), button, TRUE, TRUE, 0);
g_signal_connect(GTK_OBJECT(button), "clicked",
                 GTK_SIGNAL_FUNC(cb_gulp_keyword), (gpointer) button);
g_object_set_data(G_OBJECT(button), "key", (gpointer) MOLMEC);
g_object_set_data(G_OBJECT(button), "ptr", (gpointer) model);
if (model->gulp.coulomb == MOLMEC)
  gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(button), TRUE);

/* subsequent button... */
button = gtk_radio_button_new_with_label(gtk_radio_button_group
       (GTK_RADIO_BUTTON (button)), "Build, but retain coulomb interactions");
gtk_box_pack_start(GTK_BOX (vbox), button, TRUE, TRUE, 0);
g_signal_connect(GTK_OBJECT (button), "clicked",
                 GTK_SIGNAL_FUNC (cb_gulp_keyword), (gpointer) button);
g_object_set_data(G_OBJECT(button), "key", (gpointer) MOLQ);
g_object_set_data(G_OBJECT(button), "ptr", (gpointer) model);
if (model->gulp.coulomb == MOLQ)
  gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(button), TRUE);

/* subsequent button... */
button = gtk_radio_button_new_with_label(gtk_radio_button_group 
          (GTK_RADIO_BUTTON (button)), "Molecule building off");
gtk_box_pack_start(GTK_BOX (vbox), button, TRUE, TRUE, 0);
g_signal_connect(GTK_OBJECT (button), "clicked",
                 GTK_SIGNAL_FUNC (cb_gulp_keyword), (gpointer) button);
g_object_set_data(G_OBJECT(button), "key", (gpointer) NOBUILD);
g_object_set_data(G_OBJECT(button), "ptr", (gpointer) model);
if (model->gulp.coulomb == NOBUILD)
  gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(button), TRUE);
gui_direct_check("Fix the initial connectivity", &model->gulp.fix, NULL, NULL, vbox);
gui_direct_check("Automatic connectivity off", &model->gulp.noautobond, NULL, NULL, vbox);

/* temperature & pressure */
frame = gtk_frame_new(NULL);
gtk_box_pack_start(GTK_BOX(vbox1), frame, FALSE, FALSE, 0);
gtk_container_set_border_width(GTK_CONTAINER(frame), PANEL_SPACING/2);
vbox = gtk_vbox_new(TRUE, 0);
gtk_container_add(GTK_CONTAINER(frame),vbox);
gui_text_entry("Temperature", &model->gulp.temperature, TRUE, FALSE, vbox);
gui_text_entry("Pressure", &model->gulp.pressure, TRUE, FALSE, vbox);


/* dynamics ensemble */
frame = gtk_frame_new("Dynamics");
gtk_box_pack_start(GTK_BOX(vbox2), frame, FALSE, FALSE, 0);
gtk_container_set_border_width(GTK_CONTAINER(frame), PANEL_SPACING/2);
vbox = gtk_vbox_new(TRUE, 0);
gtk_container_add(GTK_CONTAINER(frame),vbox);

/* do the first (DEFAULT MODE) radio button */
button = gtk_radio_button_new_with_label(NULL, "NVE");
gtk_box_pack_start(GTK_BOX(vbox), button, TRUE, TRUE, 0);
g_signal_connect(GTK_OBJECT(button), "clicked",
                 GTK_SIGNAL_FUNC (cb_gulp_keyword), (gpointer) button);
g_object_set_data(G_OBJECT(button), "key", (gpointer) NVE);
g_object_set_data(G_OBJECT(button), "ptr", (gpointer) model);
if (model->gulp.ensemble == NVE)
  gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(button), TRUE);

/* make a radio group */
keyword[1] = gtk_radio_button_group(GTK_RADIO_BUTTON (button));
/* do the next button */
button = gtk_radio_button_new_with_label (keyword[1], "NVT");
gtk_box_pack_start(GTK_BOX(vbox), button, TRUE, TRUE, 0);
g_signal_connect(GTK_OBJECT(button), "clicked",
                 GTK_SIGNAL_FUNC(cb_gulp_keyword), (gpointer) button);
g_object_set_data(G_OBJECT(button), "key", (gpointer) NVT);
g_object_set_data(G_OBJECT(button), "ptr", (gpointer) model);
/* NEW */
if (model->gulp.ensemble == NVT)
  gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(button), TRUE);

/* subsequent button... */
button = gtk_radio_button_new_with_label(gtk_radio_button_group
       (GTK_RADIO_BUTTON (button)), "NPT");
gtk_box_pack_start(GTK_BOX(vbox), button, TRUE, TRUE, 0);
g_signal_connect(GTK_OBJECT(button), "clicked",
                 GTK_SIGNAL_FUNC (cb_gulp_keyword), (gpointer) button);
g_object_set_data(G_OBJECT(button), "key", (gpointer) NPT);
g_object_set_data(G_OBJECT(button), "ptr", (gpointer) model);
if (model->gulp.ensemble == NPT)
  gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(button), TRUE);

/* dynamics sampling */
gtk_box_pack_start(GTK_BOX(vbox), gtk_hseparator_new(), FALSE, FALSE, 0);

gui_text_entry("Timestep", &model->gulp.timestep, TRUE, FALSE, vbox);
gui_text_entry("Equilibration  ", &model->gulp.equilibration, TRUE, FALSE, vbox);
gui_text_entry("Production", &model->gulp.production, TRUE, FALSE, vbox);
gui_text_entry("Sample", &model->gulp.sample, TRUE, FALSE, vbox);
gui_text_entry("Write", &model->gulp.write, TRUE, FALSE, vbox);

/* frame for keyword toggles */
frame = gtk_frame_new(NULL);
gtk_box_pack_start(GTK_BOX(vbox2), frame, FALSE, FALSE, 0);
gtk_container_set_border_width(GTK_CONTAINER(frame), PANEL_SPACING/2);
vbox = gtk_vbox_new(TRUE, 0);
gtk_container_add(GTK_CONTAINER(frame), vbox);

gui_direct_check("Create input file then stop", &model->gulp.no_exec, NULL, NULL, vbox);
gui_direct_check("Build cell, then discard symmetry", &model->gulp.nosym, NULL, NULL, vbox);

gui_direct_check("No attachment energy calculation", &model->gulp.no_eatt, NULL, NULL, vbox);
gui_direct_check("QEq Electronegativity equalisation", &model->gulp.qeq, NULL, NULL, vbox);


/* NEW - text box for GULP control files */
hbox = gtk_hbox_new(FALSE,0);
label = gtk_label_new(" Files ");
gtk_notebook_append_page(GTK_NOTEBOOK(notebook), hbox, label);
vbox = gtk_vbox_new(FALSE,0);
gtk_box_pack_start(GTK_BOX(hbox), vbox, FALSE, FALSE, 0);
frame = gtk_frame_new(NULL);
gtk_box_pack_start(GTK_BOX(vbox), frame, FALSE, FALSE, 0);
gtk_container_set_border_width(GTK_CONTAINER(frame), PANEL_SPACING);
vbox = gtk_vbox_new(FALSE, PANEL_SPACING/2);
gtk_container_add(GTK_CONTAINER(frame), vbox);

/* GULP file labels */
gui_text_entry("Input file", &model->gulp.temp_file, FALSE, FALSE, vbox);
gui_text_entry("Dump file", &model->gulp.dump_file, FALSE, FALSE, vbox);
gui_text_entry("Trajectory file  ", &model->gulp.trj_file, FALSE, FALSE, vbox);

/*
gtk_box_pack_start(GTK_BOX(vbox), gtk_hseparator_new(), FALSE, FALSE, 0);

hbox = gtk_hbox_new(FALSE, 0);
gtk_box_pack_start(GTK_BOX(vbox), hbox, FALSE, FALSE, 0);
label = gtk_label_new("Customize job name  ");
gtk_box_pack_start(GTK_BOX(hbox), label, FALSE, FALSE, 0);
entry = gtk_entry_new();
gtk_box_pack_end(GTK_BOX(hbox), entry, FALSE, FALSE, 0);
*/

/* FIXME - stopped using basename, because it may have spaces in it - which */
/* cause problems reading/writing filenames (ie not quoted) */
//gtk_entry_set_text(GTK_ENTRY(entry), model->basename);
/*
gtk_entry_set_text(GTK_ENTRY(entry), "");
g_signal_connect(GTK_OBJECT(entry), "changed",
                 GTK_SIGNAL_FUNC(gulp_jobname_changed), (gpointer) model);

*/

/* CURRENT - advanced options */
frame = gtk_frame_new("Options");
gtk_box_pack_start(GTK_BOX(vbox), frame, FALSE, FALSE, 0);
gtk_container_set_border_width(GTK_CONTAINER(frame), PANEL_SPACING);
vbox = gtk_vbox_new(FALSE, PANEL_SPACING/2);
gtk_container_add(GTK_CONTAINER(frame), vbox);


gui_direct_check("Print charge with coordinates ", &model->gulp.print_charge, NULL, NULL, vbox);


/* minimize page */
page = gtk_vbox_new(FALSE,0);
label = gtk_label_new(" Optimisation ");
gtk_notebook_append_page(GTK_NOTEBOOK(notebook), page, label);
hbox = gtk_hbox_new(FALSE, 0);
gtk_box_pack_start(GTK_BOX(page), hbox, FALSE, FALSE, 0);

/* optimiser to use */
frame = gtk_frame_new("Primary optimiser");
gtk_box_pack_start(GTK_BOX(hbox), frame, TRUE, TRUE, 0);
gtk_container_set_border_width(GTK_CONTAINER(frame), PANEL_SPACING);
vbox = gtk_vbox_new(TRUE, 0);
gtk_container_add(GTK_CONTAINER(frame), vbox);

/* do the first (DEFAULT MODE) radio button */
button = gtk_radio_button_new_with_label(NULL, "bfgs");
gtk_box_pack_start(GTK_BOX(vbox), button, TRUE, TRUE, 0);
g_signal_connect(GTK_OBJECT(button), "clicked",
                 GTK_SIGNAL_FUNC(cb_gulp_keyword), (gpointer) button);
g_object_set_data(G_OBJECT(button), "key", (gpointer) BFGS_OPT);
g_object_set_data(G_OBJECT(button), "ptr", (gpointer) model);
if (model->gulp.optimiser == BFGS_OPT || model->gulp.optimiser == -1)
  gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(button), TRUE);

/* make a radio group */
keyword[1] = gtk_radio_button_group(GTK_RADIO_BUTTON (button));
/* do the next button */
button = gtk_radio_button_new_with_label (keyword[1], "conj");
gtk_box_pack_start(GTK_BOX(vbox), button, TRUE, TRUE, 0);
g_signal_connect(GTK_OBJECT(button), "clicked",
                 GTK_SIGNAL_FUNC(cb_gulp_keyword), (gpointer) button);
g_object_set_data(G_OBJECT(button), "key", (gpointer) CONJ_OPT);
g_object_set_data(G_OBJECT(button), "ptr", (gpointer) model);
if (model->gulp.optimiser == CONJ_OPT)
  gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(button), TRUE);

/* subsequent button... */
button = gtk_radio_button_new_with_label(gtk_radio_button_group
       (GTK_RADIO_BUTTON(button)), "rfo");
gtk_box_pack_start(GTK_BOX(vbox), button, TRUE, TRUE, 0);
g_signal_connect(GTK_OBJECT(button), "clicked",
                 GTK_SIGNAL_FUNC (cb_gulp_keyword), (gpointer) button);
g_object_set_data(G_OBJECT(button), "key", (gpointer) RFO_OPT);
g_object_set_data(G_OBJECT(button), "ptr", (gpointer) model);
if (model->gulp.optimiser == RFO_OPT)
  gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(button), TRUE);


/* optimiser to use */
frame = gtk_frame_new("Secondary optimiser");
gtk_box_pack_start(GTK_BOX(hbox), frame, TRUE, TRUE, 0);
gtk_container_set_border_width(GTK_CONTAINER(frame), PANEL_SPACING);
vbox = gtk_vbox_new(TRUE, 0);
gtk_container_add(GTK_CONTAINER(frame), vbox);

/* do the first (DEFAULT MODE) radio button */
button = gtk_radio_button_new_with_label(NULL, "none");
gtk_box_pack_start(GTK_BOX(vbox), button, FALSE, FALSE, 0);
g_signal_connect(GTK_OBJECT(button), "clicked",
                 GTK_SIGNAL_FUNC (switch_keyword), (gpointer) button);
g_object_set_data(G_OBJECT(button), "key", (gpointer) SWITCH_OFF);
g_object_set_data(G_OBJECT(button), "ptr", (gpointer) model);
if (model->gulp.optimiser2 == SWITCH_OFF || model->gulp.optimiser2 == -1)
  gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(button), TRUE);

/* make a radio group */
keyword[1] = gtk_radio_button_group(GTK_RADIO_BUTTON(button));
/* do the next button */
button = gtk_radio_button_new_with_label(keyword[1], "bfgs");
gtk_box_pack_start(GTK_BOX(vbox), button, FALSE, FALSE, 0);
g_signal_connect(GTK_OBJECT(button), "clicked",
                 GTK_SIGNAL_FUNC(switch_keyword), (gpointer) button);
g_object_set_data(G_OBJECT(button), "key", (gpointer) BFGS_OPT);
g_object_set_data(G_OBJECT(button), "ptr", (gpointer) model);
if (model->gulp.optimiser2 == BFGS_OPT)
  gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(button), TRUE);

/* subsequent button... */
button = gtk_radio_button_new_with_label(gtk_radio_button_group
                                        (GTK_RADIO_BUTTON(button)), "conj");
gtk_box_pack_start(GTK_BOX(vbox), button, FALSE, FALSE, 0);
g_signal_connect(GTK_OBJECT(button), "clicked",
                 GTK_SIGNAL_FUNC(switch_keyword), (gpointer) button);
g_object_set_data(G_OBJECT(button), "key", (gpointer) CONJ_OPT);
g_object_set_data(G_OBJECT(button), "ptr", (gpointer) model);
if (model->gulp.optimiser2 == CONJ_OPT)
  gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(button), TRUE);

/* subsequent button... */
button = gtk_radio_button_new_with_label(gtk_radio_button_group
                                        (GTK_RADIO_BUTTON(button)), "rfo");
gtk_box_pack_start(GTK_BOX(vbox), button, FALSE, FALSE, 0);
g_signal_connect(GTK_OBJECT(button), "clicked",
                 GTK_SIGNAL_FUNC(switch_keyword), (gpointer) button);
g_object_set_data(G_OBJECT(button), "key", (gpointer) RFO_OPT);
g_object_set_data(G_OBJECT(button), "ptr", (gpointer) model);
if (model->gulp.optimiser2 == RFO_OPT)
  gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(button), TRUE);


/* optimiser to use */
vbox = gui_frame_vbox("Switching criteria", TRUE, TRUE, hbox);

/* do the first (DEFAULT MODE) radio button */
button = gtk_radio_button_new_with_label(NULL, "cycle");
gtk_box_pack_start(GTK_BOX(vbox), button, FALSE, FALSE, 0);
g_signal_connect(GTK_OBJECT(button), "clicked",
                 GTK_SIGNAL_FUNC(switch_keyword), (gpointer) button);
g_object_set_data(G_OBJECT(button), "key", (gpointer) CYCLE);
g_object_set_data(G_OBJECT(button), "ptr", (gpointer) model);
if (model->gulp.switch_type == CYCLE)
  gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(button), TRUE);

/* make a radio group */
keyword[1] = gtk_radio_button_group(GTK_RADIO_BUTTON(button));
/* do the next button */
button = gtk_radio_button_new_with_label(keyword[1], "gnorm");
gtk_box_pack_start(GTK_BOX(vbox), button, FALSE, FALSE, 0);
g_signal_connect(GTK_OBJECT(button), "clicked",
                 GTK_SIGNAL_FUNC(switch_keyword), (gpointer) button);
g_object_set_data(G_OBJECT(button), "key", (gpointer) GNORM);
g_object_set_data(G_OBJECT(button), "ptr", (gpointer) model);
if (model->gulp.switch_type == GNORM || model->gulp.switch_type < 0)
  gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(button), TRUE);

entry = gtk_entry_new();
text = g_strdup_printf("%f", model->gulp.switch_value);
gtk_entry_set_text(GTK_ENTRY(entry), text);
g_free(text);
gtk_box_pack_end(GTK_BOX(vbox), entry, TRUE, FALSE, 0);
g_signal_connect(GTK_OBJECT(entry), "changed",
                 GTK_SIGNAL_FUNC(switch_value_changed), (gpointer) entry);
g_object_set_data(G_OBJECT(entry), "ptr", (gpointer) model);

/* frame */
hbox = gui_frame_hbox(NULL, FALSE, FALSE, page);

/* optimization cycles */
label = gtk_label_new(" Maximum cycles ");
gtk_box_pack_start(GTK_BOX(hbox), label, FALSE, FALSE, 0);
gui_direct_spin(NULL, &model->gulp.maxcyc, 0, 500, 10, NULL, NULL, hbox);

/* NEW - text box for potentials */
vbox = gtk_vbox_new(FALSE,0);
label = gtk_label_new(" Potentials ");
gtk_notebook_append_page(GTK_NOTEBOOK(notebook), vbox, label);

frame = gtk_frame_new(NULL);
gtk_box_pack_start(GTK_BOX(vbox), frame, TRUE, TRUE, 0);
swin = gui_text_window(&model->gulp.potentials, TRUE);
gtk_container_add(GTK_CONTAINER(frame), swin);

gui_text_entry(" Potential Library", &model->gulp.libfile, TRUE, FALSE, vbox);

/* NEW - text box for element and species data */
hbox = gtk_hbox_new(TRUE, 0);
label = gtk_label_new(" Elements ");
gtk_notebook_append_page(GTK_NOTEBOOK(notebook), hbox, label);

frame = gtk_frame_new("Element");
gtk_box_pack_start(GTK_BOX(hbox), frame, TRUE, TRUE, 0);
swin = gui_text_window(&model->gulp.elements, TRUE);
gtk_container_add(GTK_CONTAINER(frame), swin);

frame = gtk_frame_new("Species");
gtk_box_pack_start(GTK_BOX(hbox), frame, TRUE, TRUE, 0);
swin = gui_text_window(&model->gulp.species, TRUE);
gtk_container_add(GTK_CONTAINER(frame), swin);

/* text box for "pass through" GULP data */
vbox = gtk_vbox_new(FALSE,0);
label = gtk_label_new(" Unprocessed ");
gtk_notebook_append_page(GTK_NOTEBOOK(notebook), vbox, label);

/* unprocessed keywords */
vbox1 = gui_frame_vbox(NULL, TRUE, TRUE, vbox);

hbox = gtk_hbox_new(TRUE, 0);
gtk_box_pack_start(GTK_BOX(vbox1), hbox, FALSE, FALSE, 0);
gui_direct_check("Pass keywords through to GULP", &model->gulp.output_extra_keywords, NULL, NULL, hbox);
gui_text_entry("Keywords  ", &model->gulp.extra_keywords, TRUE, TRUE, vbox1);

/* unprocessed options */
hbox = gtk_hbox_new(TRUE, 0);
gtk_box_pack_start(GTK_BOX(vbox1), hbox, FALSE, FALSE, 0);
gui_direct_check("Pass options through to GULP", &model->gulp.output_extra, NULL, NULL, hbox);
swin = gui_text_window(&model->gulp.extra, TRUE);
gtk_box_pack_start(GTK_BOX(vbox1), swin, TRUE, TRUE, 0);

/* vibrational eigenvector display page */
page = gtk_vbox_new(FALSE,0);
label = gtk_label_new (" Vibrational ");
gtk_notebook_append_page(GTK_NOTEBOOK(notebook), page, label);
gulp_phonon_box(page, model);


/* CURRENT - solvation (cosmic) */
/* TODO - only available on some versions of GULP */
/* FIXME - windows version generates GTK CRITICALS when trying to invoke this */
#ifndef _WIN32
page = gtk_vbox_new(FALSE,0);
label = gtk_label_new (" Solvation ");
gtk_notebook_append_page(GTK_NOTEBOOK(notebook), page, label);
gulp_solvation_box(page, dialog);
#endif


/* job input filename */
hbox = gui_frame_hbox("Files", FALSE, FALSE, box);

label = gtk_label_new("Job file name  ");
gtk_box_pack_start(GTK_BOX(hbox), label, TRUE, TRUE, 0);
entry = gtk_entry_new();
gtk_box_pack_end(GTK_BOX(hbox), entry, TRUE, TRUE, 0);

//gtk_entry_set_text(GTK_ENTRY(entry), model->gulp.temp_file);
//gtk_entry_set_text(GTK_ENTRY(entry), "");

/* NB: GULP job file name entry has extension removed */
text = parse_strip_extension(model->gulp.temp_file);
gtk_entry_set_text(GTK_ENTRY(entry), text);
g_free(text);

g_signal_connect(GTK_OBJECT(entry), "changed",
                 GTK_SIGNAL_FUNC(gulp_jobname_changed), (gpointer) model);

/* summary details for the model */
hbox = gui_frame_hbox("Details", FALSE, FALSE, box);

/* left vbox */
vbox = gtk_vbox_new(TRUE, 2);
gtk_box_pack_start (GTK_BOX (hbox), vbox, TRUE, TRUE, 0);

label = gtk_label_new("Structure name");
gtk_box_pack_start (GTK_BOX (vbox), label, TRUE, TRUE, 0);

label = gtk_label_new("Total energy (eV)");
gtk_box_pack_start (GTK_BOX (vbox), label, TRUE, TRUE, 0);
if (model->periodic == 2)
  {
  label = gtk_label_new("Surface bulk energy (eV)");
  gtk_box_pack_start (GTK_BOX (vbox), label, TRUE, TRUE, 0);
  label = gtk_label_new("Surface dipole");
  gtk_box_pack_start (GTK_BOX (vbox), label, TRUE, TRUE, 0);
  label = gtk_label_new("Surface energy");
  gtk_box_pack_start (GTK_BOX (vbox), label, TRUE, TRUE, 0);
  label = gtk_label_new("Attachment energy");
  gtk_box_pack_start (GTK_BOX (vbox), label, TRUE, TRUE, 0);
  }

/* right vbox */
vbox = gtk_vbox_new(TRUE, 2);
gtk_box_pack_start(GTK_BOX(hbox), vbox, TRUE, TRUE, 0);


entry = gtk_entry_new();
gtk_entry_set_text(GTK_ENTRY(entry), model->basename);
gtk_box_pack_start(GTK_BOX(vbox), entry, FALSE, FALSE, 0);
gtk_entry_set_editable(GTK_ENTRY(entry), FALSE);


model->gulp.energy_entry = gtk_entry_new_with_max_length(LINELEN);
if (model->gulp.free)
  g_string_sprintf(line,"%f (free energy)",model->gulp.energy);
else
  g_string_sprintf(line,"%f",model->gulp.energy);

gtk_entry_set_text(GTK_ENTRY(model->gulp.energy_entry),line->str);
gtk_box_pack_start (GTK_BOX (vbox), model->gulp.energy_entry, TRUE, TRUE, 0);
gtk_entry_set_editable(GTK_ENTRY(model->gulp.energy_entry), FALSE);

if (model->periodic == 2)
  {
  model->gulp.sbe_entry = gtk_entry_new_with_max_length(LINELEN);
  g_string_sprintf(line,"%f",model->gulp.sbulkenergy);
  gtk_entry_set_text(GTK_ENTRY(model->gulp.sbe_entry),line->str);
  gtk_box_pack_start(GTK_BOX(vbox), model->gulp.sbe_entry, TRUE, TRUE, 0);

/* NEW - sbulkenergy editable */
gtk_entry_set_editable(GTK_ENTRY(model->gulp.sbe_entry), TRUE);
g_signal_connect(GTK_OBJECT(model->gulp.sbe_entry), "changed",
                 GTK_SIGNAL_FUNC(cb_gulp_keyword), (gpointer) entry);
g_object_set_data(G_OBJECT(model->gulp.sbe_entry), "key", (gpointer) SBULK_ENERGY);
g_object_set_data(G_OBJECT(model->gulp.sbe_entry), "ptr", (gpointer) model);

/* surface dipole */
  model->gulp.sdipole_entry = gtk_entry_new_with_max_length(LINELEN);
  g_string_sprintf(line,"%f e.Angs", model->gulp.sdipole);
  gtk_entry_set_text(GTK_ENTRY(model->gulp.sdipole_entry), line->str);
  gtk_box_pack_start(GTK_BOX(vbox), model->gulp.sdipole_entry, TRUE, TRUE, 0);
  gtk_entry_set_editable(GTK_ENTRY(model->gulp.sdipole_entry), FALSE);
/* surface energy */
  model->gulp.esurf_entry = gtk_entry_new_with_max_length(LINELEN);
  if (model->gulp.esurf[1] == 0.0)
     g_string_sprintf(line,"%f    %s",
                            model->gulp.esurf[0], model->gulp.esurf_units);
  else
     g_string_sprintf(line,"%f    %s",
                            model->gulp.esurf[1], model->gulp.esurf_units);

  gtk_entry_set_text(GTK_ENTRY(model->gulp.esurf_entry), line->str);
  gtk_box_pack_start(GTK_BOX(vbox), model->gulp.esurf_entry, TRUE, TRUE, 0);
  gtk_entry_set_editable(GTK_ENTRY(model->gulp.esurf_entry), FALSE);
/* attachment energy */
  model->gulp.eatt_entry = gtk_entry_new_with_max_length(LINELEN);
  if (model->gulp.eatt[1] == 0.0)
     g_string_sprintf(line,"%f    %s",
                            model->gulp.eatt[0], model->gulp.eatt_units);
  else
     g_string_sprintf(line,"%f    %s",
                            model->gulp.eatt[1], model->gulp.eatt_units);
  gtk_entry_set_text(GTK_ENTRY(model->gulp.eatt_entry), line->str);
  gtk_box_pack_start(GTK_BOX(vbox), model->gulp.eatt_entry, TRUE, TRUE, 0);
  gtk_entry_set_editable(GTK_ENTRY(model->gulp.eatt_entry), FALSE);
  }

/* done */
gtk_widget_show_all(box);

g_string_free(line, TRUE);
gui_model_select(model);
}

/****************************************/
/* run the approriate energetics method */
/****************************************/
void gulp_execute(GtkWidget *w, gpointer data)
{
struct model_pak *model = data;

g_assert(model != NULL);

gui_gulp_task(NULL, model);
}

/***********************/
/* grid job submission */
/***********************/
void gui_gulp_grid_submit(GtkWidget *w, gpointer data)
{
#ifdef WITH_GRISU
gchar *tmp;
struct model_pak *model = data;
struct grid_pak *grid;

g_assert(model != NULL);

grid = grid_new();
grid->local_cwd = g_strdup(sysenv.cwd);
grid->local_input = g_strdup(model->gulp.temp_file);

/* immediately write to disk so job is preserved even if user deletes model */
tmp = g_build_filename(grid->local_cwd, grid->local_input, NULL);
write_gulp(tmp, model);
g_free(tmp);

if (!grid_job_submit("gulp", grid))
  {
/* TODO - failed message */
  grid_free(grid);
  }
#endif
}

/*********************************/
/* cleanup for energetics dialog */
/*********************************/
void gulp_cleanup(struct model_pak *model)
{
g_assert(model != NULL);

model->gulp.energy_entry = NULL;
model->gulp.sdipole_entry = NULL;
model->gulp.esurf_entry = NULL;
model->gulp.eatt_entry = NULL;
}

/**************************/
/* GULP energetics dialog */
/**************************/
void gulp_dialog(void)
{
gpointer dialog;
GtkWidget *window, *vbox;
struct model_pak *model;

model = sysenv.active_model;
if (!model)
  return;

gulp_files_init(model);

/* request an energetics dialog */
dialog = dialog_request(GULP, "GULP configuration", NULL, gulp_cleanup, model);
if (!dialog)
  return;
window = dialog_window(dialog);

vbox = gui_frame_vbox(NULL, FALSE, FALSE, GTK_DIALOG(window)->vbox);

gui_gulp_widget(vbox, dialog);

/* terminating button */
gui_stock_button(GTK_STOCK_EXECUTE, gulp_execute, model,
                   GTK_DIALOG(window)->action_area);

#ifdef WITH_GRISU
gui_icon_button(GTK_STOCK_APPLY, "Grid Submit", gui_gulp_grid_submit, model, GTK_DIALOG(window)->action_area);
#endif

gui_stock_button(GTK_STOCK_CLOSE, dialog_destroy, dialog,
                   GTK_DIALOG(window)->action_area);

gtk_widget_show_all(window);
}

