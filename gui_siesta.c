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
#include "numeric.h"
#include "parse.h"
#include "project.h"
#include "render.h"
#include "spatial.h"
#include "quaternion.h"
#include "task.h"
#include "interface.h"
#include "dialog.h"
#include "gui_shorts.h"
/*
#include "meschach/matrix-meschach.h"
*/
#include "mesch.h"
#include "gui_siesta.h"

extern struct sysenv_pak sysenv;
extern struct elem_pak elements[];

static gboolean siestafileWRITE;

/******************************/
/* toggle eigenvector display */
/******************************/
void gui_siesta_mode_show(GtkWidget *w, gpointer dialog)
{
gint i, atom, state;
gdouble scale, x1[3], x2[3], colour[3];
gpointer button, spin;
GSList *list;
struct core_pak *core;
struct spatial_pak *spatial;
struct model_pak *model;

g_assert(dialog != NULL);

model = dialog_model(dialog);
g_assert(model != NULL);

button = dialog_child_get(dialog, "phonon_toggle");
g_assert(button != NULL);

spatial_destroy_by_label("siesta_phonons", model);

state = gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(button));
if (!state)
  {
  redraw_canvas(SINGLE);
  return;
  }

spin = dialog_child_get(dialog, "phonon_scaling");
scale = SPIN_FVAL(spin);

/* create & init the spatial object */
spatial = spatial_new("siesta_phonons", SPATIAL_VECTOR, 2, TRUE, model);

atom = 0;
/* get eigenvectors from all atoms */
for (list=model->cores ; list; list=g_slist_next(list))
  {
  core = list->data;
  ARR3SET(x1, core->x);

/* get current eigenvector */
  i = model->siesta.sorted_eig_values[model->siesta.current_animation];
  x2[0] = mesch_me_get(model->siesta.eigen_xyz_atom_mat, 3*atom, i);
  x2[1] = mesch_me_get(model->siesta.eigen_xyz_atom_mat, 3*atom+1, i);
  x2[2] = mesch_me_get(model->siesta.eigen_xyz_atom_mat, 3*atom+2, i);
  atom++;

/* compute coords */
  VEC3MUL(x2, scale);
  vecmat(model->ilatmat, x2);
  ARR3ADD(x2, x1);
/* add to spatial */
  spatial_vertex_add(x2, colour, spatial);
  spatial_vertex_add(x1, colour, spatial);
  }
    
/* drawing update */
coords_compute(model);
redraw_canvas(SINGLE);
}

/*******************************/
/* update frequency entry text */
/*******************************/
void gui_siesta_frequency_set(gpointer dialog)
{
gint m;
gchar *text;
gdouble f;
gpointer entry;
struct model_pak *model;

model = dialog_model(dialog);
g_assert(model != NULL);
entry = dialog_child_get(dialog, "phonon_entry");
g_assert(entry != NULL);

m = model->siesta.current_animation;

g_assert(m >= 0);
g_assert(m < model->siesta.num_animations);

f = make_eigvec_freq(mesch_ve_get(model->siesta.eigen_values, model->siesta.sorted_eig_values[m]));

text = g_strdup_printf("%f", f);
gtk_entry_set_text(GTK_ENTRY(entry), text);
g_free(text);
}

/***************************/
/* slider bar event change */
/***************************/
void gui_siesta_slider_changed(GtkWidget *w, gpointer dialog)
{
gint m;
struct model_pak *model;

g_assert(w != NULL);
g_assert(dialog != NULL);

model = dialog_model(dialog);

m = nearest_int(gtk_range_get_value(GTK_RANGE(w)));
model->siesta.current_animation = m-1;   /* 1..n range -> 0..n-1 */

gui_siesta_frequency_set(dialog);
gui_siesta_mode_show(NULL, dialog);
}

/*********************/
/* move to next mode */
/*********************/
void gui_siesta_mode_next(GtkWidget *w, gpointer dialog)
{
gdouble m;
gpointer slider;

slider = dialog_child_get(dialog, "phonon_slider");

m = gtk_range_get_value(GTK_RANGE(slider));
m++;
gtk_range_set_value(GTK_RANGE(slider), m);
}

/*************************/
/* move to previous mode */
/*************************/
void gui_siesta_mode_prev(GtkWidget *w, gpointer dialog)
{
gdouble m;
gpointer slider;

slider = dialog_child_get(dialog, "phonon_slider");

m = gtk_range_get_value(GTK_RANGE(slider));
m--;
gtk_range_set_value(GTK_RANGE(slider), m);
}

/***********************************/
/* cleanup a vibrational animation */
/***********************************/
void siesta_phonon_cleanup(struct model_pak *model)
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
gint siesta_phonon_timer(gpointer dialog)
{
static gint count=0;
gint atom, mode;
gchar *text, *name;
gdouble f, x1[3];
gpointer spin;
GSList *list;
struct core_pak *core;
struct model_pak *model;

/* checks */
if (!dialog_valid(dialog))
  return(FALSE);

model = dialog_model(dialog);
g_assert(model != NULL);

/* stop animation? */
if (!model->pulse_direction)
  {
  siesta_phonon_cleanup(model);
  coords_compute(model);
  redraw_canvas(SINGLE);
  return(FALSE);
  }

/* setup animation resolution */
spin = dialog_child_get(dialog, "phonon_resolution");
model->pulse_max = SPIN_FVAL(spin);

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

spin = dialog_child_get(dialog, "phonon_scaling");
f = SPIN_FVAL(spin);
f *= (gdouble) model->pulse_count;
f /= model->pulse_max;

atom = 0;

mode = model->siesta.sorted_eig_values[model->siesta.current_animation];

/* get eigenvectors from all atoms */
for (list=model->cores ; list; list=g_slist_next(list))
  {
  core = list->data;
  ARR3SET(x1, core->x);

//get x,y,z for given freq.
  core->offset[0] = mesch_me_get(model->siesta.eigen_xyz_atom_mat, 3*atom, mode);
  core->offset[1] = mesch_me_get(model->siesta.eigen_xyz_atom_mat, 3*atom+1, mode);
  core->offset[2] = mesch_me_get(model->siesta.eigen_xyz_atom_mat, 3*atom+2, mode);
  atom++;

/* pulse offset scaling */
  VEC3MUL(core->offset, f);
  vecmat(model->ilatmat, core->offset);
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

    text = g_strdup_printf("%s -delay 20 %s_*.tga %s.%s",
                           sysenv.convert_path, model->phonon_movie_name,
                           model->phonon_movie_name, model->phonon_movie_type);

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
void siesta_phonon_start(GtkWidget *w, gpointer dialog)
{
struct model_pak *model;

g_assert(dialog != NULL);
model = dialog_model(dialog);
g_assert(model != NULL);

model->pulse_count = 0;
model->pulse_direction = 1;

g_timeout_add(100, (GSourceFunc) siesta_phonon_timer, dialog);
}

/******************************/
/* stop phonon mode animation */
/*****************************/
void siesta_phonon_stop(GtkWidget *w, gpointer dialog)
{
struct model_pak *model;

g_assert(dialog != NULL);
model = dialog_model(dialog);
g_assert(model != NULL);

/* reset */
model->pulse_direction = 0;
model->pulse_count = 0;
model->phonon_movie = FALSE;
}

/**************************************/
/* create a movie of the current mode */
/**************************************/
void siesta_phonon_movie(GtkWidget *w, gpointer dialog)
{
gpointer entry;
struct model_pak *model;

model = dialog_model(dialog);
g_assert(model != NULL);

entry = dialog_child_get(dialog, "phonon_movie_name");
g_assert(entry != NULL);
/* FIXME - GULP phonon movie uses the same variable */
g_free(model->phonon_movie_name);
model->phonon_movie_name = g_strdup(gtk_entry_get_text(GTK_ENTRY(entry)));

entry = dialog_child_get(dialog, "phonon_movie_type");
g_assert(entry != NULL);
g_free(model->phonon_movie_type);
model->phonon_movie_type = g_strdup(gtk_entry_get_text(GTK_ENTRY(entry)));

model->phonon_movie = TRUE;
siesta_phonon_start(NULL, dialog);
}

/*****************************/
/* SIESTA phonon calculation */
/*****************************/
/* TODO - relocate */
gint siesta_phonon_calc(struct model_pak *model)
{
gint num_tokens, num_atoms, num_lines, i, j, k, l;
gchar ** buff;
gchar * modelFCname, *modelFCnameCSV;
gdouble wi, wj, value;
gpointer bigmat, correction_mat, mini_correction_zerooooo_mat;
GSList *list_i, *list_j;
struct core_pak *core_i;
struct core_pak *core_j;
FILE *fp, *matout=NULL;

g_assert(model != NULL);

/* check atom labels (since we must be able to do a valid weight lookup) */
for (list_i=model->cores ; list_i ; list_i=g_slist_next(list_i))
  {
  core_i = list_i->data;
  if (core_i->atom_code == 0)
    {
    gchar *text;

    text = g_strdup_printf("Unknown atom label: [%s]\n", core_i->atom_label); 
    gui_text_show(ERROR, text);
    g_free(text);

    return(1);
    }
  }

num_atoms = model->num_atoms;

//lines = 3 *    N * N     * 2;
//       xyz, each atom, back/forward
//
num_lines = 3*2*num_atoms*num_atoms;

modelFCname = g_strdup_printf("%s/%s.FC", sysenv.cwd, model->basename);
modelFCnameCSV = g_strdup_printf("%s.csv", modelFCname);

fp = fopen(modelFCname, "rt");

if (siestafileWRITE)
  {
  matout = fopen(modelFCnameCSV, "w");
  if (!matout)
    {
    gui_text_show(ERROR, "bugger - no save files\n");
    return(2);
    }
  }

if (!fp)
  {
  gchar * text;
  text = g_strdup_printf("*ERROR* - modelFCname file not opened\n");
  gui_text_show(ERROR, text);
  gui_text_show(ERROR, modelFCname);
  gui_text_show(ERROR, "\n");
  g_free(text);
  return(3);
  }

//no need for names anymore
g_free(modelFCname);

//initalise bigmat
bigmat = mesch_mat_new(3*num_atoms, 3*num_atoms);

//first line is crap.
buff = get_tokenized_line(fp, &num_tokens);

//FILE reading into bigmat
for (i = 0; i<num_lines; i++)
  {
  g_strfreev(buff);
  buff = get_tokenized_line(fp, &num_tokens);
  if (!buff)
    {
//error not enough lines in file...
//matrix_2d_free(bigmat);
    }

  if (num_tokens > 2)
    {
    if ( (i/num_atoms)%2 == 0)
      {
//first pass at row?

mesch_me_set(bigmat, i/(2*num_atoms), (3*i)%(3*num_atoms), str_to_float(*buff));
mesch_me_set(bigmat, i/(2*num_atoms), ((3*i)+1)%(3*num_atoms), str_to_float(*(buff+1)));
mesch_me_set(bigmat, i/(2*num_atoms), ((3*i)+2)%(3*num_atoms), str_to_float(*(buff+2)));

       }
     else
       {
//second pass - do the average

value = 0.5*(mesch_me_get(bigmat, i/(2*num_atoms), (3*i)%(3*num_atoms)) + str_to_float(*buff));
mesch_me_set(bigmat, i/(2*num_atoms), (3*i)%(3*num_atoms), value);

value = 0.5*(mesch_me_get(bigmat, i/(2*num_atoms), ((3*i)+1)%(3*num_atoms)) + str_to_float(*(buff+1)));
mesch_me_set(bigmat, i/(2*num_atoms), ((3*i)+1)%(3*num_atoms), value);

value = 0.5*(mesch_me_get(bigmat, i/(2*num_atoms), ((3*i)+2)%(3*num_atoms)) + str_to_float(*(buff+2)));
mesch_me_set(bigmat, i/(2*num_atoms), ((3*i)+2)%(3*num_atoms), value);

        }
      }
    else
      {
//what happened - why is there not 3 things on the line?
//matrix_2d_free(bigmat);
       }
//next line?
     }

//Symmetricalise? -> to make symmetric
for (i=0; i<(3*num_atoms); i++)
  {
  for (j=0; j<(3*num_atoms); j++)
    {
    value = mesch_me_get(bigmat, i, j) + mesch_me_get(bigmat, j, i);
    value *= 0.5;
    mesch_me_set(bigmat, i, j, value);
    mesch_me_set(bigmat, j, i, value);
    }
  }

correction_mat = mesch_mat_new(3*num_atoms,3);
mini_correction_zerooooo_mat = mesch_mat_new(3,3);

mesch_m_zero(correction_mat);
mesch_m_zero(mini_correction_zerooooo_mat);


//build the correction_matrix

for (i=0; i<mesch_rows_get(bigmat); i++)
  {
  for (j=0; j<mesch_cols_get(bigmat)/3; j++)
    {
//[3n][3] -> [i][0], [i][1], [i][2]

    value = mesch_me_get(bigmat, i, 3*j);
    mesch_me_add(correction_mat, i, 0, value);

    value = mesch_me_get(bigmat, i, (3*j)+1);
    mesch_me_add(correction_mat, i, 1, value);

    value = mesch_me_get(bigmat, i, (3*j)+2);
    mesch_me_add(correction_mat, i, 2, value);
    }

 //average each cell per row in the correction matrix
  value = 1.0 / (gdouble) num_atoms;

  mesch_me_mul(correction_mat, i, 0, value);
  mesch_me_mul(correction_mat, i, 1, value);
  mesch_me_mul(correction_mat, i, 2, value);
  }


//built mini matrix - [3][3]
for (i=0; i<mesch_rows_get(correction_mat); i++)
  {
  for (j=0; j<mesch_cols_get(correction_mat); j++)
    {
    value = mesch_me_get(correction_mat, i, j);
    mesch_me_add(mini_correction_zerooooo_mat, i%3, j, value);
    }
  }

//average the cells in mini_correction_zerooooo_mat

value = 1.0 / (gdouble) num_atoms;

for (i=0; i<mesch_rows_get(mini_correction_zerooooo_mat); i++)
  {
  for (j=0; j<mesch_cols_get(mini_correction_zerooooo_mat); j++)
    {
    mesch_me_mul(mini_correction_zerooooo_mat, i, j, value);
    }
  }

//zero point correction
//    crappy crappy fortran loop that i dont understand
//       do i=1,natoms
//         do j=1,natoms
//           do ii=1,3
//             do ij=1,3
//               correct = (zeroo(ii,ij)+zeroo(ij,ii))/2.0d0 -
//                         (zero(ii,ij,j)+zero(ij,ii,i))
//               do lx=-lxmax,lxmax
//               do ly=-lymax,lymax
//               do lz=-lzmax,lzmax
//                 phi(ii,i,ij,j,lx,ly,lz) = phi(ii,i,ij,j,lx,ly,lz) +
//                                           correct
//               enddo
//               enddo
//               enddo
//             enddo
//           enddo
//         enddo
//       enddo

        gdouble correction;
for (i=0; i<num_atoms; i++)
  {
  for (j=0; j<num_atoms; j++)
    {
    for (k=0; k<3; k++) //(ii)
      {
      for (l=0; l<3; l++) //(ij)
        {
// THIS WORKS - I HAVE TESTED IT - MANY TIMES......
        correction = mesch_me_get(mini_correction_zerooooo_mat, k, l);
        correction += mesch_me_get(mini_correction_zerooooo_mat, l, k);
        correction *= 0.5;
        correction -= mesch_me_get(correction_mat, 3*i+k, l);
        correction -= mesch_me_get(correction_mat, 3*j+l, k);

        mesch_me_add(bigmat, 3*i+k, 3*j+l, correction);
        }
      }
    }
  }

i=j=0;

for (list_i=model->cores ; list_i ; list_i=g_slist_next(list_i))
  {
  core_i = list_i->data;
  if (core_i->status & DELETED)
    {
    i++;
    continue;
    }

  wi = elements[core_i->atom_code].weight;
  g_assert(wi > G_MINDOUBLE);

  for (list_j=model->cores ; list_j ; list_j=g_slist_next(list_j))
    {
    core_j = list_j->data;
    if (core_j->status & DELETED)
      {
      j++;
      continue;
      }
//multi i rows.... 3 of them....

    wj = elements[core_j->atom_code].weight;
    g_assert(wj > G_MINDOUBLE);
    value  = 1.0 / sqrt(wi * wj);

    for (k=0; k<3; k++)
      {
      mesch_me_mul(bigmat, (3*i)+k, 3*j, value);
      mesch_me_mul(bigmat, (3*i)+k, (3*j)+1, value);
      mesch_me_mul(bigmat, (3*i)+k, (3*j)+2, value);
      }

    j++;
    }
  i++;
  j=0;
  }

model->siesta.eigen_xyz_atom_mat = mesch_mat_new(3*num_atoms, 3*num_atoms);
model->siesta.eigen_values = mesch_vec_new(3*num_atoms);

//external library call.

mesch_sev_compute(bigmat, model->siesta.eigen_xyz_atom_mat, model->siesta.eigen_values);

// stupid sort routine -> this is going to need a rewrite - its a bubble sort - O(n^2).
model->siesta.sorted_eig_values = g_malloc(sizeof(int[mesch_dim_get(model->siesta.eigen_values)]));

for (i=0; i<mesch_dim_get(model->siesta.eigen_values); i++)
  {
  model->siesta.sorted_eig_values[i] = i;
  }

gint temp_int;
gdouble freq_i, freq_ii;

gint sizeofeig = mesch_dim_get(model->siesta.eigen_values);

for (j=sizeofeig-1; j>1; j--)
  {
  for (i=0; i < j; i++)
    {
    freq_i = make_eigvec_freq(mesch_ve_get(model->siesta.eigen_values, model->siesta.sorted_eig_values[i]));

    freq_ii = make_eigvec_freq(mesch_ve_get(model->siesta.eigen_values, model->siesta.sorted_eig_values[i+1]));


    if (freq_i > freq_ii )
      {
      temp_int = model->siesta.sorted_eig_values[i];
      model->siesta.sorted_eig_values[i] = model->siesta.sorted_eig_values[i+1];
      model->siesta.sorted_eig_values[i+1] = temp_int;
      }
    }
  }

//PRINT METHOD FOR UWA VISIT
if (siestafileWRITE && matout)
  {
  fprintf(matout, "eig vectors 3N*3N\n-------------\n");
  for (j=0;j<(3*num_atoms); j++)
    {
    for (i=0; i<(3*num_atoms); i++)
      {
      fprintf(matout, "%f, ", mesch_me_get(bigmat, j, i));
      }
    fprintf(matout, "\n");
    }

  fprintf(matout, "\n\neig_vals\n-------------\n");
  for (i=0; i<mesch_dim_get(model->siesta.eigen_values); i++)
    {
    fprintf(matout, "%f, ", mesch_ve_get(model->siesta.eigen_values, i));
    }

  fprintf(matout, "\n\nfreqs\n-------------\n");
  for (i=0; i<mesch_dim_get(model->siesta.eigen_values); i++)
    {
    fprintf(matout, "%f, ", make_eigvec_freq(mesch_ve_get(model->siesta.eigen_values, i)));
    }

  fclose(matout);
  }

fclose(fp);

mesch_m_free(bigmat);
mesch_m_free(correction_mat);
mesch_m_free(mini_correction_zerooooo_mat);

model->siesta.current_animation = 0;

//Lookup using index array.
model->siesta.current_frequency = make_eigvec_freq(mesch_ve_get(model->siesta.eigen_values, model->siesta.sorted_eig_values[0]));
model->siesta.freq_disp_str = g_strdup_printf("%.2f", model->siesta.current_frequency);
model->siesta.num_animations = mesch_dim_get(model->siesta.eigen_values);
model->siesta.vibration_calc_complete = TRUE;

return(0);
}

/********************************/
/* SIESTA phonon display dialog */
/********************************/
void siesta_animation_dialog(GtkWidget *w, struct model_pak *model)
{
GList *list;
GtkWidget *window, *box, *hbox, *hbox2, *vbox, *label, *hscale, *entry, *button, *spin;
gpointer dialog;

g_assert(model != NULL);

/* don't recalculate modes if already done */
if (!model->siesta.vibration_calc_complete)
  if (siesta_phonon_calc(model))
    return;

/* request a dialog */
dialog = dialog_request(100, "Vibrational viewer", NULL, NULL, model);
if (!dialog)
  return;

window = dialog_window(dialog);
box = gtk_vbox_new(TRUE, 0);
gtk_container_add(GTK_CONTAINER(GTK_DIALOG(window)->vbox), box);

/* phonon selection */
vbox = gui_frame_vbox(NULL, FALSE, FALSE, box);
hbox = gtk_hbox_new(FALSE, 0);
gtk_box_pack_start(GTK_BOX(vbox), hbox, FALSE, FALSE, 0);

/* phonon frequency */
label = gtk_label_new("Frequency  ");
gtk_box_pack_start(GTK_BOX(hbox), label, FALSE, FALSE, 0);
entry = gtk_entry_new_with_max_length(LINELEN);
gtk_entry_set_text(GTK_ENTRY(entry), " ");
gtk_entry_set_editable(GTK_ENTRY(entry), FALSE);
gtk_box_pack_start(GTK_BOX(hbox), entry, FALSE, FALSE, 0);
dialog_child_set(dialog, "phonon_entry", entry);
gui_button(" < ", gui_siesta_mode_prev, dialog, hbox, FF);
gui_button(" > ", gui_siesta_mode_next, dialog, hbox, FF);

/* phonon slider */
hbox = gtk_hbox_new(FALSE, 0);
gtk_box_pack_start(GTK_BOX(vbox), hbox, FALSE, FALSE, 0);
hscale = gtk_hscale_new_with_range(1.0, model->siesta.num_animations, 1.0);
gtk_box_pack_start(GTK_BOX(hbox), hscale, TRUE, TRUE, 0);
dialog_child_set(dialog, "phonon_slider", hscale);
g_signal_connect(GTK_OBJECT(hscale), "value_changed",
                 GTK_SIGNAL_FUNC(gui_siesta_slider_changed), dialog);

/* phonon display options */
vbox = gui_frame_vbox(NULL, FALSE, FALSE, box);
hbox = gtk_hbox_new(FALSE, 0);
gtk_box_pack_start(GTK_BOX(vbox), hbox, FALSE, FALSE, 0);
button = new_check_button("Show eigenvectors", gui_siesta_mode_show, dialog, FALSE, hbox);
dialog_child_set(dialog, "phonon_toggle", button);

spin = new_spinner("Eigenvector scaling", 0.1, 9.9, 0.1, gui_siesta_mode_show, dialog, vbox);
gtk_spin_button_set_value(GTK_SPIN_BUTTON(spin), 4.0);
dialog_child_set(dialog, "phonon_scaling", spin);

spin = new_spinner("Animation resolution", 5.0, 50.0, 1.0, NULL, NULL, vbox);
dialog_child_set(dialog, "phonon_resolution", spin);

/* phonon mode animation */
vbox = gui_frame_vbox(NULL, FALSE, FALSE, box);

hbox = gtk_hbox_new(FALSE, 0);
gtk_box_pack_start(GTK_BOX(vbox), hbox, FALSE, FALSE, 0);
label = gtk_label_new("Animate mode ");
gtk_box_pack_start(GTK_BOX(hbox), label, FALSE, FALSE, 0);
hbox2 = gtk_hbox_new(FALSE, 0);
gtk_box_pack_end(GTK_BOX(hbox), hbox2, FALSE, FALSE, 0);
gui_icon_button("GDIS_PLAY", NULL, siesta_phonon_start, dialog, hbox2);
gui_icon_button("GDIS_STOP", NULL, siesta_phonon_stop, dialog, hbox2);

/* phonon mode movie generation */
hbox = gtk_hbox_new(FALSE, 0);
gtk_box_pack_start(GTK_BOX(vbox), hbox, FALSE, FALSE, 0);
label = gtk_label_new("Create movie ");
gtk_box_pack_start(GTK_BOX(hbox), label, FALSE, FALSE, 0);
entry = gtk_entry_new_with_max_length(LINELEN);
gtk_entry_set_text(GTK_ENTRY(entry), "movie_name");
gtk_entry_set_editable(GTK_ENTRY(entry), TRUE);
gtk_box_pack_start(GTK_BOX(hbox), entry, FALSE, FALSE, 0);
dialog_child_set(dialog, "phonon_movie_name", entry);

/* movie type */
list = NULL;
list = g_list_append(list, "gif");
list = g_list_append(list, "mpg");
entry = gui_pulldown_new(NULL, list, FALSE, hbox);
dialog_child_set(dialog, "phonon_movie_type", entry);
gui_button_x(NULL, siesta_phonon_movie, dialog, hbox);

/* init and display */
gui_siesta_frequency_set(dialog);

gtk_widget_show_all(window);
}

/*************************/
/* main SIESTA interface */
/*************************/
void gui_siesta_dialog(void)
{
siestafileWRITE = FALSE;
GtkWidget *window;
gpointer dialog;
struct model_pak *model;

model = sysenv.active_model;

if (!model)
  return;

/* request a dialog */
dialog = dialog_request(200, "Siesta Setup", NULL, NULL, model);

if (!dialog)
    return;

window = dialog_window(dialog);

file_handler_page_generator(dialog, window, model);

siesta_gui_page_generator(dialog, window, model);

gui_button("Animation", siesta_animation_dialog, model, GTK_DIALOG(window)->action_area, TT);
gui_button("FileGenerate", siesta_file_dialog, model, GTK_DIALOG(window)->action_area, TT);
gui_stock_button(GTK_STOCK_CLOSE, dialog_destroy, dialog, GTK_DIALOG(window)->action_area);

gtk_widget_show_all(window);
}

void file_handler_page_generator(gpointer * dialog, GtkWidget * window, struct model_pak *model)
{
    GtkWidget *frame, *vbox;

    vbox = gtk_vbox_new(FALSE,0);
    gtk_box_pack_start(GTK_BOX(GTK_DIALOG(window)->vbox), vbox, TRUE, TRUE, 10);

    gui_direct_check("Write csv", &siestafileWRITE, random_action, NULL, vbox);

    frame = gtk_frame_new(" File Specifics " );
    gtk_box_pack_start(GTK_BOX(vbox), frame, FALSE, FALSE, 0);
    vbox = gtk_vbox_new(FALSE,0);
    gtk_container_add(GTK_CONTAINER(frame), vbox);

     gui_text_entry("Filename ", &model->siesta.modelfilename, TRUE, TRUE, vbox);
}

void siesta_gui_page_generator(gpointer * dialog, GtkWidget * window, struct model_pak *model)
{
    GtkWidget *frame, *radio_vbox, *vbox, *vbox1, *hbox, *label, *notebook;
    GtkWidget *geomhbox, *geomhbox2, *geomvbox, *geombook, *combo;
    GtkWidget *page, *geompage, *button, *vbox2;
    //more sensitive boxes
    GtkWidget *target_pressure_sens_box, *zeta_sensitive_box, *sensitive_box;

    GList *list;

    /* create notebook */
    hbox = gtk_hbox_new(FALSE,0);
    gtk_box_pack_start(GTK_BOX(GTK_DIALOG(window)->vbox), hbox, TRUE, TRUE, 10);
    vbox = gtk_vbox_new(FALSE, 0);
    gtk_box_pack_start (GTK_BOX(hbox), vbox, FALSE, FALSE, 0);


    // Specifics Page
    frame = gtk_frame_new(" System Params " );
    gtk_box_pack_start(GTK_BOX(hbox), frame, FALSE, FALSE, 0);
    vbox = gtk_vbox_new(FALSE, 0);
    gtk_container_add(GTK_CONTAINER(frame), vbox);


    hbox = gtk_hbox_new (FALSE, 0);
    gtk_box_pack_start(GTK_BOX(vbox),hbox,TRUE,TRUE,0);
    notebook = gtk_notebook_new();
    gtk_widget_show (notebook);
    gtk_container_add(GTK_CONTAINER(hbox), notebook);

    page = gtk_vbox_new(FALSE, 0);
    label = gtk_label_new (" Electronic Structure ");
    gtk_notebook_append_page(GTK_NOTEBOOK(notebook), page, label);


    /* title display */
    frame = gtk_frame_new(" Basis Set ");
    gtk_box_pack_start(GTK_BOX(page), frame, FALSE, FALSE, 0);
    radio_vbox = gtk_vbox_new(FALSE, 0);
    gtk_container_add(GTK_CONTAINER(frame), radio_vbox);

    /* do the radio buttons */
    new_radio_group(0, radio_vbox, TT);
    button = add_radio_button("Single zeta", (gpointer) set_basis_set_sz, model);
    if (model->siesta.basis_set == SZ_ZETA)
        gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(button), TRUE);
    button = add_radio_button("Double zeta", (gpointer) set_basis_set_dz, model);
    if (model->siesta.basis_set == DZ_ZETA)
        gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(button), TRUE);
    button = add_radio_button("Single zeta polarised", (gpointer) set_basis_set_szp, model);
    if (model->siesta.basis_set == SZP_ZETA)
        gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(button), TRUE);
    button = add_radio_button("Double zeta polarised", (gpointer) set_basis_set_dzp, model);
    if (model->siesta.basis_set == DZP_ZETA)
        gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(button), TRUE);

    /*
     * Redundant code.
    button = add_radio_button("Triple zeta", (gpointer) set_basis_set_tz, model);
    if (model->siesta.basis_set == TZ)
        gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(button), TRUE);
    button = add_radio_button("Triple zeta polarised", (gpointer) set_basis_set_tzp, model);
    if (model->siesta.basis_set == TZP)
        gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(button), TRUE);
    */

    button = add_radio_button("Custom Zeta", (gpointer) set_basis_set_custom, model);


    // Set sensitive frame? how do i do that?
    frame = gtk_frame_new(" Custom Zeta Config ");
    gtk_box_pack_start(GTK_BOX(radio_vbox), frame, FALSE, FALSE, 0);
    zeta_sensitive_box = gtk_vbox_new(FALSE, 0);

    model->siesta.custom_zeta_frame = zeta_sensitive_box;

    gtk_container_add(GTK_CONTAINER(frame), zeta_sensitive_box);

    gui_direct_spin("Zeta level ",
            &model->siesta.custom_zeta, 1, 5, 1,
            zeta_warning, model, zeta_sensitive_box);
    gui_direct_spin("Polarisation level ",
            &model->siesta.custom_zeta_polarisation, 0, 5, 1,
            zeta_warning, model, zeta_sensitive_box);

    if (model->siesta.basis_set != CUSTOM_ZETA)
    {
        gtk_widget_set_sensitive(GTK_WIDGET(model->siesta.custom_zeta_frame), FALSE);
    }
    else
    {
        gtk_widget_set_sensitive(GTK_WIDGET(model->siesta.custom_zeta_frame), TRUE);
        gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(button), TRUE);
    }


    frame = gtk_frame_new("Special inputs");
    gtk_box_pack_start(GTK_BOX(page), frame, FALSE, FALSE, 0);
    vbox = gtk_vbox_new(FALSE, 0);
    gtk_container_add(GTK_CONTAINER(frame), vbox);

    gui_direct_spin("Split Zeta norm ",
            &model->siesta.split_zeta_norm, 0.01, 1.00, 0.01,
            dud_action, model, vbox);
    gui_direct_spin("Energy Shift (Ryd) ",
            &model->siesta.energy_shift, 0.001, 0.05, 0.001,
            dud_action, model, vbox);
    gui_direct_spin("Mesh Cutoff (Ryd) ",
            &model->siesta.mesh_cutoff, 40.0, 1000.0, 5,
            dud_action, model, vbox);
    gui_direct_spin("Electronic temp (K) ",
            &model->siesta.electronic_temperature, 0.0,
            500.0, 1.0, dud_action, model,
            vbox);
    gui_direct_check("Spin Polarised",
             &model->siesta.spin_polarised,
             random_action, page, vbox);

/*
 *  Periodic model checking....
 */
    if (model->periodic == TRUE)
    {
    gui_direct_check("Is it Periodic?",
             &model->siesta.is_periodic,
             kgrid_action, page, vbox);

    frame = gtk_frame_new("Periodic");
    gtk_box_pack_start(GTK_BOX(page), frame, FALSE, FALSE, 0);
    vbox1 = gtk_vbox_new(FALSE, 0);
    gtk_container_add(GTK_CONTAINER(frame), vbox1);
    gui_direct_spin("kgrid Cutoff ",
        &model->siesta.kgrid_cutoff, 1.0,
        50.0, 1.0, dud_action, model, vbox1);
    }


    //Next Page

    page = gtk_vbox_new(FALSE, 0);
    label = gtk_label_new (" SCF ");
    gtk_notebook_append_page(GTK_NOTEBOOK(notebook), page, label);

    frame = gtk_frame_new("Method");
    gtk_box_pack_start(GTK_BOX(page), frame, FALSE, FALSE, 0);
    vbox = gtk_vbox_new(FALSE, 0);
    gtk_container_add(GTK_CONTAINER(frame), vbox);

    //Common to both
    gui_direct_spin("Number of Cycles ",
            &model->siesta.no_of_cycles, 1.0, 100.0, 1.0,
            dud_action, model, vbox);
    gui_direct_spin("Mixing Weight ",
            &model->siesta.mixing_weight, 0.01, 0.5, 0.01,
            dud_action, model, vbox);


    sensitive_box = gtk_vbox_new(TRUE, 0);
    gui_direct_check("Pulay Mixing", &model->siesta.pulay_mixing,
                    set_pulay_sensitive, sensitive_box, vbox);
    gtk_box_pack_end(GTK_BOX(vbox), sensitive_box, TRUE, TRUE, 0);

    /* sensitive box */
    vbox = gtk_vbox_new(FALSE,0);
    gtk_box_pack_start(GTK_BOX(sensitive_box), vbox, TRUE, TRUE, 0);
    gtk_container_set_border_width(GTK_CONTAINER(vbox), PANEL_SPACING);

    gui_direct_spin("Number of Pulay Matrices ",
            &model->siesta.no_of_pulay_matrices,
            1.0, 1000.0, 5, dud_action, model,
            vbox);

    //inital call to set sensitive.
    set_pulay_sensitive(vbox, sensitive_box);

    frame = gtk_frame_new("Speed Hacks?");
    gtk_box_pack_start(GTK_BOX(page), frame, FALSE, FALSE, 0);
    vbox = gtk_vbox_new(FALSE, 0);
    gtk_container_add(GTK_CONTAINER(frame), vbox);

    gui_direct_check("Divide and Conquer", &model->siesta.diag_divide_and_conquer,
                    random_action, vbox, vbox);

    //Next Page

    page = gtk_vbox_new(FALSE, 0);
    label = gtk_label_new (" Geometry ");
    gtk_notebook_append_page(GTK_NOTEBOOK(notebook), page, label);

    frame = gtk_frame_new("Method");
    gtk_box_pack_start(GTK_BOX(page), frame, FALSE, FALSE, 0);
    vbox = gtk_vbox_new(FALSE, 0);
    gtk_container_add(GTK_CONTAINER(frame), vbox);


    /* do the radio buttons */
    new_radio_group(0, vbox, TT);
    button = add_radio_button("Single Point", (gpointer) set_geom_runtype_sp, model);
    if (model->siesta.run_type == SINGLE_POINT)
        gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(button), TRUE);
    button = add_radio_button("Optimisation", (gpointer) set_geom_runtype_opt, model);
    if (model->siesta.run_type == OPTIMISATION)
        gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(button), TRUE);
    button = add_radio_button("Molecular Dynamics", (gpointer) set_geom_runtype_md, model);
    if (model->siesta.run_type == MOLECULAR_DYNAMICS)
        gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(button), TRUE);
    button = add_radio_button("Phonon Calculation", (gpointer) set_geom_runtype_pc, model);
    if (model->siesta.run_type == PHONON_CALCULATION)
        gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(button), TRUE);

//options frame

    frame = gtk_frame_new("Options");
    gtk_box_pack_start(GTK_BOX(page), frame, FALSE, FALSE, 0);
    vbox = gtk_vbox_new(FALSE, 0);
    gtk_container_add(GTK_CONTAINER(frame), vbox);

    geombook = gtk_notebook_new();
    gtk_container_add(GTK_CONTAINER(vbox), geombook);

    //HACK HACK HACK
    model->siesta.geom_notebook_hack = geombook;


//  Single Point Page
    geompage = gtk_vbox_new(FALSE, 0);
    label = gtk_label_new (" Single Point ");
    gtk_notebook_append_page(GTK_NOTEBOOK(geombook), geompage, label);
//  Single Point Options

    label = gtk_label_new("Number of Steps - locked to Zero\n");
    gtk_box_pack_start(GTK_BOX(geompage), label, FALSE, FALSE, 0);
    label = gtk_label_new("(single point)\n");
    gtk_box_pack_start(GTK_BOX(geompage), label, FALSE, FALSE, 0);


    geompage = gtk_vbox_new(FALSE, 0);
    label = gtk_label_new (" Optimisation ");
    gtk_notebook_append_page(GTK_NOTEBOOK(geombook), geompage, label);

    vbox = gtk_vbox_new(FALSE,0);
    gtk_box_pack_start(GTK_BOX(geompage), vbox, FALSE, FALSE, 0);

    gui_direct_check("Optimise cell",
            &model->siesta.md_variable_cell,
            optimise_cell_action, model, vbox);

    gui_direct_spin("Number of Steps ",
            &model->siesta.number_of_steps, 2.0, 1000.0, 1.0,
            dud_action, NULL, vbox);

    //Target pressure - only if optimise cell
    //meant to be sensitive
    target_pressure_sens_box = gtk_vbox_new(FALSE,0);
    gtk_box_pack_start(GTK_BOX(geompage), target_pressure_sens_box, FALSE, FALSE, 0);

    gui_direct_spin("Target Pressure",
            &model->siesta.md_target_pressure, -5.0, 5.0, 0.1,
            dud_action, NULL, vbox);

    //Target stress

    frame = gtk_frame_new("Stress Tensors");
    gtk_box_pack_start(GTK_BOX(geompage), frame, FALSE, FALSE, 0);

    //int labled stress vectors.
    geomvbox = gtk_hbox_new(FALSE, 0);
    gtk_container_add(GTK_CONTAINER(frame), geomvbox);

    geomhbox = gtk_vbox_new(FALSE, 0);
    gtk_container_add(GTK_CONTAINER(geomvbox), geomhbox);
    geomhbox2 = gtk_vbox_new(FALSE, 0);
    gtk_container_add(GTK_CONTAINER(geomvbox), geomhbox2);

    gui_direct_spin("   xx =",
            &model->siesta.md_target_stress_xx, -5.0, 5.0, 0.1,
            dud_action, model, geomhbox);
    gui_direct_spin("   yy =",
            &model->siesta.md_target_stress_yy, -5.0, 5.0, 0.1,
            dud_action, model, geomhbox);
    gui_direct_spin("   zz =",
            &model->siesta.md_target_stress_zz, -5.0, 5.0, 0.1,
            dud_action, model, geomhbox);
    gui_direct_spin("   xy =",
            &model->siesta.md_target_stress_xy, -5.0, 5.0, 0.1,
            dud_action, model, geomhbox2);
    gui_direct_spin("   xz =",
            &model->siesta.md_target_stress_xz, -5.0, 5.0, 0.1,
            dud_action, model, geomhbox2);
    gui_direct_spin("   yz =",
            &model->siesta.md_target_stress_yz, -5.0, 5.0, 0.1,
            dud_action, model, geomhbox2);


    frame = gtk_frame_new("Termination Options");
    gtk_box_pack_start(GTK_BOX(geompage), frame, FALSE, FALSE, 0);
    vbox = gtk_vbox_new(FALSE, 0);
    gtk_container_add(GTK_CONTAINER(frame), vbox);

    gui_direct_spin("Max CG displacement",
            &model->siesta.md_max_cg_displacement, 0.0, 2.0, 0.01,
            dud_action, NULL, vbox);
    gui_direct_spin("Max Force Tolerance",
            &model->siesta.md_max_force_tol, 0.0, 2.0, 0.01,
            dud_action, NULL, vbox);
    gui_direct_spin("Max Stress Tolerance",
            &model->siesta.md_max_stress_tol, 0.0, 2.0, 0.1,
            dud_action, NULL, vbox);

    frame = gtk_frame_new("Job Options");
    gtk_box_pack_start(GTK_BOX(geompage), frame, FALSE, FALSE, 0);
    vbox = gtk_vbox_new(FALSE, 0);
    gtk_container_add(GTK_CONTAINER(frame), vbox);

    gui_direct_check("Restart (Use saved data)",
            &model->siesta.use_saved_data,
            random_action, page, vbox);



    geompage = gtk_vbox_new(FALSE, 0);
    label = gtk_label_new (" Molecular Dynamics ");
    gtk_notebook_append_page(GTK_NOTEBOOK(geombook), geompage, label);

    frame = gtk_frame_new("Run Type");
    gtk_box_pack_start(GTK_BOX(geompage), frame, FALSE, FALSE, 0);
    vbox = gtk_vbox_new(FALSE, 0);
    gtk_container_add(GTK_CONTAINER(frame), vbox);

    list = NULL;
    list = g_list_append(list, "Verlet");
    list = g_list_append(list, "Nose");
    list = g_list_append(list, "Parrinello-Rahman");
    list = g_list_append(list, "Nose-Parrinello-Rahman");
    list = g_list_append(list, "Anneal");
    list = g_list_append(list, "Phonon");

    combo = gtk_combo_new();
    gtk_entry_set_editable(GTK_ENTRY(GTK_COMBO(combo)->entry), FALSE);
    gtk_combo_set_popdown_strings(GTK_COMBO(combo), list);
    gtk_box_pack_end(GTK_BOX(vbox), combo, FALSE, FALSE, 0);
    g_signal_connect(GTK_OBJECT(GTK_COMBO(combo)->entry), "changed",
                    GTK_SIGNAL_FUNC(set_md_run_type), (gpointer *) combo);

    frame = gtk_frame_new("Temperature");
    gtk_box_pack_start(GTK_BOX(geompage), frame, FALSE, FALSE, 0);
    vbox = gtk_vbox_new(FALSE, 0);
    gtk_container_add(GTK_CONTAINER(frame), vbox);

    gui_direct_spin("Inital Temperature ",
            &model->siesta.md_inital_temperature, 0.0, 500.0, 1.0,
            dud_action, NULL, vbox);
    gui_direct_spin("Target Temperature ",
            &model->siesta.md_target_temperature, 0.0, 500.0, 1.0,
            dud_action, NULL, vbox);

    frame = gtk_frame_new("Pressure");
    gtk_box_pack_start(GTK_BOX(geompage), frame, FALSE, FALSE, 0);
    vbox = gtk_vbox_new(FALSE, 0);
    gtk_container_add(GTK_CONTAINER(frame), vbox);

    gui_direct_spin("Target Pressure ",
            &model->siesta.pressure, -10.0, 10.0, 0.01,
            dud_action, model, vbox);



    frame = gtk_frame_new("Time");
    gtk_box_pack_start(GTK_BOX(geompage), frame, FALSE, FALSE, 0);
    vbox = gtk_vbox_new(FALSE, 0);
    gtk_container_add(GTK_CONTAINER(frame), vbox);

    gui_direct_spin("Initial Timestep ",
            &model->siesta.md_inital_time_step, 0.0, 500.0, 1.0,
            dud_action, NULL, vbox);
    gui_direct_spin("Final Timestep ",
            &model->siesta.md_final_time_step, 0.0, 500.0, 1.0,
            dud_action, NULL, vbox);
    gui_direct_spin("Length of Timestep ",
            &model->siesta.md_length_time_step, 0.1, 200.0, 0.1,
            dud_action, NULL, vbox);

    frame = gtk_frame_new("Job Options");
    gtk_box_pack_start(GTK_BOX(geompage), frame, FALSE, FALSE, 0);
    vbox = gtk_vbox_new(FALSE, 0);
    gtk_container_add(GTK_CONTAINER(frame), vbox);

    gui_direct_check("Restart (Use saved data)",
            &model->siesta.use_saved_data,
            random_action, page, vbox);


    geompage = gtk_vbox_new(FALSE, 0);
    label = gtk_label_new (" Phonon ");
    gtk_notebook_append_page(GTK_NOTEBOOK(geombook), geompage, label);

    frame = gtk_frame_new("Differencing");
    gtk_box_pack_start(GTK_BOX(geompage), frame, FALSE, FALSE, 0);
    vbox = gtk_vbox_new(FALSE, 0);
    gtk_container_add(GTK_CONTAINER(frame), vbox);

    gui_direct_spin("Finite Difference step size ",
            &model->siesta.finite_diff_step_size, 0.1, 10.0, 0.1,
            dud_action, NULL, vbox);

    frame = gtk_frame_new("Job Options");
    gtk_box_pack_start(GTK_BOX(geompage), frame, FALSE, FALSE, 0);
    vbox = gtk_vbox_new(FALSE, 0);
    gtk_container_add(GTK_CONTAINER(frame), vbox);

    gui_direct_check("Restart (Use saved data)",
            &model->siesta.use_saved_data,
            random_action, page, vbox);

    gtk_notebook_set_current_page (GTK_NOTEBOOK(geombook), 3); // BROKEN
    gtk_notebook_set_show_tabs (GTK_NOTEBOOK(geombook), FALSE);
    gtk_notebook_set_show_border (GTK_NOTEBOOK(geombook), FALSE);


    //Next Page

    page = gtk_vbox_new(FALSE, 0);
    label = gtk_label_new (" File I/O ");
    gtk_notebook_append_page(GTK_NOTEBOOK(notebook), page, label);

    frame = gtk_frame_new("Files");
    gtk_box_pack_start(GTK_BOX(page), frame, FALSE, FALSE, 0);
    vbox = gtk_vbox_new(FALSE, 0);
    gtk_container_add(GTK_CONTAINER(frame), vbox);


    gui_direct_check("Long output ",
            &model->siesta.long_output,
            long_output_click, model, vbox);


    frame = gtk_frame_new("Mesh potential");
    gtk_box_pack_start(GTK_BOX(page), frame, FALSE, FALSE, 0);
    vbox = gtk_vbox_new(FALSE, 0);
    gtk_container_add(GTK_CONTAINER(frame), vbox);

    gui_direct_spin("Density of states ",
            &model->siesta.density_of_states, 0.1, 10.0, 0.1,
            dud_action, model, vbox);
    gui_direct_spin("Density on mesh ",
            &model->siesta.density_on_mesh, 0.1, 10.0, 0.1,
            dud_action, model, vbox);
    gui_direct_spin("Electrostatic pot on mesh ",
            &model->siesta.electrostatic_pot_on_mesh, 0.1, 10.0, 0.1,
            dud_action, model, vbox);

    frame = gtk_frame_new("Output Options");
    gtk_box_pack_start(GTK_BOX(page), frame, FALSE, FALSE, 0);

    hbox = gtk_hbox_new(FALSE, 0);
    gtk_container_add(GTK_CONTAINER(frame), hbox);


    vbox = gtk_vbox_new(FALSE, 0);
    vbox2 = gtk_vbox_new(FALSE, 0);
    gtk_box_pack_start(GTK_BOX(hbox), vbox, FALSE, FALSE, 0);
    gtk_box_pack_start(GTK_BOX(hbox), vbox2, FALSE, FALSE, 0);

    model->siesta.long_output_widget = hbox;

    gui_direct_check("WriteCoorStep ",
            &model->siesta.file_output_write_coor_step,
            random_action, page, vbox);
    gui_direct_check("WriteForces ",
            &model->siesta.file_output_write_forces,
            random_action, page, vbox);
    gui_direct_check("WriteKpoints ",
            &model->siesta.file_output_write_kpoints,
            random_action, page, vbox);
    gui_direct_check("WriteEigenvalues ",
            &model->siesta.file_output_write_eigenvalues,
            random_action, page, vbox);
    gui_direct_check("WriteKbands ",
            &model->siesta.file_output_write_kbands,
            random_action, page, vbox2);
    gui_direct_check("WriteBands ",
            &model->siesta.file_output_write_bands,
            random_action, page, vbox2);
    gui_direct_check("WriteWaveFunctions ",
            &model->siesta.file_output_write_wavefunctions,
            random_action, page, vbox2);

    gui_direct_spin("WriteMullikenPop",
            &model->siesta.file_output_write_mullikenpop, 0.0, 3.0, 1.0,
            dud_action, page, vbox2);

    frame = gtk_frame_new("Extra Output Options");
    gtk_box_pack_start(GTK_BOX(page), frame, FALSE, FALSE, 0);
    vbox = gtk_vbox_new(FALSE, 0);
    gtk_container_add(GTK_CONTAINER(frame), vbox);

    gui_direct_check("Write density matrix",
            &model->siesta.file_output_write_dm,
            random_action, page, vbox);
    gui_direct_check("Write Xmol coordinates",
            &model->siesta.file_output_write_coor_xmol,
            random_action, page, vbox);
    gui_direct_check("Write cerius coordinates",
            &model->siesta.file_output_write_coor_cerius,
            random_action, page, vbox);
    gui_direct_check("Write MD xmol",
            &model->siesta.file_output_write_md_xmol,
            random_action, page, vbox);
    gui_direct_check("Write MD history",
            &model->siesta.file_output_write_md_history,
            random_action, page, vbox);

    gui_relation_update(NULL);
    g_list_free(list);
}

/*-------------------
 * Modifiers for the radio buttons
 * Electronic Structure PAGE
 *-------------------
 */
void set_basis_set_sz(struct model_pak * model)
{
    model->siesta.basis_set = SZ_ZETA;
    gtk_widget_set_sensitive(GTK_WIDGET(model->siesta.custom_zeta_frame), FALSE);
}
void set_basis_set_dz(struct model_pak * model)
{
    model->siesta.basis_set = DZ_ZETA;
    gtk_widget_set_sensitive(GTK_WIDGET(model->siesta.custom_zeta_frame), FALSE);
}
void set_basis_set_szp(struct model_pak * model)
{
    model->siesta.basis_set = SZP_ZETA;
    gtk_widget_set_sensitive(GTK_WIDGET(model->siesta.custom_zeta_frame), FALSE);
}
void set_basis_set_dzp(struct model_pak * model)
{
    model->siesta.basis_set = DZP_ZETA;
    gtk_widget_set_sensitive(GTK_WIDGET(model->siesta.custom_zeta_frame), FALSE);
}

void set_basis_set_custom(struct model_pak * model)
{
    model->siesta.basis_set = CUSTOM_ZETA;
    gtk_widget_set_sensitive(GTK_WIDGET(model->siesta.custom_zeta_frame), TRUE);
}


/*-----------------------------
 * More radio buttons - GEOM PAGE
 */
void set_geom_runtype_sp(struct model_pak * data)
{
    data->siesta.run_type = SINGLE_POINT;
    //change the notebook - HACK HACK HACK
    gtk_notebook_set_current_page(GTK_NOTEBOOK(data->siesta.geom_notebook_hack), 0);
}
void set_geom_runtype_opt(struct model_pak * data)
{
    data->siesta.run_type = OPTIMISATION;
    //change the notebook - HACK HACK HACK
    gtk_notebook_set_current_page(GTK_NOTEBOOK(data->siesta.geom_notebook_hack), 1);

}
void set_geom_runtype_md(struct model_pak * data)
{
    data->siesta.run_type = MOLECULAR_DYNAMICS;
    //change the notebook - HACK HACK HACK
    gtk_notebook_set_current_page( GTK_NOTEBOOK(data->siesta.geom_notebook_hack), 2);
}
void set_geom_runtype_pc(struct model_pak * data)
{
    data->siesta.run_type = PHONON_CALCULATION;
    data->siesta.md_type_of_run = FC_MDRUN;
    //change the notebook - HACK HACK HACK
    gtk_notebook_set_current_page( GTK_NOTEBOOK(data->siesta.geom_notebook_hack), 3);
}

void set_geom_constantcomp_cp(struct model_pak * data)
{
    data->siesta.constant_component = CONSTANT_PRESSURE;
}
void set_geom_constantcomp_cv(struct model_pak * data)
{
    data->siesta.constant_component = CONSTANT_VOLUME;
}


void kgrid_action (GtkWidget * mywidget, gboolean checkit)
{
    if (checkit == TRUE)
    {
        gtk_widget_set_sensitive(mywidget, TRUE);
    }
    else
    {
        gtk_widget_set_sensitive(mywidget, FALSE);
    }
}

void dud_action(GtkWidget *w, struct model_pak * data)
{
//Do i need this?
//maybe i dooo
}

void set_pulay_sensitive(GtkWidget *w, GtkWidget *a_frame)
{
struct model_pak * model;
model = sysenv.active_model;

if (model->siesta.pulay_mixing)
  gtk_widget_set_sensitive(GTK_WIDGET(a_frame), TRUE);
else
  gtk_widget_set_sensitive(GTK_WIDGET(a_frame), FALSE);
}


void random_action(GtkWidget *w, GtkWidget *box)
{
/* do nothing
 */
}

void optimise_cell_action(GtkWidget *w, struct model_pak * model)
{
    // Must set it to constant pressure if volume is changing?
}


void long_output_click(GtkWidget *w, struct model_pak * model)
{
    if (model->siesta.long_output)
    {
        model->siesta.file_output_write_coor_step = TRUE;
        model->siesta.file_output_write_forces = TRUE;
        model->siesta.file_output_write_kpoints = TRUE;
        model->siesta.file_output_write_eigenvalues = TRUE;
        model->siesta.file_output_write_bands = TRUE;
        model->siesta.file_output_write_kbands = TRUE;
        model->siesta.file_output_write_wavefunctions = TRUE;
        model->siesta.file_output_write_mullikenpop = 1;
        gui_relation_update(NULL);
    }
    else
    {
        model->siesta.file_output_write_coor_step = FALSE;
        model->siesta.file_output_write_forces = FALSE;
        model->siesta.file_output_write_kpoints = FALSE;
        model->siesta.file_output_write_eigenvalues = FALSE;
        model->siesta.file_output_write_bands = FALSE;
        model->siesta.file_output_write_kbands = FALSE;
        model->siesta.file_output_write_wavefunctions = FALSE;
        model->siesta.file_output_write_mullikenpop = 0;
        gui_relation_update(NULL);
    }
}



void zeta_warning (GtkWidget *w, struct model_pak * model)
{
    if ((int) model->siesta.custom_zeta > 4)
    {
        gchar * message;

        message = g_strdup("What are you doing?\nMore than 4 zeta levels\nGrind Grind Grind?");
        gui_text_show(ERROR, message);
        g_free(message);
    }
    if ((int) model->siesta.custom_zeta_polarisation > 3)
    {

        gchar * message;
        message = g_strdup("What are you doing?\nMore than 3 polarisations?");
        gui_text_show(ERROR, message);
        g_free(message);
    }
}

/********************************/
/* SIESTA phonon display dialog */
/********************************/

void set_md_run_type(GtkWidget *w, gpointer *data)
{
    //evil active model call - evil evil evil
    struct model_pak *model;
    model = sysenv.active_model;

    gchar *buff;

    buff = g_strdup(gtk_entry_get_text(GTK_ENTRY(GTK_COMBO(data)->entry)));
    if (g_ascii_strncasecmp("Verlet", buff, 6) == 0)
        model->siesta.md_type_of_run = VERLET_MDRUN;
    else if (g_ascii_strncasecmp("Nose", buff, 4) == 0)
        model->siesta.md_type_of_run = NOSE_MDRUN;
    else if (g_ascii_strncasecmp("Parrinello-Rahman", buff, 17) == 0)
        model->siesta.md_type_of_run = PARRINELLOPAHMAN_MDRUN;
    else if (g_ascii_strncasecmp("Nose-Parrinello-Rahman", buff, 22) == 0)
        model->siesta.md_type_of_run = NOSEPARRINELLOPAHMAN_MDRUN;
    else if (g_ascii_strncasecmp("Anneal", buff, 6) == 0)
        model->siesta.md_type_of_run = ANNEAL_MDRUN;


    g_free(buff);
}

gdouble make_eigvec_freq(gdouble value)
{
    gdouble xmagic;
    xmagic = 519.6;
    if (value > 0.0)
    {
        return sqrt(xmagic*xmagic*value);
    }
    else    // fake number!
    {
        return -sqrt(-xmagic*xmagic*value);
    }
}

void siesta_file_dialog(GtkWidget *w, struct model_pak * model)
{
    GtkWidget *hbox, *vbox, *vbox2, *window;
    GtkWidget *dialog, *frame, *label;
    /* Create the widgets */

    /* request a dialog */
    dialog = dialog_request(100, "Multi File Writer", NULL, NULL, model);

    if (!dialog)
        return;

    window = dialog_window(dialog);

    vbox = gtk_vbox_new(TRUE,5);
    gtk_container_add (GTK_CONTAINER (GTK_DIALOG(window)->vbox), vbox);

    frame = gtk_frame_new(" Atoms in system ");
    gtk_box_pack_start(GTK_BOX(vbox), frame, TRUE, FALSE, 0);
    hbox = gtk_vbox_new(TRUE, 0);
    gtk_container_add(GTK_CONTAINER(frame), hbox);

    label = gtk_label_new(g_strdup_printf("Atoms seen by GDIS = %d\n", model->num_atoms));
    gtk_box_pack_start(GTK_BOX(hbox), label, TRUE, FALSE, 0);

    //how to round it?
    model->siesta.atoms_per_job = (gdouble) ((int) 100/(model->num_atoms));
    label = gtk_label_new(g_strdup_printf("Recommended atoms per job = %d\n", (int) 100/(model->num_atoms)));
    gtk_box_pack_start(GTK_BOX(hbox), label, TRUE, FALSE, 0);

    frame = gtk_frame_new(" Job Totals ");
    gtk_box_pack_start(GTK_BOX(vbox), frame, TRUE, FALSE, 0);
    hbox = gtk_hbox_new(TRUE, 0);
    gtk_container_add(GTK_CONTAINER(frame), hbox);

    vbox2 = gtk_vbox_new(TRUE, 0);
    gtk_box_pack_start(GTK_BOX(hbox), vbox2, TRUE, FALSE, 0);

    gui_direct_spin("Atoms per Job   ", &model->siesta.atoms_per_job, 1.0, model->num_atoms, 1.0, dud_action, model, vbox2);

    frame = gtk_frame_new(" Job Setup ");
    gtk_box_pack_start(GTK_BOX(vbox), frame, TRUE, FALSE, 0);
    hbox = gtk_hbox_new(TRUE, 0);
    gtk_container_add(GTK_CONTAINER(frame), hbox);


    vbox = gtk_vbox_new(TRUE, 0);
    gtk_box_pack_start(GTK_BOX(hbox), vbox, TRUE, FALSE, 0);

    // file location

    gui_stock_button(GTK_STOCK_SAVE, siesta_file_save_loop, model, GTK_DIALOG(window)->action_area);
    gui_button("Save and Queue", siesta_save_n_queue, model, GTK_DIALOG(window)->action_area, TT);
    gui_stock_button(GTK_STOCK_CLOSE, dialog_destroy, dialog, GTK_DIALOG(window)->action_area);

    gtk_widget_show_all (window);
}

void siesta_file_save_loop(GtkWidget *w, struct model_pak * model)
{
    int i;
    int total_files;
    gchar * filename, * orig_basename;

    total_files = model->num_atoms / model->siesta.atoms_per_job +1;

    orig_basename = g_strdup(model->basename);

    for (i=0; i<total_files; i++)
    {
        //set the range for md_fc_first, and md_fc_last
        model->siesta.md_fc_first = i* model->siesta.atoms_per_job + 1;
        model->siesta.md_fc_last = i * model->siesta.atoms_per_job + model->siesta.atoms_per_job;
        filename = g_strdup_printf("%s.%d.%d", orig_basename, i+1, total_files);
        if (model->siesta.md_fc_last >= model->num_atoms)
        {
            i++;
            model->siesta.md_fc_last = model->num_atoms;
        }
        model->basename = g_strdup(filename);

        filename = g_strdup_printf("%s.fdf", filename);

        // call the fdf save bit....
        file_save_as(filename, model);

        //grab current grid_config_pak
        struct grid_config_pak * grid_pak=NULL;

        //get job directory
        siesta_make_runscript(model->basename, sysenv.cwd, grid_pak);

        gui_text_show(OUTPUT, g_strdup_printf("fc-first:\t%d\tfc-last:\t%d\t\t", model->siesta.md_fc_first, model->siesta.md_fc_last));
        gui_text_show(OUTPUT, model->basename);
        gui_text_show(OUTPUT, "\n");
        g_free(filename);
    }

    model->basename = g_strdup(orig_basename);

}

void siesta_make_runscript(gchar * model_name, gchar * directory, struct grid_config_pak * grid)
{
    FILE *script;

    gchar *script_name;

    script_name = g_strdup_printf("run.script.%s", model_name);

    g_build_filename(directory, script_name, NULL);

    script = fopen (script_name, "w");
    if (!script_name)
        return;

    fprintf(script, "#!/bin/bash\n");
    fprintf(script, "#PBS -P d64\n");
    fprintf(script, "#PBS -l ncpus=8\n");
    fprintf(script, "#PBS -l vmem=2800mb\n");
    fprintf(script, "#PBS -l walltime=2:00:00\n");
    fprintf(script, "#PBS -l software=siesta\n");
    fprintf(script, "#PBS -l other=mpi\n");
    fprintf(script, "#PBS -wd\n");

    fprintf(script, "module load siesta\n");
    fprintf(script, "mpirun siesta < %s.fdf > %s.sot\n", model_name, model_name);

    fclose (script);
    g_free (script_name);

}

void siesta_load_animation(GtkWidget *w, struct model_pak *model)
{
    //change the current animation to the slider value
    //slider value auto updated due to shortcut.


if (!model->siesta.eigen_values)
  return;


    g_free(model->siesta.freq_disp_str);

/*
    model->siesta.current_frequency = make_eigvec_freq(model->siesta.eigen_values->ve[model->siesta.sorted_eig_values[model->siesta.current_animation]]);
*/

    model->siesta.current_frequency = make_eigvec_freq(mesch_ve_get(model->siesta.eigen_values, model->siesta.sorted_eig_values[model->siesta.current_animation]));

    model->siesta.freq_disp_str = g_strdup_printf("%.2f", model->siesta.current_frequency);

    gtk_entry_set_text((GtkEntry *) model->siesta.freq_text_box, model->siesta.freq_disp_str);
}

void siesta_animation_dialog_destroy(GtkWidget *w, gpointer *dialog)
{
    struct model_pak * model;

    model = sysenv.active_model;

    //stop the animation timer
/*
    siesta_destroy_timer(model);
*/

    //kill the dialog
    dialog_destroy(w, dialog);
}

void siesta_save_n_queue(GtkWidget *w, struct model_pak * model)
{
    int i, jobID=0;
    int total_files;
    gchar * filename, * orig_basename, *queuepath, *jobID_string;

    total_files = model->num_atoms / model->siesta.atoms_per_job +1;

    orig_basename = g_strdup(model->basename);

    struct grid_config_pak * grid_pak=NULL;

    //default jobstorage dir??
    //"$homedir/.gdis_jobs/  ?
    queuepath = g_build_filename(g_get_home_dir(), ".gdis_jobman", NULL);

    jobID_string = g_strdup_printf("%d", jobID);

    if (g_file_test (queuepath, G_FILE_TEST_IS_DIR))
    {
        // yay - it exists.
        // make new directory
//        grid_new_job();

        //change to dir.

        for (i=0; i<total_files; i++)
        {
            //set the range for md_fc_first, and md_fc_last
            model->siesta.md_fc_first = i * model->siesta.atoms_per_job + 1;
            model->siesta.md_fc_last = i * model->siesta.atoms_per_job + model->siesta.atoms_per_job;
            filename = g_strdup_printf("%s.%d.%d", orig_basename, i+1, total_files);
            if (model->siesta.md_fc_last >= model->num_atoms)
            {
                i++;
                model->siesta.md_fc_last = model->num_atoms;
            }
            model->basename = g_strdup(filename);

            filename = g_strdup_printf("%s.fdf", filename);

            // call the fdf save bit....
            file_save_as(filename, model);

            //grab current grid_config_pak

            //get job directory
            siesta_make_runscript(model->basename, sysenv.cwd, grid_pak);

            gui_text_show(OUTPUT, g_strdup_printf("fc-first:\t%d\tfc-last:\t%d\t\t", model->siesta.md_fc_first, model->siesta.md_fc_last));
            gui_text_show(OUTPUT, model->basename);
            gui_text_show(OUTPUT, "\n");
            g_free(filename);
        }

        model->basename = g_strdup(orig_basename);


    }
    else
    {
        //no queue directory

    }




}
