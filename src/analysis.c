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
#include <math.h>
#include <string.h>
#include <time.h>

#include "gdis.h"
#include "coords.h"
#include "model.h"
#include "edit.h"
#include "file.h"
#include "matrix.h"
#include "measure.h"
#include "parse.h"
#include "graph.h"
#include "gui_shorts.h"
#include "interface.h"
#include "analysis.h"
#include "numeric.h"
#include "count.h"

extern struct elem_pak elements[];
extern struct sysenv_pak sysenv;

#define MAX_PIE_SHELL 0


/* determine number of vector repeats to satisfy a given distance requirement in a given direction */
/* NB: direction should be given in normalize form */
gint calc_rep_3d(gdouble *vec, gdouble *normal, gdouble max)
{
gdouble len, x[3];

/* project vector and compute length */
ARR3SET(x, vec);
ARR3MUL(x, normal);
len = VEC3MAG(x);

/* orthogonal */
if (len <  G_MINDOUBLE)
  return(0);
/* one vector is all we need */
if (len > max)
  return(1);

/* return multiple required (NB: round up for safety) */
return(0.999 + max/len);
}

/* separation + lattice => compute repetitions */

void test(gint *images, gdouble *latmat, gdouble max)
{
gint i, j, m;
gdouble n[3], v[3];

/* compute minimum lattice images needed to ensure the 3 cartesian directions */
/* will be filled out to a specified minimum distance */
VEC3SET(images, G_MAXINT, G_MAXINT, G_MAXINT);
for (i=0 ; i<3 ; i++)
  {
/* normal */
  VEC3SET(n, 0.0, 0.0, 0.0);
  n[i] = 1.0;

  for (j=0 ; j<3 ; j++)
    {
/* iterate through lattice vectors */
    VEC3SET(v, latmat[j], latmat[j+3], latmat[j+6]);

/* record minimum non-zero value for each lattice vector */
    m = calc_rep_3d(v, n, max);
    if (m > 0) 
      {
      if (m < images[j])
        images[j] = m; 
      }
    }
  }

/*
for (i=0 ; i<3 ; i++)
  {
  printf("[%d] ", min[i]);
  }
printf("\n");
*/

}


/* fill out bins (NB: +) for 2 points, given periodic boundrary conditions */
void count_pie_3d(gdouble *xsep, gdouble *latmat, gpointer count)
{
gint i, j, k, images[3];
gdouble r, x[3], xfa[3], xfb[3], x1[3], x2[3], offset[3];

/* relocate so centroid is at the middle of the cell */
ARR3SET(x, xsep);
VEC3MUL(x, 0.5);
VEC3SET(xfa, 0.5, 0.5, 0.5);
ARR3SUB(xfa, x);
VEC3SET(xfb, 0.5, 0.5, 0.5);
ARR3ADD(xfb, x);

/* compute images needed to satsify the count max distance */
test(images, latmat, count_stop(count));

/* FIXME - need to correct the volume */

printf("using extents: [%d][%d][%d]\n", images[0], images[1], images[2]);

for (k=-(images[2]-1) ; k<images[2] ; k++)
  {
  for (j=-(images[1]-1) ; j<images[1] ; j++)
    {
    for (i=-(images[0]-1) ; i<images[0] ; i++)
      {

/*
printf("computing distance for: [%d][%d][%d]\n", i, j, k);
*/

/* cartesian offset */
      VEC3SET(offset, i, j, k);

      ARR3SET(x1, xfa);
      ARR3SET(x2, xfb);
      ARR3ADD(x2, offset);

      ARR3SUB(x2, x1);
      vecmat(latmat, x2);
      r = VEC3MAG(x2);

/* avoid adding 0.0 - ie distance between atom and itself */
if (r > G_MINDOUBLE)
  {
      count_insert(r, count);
  }

      }
    }
  }

}

/****************************/
/* display analysis results */
/****************************/
void analysis_show(struct model_pak *model)
{
g_assert(model != NULL);

/* clear any other special objects displayed */
model->picture_active = NULL;

/* display the new plot */
/*
tree_model_refresh(model);
*/

/* TODO - refresh_tree/gui/interface ? */
sysenv.refresh_dialog = TRUE;
redraw_canvas(SINGLE);
}

/**************************/
/* read in all frame data */
/**************************/
#define DEBUG_LOAD_ANALYSIS 1
gint analysis_load(struct analysis_pak *analysis, struct model_pak *model, struct task_pak *task)
{
gint i, j, k, m, n;
GSList *list;
struct core_pak *core;
struct model_pak temp;
FILE *fp;

/* checks */
g_assert(analysis != NULL);
g_assert(model != NULL);
g_assert(task != NULL);

m = analysis->num_frames;
n = analysis->num_atoms;

#if DEBUG_LOAD_ANALYSIS
printf("Allocating for %d frames and %d atoms.\n", m, n);
#endif

analysis->latmat = g_malloc((m+1)*sizeof(struct gd9_pak));

analysis->time = g_malloc((m+1)*sizeof(gdouble));
analysis->ke = g_malloc((m+1)*sizeof(gdouble));
analysis->pe = g_malloc((m+1)*sizeof(gdouble));
analysis->temp = g_malloc((m+1)*sizeof(gdouble));

analysis->position = g_malloc((m+1)*n*sizeof(struct gd3_pak));
analysis->velocity = g_malloc((m+1)*n*sizeof(struct gd3_pak));

/* init temp model */
model_init(&temp);
temp.frame_list = g_list_copy(model->frame_list);
temp.id = model->id;
temp.periodic = model->periodic;
temp.fractional = model->fractional;
/* init to use external trj frames */
if (model->id == GULP)
  {
  temp.gulp.trj_file = g_strdup(model->gulp.trj_file);
  temp.header_size = model->header_size;
  temp.frame_size = model->frame_size;
  temp.file_size = model->file_size;
  temp.expected_cores = model->expected_cores;
  temp.expected_shells = model->expected_shells;
  temp.trj_swap = model->trj_swap;
  memcpy(temp.latmat, model->latmat, 9*sizeof(gdouble));
  memcpy(temp.ilatmat, model->ilatmat, 9*sizeof(gdouble));
  temp.construct_pbc = TRUE;
  }

fp = fopen(model->filename, "r");
if (!fp)
  {
  printf("Could not open source file.");
  return(1);
  }
else
  {
#if DEBUG_LOAD_ANALYSIS
printf("Reading frames from: %s\n", model->filename);
#endif
  }

/* read frames */
k=0;
for (i=0 ; i<m ; i++)
  {
/* NB: frames start at 1 */
  read_raw_frame(fp, i+1, &temp);

/* NB: this causes problems - i guess model_prep() must have some gui update */
/* stuff in it that screws up as we're running this in a thread */
/*
  model_prep(&temp);
*/

  memcpy((analysis->latmat+i)->m, temp.latmat, 9*sizeof(gdouble));

/* fill out analysis frame with data from file load */

  if (!i)
    analysis->time_start = temp.frame_time;
  *(analysis->time+i) = temp.frame_time;
  *(analysis->ke+i) = temp.frame_ke;
  *(analysis->pe+i) = temp.frame_pe;
  *(analysis->temp+i) = temp.frame_temp;

  j = g_slist_length(temp.cores);

/* fill out coordinates */
  if (j == n)
    {
    j=0;
    for (list=temp.cores ; list ; list=g_slist_next(list))
      {
      core = list->data;
      ARR3SET((analysis->position+i*n+j)->x, core->x);
      ARR3SET((analysis->velocity+i*n+j)->x, core->v);
      j++;
      }
    k++;
    }
  else
    printf("WARNING: inconsistent frame %d has %d/%d cores.\n", i, j, n);

/* progress */
  task->progress += 50.0/analysis->num_frames;
  }

analysis->num_frames = k;
analysis->time_stop = temp.frame_time;

#if DEBUG_LOAD_ANALYSIS
printf("Read in %d complete frames.\n", k);
#endif

/* NB: we only copied the list from source model */
/* so make sure we don't free the elements */
g_list_free(temp.frame_list);
temp.frame_list = NULL;
model_free(&temp);
return(0);
}

/***********************/
/* free all frame data */
/***********************/
void analysis_free(struct analysis_pak *analysis)
{
if (!analysis)
  return;

if (analysis->latmat)
  g_free(analysis->latmat);
if (analysis->time)
  g_free(analysis->time);
if (analysis->ke)
  g_free(analysis->ke);
if (analysis->pe)
  g_free(analysis->pe);
if (analysis->temp)
  g_free(analysis->temp);
if (analysis->position)
  g_free(analysis->position);
if (analysis->velocity)
  g_free(analysis->velocity);

g_free(analysis);
}

/***********************************/
/* allocate the analysis structure */
/***********************************/
gpointer analysis_new(void)
{
struct analysis_pak *analysis;

analysis = g_malloc(sizeof(struct analysis_pak));

analysis->latmat = NULL;
analysis->time = NULL;
analysis->ke = NULL;
analysis->pe = NULL;
analysis->temp = NULL;
analysis->position = NULL;
analysis->velocity = NULL;

return(analysis);
}

/***********************************/
/* duplicate an analysis structure */
/***********************************/
/* intended to preserve the dialog values at the time a job is requested */
gpointer analysis_dup(gpointer data)
{
struct analysis_pak *src=data, *dest;

g_assert(src != NULL);

dest = analysis_new();

dest->start = src->start;
dest->stop = src->stop;
dest->step = src->step;

dest->num_atoms = src->num_atoms;
dest->num_frames = src->num_frames;

dest->atom1 = g_strdup(src->atom1);
dest->atom2 = g_strdup(src->atom2);

dest->rdf_normalize = src->rdf_normalize;

return(dest);
}

/**********************/
/* setup the analysis */
/**********************/
#define DEBUG_INIT_ANALYSIS 0
gint analysis_init(struct model_pak *model)
{
struct analysis_pak *analysis;

/* checks */
g_assert(model != NULL);
g_assert(model->analysis != NULL);

analysis = model->analysis;

analysis->num_atoms = g_slist_length(model->cores);
analysis->atom1 = NULL;
analysis->atom2 = NULL;
analysis->num_frames = model->num_frames;

/* TODO - should really be renamed eg rdf_start/stop/step */
analysis->start = 0.0;
/* make integer - looks nicer */
analysis->stop = (gint) model->rmax;
analysis->step = 0.1;

analysis->rdf_normalize = TRUE;

return(0);
}

/*******************/
/* RDF calculation */
/*******************/
#define DEBUG_CALC_RDF 0
void analysis_plot_rdf(struct analysis_pak *analysis, struct task_pak *task)
{
gint i, j, n;
gint size;
gint *bins;
gdouble *cz, slice;
gdouble volume, c, r1, r2, xi[3], xj[3];
gdouble *latmat;
/*
gdouble latmat[9];
*/
gpointer graph, count1;
gchar *label;
struct core_pak *corei, *corej;
struct model_pak *model;

/* checks */
g_assert(analysis != NULL);
g_assert(task != NULL);
model = task->locked_model;
g_assert(model != NULL);

/* load frames? */
if (!analysis->position)
  {
/* time slice - if load_analysis() needed assume it consumes half ie 50.0 */
  slice = 50.0;
  if (analysis_load(analysis, model, task))
    return;
  }
else
  slice = 100.0;


count1 = count_new(analysis->start, analysis->stop, analysis->step);


#if DEBUG_CALC_RDF
printf("%f - %f : %f\n", analysis->start, analysis->stop, analysis->step);
printf("%s - %s\n", analysis->atom1, analysis->atom2);
#endif

/* loop over frames */
volume = 0.0;
for (n=0 ; n<analysis->num_frames ; n++)
  {
/*
  memcpy(latmat, (analysis->latmat+n)->m, 9*sizeof(gdouble));
*/

  latmat = (analysis->latmat+n)->m;

  volume += calc_volume(latmat);

/* loop over unique atom pairs */
/*
  for (i=0 ; i<analysis->num_atoms-1 ; i++)
*/

/* CURRENT - loop over all atoms - even duplicates so that we can record */
/* distance between an atom and it's periodic images (NB: correct for double counting) */
  for (i=0 ; i<analysis->num_atoms ; i++)
    {
    corei = g_slist_nth_data(model->cores, i);

    ARR3SET(xi, (analysis->position+n*analysis->num_atoms+i)->x);

    for (j=0 ; j<analysis->num_atoms ; j++)
      {
      corej = g_slist_nth_data(model->cores, j);

      if (!pair_match(analysis->atom1, analysis->atom2, corei, corej))
        continue;

/* calculate closest distance */
      ARR3SET(xj, (analysis->position+n*analysis->num_atoms+j)->x);
      ARR3SUB(xj, xi);
      fractional_min(xj, model->periodic);

#define OLD 0
#if OLD
      vecmat(latmat, xj);
      r1 = VEC3MAG(xj); 
      count_insert(r1, count1);
#else
      count_pie_3d(xj, latmat, count1);
#endif

      }
    }

/* update progess */
  task->progress += slice/analysis->num_frames;
  }

/* convert count integer array into gdouble array */
size = count_size(count1);
bins = count_bins(count1);
cz = g_malloc(size * sizeof(gdouble));
for (i=size ; i-- ; )
  {
  *(cz+i) = (gdouble) *(bins+i);
  
/* correct for double counting */
  *(cz+i) *= 0.5;
  }

count_free(count1);

/* normalize against ideal gas with same AVERAGE volume */
if (analysis->rdf_normalize)
  {
  volume /= analysis->num_frames;
  c = 1.333333333*PI*analysis->num_atoms/volume;
  c *= analysis->num_atoms*analysis->num_frames;
  r1 = r2 = analysis->start;
  r2 += analysis->step;
  for (i=0 ; i<size ; i++)
    {
    *(cz+i) /= (c * (r2*r2*r2 - r1*r1*r1));
    r1 += analysis->step;
    r2 += analysis->step;
    }

  label = g_strdup_printf("RDF_%s_%s", analysis->atom1, analysis->atom2);
  }
else
  label = g_strdup_printf("Pair_%s_%s", analysis->atom1, analysis->atom2);


/* FIXME - alt method? (A&T p54) */
/*
ntot = analysis->num_frames * analysis->num_atoms;
c = volume / (ntot*ntot);
for (i=0 ; i<nz ; i++)
  *(cz+i) *= c;
*/

/* NB: add graph to the model's data structure, but don't update */
/* the model tree as we're still in a thread */
graph = graph_new(label, model);
g_free(label);

graph_add_data(size, cz, analysis->start, analysis->stop, graph);
graph_set_yticks(TRUE, 2, graph);
g_free(cz);

analysis_free(analysis);
}

/**************************************/
/* plot the evolution of measurements */
/**************************************/
/* CURRENT - only do the 1st, until we figure out if we should */
/* distinguish bonds angles etc ans selectable issues etc */
void analysis_plot_meas(struct analysis_pak *analysis, struct task_pak *task)
{
gint i, j, n, type, index[4];
gdouble value, *data, x[4][3];
gpointer m, graph;
GSList *list, *clist;
struct model_pak *model;

/* checks */
g_assert(analysis != NULL);
g_assert(task != NULL);
model = task->locked_model;
g_assert(model != NULL);

/* load frames? */
if (!analysis->position)
  if (analysis_load(analysis, model, task))
    return;

/* CURRENT - only dealing with 1st measurement */
m = g_slist_nth_data(model->measure_list, 0);
if (!m)
  {
  printf("No measurements to plot.\n");
  return;
  }

type = measure_type_get(m);
clist = measure_cores_get(m);

/*
printf("measurement cores: %d\n", g_slist_length(clist));
*/

/* record the position of each core */
n=0;
for (list=clist ; list ; list=g_slist_next(list))
  {
  if (n<4)
    {
    index[n] = g_slist_index(model->cores, list->data);
/*
printf(" - %p (%d)\n", list->data, index[n]);
*/
    n++;
    }
  else
    break;
  }

/* allocation for graph plot */
data = g_malloc(model->num_frames*sizeof(gdouble));

/* loop over each frame */
for (i=0 ; i<analysis->num_frames ; i++)
  {

/*
printf("frame: %d ", i);
*/

/* get cartesian coords for constituent atoms */
  for (j=0 ; j<n ; j++)
    {
    ARR3SET(&x[j][0], (analysis->position+i*analysis->num_atoms+index[j])->x);
    vecmat(model->latmat, x[j]);
    }

/* CURRENT - make measurement -> data */
  switch (type)
    {
    case MEASURE_TORSION:
/* FIXME - check this */
      value = measure_dihedral(x[0], x[1], x[2], x[3]);
      break;

    case MEASURE_ANGLE:
      value = measure_angle(x[0], x[1], x[2]);
      break;

    default:
      value = measure_distance(x[0], x[1]);
    }

/*
printf(" (value = %f)\n", value);
*/

  *(data+i) = value;
  }

/* graph */
graph = graph_new("Measurement", model);
graph_add_data(model->num_frames, data, 1, model->num_frames, graph);

/* done */
g_free(data);
g_slist_free(clist);

analysis_free(analysis);
}

/*********************************************/
/* plot the FFT of the VACF (power spectrum) */
/*********************************************/
void analysis_plot_power(gint m, gdouble *vacf, gdouble step, struct model_pak *model)
{
gint i, e, n;
gdouble r_power, s_power, *x, *y;
struct graph_pak *graph;

/* checks */
g_assert(vacf != NULL);
g_assert(model != NULL);

/* get the best n = 2^e >= m */
e = log((gdouble) m) / log(2.0);
n = pow(2.0, ++e);

/*
printf("input size: %d, fft size: %d\n", m, n);
*/

g_assert(n >= m);
x = g_malloc(2*n*sizeof(gdouble));

/* assign data */
r_power = 0.0;
for (i=0 ; i<m ; i++)
  {
  x[2*i] = vacf[i];
  x[2*i+1] = 0.0; 
  r_power += vacf[i] * vacf[i];
  }

/* pad the rest */
for (i=2*m ; i<2*n ; i++)
  x[i] = 0.0; 

/* compute fourier transform */
fft(x, n, 1);

/* compute power spectrum */
y = g_malloc(n*sizeof(gdouble));
s_power = 0.0;
for (i=0 ; i<n ; i++)
  {
/* see numerical recipies in C p401 */
  y[i] = x[2*i]*x[2*i] + x[2*i+1]*x[2*i+1];
  s_power += y[i];
  }

/* see num recipics - p407 (12.1.10) */
printf("total power (real)    : %f\n", r_power);
printf("total power (spectral): %f\n", s_power / (gdouble) m);

/* plot */
graph = graph_new("POWER", model);
/* only plot +ve frequencies (symmetrical anyway, since input is real) */
graph_add_data(n/2, y, 0.0, 0.5/step, graph);
graph_set_yticks(FALSE, 2, graph);

g_free(x);
g_free(y);
}

/********************/
/* VACF calculation */
/********************/
#define DEBUG_CALC_VACF 0
void analysis_plot_vacf(struct analysis_pak *analysis, struct task_pak *task)
{
gint i, n;
gdouble slice, step=1.0;
gdouble *vacf, v0[3], vi[3];
gpointer graph;
struct model_pak *model;

/* checks */
g_assert(analysis != NULL);
g_assert(task != NULL);
model = task->locked_model;
g_assert(model != NULL);

/* vacf setup */
vacf = g_malloc(analysis->num_frames * sizeof(gdouble));

/* load frames? */
if (!analysis->position)
  {
/* time slice - if load_analysis() needed assume it consumes half ie 50.0 */
  slice = 50.0;
  if (analysis_load(analysis, model, task))
    return;
  }
else
  slice = 100.0;

for (n=0 ; n<analysis->num_frames ; n++)
  {
#if DEBUG_CALC_VACF
printf(" --- frame %d\n", n);
#endif
  vacf[n] = 0.0;
  for (i=0 ; i<analysis->num_atoms ; i++)
    {
    ARR3SET(v0, (analysis->velocity+i)->x);
    ARR3SET(vi, (analysis->velocity+n*analysis->num_atoms+i)->x);

    ARR3MUL(vi, v0);
    vacf[n] += vi[0]+vi[1]+vi[2];

#if DEBUG_CALC_VACF
printf("atom %d, ", i);
P3VEC("vi : ", vi);
#endif
    }
  vacf[n] /= (gdouble) analysis->num_atoms;
/* update progess */
  task->progress += slice/analysis->num_frames;
  }

step = analysis->time_stop - analysis->time_start;
step /= (gdouble) (analysis->num_frames-1);

#if DEBUG_CALC_VACF
printf("t : [%f - %f](%d x %f)\n", analysis->time_start, analysis->time_stop,
                                   analysis->num_frames, step);
#endif

/* NB: add graph to the model's data structure, but don't update */
/* the model tree as we're still in a thread */
graph = graph_new("VACF", model);
graph_add_data(analysis->num_frames, vacf,
                  analysis->time_start, analysis->time_stop, graph);
graph_set_yticks(FALSE, 2, graph);

/* NEW */
/* TODO - make optional */
analysis_plot_power(analysis->num_frames, vacf, step, model);

g_free(vacf);

analysis_free(analysis);
}

/********************/
/* temperature plot */
/********************/
void analysis_plot_temp(struct analysis_pak *analysis, struct task_pak *task)
{
gpointer graph;
struct model_pak *model;

/* checks */
g_assert(analysis != NULL);
g_assert(task != NULL);
model = task->locked_model;
g_assert(model != NULL);

/* load frames? */
if (!analysis->position)
  if (analysis_load(analysis, model, task))
    return;

graph = graph_new("Temp", model);
graph_add_data(analysis->num_frames, analysis->temp,
                  analysis->time_start, analysis->time_stop, graph);

analysis_free(analysis);
}

/***********************/
/* kinetic energy plot */
/***********************/
void analysis_plot_ke(struct analysis_pak *analysis, struct task_pak *task)
{
gpointer graph;
struct model_pak *model;

/* checks */
g_assert(analysis != NULL);
g_assert(task != NULL);
model = task->locked_model;
g_assert(model != NULL);

/* load frames? */
if (!analysis->position)
  if (analysis_load(analysis, model, task))
    return;

graph = graph_new("KE", model);
graph_add_data(analysis->num_frames, analysis->ke,
                  analysis->time_start, analysis->time_stop, graph);

analysis_free(analysis);
}

/*************************/
/* potential energy plot */
/*************************/
void analysis_plot_pe(struct analysis_pak *analysis, struct task_pak *task)
{
gpointer graph;
struct model_pak *model;

/* checks */
g_assert(analysis != NULL);
g_assert(task != NULL);
model = task->locked_model;
g_assert(model != NULL);

/* load frames? */
if (!analysis->position)
  if (analysis_load(analysis, model, task))
    return;

graph = graph_new("PE", model);
graph_add_data(analysis->num_frames, analysis->pe,
                  analysis->time_start, analysis->time_stop, graph);

analysis_free(analysis);
}
