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
#include <unistd.h>

#include "gdis.h"
#include "coords.h"
#include "edit.h"
#include "error.h"
#include "file.h"
#include "gulp_keyword.h"
#include "parse.h"
#include "scan.h"
#include "space.h"
#include "matrix.h"
#include "model.h"
#include "numeric.h"
#include "interface.h"
#include "ff.h"
#include "task.h"

#define DEBUG_MORE 0
#define MAX_KEYS 15

/* main structures */
extern struct sysenv_pak sysenv;
extern struct elem_pak elements[];

/* external function from ff-gulp.c */
gint ff_gulp_write(FILE *, struct model_pak *);

/*********************************/
/* GULP structure initialization */
/*********************************/
void gulp_init(gpointer data)
{
struct gulp_pak *gulp = data;

g_assert(gulp != NULL);

/* gulp template init */
gulp->maxcyc=0;
gulp->energy = 0.0;
gulp->no_esurf=FALSE;
gulp->esurf[0]=0.0;
gulp->esurf[1]=0.0;
gulp->esurf_units=g_strdup(" ");
gulp->no_eatt=FALSE;
gulp->eatt[0]=0.0;
gulp->eatt[1]=0.0;
gulp->eatt_units=g_strdup(" ");
gulp->sbulkenergy=0.0;
gulp->sdipole=999999999.0;
gulp->gnorm = -1.0;
gulp->qsum = 0.0;
gulp->no_exec = FALSE;
gulp->run = E_SINGLE;
gulp->rigid = FALSE;
gulp->rigid_move = g_strdup("z");
gulp->qeq = FALSE;
gulp->free = FALSE;
gulp->zsisa = FALSE;
gulp->compare = FALSE;
gulp->prop = FALSE;
gulp->nosym = FALSE;
gulp->fix = FALSE;
gulp->noautobond = FALSE;
gulp->method = CONP;
/* dynamics info */
/*
gulp->frame_time = 0.0;
gulp->frame_ke = 0.0;
gulp->frame_pe = 0.0;
gulp->frame_temp = 0.0;
*/
/* NEW - original frame data (stored since reading in trajectory overwrites) */
gulp->orig_fractional = FALSE;
gulp->orig_construct_pbc = FALSE;
gulp->orig_cores = NULL;
gulp->orig_shells = NULL;
gulp->potentials = NULL;
gulp->elements = NULL;
gulp->species = NULL;
/* potgrid replacement */
gulp->epot_vecs = NULL;
gulp->epot_vals = NULL;
gulp->epot_min = 999999999.9;
gulp->epot_max = -999999999.9;

gulp->optimiser = -1;
gulp->optimiser2 = SWITCH_OFF;
gulp->switch_type = -1;
gulp->switch_value = 0.1;
gulp->unit_hessian = FALSE;
gulp->ensemble = -1;
gulp->pressure = NULL;
gulp->temperature = NULL;
VEC3SET(gulp->super, 1, 1, 1);
gulp->coulomb = NOBUILD;
gulp->libfile = NULL;
gulp->output_extra_keywords = TRUE;
gulp->output_extra = FALSE;
gulp->extra_keywords = NULL;
gulp->extra = NULL;

gulp->timestep = NULL;
gulp->equilibration = NULL;
gulp->production = NULL;
gulp->sample = NULL;
gulp->write = NULL;

/* io filenames */
gulp->temp_file = NULL;
gulp->dump_file = NULL;
gulp->out_file = NULL;
gulp->trj_file = NULL;
gulp->mov_file = NULL;

/* io options */
gulp->print_charge = TRUE;

gulp->solvation_model = GULP_SOLVATION_NONE;
gulp->cosmo_sas = TRUE;
gulp->cosmo_shape = 0;
gulp->cosmo_shape_index[0] = 5;
gulp->cosmo_shape_index[1] = 0;
gulp->cosmo_points = 974;
gulp->cosmo_segments = 110;
gulp->cosmo_smoothing = 0.0;
gulp->cosmo_solvent_epsilon = 1.0;
gulp->cosmo_solvent_radius = 1.0;
gulp->cosmo_solvent_delta = 0.1;
gulp->cosmo_solvent_rmax = 10.0;
}

/***********************************************************/
/* use number of points and shape to compute shape indices */
/***********************************************************/
#define COSMO_SHAPE_INDEX_MAX 10
gint gulp_shape_indices_compute(struct gulp_pak *gulp)
{
gint i, j, n;
gdouble a, b, c;

g_assert(gulp != NULL);

switch (gulp->cosmo_shape)
  {
  case 0:
    c = 4.0;
    break;
  default:
    c = 10.0;
    break;
  }

a = gulp->cosmo_points;
a = log((a - 2.0) / c) / log(3.0);
b = log(4.0) / log(3.0);

/* simple for loop search for suitable index values */
for (j=0 ; j<COSMO_SHAPE_INDEX_MAX ; j++)
  {
  i = nearest_int(a - j*b);

  n = nearest_int (c * pow(3.0, (gdouble) i) * pow(4.0, (gdouble) j) + 2.0);

/*
printf("[%d,%d] : %d / %d\n", i, j, n, gulp->cosmo_points);
*/

  if (!(n - gulp->cosmo_points))
    {
    gulp->cosmo_shape_index[0] = i;
    gulp->cosmo_shape_index[1] = j;

    return(TRUE);
    }
  }
return(FALSE);
}

/*****************************/
/* core file print primitive */
/*****************************/
void fprint_core(FILE *fp, struct core_pak *core, struct model_pak *model)
{
gint i;
gdouble vec[3];

/* get cartesian coords */
ARR3SET(vec, core->x);
vecmat(model->latmat, vec); 

/* set fractional part */
for (i=0 ; i<model->periodic ; i++)
  if (model->fractional)
    vec[i] = core->x[i];

/* only print specific atom charges if we read them in */
/* in that fashion (ie lookup_charge is false) */
if (core->breathe)
  fprintf(fp,"%-4s bcore %11.6f %11.6f %11.6f", core->atom_label,
                                                vec[0], vec[1], vec[2]);
else
  fprintf(fp,"%-4s core  %11.6f %11.6f %11.6f", core->atom_label,
                                                vec[0], vec[1], vec[2]);

/* NB: 7 dp for charges - 8 occasionally introduced small errors (cf GULP) */
if (model->gulp.print_charge)
  if (!core->lookup_charge)
    fprintf(fp,"  %11.7lf", core->charge);

/* 7dp for occupancy - to cover cases like 1/3 */
if (core->has_sof)
  fprintf(fp,"  %9.7f", core->sof);

if (core->radius != 0.0)
  fprintf(fp,"  %7.5f", core->radius);

if (core->flags)
  fprintf(fp," %s", core->flags);

if (model->periodic == 2)
  if (core->growth)
    fprintf(fp, " %%");

if (core->translate)
  fprintf(fp, " T");

fprintf(fp,"\n");
}

/******************************/
/* shell file print primitive */
/******************************/
void fprint_shell(FILE *fp, struct shel_pak *shell, struct model_pak *model)
{
gint i;
gdouble vec[3];

/* get cartesian coordinates */
ARR3SET(vec, shell->x);
vecmat(model->latmat, vec); 

/* set fractional part */
for (i=0 ; i<model->periodic ; i++)
  if (model->fractional)
    vec[i] = shell->x[i];

if (shell->breathe)
  fprintf(fp,"%-4s bshel %11.6f %11.6f %11.6f", shell->shell_label,
                                                vec[0], vec[1], vec[2]);
else
  fprintf(fp,"%-4s shel  %11.6f %11.6f %11.6f", shell->shell_label,
                                                vec[0], vec[1], vec[2]);

/* only print specific atom charges if we read them in */
if (model->gulp.print_charge)
  if (!shell->lookup_charge)
    fprintf(fp,"  %11.7lf", shell->charge);

/* 7dp for occupancy - to cover cases like 1/3 */
if (shell->has_sof)
  fprintf(fp,"  %9.7f", shell->sof);

if (shell->radius != 0.0)
  fprintf(fp,"  %7.5f", shell->radius);

if (shell->flags)
  fprintf(fp," %s", shell->flags);

if (model->periodic == 2)
  if ((shell->core)->growth)
    fprintf(fp, " %%");

if (shell->translate)
  fprintf(fp, " T");

fprintf(fp,"\n");
}

/*********************************************/
/* enforce valid input/output GULP filenames */
/*********************************************/
void gulp_files_init(struct model_pak *model)
{
gchar *basename;

g_assert(model != NULL);
g_assert(model->basename != NULL);

/* remove spaces from basename */
basename = g_strdup(model->basename);
parse_space_replace(basename, '_');

/* must have input and ouptut filenames */
if (!model->gulp.temp_file)
  {
  model->gulp.temp_file = g_strdup_printf("%s.gin", model->basename);
  parse_space_replace(model->gulp.temp_file, '_');
  }
if (!model->gulp.out_file)
  {
  model->gulp.out_file = g_strdup_printf("%s.got", model->basename);
  parse_space_replace(model->gulp.out_file, '_');
  }
/* default dump file - gulp opti ASSUMES a dump file - see proc_gulp_task() */
if (!model->gulp.dump_file)
  {
  model->gulp.dump_file = g_strdup_printf("%s.res", model->basename);
  parse_space_replace(model->gulp.dump_file, '_');
  }
}

/*********************************/
/* cleanup for GULP file control */
/*********************************/
void gulp_files_free(struct model_pak *model)
{
g_assert(model != NULL);

g_free(model->gulp.temp_file);
g_free(model->gulp.dump_file);
g_free(model->gulp.out_file);
}

/**********************************/
/* free memory used for gulp data */
/**********************************/
#define DEBUG_FREE_GULP_DATA 0
void gulp_data_free(struct model_pak *model)
{
g_assert(model != NULL);

/* NEW - GULP files with associated trajectory only */
if (model->gulp.orig_cores)
  free_slist(model->gulp.orig_cores);
if (model->gulp.orig_shells)
  free_slist(model->gulp.orig_shells);

g_free(model->gulp.potentials);
g_free(model->gulp.elements);
g_free(model->gulp.species);
g_free(model->gulp.kpoints);
model->gulp.potentials=NULL;
model->gulp.elements=NULL;
model->gulp.species=NULL;
model->gulp.kpoints=NULL;

if (model->gulp.epot_vecs)
  free_slist(model->gulp.epot_vecs);
if (model->gulp.epot_vals)
  free_slist(model->gulp.epot_vals);
model->gulp.epot_vecs=NULL;
model->gulp.epot_vals=NULL;

g_free(model->gulp.libfile);
g_free(model->gulp.extra);
g_free(model->gulp.eatt_units);
g_free(model->gulp.esurf_units);
g_free(model->gulp.pressure);
g_free(model->gulp.temperature);
model->gulp.libfile=NULL;
model->gulp.extra=NULL;
model->gulp.eatt_units=NULL;
model->gulp.esurf_units=NULL;
model->gulp.pressure=NULL;
model->gulp.temperature=NULL;
}

/***********************************************/
/* copy unrecognized gulp keywords and options */
/***********************************************/
void gulp_extra_copy(struct model_pak *src, struct model_pak *dest)
{
if (src->gulp.extra)
  {
  g_free(dest->gulp.extra);
  dest->gulp.extra = g_strdup(src->gulp.extra);
  }
if (src->gulp.extra_keywords)
  {
  g_free(dest->gulp.extra_keywords);
  dest->gulp.extra_keywords = g_strdup(src->gulp.extra_keywords);
  }
}

/**********************/
/* transfer gulp data */
/**********************/
void gulp_data_copy(struct model_pak *src, struct model_pak *dest)
{
/* free any old data in the destination model */
gulp_data_free(dest);

/* copy keywords */
dest->gulp.no_exec = src->gulp.no_exec;
dest->gulp.run = src->gulp.run;
dest->gulp.method = src->gulp.method;
dest->gulp.optimiser = src->gulp.optimiser;
dest->gulp.optimiser2 = src->gulp.optimiser2;
dest->gulp.switch_type = src->gulp.switch_type;
dest->gulp.switch_value = src->gulp.switch_value;
dest->gulp.free = src->gulp.free;
dest->gulp.qeq = src->gulp.qeq;
dest->gulp.zsisa = src->gulp.zsisa;
dest->gulp.compare = src->gulp.compare;
dest->gulp.nosym = src->gulp.nosym;
dest->gulp.fix = src->gulp.fix;
dest->gulp.cycles = src->gulp.cycles;
dest->gulp.ensemble = src->gulp.ensemble;
dest->gulp.coulomb = src->gulp.coulomb;
dest->gulp.maxcyc = src->gulp.maxcyc;

/* FIXME AR change - copy over cosmo stuff - may cause problems with surfaces (but thats what I want) */
/* Note to SEAN would be a lot easier if were a separate struct as I could then replace with a single memcopy */
dest->gulp.solvation_model = src->gulp.solvation_model;
dest->gulp.cosmo_sas = src->gulp.cosmo_sas; 
dest->gulp.cosmo_shape = src->gulp.cosmo_shape; 
dest->gulp.cosmo_shape_index[0] = src->gulp.cosmo_shape_index[0]; 
dest->gulp.cosmo_shape_index[1] = src->gulp.cosmo_shape_index[1];
dest->gulp.cosmo_points = src->gulp.cosmo_points;
dest->gulp.cosmo_segments = src->gulp.cosmo_segments;
dest->gulp.cosmo_smoothing = src->gulp.cosmo_smoothing; 
dest->gulp.cosmo_solvent_epsilon = src->gulp.cosmo_solvent_epsilon;
dest->gulp.cosmo_solvent_radius = src->gulp.cosmo_solvent_radius;
dest->gulp.cosmo_solvent_delta = src->gulp.cosmo_solvent_delta;
dest->gulp.cosmo_solvent_rmax = src->gulp.cosmo_solvent_rmax;

/* FIXME - may cause problems when creating surfaces from unit cells */
dest->gulp.output_extra_keywords = src->gulp.output_extra_keywords;
dest->gulp.output_extra = src->gulp.output_extra;
gulp_extra_copy(src, dest);

/* copy surface related data */
if (dest->periodic == 2)
  dest->gulp.sbulkenergy = src->gulp.sbulkenergy;
if (dest->periodic == 2 && src->periodic == 2)
  dest->surface.dspacing = src->surface.dspacing;

/* duplicate gulp character data */
dest->gulp.libfile = g_strdup(src->gulp.libfile);
dest->gulp.potentials = g_strdup(src->gulp.potentials);
dest->gulp.elements = g_strdup(src->gulp.elements);
dest->gulp.species = g_strdup(src->gulp.species);
dest->gulp.pressure = g_strdup(src->gulp.pressure);
dest->gulp.temperature = g_strdup(src->gulp.temperature);
}

/**********************************/
/* debugging - dump a phonon mode */
/**********************************/
/* NB: modes are from 1 - num_phonons */
void dump_phonon(gint atom, gint mode, struct model_pak *model)
{
gchar *freq;
gdouble v[3];
GSList *xlist, *ylist, *zlist;
struct core_pak *core;

freq = g_slist_nth_data(model->phonons, mode-1);
core = g_slist_nth_data(model->cores, atom-1);

xlist = core->vibx_list; 
ylist = core->viby_list; 
zlist = core->vibz_list; 

v[0] = *((gdouble *) g_slist_nth_data(xlist, mode-1));
v[1] = *((gdouble *) g_slist_nth_data(ylist, mode-1));
v[2] = *((gdouble *) g_slist_nth_data(zlist, mode-1));

printf("[%s] : ", freq);
printf("[%s]", core->atom_label);
P3VEC(" : ", v);
}

/************************************************/
/* debugging - dump the vibrational frequencies */
/************************************************/
void dump_phonons(struct model_pak *model)
{
GSList *list;

for (list=model->phonons ; list ; list=g_slist_next(list))
  {
  printf("[%s]\n", (gchar *) list->data);
  }
}

/**********************************/
/* print appropriate region label */
/**********************************/
void write_region_header(FILE *fp, gint region, struct model_pak *model)
{
/* check periodicity */
switch (model->periodic)
  {
  case 1:
    if (model->fractional)
      fprintf(fp,"pfractional region %d", region+1);
    else
      fprintf(fp,"cartesian region %d", region+1);
    break;

  case 2:
    if (model->fractional)
      fprintf(fp,"sfractional region %d", region+1);
    else
      fprintf(fp,"cartesian region %d", region+1);
    break;

  case 3:
    if (model->fractional)
      fprintf(fp,"fractional");
    else
      fprintf(fp,"cartesian");
    break;

  default:
    fprintf(fp,"cartesian");
  }

/* NB: rigid xyz allows all translation as a rigid body */
/* rigid z allows vertical movement only -> mandatory for grid based dock */
if (model->gulp.rigid && region == 2)
  {
  fprintf(fp," rigid ");
  if (model->gulp.rigid_move)
    fprintf(fp,"%s ", model->gulp.rigid_move);
  }

fprintf(fp,"\n");
}

/*********************/
/* write a GULP file */
/*********************/
#define DEBUG_WRITE_GULP 0
gint write_gulp(gchar *filename, struct model_pak *data)
{
gint i, j, libflag, flag, region;
gdouble x[3];
GString *line;
GSList *list;
struct core_pak *core, *core1, *core2;
struct shel_pak *shell;
struct bond_pak *bond;
struct vec_pak *vec;
FILE *fp;

g_return_val_if_fail(data != NULL, 1);

/* is this really necessary??? */
unlink(filename);

/* is a valid potential library specified? */
libflag = 0;
if (data->gulp.libfile)
  {
  if (strlen(data->gulp.libfile))
    {
/* test for library existance */
    if (!access(data->gulp.libfile, R_OK))
      libflag++;
    else
      printf("Warning: library %s not found.\n", data->gulp.libfile);
    }
  }

/* off we go */
fp = fopen(filename, "wt");
if (!fp)
  {
  printf("Failed to open %s.\n", filename);
  return(2);
  }

/* run type */
line = g_string_new(NULL);
switch(data->gulp.run)
  {
  case E_SINGLE:
    g_string_assign(line,"single");
    break;
  case E_OPTIMIZE:
    g_string_assign(line,"opti");
    break;
  case MD:
    g_string_assign(line,"md");
    break;
/* no run -> single pt (valid default?) */
  default:
    g_string_assign(line,"single");
    break;
  }

switch(data->gulp.method)
  {
  case CONP:
    g_string_append(line," conp");
    break;
  case CONV:
    g_string_append(line," conv");
    break;
  default:
    break;
  }

switch(data->gulp.optimiser)
  {
  case BFGS_OPT:
    g_string_append(line," bfgs");
    break;
  case CONJ_OPT:
    g_string_append(line," conjugate");
    break;
  case RFO_OPT:
    g_string_append(line," rfo");
    break;
  }

if (data->gulp.unit_hessian)
    g_string_append(line," unit");

switch(data->gulp.coulomb)
  {
  case MOLE:
    g_string_append(line," mole");
    break;
  case MOLMEC:
    g_string_append(line," molmec");
    break;
  case MOLQ:
    g_string_append(line," molq");
    break;
  default:
    break;
  }

/* cosmo calc */
/* TODO - change variable name to solvation model */
switch (data->gulp.solvation_model)
  {
  case GULP_SOLVATION_COSMIC:
    g_string_append(line," cosmic");
    break;
  case GULP_SOLVATION_COSMO:
    g_string_append(line," cosmo");
    break;
  }

/* phonon calc */
if (data->gulp.phonon)
  {
  g_string_append(line," phonon");
  if (data->gulp.eigen)
    g_string_append(line," eigen");
  }

/* NEW - this makes it easy to use GULP to compute properties */
if (data->gulp.prop)
  g_string_append(line," eigen");

/* additional keywords */
if (data->gulp.qeq)
  g_string_append(line," qeq");
if (data->gulp.free)
  g_string_append(line," free");
if (data->gulp.zsisa)
  g_string_append(line," zsisa");
if (data->gulp.compare)
  g_string_append(line," comp");
if (data->gulp.fix)
  g_string_append(line," fix");
if (data->gulp.noautobond)
  g_string_append(line," noautobond");

/* if nosym specified, force output of non-prim cell for consistency */
/* FIXME - this will generate 2 x full if full is already specified, */
/* not a problem currently since gulp doesnt complain about duplication */
if (data->gulp.nosym)
  g_string_append(line," nosym full");
/* electrostatic potential calcs */
if (g_slist_length(data->gulp.epot_vecs))
  g_string_append(line," pot");

/* print out the keywords line */
if (data->gulp.output_extra_keywords && data->gulp.extra_keywords)
  fprintf(fp,"%s %s\n\n", line->str, data->gulp.extra_keywords);
else
  fprintf(fp,"%s\n\n",line->str);

fprintf(fp, "# ");
gdis_blurb(fp);
fprintf(fp, "#\n\n");

if (data->periodic == 2)
  {
/* NB: this surface data is only meaningful if generate from within GDIS */
  if (VEC3MAGSQ(data->surface.depth_vec) > FRACTION_TOLERANCE)
    {
    fprintf(fp, "#       miller: %d %d %d\n",
                     (gint) data->surface.miller[0],
                     (gint) data->surface.miller[1],
                     (gint) data->surface.miller[2]);
    fprintf(fp, "#        shift: %lf\n", data->surface.shift);
    fprintf(fp, "# region sizes: %d %d\n\n",
                    (gint) data->surface.region[0],
                    (gint) data->surface.region[1]);
    }
  }

/* optional data */
if (data->gulp.maxcyc > 0)
  fprintf(fp,"maxcyc %d\n\n", (gint) data->gulp.maxcyc);

if (data->gulp.optimiser2 != SWITCH_OFF)
  {
  fprintf(fp,"switch ");
  switch(data->gulp.optimiser2)
    {
    case CONJ_OPT:
      fprintf(fp,"conj ");
      break;
    case RFO_OPT:
      fprintf(fp,"rfo ");
      break;
    case BFGS_OPT:
    default:
      fprintf(fp,"bfgs ");
      break;
    }

  switch(data->gulp.switch_type)
    {
    case CYCLE:
      fprintf(fp,"cycle %d\n\n", (gint) data->gulp.switch_value);
      break;
    case GNORM:
    default:
      fprintf(fp,"gnorm %lf\n\n", data->gulp.switch_value);
      break;
    }
  }

/* only write if we have a valid supplied library */
if (libflag)
  fprintf(fp,"library %s\n\n", data->gulp.libfile);

fprintf(fp,"name %s\n\n", data->basename);

if (data->gulp.dump_file)
  if (strlen(data->gulp.dump_file))
    fprintf(fp,"\ndump %s\n\n", data->gulp.dump_file);

if (data->gulp.temperature)
  if (strlen(data->gulp.temperature))
    fprintf(fp,"temperature %s\n\n", data->gulp.temperature);

if (data->gulp.pressure)
  if (strlen(data->gulp.pressure))
    fprintf(fp,"pressure %s\n\n", data->gulp.pressure);

/* dynamics stuff */
if (data->gulp.run == MD)
  {
  switch (data->gulp.ensemble)
    {
    case NVE:
      fprintf(fp,"ensemble nve\n\n");
      break;
    case NVT:
      fprintf(fp,"ensemble nvt\n\n");
      break;
    case NPT:
      fprintf(fp,"ensemble npt\n\n");
      break;
    }

  fprintf(fp, "timestep %s\n", data->gulp.timestep);
  fprintf(fp, "equilibration %s\n", data->gulp.equilibration);
  fprintf(fp, "production %s\n", data->gulp.production);
  fprintf(fp, "sample %s\n", data->gulp.sample);
  fprintf(fp, "write %s\n", data->gulp.write);
  fprintf(fp, "\n");

  if (data->gulp.trj_file)
    fprintf(fp, "output trajectory %s\n", data->gulp.trj_file);
  }

if (data->gulp.mov_file)
  {
  fprintf(fp, "output movie arc %s\n", data->gulp.mov_file);
  if (data->gulp.run == MD)
    fprintf(fp, "mdarchive %s\n", data->gulp.mov_file);
  }


/* NEW - cosmo stuff */
if (data->gulp.solvation_model != GULP_SOLVATION_NONE)
  {
  switch (data->gulp.cosmo_shape)
    {
    case 0:
      fprintf(fp, "cosmoshape octahedron\n");
      break;
    default:
      fprintf(fp, "cosmoshape dodecahedron\n");
    }
  fprintf(fp, "pointsperatom %d\n", data->gulp.cosmo_points);
  fprintf(fp, "segmentsperatom %d\n", data->gulp.cosmo_segments);

  fprintf(fp, "solventepsilon %f\n", data->gulp.cosmo_solvent_epsilon);
  fprintf(fp, "solventradius %f %f\n", data->gulp.cosmo_solvent_radius, data->gulp.cosmo_solvent_delta);
  fprintf(fp, "solventrmax %f\n", data->gulp.cosmo_solvent_rmax);
  fprintf(fp, "rangeforsmooth %f\n", data->gulp.cosmo_smoothing);

  fprintf(fp, "\n");

  if (data->gulp.cosmo_sas)
    fprintf(fp, "output sas %s.sas\n", data->basename);
  }


#if DEBUG_WRITE_GULP
printf("  periodic: %d\n", data->periodic);
printf("fractional: %d\n", data->fractional);
#endif

switch(data->periodic)
  {
  case 3:
    if (data->construct_pbc)
      {
/* NB: gulp matrices are transposed wrt gdis */
      fprintf(fp,"vectors\n");
      fprintf(fp,"%15.10f %15.10f %15.10f\n",
                  data->latmat[0], data->latmat[3], data->latmat[6]);
      fprintf(fp,"%15.10f %15.10f %15.10f\n",
                  data->latmat[1], data->latmat[4], data->latmat[7]);
      fprintf(fp,"%15.10f %15.10f %15.10f\n",
                  data->latmat[2], data->latmat[5], data->latmat[8]);
      }
    else
      {
      fprintf(fp,"cell\n");
      fprintf(fp,"%lf %lf %lf  %lf %lf %lf\n",data->pbc[0],data->pbc[1],data->pbc[2],
                data->pbc[3]*180/PI, data->pbc[4]*180/PI, data->pbc[5]*180/PI);
      }
    break;

  case 2:
    if (data->construct_pbc)
      {
      fprintf(fp,"svectors\n");
      fprintf(fp,"%lf  %lf\n%lf  %lf\n",data->latmat[0], data->latmat[3],
                                    data->latmat[1], data->latmat[4]);
      }
    else
      {
      fprintf(fp,"scell\n");
      fprintf(fp,"%lf %lf %lf\n",data->pbc[0], data->pbc[1],
                              180.0*data->pbc[5]/PI);
      }
    break;

  case 1:
    if (data->construct_pbc)
      {
      fprintf(fp,"pvector\n");
      fprintf(fp,"%lf\n",data->latmat[0]);
      }
    else
      {
      fprintf(fp,"pcell\n");
      fprintf(fp,"%lf\n",data->pbc[0]);
      }
    break;
  }

/* coord section write */
flag = 0;

/* CURRENT - no more than 3 regions allowed */
for (region=0 ; region<3 ; region++)
  {
/* cores */
  for (list=data->cores ; list ; list=g_slist_next(list))
    {
    core = list->data;

    if (core->status & DELETED)
      continue;
    if (!core->primary)
      continue;
    if (core->region != region)
       continue;

/* write suitable header if it hasn't yet been written */
    if (!flag)
      {
      write_region_header(fp, region, data);
      flag++;
      }
    fprint_core(fp, core, data);
    }

/* shells */
  for (list=data->shels ; list ; list=g_slist_next(list))
    {
    shell = list->data;

    if (shell->status & DELETED)
      continue;
    if (!shell->primary)
      continue;
    if (shell->region != region)
      continue;

    fprint_shell(fp, shell, data);
    }
  flag = 0;
  }

/* post coord data */
switch(data->periodic)
  {
  case 2:
/* write dhkl for Eatt calcs */
    if (data->surface.dspacing > 0.1 && !data->gulp.no_eatt)
      fprintf(fp, "\ndhkl %lf\n",data->surface.dspacing);

/* write surface bulk energy for Esurf calcs */
    if (!data->gulp.no_esurf)
      fprintf(fp, "\nsbulkenergy %lf\n",data->gulp.sbulkenergy);
    break;

/* SG info */
  case 3:
    if (data->sginfo.spacename)
      {
      fprintf(fp,"\nspace\n%s\n",g_strstrip(data->sginfo.spacename));
      if (data->sginfo.cellchoice)
        fprintf(fp,"\norigin %d\n",data->sginfo.cellchoice);
      }
    else
      {
      fprintf(fp,"\nspace\nP 1\n");
      }
    break;
  }
fprintf(fp,"\n");

/* NEW - output connectivity */
if (data->gulp.noautobond)
  {
  for (list=data->bonds ; list ; list=g_slist_next(list))
    {
    bond = list->data;

    core1 = bond->atom1;
    core2 = bond->atom2;

    i = 1+g_slist_index(data->cores, core1);
    j = 1+g_slist_index(data->cores, core2);

/* connect i j <type> ix iy iz */
    fprintf(fp, "connect %4d %4d", i, j);

    switch(bond->type)
      {
      case BOND_TRIPLE:
        fprintf(fp, "  triple");
        break;
      case BOND_DOUBLE:
        fprintf(fp, "  double");
        break;
      case BOND_SINGLE:
      default:
        fprintf(fp, "  single");
        break;
      }

/* image <ix iy iz> : the translation of atom j needed to make the bond */
    ARR3SET(x, bond->offset);
    VEC3MUL(x, 2.0);
    ARR3ADD(x, core1->x);
    ARR3SUB(x, core2->x);
    fprintf(fp, "  %d %d %d\n", nearest_int(x[0]), nearest_int(x[1]), nearest_int(x[2]));

    fprintf(fp, "\n");
    }
  fprintf(fp, "\n");
  }

/* output species data */
if (data->gulp.species)
  {
  strip_extra(data->gulp.species);
  if (strlen(data->gulp.species))
    {
    fprintf(fp,"species\n");
    fprintf(fp,"%s\n", data->gulp.species);
    fprintf(fp,"end\n\n");
    }
  }

/* output cova data */
/* FIXME - what if other element properties have been modified (within GDIS) */
if (data->gulp.elements)
  {
  strip_extra(data->gulp.elements);
  if (strlen(data->gulp.elements))
    {
    fprintf(fp,"elements\n");
    fprintf(fp,"%s\n", data->gulp.elements);
    fprintf(fp,"end\n\n");
    }
  }

/* potentials */
if (data->gulp.potentials)
  {
  if (strlen(data->gulp.potentials))
    fprintf(fp,"%s\n", data->gulp.potentials);
  }
else
  ff_gulp_write(fp, data);

/* electrostatic potential sites */
if (data->gulp.epot_vecs)
  {
  fprintf(fp, "potsites cart\n");

  for (list=data->gulp.epot_vecs ; list ; list=g_slist_next(list))
    {
    vec = list->data;

    fprintf(fp, "%f  %f  %f\n", vec->rx[0], vec->rx[1], vec->rx[2]);
    }
  fprintf(fp,"\n");
  }

/* specified kpoints? */
if (data->periodic && data->gulp.kpoints)
  {
  strip_extra(data->gulp.kpoints);
  if (strlen(data->gulp.kpoints))
    {
    fprintf(fp, "kpoints %d\n", 1+char_count(data->gulp.kpoints, '\n'));
    fprintf(fp, "%s\n\n", data->gulp.kpoints);
    }
  }

/* unrecognized stuff */
if (data->gulp.output_extra)
  {
  if (data->gulp.extra)
    {
    fprintf(fp, "\n#--- Additional options/data\n\n");
    fprintf(fp, "%s\n\n", data->gulp.extra);
    }
  }
fprintf(fp, "print 1\n");

/* done */
fclose(fp);
g_string_free(line, TRUE);
return(0);
}

/*****************/
/* debugging aid */
/*****************/
void print_gulp_flags(gint *list, gint size)
{
gint i;

printf("Flags:");
for (i=0 ; i<size ; i++)
  {
/* active? */
  if (*(list+i))
    {
    switch(i)
      {
      case NAME:
        printf(" name");
        break;
      case SURFACE_CELL:
        printf(" scell");
        break;
      case CELL:
        printf(" cell");
        break;
      case FRAC:
        printf(" frac");
        break;
      case SFRAC:
        printf(" sfrac");
        break;
      case PFRAC:
        printf(" pfrac");
        break;
      case CART:
        printf(" cart");
        break;
      case DUMP_FILE:
        printf(" dump");
        break;
      case SURFACE_VECTORS:
        printf(" svec");
        break;
      case LATTICE_VECTORS:
        printf(" lvec");
        break;
      }
    }
  }
printf("\n");
}

/***********************/
/* gulp keyword parser */
/***********************/
gint gulp_is_keyword(gchar *text)
{
gint i;
size_t n;

if (!text)
return(-1);

i=0;
while (gulp_keyword[i].label)
  {
  n = strlen(gulp_keyword[i].label);
  if (g_ascii_strncasecmp(text, gulp_keyword[i].label, n) == 0)
    return(gulp_keyword[i].code);
  i++;
  }
return(-1);
}

/***********************/
/* gulp options parser */
/***********************/
gint gulp_is_option(gchar *text)
{
gint i;
size_t n;

if (!text)
  return(-1);

i=0;
while (gulp_option[i].label)
  {
  n = strlen(gulp_option[i].label);
  if (g_ascii_strncasecmp(text, gulp_option[i].label, n) == 0)
    return(gulp_option[i].code);
  i++;
  }
return(-1);
}

enum
{
GULP_NEW_NAME, GULP_NEW_LATTICE, GULP_NEW_COORDS, GULP_NEW_FLAGS
};

/************************************************/
/* determines if we should allocate a new model */
/************************************************/
/* FIXME - need a better way of handling this so that keyword line directives are inherited */
struct model_pak *gulp_change_model(gint trigger, struct model_pak *model)
{
gint i, flag_new=FALSE;
static gint flag_level=1, flag_table[GULP_NEW_FLAGS];
struct model_pak *new;

if (!model)
  {
/* init flag table */
  flag_level = 1;
  for (i=GULP_NEW_FLAGS ; i-- ; )
    flag_table[i] = 0;
  return(NULL);
  }

/* handle possible new model triggers (NB: could occur in any order) */
switch (trigger)
  {
  case GULP_NEW_COORDS:
  case GULP_NEW_NAME:
  case GULP_NEW_LATTICE:
    flag_table[trigger]++;
/* only new flag triggers above the current level trigger a new model */
    if (flag_table[trigger] > flag_level)
      {
      flag_new = TRUE;
      flag_level = flag_table[trigger];
      }
    break;
  }

/* once the new coords flag is tripped - any subsequent new flag */
/* triggers should generate a new model. */
/* NB: name/cell options must come before the coords (GULP incompatibility?) */
if (trigger == GULP_NEW_COORDS)
  {
  for (i=GULP_NEW_FLAGS ; i-- ; )
    flag_table[i] = flag_level;
  }

/*
printf("%d : %d : %d\n", flag_table[0], flag_table[1], flag_table[2]); 
*/

/* return new model if flagged, else the old model */
if (flag_new)
  {
/*
printf(" + new model triggered!\n");
*/
  new = model_new();
  new->id = -1;
  return(new);
  }
return(model);
}

/*********************************/
/* yet another gulp read rewrite */
/*********************************/
#define DEBUG_READ_GULP 0
gint read_gulp(gchar *filename, struct model_pak *model)
{
gint i, j, nc, ns;
gint keywords_found, keywords_expected, num_tokens, remaining;
gint code, context, run, method, optimiser, optimiser2, unit, fix, switch_type;
gint region, breathe, coulomb, maxcyc, qeq, free, zsisa, compare, nosym, phonon, eigen;
gint type, gen_flag, temp_code, noflags, growth, translate, noautobond;
gchar *tmp, *line, **buff;
gdouble switch_value;
gdouble x[3], vec1[3], vec2[3], vec3[3];
gdouble strainmat[9], svec[4];
gpointer scan;
GSList *item, *clist, *slist;
GString *key_buff, *elem_buff, *ex_buff, *ff_buff, *species_buff, *kpoints_buff;
struct gulp_pak gulp;
struct elem_pak elem_data;
struct core_pak *core, *core1, *core2;
struct shel_pak *shell;
struct bond_pak *bond;

scan = scan_new(filename);
if (!scan)
  return(1);

/* init */
model->id = -1;
keywords_found = keywords_expected = 0;
/* deprec - use gulp_pak structure/gulp_init() instead */
/*code = */context = run = method = optimiser = optimiser2 = switch_type = switch_value = -1;/*FIX 289935*/
qeq = free = zsisa = compare = nosym = phonon = eigen = unit = fix = FALSE;

/* NEW*/
noautobond = FALSE;

gulp_init(&gulp);

/* NEW - noflags assumed false now */
noflags = FALSE;
maxcyc = 0;
region = REGION1A;
coulomb = NOBUILD;
gulp_change_model(0, NULL);
key_buff = g_string_new(NULL);
elem_buff = g_string_new(NULL);
ex_buff = g_string_new(NULL);
ff_buff = g_string_new(NULL);
species_buff = g_string_new(NULL);
kpoints_buff = g_string_new(NULL);

/* main parse loop */
buff = scan_get_tokens(scan, &num_tokens);
while (!scan_complete(scan))
  {

g_assert(buff != NULL);

  if (!keywords_found)
    {
/* keyword only search */
    for (i=num_tokens ; i-- ; )
      {
      code = gulp_is_keyword(*(buff+i));
      switch (code)
        {
        case E_SINGLE:
        case E_OPTIMIZE:
        case MD:
          run = code;
          keywords_found++;
          keywords_expected = num_tokens;
          break;
        case RFO_OPT:
        case BFGS_OPT:
        case CONJ_OPT:
          optimiser = code;
          keywords_found++;
          keywords_expected = num_tokens;
          break;
        case UNIT_HESSIAN:
          unit = TRUE;
          keywords_found++;
          keywords_expected = num_tokens;
          break;

        case MOLE:
        case MOLMEC:
        case MOLQ:
          coulomb = code;
          keywords_found++;
          keywords_expected = num_tokens;
          break;
        case FIX:
          fix = TRUE;
          keywords_found++;
          break;
        case GULP_NOAUTOBOND:
/* set flag incase we have more than one structure in the file */
          noautobond = TRUE;
          keywords_found++;
          break;

        case QEQ:
          qeq = TRUE;
          keywords_found++;
          keywords_expected = num_tokens;
          break;
        case FREE_ENERGY:
          free = TRUE;
          keywords_found++;
          keywords_expected = num_tokens;
          break;
        case ZSISA:
          zsisa = TRUE;
          keywords_found++;
          keywords_expected = num_tokens;
          break;
        case COMPARE:
          compare = TRUE;
          keywords_found++;
          keywords_expected = num_tokens;
          break;
        case NOSYM:
          nosym = TRUE;
          keywords_found++;
          keywords_expected = num_tokens;
          break;
        case CONP:
        case CONV:
          method = code;
          keywords_found++;
          keywords_expected = num_tokens;
          noflags = TRUE;
          break;

        case CELL:
        case NOFLAGS:
          noflags = TRUE;
          break;

        case GULP_SOLVATION_COSMIC:
        case GULP_SOLVATION_COSMO:
          gulp.solvation_model = code;
          keywords_found++;
          keywords_expected = num_tokens;
          break;

        case PHONON:
          phonon = TRUE;
          keywords_found++;
          keywords_expected = num_tokens;
          break;
        case EIGEN:
          eigen = TRUE;
          keywords_found++;
          keywords_expected = num_tokens;
          break;

        default:
          g_string_sprintfa(key_buff, "%s ", *(buff+i));

        }
      }

/* FIXME - check if noflags should be made false */

    }
  else
    {
/* main option/setup parse */
/* acquire next option context? */
    code = gulp_is_option(*buff);
    if (code != -1)
      context = code;

/* pass everything we don't recognize to the extra string */
if (code == -1)
  {
/* FIXME - should have finished all contextual stuff */
  if (context != -1)
    {
#if DEBUG_READ_GULP
printf("Warning: context = %d\n", context);
printf(" > %s", scan_cur_line(scan));
#endif
    }
  g_string_sprintfa(ex_buff, "%s", scan_cur_line(scan));
  }

/*
printf("code: %d, context: %d\n", code, context);
*/
/* parse current option */
      switch (code)
        {
        case GULP_CONNECT:

          if (num_tokens < 6)
            {
/* FIXME - why isnt the error table working? */
            error_table_entry("file_gulp(): bad connect line.\n");
            break;
            }

          i = str_to_float(*(buff+1));
          j = str_to_float(*(buff+2));
          nc = g_slist_length(model->cores);

/* TODO - print warning (bad connect entry) */
          if (i > nc || j > nc)
            break;

          if (num_tokens == 6)
            {
/* exactly 6 -> assume no bond type specified */
            vec1[0] = str_to_float(*(buff+3));
            vec1[1] = str_to_float(*(buff+4));
            vec1[2] = str_to_float(*(buff+5));
            } 
          else
            {
/* more than 6 -> assume bond type is specified */
            vec1[0] = str_to_float(*(buff+4));
            vec1[1] = str_to_float(*(buff+5));
            vec1[2] = str_to_float(*(buff+6));
            }

/* reverse order (cores prepended) */
          core1 = g_slist_nth_data(model->cores, nc - i);
          core2 = g_slist_nth_data(model->cores, nc - j);

          g_assert(core1 != NULL);
          g_assert(core2 != NULL);

          bond = g_malloc(sizeof(struct bond_pak));

          model->bonds = g_slist_prepend(model->bonds, bond);
          core1->bonds = g_slist_prepend(core1->bonds, bond);
          core2->bonds = g_slist_prepend(core2->bonds, bond);

          bond->type = BOND_SINGLE;
          bond->status = NORMAL;
          bond->atom1 = core1;
          bond->atom2 = core2;

          ARR3SET(bond->offset, core2->x);
          ARR3SUB(bond->offset, core1->x);

/* bond midpoint offset must be in fractionals */
          if (!model->fractional)
            {
            matrix_lattice_init(model);
            vecmat(model->ilatmat, bond->offset);
            }

/* compute the bond offset (midpoint vector in fractionals) */
          ARR3ADD(bond->offset, vec1);
          VEC3MUL(bond->offset, 0.5);
          break;

        case GULP_COSMO_SHAPE:
          if (num_tokens > 1)
            {
            if (g_ascii_strncasecmp(*(buff+1), "dode", 4) == 0)
              gulp.cosmo_shape = 1;
            else
              gulp.cosmo_shape = 0;
            }
          break;
        case GULP_COSMO_POINTS:
          if (num_tokens > 1)
            gulp.cosmo_points = str_to_float(*(buff+1));
          break;
        case GULP_COSMO_SEGMENTS:
          if (num_tokens > 1)
            gulp.cosmo_segments = str_to_float(*(buff+1));
          break;
        case GULP_COSMO_SMOOTHING:
          if (num_tokens > 1)
            gulp.cosmo_smoothing = str_to_float(*(buff+1));
          break;
        case GULP_COSMO_SOLVENT_EPSILON:
          if (num_tokens > 1)
            gulp.cosmo_solvent_epsilon = str_to_float(*(buff+1));
          break;
        case GULP_COSMO_SOLVENT_RADIUS:
          if (num_tokens > 1)
            gulp.cosmo_solvent_radius = str_to_float(*(buff+1));
          if (num_tokens > 2)
            gulp.cosmo_solvent_delta = str_to_float(*(buff+2));
          break;
        case GULP_COSMO_SOLVENT_RMAX:
          if (num_tokens > 1)
            gulp.cosmo_solvent_rmax = str_to_float(*(buff+1));
          break;

        case ENDFORCE:
          g_string_sprintfa(ex_buff, "%s", scan_cur_line(scan));
          break;

        case PRINT:
/* skip, don't store the print option as gdis always appends it */
          if (num_tokens < 2)
            {
            while (!scan_complete(scan))
              {
              g_strfreev(buff);
              buff = scan_get_tokens(scan, &num_tokens);
              if (num_tokens)
                break;
              }
            }
          break;

        case MAXCYC:
          if (num_tokens > 1)
            maxcyc = str_to_float(*(buff+1));
          else
            {
            while (!scan_complete(scan))
              {
              g_strfreev(buff);
              buff = scan_get_tokens(scan, &num_tokens);
              if (num_tokens)
                {
                maxcyc = str_to_float(*buff);
                break;
                }
              }
            }
          break;

        case SWITCH_ON:
          for (j=1 ; j<num_tokens ; j++)
            {
            code = gulp_is_option(*(buff+j));
            switch (code)
              {
              case BFGS_OPT:
              case CONJ_OPT:
              case RFO_OPT:
                optimiser2 = code;
                break;

              case GNORM:
              case CYCLE:
                switch_type = code;
                break;
 
              default:
                if (str_is_float(*(buff+j)))
                  switch_value = str_to_float(*(buff+j));
              }
            }
          break;

/* NB: these keywords can initiate a new model (ie if repeat occurs) */
        case PFRAC:
        case SFRAC:
        case FRAC:
        case CART:
          if (num_tokens > 2)
            {
            if (g_ascii_strncasecmp("region",*(buff+1),6) == 0)
              region = str_to_float(*(buff+2)) - 1;
            model->region_empty[region] = FALSE; 
            if (region == REGION1A)
              model = gulp_change_model(GULP_NEW_COORDS, model);
            }
          else
            model = gulp_change_model(GULP_NEW_COORDS, model);
          nc = ns = 0;
          break;

/* another keyword that can initiate a new model (on repeat) */
        case NAME:
          model = gulp_change_model(GULP_NEW_NAME, model);
          break;

/* other keywords that can initiate a new model (on repeat) */
        case CELL:
        case POLYMER_CELL:
        case POLYMER_VECTOR:
        case SURFACE_CELL:
        case SURFACE_VECTORS:
        case LATTICE_VECTORS:
          model = gulp_change_model(GULP_NEW_LATTICE, model);
          break;
        }
    }

/* main option/data parse */
/* NB: after each case - set context = -1 if parsing for that option is completed */
    switch (context)
      {
/* pass through everything in these blocks */
      case TITLE:
      case WEIGHT:
      case GULP_OBSERVABLES:
      case GULP_VARIABLES:
        g_string_sprintfa(ex_buff, "%s", scan_cur_line(scan));
        while (!scan_complete(scan))
          {
          g_strfreev(buff);
          buff = scan_get_tokens(scan, &num_tokens);
          if (num_tokens)
            {
/* check first item only */
            code = gulp_is_option(*(buff));
            if (code == -1)
              g_string_sprintfa(ex_buff, "%s", scan_cur_line(scan));
            else
              {
/* terminate when we hit a new option */
/* NB: if option was not an "end" - buffer it for the next parse cycle */
              if (code != END)
                {
                scan_put_line(scan);
                }
              else
                g_string_sprintfa(ex_buff, "%s", scan_cur_line(scan));
              break;
              }
            }
          }
        break;

/* hack for getting an associated .trg file */
      case OUTPUT:
        switch (num_tokens)
          {
          case 4:
            if (g_ascii_strncasecmp(*(buff+1),"mov",3) == 0)
              {
              if (model->gulp.mov_file)
                g_free(model->gulp.mov_file);
              model->gulp.mov_file = g_strdup(*(buff+3));
              }
            break;
          case 3:
            if (g_ascii_strncasecmp(*(buff+1),"traj",4) == 0)
              {
              g_free(model->gulp.trj_file);
              model->gulp.trj_file = g_strdup(*(buff+2));
              model->animation = TRUE;
              }
            break;
          }
        context = -1;
        break;

        case ENSEMBLE:
          if (num_tokens > 1)
            {
            if (g_ascii_strncasecmp(*(buff+1),"nve",3) == 0)
              model->gulp.ensemble = NVE;
            if (g_ascii_strncasecmp(*(buff+1),"nvt",3) == 0)
              model->gulp.ensemble = NVT;
            if (g_ascii_strncasecmp(*(buff+1),"npt",3) == 0)
              model->gulp.ensemble = NPT;
            }
          context = -1;
          break;

        case TEMPERATURE:
          if (num_tokens > 1)
            {
            g_free(model->gulp.temperature);
            model->gulp.temperature = g_strjoinv(" ", buff+1);
            }
          else
            {
/* get next (non empty) line */
            g_strfreev(buff);
            buff = scan_get_tokens(scan, &num_tokens);
            if (buff)
              {
              g_free(model->gulp.temperature);
              model->gulp.temperature = g_strjoinv(" ", buff);
              }
            }
          context = -1;
          break;

        case PRESSURE:
          if (num_tokens > 1)
            {
            g_free(model->gulp.pressure);
            model->gulp.pressure = g_strjoinv(" ", buff+1);
            }
          else
            {
/* get next (non empty) line */
            g_strfreev(buff);
            buff = scan_get_tokens(scan, &num_tokens);
            if (buff)
              {
              g_free(model->gulp.pressure);
              model->gulp.pressure = g_strjoinv(" ", buff);
              }
            }
          context = -1;
          break;

       case TIMESTEP:
         g_free(model->gulp.timestep);
         model->gulp.timestep = g_strdup(get_token_pos(scan_cur_line(scan), 1));
         g_strstrip(model->gulp.timestep);
         break;
       case EQUILIBRATION:
         g_free(model->gulp.equilibration);
         model->gulp.equilibration = g_strdup(get_token_pos(scan_cur_line(scan), 1));
         g_strstrip(model->gulp.equilibration);
         break;
       case PRODUCTION:
         g_free(model->gulp.production);
         model->gulp.production = g_strdup(get_token_pos(scan_cur_line(scan), 1));
         g_strstrip(model->gulp.production);
         break;
       case SAMPLE:
         g_free(model->gulp.sample);
         model->gulp.sample = g_strdup(get_token_pos(scan_cur_line(scan), 1));
         g_strstrip(model->gulp.sample);
         break;
       case WRITE:
         g_free(model->gulp.write);
         model->gulp.write = g_strdup(get_token_pos(scan_cur_line(scan), 1));
         g_strstrip(model->gulp.write);
         break;

        case SUPER_CELL:
          if (num_tokens > 3)
            {
            model->gulp.super[0] = str_to_float(*(buff+1));
            model->gulp.super[1] = str_to_float(*(buff+2));
            model->gulp.super[2] = str_to_float(*(buff+3));
            }
          break;

        case NAME:
          if (num_tokens > 1)
            {
            g_free(model->basename);
            model->basename = g_strjoinv(" ", buff+1);
            g_strstrip(model->basename);
            strcpy(model->filename, model->basename);
            }
          break;

        case GULP_LIBRARY:
          if (num_tokens > 1)
            {
            g_free(model->gulp.libfile);
            model->gulp.libfile = g_strdup(*(buff+1));
            }
          break;

        case SPACE:
          if (num_tokens > 1)
            {
            line = scan_cur_line(scan);
            tmp = get_token_pos(line, 1);
            }
          else
            tmp = scan_get_line(scan);
 
/* process for spacgroup symbol or number */
          gen_flag=0;
          for (i=0 ; i<strlen(tmp) ; i++)
            {
/* any alphabetic chars => space *name* */
            if (g_ascii_isalpha(*(tmp+i)))
              {
              model->sginfo.spacename = g_strdup(tmp);
/* indicate that name should used in lookup */
              model->sginfo.spacenum = -1;
              gen_flag++;
              break;
              }
            }
/* no alphabetic chars, assume space group number */
          if (!gen_flag)
            model->sginfo.spacenum = str_to_float(tmp);
          context = -1;
          break;

        case ORIGIN:
          switch (num_tokens)
            {
            case 2:
              model->sginfo.cellchoice = str_to_float(*(buff+1));
              break;
            default:
              printf("Warning: unsupported origin specification.\n");
            }
          break;

        case POLYMER_CELL:
          g_strfreev(buff);
          buff = scan_get_tokens(scan, &num_tokens);
          if (num_tokens)
            {
            model->pbc[0] = str_to_float(*buff);
            model->periodic = 1;
            }
          break;

        case POLYMER_VECTOR:
          g_strfreev(buff);
          buff = scan_get_tokens(scan, &num_tokens);
          if (num_tokens)
            {
            model->latmat[0] = str_to_float(*(buff+0)); 
            model->latmat[3] = 0.0;
            model->latmat[6] = 0.0;
            model->construct_pbc = TRUE;
            model->periodic = 1;
            matrix_lattice_init(model);
            g_strfreev(buff);
/* check for optimisation flags */
/* At the moment just skips so input file can be read successfully */
            buff = scan_get_tokens(scan, &num_tokens);
            if ((buff[0][0] != '0' || buff[0][0] != '1' ) && num_tokens != 1)
              scan_put_line(scan);
            }
          break;

        case SURFACE_CELL:
          g_strfreev(buff);
          buff = scan_get_tokens(scan, &num_tokens);
          if (num_tokens > 2)
            {
            model->pbc[0] = str_to_float(*(buff+0));
            model->pbc[1] = str_to_float(*(buff+1));
            model->pbc[5] = str_to_float(*(buff+2));
            model->pbc[5] *= D2R;
            model->periodic = 2;
            }
          break;

        case SURFACE_VECTORS:
          g_strfreev(buff);
          buff = scan_get_tokens(scan, &num_tokens);
          if (num_tokens > 1)
            {
/* FIXME - a is not necessarily colinear with x */
            vec1[0] = str_to_float(*(buff+0));
            vec2[0] = str_to_float(*(buff+1));
            g_strfreev(buff);
            buff = scan_get_tokens(scan, &num_tokens);
            if (num_tokens > 1)
              {
              vec1[1] = str_to_float(*(buff+0));
              vec2[1] = str_to_float(*(buff+1));
/* NB: gdis wants transposed gulp matrix */
              VEC3SET(&(model->latmat[0]), vec1[0], vec1[1], 0.0);
              VEC3SET(&(model->latmat[3]), vec2[0], vec2[1], 0.0);
              VEC3SET(&(model->latmat[6]), 0.0, 0.0, 1.0);
              model->construct_pbc = TRUE;
              model->periodic = 2;
              matrix_lattice_init(model);
              }
            g_strfreev(buff);
/* check for optimisation flags */
/* At the moment just skips so input file can be read successfully */
            buff = scan_get_tokens(scan, &num_tokens);
            if ((buff[0][0] != '0' || buff[0][0] != '1' ) && num_tokens != 2)
              scan_put_line(scan);
            }
         break;

        case CELL:
/* get next line if insufficient tokens */
          if (num_tokens < 2)
            {
            g_strfreev(buff);
            buff = scan_get_tokens(scan, &num_tokens);
            i=0;
            }
          else
            i=1;
/* parse for cell data */
          if (num_tokens > 5)
            {
            model->pbc[0] = str_to_float(*(buff+i+0));
            model->pbc[1] = str_to_float(*(buff+i+1));
            model->pbc[2] = str_to_float(*(buff+i+2));
            model->pbc[3] = str_to_float(*(buff+i+3));
            model->pbc[4] = str_to_float(*(buff+i+4));
            model->pbc[5] = str_to_float(*(buff+i+5));
/* hack to determine if 2D or 3D */
/* FIXME - possible for pbc[3] or pbc[4] to be 0 instead */
            if (fabs(model->pbc[2]) < FRACTION_TOLERANCE)
              model->periodic = 2;
            else
              model->periodic = 3;
            model->pbc[3] *= D2R;
            model->pbc[4] *= D2R;
            model->pbc[5] *= D2R;
            }
          break;

        case LATTICE_VECTORS:
          VEC3SET(vec1, 1.0, 0.0, 0.0);
          VEC3SET(vec2, 0.0, 1.0, 0.0);
          VEC3SET(vec3, 0.0, 0.0, 1.0);
          model->periodic = 0;
          g_strfreev(buff);
          buff = scan_get_tokens(scan, &num_tokens);
          if (num_tokens > 2)
            {
            vec1[0] = str_to_float(*(buff+0));
            vec2[0] = str_to_float(*(buff+1));
            vec3[0] = str_to_float(*(buff+2));
            model->periodic++;
            }
          g_strfreev(buff);
          buff = scan_get_tokens(scan, &num_tokens);
          if (num_tokens > 2)
            {
            vec1[1] = str_to_float(*(buff+0));
            vec2[1] = str_to_float(*(buff+1));
            vec3[1] = str_to_float(*(buff+2));
            model->periodic++;
            }
          g_strfreev(buff);
          buff = scan_get_tokens(scan, &num_tokens);
          if (num_tokens > 2)
            {
            vec1[2] = str_to_float(*(buff+0));
            vec2[2] = str_to_float(*(buff+1));
            vec3[2] = str_to_float(*(buff+2));
            model->periodic++;
            }
          g_strfreev(buff);
          if (model->periodic == 3)
            {
/* use the supplied lattice matrix */
            ARR3SET(&(model->latmat[0]), vec1);
            ARR3SET(&(model->latmat[3]), vec2);
            ARR3SET(&(model->latmat[6]), vec3);
            model->construct_pbc = TRUE;
            matrix_lattice_init(model);
            }
/* check for optimisation flags */
/* At the moment just skips so input file can be read successfully */
          buff = scan_get_tokens(scan, &num_tokens);
          if ((buff[0][0] != '0' || buff[0][0] != '1' ) && num_tokens != 6)
            scan_put_line(scan);
          break;

        case CELL_STRAIN:
          if (model->periodic == 3 && num_tokens == 7)
            {
            model->cell_strains[0] = str_to_float(*(buff+1));
            model->cell_strains[1] = str_to_float(*(buff+2));
            model->cell_strains[2] = str_to_float(*(buff+3));
            model->cell_strains[3] = str_to_float(*(buff+4));
            model->cell_strains[4] = str_to_float(*(buff+5));
            model->cell_strains[5] = str_to_float(*(buff+6));
              
            strainmat[0] = model->cell_strains[0] + 1.0;
            strainmat[4] = model->cell_strains[1] + 1.0;
            strainmat[8] = model->cell_strains[2] + 1.0;
            strainmat[1] = strainmat[3] = model->cell_strains[5]*0.5;
            strainmat[2] = strainmat[6] = model->cell_strains[4]*0.5;
            strainmat[5] = strainmat[7] = model->cell_strains[3]*0.5;
#if DEBUG_READ_GULP
P3MAT("lattice matrix: ", model->latmat);
P3MAT("strain matrix: ", strainmat);
#endif
            matmat(strainmat, model->latmat);
#if DEBUG_READ_GULP
P3MAT("new lattice matrix: ", model->latmat);
#endif
            }
          else if (model->periodic == 2 && num_tokens == 4)
            {
            model->cell_strains[0] = str_to_float(*(buff+1));
            model->cell_strains[1] = str_to_float(*(buff+2));
            model->cell_strains[2] = str_to_float(*(buff+3));
              
            strainmat[0] = model->cell_strains[0] + 1.0;
            strainmat[3] = model->cell_strains[1] + 1.0;
            strainmat[1] = strainmat[2] = model->cell_strains[2]*0.5;
/* GDIS stores the latmat as a 3x3 matrix; pullout 2 x 2 do matrix ops and put back */
            svec[0] = model->latmat[0];
            svec[1] = model->latmat[1];
            svec[2] = model->latmat[3];
            svec[3] = model->latmat[4];
#if DEBUG_READ_GULP
P2MAT("lattice matrix: ", svec);
P2MAT("strain matrix: ", strainmat);
#endif
            mat2mat(strainmat, svec);
#if DEBUG_READ_GULP
P2MAT("new lattice matrix: ", svec);
#endif
            model->latmat[0] = svec[0];
            model->latmat[1] = svec[1];
            model->latmat[3] = svec[2];
            model->latmat[4] = svec[3];
            }
          else if (model->periodic == 1 && num_tokens == 2)
            {
            model->cell_strains[0] = str_to_float(*(buff+1));
              
            strainmat[0] = model->cell_strains[0] + 1.0;
            model->latmat[0] *= strainmat[0];
            }
          else
            {
            gui_text_show(WARNING, "Invalid use of cstrain. Ignoring...\n");
            }
          break;
/* NB: errors -> break (not return) so we cope with empty regions */
        case PFRAC:
        case SFRAC:
        case FRAC:
          model->fractional = TRUE;
        case CART:
/* au check */
          gen_flag = 0;
          if (num_tokens > 1)
            if (g_ascii_strncasecmp("au",*(buff+1),2) == 0)
              gen_flag++;

/* get new line */
/* insufficient items => stupid gulp frac/cart \n number \n coords format */
          g_strfreev(buff);
          buff = scan_get_tokens(scan, &num_tokens);
          if (num_tokens < 4)
            {
            g_strfreev(buff);
            buff = scan_get_tokens(scan, &num_tokens);
            }

/* core, shell indices */
          while (!scan_complete(scan))
            {
            if (gulp_is_option(*buff) > -1)
              {
              scan_put_line(scan);
              break;
              }

/* NEW - translation and growth slice check */
            growth = translate = FALSE;
            switch (*(*(buff+num_tokens-1)))
              {
              case '%':
                growth = TRUE;
                break;

              case 'T':
                translate = TRUE;
                break;
              }

/* sanity check for coordinate line, atom + x,y,z */
            if (num_tokens < 4)
              break;

/* if an atom type (ie core/shell) is specified - adjust data positions */
            type = 0;
/* next data token position */
            i = 1;
            breathe = FALSE;

            if (*(buff+1)[0] == 'b' || *(buff+1)[0] == 'B')
              {
              breathe = TRUE;
              i++;
              if (*(buff+1)[1] == 'c' || *(buff+1)[1] == 'C')
                type = 1;
              else
                type = 2;
              }
            else
              {
              if (*(buff+1)[0] == 'c' || *(buff+1)[0] == 'C')
                {
                type = 1;
                i++;
                }
              if (*(buff+1)[0] == 's' || *(buff+1)[0] == 'S')
                {
                type = 2;
                i++;
                }
              }

/* x y z for core and shell */
            x[0] = str_to_float(*(buff+i));
            x[1] = str_to_float(*(buff+i+1));
            x[2] = str_to_float(*(buff+i+2));
            if (gen_flag)
              {
              VEC3MUL(x, AU2ANG);
              }
            i+=3;

            if (growth || translate)
              i++; 

/* number of items after x y z */
            remaining = num_tokens - i;
/* do this so the switch still works (ie ignores anything more than what it expects) */
            remaining = CLAMP(remaining, 0, 6);

            if (type == 2)
              {
/* new shell */
              shell = new_shell(*buff, model);
              ARR3SET(shell->x, x);

              shell->breathe = breathe;
              shell->region = region;
              shell->translate = translate;

              switch (remaining)
                {
                case 6:
                  shell->lookup_charge = FALSE;
                  shell->charge = str_to_float(*(buff+i));
                  shell->has_sof = TRUE;
                  shell->sof = str_to_float(*(buff+i+1));
                  shell->radius = str_to_float(*(buff+i+2));
                  shell->flags = g_strdup_printf("%s %s %s",
                                  *(buff+i+3), *(buff+i+4), *(buff+i+5));
                  break;

                case 5:
                  shell->lookup_charge = FALSE;
                  shell->charge = str_to_float(*(buff+i));
                  shell->has_sof = TRUE;
                  shell->sof = str_to_float(*(buff+i+1));
                  shell->flags = g_strdup_printf("%s %s %s",
                                  *(buff+i+2), *(buff+i+3), *(buff+i+4));
                  break;

                case 4:
                  shell->lookup_charge = FALSE;
                  shell->charge = str_to_float(*(buff+i));

                  shell->flags = g_strdup_printf("%s %s %s",
                                  *(buff+i+1), *(buff+i+2), *(buff+i+3));
                  break;

                case 3:
                  if (noflags)
                    {
                    shell->lookup_charge = FALSE;
                    shell->charge = str_to_float(*(buff+i));
                    shell->has_sof = TRUE;
                    shell->sof = str_to_float(*(buff+i+1));
                    shell->radius = str_to_float(*(buff+i+2));
                    }
                  else
                    shell->flags = g_strdup_printf("%s %s %s",
                                    *(buff+i), *(buff+i+1), *(buff+i+2));
                  break;

                case 2:
                  shell->has_sof = TRUE;
                  shell->sof = str_to_float(*(buff+i+1));
                case 1:
                  shell->lookup_charge = FALSE;
                  shell->charge = str_to_float(*(buff+i));
                  break;
                }
              model->shels = g_slist_prepend(model->shels, shell);
              ns++;
              }
            else
              {
/* new core */
              core = new_core(*buff, model);
              ARR3SET(core->x, x);

              core->breathe = breathe;
              core->region = region;
              core->growth = growth;
              core->translate = translate;

              switch (remaining)
                {
                case 6:
                  core->lookup_charge = FALSE;
                  core->charge = str_to_float(*(buff+i));
                  core->has_sof = TRUE;
                  core->sof = str_to_float(*(buff+i+1));
                  core->radius = str_to_float(*(buff+i+2));
                  core->flags = g_strdup_printf("%s %s %s",
                                  *(buff+i+3), *(buff+i+4), *(buff+i+5));
                  break;

                case 5:
                  core->lookup_charge = FALSE;
                  core->charge = str_to_float(*(buff+i));
                  core->has_sof = TRUE;
                  core->sof = str_to_float(*(buff+i+1));
                  core->flags = g_strdup_printf("%s %s %s",
                                  *(buff+i+2), *(buff+i+3), *(buff+i+4));
                  break;

                case 4:
                  core->lookup_charge = FALSE;
                  core->charge = str_to_float(*(buff+i));

                  core->flags = g_strdup_printf("%s %s %s",
                                  *(buff+i+1), *(buff+i+2), *(buff+i+3));
                  break;

                case 3:
                  if (noflags)
                    {
                    core->lookup_charge = FALSE;
                    core->charge = str_to_float(*(buff+i));
                    core->has_sof = TRUE;
                    core->sof = str_to_float(*(buff+i+1));
                    core->radius = str_to_float(*(buff+i+2));
                    }
                  else
                    core->flags = g_strdup_printf("%s %s %s",
                                    *(buff+i), *(buff+i+1), *(buff+i+2));
                  break;

                case 2:
                  core->has_sof = TRUE;
                  core->sof = str_to_float(*(buff+i+1));
                case 1:
                  core->lookup_charge = FALSE;
                  core->charge = str_to_float(*(buff+i));
                  break;
                }

              model->cores = g_slist_prepend(model->cores, core);
              nc++;
              }



#if OLD_WAY
                if (breathe)
                  {
                  if (num_tokens > 5)
                    core->radius = str_to_float(*(buff+5));
/* FIXME - what about charge? */
                  }
                else
                  {
                  i = 4+type; 
                  remaining = num_tokens - i;

                  remaining = CLAMP(remaining, 0, 5);

                  if (noflags)
                    remaining = CLAMP(remaining, 0, 2); 

                  switch (remaining)
                    {
                    case 5:
                      core->lookup_charge = FALSE;
                      core->charge = str_to_float(*(buff+i));
                      core->has_sof = TRUE;
                      core->sof = str_to_float(*(buff+i+1));
                      core->flags = g_strdup_printf("%s %s %s", *(buff+i+2), *(buff+i+3), *(buff+i+4));
                      break;

                    case 4:
                      core->lookup_charge = FALSE;
                      core->charge = str_to_float(*(buff+i));
                      i++;
                    case 3:
                      core->flags = g_strdup_printf("%s %s %s", *(buff+i), *(buff+i+1), *(buff+i+2));
                      break;

                    case 2:
                      core->has_sof = TRUE;
                      core->sof = str_to_float(*(buff+i+1));
                    case 1:
                      core->lookup_charge = FALSE;
                      core->charge = str_to_float(*(buff+i));
                      break;
                    }
                  }



                if (gen_flag)
                  {
                  VEC3MUL(core->x, AU2ANG);
                  }
                nc++;
                break;

              case 2:
                if (num_tokens > 4)
                  {
                  shell = new_shell(*buff, model);
                  model->shels = g_slist_prepend(model->shels, shell);

                  shell->x[0] = str_to_float(*(buff+2));
                  shell->x[1] = str_to_float(*(buff+3));
                  shell->x[2] = str_to_float(*(buff+4));
                  }
                else
                  break;

                shell->breathe = breathe;
                shell->region = region;
                shell->translate = translate;

                if (breathe)
                  {
                  if (num_tokens > 5)
                    shell->radius = str_to_float(*(buff+5));
/* FIXME - what about charge? */
                  }
                else
                  {
                  i = 5;
                  remaining = num_tokens - i;
                  remaining = CLAMP(remaining, 0, 5);

                  if (noflags)
                    remaining = CLAMP(remaining, 0, 2); 

                  switch (remaining)
                    {
                    case 5:
                      shell->lookup_charge = FALSE;
                      shell->charge = str_to_float(*(buff+i));
                      shell->has_sof = TRUE;
                      shell->sof = str_to_float(*(buff+i+1));
                      shell->flags = g_strdup_printf("%s %s %s", *(buff+i+2), *(buff+i+3), *(buff+i+4));
                      break;

                    case 4:
                      shell->lookup_charge = FALSE;
                      shell->charge = str_to_float(*(buff+i));
                      i++;
                    case 3:
                      shell->flags = g_strdup_printf("%s %s %s", *(buff+i), *(buff+i+1), *(buff+i+2));
                      break;

                    case 2:
                      shell->has_sof = TRUE;
                      shell->sof = str_to_float(*(buff+i+1));
                    case 1:
                      shell->lookup_charge = FALSE;
                      shell->charge = str_to_float(*(buff+i));
                      break;
                    }
                  }
                if (gen_flag)
                  {
                  VEC3MUL(shell->x, AU2ANG);
                  }
                ns++;
                break;
              }
#endif

            g_strfreev(buff);
            buff = scan_get_tokens(scan, &num_tokens);
            }
          break;

        case D_HKL:
          if (num_tokens > 1)
            model->surface.dspacing = str_to_float(*(buff+1));
          context = -1;
          break;

        case SBULK_ENERGY:
          if (num_tokens > 1)
            model->gulp.sbulkenergy = str_to_float(*(buff+1));
          context = -1;
          break;

        case TOTAL_ENERGY:
          if (num_tokens > 1)
            model->gulp.energy = str_to_float(*(buff+1));
          context = -1;
          break;

        case SPECIES:
          g_strfreev(buff);
          buff = scan_get_tokens(scan, &num_tokens);
          while (!scan_complete(scan))
            {
/* end on an option */
            if (num_tokens)
              {
              if (gulp_is_option(*buff) > -1)
                {
                scan_put_line(scan);
                break;
                }

/* prevent core dump on comment */
              if (**buff == '#')
                break;

/* add new species data if the 1st item an element, otherwise end */
              if (elem_test(*buff))
                g_string_sprintfa(species_buff, "%s", scan_cur_line(scan));
              else
                break;
              }
            g_strfreev(buff);
            buff = scan_get_tokens(scan, &num_tokens);
            }
          context = -1;
          break;

        case ELEMENT:
          g_strfreev(buff);
          buff = scan_get_tokens(scan, &num_tokens);
          while (!scan_complete(scan))
            {
            if (num_tokens)
              {
/* end on an option */
              temp_code = gulp_is_option(*buff);
              if (temp_code != COVALENT && temp_code != IONIC && temp_code != VDW)
                {
                scan_put_line(scan);
                break;
                }
              else
                g_string_sprintfa(elem_buff, "%s", scan_cur_line(scan));
              }
            g_strfreev(buff);
            buff = scan_get_tokens(scan, &num_tokens);
            }
          context = -1;
          break;

        case KPOINTS:
/* get all subsequent non-empty & numeric lines */
          while (!scan_complete(scan))
            {
            g_strfreev(buff);
            buff = scan_get_tokens(scan, &num_tokens);
            if (num_tokens)
              {
              if (str_is_float(*buff))
                g_string_sprintfa(kpoints_buff, "%s", scan_cur_line(scan));
              }
            else
              break;
            }
          context = -1;
          break;

        case DUMP_FILE:
          if (num_tokens > 1)
            {
/* NEW - grab everything except the dump option itself */
            g_free(model->gulp.dump_file);
            line = scan_cur_line(scan);
            model->gulp.dump_file = g_strdup(get_token_pos(line,1));
            strip_extra(model->gulp.dump_file);
            }
          context = -1;
          break;

        case IGNORE:
/* read until hit an ERONGI */
          while (!scan_complete(scan))
            {
            line = scan_get_line(scan);
            if (gulp_is_option(line) == ERONGI)
              break;
            }
          context = -1;
          break;

/* special case - only one line expected */
        case GULP_SINGLE_LINE_POTENTIAL:
          g_string_sprintfa(ff_buff, "%s", scan_cur_line(scan));
          context = -1;
          break;

/* multi-line potentials */
	case GULP_POTENTIAL:
/* store potential header */
          g_string_sprintfa(ff_buff, "%s", scan_cur_line(scan));

          while (!scan_complete(scan))
            {
            line = scan_get_line(scan);

/* ignore (but keep reading) empty/blank lines */
/* NB - keep both tests below (NULL and \n tests) */
            if (!line)
              continue;
            if (strlen(line) < 2)
              continue;
/* terminate on option */
            if (gulp_is_option(line) > -1)
              {
              context = -1;
              scan_put_line(scan);
              break;
              }
/* terminate if some chars of 1st token are alphabetic and not an element */
            g_strfreev(buff);
            buff = tokenize(line, &num_tokens);
            if (num_tokens)
              {
              tmp = NULL;
              for (i=(gint) strlen(*buff) ; i-- ; )
                {
                if (g_ascii_isalpha((*buff)[i]))
                  {
                  tmp = *buff;
                  break;
                  }
                }
              if (tmp)
                {
/* NEW - allow GULP wildcard element */
                if (g_ascii_strncasecmp(tmp, "x\0", 2) != 0)
                  {
                  if (!elem_symbol_test(tmp))
                    {
                    context = -1;
                    scan_put_line(scan);
                    break;
                    }
                  }
                } 
              }
/* store potential data */
            g_string_sprintfa(ff_buff, "%s", line);
            }
          break;
        }

/* get next non-empty line */
/* TODO - encapsulate in function call? (buff & line) */
  g_strfreev(buff);
  buff = scan_get_tokens(scan, &num_tokens);
  }
g_strfreev(buff);

/* TODO - if keywords found > 0, but doesn't match num tokens */
/* then we've found an unrecognized keyword -> pass through?/print warning? */
/*
printf("Keywords processed: %d/%d\n", keywords_found, keywords_expected);
*/

/* NEW */
if (gulp.solvation_model != GULP_SOLVATION_NONE)
  if (!gulp_shape_indices_compute(&gulp))
    printf("WARNING: input points may be invalid.\n");

/* setup newly loaded models */
for (item=sysenv.mal ; item ; item=g_slist_next(item))
  {
  model = item->data;

  if (model->id == -1)
    {
    model->id = GULP;
    strcpy(model->filename, filename);

/* NEW - unrecognized keyword transfer */
    if (keywords_found != keywords_expected)
      {
      g_free(model->gulp.extra_keywords);
      model->gulp.extra_keywords = g_strdup(key_buff->str);
      }

/* no depth info for loaded surfaces, so pretend they're MARVIN style */
    if (model->periodic == 2)
      model->surface.true_cell = FALSE;

/* fill out current model with GULP data */

/* NEW - connectivity */
model->gulp.noautobond = noautobond;
if (noautobond)
  model->build_molecules = FALSE;

/* solvation */
model->gulp.solvation_model = gulp.solvation_model;
model->gulp.cosmo_shape = gulp.cosmo_shape;
model->gulp.cosmo_shape_index[0] = gulp.cosmo_shape_index[0];
model->gulp.cosmo_shape_index[1] = gulp.cosmo_shape_index[1];
model->gulp.cosmo_points = gulp.cosmo_points;
model->gulp.cosmo_segments = gulp.cosmo_segments;
model->gulp.cosmo_smoothing = gulp.cosmo_smoothing;
model->gulp.cosmo_solvent_epsilon = gulp.cosmo_solvent_epsilon;
model->gulp.cosmo_solvent_radius = gulp.cosmo_solvent_radius;
model->gulp.cosmo_solvent_delta = gulp.cosmo_solvent_delta;
model->gulp.cosmo_solvent_rmax = gulp.cosmo_solvent_rmax;


/* transfer universal gulp data */
/* FIXME - deprec - use gulp_pak instead (see above cosmo stuff) */
    model->gulp.run = run;
    model->gulp.method = method;
    model->gulp.optimiser = optimiser;
    if (optimiser2 > 0)
      {
      model->gulp.optimiser2 = optimiser2;
      model->gulp.switch_type = switch_type;
      model->gulp.switch_value = switch_value;
      }
    model->gulp.unit_hessian = unit;
    model->gulp.coulomb = coulomb;
    model->gulp.maxcyc = maxcyc;
    model->gulp.qeq = qeq;
    model->gulp.free = free;
    model->gulp.zsisa = zsisa;
    model->gulp.compare = compare;
    model->gulp.nosym = nosym;
    model->gulp.fix = fix;
    model->gulp.phonon = phonon;
    model->gulp.eigen = eigen;

    model->gulp.kpoints = g_strdup(kpoints_buff->str);
    model->gulp.potentials = g_strdup(ff_buff->str);
/* NEW */
model->ff_list =  ff_gulp_parse(ff_buff->str);

/* duplicate species data */
    model->gulp.species = g_strdup(species_buff->str);

    buff = tokenize(model->gulp.species, &num_tokens);
    i = 0;
    while (i < num_tokens)
      {
      type = elem_test(*(buff+i));
      j = i;

      get_elem_data(type, &elem_data, model);
      if (type && ++i < num_tokens)
        {
        if (g_ascii_strncasecmp(*(buff+i), "shel", 4) == 0)
          {
          if (++i < num_tokens)
            {
            elem_data.shell_charge = str_to_float(*(buff+i));
            for (slist=model->shels ; slist ; slist=g_slist_next(slist))
              {
              shell = slist->data;
              if (shel_match(*(buff+j), shell))
                {
/* NEW - only assign species data if explict charges were absent */
                if (shell->lookup_charge)
                  {
                  shell->charge = elem_data.shell_charge;
                  shell->lookup_charge = FALSE;
                  }
                }
              }
            }
          }
        else
          {
          if (g_ascii_strncasecmp(*(buff+i), "core", 4) == 0)
            {
            if (++i < num_tokens)
              elem_data.charge = str_to_float(*(buff+i));
            else
              elem_data.charge = 0.0;
            }
          else
            elem_data.charge = str_to_float(*(buff+i));
          for (clist=model->cores ; clist ; clist=g_slist_next(clist))
            {
            core = clist->data;
            if (core_match(*(buff+j), core))
              {
/* NEW - only assign species data if explict charges were absent */
              if (core->lookup_charge)
                {
                core->charge = elem_data.charge;
                core->lookup_charge = FALSE;
                }
              }
            }
          }
        } 
      i++;
      }
    g_strfreev(buff);

/* process element data */
    model->gulp.elements = g_strdup(elem_buff->str);
    buff = tokenize(model->gulp.elements, &num_tokens);
    i = 0;
    type = -1;
    while (i < num_tokens)
      {
      switch (type)
        {
        case COVALENT:
          code = elem_test(*(buff+i));
          i++;
          if (!get_elem_data(code, &elem_data, model) && i<num_tokens)
            {
            elem_data.cova = str_to_float(*(buff+i));
            put_elem_data(&elem_data, model);
            }
          type = -1;
          break;

        case VDW:
          code = elem_test(*(buff+i));
          i++;
          if (!get_elem_data(code, &elem_data, model) && i<num_tokens)
            {
            elem_data.vdw = str_to_float(*(buff+i));
            put_elem_data(&elem_data, model);
            }
          type = -1;
          break;

        default:
          if (g_ascii_strncasecmp(*(buff+i), "cova", 4) == 0)
            type = COVALENT;
          if (g_ascii_strncasecmp(*(buff+i), "vdw", 3) == 0)
            type = VDW;
        }
      i++;
      }
    g_strfreev(buff);

/* extra (unprocessed lines) */
    if (strlen(ex_buff->str))
      model->gulp.extra = g_strdup(ex_buff->str);

/* NEW - store initial frame coord type, as trajectory animation will overwrite */
    model->gulp.orig_fractional = model->fractional;
    model->gulp.orig_construct_pbc = model->construct_pbc;

/* ensure the ordering is correct, as GULP trajectory files depend on it */
    model->cores = g_slist_reverse(model->cores);
    model->shels = g_slist_reverse(model->shels);

/* display init */
    model_prep(model);
    }
  }

/* FIXME - duplicate ff list for all models */
/*
if (model)
  {
printf("FIXME - adding forcefields [%d] to model %p\n", g_slist_length(ff_list), model);
ff_dump_all(ff_list);

  model->ff_list = ff_list;
  }
*/

/* cleanup */
g_string_free(key_buff, TRUE);
g_string_free(elem_buff, TRUE);
g_string_free(ex_buff, TRUE);
g_string_free(ff_buff, TRUE);
g_string_free(species_buff, TRUE);
g_string_free(kpoints_buff, TRUE);

scan_free(scan);
return(0);
}

/*******************************************************************/
/* primitive for adding two strings to generate the property value */
/*******************************************************************/
void property_add_with_units(gint rank,
                             gchar *key, gchar *value, gchar *units,
                             struct model_pak *model)
{
gchar *text;

text = g_strjoin(" ", value, units, NULL);
property_add_ranked(rank, key, text, model);
g_free(text);
}

/*****************************/
/* generic gulp output parse */
/*****************************/
/* NB: new model is not made, as only energy etc. values */
/* are of interest if we're parsing the results of a gulp run */
#define DEBUG_READ_GULP_OUTPUT 0
gint read_gulp_output(gchar *filename, struct model_pak *data)
{
gint i, j, m, n, flag, region, primitive=FALSE, shift=0 /*, comp=0,*/;
gint /*status,*/ count=0, dflag=0, fflag=0, pflag=0, read_pass=0, num_tokens;/*FIX 68c4ec e1cc11*/
gint num_structures = 0;
gchar line[LINELEN], **buff, *text, coord;
gdouble *value;
GSList *list, *flist=NULL, *ir_list=NULL, *raman_list=NULL, *model_list=NULL;
GString *build_text;
GHashTable *types=NULL;
struct model_pak *temp;
struct core_pak *core;
struct shel_pak *shell;
FILE *fp;

/* checks */
g_assert(data != NULL);
g_assert(filename != NULL);

fp = fopen(filename, "rt");
if (!fp)
  return(1);

/* init */
data->region_empty[REGION1A] = FALSE;
data->gulp.esurf[0] = 0.0;
data->gulp.esurf[1] = 0.0;
data->gulp.eatt[0] = 0.0;
data->gulp.eatt[1] = 0.0;

types = g_hash_table_new_full(g_str_hash, g_str_equal, g_free, g_free);

#if DEBUG_READ_GULP_OUTPUT
printf("Grafted: %d\n", data->grafted);
#endif

/* init for new model (ie not on the tree) */
if (!data->grafted)
  {
  data->id = GULPOUT;
  g_free(data->basename);
  data->basename = parse_strip(filename);
  strcpy(data->filename, filename);
  }

/* append model to list irrespective of whether its grafted or not */
  model_list = g_slist_append(model_list, data);

/* look for interesting output */
flag=0;//status=0;/*FIX 68c4ec e1cc11*/
while(!fgetline(fp,line))
  {
  buff = tokenize(line, &num_tokens);

/* GULP trapping */
  if (g_ascii_strncasecmp("!! ERROR",line,8) == 0)
    {
    data->error_file_read = g_string_append(data->error_file_read, line);
    g_strfreev(buff);
//  status = 2;/*FIX e1cc11*/
    break;
    }
  if (g_ascii_strncasecmp("  **** Unit cell is not charge neutral",line,38) == 0)
    {
    data->error_file_read = g_string_append(data->error_file_read, line);
    if (!fgetline(fp,line))
      data->error_file_read = g_string_append(data->error_file_read, line);
    g_strfreev(buff);
//  status = 2;/*FIX 68c4ec*/
    break;
    }
/*
  if (g_ascii_strncasecmp("  **** Warning", line, 14) == 0)
    gui_text_show(WARNING, line);
*/
  
  /* make new model if new structure number */
  if (g_ascii_strncasecmp("*  Input for Configuration", line, 26) == 0)
    {
    if (num_tokens > 5)
      {
      if ((gint) str_to_float(*(buff+5)) == (num_structures + 1))
        num_structures++;
      else
        {
        build_text = g_string_new("");
        g_string_printf(build_text, "Inconsistent configuration numbers in GULP: Expecing %d, got %d\n", num_structures+1, (gint)str_to_float(*(buff+5)));
        gui_text_show(ERROR, build_text->str);
        g_string_free(build_text, TRUE);
        return(1);
        }
      if (num_structures > 1)
        {
        /* create new model */
        temp = model_new();
        temp->id = GULPOUT;
        /* this is now our new model */
        data = temp;
        model_list = g_slist_append(model_list, data);
        }
      }
    /* have name? */
    build_text = g_string_new("");
    if (num_tokens > 7)
      {
      g_free(data->basename);
      data->basename = g_strdup(buff[7]);
      g_string_free(build_text, TRUE);
      }
    else
      {
      g_string_printf(build_text, "%s_%d", data->basename, (gint)str_to_float(*(buff+5)));
      g_free(data->basename);
      data->basename = build_text->str;
      }
    }
  
/* Swap structure if neccessary to match output with input */
  if (g_ascii_strncasecmp("*  Output for Configuration", line, 17) == 0)
    {
    /* put energy property into model before switching */

    text = g_strdup_printf("%.5f eV", data->gulp.energy);
    property_add_ranked(2, "Energy", text, data);
    g_free(text);

#if DEBUG_READ_GULP_OUTPUT
    printf("retrieving structure %d\n", (gint) str_to_float(*(buff+4)) - 1);
#endif
    data = g_slist_nth_data (model_list, (gint) str_to_float(*(buff+4)) - 1);
    }
  
/* periodicity */
  if (g_ascii_strncasecmp("  Dimensionality", line, 16) == 0)
    {
    if (num_tokens > 2)
      data->periodic = str_to_float(*(buff+2));
    }

/* space group */
  if (g_ascii_strncasecmp("  Space group ", line, 14) == 0)
    {
    if (!data->grafted)
      {
      for (i=(gint) strlen(line) ; i-- ; )
        {
        if (line[i] == ':')
          {
          data->sginfo.spacenum = -1;
          data->sginfo.spacename = g_strstrip(g_strdup(&line[i+1]));

/*
          if (data->sginfo.spacename[0] == 'R')
            primitive = TRUE;
      if (g_strrstr(data->sginfo.spacename, ":R") )
        {
printf("1 prim\n");
        primitive = TRUE;
        }
*/

#if DEBUG_READ_GULP_OUTPUT
printf("space group: %s\n", data->sginfo.spacename);
#endif
          }
        }
      }
    else
      {
/* NEW - check if we should be looking for the primitve cell */
/*
      if (g_strrstr(data->sginfo.spacename, ":R") )
        {
printf("2 prim\n");
        primitive = TRUE;
        }
*/

      }
    }

/* NEW - process GULP's library symbols as FF types */
  if (g_ascii_strncasecmp("  Species output", line, 16) == 0)
    {
    for (i=4 ; i-- ; )
      {
      g_strfreev(buff);
      buff = get_tokenized_line(fp, &num_tokens);
      }

    while (*buff)
      {
      g_strfreev(buff);
      buff = get_tokenized_line(fp, &num_tokens);

      if (elem_symbol_test(*buff))
        {
        if (num_tokens > 8)
          {
          g_hash_table_insert(types, g_strdup(*buff), g_strdup(*(buff+8)));
          }
        }
      else
        break; 
      }
    }

/* create new model if the calculation was a defect calc */
  if (g_ascii_strncasecmp("*  Defect calculation for configuration", line, 39) == 0)
    {
/* create new model */
    temp = model_new();
    temp->id = GULPOUT;
    g_free(temp->basename);
    temp->basename = g_strdup_printf("%s_defect_r1", data->basename);
/* this is now our new model */
    data = temp;
    model_list = g_slist_append(model_list, data);
    }

/* create new model on multiple instances of output data */
/* introducted to cope with output from the translate option */
//  if (g_ascii_strncasecmp("  Components of energy", line, 22) == 0)
//    {
//    if (comp && !(comp % 2))
//      {
///* create new model */
//      temp = model_new();
//      temp->id = GULPOUT;
///* pass model data through */
//      temp->fractional = data->fractional;
//      temp->periodic = data->periodic;
//      temp->sginfo.spacename = g_strdup(data->sginfo.spacename);
//      temp->sginfo.spacenum = -1;
//      temp->construct_pbc = data->construct_pbc;
//      if (data->construct_pbc)
//        memcpy(temp->latmat, data->latmat, 9*sizeof(gdouble));
//      else
//        memcpy(temp->pbc, data->pbc, 6*sizeof(gdouble));
///* this is now our new model */
//      data = temp;
//      model_list = g_slist_append(model_list, data);
//      }
//    comp++;
//    }

/* minimized energy search */
  if (g_ascii_strncasecmp("  Final energy =",line,16) == 0)
    {
    if (num_tokens > 3)
      {
/* record total energy */
      data->gulp.energy = str_to_float(*(buff+3));
      flag++;
      }
    }

/* surface energy search */
  if (g_ascii_strncasecmp("  Surface energy ",line,17) == 0)
    {
    if (num_tokens > 6)
      {
      if (data->gulp.esurf[0] == 0.0)
        {
        data->gulp.esurf[0] = str_to_float(*(buff+5));
        property_add_with_units(3, "Esurf (u)", *(buff+5), *(buff+6), data);
        }
      else
        {
        data->gulp.esurf[1] = str_to_float(*(buff+5));
        property_add_with_units(4, "Esurf (r)", *(buff+5), *(buff+6), data);
        }
      g_free(data->gulp.esurf_units);
      data->gulp.esurf_units = g_strdup(*(buff+6));
      }
    }

/* attachment energy search 1 */
  if (g_ascii_strncasecmp("  Attachment energy ",line,20) == 0)
    {
    if (num_tokens > 4)
      {
      if (data->gulp.eatt[0] == 0.0)
        {
        data->gulp.eatt[0] = str_to_float(*(buff+3));
        property_add_with_units(5, "Eatt (u)", *(buff+3), *(buff+4), data);
        }
      else
        {
        data->gulp.eatt[1] = str_to_float(*(buff+3));
        property_add_with_units(6, "Eatt (r)", *(buff+3), *(buff+4), data);
        }
      g_free(data->gulp.eatt_units);
      data->gulp.eatt_units = g_strdup(*(buff+4));
      }
    }

/* attachment energy search 2 */
  if (g_ascii_strncasecmp("  Attachment energy/",line,20) == 0)
    {
    if (num_tokens > 3)
      {
/* NB: assumes this match occurs just after the previous match */
      if (data->gulp.eatt[1] == 0.0)
        {
        data->gulp.eatt[0] = str_to_float(*(buff+3));
        property_add_with_units(5, "Eatt (u)", *(buff+3), "eV/unit", data);
        }
      else
        {
        data->gulp.eatt[1] = str_to_float(*(buff+3));
        property_add_with_units(6, "Eatt (r)", *(buff+3), "eV/unit", data);
        }
      g_free(data->gulp.eatt_units);
      data->gulp.eatt_units = g_strdup("eV/unit");
      }
    }

/* search for gnorm */
  if (g_ascii_strncasecmp("  Final Gnorm  =",line,16) == 0)
    {
    if (num_tokens > 3)
      {
      data->gulp.gnorm = str_to_float(*(buff+3));
      property_add_ranked(9, "Gnorm", *(buff+3), data);
      flag++;
      }
    }

/* single point energy search */
  if (g_ascii_strncasecmp("  Total lattice energy",line,22) == 0)
    {
/* NB: assumes we have the correct cell parameters */
    primitive = space_primitive_cell(data);

/* if nothing on the rest of this line, then should be on the next */
    switch (num_tokens)
      {
      case 4:
      while(!fgetline(fp,line))
        {
        if (primitive)
          {
          if (g_ascii_strncasecmp("    Primitive unit cell",line,23) == 0)
            {
            g_strfreev(buff);
            buff = tokenize(line, &num_tokens);
            if (num_tokens > 5)
              data->gulp.energy = str_to_float(*(buff+4));
            break;
            }
          }
        else
          {
          if (g_ascii_strncasecmp("    Non-primitive unit cell",line,27) == 0)
            {
            g_strfreev(buff);
            buff = tokenize(line, &num_tokens);
            if (num_tokens > 5)
              data->gulp.energy = str_to_float(*(buff+4));
            break;
            }
          }
        }
      break;

      case 6:
        if (g_ascii_strncasecmp("eV",*(buff+5),2) == 0)
          data->gulp.energy = str_to_float(*(buff+4));
      break;

      }
    }

/* single point free energy search */
  if (g_ascii_strncasecmp("  Initial surface dipole ",line,25) == 0)
    {
    if (num_tokens > 4)
      data->gulp.sdipole = str_to_float(*(buff+4));
    }

/* single point free energy search */
  if (g_ascii_strncasecmp("  Total free energy",line,19) == 0)
    {
    if (num_tokens > 5)
      if (g_ascii_strncasecmp("eV",*(buff+5),2) == 0)
        data->gulp.energy = str_to_float(*(buff+4));
    }

/* NEW - vibrational info search */
  if (g_ascii_strncasecmp(" Frequencies (cm-1) and Eigenvectors",line,36) == 0)
    fflag++;
  if (g_ascii_strncasecmp("  Zero point energy",line,19) == 0)
    fflag=0;

  if (fflag && g_ascii_strncasecmp(" Frequency",line,10) == 0)
    {
/* space group processing here, since we want to attach */
/* the vibration lists to ALL atoms in the full cell, */
/* not just the asymmetric ones. */
    if (!data->grafted && !pflag)
      {
      model_prep(data);
      pflag++;
      }

    for (i=1 ; i<num_tokens ; i++)
      flist = g_slist_prepend(flist, g_strdup(*(buff+i)));

/* skip to data */
    while (!fgetline(fp,line))
      {
      g_strfreev(buff);
      buff = tokenize(line, &n);

/* NEW - store IR intensities */
      if (g_ascii_strncasecmp(" IR Intensity", line, 13) == 0)
        {
//      i = 2;/*FIX 2b7a7b*/
        for (i=2 ; i<n ; i++)
          ir_list = g_slist_prepend(ir_list, g_strdup(*(buff+i)));
        }

/* NEW - store Raman intensities */
      if (g_ascii_strncasecmp(" Raman", line, 6) == 0)
        {
//      i = 2;/*FIX b35039*/
        for (i=2 ; i<n ; i++)
          raman_list = g_slist_prepend(raman_list, g_strdup(*(buff+i)));
        }

      if (n > 1)
        coord = **(buff+1);
      else
        coord = '?';

      if (coord == 'x')
        break;
      }

/* read the eigenvectors in */
    m=0;
    while (m<3*g_slist_length(data->cores))
      {
      g_strfreev(buff);
      buff = tokenize(line, &n);
/* blank line termination */
      if (!n)
        break;
/* get reference atom number (NB: GULP starts at 1, list starts at 0) */
      i = str_to_float(*buff);
      i--;

/* FIXME - if we're reading in cores from this gulp file - order will be reverse */
/* due to prepending to core list */
      core = g_slist_nth_data(data->cores, i);
      g_assert(core != NULL);

/*
      if (i<0 || i>data->num_atoms-1)
        {
        printf("Bad line: %s\n", line);
        printf("ref atom: %d, num atoms: %d\n", i, data->num_atoms);
        g_assert_not_reached();
        }
*/
      for (j=2 ; j<n ; j++)
        {
/* get the value */
        value = g_malloc(sizeof(gdouble));
        *value = str_to_float(*(buff+j));
/* get reference coordinate */
        coord = **(buff+1);
        switch(coord)
          {
          case 'x':
            core->vibx_list = g_slist_append(core->vibx_list,
                                         (gpointer *) value);
            break;
          case 'y':
            core->viby_list = g_slist_append(core->viby_list,
                                         (gpointer *) value);
            break;
          case 'z':
            core->vibz_list = g_slist_append(core->vibz_list,
                                         (gpointer *) value);
            break;
          default:
            g_assert_not_reached();
          }
        }

/* next line (if available) */
      if (fgetline(fp,line))
        break;
      m++;
      }
    }

/* minimized coords search */
/* TODO - unminimized coords? */
/* added 'of surface' to avoid reading in growth slice */
  dflag=0;
  if (g_ascii_strncasecmp("  Final asymmetric unit coordinates",line,35) == 0
   || g_ascii_strncasecmp("  Final fractional coordinates",line,30) == 0
   || g_ascii_strncasecmp("  Fractional coordinates of asymmetric",line,38) == 0
   || g_ascii_strncasecmp("  Final fractional/Cartesian coordinates", line, 40) == 0
   || g_ascii_strncasecmp("  Mixed fractional/Cartesian coordinates", line, 40) == 0)
    {
    data->fractional = TRUE;
    dflag++;
    }
  else if (g_ascii_strncasecmp("  Final coordinates of region", line, 29) == 0)
    {
    dflag++;
    shift=1;
    }


/* don't read in coords if already on the tree (eg single pt) */
  if (dflag && !data->grafted)
    {
/* enforce empty core/shell lists */
    free_core_list(data);
/* skip the two (minimum number) ----- dividers */
    region = REGION1A;
    count=0;
    while (count < 2)
      {
      if (fgetline(fp,line))
        return(4);
      if (g_ascii_strncasecmp("-------", line, 7) == 0)
        count++;
      }

/* loop until we run out of atomic data lines */
    i=j=n=0;
    flag=0;
/* read as many atoms as we can */
    while(!fgetline(fp,line))
      {
      g_strfreev(buff);
      buff = tokenize(line, &num_tokens);

/* TODO - dont exit if ----- (could be region divider) */
/*
      if (g_ascii_strncasecmp("-------", line, 7) == 0)
        continue;
*/

/* blank line termination */
if (!num_tokens)
  break;

/* determine region */
      if (g_ascii_strncasecmp("  Region ", line, 9) == 0)
        {
        if (num_tokens > 1)
        switch(*(*(buff+1)))
          {
          case '1':
#if DEBUG_READ_GULP_OUTPUT
printf("Labelling region 1...\n");
#endif
            region = REGION1A;
            break;
          case '2':
#if DEBUG_READ_GULP_OUTPUT
printf("Labelling region 2...\n");
#endif
            region = REGION2A;
            break;
          default:
            printf("WARNING: unknown region specification\n");
            region = REGION1A;
          }
        data->region_empty[region] = FALSE;
        continue;
        }

/* HACK - GULP sometimes puts *'s after coords (region 1 only?) */
      if (num_tokens > 7)
      if (g_ascii_strncasecmp("*", *(buff+4), 1) == 0)
        {
        g_free(*(buff+4)); 
        *(buff+4) = g_strdup(*(buff+5));
        g_free(*(buff+5)); 
        *(buff+5) = g_strdup(*(buff+7));
        }

/* exit only if insufficient data items (NB: changes b/w sections, min=6) */
      if (num_tokens > 6)
        {
/* read the data in - assume the same order */
        if (g_ascii_strncasecmp("c", *(buff+2), 1) == 0)
          {
          core = new_core(*(buff+1), data);
          data->cores = g_slist_prepend(data->cores, core);

          if (!read_pass)
            core->region = region;
          core->x[0] = str_to_float(*(buff+3));
          core->x[1] = str_to_float(*(buff+4));
          core->x[2] = str_to_float(*(buff+5));
          core->charge = str_to_float(*(buff+6+shift));
          core->lookup_charge = FALSE;
/*
          core->sof = str_to_float(*(buff+7));
*/
          i++;

#if DEBUG_READ_GULP_OUTPUT
printf("core coords: %lf %lf %lf\n", core->x[0], core->x[1], core->x[2]);
#endif
          }

        if (g_ascii_strncasecmp("s", *(buff+2), 1) == 0)
          {
          shell = new_shell(*(buff+1), data);
          data->shels = g_slist_prepend(data->shels, shell);

          if (!read_pass)
            shell->region = region;
          shell->x[0] = str_to_float(*(buff+3));
          shell->x[1] = str_to_float(*(buff+4));
          shell->x[2] = str_to_float(*(buff+5));
          shell->charge = str_to_float(*(buff+6+shift));
          shell->lookup_charge = FALSE;

#if DEBUG_READ_GULP_OUTPUT
printf("shel coords: %lf %lf %lf\n", shell->x[0], shell->x[1], shell->x[2]);
#endif

          j++;
          }
        }

      }

data->cores = g_slist_reverse(data->cores);
data->shels = g_slist_reverse(data->shels);

/* done one pass => read in unrelaxed coords if reading optimization */
    read_pass++;

#if DEBUG_READ_GULP_OUTPUT
printf("Retrieved %d atoms & %d shells (periodic).\n",i,j);
#endif
    }

/* cell parameters - search 1 - read from Final cell parameters block */
if (g_ascii_strncasecmp("  Final cell parameters and derivatives", line, 39) == 0 && !data->grafted)
  {
  
  data->periodic = 3;
  data->construct_pbc = FALSE;
  
  /* skip to 1st line of data */
  fgetline(fp,line);
  fgetline(fp,line);
  
  /* read cell lengths */
  for (i=0; i<3; i++)
    {
    g_strfreev(buff);
    buff = get_tokenized_line(fp, &num_tokens);
    if (num_tokens > 1)
      data->pbc[i] = str_to_float(*(buff+1));
    }

  /* read cell angles */

  for (i=3; i<6; i++)
    {
    g_strfreev(buff);
    buff = get_tokenized_line(fp, &num_tokens);
    if (num_tokens > 1)
      data->pbc[i] = D2R*str_to_float(*(buff+1));
    }
  }
  
///* cell parameters - search 1 */
//  if (g_ascii_strncasecmp("  Final cell parameters",line,22) == 0 && !data->grafted)
//    {
///* skip to data */
//    fgetline(fp,line);
//    fgetline(fp,line);
///* get cell lengths */
//    for (i=0 ; i<6 ; i++)
//      {
//      if (fgetline(fp,line))
//        break;
//      g_strfreev(buff);
//      buff = tokenize(line, &num_tokens);
//      if (num_tokens > 1)
//        data->pbc[i] = str_to_float(*(buff+1));
//      }
///* convert to radians */
//    data->pbc[3] *= D2R;
//    data->pbc[4] *= D2R;
//    data->pbc[5] *= D2R;
//    }
//

/* cell parameters - search 2 */
  if (g_ascii_strncasecmp("  Non-primitive lattice parameters",line,34) == 0 && !data->grafted)
    {
/* skip to data */
    if (fgetline(fp,line))
      break;
    if (fgetline(fp,line))
      break;

/* get cell lengths */
    g_strfreev(buff);
    buff = tokenize(line, &num_tokens);
    if (num_tokens > 8)
      {
      data->pbc[0] = str_to_float(*(buff+2));
      data->pbc[1] = str_to_float(*(buff+5));
      data->pbc[2] = str_to_float(*(buff+8));
      }   
/* get cell angles */
    if (fgetline(fp,line))
      break;
    g_strfreev(buff);
    buff = tokenize(line, &num_tokens);
    if (num_tokens > 5)
      {
      data->pbc[3] = D2R*str_to_float(*(buff+1));
      data->pbc[4] = D2R*str_to_float(*(buff+3));
      data->pbc[5] = D2R*str_to_float(*(buff+5));
      }
    }

/* cell parameters - search 3 */
  if (g_ascii_strncasecmp("  Final surface cell parameters ", line, 32) == 0 && !data->grafted)
    {
    data->periodic = 2;
    data->construct_pbc = FALSE;
    data->fractional = TRUE;
    data->pbc[2] = 0.0;
    data->pbc[3] = PI/2.0;
    data->pbc[4] = PI/2.0;

/* skip to 1st line of data */
    fgetline(fp,line);
    fgetline(fp,line);

    g_strfreev(buff);
    buff = get_tokenized_line(fp, &num_tokens);
    if (num_tokens > 1)
      data->pbc[0] = str_to_float(*(buff+1));

    g_strfreev(buff);
    buff = get_tokenized_line(fp, &num_tokens);
    if (num_tokens > 1)
      data->pbc[1] = str_to_float(*(buff+1));

    g_strfreev(buff);
    buff = get_tokenized_line(fp, &num_tokens);
    if (num_tokens > 1)
      data->pbc[5] = D2R*str_to_float(*(buff+1));
    }

/* cell parameters - search 4 (eg a CONV surface calc) */
  if (g_ascii_strncasecmp("  Surface cell parameters ",line, 26) == 0 && !data->grafted)
    {
    data->periodic = 2;
    data->construct_pbc = FALSE;
    data->fractional = TRUE;
    data->pbc[2] = 0.0;
    data->pbc[3] = PI/2.0;
    data->pbc[4] = PI/2.0;
/* skip to 1st line of data */
    fgetline(fp,line);

    g_strfreev(buff);
    buff = get_tokenized_line(fp, &num_tokens);
    if (num_tokens > 5)
      {
      data->pbc[0] = str_to_float(*(buff+2));
      data->pbc[5] = D2R*str_to_float(*(buff+5));
      }
    g_strfreev(buff);
    buff = get_tokenized_line(fp, &num_tokens);
    if (num_tokens > 1)
      data->pbc[1] = str_to_float(*(buff+2));
#if DEBUG_READ_GULP_OUTPUT
    printf("Search 4\n");
    P3VEC("cell lengths ",&data->pbc[0]);
    P3VEC("cell angles ",&data->pbc[3]);
#endif
    }
  

/* cell parameters - search 5 (unrelaxed) */
  if (g_ascii_strncasecmp("  Cell parameters (Angstroms/Degrees):",line,38) == 0 && !data->grafted)
    {
/* skip blank line */
    fgetline(fp,line);
/* get cell lengths */
    for (i=0 ; i<3 ; i++)
      {
      g_strfreev(buff);
      buff = get_tokenized_line(fp, &num_tokens);

      if (num_tokens > 5)
        {
        data->pbc[i] = str_to_float(*(buff+2));
        data->pbc[i+3] = D2R*str_to_float(*(buff+5));
        }
      }
    }

/* cell parameters - search 6 */
  if (g_ascii_strncasecmp("  Surface Cartesian vectors (Angstroms) :",line,41) == 0 && !data->grafted)
    {
    data->periodic = 2;
    data->construct_pbc = TRUE;
/* skip blank line */
    fgetline(fp,line);

/* get vec 1 */
    g_strfreev(buff);
    buff = get_tokenized_line(fp, &num_tokens);
    if (num_tokens > 1)
      {
      data->latmat[0] = str_to_float(*buff);
      data->latmat[3] = str_to_float(*(buff+1));
      data->latmat[6] = 0.0;
      }

/* get vec 2 */
    g_strfreev(buff);
    buff = get_tokenized_line(fp, &num_tokens);
    if (num_tokens > 1)
      {
      data->latmat[1] = str_to_float(*buff);
      data->latmat[4] = str_to_float(*(buff+1));
      data->latmat[7] = 0.0;
      }

/* set vec 3 */
    data->latmat[2] = 0.0;
    data->latmat[5] = 0.0;
    data->latmat[8] = 1.0;
    }

/* cell parameters - search 8 */
  if (g_ascii_strncasecmp("  Primitive cell parameters :",line,29) == 0 && !data->grafted)
    {
    data->periodic = 3;
/* skip blank line */
    fgetline(fp,line);

/* get line 1 */
    g_strfreev(buff);
    buff = get_tokenized_line(fp, &num_tokens);
    if (num_tokens > 11)
      {
      data->pbc[0] = str_to_float(*(buff+8));
      data->pbc[3] = D2R*str_to_float(*(buff+11));
      }
/* get line 2 */
    g_strfreev(buff);
    buff = get_tokenized_line(fp, &num_tokens);
    if (num_tokens > 11)
      {
      data->pbc[1] = str_to_float(*(buff+8));
      data->pbc[4] = D2R*str_to_float(*(buff+11));
      }
/* get line 3 */
    g_strfreev(buff);
    buff = get_tokenized_line(fp, &num_tokens);
    if (num_tokens > 11)
      {
      data->pbc[2] = str_to_float(*(buff+8));
      data->pbc[5] = D2R*str_to_float(*(buff+11));
      }
    }

/* cell parameters - search 9 */
  if (g_ascii_strncasecmp("  Polymer cell parameter",line,24) == 0 && !data->grafted)
    {
/* init */
    data->periodic = 1;
    data->construct_pbc = FALSE;
/* blank line skip */
    fgetline(fp,line);

    g_strfreev(buff);
    buff = get_tokenized_line(fp, &num_tokens);
    if (num_tokens > 2)
      data->pbc[0] = str_to_float(*(buff+2));
    }

/* cell parameters - search 10 */
  if (g_ascii_strncasecmp("  Final polymer cell parameter",line,30) == 0 && !data->grafted)
    {
/* init */
    data->periodic = 1;
    data->construct_pbc = FALSE;
/* blank line skip */
    fgetline(fp,line);
    fgetline(fp,line);

    g_strfreev(buff);
    buff = get_tokenized_line(fp, &num_tokens);
    if (num_tokens > 2)
      data->pbc[0] = str_to_float(*(buff+1));
    }

/* isolated coords search */
  if ((g_ascii_strncasecmp("  Final cartesian coordinates", line, 29) == 0 
   ||  g_ascii_strncasecmp("  Cartesian coordinates of cluster", line, 34) == 0
   ||  g_ascii_strncasecmp("  Region 1 (absolute coordinates)", line, 33) == 0)
   && !data->grafted)
    {
    data->periodic = 0;
    data->fractional = FALSE;

/* enforce empty core/shell lists */
    free_core_list(data);

    for (i=0 ; i<5 ; i++)
      fgetline(fp,line);
/* loop until we run out of atomic data lines */
 #if DEBUG_READ_GULP_OUTPUT
    i=j=0;/*FIX 352de5*/
#endif
    n=0;
    flag=0;
/* don't read any more than num_atoms -> memory allocation problems! */
    while(!fgetline(fp,line) && !flag)
      {
      g_strfreev(buff);
      buff = tokenize(line, &num_tokens);

/* blank line termination */
if (!num_tokens)
  break;

/* exit if atom number is not consecutive */ 
      m = str_to_float(*buff);
      if (m != ++n)
        flag=1;
      else
        {
/* read the data in - assume the same order */
        if (num_tokens > 6)
          {
          if (g_ascii_strncasecmp("c",*(buff+2),1) == 0)
            {
            core = new_core(*(buff+1), data);
            data->cores = g_slist_prepend(data->cores, core);

            core->x[0] = str_to_float(*(buff+3));
            core->x[1] = str_to_float(*(buff+4));
            core->x[2] = str_to_float(*(buff+5));
            core->charge = str_to_float(*(buff+6));
            core->lookup_charge = FALSE;
            }

          if (g_ascii_strncasecmp("s",*(buff+2),1) == 0)
            {
            shell = new_shell(*(buff+1), data);
            data->shels = g_slist_prepend(data->shels, shell);

            shell->x[0] = str_to_float(*(buff+3));
            shell->x[1] = str_to_float(*(buff+4));
            shell->x[2] = str_to_float(*(buff+5));
            shell->charge = str_to_float(*(buff+6));
            shell->lookup_charge = FALSE;
            }
          }
        }
      }

    data->cores = g_slist_reverse(data->cores);
    data->shels = g_slist_reverse(data->shels);

#if DEBUG_READ_GULP_OUTPUT
printf("Retrieved %d atoms & %d shells (cluster).\n",i,j);
#endif
    }

/* NEW */
  if (g_ascii_strncasecmp("  Electrostatic potential at input sites", line, 40) == 0)
    {
/* skip to data */
    for (i=5 ; i-- ; )
      fgetline(fp,line);

/* clean the list */
    if (data->gulp.epot_vals)
      free_slist(data->gulp.epot_vals);
    data->gulp.epot_vals = NULL;

/* acquire lines */
    while (!fgetline(fp, line))
      {
      if (g_ascii_strncasecmp("----------",line,10) == 0)
        break;

      g_strfreev(buff);
      buff = tokenize(line, &num_tokens);
      if (num_tokens > 4)
        {
/* get potential value */
        value = g_malloc(sizeof(gdouble));
        *value = str_to_float(*(buff+4));
        data->gulp.epot_vals = g_slist_prepend(data->gulp.epot_vals, value);
/* record min & max */
        if (*value > data->gulp.epot_max && data->epot_autoscale)
          data->gulp.epot_max = *value;
        if (*value < data->gulp.epot_min && data->epot_autoscale)
          data->gulp.epot_min = *value;
        }
      }

/* NB: prepend -> reverse to get matching direction with epot_vecs */
    data->gulp.epot_vals = g_slist_reverse(data->gulp.epot_vals); 

/* confirm that we got sufficient points */
#if DEBUG_READ_GULP_OUTPUT
    printf("got: %d, expected: %d\n",
           g_slist_length(data->gulp.epot_vecs),
           g_slist_length(data->gulp.epot_vals));
    printf("min: %f, max: %f\n", data->gulp.epot_min, data->gulp.epot_max);
#endif
    }

  g_strfreev(buff);
  }

#if DEBUG_READ_GULP_OUTPUT
P3VEC("cell lengths ",&data->pbc[0]);
P3VEC("cell angles ",&data->pbc[3]);
#endif

text = g_strdup_printf("%.5f eV", data->gulp.energy);
property_add_ranked(2, "Energy", text, data);
g_free(text);

/* surfaces are always const vol */
//if (data->periodic == 2)
//  data->gulp.method = CONV;

/* init lists */
if (flist)
  data->phonons = g_slist_reverse(flist);
else
  data->phonons = NULL;
data->num_phonons = g_slist_length(data->phonons);

if (ir_list)
  data->ir_list = g_slist_reverse(ir_list); 
if (raman_list)
  data->raman_list = g_slist_reverse(raman_list); 

#if DEBUG_READ_GULP_OUTPUT
printf("vibrational eigenvectors: %d\n", data->num_phonons);
if (data->num_phonons)
  dump_phonons(data);
#endif

/* initialize the models for display */
for (list=model_list ; list ; list=g_slist_next(list))
  {
  temp = list->data;

/* NEW - process the library symbols */
  for (flist=temp->cores ; flist ; flist=g_slist_next(flist))
    {
    core = flist->data;

    text = g_hash_table_lookup(types, core->atom_label); 

    if (text)
      {
      g_free(core->atom_type);
      core->atom_type = g_strdup(text);
      }
    }

/* init for display (if not already displayed) */
  if (!temp->grafted && !pflag)
    {
    model_prep(temp);
    }
  }

/* CURRENT */
/*
dump_phonon(289, 873, data);
*/

g_hash_table_destroy(types);

fclose(fp);
return(0);
}

/**************************************/
/* show the current file stream error */
/**************************************/
void print_file_error(FILE *fp)
{

printf("error = %d\n", ferror(fp));

clearerr(fp);
}

/***********************************/
/* GULP dynamics trajectory header */
/***********************************/
/* macro for fortran start/end record padding */
/*
int int_buffer=0;
#define READ_RECORD fread(&int_buffer, sizeof(int), 1, fp)
#define WRITE_RECORD fwrite(&int_buffer, sizeof(int), 1, fp)
*/

#define DEBUG_TRJ_HEADER 0
gint read_trj_header(FILE *fp, struct model_pak *model)
{
int num_atoms, periodic;
double version;

/* record 1 - version */
READ_RECORD;
fread(&version, sizeof(version), 1, fp);
if (fabs(version) > 100.0)
  {
  swap_bytes(&version, sizeof(version));
  if (fabs(version) > 100.0)
    {
    printf("Error: file seems garbled.\n");
    return(1);
    }
  model->trj_swap = TRUE;
  }
READ_RECORD;
/* record 2 */
READ_RECORD;
/* # of atoms + shells */
fread(&num_atoms, sizeof(num_atoms), 1, fp);
if (model->trj_swap)
  swap_bytes(&num_atoms, sizeof(num_atoms));
/* dimension */
fread(&periodic, sizeof(periodic), 1, fp);
if (model->trj_swap)
  swap_bytes(&periodic, sizeof(periodic));
READ_RECORD;

#if DEBUG_TRJ_HEADER
printf("trj version: %lf\n", version);
printf(" Swap bytes: %d\n", model->trj_swap);
printf("total atoms: %d\n", num_atoms);
printf("size of integer: %d\n", sizeof(int));
printf("size of double: %d\n", sizeof(double));
#endif

return(0);
}

/**********************************/
/* GULP dynamics trajectory frame */
/**********************************/
#define DEBUG_TRJ_FRAME 0
gint read_trj_frame(FILE *fp, struct model_pak *model, gint replace_coords)
{
gint i, j;
int num_atoms;
#if DEBUG_TRJ_FRAME
int periodic;
#endif
double time, ke, pe, temp, *x[6];
double cell[9], velc[6];
GSList *list;
struct core_pak *core;
struct shel_pak *shell;

/* standard frame read */
num_atoms = model->expected_cores + model->expected_shells;
#if DEBUG_TRJ_FRAME
periodic = model->periodic;
#endif

/* alloc for x,y,z & vx, vy, vz */
for (j=0 ; j<6 ; j++)
  x[j] = g_malloc(num_atoms * sizeof(double));

/* read one frame */
READ_RECORD;

if (!fread(&time, sizeof(time), 1, fp))
  print_file_error(fp);
if (model->trj_swap)
  swap_bytes(&time, sizeof(time));

if (!fread(&ke, sizeof(ke), 1, fp))
  print_file_error(fp);
if (model->trj_swap)
  swap_bytes(&ke, sizeof(ke));

if (!fread(&pe, sizeof(pe), 1, fp))
  print_file_error(fp);
if (model->trj_swap)
  swap_bytes(&pe, sizeof(pe));

if (!fread(&temp, sizeof(temp), 1, fp))
  print_file_error(fp);
if (model->trj_swap)
  swap_bytes(&temp, sizeof(temp));

READ_RECORD;

#if DEBUG_TRJ_FRAME
printf("[%1dD][%d atoms]", periodic, num_atoms);
printf("[t : %.2lf]",time);
printf("[ke : %.2lf]",ke);
printf("[pe : %.2lf]",pe);
printf("[T : %.2lf][%d]\n",temp, replace_coords);
#endif

/* record dynamics frame */
model->frame_time = time;
model->frame_ke = ke;
model->frame_pe = pe;
model->frame_temp = temp;

/* loop over x,y,z & assoc velocity components */
for (j=0 ; j<6 ; j++)
  {
/* NB: cores first, then shells */
  READ_RECORD;
  for (i=0 ; i<num_atoms ; i++)
    {
    if (!fread(x[j]+i, sizeof(double), 1, fp))
      {
      print_file_error(fp);
      break;
      }
    if (model->trj_swap)
      swap_bytes(x[j]+i, sizeof(double));
    }
  READ_RECORD;
  }

#if DEBUG_TRJ_FRAME
printf("expected: cores = %d, shells = %d\n",
        model->expected_cores, model->expected_shells);
printf("existing: cores = %d, shells = %d\n",
        g_slist_length(model->cores), g_slist_length(model->shels));
#endif

/* read core data */
list = model->cores;
if (replace_coords)
  {
  for (i=0 ; i<model->expected_cores ; i++)
    {
    if (list)
      core = list->data;
    else
      {
      core = new_core("X", model);
      model->cores = g_slist_append(model->cores, core);
      list = g_slist_last(model->cores);
      }

    core->x[0] = *(x[0]+i);
    core->x[1] = *(x[1]+i);
    core->x[2] = *(x[2]+i);

    core->v[0] = *(x[3]+i);
    core->v[1] = *(x[4]+i);
    core->v[2] = *(x[5]+i);
  
    list = g_slist_next(list);
    }
  }

/* read shells data */
list = model->shels;
if (replace_coords)
  {
  for (i=model->expected_cores ; i<num_atoms ; i++)
    {
    if (list)
      shell = list->data;
    else
      {
      shell = new_shell("X", model);
      model->shels = g_slist_append(model->shels, shell);
      list = g_slist_last(model->shels);
      }

    shell->x[0] = *(x[0]+i);
    shell->x[1] = *(x[1]+i);
    shell->x[2] = *(x[2]+i);
  
    shell->v[0] = *(x[3]+i);
    shell->v[1] = *(x[4]+i);
    shell->v[2] = *(x[5]+i);
  
    list = g_slist_next(list);
    }
  }

/* read cell info */
if (model->gulp.ensemble == NPT)
  {
/* get cell vectors */
  READ_RECORD;
  fread(cell, sizeof(double), 9, fp);
  READ_RECORD;

/* get cell velocities */
  READ_RECORD;
/* calculate nstrains [0, 6] */
  switch (model->periodic)
    {
    case 3:
      j=6;
      break;
    case 2:
      j=3;
      break;
    case 1:
      j=1;
      break;
    default:
      j=0;
    }
  fread(velc, sizeof(double), j, fp);
  READ_RECORD;

/* compute latmat/ilatmat */
  memcpy(model->latmat, cell, 9*sizeof(gdouble));
  model->construct_pbc = TRUE;
  matrix_lattice_init(model); 
  }

/* convert input cartesian to fractional */
if (replace_coords)
  if (model->fractional)
    coords_make_fractional(model);

/* clean up */
for (j=0 ; j<6 ; j++)
  g_free(x[j]);

return(0);
}

/**********************************************/
/* process a TRJ file and get frame positions */
/**********************************************/
gint mark_trj_frames(struct model_pak *model)
{
gint i;
FILE *fp;
gchar *filename;

filename = g_strdup_printf("%s/%s", sysenv.cwd, model->gulp.trj_file);

/* NB: GULP binary trajectory files - MUST specify "rb" */
/* FIXME - some GULP traj files may be ascii I think? */
fp = fopen(filename, "rb");
if (!fp)
  {
  printf("Failed to open: %s\n", filename);
  g_free(filename);
  return(1);
  }

read_trj_header(fp, model);

for (i=0 ; i<model->num_frames ; i++)
  {
  add_frame_offset(fp, model);
  read_trj_frame(fp, model, FALSE);
  }

fclose(fp);

g_free(filename);

return(0);
}

/***********************************/
/* GULP dynamics trajectory header */
/***********************************/
gint write_trj_header(FILE *fp, struct model_pak *model)
{
int num_atoms, periodic;
double version = 1.4;

WRITE_RECORD;
fwrite(&version, sizeof(version), 1, fp);
WRITE_RECORD;

WRITE_RECORD;
num_atoms = g_slist_length(model->cores) + g_slist_length(model->shels);
fwrite(&num_atoms, sizeof(num_atoms), 1, fp);
periodic = model->periodic;
fwrite(&periodic, sizeof(periodic), 1, fp);
WRITE_RECORD;

return(0);
}

/****************************************/
/* write GULP dynamics trajectory frame */
/****************************************/
gint write_trj_frame(FILE *fp, struct model_pak *model)
{
gint i, j;
double temp, x[3];
GSList *list;
struct core_pak *core;
struct shel_pak *shell;

/* fake time, ke, pe, and temp */
WRITE_RECORD;
temp = 0.0;
fwrite(&temp, sizeof(double), 4, fp);
WRITE_RECORD;

/* core/shell coordinate write */
for (i=0 ; i<3 ; i++)
  {
  WRITE_RECORD;
  for (list=model->cores ; list ; list=g_slist_next(list))
    {
    core = list->data;

    ARR3SET(x, core->x);
    vecmat(model->latmat, x);
  
    fwrite(&x[i], sizeof(double), 1, fp);
    }
  for (list=model->shels ; list ; list=g_slist_next(list))
    {
    shell = list->data;

    ARR3SET(x, shell->x);
    vecmat(model->latmat, x);
  
    fwrite(&x[i], sizeof(double), 1, fp);
    }
  WRITE_RECORD;
  }

/* core/shell velocity write */
for (i=0 ; i<3 ; i++)
  {
  WRITE_RECORD;
  for (list=model->cores ; list ; list=g_slist_next(list))
    {
    core = list->data;

    x[i] = core->v[i];
  
    fwrite(&x[i], sizeof(double), 1, fp);
    }
  for (list=model->shels ; list ; list=g_slist_next(list))
    {
    shell = list->data;

    x[i] = shell->v[i];
  
    fwrite(&x[i], sizeof(double), 1, fp);
    }
  WRITE_RECORD;
  }

/* lattice vector write */
if (model->gulp.ensemble == NPT)
  {
  WRITE_RECORD;
/* TODO - check matrix storage order is ok (ie in rows or columns?) */
  for (i=0 ; i<9 ; i++)
    {
    temp = model->latmat[i];
    fwrite(&temp, sizeof(double), 1, fp);
    }
  WRITE_RECORD;

/* calculate nstrains [0, 6] */
  switch (model->periodic)
    {
    case 3:
      j=6;
      break;
    case 2:
      j=3;
      break;
    case 1:
      j=1;
      break;
    default:
      j=0;
    }

/* fake lattice velocity write */
  WRITE_RECORD;
  temp = 0.0;
  for (i=0 ; i<j ; i++)
    fwrite(&temp, sizeof(double), 1, fp);
  WRITE_RECORD;
  }

return(0);
}

/*******************/
/* run a gulp file */
/*******************/
#define DEBUG_EXEC_GULP 0
gint exec_gulp(const gchar *input, const gchar *output)
{
gint status=0;
gchar *cmd;

/* checks */
if (!sysenv.gulp_path)
  return(-1);

/* delete the old file to be sure output is only data from current run */
chdir(sysenv.cwd);
unlink(output);

/* put the GULP path in quotes to avoid problems with spaces in pathnames etc */
/* TODO - need to do the same with input/output names? */

#ifdef _WIN32
cmd = g_strdup_printf("\"%s\" < %s > %s", sysenv.gulp_path, input, output);
#else
cmd = g_strdup_printf("%s < %s > %s", sysenv.gulp_path, input, output);
#endif


#if DEBUG_EXEC_GULP
printf("executing: [%s]\n",cmd);
#endif

task_sync(cmd);

/* done */
g_free(cmd);
return(status);
}
