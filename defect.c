/*
Copyright (C) 2003 by Sean David Fleming and Andrew Walker

sean@ivec.org
andrew.m.walker@anu.edu.au

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
#include <unistd.h>
#include <math.h>
#include <time.h>

#include "gdis.h"
#include "file.h"
#include "scan.h"
#include "parse.h"
#include "coords.h"
#include "defect.h"
#include "edit.h"
#include "graph.h"
#include "dialog.h"
#include "model.h"
#include "matrix.h"
#include "geometry.h"
#include "numeric.h"
#include "space.h"
#include "surface.h"
#include "task.h"
#include "zone.h"
#include "interface.h"
#include "gui_shorts.h"

extern struct sysenv_pak sysenv;

/*************************/
/* check for 2-fold axes */
/*************************/
void defect_check_symmetry(struct model_pak *model)
{
gint i, n;

for (i=1 ; i<model->sginfo.order ; i++)
  {

  n = matrix_is_z_rotation(*(model->sginfo.matrix+i));

  if (n)
    {

printf("matrix %d", i);
P3MAT(" ", *(model->sginfo.matrix+i));
P3VEC("+", *(model->sginfo.offset+i));
printf(" : is a %d rotation axis\n", n);

    }

  }
}

/*****************************************/
/* acquire the elastic compliance matrix */
/*****************************************/
gint defect_compliance_get(gdouble *s, struct plane_pak *plane, struct model_pak *source)
{
gint i, j, tokens;
gchar *inp, *out, *line, **buff;
gchar *full_inp, *full_out;
gpointer scan;
struct model_pak *model;

g_assert(source != NULL);

/* I thought we could use the true cell stuff, but it's too messy */
/* given that the dislocation is based on distance from cylinder center */

printf("generating transformed cell.\n");

model = model_new();
source->surface.true_cell = TRUE;
ARR3SET(model->surface.miller, plane->m);
model->surface.region[0] = 1;
model->surface.region[1] = 0;
model->surface.shift = 0.0;
generate_surface(source, model);

/* setup "surface" for GULP calc */
gulp_data_copy(source, model);

/* this can have crap in it that stuffs up the GULP run */
/*
gulp_extra_copy(source, model);
*/

/* run the required GULP job */
model->gulp.prop = TRUE;
model->gulp.run = E_OPTIMIZE;

/* this tells GULP to output lattice vectors */
model->construct_pbc = TRUE;
model->periodic = 3;

/* get a temporary name */
inp = gun("gin");
out = gun("got");
full_inp = g_build_filename(sysenv.cwd, inp, NULL);
full_out = g_build_filename(sysenv.cwd, out, NULL);
g_free(inp);
g_free(out);

printf("acquiring elastic compliance matrix...\n");

write_gulp(full_inp, model);
if (exec_gulp(full_inp, full_out))
  return(1);

/* process the output */
scan = scan_new(full_out);
if (!scan)
  return(2);
 
while (!scan_complete(scan))
  {
  line = scan_get_line(scan);
  
  if (line)
  if (g_ascii_strncasecmp(line, "  Elastic Compliance Matrix:", 28) == 0)
    {
    for (i=4 ; i-- ; )
      scan_get_line(scan);
    
    for (i=0 ; i<6 ; i++)
      {
      buff = scan_get_tokens(scan, &tokens);
      if (tokens > 6)
        {
        for (j=0 ; j<6 ; j++)
          s[6*i+j] = str_to_float(*(buff+j+1));
        }
      else
        printf("Error parsing constants.\n");

      g_strfreev(buff);
      }
    } 
  }
scan_free(scan);

unlink(full_inp);
unlink(full_out);

g_free(full_inp);
g_free(full_out);

model_free(model);
sysenv.mal = g_slist_remove(sysenv.mal, model);
g_free(model);

return(0);
}

/********************************************************/
/* create a cylinder (xy cross section) of given radius */
/********************************************************/
void defect_make_cylinder(gdouble *o, gdouble r1, gdouble r2, struct model_pak *model)
{
gint region;
gdouble r, rsq, r1sq, d2, x[3];
GSList *mlist, *clist;
struct mol_pak *mol;
struct core_pak *core;
struct shel_pak *shell;

g_assert(model != NULL);

r = r1+r2;
rsq = r*r;
r1sq = r1*r1;

/* generate enough unit cells to satisfy input radius */
/* NB: round up the image number required */
model->image_limit[0] = (gint) (1 + r / model->pbc[0]);
model->image_limit[1] = (gint) (1 + r / model->pbc[0]);
model->image_limit[2] = (gint) (1 + r / model->pbc[1]);
model->image_limit[3] = (gint) (1 + r / model->pbc[1]);
model->image_limit[4] = 0.0;
model->image_limit[5] = 1.0;
space_make_images(CREATE, model);
space_make_supercell(model);

/* convert to cartesian block */
coords_make_cartesian(model);
model->periodic = 0;
model->fractional = FALSE;
coords_compute(model);
zone_init(model);

/* if we cleave by atoms - prevent bonding anything, so each molecule is an atom */
/* NB: not sure why, but source's bonds are transferred to model */
if (model->surface.ignore_bonding)
  wipe_bonds(model);
else
  connect_bonds(model);

connect_molecules(model);

/* offset by desired amount and clip */
/* TODO - if (!cleave) else by atom */
for (mlist=model->moles ; mlist ; mlist=g_slist_next(mlist))
  {
  mol = mlist->data;

/* transform the centroid */
  ARR3SET(x, mol->centroid);

/* FIXME - o is fractional - x is not? */
  x[0] -= o[0];
  x[1] -= o[1];

  d2 = x[0]*x[0] + x[1]*x[1];

  if (d2 < rsq)
    {

    if (d2 < r1sq)
      region = 0; 
    else
      region = 1;

    for (clist=mol->cores ; clist ; clist=g_slist_next(clist))
      {
      core = clist->data;
  core->x[0] -= o[0];
  core->x[1] -= o[1];
      core->region = region;
      if (core->shell)
        {
        shell = core->shell;
  shell->x[0] -= o[0];
  shell->x[1] -= o[1];
        shell->region = region;
        }
      }
    }
  else
    {
    for (clist=mol->cores ; clist ; clist=g_slist_next(clist))
      {
      core = clist->data;
      delete_core(core);
      }
    }
  }
}

/******************************/
/* displace the supplied core */
/******************************/
gdouble defect_anisotropic_1(gdouble *x, gdouble r44, gdouble r55, gdouble *b)
{
gdouble angle;

angle = angle_x_compute(x[0], x[1]);
/* -ve ??? */
angle *= sqrt(fabs(r44)/fabs(r55));
return(0.5 * VEC3MAG(b) * angle / G_PI);
}

/************************************************/
/* do coordinate offsets for a scew dislocation */
/************************************************/
void defect_screw_dislocate(gdouble *s, gdouble *v, struct model_pak *model)
{
gdouble dz, r44, r55, x[3];
GSList *clist, *mlist;
struct mol_pak *mol;
struct core_pak *core;
struct shel_pak *shell;

g_assert(model != NULL);

/* r44 = s44 - (s34 s34 / s33) */
r44 = s[21] - ((s[15]*s[15])/s[14]);

/* r55 = s55 - (s35 s35 / s33) */
r55 = s[28] - ((s[16]*s[16])/s[14]); 

printf("\nReduced elastic compliances: S44 = %f, S55 = %f\n", r44, r55);
if (r44 == r55)
  printf("Isotropic case...\n");
 
if (fabs(r44) < FRACTION_TOLERANCE || fabs(r55) < FRACTION_TOLERANCE)
  {
  printf("Bad compliance matrix.\n"); 
  return;
  }

for (mlist=model->moles ; mlist ; mlist=g_slist_next(mlist))
  {
  mol = mlist->data;

  ARR3SET(x, mol->centroid);

  dz = defect_anisotropic_1(x, r44, r55, v);

  for (clist=mol->cores ; clist ; clist=g_slist_next(clist))
    {
    core = clist->data;
    core->x[2] += dz;
    if (core->shell)
      {
      shell = core->shell;
      shell->x[2] += dz;
      }
    }
  }
}

/*********************************************************************/
/* attempt to find an int which will multiply to make a float an int */
/*********************************************************************/
/* TODO - put in numeric */
gint numeric_make_whole(gdouble x)
{
gint i;
gdouble frac, f;

frac = fabs(x - nearest_int(x));

if (frac > FRACTION_TOLERANCE)
  {
/* crude hack - look for common extensions */
  for (i=2 ; i<9 ; i++)
    {
    f = 1.0 / (gdouble) i;
    if (fabs(frac-f) < FRACTION_TOLERANCE)
      return(i);
    }
  return(0);
  }
else
  return(nearest_int(x));
}

/********************/
/* build the defect */
/********************/
#define DEBUG_DEFECT_NEW 1
void defect_new(struct defect_pak *defect, struct model_pak *source)
{
gint height, gcd, mult1, mult2;
gdouble r, len;
gdouble o[2], v[3], b[3], m[9], s[36];
gdouble w[3], x[3], x1[3], x2[3], x3[3];
GSList *list;
struct model_pak *model;
struct plane_pak *plane;
struct core_pak *core;
struct shel_pak *shell;

g_assert(defect != NULL);
g_assert(source != NULL);

/* check for 2-fold axes */
defect_check_symmetry(source);

if (VEC3MAGSQ(defect->orient) < FRACTION_TOLERANCE)
  {
  gui_text_show(ERROR, "Please specify a non-zero dislocation vector.\n");
  return;
  }

/* setup the defect geometry */
ARR3SET(v, defect->orient);
ARR3SET(b, defect->burgers);

vecmat(source->latmat, v);
vecmat(source->latmat, b);
matrix_z_alignment(m, v);
vecmat(m, v);
vecmat(m, b);

#if DEBUG_DEFECT_NEW
P3VEC("input dislocation: ", defect->orient);
P3VEC("input burgers: ", defect->burgers);
P3VEC("transformed dislocation: ", v);
P3VEC("transformed burgers: ", b);
#endif

/* get the transformed lattice vectors, so we know how deep the cell needs to be */
plane = plane_new(defect->orient, source);

VEC3SET(x1, plane->lattice[0], plane->lattice[3], plane->lattice[6]);
VEC3SET(x2, plane->lattice[1], plane->lattice[4], plane->lattice[7]);
VEC3SET(x3, plane->lattice[2], plane->lattice[5], plane->lattice[8]);

P3MAT("cell: ", plane->lattice);

printf("%f x %f\n", VEC3MAG(x1), VEC3MAG(x2));

mult1 = mult2 = 1;

/* check if periodicity along first surface vector needs help */
ARR3SET(x, x3);
vector_v_project(x, x1);
len = VEC3MAG(x);
if (len > FRACTION_TOLERANCE)
  {
printf("non-zero projection of depth onto surface vec 1 = %f\n", len);

  r = VEC3MAG(x1)/len;
  mult1 = numeric_make_whole(r);
  if (!mult1)
    mult1 = 1;
  }

/* check if periodicity along second surface vector needs help */
ARR3SET(x, x3);
vector_v_project(x, x2);
len = VEC3MAG(x);
if (len > FRACTION_TOLERANCE)
  {
printf("non-zero projection of depth onto surface vec 2 = %f\n", len);

  r = VEC3MAG(x2)/len;
  mult2 = numeric_make_whole(r);
  if (!mult2)
    mult2 = 1;
  }

/* compute compiance tensor */
if (defect_compliance_get(s, plane, source))
  {
  printf("Failed - terminating.\n");
  return;
  }
else
  {
printf("s44 = %f, s55 = %f\n", s[21], s[28]);
  }

gcd = GCD(mult1, mult2);
height = mult1 * mult2 / gcd;

printf("mult1 = %d, mult2 = %d, gcd = %d\n", mult1, mult2, gcd);
printf("height = %d\n", height);

/* DEBUG - test that our new depth vector (in z) will hit a lattice point */
ARR3SET(w, x3);
VEC3MUL(w, height);
vector_v_project(w, x1);
printf("lattice point along surface vector 1: %f\n", VEC3MAG(w)/VEC3MAG(x1));

ARR3SET(w, x3);
VEC3MUL(w, height);
vector_v_project(w, x2);
printf("lattice point along surface vector 2: %f\n", VEC3MAG(w)/VEC3MAG(x2));

/* CURRENT - cope with reduced orientation vector (eg 220 = halved 110) */
gcd = GCD(GCD(defect->orient[0], defect->orient[1]), GCD(defect->orient[1], defect->orient[2]));
height /= gcd;

/* TODO - if threaded - do this in the gui *before* defect_new() is called */

model = model_new();

/* I thought we could use the true cell stuff, but it's too messy */
/* given that the dislocation is based on distance from cylinder center */
source->surface.true_cell = FALSE;
ARR3SET(model->surface.miller, defect->orient);
model->surface.region[0] = height;
model->surface.region[1] = 0;
model->surface.shift = 0.0;
generate_surface(source, model);

/* model is fractional - ie the transformed unit cell */
if (defect->cleave)
  model->surface.ignore_bonding = TRUE;
matrix_lattice_init(model);


/* convert fractional origin shift to cartesian */

printf("f origin: %f,%f\n", defect->origin[0], defect->origin[1]);

o[0] = plane->lattice[0] * defect->origin[0] + plane->lattice[1] * defect->origin[1];
o[1] = plane->lattice[3] * defect->origin[0] + plane->lattice[4] * defect->origin[1];

printf("c origin: %f,%f\n", o[0], o[1]);

defect_make_cylinder(o, defect->region[0], defect->region[1], model);

/* coordinate dislocation */
if (via(v, b, 3) < 0.001)
  {
printf("creating SCREW dislocation...\n");
defect_screw_dislocate(s, b, model);
  }
else
  {
printf("TODO - create EDGE dislocation...\n");
  }

/* CURRENT - rotate so that z (ie defect orientation) points along x */
VEC3SET(&m[0],  0.0, 0.0, 1.0);
VEC3SET(&m[3],  0.0, 1.0, 0.0);
VEC3SET(&m[6], -1.0, 0.0, 0.0);
for (list=model->cores ; list ; list=g_slist_next(list))
  {
  core = list->data;
  vecmat(m, core->x);
  if (core->shell)
    {
    shell = core->shell;
    vecmat(m, shell->x);
    }
  }

/* display the model - if threaded have to shunt this elsewhere */
sysenv.mal = g_slist_append(sysenv.mal, model);

/* all solated cluster for debugging purposes */
if (defect->cluster)
  {
  matrix_identity(model->latmat);
  model->periodic = 0;
  }
else
  {
/* setup the lattice matrix */
  matrix_identity(model->latmat);
  model->periodic = 1;

/* get z projected length of depth vector */
  VEC3SET(x, 0.0, 0.0, 1.0);
  vector_v_project(x3, x);
  model->latmat[0] = height * VEC3MAG(x3);
  }


plane_free(plane);


model->fractional = FALSE;
model->id = GULP;
gulp_data_copy(source, model);

model_prep(model);
sysenv.active_model = model;
tree_model_add(model);
tree_select_active();
}
