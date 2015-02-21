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

#include <ctype.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#ifdef __APPLE__
#include <OpenGL/gl.h>
#else
#include <GL/gl.h>
#endif

#include "gdis.h"
#include "coords.h"
#include "edit.h"
#include "file.h"
#include "parse.h"
#include "library.h"
#include "matrix.h"
#include "zmatrix.h"
#include "measure.h"
#include "model.h"
#include "morph.h"
#include "numeric.h"
#include "select.h"
#include "space.h"
#include "spatial.h"
#include "surface.h"
#include "gui_shorts.h"
#include "type.h"
#include "interface.h"
#include "dialog.h"
#include "opengl.h"
#include "zone.h"

#define DEBUG 0

extern struct sysenv_pak sysenv;
extern struct elem_pak elements[];
extern GtkWidget *apd_label;
extern GtkWidget *window;

/* globals */
gdouble tmat[9];
gdouble tvec[3];
gdouble edit_anim_n=0;
gdouble edit_chirality[2] = {6, 6};
gdouble edit_length = 1.44;
gpointer edit_construct;
gchar *edit_basis[2] = {NULL, NULL};
GtkWidget *elem_entry, *grid_dim, *grid_sep;
GtkWidget *axis_entry, *obj_entry, *periodicity_spin;
GtkWidget *transmat[12];

/* NEW */
gdouble edit_spatial_colour[3] = {1.0, 0.0, 0.0};
gdouble edit_label_colour[3] = {1.0, 1.0, 1.0};

enum{AT_ANY, AT_SELECTED, AT_LATTICE};

/*********************************************/
/* update widget values with internal values */
/*********************************************/
void update_transmat(void)
{
gint i, j;
gchar *text;

/* set transformation matrix widget entries */
for (j=0 ; j<3 ; j++)
  {
  for (i=0 ; i<3 ; i++)
    {
    text = g_strdup_printf("%f", tmat[3*j+i]);
/* NB: display transpose */
    gtk_entry_set_text(GTK_ENTRY(transmat[3*j+i]), text);
    g_free(text);
    }
  }

/* set translation widget entries */
for (i=0 ; i<3 ; i++)
  {
  text = g_strdup_printf("%f", tvec[i]);
  gtk_entry_set_text(GTK_ENTRY(transmat[9+i]), text);
  g_free(text);
  }
}

/********************************************/
/* restore original transformation settings */
/********************************************/
void reset_transmat(void)
{
VEC3SET(tvec, 0.0, 0.0, 0.0);
matrix_identity(tmat);
update_transmat();
}

/*************************************/
/* construct a transformation matrix */
/*************************************/
void construct_transmat(void)
{
gint n, flag=0, type=-1;
const gchar *text;
gdouble angle;
gdouble v1[3];
struct vec_pak *p[3];
struct spatial_pak *spatial=NULL;
struct model_pak *data;

/* checks */
data = sysenv.active_model;
if (!data)
  return;
if (!GTK_IS_ENTRY(edit_construct))
  return;

/* requested transformation */
text = gtk_entry_get_text(GTK_ENTRY(edit_construct));
if (g_ascii_strncasecmp(text, "z alignment", 11) == 0)
  type = ALIGNMENT;
if (g_ascii_strncasecmp(text, "reflection", 10) == 0)
  type = REFLECTION;
if (g_ascii_strncasecmp(text, "rotation", 8) == 0)
  type = PAXIS;
if (g_ascii_strncasecmp(text, "lattice", 7) == 0)
  type = LATMAT;
if (g_ascii_strncasecmp(text, "identity", 8) == 0)
  type = IDENTITY;

/* special cases */
switch(type)
  {
  case IDENTITY:
    reset_transmat();
    return;

  case LATMAT:
    memcpy(tmat, data->latmat, 9*sizeof(gdouble));
    VEC3SET(tvec, 0.0, 0.0, 0.0);
    update_transmat();
    return;
  }

/* rotation angle */
text = gtk_entry_get_text(GTK_ENTRY(axis_entry));
angle = str_to_float(text);

/* check for special reference objects */
text = gtk_entry_get_text(GTK_ENTRY(obj_entry));
if (g_ascii_strncasecmp("x", text, 1) == 0)
  flag = 1;
if (g_ascii_strncasecmp("y", text, 1) == 0)
  flag = 2;
if (g_ascii_strncasecmp("z", text, 1) == 0)
  flag = 3;
if (g_ascii_strncasecmp("a", text, 1) == 0)
  flag = 4;
if (g_ascii_strncasecmp("b", text, 1) == 0)
  flag = 5;
if (g_ascii_strncasecmp("c", text, 1) == 0)
  flag = 6;

/* construct appropriate reference vector */
switch (flag)
  {
  case 1:
    VEC3SET(v1, 1.0, 0.0, 0.0);
    break;
  case 2:
    VEC3SET(v1, 0.0, 1.0, 0.0);
    break;
  case 3:
    VEC3SET(v1, 0.0, 0.0, 1.0);
    break;
  case 4:
    VEC3SET(v1, 1.0, 0.0, 0.0);
    vecmat(data->latmat, v1);
    break;
  case 5:
    VEC3SET(v1, 0.0, 1.0, 0.0);
    vecmat(data->latmat, v1);
    break;
  case 6:
    VEC3SET(v1, 0.0, 0.0, 1.0);
    vecmat(data->latmat, v1);
    break;

  default:
/* no special reference object - assume a numerical spatial reference */
    n = str_to_float(text);
    spatial = g_slist_nth_data(data->spatial, n);
    if (!spatial)
      {
      gui_text_show(ERROR, "Undefined reference spatial object.\n");
      return;
      }

/* compute orientation vector (plane normal/vector direction) */
    switch (spatial->type)
      {
      case SPATIAL_VECTOR:
        p[0] = g_slist_nth_data(spatial->list, 0);
        p[1] = g_slist_nth_data(spatial->list, 1);
        g_assert(p[0] != NULL);
        g_assert(p[1] != NULL);
/* vector points from 1st to 2nd */
        ARR3SET(v1, p[1]->rx);
        ARR3SUB(v1, p[0]->rx);
        break;

      default:
        p[0] = g_slist_nth_data(spatial->list, 0);
        p[1] = g_slist_nth_data(spatial->list, 1);
        p[2] = g_slist_nth_data(spatial->list, 2);
        g_assert(p[0] != NULL);
        g_assert(p[1] != NULL);
        g_assert(p[2] != NULL);
        calc_norm(v1, p[0]->rx, p[1]->rx, p[2]->rx);
        break;
      }
    break;
  }

/* construct */
switch (type)
  {
  case ALIGNMENT:
    matrix_z_alignment(tmat, v1);
    break;

  case PAXIS:
    matrix_v_rotation(tmat, v1, D2R*angle);
    break;

  case REFLECTION:
    matrix_v_reflection(tmat, v1);
    break;
  }
VEC3SET(tvec, 0.0, 0.0, 0.0);

update_transmat();
}

/********************************************/
/* put text entry values into actual matrix */
/********************************************/
void change_transmat(GtkWidget *w, gint i)
{
const gchar *text;

text = gtk_entry_get_text(GTK_ENTRY(transmat[i]));
if (i < 9)
  tmat[i] = str_to_float(text);
else
  tvec[i-9] = str_to_float(text);
}

/************************************************/
/* apply the current transformation/translation */
/************************************************/
void apply_transmat(GtkWidget *w, gint mode)
{
GSList *item, *list=NULL;
struct model_pak *data;
struct core_pak *core;
struct shel_pak *shel;

data = sysenv.active_model;
if (!data)
  return;

/* application mode */
switch(mode)
  {
  case AT_ANY:
    list = data->cores;
    break;

  case AT_SELECTED:
    list = data->selection;
    break;

  case AT_LATTICE:
    matrix_lattice_new(tmat, data);
    return;

  default:
    g_assert_not_reached();
  }

for (item=list ; item ; item=g_slist_next(item))
  {
  core = item->data;

/* transform */
  vecmat(data->latmat, core->x);
  vecmat(tmat, core->x);
  ARR3ADD(core->x, tvec);
  vecmat(data->ilatmat, core->x);

/* assoc. shell? */
  if (core->shell)
    {
    shel = core->shell;

/* transform */
    vecmat(data->latmat, shel->x);
    vecmat(tmat, shel->x);
    ARR3ADD(shel->x, tvec);
    vecmat(data->ilatmat, shel->x);
    }
  }
/* update */
coords_compute(data);
zone_init(data);
connect_bonds(data);
connect_molecules(data);

redraw_canvas(SINGLE);
}

/***************************************************/
/* construct animation using multiple applications */
/***************************************************/
/* FIXME - this is all broken (due to new quat based camera) */
void apply_n_transmat(GtkWidget *w, gint mode)
{
/*
gint i;
gdouble mat[9];
struct model_pak *model;
struct transform_pak *transform;

model = sysenv.active_model;
if (!model)
  return;

dialog_destroy_single(ANIM, model);

memcpy(mat, model->rotmat, 9*sizeof(gdouble));

for (i=0 ; i<edit_anim_n ; i++)
  {
  transform = g_malloc(sizeof(struct transform_pak));
  model->transform_list = g_slist_append(model->transform_list, transform);
  transform->id = ROTATION;

transpose(mat);
  matmat(tmat, mat);
transpose(mat);
  memcpy(transform->matrix, mat, 9*sizeof(gdouble));

  memcpy(transform->matrix, tmat, 9*sizeof(gdouble));

  VEC3SET(transform->vector, model->offset[0], model->offset[1], 0.0);
  transform->scalar = model->scale;
  model->num_frames++; 
  }

if (model->num_frames)
  model->animation = TRUE;

redraw_canvas(SINGLE);
*/
}

/****************************/
/* move the whole selection */
/****************************/
#define DEBUG_REGION_MOVE 0
void region_move(GtkWidget *w, gint direction)
{
gchar txt[60];
GSList *list;
struct model_pak *data;
struct core_pak *core;

data = sysenv.active_model;

/* checks */
if (!data)
  return;
if (data->periodic != 2)
  return;
if (!data->selection)
  {
  gui_text_show(WARNING, "Empty selection.\n");
  return;
  }

#if DEBUG_REGION_MOVE 
P3MAT("depth vec:", data->surface.depth_vec);
#endif

/* quit if no depth info ie a loaded 2D file */
if (VEC3MAG(data->surface.depth_vec) < FRACTION_TOLERANCE)
  {
  gui_text_show(ERROR,
  "This is only possible for surfaces generated in the current session.\n");
  return;
  }

/* move all atoms in selection */
for (list=data->selection ; list ; list=g_slist_next(list))
  {
  core = list->data;
  region_move_atom(core, direction, data); 
  }

/* update */
coords_compute(data);
connect_bonds(data);
connect_molecules(data);

#if DEBUG_REGION_MOVE 
printf("old dipole: %f\n", data->gulp.sdipole);
#endif

/* recompute */
calc_emp(data);

#if DEBUG_REGION_MOVE 
printf("new dipole: %f\n", data->gulp.sdipole);
#endif

/* inform the user */
sprintf(txt, "New surface dipole: %f\n", data->gulp.sdipole);
gui_text_show(STANDARD, txt);

redraw_canvas(SINGLE);
}

/****************************/
/* connectivity only update */
/****************************/
void update_bonds(void)
{
struct model_pak *data;

data = sysenv.active_model;
if (!data)
  return;

coords_compute(data);
connect_bonds(data);
connect_molecules(data);
redraw_canvas(SINGLE);
}

/*******************************/
/* creates a new (blank) model */
/*******************************/
void edit_model_create(void)
{
struct model_pak *data;

/* make a slot for the new model */
data = model_new();
sysenv.active_model = data;

/* setup model parameters */
data->id = CREATOR;
strcpy(data->filename, "new_model");
g_free(data->basename);
data->basename = g_strdup("new_model");

/* initialize */
model_prep(data);
data->mode = FREE;

/* rmax has to be after coords_init() as INIT_COORDS will set it to 1.0 */
data->rmax = 5.0*RMAX_FUDGE;

/* update/redraw */
tree_model_add(data);
tree_select_model(data);
redraw_canvas(SINGLE);
}

/***********************/
/* add atom event hook */
/***********************/
void add_atom(gint x, gint y, struct model_pak *data)
{
gchar *elem;
const gchar *orig;
gdouble r[3];
struct core_pak *core;

g_assert(data != NULL);

/* get the atom type from the label */
orig = gtk_entry_get_text(GTK_ENTRY(apd_label));
elem = g_strdup(orig);

/* add the atom to the model */
core = new_core(elem, data);

/* set atom position in the plane running through the origin */
gl_project(r, x, y, canvas_find(data));

/*
ARR3ADD(r, data->centroid);
VEC3MUL(r, 1.0/data->scale);
*/

ARR3SET(core->rx, r);

ARR3SET(core->x, r);
vecmat(data->ilatmat, core->x);
ARR3ADD(core->x, data->centroid);

data->cores = g_slist_append(data->cores, core);

/* make sure this atom is in the apd box, so we can keep adding this type */
select_clear(data);
select_add_core(core, data);

/* TODO - need a more fine grained update (atoms/selection) */
/*
coords_compute(data);
*/
zone_init(data);
connect_bonds(data);
connect_molecules(data);

/* REFRESH */
gui_refresh(GUI_MODEL_PROPERTIES);

g_slist_free(data->unique_atom_list);
data->unique_atom_list = find_unique(ELEMENT, data);
init_atom_colour(core, data);
init_atom_charge(core, data);
calc_emp(data);

/* done */
g_free(elem);
}

/*****************************/
/* change the lattice matrix */
/*****************************/
void edit_transform_latmat(void)
{
GSList *list;
struct model_pak *model;
struct core_pak *core;
struct shel_pak *shell;

/* checks */
model = sysenv.active_model;
if (!model)
  return;

/* remove all symmetry and convert to cartesian */
space_make_p1(model);

/* make coordinates cartesian for the prep stqage */
for (list=model->cores ; list ; list=g_slist_next(list))
  {
  core = list->data;
  vecmat(model->latmat, core->x);
  }
for (list=model->shels ; list ; list=g_slist_next(list))
  {
  shell = list->data;
  vecmat(model->latmat, shell->x);
  }

/* create new lattice */
matmat(tmat, model->latmat);
model->fractional = FALSE;
model->construct_pbc = TRUE;
model_prep(model);

/* GUI updates */
dialog_destroy_single(GENSURF, model);
tree_model_refresh(model);
gui_active_refresh();
gui_relation_update(model);
redraw_canvas(SINGLE);
}

/***********************************/
/* change periodicity of the model */
/***********************************/
void cb_modify_periodicity(void)
{
gint i, n, button;
gdouble d, x[3];
GSList *list;
GtkWidget *dialog;
struct model_pak *model;
struct core_pak *core;
struct shel_pak *shell;

/* checks */
model = sysenv.active_model;
if (!model)
  return;

n = SPIN_IVAL(GTK_SPIN_BUTTON(periodicity_spin));

/* if we're increasing the periodicity, we need repeat vector info */
if (n > model->periodic)
  {
/* check if its the unit vector (ie likely unchanged from default) */
  for (i=model->periodic ; i<n ; i++)
    {
    VEC3SET(x, tmat[i], tmat[i+3], tmat[i+6]);
    d = fabs(VEC3MAGSQ(x) - 1.0);
    if (d < FRACTION_TOLERANCE)
      {
      dialog = gtk_message_dialog_new(GTK_WINDOW(window),
                                      GTK_DIALOG_DESTROY_WITH_PARENT,
                                      GTK_MESSAGE_WARNING,
                                      GTK_BUTTONS_YES_NO,
                                     "Are you sure the transformation matrix contains the appropriate repeat vectors for this change in periodicity?");

      button = gtk_dialog_run(GTK_DIALOG(dialog));
      gtk_widget_destroy(dialog);
      if (button != -8)
        return;
      else
        break;
      }
    }
  }

/* remove all symmetry and convert to cartesian */
space_make_p1(model);

/* make coordinates cartesian for the prep stqage */
for (list=model->cores ; list ; list=g_slist_next(list))
  {
  core = list->data;
  vecmat(model->latmat, core->x);
  }
for (list=model->shels ; list ; list=g_slist_next(list))
  {
  shell = list->data;
  vecmat(model->latmat, shell->x);
  }

/* set the (extra) new lattice matrix periodicty */
for (i=model->periodic ; i<n ; i++)
  {
  model->latmat[i+0] = tmat[i+0];
  model->latmat[i+3] = tmat[i+3];
  model->latmat[i+6] = tmat[i+6];
  }

/* init new lattice */
model->periodic = n;
model->fractional = FALSE;
model->construct_pbc = TRUE;
model_prep(model);

/* GUI updates */
/* REFRESH */
dialog_destroy_single(GENSURF, model);
tree_model_refresh(model);

gui_relation_update(model);

gui_refresh(GUI_MODEL_PROPERTIES);
gui_refresh(GUI_CANVAS);
}

/********************/
/* build a nanotube */
/********************/
#define EDIT_NANOTUBE_NEW 0
void edit_nanotube_new(void)
{
gint i, j, k, n, m, d, dr, np, mp;
gint imin, imax, jmin, jmax, dummy[3];
gdouble a1[3], a2[3], a3[3];
gdouble a, r, x[3], y[3], ch[3], t[3], v[3], p[2];
gchar *text;
struct core_pak *core;
struct model_pak *model;

/* checks */
model = sysenv.active_model;
if (!model)
  {
  edit_model_create();
  model = sysenv.active_model;
  g_assert(model != NULL);
  }

/* default to brenner if potential set is absent */
if (!model->gulp.potentials && !model->gulp.libfile)
  model->gulp.potentials = g_strdup("brenner\n");

#define ROOT3 sqrt(3)

/* TODO - test basis atom strings for validity */
#if EDIT_NANOTUBE_NEW
printf("Basis: %s, %s (%f)\n", edit_basis[0], edit_basis[1], edit_length);
#endif

/* compute lattice basis vectors */
VEC3SET(a1, 1.5*edit_length, 0.5*ROOT3*edit_length, 0.0);
VEC3SET(a2, 1.5*edit_length, -0.5*ROOT3*edit_length, 0.0);
VEC3SET(a3, edit_length, 0.0, 0.0);

/* compute indices */
n = edit_chirality[0];
m = edit_chirality[1];

d = gcd(n, m);
if ((n - m) % (3*d))
  dr = d;
else
  dr = 3.0*d;

np = 2*m + n;
np /= dr;
mp = 2*n + m;
mp /= dr;

/* chirality vector */
ARR3SET(x, a1);
ARR3SET(y, a2);
VEC3MUL(x, n);
VEC3MUL(y, m);
ARR3SET(ch, x);
ARR3ADD(ch, y);

/* tube radius */
r = 0.5 * VEC3MAG(ch) / G_PI;

/* translation vector */
ARR3SET(x, a1);
ARR3SET(y, a2);
VEC3MUL(x, np);
VEC3MUL(y, mp);
ARR3SET(t, x);
ARR3SUB(t, y);

/* loop limits */
imin = MIN(MIN(np, 0), n);
imax = MAX(MAX(n+np, n), np);
jmin = MIN(MIN(-mp, 0), m);
jmax = MAX(MAX(m-np, m), -mp);

#if EDIT_NANOTUBE_NEW
printf("chirality vector index: (%d, %d)\n", n, m);
printf("chirality xlat index: (%d, %d)\n", np, -mp);
P3VEC("chiral vector: ", ch);
P3VEC("xlat vector: ", t);
#endif

/* loop over graphite lattice */
for (i=imin ; i<=imax ; i++)
  {
  for (j=jmin ; j<=jmax ; j++)
    {
/* compute lattice vector */
    ARR3SET(x, a1);
    ARR3SET(y, a2);
    VEC3MUL(x, i);
    VEC3MUL(y, j);
    ARR3SET(v, x);
    ARR3ADD(v, y);

/* compute basis atom coord */
    for (k=0 ; k<2 ; k++)
      {
      if (k)
        {
        ARR3ADD(v, a3);
        }
      ARR3SET(x, v);
      ARR3MUL(x, ch);
      p[0] = x[0] + x[1] + x[2];
      p[0] /= VEC3MAGSQ(ch);

      ARR3SET(y, v);
      ARR3MUL(y, t);
      p[1] = y[0] + y[1] + y[2];
      p[1] /= VEC3MAGSQ(t);

/* clamp to one unit's worth of atoms */
      fractional_clamp(p, dummy, 2);

/* choose basis atom type */
      if (k)
        core = core_new(edit_basis[0], NULL, model);
      else
        core = core_new(edit_basis[1], NULL, model);
      model->cores = g_slist_prepend(model->cores, core);

/* compute 1D fractional coords */
      a = 2.0 * G_PI * p[0];
/* xlat in x */
/*
      core->x[0] = -p[1];
*/
      core->x[0] = p[1];
      core->x[1] = r*sin(a);
      core->x[2] = r*cos(a);
      }
    }
  }

/* setup */
/* TODO - what if periodicity in z is desired??? */
model->periodic = 1;
model->fractional = TRUE;
model->pbc[0] = VEC3MAG(t);
model_prep(model);

/* model info */
text = g_strdup_printf("(%d, %d)", n, m);
property_add_ranked(2, "Chirality", text, model);
g_free(text);
text = g_strdup_printf("%.4f Angs", r);
property_add_ranked(3, "Radius", text, model);
g_free(text);

tree_model_refresh(model);

gui_refresh(GUI_MODEL_PROPERTIES);
gui_refresh(GUI_CANVAS);

gui_relation_update(model);
}

/***********************************************/
/* periodicity dialog hook to make a supercell */
/***********************************************/
void edit_make_supercell(void)
{
struct model_pak *model;

model = sysenv.active_model;
if (model)
  {
  space_make_supercell(model);
  model_prep(model);
  }

/* REFRESH */
gui_refresh(GUI_MODEL_PROPERTIES);
gui_refresh(GUI_CANVAS);
}

/************************/
/* make active model P1 */
/************************/
void edit_make_p1(void)
{
struct model_pak *model;

model = sysenv.active_model;
if (!model)
  return;
if (model->periodic != 3)
  return;

space_make_p1(model);

/* REFRESH */
gui_refresh(GUI_MODEL_PROPERTIES);
gui_refresh(GUI_CANVAS);
}

/***********************/
/* start add atom mode */
/***********************/
void edit_atom_add(void)
{
if (!sysenv.active_model)
  edit_model_create();
gui_mode_switch(ATOM_ADD);
}

/********************************/
/* add shells to selected cores */
/********************************/
void edit_shells_add(void)
{
GSList *list;
struct core_pak *core;
struct shel_pak *shell;
struct model_pak *model;

model = sysenv.active_model;
if (model)
  {
  for (list=model->selection ; list ; list=g_slist_next(list))
    {
    core = list->data;

/* create new shell if none present */
    if (!core->shell)
      {
      shell = shell_new(core->atom_label, NULL, model);
      model->shels = g_slist_prepend(model->shels, shell);

/* start shell at core coords */
      ARR3SET(shell->x, core->x);
      ARR3SET(shell->rx, core->rx);

/* transfer appropriate core characteristics */
      shell->primary = core->primary;
      shell->orig = core->orig;
      shell->region = core->region;

/* do core-shell link */
      shell->core = core;
      core->shell = shell;
      }
    }
  }
}

/*****************************/
/* atom/mol cell confinement */
/*****************************/
void edit_confine(GtkWidget *w, gint mode)
{
struct model_pak *model;

model = sysenv.active_model;
if (!model)
  return;

switch (mode)
  {
  case CORE:
    coords_confine_cores(model->cores, model);
    coords_compute(model);
    connect_bonds(model);
    break;
  case MOL:
    connect_molecules(model);
    break;
  }
redraw_canvas(SINGLE);
}

/******************/
/* bonding toggle */
/******************/
void gui_connect_toggle(void)
{
struct model_pak *model = sysenv.active_model;

if (model)
  {
  model->show_bonds ^= 1;

  redraw_canvas(SINGLE);
  }
}

/**********************/
/* model builder page */
/**********************/
/* TODO - molecule fragments (library?) */
void build_page(GtkWidget *box)
{
GtkWidget *vbox, *hbox;
GtkWidget *label, *entry, *spin;

g_return_if_fail(box != NULL);

/* Frame */
vbox = gui_frame_vbox("Atoms", FALSE, FALSE, box);

gui_button_x("Add atoms", edit_atom_add, NULL, vbox);
gui_button_x("Confine atoms to cell", edit_confine, (gpointer) CORE, vbox);
gui_button_x("Confine molecules to cell", edit_confine, (gpointer) MOL, vbox);
gui_button_x("Add shells to selected cores", edit_shells_add, NULL, vbox);

/* Frame */
/* FIXME - fix up user bonds then re-introduce */
vbox = gui_frame_vbox("Bonds", FALSE, FALSE, box);

gui_button_x("Add single bonds ", gtk_mode_switch, (gpointer) BOND_SINGLE, vbox);
/*
gui_button_x("Add double bonds ", gtk_mode_switch, (gpointer) BOND_DOUBLE, vbox);
gui_button_x("Add triple bonds ", gtk_mode_switch, (gpointer) BOND_TRIPLE, vbox);
*/
gui_button_x("Delete bonds ", gtk_mode_switch, (gpointer) BOND_DELETE, vbox);
gui_button_x("Toggle bonding ", gui_connect_toggle, NULL, vbox);

/* Frame */
vbox = gui_frame_vbox("Structural", FALSE, FALSE, box);

gui_button_x("Make new model", edit_model_create, NULL, vbox);
gui_button_x("Make supercell", edit_make_supercell, NULL, vbox);
gui_button_x("Force structure to P1", edit_make_p1, NULL, vbox);

/* nanotube setup */
vbox = gui_frame_vbox("Nanotube", FALSE, FALSE, box);

/* chirality */
hbox = gtk_hbox_new(FALSE,0);
gtk_box_pack_start(GTK_BOX(vbox), hbox, FALSE, FALSE, 0);

label = gtk_label_new("Chirality  ");
gtk_box_pack_start(GTK_BOX(hbox),label,FALSE,FALSE,0); 
spin = gui_direct_spin(NULL, &edit_chirality[0], 0, 99, 1, NULL, NULL, NULL);
gtk_box_pack_start(GTK_BOX(hbox),spin,FALSE,FALSE,0); 
spin = gui_direct_spin(NULL, &edit_chirality[1], 0, 99, 1, NULL, NULL, NULL);
gtk_box_pack_start(GTK_BOX(hbox),spin,FALSE,FALSE,0); 

/* basis atoms */
/*
hbox = gtk_hbox_new(FALSE,0);
gtk_box_pack_start(GTK_BOX(vbox), hbox, FALSE, FALSE, 0);
*/
if (!edit_basis[0])
  edit_basis[0] = g_strdup("C");
if (!edit_basis[1])
  edit_basis[1] = g_strdup("C");
entry = gui_text_entry("  Basis  ", &edit_basis[0], TRUE, TRUE, hbox);
gtk_entry_set_width_chars(GTK_ENTRY(entry), 4);
entry = gui_text_entry(" - ", &edit_basis[1], TRUE, TRUE, hbox);
gtk_entry_set_width_chars(GTK_ENTRY(entry), 4);

/* characteristic length */
gui_direct_spin(" : ", &edit_length, 0.1, 5.0, 0.05, NULL, NULL, hbox);

gui_button_x("Create nanotube ", edit_nanotube_new, NULL, vbox);

gtk_widget_show_all(box);
}

/*******************/
/* spatial globals */
/*******************/
GtkListStore *spatial_list=NULL;
GtkWidget *spatial_tree=NULL;
gpointer spatial_selected=NULL;

/**************************************/
/* delete all spatials of given label */
/**************************************/
void gui_spatial_delete(GtkWidget *w, gpointer data)
{
const gchar *label = data;

if (sysenv.active_model)
  spatial_destroy_by_label(label, sysenv.active_model);

sysenv.refresh_dialog=TRUE;

redraw_canvas(SINGLE);
}

/***********************/
/* delete all spatials */
/***********************/
void gui_spatial_delete_all(GtkWidget *w, gpointer data)
{
struct model_pak *model = sysenv.active_model;

if (model)
  {
  spatial_destroy_all(model);
  sysenv.refresh_dialog=TRUE;
  redraw_canvas(SINGLE);
  }
}

/*****************************************/
/* delete the currently selected spatial */
/*****************************************/
void gui_spatial_delete_selected(GtkWidget *w, gpointer data)
{
spatial_destroy(spatial_selected, sysenv.active_model);
sysenv.refresh_dialog=TRUE;
redraw_canvas(SINGLE);
}

/*****************************/
/* fill out the spatial list */
/*****************************/
void gui_spatial_populate(void)
{
gint n;
gchar *size, *type;
GSList *list;
GtkTreeIter iter;
struct spatial_pak *spatial;
struct model_pak *model;

gtk_list_store_clear(spatial_list);

model = sysenv.active_model;
if (!model)
  return;

for (list=model->spatial ; list ; list=g_slist_next(list))
  {
  spatial = list->data;

/* compute number of primitives */
  n = g_slist_length(spatial->list);
  if (spatial->size)
    n /= spatial->size;
  size = g_strdup_printf("%d", n);

/* primitive type */
  switch (spatial->size)
    {
    case 1:
      type = g_strdup("points");
      break;
    case 2:
      type = g_strdup("vectors");
      break;
    case 3:
      type = g_strdup("triangles");
      break;
    case 4:
      type = g_strdup("quads");
      break;
    default:
      type = g_strdup("vertices");
      break;
    }

/* add the spatial and descriptors */
  gtk_list_store_append(spatial_list, &iter);
  gtk_list_store_set(spatial_list, &iter, 0, spatial->label, 
                                          1, size,
                                          2, type,
                                          3, spatial,
                                         -1);

/* NB: the list will make its own copy */
  g_free(size);
  g_free(type);
  }
}

/*******************************/
/* select a particular spatial */
/*******************************/
void gui_spatial_select(GtkTreeSelection *selection, gpointer data)
{
GtkTreeIter iter;
GtkTreeModel *treemodel;
gpointer spatial;
struct model_pak *model;

/* checks */
model = sysenv.active_model;
if (!model)
  return;

/* record selection as active */
if (gtk_tree_selection_get_selected(selection, &treemodel, &iter))
  {
  gtk_tree_model_get(treemodel, &iter, 3, &spatial, -1);
  spatial_selected = spatial;
  }
}

/**************************************/
/* callback for spatial colour change */
/**************************************/
void gui_spatial_colour_all(GtkWidget *w, gdouble *colour)
{
GSList *list1, *list2;
struct model_pak *model;
struct spatial_pak *spatial;
struct vec_pak *vertex;

model = sysenv.active_model;

if (model)
  {
  for (list1=model->spatial ; list1 ; list1=g_slist_next(list1))
    {
    spatial = list1->data;

    for (list2=spatial->list ; list2 ; list2=g_slist_next(list2))
      {
      vertex = list2->data;
      ARR3SET(vertex->colour, colour);
      } 
    }
  }
gui_refresh(GUI_CANVAS);
}

/***********************************************/
/* callback for selected spatial colour change */
/***********************************************/
void gui_spatial_colour_select(GtkWidget *w, gdouble *colour)
{
GSList *list;
struct vec_pak *vertex;
struct spatial_pak *spatial = spatial_selected;

if (spatial)
  {
  for (list=spatial->list ; list ; list=g_slist_next(list))
    {
    vertex = list->data;
    ARR3SET(vertex->colour, colour);
    }
  gui_refresh(GUI_CANVAS);
  }
}

/********************************************/
/* callback for spatial label colour change */
/********************************************/
void gui_spatial_colour_label_all(GtkWidget *w, gdouble *colour)
{
GSList *list;
struct model_pak *model;
struct spatial_pak *spatial;

model = sysenv.active_model;

if (model)
  {
  for (list=model->spatial ; list ; list=g_slist_next(list))
    {
    spatial = list->data;
    ARR3SET(spatial->c, colour);
    }
  }
gui_refresh(GUI_CANVAS);
}

/*****************************************************/
/* callback for selected spatial label colour change */
/*****************************************************/
void gui_spatial_colour_label_select(GtkWidget *w, gdouble *colour)
{
struct spatial_pak *spatial = spatial_selected;

if (spatial)
  {
  ARR3SET(spatial->c, colour);
  gui_refresh(GUI_CANVAS);
  }
}

/************************/
/* spatial objects page */
/************************/
void spatial_page(GtkWidget *box)
{
gint i;
gchar *titles[] = {" Label ", " Size ", " Primitive "};
GtkCellRenderer *r;
GtkTreeViewColumn *c;
GtkTreeSelection *select;
GtkWidget *swin, *hbox, *vbox1, *vbox2, *vbox;

/* initialize - nothing selected */
spatial_selected = NULL;

/* left & right pane split */
hbox = gtk_hbox_new(FALSE, PANEL_SPACING);
gtk_container_add(GTK_CONTAINER(box), hbox);
vbox1 = gtk_vbox_new(FALSE, PANEL_SPACING);
gtk_box_pack_start(GTK_BOX(hbox), vbox1, FALSE, FALSE, 0);
vbox2 = gtk_vbox_new(FALSE, PANEL_SPACING);
gtk_box_pack_start(GTK_BOX(hbox), vbox2, TRUE, TRUE, 0);

/* Frame */
vbox = gui_frame_vbox("Adding", FALSE, FALSE, vbox1);

/* mode buttons */
gui_button_x("Add vectors ", 
               gtk_mode_switch, GINT_TO_POINTER(DEFINE_VECTOR),
               vbox);
gui_button_x("Add planes ",  
               gtk_mode_switch, GINT_TO_POINTER(DEFINE_PLANE),
               vbox);

/* FIXME - currently, ribbons are not spatial objects */
gui_button_x("Add ribbons ", 
               gtk_mode_switch, GINT_TO_POINTER(DEFINE_RIBBON),
               vbox);

/* Frame */
vbox = gui_frame_vbox("Deleting", FALSE, FALSE, vbox1);

gui_button_x("Delete all vectors ", gui_spatial_delete, "vector", vbox);
gui_button_x("Delete all planes ", gui_spatial_delete, "plane", vbox);
gui_button_x("Delete all ", gui_spatial_delete_all, NULL, vbox);
gui_button_x("Delete selected ", gui_spatial_delete_selected, NULL, vbox);

/* Frame */
vbox = gui_frame_vbox(NULL, FALSE, FALSE, vbox1);

/* TODO - apply to selection */
gui_colour_box("Spatial fill colour  ", edit_spatial_colour, vbox);
gui_button_x("Apply to all ",
               gui_spatial_colour_all, edit_spatial_colour, vbox);
gui_button_x("Apply to selected ",
               gui_spatial_colour_select, edit_spatial_colour, vbox);

/* Frame */
vbox = gui_frame_vbox(NULL, FALSE, FALSE, vbox1);

/* TODO - apply to selection */
gui_colour_box("Spatial label colour ", edit_label_colour, vbox);
gui_button_x("Apply to all ",
               gui_spatial_colour_label_all, edit_label_colour, vbox);
gui_button_x("Apply to selected ",
               gui_spatial_colour_label_select, edit_label_colour, vbox);

/* spatial list */
vbox = gui_frame_vbox(NULL, TRUE, TRUE, vbox2);

/* scrolled window */
swin = gtk_scrolled_window_new(NULL, NULL);
gtk_scrolled_window_set_policy(GTK_SCROLLED_WINDOW(swin),
                               GTK_POLICY_AUTOMATIC, GTK_POLICY_AUTOMATIC);
gtk_box_pack_start(GTK_BOX(vbox), swin, TRUE, TRUE, 0);

/* list */
spatial_list = gtk_list_store_new(4, G_TYPE_STRING, G_TYPE_STRING, G_TYPE_STRING, G_TYPE_POINTER);
spatial_tree = gtk_tree_view_new_with_model(GTK_TREE_MODEL(spatial_list));
gtk_scrolled_window_add_with_viewport(GTK_SCROLLED_WINDOW(swin), spatial_tree);
for (i=0 ; i<3 ; i++)
  {
  r = gtk_cell_renderer_text_new();
  c = gtk_tree_view_column_new_with_attributes(titles[i], r, "text", i, NULL);
  gtk_tree_view_append_column(GTK_TREE_VIEW(spatial_tree), c);
  gtk_tree_view_set_headers_visible(GTK_TREE_VIEW(spatial_tree), FALSE);
  }
gtk_tree_view_set_headers_visible(GTK_TREE_VIEW(spatial_tree), TRUE);

/* selection handler */
select = gtk_tree_view_get_selection(GTK_TREE_VIEW(spatial_tree));
gtk_tree_selection_set_mode(select, GTK_SELECTION_SINGLE);
g_signal_connect(G_OBJECT(select), "changed",
                 G_CALLBACK(gui_spatial_select),
                 NULL);

gui_spatial_populate();

gtk_widget_show_all(box);
}

/************************/
/* transformations page */
/************************/
void trans_page(GtkWidget *box)
{
gint i, j;
GList *list;
GtkWidget *table, *hbox, *vbox, *label;

g_assert(box != NULL);

/* transformation construction */
vbox = gui_frame_vbox("Construction", FALSE, FALSE, box);

/* TODO - use spinners instead of gtk entries for these numbers */
/* axes order */
hbox = gtk_hbox_new(FALSE, PANEL_SPACING);
gtk_box_pack_start(GTK_BOX(vbox), hbox, TRUE, FALSE, 0);

label = gtk_label_new("Rotation angle (degrees):");
gtk_box_pack_start(GTK_BOX(hbox), label, FALSE, FALSE, 0);
axis_entry = gtk_entry_new();
gtk_box_pack_end(GTK_BOX(hbox), axis_entry, FALSE, FALSE, 0);
/*
gtk_entry_set_text(GTK_ENTRY(axis_entry), "90");
*/

/* spatial object */
hbox = gtk_hbox_new(FALSE, PANEL_SPACING);
gtk_box_pack_start(GTK_BOX(vbox), hbox, TRUE, FALSE, 0);

label = gtk_label_new("Reference spatial object: ");
gtk_box_pack_start(GTK_BOX(hbox), label, FALSE, FALSE, 0);
obj_entry = gtk_entry_new();
gtk_box_pack_end(GTK_BOX(hbox), obj_entry, FALSE, FALSE, 0);
/*
gtk_entry_set_text(GTK_ENTRY(obj_entry), "0");
*/

/* construction buttons */
hbox = gtk_hbox_new(FALSE, PANEL_SPACING);
gtk_box_pack_start(GTK_BOX(vbox), hbox, FALSE, FALSE, 0);

list = NULL;
list = g_list_prepend(list, "identity matrix");
list = g_list_prepend(list, "lattice matrix");
list = g_list_prepend(list, "reflection matrix");
list = g_list_prepend(list, "rotation matrix");
list = g_list_prepend(list, "z alignment matrix");
list = g_list_reverse(list);

edit_construct = gui_pulldown_new("Construct", list, FALSE, hbox);

gui_button_x(NULL, construct_transmat, NULL, hbox);


/* table */
table = gtk_table_new(5, 4, FALSE);
gtk_container_add(GTK_CONTAINER(GTK_BOX(vbox)),table);

label = gtk_label_new("a*");
gtk_table_attach_defaults(GTK_TABLE(table),label,0,1,1,2);
label = gtk_label_new("b*");
gtk_table_attach_defaults(GTK_TABLE(table),label,0,1,2,3);
label = gtk_label_new("c*");
gtk_table_attach_defaults(GTK_TABLE(table),label,0,1,3,4);

label = gtk_label_new("a");
gtk_table_attach_defaults(GTK_TABLE(table),label,1,2,0,1);
label = gtk_label_new("b");
gtk_table_attach_defaults(GTK_TABLE(table),label,2,3,0,1);
label = gtk_label_new("c");
gtk_table_attach_defaults(GTK_TABLE(table),label,3,4,0,1);
label = gtk_label_new("t");
gtk_table_attach_defaults(GTK_TABLE(table),label,4,5,0,1);

/* column loop */
for (j=0 ; j<3 ; j++)
  {
/* row loop */
  for (i=0 ; i<3 ; i++)
    {
    transmat[3*j+i] = gtk_entry_new();
    gtk_table_attach_defaults(GTK_TABLE(table),transmat[3*j+i],i+1,i+2,j+1,j+2);
    gtk_widget_set_size_request(transmat[3*j+i], 7*sysenv.gtk_fontsize, -1);
    }
  }
/* translation */
for (j=0 ; j<3 ; j++)
  {
  transmat[9+j] = gtk_entry_new();
  gtk_table_attach_defaults(GTK_TABLE(table),transmat[9+j],4,5,j+1,j+2);
  gtk_widget_set_size_request(transmat[9+j], 7*sysenv.gtk_fontsize, -1);
  }

/* control buttons */
vbox = gui_frame_vbox("Coordinate transformations", FALSE, FALSE, box);

gui_button_x("Apply to all atoms", 
               apply_transmat, GINT_TO_POINTER(AT_ANY),
               vbox);
gui_button_x("Apply to selected atoms",
               apply_transmat, GINT_TO_POINTER(AT_SELECTED),
               vbox);

/* animation construction */
/* FIXME - broken by new quaternion camera */
/*
hbox = gtk_hbox_new(FALSE, 0);
gtk_box_pack_start(GTK_BOX(vbox), hbox, FALSE, FALSE, 0);
label = gtk_label_new("Generate animation frames ");
gtk_box_pack_start(GTK_BOX(hbox), label, FALSE, FALSE, 0);
gui_direct_spin(NULL, &edit_anim_n, 0.0, 1000.0, 1.0, NULL, NULL, hbox);
gui_button_x(NULL, apply_n_transmat, GINT_TO_POINTER(AT_ANY), hbox);
*/

/* control buttons */
vbox = gui_frame_vbox("Lattice transformations", FALSE, FALSE, box);

gui_button_x("Apply to lattice matrix", edit_transform_latmat, NULL, vbox);

hbox = gtk_hbox_new(FALSE, PANEL_SPACING);
gtk_box_pack_start(GTK_BOX(vbox), hbox, FALSE, FALSE, 0);
label = gtk_label_new("Alter lattice periodicity ");
gtk_box_pack_start(GTK_BOX(hbox), label, FALSE, FALSE, 0);
periodicity_spin = gtk_spin_button_new_with_range(0, 3, 1);
gtk_box_pack_start(GTK_BOX(hbox), periodicity_spin, FALSE, FALSE, 0);
gui_button_x(NULL, cb_modify_periodicity, NULL, hbox);

gui_button_x("Create new lattice model from linear combination",
               apply_transmat, GINT_TO_POINTER(AT_LATTICE),
               vbox);
}

/*************************************/
/* save a diffraction (DIFFAX) setup */
/*************************************/
struct model_pak *diffract_model=NULL;

void diffract_save(gchar *name)
{
write_diffax(name, diffract_model);

dialog_destroy_type(FILE_SELECT);
}

/******************************************/
/* create a file dialog for a DIFFAX save */
/******************************************/
void diffract_save_dialog(GtkWidget *w, struct model_pak *model)
{
diffract_model = model;
file_dialog("Save DIFFAX file", model->basename, FILE_SAVE,
            (gpointer) diffract_save, DIFFAX_INP);
}

/*******************/
/* add a new layer */
/*******************/
GtkWidget *diffract_layer_total;

#define DEBUG_DIFFRACT_LAYER_SETUP 0
void diffract_layer_setup(struct model_pak *data)
{
gint n, region;
gint num_layer, tot_layer;
gdouble start, stop;
GSList *clist, *list1=NULL;
struct layer_pak *layer;
struct core_pak *core;

if (!data)
  return;

if (GTK_IS_SPIN_BUTTON(diffract_layer_total))
  tot_layer = SPIN_IVAL(GTK_SPIN_BUTTON(diffract_layer_total));
else
  return;

/* get rid of old list */
free_slist(data->layer_list);
data->layer_list = NULL;

for (num_layer=0 ; num_layer<tot_layer ; num_layer++)
  {

/* create the new layer */
layer = g_malloc(sizeof(struct layer_pak));
layer->width = 1.0/(gdouble) tot_layer;
VEC3SET(layer->centroid, 0.0, 0.0, 0.0);

start = (gdouble) num_layer / (gdouble) tot_layer;
stop = (gdouble) (num_layer+1) / (gdouble) tot_layer;

#if DEBUG_DIFFRACT_LAYER_SETUP
printf("new layer: %f-%f (%f)\n", start, stop, layer->width);
#endif

region = num_layer;

/* add cores that satisfy the start/stop condition */
list1 = NULL;
for (clist=data->cores ; clist ; clist=g_slist_next(clist))
  {
  core = (struct core_pak *) clist->data;

/* FIXME - floating point boundary problems? */
  if (core->x[2] >= start && core->x[2] < stop)
    {
    list1 = g_slist_prepend(list1, core); 
    ARR3ADD(layer->centroid, core->x);
    core->region = region;
/* update region colouring (if necessary) */
    if (data->colour_scheme == REGION)
      atom_colour_scheme(REGION, core, data);
    }
  }

/* centroid calc. */
n = g_slist_length(list1);
if (n)
  {
#if DEBUG_DIFFRACT_LAYER_SETUP
printf("[%d] added layer with %d cores.\n", num_layer, n);
#endif

  layer->cores = list1;
  VEC3MUL(layer->centroid, 1.0/(gdouble) n);
  }

/* always save the layer - may *want* to have an empty layer */
/* eg to simulate a vacuum gap */
data->layer_list = g_slist_prepend(data->layer_list, layer);
  }

#if DEBUG_DIFFRACT_LAYER_SETUP
printf("Total layers: %d\n", g_slist_length(data->layer_list));
#endif

redraw_canvas(SINGLE);
}

/***********************************************/
/* create a defect structure from input layers */
/***********************************************/
GtkWidget *diffract_layer_order;

#define DEBUG_DIFFRACT_MODEL_CREATE 0
void diffract_model_create(GtkWidget *w, struct model_pak *model)
{
gint i, n, num_layer, tot_layer;
gdouble new_width, old_width, offset;
const gchar *text;
GSList *list;
struct model_pak *dest_model;
struct layer_pak *layer;
struct core_pak *core;

if (!model)
  return;

diffract_layer_setup(model);
text = gtk_entry_get_text(GTK_ENTRY(diffract_layer_order));

/* create a new model */
dest_model = model_new();
g_return_if_fail(dest_model != NULL);
/* NB: not DIFFAX_INP, since this is demonstrating a particular layer packing */
dest_model->id = CREATOR;

strcpy(dest_model->filename, "new_model");
g_free(dest_model->basename);
dest_model->basename = g_strdup("new model");

/* find how many valid layers we will create */
tot_layer = 0;
for (i=0 ; i<strlen(text) ; i++)
  {
  if (g_ascii_isdigit(text[i]))
    {
    n = text[i] - '0';
    layer = g_slist_nth_data(model->layer_list, n-1);
    if (layer)
      tot_layer++;
    }
  }

#if DEBUG_DIFFRACT_MODEL_CREATE
printf("Requested layers: %d\n", tot_layer);
#endif

/* build the new model from the input layer string */
num_layer=0;
new_width=1.0/tot_layer;
old_width=1.0;

for (i=strlen(text) ; i-- ; )
  {
  if (g_ascii_isdigit(text[i]))
    {
    n = text[i] - '0';
/* NB: numbering starts from 0 */
    layer = g_slist_nth_data(model->layer_list, n-1);
    if (layer)
      {
      old_width = layer->width;

/* add j based z offset */
      offset = num_layer+0.5;
      offset *= new_width;

#if DEBUG_DIFFRACT_MODEL_CREATE
printf("[%d] inserting layer: %d (scale: %f) (offset: %f)\n",
        num_layer, n, new_width/old_width, offset);
#endif

/* add the layer's cores */
      for (list=layer->cores ; list ; list=g_slist_next(list))
        {
        core = dup_core(list->data);

        core->region = tot_layer - i - 1;

/* remove z centroid */
        core->x[2] -= layer->centroid[2];

/* z values need to be SCALED (to meet the new layer width) */
        core->x[2] *= new_width/old_width;

        core->x[2] += offset;

        dest_model->cores = g_slist_prepend(dest_model->cores, core);
        }
      num_layer++;
      }
    else
      printf("Layer number %d not defined.\n", n);
    }
  }

/* test if anything was created */
if (!num_layer)
  {
  model_delete(dest_model);
  return;
  }

dest_model->cores = g_slist_reverse(dest_model->cores);

/* setup */
dest_model->periodic = 3;
dest_model->fractional = TRUE;
ARR3SET(&dest_model->pbc[0], &model->pbc[0]);
ARR3SET(&dest_model->pbc[3], &model->pbc[3]);
dest_model->pbc[2] *= num_layer*old_width;

#if DEBUG_DIFFRACT_MODEL_CREATE
printf("Setting c to: %f\n", dest_model->pbc[2]);
#endif

model_prep(dest_model);

model_colour_scheme(model->colour_scheme, dest_model);

tree_model_add(dest_model);
}

/**********************************/
/* callback for forcefield typing */
/**********************************/
void cb_type_model(GtkWidget *w, gpointer data)
{
struct model_pak *model;

model = sysenv.active_model;
if (!model)
  return;

type_model(gtk_entry_get_text(GTK_ENTRY(data)), model);

/* possible label update */
redraw_canvas(ALL);
}

/**************************/
/* cell sculpting globals */
/**************************/
gdouble sculpt_length = 20.0;

/***********************************************************/
/* tests if a given image t[] is inside the list of planes */
/***********************************************************/
gint sculpt_image_test(gint *t, gdouble scale, struct model_pak *model)
{
gdouble r[3], n[3];
GSList *list;
struct plane_pak *plane;

/* convert periodic image to cartesian point */
/* NB: ensure we get the closest lattice point to the origin */
if (t[0] < 0)
  r[0] = t[0]+1;
else
  r[0] = t[0];

if (t[1] < 0)
  r[1] = t[1]+1;
else
  r[1] = t[1];

if (t[2] < 0)
  r[2] = t[2]+1;
else
  r[2] = t[2];

vecmat(model->latmat, r);

/* compare against planes to see if this image is excluded */
for (list=model->planes ; list ; list=g_slist_next(list))
  {
  plane = list->data;

/* get cartesian normal */
  ARR3SET(n, plane->index);
  vecmat(model->rlatmat, n);
  normalize(n, 3);

/* test dot product against required distance to plane */
  ARR3MUL(n, r);
  if ((n[0]+n[1]+n[2]) > scale*plane->f[0])
    {
    return(FALSE);
    }
  }

return(TRUE);
}

/***************************/
/* cell sculpting callback */
/***************************/
#define DEBUG_SCULPT_CREATE 0
void sculpt_model_create(struct model_pak *model)
{
gint i, t[3], limit[3];
gdouble scale, r=1.0, s, rmin, x[3], n[3];
GSList *list, *clist, *slist, *plist;
struct model_pak *dest;
struct core_pak *core;
struct shel_pak *shell;
struct mol_pak *mol;
struct plane_pak *plane;
struct spatial_pak *spatial;

/* checks */
g_assert(model != NULL);
if (model->periodic != 3)
  {
  gui_text_show(ERROR, "Source model is not 3D periodic.\n");
  return;
  }
if (!model->planes)
  {
  gui_text_show(ERROR, "No cleavage planes supplied.\n");
  return;
  }
#if DEBUG_SCULPT_CREATE
else
  printf("Planes: %d\n", g_slist_length(model->planes));
#endif

/* compute required periodic images (assumed symmetric) */
/* FIXME - most of the nuclei shape problems come from an insufficient number of repeats */
VEC3SET(limit, 0, 0, 0);
for (i=0 ; i<model->periodic ; i++)
  limit[i] = 1 + sculpt_length/model->pbc[i];

/* init destination model for sculpture */
dest = model_new();
if (!dest)
  {
  gui_text_show(ERROR, "Failed to allocate for new model.\n");
  return;
  }
model_init(dest);
gulp_data_copy(model, dest);

#if DEBUG_SCULPT_CREATE 
printf("maximum length: %f\n", sculpt_length);
printf("periodic images: %d %d %d\n", limit[0], limit[1], limit[2]);
#endif

/* TODO - going to have to put this ABOVE the sculpt_image_test() call, so we get rmin */
/* morphology computation */
morph_build(model);
/* compute the facet closest to the center */
rmin = G_MAXDOUBLE;
for (plist = model->planes ; plist ; plist=g_slist_next(plist))
  {
  plane = plist->data;

/* ensure we've got the best shift value data sitting in the plane structure */
  update_plane_energy(plane, model);

/* compute distance to plane facet */
/* NEW - a little naughty, but using the f[] (structure factor) to store the length/shift */
/* value of the plane in the nuclei sculpting so we dont keep recalculating */
  switch (model->morph_type)
    {
    case EQUIL_UN:
      plane->f[0] = plane->esurf[0];
      plane->f[1] = plane->esurf_shift;
      break;
    case GROWTH_UN:
      plane->f[0] = fabs(plane->eatt[0]);
      plane->f[1] = plane->eatt_shift;
      break;
    case EQUIL_RE:
      plane->f[0] = plane->esurf[1];
      plane->f[1] = plane->esurf_shift;
      break;
    case GROWTH_RE:
      plane->f[0] = fabs(plane->eatt[1]);
      plane->f[1] = plane->eatt_shift;
      break;
    case DHKL:
    default:
      plane->f[0] = 1.0/plane->dhkl;
      plane->f[1] = 0.0;
      break;
    }
/* NEW - stop zero result (failed/bad calc) from messing things up */
  if (plane->f[0] != 0.0 && plane->f[0] < rmin)
    rmin = plane->f[0];
  }

#if DEBUG_SCULPT_CREATE 
printf("rmin = %f\n", rmin);
#endif

scale = 0.5*sculpt_length/rmin;

/* periodic image creation */
for (t[0] = -limit[0] ; t[0] <= limit[0] ; t[0]++)
  {
  for (t[1] = -limit[1] ; t[1] <= limit[1] ; t[1]++)
    {
    for (t[2] = -limit[2] ; t[2] <= limit[2] ; t[2]++)
      {
/* test vertices of the t[] image against planes */
      if (!sculpt_image_test(t, scale, model))
        continue;

/* duplicate cores, add periodic image offset & convert to cartesian */
      clist = dup_core_list(model->cores);
      slist = dup_shell_list(model->shels);
      dest->cores = g_slist_concat(dest->cores, clist);
      dest->shels = g_slist_concat(dest->shels, slist);
      for (list=clist ; list ; list=g_slist_next(list))
        {
        core = list->data;
        core->primary = TRUE;
        ARR3ADD(core->x, t);
        vecmat(model->latmat, core->x);
        }
      for (list=slist ; list ; list=g_slist_next(list))
        {
        shell = list->data;
        shell->primary = TRUE;
        ARR3ADD(shell->x, t);
        vecmat(model->latmat, shell->x);
        }
      }
    }
  }

/* initialize connectivity and core-shell links for the new model */
zone_init(dest);
connect_bonds(dest);
connect_molecules(dest);
shell_make_links(dest);

/* plane cutoff tests */
for (plist=model->planes ; plist ; plist=g_slist_next(plist))
  {
  plane = plist->data;

/* get cartesian distance to plane facet */
  r = plane->f[0] * scale;
  s = plane->f[1];

/* attempt to create the correct (ie surface shift) termination */
if (model->sculpt_shift_use)
  {
/*
  g = GCD(GCD(plane->index[0], plane->index[1]), GCD(plane->index[1], plane->index[2]));
*/
  r /= plane->dhkl;

/* TODO - if nearest_int(r) == 0 -> set to 1 */
  s += nearest_int(r);
  r = s * plane->dhkl;
  }

/* get cartesian normal */
  ARR3SET(n, plane->index);
  vecmat(model->rlatmat, n);
  normalize(n, 3);

/* project coords onto plane normal & remove if greater than distance to facet */
/* TODO - do cutoff by molecule centroid when combining with periodic image loop */
  for (list=dest->moles ; list ; list=g_slist_next(list))
    {
    mol = list->data;

    ARR3SET(x, mol->centroid);
    ARR3MUL(x, n);
    if ((x[0]+x[1]+x[2]) > r)
      {
      for (clist=mol->cores ; clist ; clist=g_slist_next(clist))
        {
        core = clist->data;
        core->status |= DELETED;
        if (core->shell)
          {
          shell = core->shell;
          shell->status |= DELETED;
          }
        }
      }
    }
  }

/* transfer the spatials from the source model morphology to the nuclei */
list = model->spatial;
while (list)
  {
  spatial = list->data;
  list = g_slist_next(list);

/* search for all morphology related spatials */
/* a bit crude... */
  if (g_strrstr(spatial->label, "(") && g_strrstr(spatial->label, ")"))
    {
    dest->spatial = g_slist_prepend(dest->spatial, spatial);
    model->spatial = g_slist_remove(model->spatial, spatial);
    spatial->method = GL_LINE_LOOP;

/* mult vertices by len */
    for (plist=spatial->list ; plist ; plist=g_slist_next(plist))
      {
      struct vec_pak *vec = plist->data;

      ARR3SET(vec->colour, sysenv.render.fg_colour);
/* convert fractional vertices to cartesian */
      vecmat(model->latmat, vec->x);
/* scale the morphology to match the constructed nuclei */
      VEC3MUL(vec->x, scale);
      }
    }
  }

/* init for display */
delete_commit(dest);
model_prep(dest);
tree_model_add(dest);
redraw_canvas(ALL);
}

/********************************/
/* callback for nuclei creation */
/********************************/
void morph_sculpt(GtkWidget *w, gpointer data)
{
gpointer model;

sculpt_length = SPIN_FVAL(GTK_SPIN_BUTTON(data));

model = g_object_get_data(data, "model");

sculpt_model_create(model);
}

/*******************************************/
/* unit cell sculpting (eg via morphology) */
/*******************************************/
void sculpting_page(GtkWidget *box)
{
GtkWidget *hbox, *hbox2;

hbox = gtk_hbox_new(FALSE, 0);
gtk_box_pack_start(GTK_BOX(box), hbox, FALSE, FALSE, 0);
gui_direct_spin("Maximum length", &sculpt_length, 10.0, 100.0, 10.0, NULL, NULL, hbox);

/* action buttons */
hbox2 = gtk_hbox_new(FALSE, 0);
gtk_box_pack_start(GTK_BOX(box), hbox2, FALSE, FALSE, PANEL_SPACING);
hbox = gtk_hbox_new(FALSE, PANEL_SPACING);
gtk_box_pack_end(GTK_BOX(hbox2), hbox, TRUE, FALSE, 0);

gui_icon_button(GTK_STOCK_APPLY, "Create ",
                  sculpt_model_create, NULL,
                  hbox);
}

/************************************/
/* build a new rule from the dialog */
/************************************/
void gui_make_rule(GtkWidget *w, gpointer dialog)
{
gint count, level;
const gchar *ff_type, *ff_elem;
gpointer type;
GtkWidget *obj;
struct model_pak *model;

model = sysenv.active_model;
if (!model)
  return;

obj = dialog_child_get(dialog, "FF_label");
ff_type = gtk_entry_get_text(GTK_ENTRY(obj));

obj = dialog_child_get(dialog, "FF_level");
level = str_to_float(gtk_entry_get_text(GTK_ENTRY(obj)));

obj = dialog_child_get(dialog, "FF_element");
ff_elem = gtk_entry_get_text(GTK_ENTRY(obj));

obj = dialog_child_get(dialog, "FF_count");
count = str_to_float(gtk_entry_get_text(GTK_ENTRY(obj)));

type = type_new();
type_ff_set(TRUE, ff_type, type);
type_rule_add(level, count, ff_elem, type);

type_apply(type, model->cores);

type_free(type);
}

/**********************************/
/* read in a set of FF parameters */
/**********************************/
void gui_import_ff(GtkWidget *w, gpointer dialog)
{
const gchar *name;
GSList *list;
GtkWidget *obj;
struct model_pak *model;

model = sysenv.active_model;
if (!model)
  return;

obj = dialog_child_get(dialog, "FF_import");
name = gtk_entry_get_text(GTK_ENTRY(obj));

list = gromacs_read_ff(name);

if (list)
  {
  gui_text_show(NORMAL, "Imported Forcefield.\n");

/* remove any previous GULP FF set */
  if (model->gulp.potentials)
    {
    g_free(model->gulp.potentials);
    model->gulp.potentials = NULL;
    }
/* remove any previous internal FF set */
  if (model->ff_list)
    free_slist(model->ff_list);

  model->ff_list = list;
  }
else
  gui_text_show(WARNING, "Empty or missing Forcefield.\n");
}

/*******************************/
/* core labelling manipulation */
/*******************************/
void labelling_page(GtkWidget *box, gpointer dialog)
{
GtkWidget *w, *vbox, *hbox;

/* Surface options */
vbox = gui_frame_vbox("Regions", FALSE, FALSE, box);

gui_button_x("Move selection up", region_move, (gpointer) UP, vbox);
gui_button_x("Move selection down", region_move, (gpointer) DOWN, vbox);

#define FF_TYPING 1
#if FF_TYPING
{
GList *list;
GtkWidget *combo, *hbox, *label;

/* EXP - forcefield assignment */
vbox = gui_frame_vbox("Typing", FALSE, FALSE, box);

/* typing setup */
hbox = gtk_hbox_new(FALSE, PANEL_SPACING);
gtk_box_pack_start(GTK_BOX(vbox), hbox, FALSE, FALSE, 0);

label = gtk_label_new("Assign atom: ");
gtk_box_pack_start(GTK_BOX(hbox), label, FALSE, FALSE, 0);

list = NULL;
list = g_list_prepend(list, "QEq charges");
list = g_list_prepend(list, "Gasteiger charges");
list = g_list_prepend(list, "Dreiding labels");

/* experimental */
/*
list = g_list_prepend(list, "CVFF labels");
*/

combo = gtk_combo_new();
gtk_entry_set_editable(GTK_ENTRY(GTK_COMBO(combo)->entry), FALSE);
gtk_combo_set_popdown_strings(GTK_COMBO(combo), list);
gtk_box_pack_start(GTK_BOX(hbox), combo, FALSE, FALSE, PANEL_SPACING);

gui_button_x(NULL, cb_type_model, GTK_COMBO(combo)->entry, hbox);
}
#endif

#if DIFFAX
/* DIFFAX stuff - shunted here */
/* frame */
vbox = gui_frame_vbox("DIFFAX", FALSE, FALSE, box);

/* source */
hbox = gtk_hbox_new(FALSE, 0);
gtk_box_pack_start(GTK_BOX(vbox), hbox, FALSE, FALSE, 0);
label = gtk_label_new(g_strdup_printf("Source model: %s", data->basename));
gtk_box_pack_start(GTK_BOX(hbox), label, FALSE, FALSE, 0);

/* layer subdivision */
hbox = gtk_hbox_new(FALSE, 0);
gtk_box_pack_start(GTK_BOX(vbox), hbox, FALSE, FALSE, 0);
label = gtk_label_new("Number of layers ");
gtk_box_pack_start(GTK_BOX(hbox), label, FALSE, FALSE, 0);

diffract_layer_total = gtk_spin_button_new_with_range(0, 10, 1);
gtk_spin_button_set_value(GTK_SPIN_BUTTON(diffract_layer_total), 0);
gtk_box_pack_end(GTK_BOX(hbox), diffract_layer_total, FALSE, FALSE, 0);

/* layer stacking */
hbox = gtk_hbox_new(FALSE, 0);
gtk_box_pack_start(GTK_BOX(vbox), hbox, TRUE, TRUE, PANEL_SPACING);

label = gtk_label_new("Stacking sequence ");
gtk_box_pack_start(GTK_BOX(hbox), label, FALSE, FALSE, 0);
diffract_layer_order = gtk_entry_new();
gtk_box_pack_end(GTK_BOX(hbox), diffract_layer_order, FALSE, FALSE, 0);

/* action buttons */
hbox2 = gtk_hbox_new(FALSE, 0);
gtk_box_pack_start(GTK_BOX(vbox), hbox2, FALSE, FALSE, PANEL_SPACING);
hbox = gtk_hbox_new(FALSE, PANEL_SPACING);
gtk_box_pack_end(GTK_BOX(hbox2), hbox, TRUE, FALSE, 0);

gui_icon_button(GTK_STOCK_APPLY, "Create ",
                  diffract_model_create, data,
                  hbox);
gui_icon_button(GTK_STOCK_SAVE, "Save  ",
                  diffract_save_dialog, data,
                  hbox);
#endif

/* frame */
vbox = gui_frame_vbox("Experimental typing", FALSE, FALSE, box);

w = gui_text_entry("FF label ", NULL, TRUE, FALSE, vbox);
dialog_child_set(dialog, "FF_label", w);

w = gui_text_entry("Neighbour Element ", NULL, TRUE, FALSE, vbox);
dialog_child_set(dialog, "FF_element", w);

w = gui_text_entry("Neighbour Distance ", NULL, TRUE, FALSE, vbox);
dialog_child_set(dialog, "FF_level", w);

w = gui_text_entry("Neighbour Count ", NULL, TRUE, FALSE, vbox);
dialog_child_set(dialog, "FF_count", w);

/* CURRENT - experimental */
/*
gtk_widget_set_sensitive(GTK_WIDGET(vbox), FALSE);
*/

/* apply single rule */
gui_button_x("Apply rule ", gui_make_rule, dialog, vbox);

vbox = gui_frame_vbox("Experimental - Import Forcefield", FALSE, FALSE, box);

hbox = gtk_hbox_new(FALSE, 0);
gtk_box_pack_start(GTK_BOX(vbox), hbox, FALSE, FALSE, 0);

w = gui_text_entry("Import GROMACS FF ", NULL, TRUE, FALSE, hbox);
dialog_child_set(dialog, "FF_import", w);

/* FIXME - for convenience during testing */
gtk_entry_set_text(GTK_ENTRY(w), "ffoplsaabon.itp");

gui_button_x(NULL, gui_import_ff, dialog, hbox);
}

/**************************/
/* dialog update function */
/**************************/
void gui_edit_refresh(void)
{
/* may be other updates */
gui_spatial_populate();
}

/****************************/
/* the model editing dialog */
/****************************/
void gui_edit_dialog(void)
{
gint i;
gpointer dialog;
GtkWidget *window, *frame, *label;
GtkWidget *notebook, *page;

/* request a new dialog */
dialog = dialog_request(CREATOR, "Model editing", gui_edit_refresh, NULL, NULL);
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
gtk_notebook_set_show_border(GTK_NOTEBOOK(notebook), FALSE);

/* add page */
page = gtk_vbox_new(FALSE,0);
label = gtk_label_new("Builder");
gtk_notebook_append_page(GTK_NOTEBOOK(notebook), page, label);
build_page(page);

/* add page */
page = gtk_vbox_new(FALSE,0);
label = gtk_label_new("Spatials");
gtk_notebook_append_page(GTK_NOTEBOOK(notebook), page, label);
spatial_page(page);

/* add page */
page = gtk_vbox_new(FALSE,0);
label = gtk_label_new("Transformations");
gtk_notebook_append_page(GTK_NOTEBOOK(notebook), page, label);
trans_page(page);

/* add page */
page = gtk_vbox_new(FALSE,0);
label = gtk_label_new("Labelling");
gtk_notebook_append_page(GTK_NOTEBOOK(notebook), page, label);
labelling_page(page, dialog);

/* add page */
page = gtk_vbox_new(FALSE,0);
label = gtk_label_new("Library");
gtk_notebook_append_page(GTK_NOTEBOOK(notebook), page, label);
gui_library_window(page);

/* terminating buttons */
gui_stock_button(GTK_STOCK_CLOSE, dialog_destroy, dialog,
                   GTK_DIALOG(window)->action_area);

gtk_widget_show_all(window);
 
/* init the transformation values */
reset_transmat();
for (i=0 ; i<12; i++)
  g_signal_connect(GTK_OBJECT(transmat[i]), "changed", 
                   GTK_SIGNAL_FUNC(change_transmat), (gpointer) i);
}

/* globals for the atom properties dialog */
GtkWidget *apd_label, *apd_type, *apd_charge, *apd_x, *apd_y, *apd_z;
GtkWidget *apd_growth, *apd_region, *apd_translate;
struct model_pak *apd_data=NULL;
struct core_pak *apd_core=NULL;

/***********************************************/
/* change the properties of all selected atoms */
/***********************************************/
/* NB: this is a bit dangerous as you change everything in the */
/* selection to the specified value - even (eg) incompatible elements */
void selection_properties_change(gint type)
{
gint n, growth, region, translate;
gdouble charge;
const gchar *text;
GSList *list;
struct elem_pak edata;
struct model_pak *model;
struct core_pak *core;

model = sysenv.active_model;
if (!model)
  return;
if (!model->selection)
  return;

switch (type)
  {
  case NAME:
    text = gtk_entry_get_text(GTK_ENTRY(apd_label));
    n = elem_symbol_test(text);
    if (n)
      {
      get_elem_data(n, &edata, model);
/* make sure we alow enough space for the string and the \0 */
      n = (LABEL_SIZE-1 > strlen(text)) ? strlen(text) : LABEL_SIZE-1;
      for (list=model->selection ; list ; list=g_slist_next(list))
        {
        core = list->data;
        g_free(core->atom_label);
        core->atom_label = g_strdup(text);
/* update atttached shell */
        if (core->shell)
          {
          struct shel_pak *shell = core->shell;
          g_free(shell->shell_label);
          shell->shell_label = g_strdup(text);
          }
/* NEW - don't update element specific data if the element type was not */
/* changed - ie the user has just made a labelling change (eg C -> C1) */
        if (edata.number != core->atom_code)
          {
          core->atom_code = edata.number;
          core->bond_cutoff = edata.cova;
          }
        init_atom_colour(core, model);
        init_atom_charge(core, model);
        }
/* model updates */
      g_slist_free(model->unique_atom_list);
      model->unique_atom_list = find_unique(ELEMENT, model);
      calc_emp(model);
      }
    break;

  case CHARGE:
    charge = str_to_float(gtk_entry_get_text(GTK_ENTRY(apd_charge)));
    for (list=model->selection ; list ; list=g_slist_next(list))
      {
      core = list->data;
/* core updates */
      core->charge = charge;
      core->lookup_charge = FALSE;
      }
    calc_emp(model);
    break;

  case CORE_GROWTH_SLICE:
    growth = str_to_float(gtk_entry_get_text(GTK_ENTRY(apd_growth)));
    growth = CLAMP(growth, 0, 1);
    for (list=model->selection ; list ; list=g_slist_next(list))
      {
      core = list->data;
      core->growth = growth;
      if (model->colour_scheme == GROWTH_SLICE)
        atom_colour_scheme(GROWTH_SLICE, core, model);
      }
    break;

  case CORE_REGION:
    region = str_to_float(gtk_entry_get_text(GTK_ENTRY(apd_region)));
    if (region > model->region_max)
      model->region_max = region;

    for (list=model->selection ; list ; list=g_slist_next(list))
      {
      core = list->data;
      core->region = region;
      if (core->shell)
        (core->shell)->region = region;

      if (model->colour_scheme == REGION)
        atom_colour_scheme(REGION, core, model);
      }
    break;

  case CORE_TRANSLATE:
    translate = str_to_float(gtk_entry_get_text(GTK_ENTRY(apd_translate)));
    translate = CLAMP(translate, 0, 1);
    for (list=model->selection ; list ; list=g_slist_next(list))
      {
      core = list->data;
      core->translate = translate;
      if (core->shell)
        (core->shell)->translate = translate;

      if (model->colour_scheme == TRANSLATE)
        atom_colour_scheme(TRANSLATE, core, model);
      }
    break;

  case CORE_FF:
    text = gtk_entry_get_text(GTK_ENTRY(apd_type));
    for (list=model->selection ; list ; list=g_slist_next(list))
      {
      core = list->data;
      if (core->atom_type)
        g_free(core->atom_type);
      core->atom_type = g_strdup(text);
      }
    break;
  }

gui_refresh(GUI_MODEL_PROPERTIES);
gui_refresh(GUI_CANVAS);
}

/*****************************/
/* commit changes to an atom */
/*****************************/
void atom_properties_change(GtkWidget *w, gint type)
{
gint n, growth, region, translate;
const gchar *text;
struct elem_pak edata;
struct model_pak *model;

model = sysenv.active_model;
if (!model)
  return;

/* act on multiple atoms? */
if (g_slist_length(model->selection) > 1)
  {
  selection_properties_change(type);
  return;
  }

if (!apd_core)
  return;

switch(type)
  {
  case NAME:
    text = gtk_entry_get_text(GTK_ENTRY(apd_label));
    g_free(apd_core->atom_label);
    apd_core->atom_label = g_strdup(text);
    n = elem_symbol_test(text);

/* update atttached shell */
if (apd_core->shell)
  {
  struct shel_pak *shell = apd_core->shell;
  g_free(shell->shell_label);
  shell->shell_label = g_strdup(text);
  }

/* if recognized -> update */
    if (n)
      {
      get_elem_data(n, &edata, model);

/* NEW - don't update element specific data if the element type was not */
/* changed - ie the user has just made a labelling change (eg C -> C1) */
      if (n != apd_core->atom_code)
        {
        apd_core->atom_code = n;
        apd_core->bond_cutoff = edata.cova;
        }
      init_atom_colour(apd_core, model);
      init_atom_charge(apd_core, model);

      g_slist_free(model->unique_atom_list);
      model->unique_atom_list = find_unique(ELEMENT, model);
      calc_emp(model);

/* REFRESH */
      gui_refresh(GUI_MODEL_PROPERTIES);
      }
    break;

  case CORE_FF:
    text = gtk_entry_get_text(GTK_ENTRY(apd_type));
    if (apd_core->atom_type)
      g_free(apd_core->atom_type);
    apd_core->atom_type = g_strdup(text);
    break;

  case CHARGE:
    text = gtk_entry_get_text(GTK_ENTRY(apd_charge));
    apd_core->charge = str_to_float(text);
    apd_core->lookup_charge = FALSE;
    calc_emp(model);
    break;

  case COORD_X:
    text = gtk_entry_get_text(GTK_ENTRY(apd_x));
    apd_core->x[0] = str_to_float(text);
    coords_compute(model);
    break;

  case COORD_Y:
    text = gtk_entry_get_text(GTK_ENTRY(apd_y));
    apd_core->x[1] = str_to_float(text);
    coords_compute(model);
    break;

  case COORD_Z:
    text = gtk_entry_get_text(GTK_ENTRY(apd_z));
    apd_core->x[2] = str_to_float(text);
    coords_compute(model);
    break;

  case CORE_GROWTH_SLICE:
    text = gtk_entry_get_text(GTK_ENTRY(apd_growth));
    growth = CLAMP(str_to_float(text), 0, 1);
    apd_core->growth = growth;
    if (model->colour_scheme == GROWTH_SLICE)
      atom_colour_scheme(GROWTH_SLICE, apd_core, model);
    break;

  case CORE_REGION:
    text = gtk_entry_get_text(GTK_ENTRY(apd_region));
    region = str_to_float(text);
    if (region > model->region_max)
      model->region_max = region;
    apd_core->region = region;
    if (apd_core->shell)
      (apd_core->shell)->region = region;
    if (model->colour_scheme == REGION)
      atom_colour_scheme(REGION, apd_core, model);
    break;

  case CORE_TRANSLATE:
    text = gtk_entry_get_text(GTK_ENTRY(apd_translate));
    translate = CLAMP(str_to_float(text), 0, 1);
    apd_core->translate = translate;
    if (apd_core->shell)
      (apd_core->shell)->translate = translate;
    if (model->colour_scheme == TRANSLATE)
      atom_colour_scheme(TRANSLATE, apd_core, model);
    break;

  default:
    printf("Not yet modifiable...\n");
  }
gui_refresh(GUI_CANVAS);
}

/*************************************/
/* updates the dialog for a new atom */
/*************************************/
void gui_refresh_selection(void)
{
gint n, cflag=FALSE;
gdouble q, centroid[3];
gchar *element, *label, *type, *charge, *x, *y, *z, *growth, *region, *translate;
struct core_pak *core;
struct model_pak *model;
GSList *list;

model = sysenv.active_model;
core = NULL;
if (model)
  {
  list = model->selection;
  n = g_slist_length(list);
  switch (n)
    {
    case -1:
      g_assert_not_reached();
    case 0:
      break;
    case 1:
      core = list->data;
      break;

    default:
      VEC3SET(centroid, 0.0, 0.0, 0.0);
      for (list=model->selection ; list ; list=g_slist_next(list))
        {
        core = list->data;
        ARR3ADD(centroid, core->x);
        }
      VEC3MUL(centroid, 1.0 / (gdouble) n);
/* special print case - centroid display */
      cflag = TRUE;
      core = NULL;
    }
  }

if (core && model)
  {
/* data available */
  element = g_strdup(elements[core->atom_code].symbol);
  label = g_strdup(core->atom_label);

  if (core->atom_type)
    type = g_strdup(core->atom_type);
  else
    type = g_strdup("");

  q = atom_charge(core); /* Replaced by C. Fisher 2004 */
  charge = g_strdup_printf("%9.4f", q);

  x = g_strdup_printf("%9.4f", core->x[0]);
  y = g_strdup_printf("%9.4f", core->x[1]);
  z = g_strdup_printf("%9.4f", core->x[2]);

  growth = g_strdup_printf("%d", core->growth);
  region = g_strdup_printf("%d", core->region);
  translate = g_strdup_printf("%d", core->translate);

  apd_core = core;
  }
else
  {
/* otherwise defaults */
  element = g_strdup("");
  type = g_strdup("");
  charge = g_strdup("");

  if (cflag)
    {
    label = g_strdup("centroid");
    x = g_strdup_printf("%9.4f", centroid[0]);
    y = g_strdup_printf("%9.4f", centroid[1]);
    z = g_strdup_printf("%9.4f", centroid[2]);
    }
  else
    {
    label = g_strdup("");
    x = g_strdup("");
    y = g_strdup("");
    z = g_strdup("");
    }

  growth = g_strdup("");
  region = g_strdup("");
  translate = g_strdup("");
  }

/* prevent changes from messing up the atom_properties_change() callback */
apd_data = NULL;

/* entry updates */
gtk_entry_set_text(GTK_ENTRY(apd_label), label);
gtk_entry_set_text(GTK_ENTRY(apd_type), type);
gtk_entry_set_text(GTK_ENTRY(apd_charge), charge);
gtk_entry_set_text(GTK_ENTRY(apd_x), x);
gtk_entry_set_text(GTK_ENTRY(apd_y), y);
gtk_entry_set_text(GTK_ENTRY(apd_z), z);
gtk_entry_set_text(GTK_ENTRY(apd_growth), growth);
gtk_entry_set_text(GTK_ENTRY(apd_region), region);
gtk_entry_set_text(GTK_ENTRY(apd_translate), translate);

apd_data = model;

/* cleanup */
g_free(element);
g_free(label);
g_free(type);
g_free(charge);
g_free(x);
g_free(y);
g_free(z);
g_free(growth);
g_free(region);
g_free(translate);
}

/*******************************************/
/* display the properties of a single atom */
/*******************************************/
void gui_edit_widget(GtkWidget *box)
{
GtkWidget *frame, *hbox, *vbox, *entry;

/* checks */
g_return_if_fail(box != NULL);

/* two column element data display */
hbox = gtk_hbox_new(FALSE, 0);
gtk_box_pack_start(GTK_BOX(box), hbox, TRUE, TRUE, 0);

/* left vbox - titles */
vbox = gtk_vbox_new(TRUE, 0);
gtk_box_pack_start(GTK_BOX(hbox), vbox, TRUE, TRUE, 0);

/* TODO - put in a for loop? */

entry = gtk_entry_new();
gtk_entry_set_text(GTK_ENTRY(entry), "Label");
gtk_entry_set_editable(GTK_ENTRY(entry), FALSE);
gtk_box_pack_start(GTK_BOX(vbox), entry, FALSE, FALSE, 0);

entry = gtk_entry_new();
gtk_entry_set_text(GTK_ENTRY(entry), "FF Type");
gtk_entry_set_editable(GTK_ENTRY(entry), FALSE);
gtk_box_pack_start(GTK_BOX(vbox), entry, FALSE, FALSE, 0);

entry = gtk_entry_new();
gtk_entry_set_text(GTK_ENTRY(entry), "X");
gtk_entry_set_editable(GTK_ENTRY(entry), FALSE);
gtk_box_pack_start(GTK_BOX(vbox), entry, FALSE, FALSE, 0);

entry = gtk_entry_new();
gtk_entry_set_text(GTK_ENTRY(entry), "Y");
gtk_entry_set_editable(GTK_ENTRY(entry), FALSE);
gtk_box_pack_start(GTK_BOX(vbox), entry, FALSE, FALSE, 0);

entry = gtk_entry_new();
gtk_entry_set_text(GTK_ENTRY(entry), "Z");
gtk_entry_set_editable(GTK_ENTRY(entry), FALSE);
gtk_box_pack_start(GTK_BOX(vbox), entry, FALSE, FALSE, 0);

entry = gtk_entry_new();
gtk_entry_set_text(GTK_ENTRY(entry), "Charge");
gtk_entry_set_editable(GTK_ENTRY(entry), FALSE);
gtk_box_pack_start(GTK_BOX(vbox), entry, FALSE, FALSE, 0);

entry = gtk_entry_new();
gtk_entry_set_text(GTK_ENTRY(entry), "Growth");
gtk_entry_set_editable(GTK_ENTRY(entry), FALSE);
gtk_box_pack_start(GTK_BOX(vbox), entry, FALSE, FALSE, 0);

entry = gtk_entry_new();
gtk_entry_set_text(GTK_ENTRY(entry), "Region");
gtk_entry_set_editable(GTK_ENTRY(entry), FALSE);
gtk_box_pack_start(GTK_BOX(vbox), entry, FALSE, FALSE, 0);

entry = gtk_entry_new();
gtk_entry_set_text(GTK_ENTRY(entry), "Translate");
gtk_entry_set_editable(GTK_ENTRY(entry), FALSE);
gtk_box_pack_start(GTK_BOX(vbox), entry, FALSE, FALSE, 0);

/* right vbox - data */
vbox = gtk_vbox_new(TRUE, 0);
gtk_box_pack_end(GTK_BOX(hbox), vbox, TRUE, TRUE, 0);

apd_label = gtk_entry_new();
gtk_box_pack_start(GTK_BOX(vbox), apd_label, FALSE, FALSE, 0);

apd_type = gtk_entry_new();
gtk_box_pack_start(GTK_BOX(vbox), apd_type, FALSE, FALSE, 0);

apd_x = gtk_entry_new();
gtk_box_pack_start(GTK_BOX(vbox), apd_x, FALSE, FALSE, 0);

apd_y = gtk_entry_new();
gtk_box_pack_start(GTK_BOX(vbox), apd_y, FALSE, FALSE, 0);

apd_z = gtk_entry_new();
gtk_box_pack_start(GTK_BOX(vbox), apd_z, FALSE, FALSE, 0);

apd_charge = gtk_entry_new();
gtk_box_pack_start(GTK_BOX(vbox), apd_charge, FALSE, FALSE, 0);

apd_growth = gtk_entry_new();
gtk_box_pack_start(GTK_BOX(vbox), apd_growth, FALSE, FALSE, 0);

apd_region = gtk_entry_new();
gtk_box_pack_start(GTK_BOX(vbox), apd_region, FALSE, FALSE, 0);

apd_translate = gtk_entry_new();
gtk_box_pack_start(GTK_BOX(vbox), apd_translate, FALSE, FALSE, 0);

/* attach callbacks (NB: set initial data first) */
g_signal_connect(GTK_OBJECT(apd_label), "activate",
                 GTK_SIGNAL_FUNC(atom_properties_change), GINT_TO_POINTER(NAME));
g_signal_connect(GTK_OBJECT(apd_type), "activate",
                 GTK_SIGNAL_FUNC(atom_properties_change), GINT_TO_POINTER(CORE_FF));

g_signal_connect(GTK_OBJECT(apd_x), "activate",
                 GTK_SIGNAL_FUNC(atom_properties_change), GINT_TO_POINTER(COORD_X));
g_signal_connect(GTK_OBJECT(apd_y), "activate",
                 GTK_SIGNAL_FUNC(atom_properties_change), GINT_TO_POINTER(COORD_Y));
g_signal_connect(GTK_OBJECT(apd_z), "activate",
                 GTK_SIGNAL_FUNC(atom_properties_change), GINT_TO_POINTER(COORD_Z));

g_signal_connect(GTK_OBJECT(apd_charge), "activate",
                 GTK_SIGNAL_FUNC(atom_properties_change), GINT_TO_POINTER(CHARGE));

g_signal_connect(GTK_OBJECT(apd_growth), "activate",
                 GTK_SIGNAL_FUNC(atom_properties_change),
                 GINT_TO_POINTER(CORE_GROWTH_SLICE));

g_signal_connect(GTK_OBJECT(apd_region), "activate",
                 GTK_SIGNAL_FUNC(atom_properties_change), GINT_TO_POINTER(CORE_REGION));
g_signal_connect(GTK_OBJECT(apd_translate), "activate",
                 GTK_SIGNAL_FUNC(atom_properties_change), GINT_TO_POINTER(CORE_TRANSLATE));

/* CURRENT */
frame = gtk_frame_new(NULL);
gtk_box_pack_start(GTK_BOX(box), frame, FALSE, FALSE, 0);
gtk_container_set_border_width(GTK_CONTAINER(frame), PANEL_SPACING);
vbox = gtk_vbox_new(TRUE, 1);
gtk_container_add(GTK_CONTAINER(frame), vbox);

gui_button_x("Add atoms", gtk_mode_switch, (gpointer) ATOM_ADD, vbox);
gui_button_x("Add bonds", gtk_mode_switch, (gpointer) BOND_SINGLE, vbox);
gui_button_x("Delete bonds", gtk_mode_switch, (gpointer) BOND_DELETE, vbox);
gui_button_x("Normal mode", gtk_mode_switch, (gpointer) FREE, vbox);

frame = gtk_frame_new(NULL);
gtk_box_pack_start(GTK_BOX(box), frame, FALSE, FALSE, 0);
gtk_container_set_border_width(GTK_CONTAINER(frame), PANEL_SPACING);
vbox = gtk_vbox_new(TRUE, 1);
gtk_container_add(GTK_CONTAINER(frame), vbox);

gui_button_x("Mark as ghost", select_flag_ghost, NULL, vbox);
gui_button_x("Mark as normal", select_flag_normal, NULL, vbox);
}
