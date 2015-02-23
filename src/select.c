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
#include <string.h>

#include "gdis.h"
#include "coords.h"
#include "edit.h"
#include "matrix.h"
#include "quaternion.h"
#include "measure.h"
#include "spatial.h"
#include "opengl.h"
#include "render.h"
#include "select.h"
#include "zone.h"
#include "interface.h"
#include "gui_shorts.h"

extern struct sysenv_pak sysenv;
extern struct elem_pak elements[];

gint backbone = 0;

/*********************************/
/* update selection connectivity */
/*********************************/
void select_connect_update(struct model_pak *model)
{
GSList *list;

/* remove all old bonds first */
for (list=model->selection ; list ; list=g_slist_next(list))
  connect_atom_clear(list->data, model);

/* now create the new connectivity */
for (list=model->selection ; list ; list=g_slist_next(list))
  connect_atom_compute(list->data, model);

connect_molecules(model);
}

/****************************/
/* mark selection as ghosts */
/****************************/
void select_flag_ghost(void)
{
GSList *list;
struct core_pak *core;
struct model_pak *data;

data = sysenv.active_model;
if (!data)
  return;

/* clean up selection highlighting */
for (list=data->selection ; list ; list=g_slist_next(list))
  {
  core = list->data;
  core->ghost = TRUE;
  init_atom_colour(core, data);
  }
redraw_canvas(SINGLE);
}

/****************************/
/* mark selection as normal */
/****************************/
void select_flag_normal(void)
{
GSList *list;
struct core_pak *core;
struct model_pak *data;

data = sysenv.active_model;
if (!data)
  return;

/* clean up selection highlighting */
for (list=data->selection ; list ; list=g_slist_next(list))
  {
  core = list->data;
  core->ghost = FALSE;
  init_atom_colour(core, data);
  }
redraw_canvas(SINGLE);
}

/***********************************/
/* change selection rendering mode */
/***********************************/
void core_render_mode_set(gint mode, GSList *cores)
{
GSList *list;
struct core_pak *core;

for (list=cores ; list ; list=g_slist_next(list))
  {
  core = list->data;
  core->render_mode = mode;
  }
}

void core_render_wire_set(gint wire, GSList *cores)
{
GSList *list;
struct core_pak *core;

for (list=cores ; list ; list=g_slist_next(list))
  {
  core = list->data;
  core->render_wire = wire;
  }
}

/*******************/
/* clear selection */
/*******************/
void select_clear(struct model_pak *model)
{
GSList *list;
struct core_pak *core;

g_assert(model != NULL);

/* clean up selection highlighting */
for (list=model->selection ; list ; list=g_slist_next(list))
  {
  core = list->data;
  core->status &= ~SELECT;
  }
g_slist_free(model->selection);
model->selection=NULL;
}

/******************/
/* copy selection */
/******************/
void select_copy(void)
{
struct model_pak *model;

/* setup & check */
model = sysenv.active_model;
if (!model)
  return;
if (!model->selection)
  return;

sysenv.select_source = model;
}

/*********************/
/* paste a selection */
/*********************/
void select_paste(void)
{
struct model_pak *data, *src;
struct core_pak *core, *orig;
GSList *clist, *slist;

/* checks */
data = sysenv.active_model;
src = sysenv.select_source;
if (data == NULL || src == NULL)
  return;
if (data->id == NODATA || !src->selection)
  return;

/* copy the selection - just in case src == data */
slist = g_slist_copy(src->selection);
/* ensure pasted atoms will be the only ones that are selected */
select_clear(data);

/* duplicate cores */
for (clist=slist ; clist ; clist=g_slist_next(clist))
  {
  orig = clist->data;

/* duplicate the core & transform into the new model's lattice space */
  core = copy_core(orig, src, data);

/* NEW - if destination model is 2D, force to region 1 */
  if (src->periodic == 2)
    core->region = REGION1A;

/* add core to the selection */
  select_add_core(core, data);
  }

g_slist_free(slist);

/* updates */
coords_compute(data);
zone_init(data);
connect_bonds(data);
connect_molecules(data);
redraw_canvas(SINGLE);
}

/*****************************/
/* colour selection callback */
/*****************************/
void cb_select_colour(gpointer *csd)
{
gdouble colour[3];
GSList *list;
struct core_pak *core;
struct model_pak *model;

model = sysenv.active_model;
if (!model)
  return;

gtk_color_selection_get_color((GtkColorSelection *) csd, colour);

/*
P3VEC("colour: ", colour);
*/

for (list=model->selection ; list ; list=g_slist_next(list))
  {
  core = list->data;

  ARR3SET(core->colour, colour);
  VEC3MUL(core->colour, 65535.0);
  }

/* TODO - close the csd */
redraw_canvas(SINGLE);
}

/********************/
/* colour selection */
/********************/
void select_colour(void)
{
#ifdef WITH_GUI
/* no action if nothing is loaded */
if (sysenv.active_model)
  new_csd("Selection colour", (gpointer) cb_select_colour);
#endif
}

/********************/
/* delete selection */
/********************/
void select_delete(void)
{
GSList *list;
struct model_pak *data;
struct core_pak *core;

/* deletion for the active model only */
data = sysenv.active_model;
if (!data)
  return;

/* delete */
list = data->selection;
while (list)
  {
  core = list->data;
  list = g_slist_next(list);
/*
  delete_core(core);
*/

/* NEW */
  core_delete_single(core, data);

  }

/*
delete_commit(data);
*/

data->selection = NULL;
sysenv.select_source = NULL;

redraw_canvas(SINGLE);
}

/******************/
/* hide selection */
/******************/
void select_hide(void)
{
GSList *list;
struct model_pak *data;
struct core_pak *core;

/* deletion for the active model only */
data = sysenv.active_model;
if (!data)
  return;
/* exit if nothing selected */
if (!data->selection)
  return;

/* hide */
for (list=data->selection ; list ; list=g_slist_next(list))
  {
  core = list->data;
  core->status |= HIDDEN;
  core->status &= ~SELECT;
  }
sysenv.select_source = NULL;

/* clear the selection, so subsequent selection operations  */
/* (ctrl+whatever) don't move invisible things around & cause havoc */
g_slist_free(data->selection);
data->selection = NULL;

/* update */
redraw_canvas(SINGLE);
}

/********************/
/* select all cores */
/********************/
void select_all(void)
{
GSList *list;
struct core_pak *core;
struct model_pak *data;

/* deletion for the active model only */
data = sysenv.active_model;
if (!data)
  return;

/* assign the new selection */
g_slist_free(data->selection);
data->selection = g_slist_copy(data->cores);

/* update the highlighting */
for (list=data->selection ; list ; list=g_slist_next(list))
  {
  core = list->data;
  core->status |= SELECT;
  }
redraw_canvas(SINGLE);
}

/************************/
/* invert the selection */
/************************/
void select_invert(void)
{
GSList *new, *list;
struct core_pak *core;
struct model_pak *data;

/* deletion for the active model only */
data = sysenv.active_model;
if (!data)
  return;

/* copy the core list */
new = g_slist_copy(data->cores);

/* remove the old selection */
for (list=data->selection ; list ; list=g_slist_next(list))
  {
  core = list->data;

  core->status &= ~SELECT;
  new = g_slist_remove(new, core);
  }

/* assign the new selection */
g_slist_free(data->selection);
data->selection = new;

/* update the highlighting */
for (list=data->selection ; list ; list=g_slist_next(list))
  {
  core = list->data;

  core->status |= SELECT;
  }
redraw_canvas(SINGLE);
}

/*************************************************/
/* request that a core be added to the selection */
/*************************************************/
#define DEBUG_ADD_SELECT 0
gint select_add_core(struct core_pak *core, struct model_pak *data)
{
g_return_val_if_fail(core != NULL, 0);
g_return_val_if_fail(data != NULL, 0);

/* should we add? */
if (g_slist_find(data->selection, core))
  return(1);
if (core->status & HIDDEN)
  return(0);

data->selection = g_slist_append(data->selection, core);
core->status |= SELECT;

gui_refresh(GUI_MODEL_PROPERTIES);
return(0);
}

/***********************************************/
/* select all cores of a given atom label type */
/***********************************************/
void select_all_types(struct core_pak *core, struct model_pak *model)
{
GSList *list;
struct core_pak *comp;

for (list=model->cores ; list ; list=g_slist_next(list))
  {
  comp = list->data;

  if (core->atom_type)
    if (g_ascii_strcasecmp(comp->atom_type, core->atom_type) == 0)
      select_add_core(comp, model);
  }
}

/***********************************************/
/* select all cores of a given atom label type */
/***********************************************/
void select_all_labels(struct core_pak *core, struct model_pak *model)
{
GSList *list;
struct core_pak *comp;

for (list=model->cores ; list ; list=g_slist_next(list))
  {
  comp = list->data;
 
  if (g_ascii_strcasecmp(comp->atom_label, core->atom_label) == 0)
    select_add_core(comp, model);
  }
}

/**********************************************************/
/* select all elements of given type from a list of cores */
/**********************************************************/
void select_all_elem(gint code, struct model_pak *model)
{
GSList *list;
struct core_pak *core;

g_assert(model != NULL);

for (list=model->cores ; list ; list=g_slist_next(list))
  {
  core = list->data;
 
  if (core->atom_code == code)
    select_add_core(core, model);
  }
}

/**********************************************************/
/* select all elements of given type from a list of cores */
/**********************************************************/
void select_all_elem_in_mol(struct core_pak *core, struct model_pak *model)
{
GSList *list;
struct core_pak *comp;

g_assert(model != NULL);
g_assert(core != NULL);

for (list=model->cores ; list ; list=g_slist_next(list))
  {
  comp = list->data;
  if (core->atom_code == comp->atom_code && core->molecule == comp->molecule)
    select_add_core(comp, model);
  }
}

/**********************************************/
/* add a particular molecule to the selection */
/**********************************************/
void select_add_mol(struct mol_pak *mol, struct model_pak *model)
{
GSList *list;
struct core_pak *core;

g_assert(model != NULL);
if (!mol)
  return;

for (list=mol->cores ; list ; list=g_slist_next(list))
  {
  core = list->data;
  select_add_core(core, model);
  }
}

/**********************************************/
/* add a particular molecule to the selection */
/**********************************************/
#define DEBUG_ADD_FRAGMENT 0
void select_add_fragment(struct core_pak *core, struct model_pak *model)
{
GSList *list, *cores;
static struct core_pak *first=NULL;

/* test for first or second call */
if (model->state)
  {
  g_assert(first != NULL);

/* aqcuire the fragment */
  connect_fragment_init(model);
  cores = connect_fragment_get(first, core, model);

/* add all cores to the selection */
  select_add_core(first, model);
  for (list=cores ; list ; list=g_slist_next(list))
    select_add_core(list->data, model);

  g_slist_free(list);
  model->state = 0;
  first = NULL;
  }
else
  {
  first = core;
  select_add_core(core, model);
  model->state++;
  }
}

/******************************************/
/* add all cores of the same region label */
/******************************************/
void select_add_region(struct core_pak *core, struct model_pak *model)
{
GSList *list;
struct core_pak *comp;

g_assert(model != NULL);
g_assert(core != NULL);

for (list=model->cores ; list ; list=g_slist_next(list))
  {
  comp = list->data;
  if (core->region == comp->region)
    select_add_core(comp, model);
  }
}

/******************************************/
/* remove an atom from the selection list */
/******************************************/
void select_del_core(struct core_pak *core, struct model_pak *model)
{
g_assert(model != NULL);
g_assert(core != NULL);

core->status &= ~SELECT;
model->selection = g_slist_remove(model->selection, core);
}

/******************************************************/
/* select atom(s) according to current selection mode */
/******************************************************/
void select_core(struct core_pak *core, gint toggle, struct model_pak *model)
{
switch (sysenv.select_mode)
  {
  case ATOM_LABEL:
    select_all_labels(core, model);
    break;
  case ATOM_TYPE:
    select_all_types(core, model);
    break;
  case CORE:
/* toggle = true -> if selected already; deselect */
/* toggle = false -> if selected already; remain selected */
    if (toggle)
      {
      if (select_add_core(core, model))
        select_del_core(core, model);
      }
    else
      select_add_core(core, model);
    break;
  case ELEM:
    select_all_elem(core->atom_code, model);
    break;
  case ELEM_MOL:
    select_all_elem_in_mol(core, model);
    break;
  case MOL:
    select_add_mol(core->mol, model);
    break;
  case FRAGMENT:
    select_add_fragment(core, model);
    break;
  case REGION:
    select_add_region(core, model);
    break;
  }
}

/********************************/
/* rotate the current selection */
/********************************/
#define DEBUG_ROTATE_SELECT 0
void rotate_select(struct model_pak *data, gdouble *mat)
{
gint n;
gdouble scent[3], vec[3];
GSList *list;
struct core_pak *core;
struct shel_pak *shell;

g_assert(data != NULL);

n = g_slist_length(data->selection);
if (!n)
  return;

#if DEBUG_ROTATE_SELECT
P3MAT("rot: ", mat);
#endif

/* calculate selection's centroid */
VEC3SET(scent, 0.0, 0.0, 0.0);
for (list=data->selection ; list ; list=g_slist_next(list))
  {
  core = list->data;
  ARR3ADD(scent, core->x);
  }
VEC3MUL(scent, 1.0 / (gdouble) n);

/* FIXME - why is this different for INITIAL orientation */
#if DEBUG_ROTATE_SELECT
P3VEC("scent: ", scent);
#endif

for (list=data->selection ; list ; list=g_slist_next(list))
  {
  core = list->data;

  ARR3SET(vec, core->x);
  ARR3SUB(vec, scent);
  vecmat(mat, vec);
  ARR3ADD(vec, scent);
  ARR3SET(core->x, vec);

/* shells */
  if (core->shell)
    {
    shell = core->shell;

    ARR3SET(vec, shell->x);
    ARR3SUB(vec, scent);

/* rotate to get on screen positions, and store the values */
    vecmat(mat, vec);
    ARR3ADD(vec, scent);
    ARR3SET(shell->x, vec);
    }
  }

/* updates */
coords_compute(data);
zone_init(data);

/* FIXME - fine grained approach sometimes misses bonds */
/*
select_connect_update(data);
*/
connect_bonds(data);
connect_molecules(data);

meas_graft_model(data);
}

/*************************/
/* translate a selection */
/*************************/
void select_translate(gint px, gint py, struct model_pak *model)
{
gdouble x[3], r[3], qi[9], tm[9], ti[9];
GSList *list;
struct core_pak *core;
struct shel_pak *shel;

g_assert(model != NULL);

/* translation */
VEC3SET(x, px, 0.0, -py);
VEC3MUL(x, 0.01);

/* lattice transformation */
memcpy(ti, model->latmat, 9*sizeof(gdouble));

/* camera transformation */
quat_matrix(tm, camera_q(model->camera));
quat_matrix(qi, camera_q(model->camera));
matrix_invert(qi);

/* build transform matrices */
matmat(qi, ti);
matmat(model->ilatmat, tm);

for (list=model->selection ; list ; list=g_slist_next(list))
  {
  core = list->data;

/* transform */
  ARR3SET(r, core->x);
  vecmat(ti, r);
  ARR3ADD(r, x);
  vecmat(tm, r);
  ARR3SET(core->x, r);

  if (core->shell)
    {
    shel = core->shell;

/* transform */
    ARR3SET(r, shel->x);
    vecmat(ti, r);
    ARR3ADD(r, x);
    vecmat(tm, r);
    ARR3SET(shel->x, r);
    }
  }

/* updates */
coords_compute(model);
zone_init(model);
/* FIXME - fine grained approach sometimes misses bonds */
/*
select_connect_update(model);
*/
connect_bonds(model);
connect_molecules(model);

meas_graft_model(model);
}

/************************************/
/* total pathways from a given atom */
/************************************/
gint total_branches(struct core_pak *core)
{
gint branches=0;
GSList *list;
struct bond_pak *bond;
struct core_pak *exit;

for (list=core->bonds ; list ; list=g_slist_next(list))
  {
  bond = list->data;

  if (core == bond->atom1)
    exit = bond->atom2;
  else
    exit = bond->atom1;

/* NEW - preprocess has removed all the dead end branches */
  if (exit->status & PRUNED)
    continue;

  if (backbone)
    if (exit->atom_code == backbone)
      branches++;
  }
return(branches);
}

/****************************************/
/* exit pathways from atom2, given 1->2 */
/****************************************/
GSList *exit_branches(struct core_pak *core1, struct core_pak *core2)
{
GSList *list1, *exits=NULL;
struct bond_pak *bond;
struct core_pak *exit;

/* checks */
if (backbone)
  if (core2->atom_code != backbone)
    return(NULL);

for (list1=core2->bonds ; list1 ; list1=g_slist_next(list1))
  {
  bond = list1->data;
  if (bond->type == BOND_HBOND)
    continue;

/* get exits */
  if (core2 == bond->atom1)
    exit = bond->atom2;
  else
    exit = bond->atom1;

/* ignore entry point */
  if (exit == core1)
    continue;

/* ignore marked (pruned) branches */
  if (exit->status & PRUNED)
    continue;

/* if constructing ribbon first atom is backbone type */
  if (backbone)
    {
    if (exit->atom_code != backbone)
       continue;
    }

  exits = g_slist_prepend(exits, exit);
  }

return(exits);
}

/**************************/
/* node creation and init */
/**************************/
#define DEBUG_CREATE_NODE 0
struct node_pak *create_node(struct core_pak *core)
{
struct node_pak *node_data;

node_data = g_malloc(sizeof(struct node_pak));
node_data->core = core;
ARR3SET(node_data->x, core->x);
node_data->num_branches = total_branches(core);

#if DEBUG_CREATE_NODE
printf(" >>> NEW NODE @ atom: %s [%p], order: %d\n", core->label, core,
                                              node_data->num_branches);
#endif

/* number of exits traversed */
node_data->num_explored = 0;

/* these are filled out by the traverse function */
node_data->exit_list=NULL;
node_data->link_list=NULL;

return(node_data);
}

/***************/
/* find a node */
/***************/
struct node_pak *get_node_by_atom(struct core_pak *core, GList *node_list)
{
GList *list;
struct node_pak *node_data;

for (list=node_list ; list ; list=g_list_next(list))
  {
  node_data = list->data;

  if (node_data->core == core)
    return(node_data);
  }

return(NULL);
}

/**************************/
/* link creation and init */
/**************************/
#define DEBUG_CREATE_LINK 0
struct link_pak *create_link(GList *chain)
{
GList *list;
struct link_pak *link;
static gint link_id=0;

/* first time init */
if (!chain)
  {
  link_id = 0;
  return(NULL);
  }

link = g_malloc(sizeof(struct link_pak));
link->chain = NULL;
link->size = 0;
link->id = link_id++;

list = g_list_last(chain);
link->last = list->data;
list = chain;
link->first = list->data;

while (list)
  {
/* need to do this? */

  link->chain = g_list_append(link->chain, list->data);
  link->size++;
  list = g_list_next(list); 
  }

#if DEBUG_CREATE_LINK
printf("Creating link [%p:%p], size = %d\n", link->first, link->last, link->size);
#endif

return(link);
}

/*****************************************************/
/* return the atom number for the next in a sequence */
/*****************************************************/
#define DEBUG_TRAVERSE 0
struct core_pak *traverse(struct node_pak *node, GList *chain)
{
GList *list;
GSList *slist, *slist2;
struct link_pak *last_link;
struct core_pak *curr=NULL, *prev=NULL, *next=NULL;

/* FIXME - if chain is NULL -> previous list on current node? */

/* current atom */
list = g_list_last(chain);
curr = list->data;

/* attempt to get the atom we last left */
list = g_list_previous(list);

if (list)
  prev = list->data;
else
  {
/* last submitted link of current node */
  list = g_list_last(node->link_list);
/* last atom in this list should be the current one - get previous */
  if (list)
    {
    last_link = list->data;
    list = g_list_last(last_link->chain);
    list = g_list_previous(list);

    g_assert(list != NULL);
    prev = list->data;
    }
  }

/* get available pathways */
slist = exit_branches(prev, curr);

/* are we at a node? */
if (node->core == curr)
  {
/* go through candidate pathways */
  while (slist)
    {
    next = slist->data;
/* if explored previously, then skip */
    slist2 = node->exit_list;
    while (slist2)
      {
      if (next == slist2->data)
        break;
      slist2 = g_slist_next(slist2);
      }

/* no match found - this pathways is unexplored */
    if (!slist2)
      break;

    slist = g_slist_next(slist);
    }

/* that's one more branch done */
  if (slist)
    {
    node->exit_list = g_slist_append(node->exit_list, next);
    node->num_explored++;
    }
  else
    g_assert_not_reached();

#if DEBUG_TRAVERSE
printf("traverse() node @ atom: %p\n", node->core);
printf("traverse() branches: %d\n", node->num_branches);
printf("traverse() explored: %d\n", node->num_explored);
#endif

  }
else
  if (slist)
    next = slist->data;

#if DEBUG_TRAVERSE
printf("traverse() next atom: %p\n", next);
#endif

return(next);
}

/***********************************/
enum { OPEN, CLOSED };
/***********************************/
#define DEBUG_SUBMIT_LINK 0
void submit_link(struct node_pak *rnode, struct node_pak *cnode, 
                  GList *chain, struct model_pak *data)
{
gint type; /* only really important to flag cyclic/noncyclic ie 4 drawing */
struct link_pak *dest_link=NULL, *temp_link=NULL;
GList *list;
GList *dest_list;
struct core_pak *core;

g_return_if_fail(cnode != NULL);

type = OPEN;

if (rnode)
  {
  rnode->num_explored++;
/* get atom before the repeat node */
  list = chain;
  list = g_list_last(list);
  list = g_list_previous(list);
  g_assert(list != NULL);
  core = list->data;

/* CURRENT - add to node's explored list */
  rnode->exit_list = g_slist_append(rnode->exit_list, core);

/* one node in cycle case... */
  if (cnode->core == rnode->core)
    {
    type = CLOSED;
#if DEBUG_SUBMIT_LINK
printf("terminating cyclic group.\n");
#endif
    }
/* two node in cyclic group case... */
/* TODO - general case of multiple nodes */
  list = cnode->link_list;
  while(list)
    {
    temp_link = list->data;
 
/* TODO - other cases where only one of these match - eg type = JOINED */
/* loop: A-C-B-D */
/* A-C-B + A-D-B case */
    if (temp_link->first == cnode->core && temp_link->last == rnode->core)
      {
#if DEBUG_SUBMIT_LINK
printf("terminating cyclic group.\n");
#endif
      dest_link = temp_link;
      type = CLOSED;
/* reverse input chain to make the sequence continuous */
      chain = g_list_reverse(chain);
      break;
      }
/* loop: A-C-B-D */
/* A-C-B + B-D-A case */
    if (temp_link->first == rnode->core && temp_link->last == cnode->core)
      {
/*if this happens - same as above, except chain won't need to be reversed */
      g_assert_not_reached();
      }

    list = g_list_next(list);
    }
  }

/* NB: remember the supplied chain will be destroyed by the calling routine */
/* so we must duplicate, rather than copy the list of atoms */

if (dest_link)
  {
/* modify existing link to indicate that it has common first & last */
/* nodes with another link - ie it's a cyclic link */
  dest_link->type = type;
/* NB: for cyclic chains - first & last become entry & exit point */
/* for the cyclic group (used for drawing?) */
  list = chain;

/* skip 1st element of input chain - as we only want a repeat of the ends */
/* not a repeat of the two in the middle */
  list = g_list_next(list);

/* append the supplied chain to ALREADY EXISTING link list */
  dest_list = dest_link->chain;
/* loop through the input chain & store at destination */
  while (list)
    {
    core = list->data;
    dest_list = g_list_append(dest_list, core);
    dest_link->size++;
    list = g_list_next(list);
    }
  }
else
  {
/* create new link */
  dest_link = create_link(chain);
  dest_link->type = type;
/* input chain */
  list = chain;
/* store the completed chain as a link */
  cnode->link_list = g_list_append(cnode->link_list, dest_link);
  }
}

/*********************************************/
/* hides a carbon and any attached hydrogens */
/*********************************************/
void hide_carbon(struct core_pak *core)
{
GSList *list;
struct bond_pak *bond;

g_return_if_fail(core->atom_code == 6);

core->status |= HIDDEN;

for (list=core->bonds ; list ; list=g_slist_next(list))
  {
  bond = (struct bond_pak *) list->data;

  if ((bond->atom1)->atom_code == 1)
    (bond->atom1)->status |= HIDDEN;
  if ((bond->atom2)->atom_code == 1)
    (bond->atom2)->status |= HIDDEN;
  }
}

/*********************************************/
/* eliminate (prune) all dead end type links */
/*********************************************/
#define DEBUG_PREPROCESS 0
void preprocess(struct core_pak *start, struct model_pak *data)
{
gint flag;
GSList *list;
struct core_pak *core;

g_assert (data != NULL);

/* set global - yuck */
if (start)
  backbone = start->atom_code;
else
  backbone = 0;

/* clear the pruned flag for backbone & prune all others */
for (list=data->cores ; list ; list=g_slist_next(list))
  {
  core = list->data;

  if (backbone)
    {
  if (core->atom_code == backbone)
    core->status &= ~PRUNED;
  else
    core->status |= PRUNED;
    }
  else
    core->status &= ~PRUNED;
  }

/* prune "dangling" backbone type links */
if (backbone)
  {
flag=1;
while (flag)
  {
  flag=0;
  for (list=data->cores ; list ; list=g_slist_next(list))
    {
    core = list->data;

/* only prune once */
    if (core->status & PRUNED)
      continue;

/* prune off dead-ends */
    if (total_branches(core) < 2)
      {
      core->status |= PRUNED;
      flag++;
      }
/* hack's for nasty disordered structures */
/* hack 1 - discard bonds if more than 4 */
    if (total_branches(core) > 4)
      {
      core->status |= PRUNED;
      flag++;
      }
/* hack 2 - prune off partial occupancies */
    if (core->sof < 1.0)
      {
      core->status |= PRUNED;
      flag++;
      }
    }
  }
  }

#if DEBUG_PREPROCESS
for (list=data->cores ; list ; list=g_slist_next(list))
  {
  core = list->data;
  if (core->status & PRUNED)
    core->status |= HIDDEN;
  }
#endif
}

/****************************************/
/* select connected atom types for the  */
/* carbon backbone/ribbon drawing stuff */
/****************************************/
#define DEBUG_CONSTRUCT_BACKBONE 0
void construct_backbone(struct core_pak *start, struct model_pak *data)
{
gint match, num_exits;
GList *list, *list2, *node_list=NULL, *link_list;
GList *current_chain=NULL;
GSList *tmp_list;
struct node_pak *current_node, *repeat_node, *temp_node;
struct link_pak *current_link;
struct core_pak *core, *path;

/* checks */
g_assert(start != NULL);
g_assert(data != NULL);

/* clicked atom */
core = start;

/* prune the model first (ie remove dangling chains) */
/* to make tracing the connectivity easier */
preprocess(start, data);

/* attempt to set up initial node */
switch(total_branches(core))
  {
/* no pathways out of this atom - we're done already */
  case 0:
    current_node = NULL;
    break;

/* should have already been pruned */
  case 1:
#if DEBUG_CONSTRUCT_BACKBONE
printf("Isolated chain.\n");
#endif
    return;
   
/* if first atom is not a node - traverse until find one */
  case 2:
#if DEBUG_CONSTRUCT_BACKBONE
printf("atom: %p is not a node.\n", core);
#endif
/* get next atom & append to current chain */
    tmp_list = exit_branches(NULL, core);

if (!tmp_list)
  return;

    path = tmp_list->data; 
    g_slist_free(tmp_list);

    while (total_branches(path) < 3)
      {
      tmp_list = exit_branches(core, path);
      core = path;
      path = (struct core_pak *) tmp_list->data; 
      g_slist_free(tmp_list);

/* failsafe */
if (path == start)
  {
#if DEBUG_CONSTRUCT_BACKBONE
printf("Lone cyclic group.\n");
#endif
  return;
  }

      }
#if DEBUG_CONSTRUCT_BACKBONE
printf("starting at atom: %p\n", path);
#endif
    core = path;

/* at a node */
  default:
    current_node = create_node(core);
    node_list = g_list_append(node_list, current_node);
    current_chain = g_list_append(current_chain, core);

    hide_carbon(core);
  }

/* init linkages */
create_link(NULL);
list=list2=link_list=NULL;

/* while we have a node to explore... */
while (current_node)
  {
/* get next atom & append to current chain */
  core = traverse(current_node, current_chain);
  g_assert(core != NULL);

/* NEW - hide all atoms traversed */
/* FIXME - make more sophisticated? eg hide "danglers" as well */
  hide_carbon(core);

  current_chain = g_list_append(current_chain, core);

/* NB: num_exits is one less than number of bonds - since we ignore the */
/* bond that we came in by when traversing the connectivity tree */
  num_exits = total_branches(core) - 1;

/* FIXME - hack - the one possiblity we've missed: */
/* if current atom is a node BUT, only one exit branch */
/* ie it's the first atom clicked - which is automatically */
/* made a node */
/* make the switch statement trigger the node case */
  if (core == current_node->core)
    num_exits = 2;

  switch (num_exits)
    {
/* pre-process should remove all such possibilities */
    case 0:
      g_assert_not_reached();
      break;

/* only 1 possible exit - continue to traverse the current chain */
    case 1:
      break;

/* >= 2 : this atom is a node */
    default:
/* is it a pre-existing node? */
      match=0;

#if GOK
      repeat_node = NULL;
      list = node_list;
      while (list)
        {
        temp_node = list->data;
        if (core == temp_node->core)
          {
/* repeated node - join incomplete chains of the current node */
/* NB: this is the part that will fail on more complex chain structures */
/* such as a cyclic structure that has more than two nodes */
/* TODO - maybe repeat_node could be replaced by a list of common nodes */
/* that build up the cyclic structure */
          repeat_node = temp_node;
          break;
          }
        list = g_list_next(list);
        }
#endif

repeat_node = get_node_by_atom(core, node_list);

/* store the current chain (repeat_node is NULL if node is new) */
/* this routine is the core of the whole thing, it should match */
/* up uncompleted chains to form cyclic chains, update number */
/* of explored pathways etc. etc. */
/* NB: the current chain must be duplicated before it's stored, */
/* since this routine will re-used the pointer for the next */
/* atom traversing loop */
      submit_link(repeat_node, current_node, current_chain, data);

/* if the node is not a repeat (ie new) then create */
      if (!repeat_node)
        {
        temp_node = create_node(core);

/* add the path we came in by */
        list = g_list_last(current_chain);
        list = g_list_previous(list);
        g_assert(list != NULL);

        path = list->data;

        temp_node->exit_list = g_slist_append(temp_node->exit_list, path);
        temp_node->num_explored++;

/* add to main list */
        node_list = g_list_append(node_list, temp_node);
        }

/* scan previous nodes for unexplored pathways  */
      current_node = NULL;
      list = node_list;
      while (list)
        {
        current_node = list->data;
        if (current_node->num_explored < current_node->num_branches)
          break;
        list = g_list_next(list);
        current_node = NULL;
        }

/* init for the next traverse call */
      if (current_node)
        {
#if DEBUG_CONSTRUCT_BACKBONE
printf("init traverse with node [%p] (%d/%d)\n", current_node->core,
            current_node->num_explored, current_node->num_branches);
#endif

        core = current_node->core;
        }

/* new pathway => new chain to explore */
      g_list_free(current_chain);
      current_chain = NULL;
      current_chain = g_list_append(current_chain, core);
      break;
   
    }
  }

#if DEBUG_CONSTRUCT_BACKBONE
printf("node list traversed.\n");
printf("-------------------------\n");
#endif

/* summary (stop compiler warnings) */
current_link=NULL;
#if DEBUG_CONSTRUCT_BACKBONE
list = node_list;
while (list)
  {
  current_node = list->data;
  printf("node: %p, branches: %d, explored: %d\n", 
  current_node->core, current_node->num_branches, current_node->num_explored);
  link_list = current_node->link_list; 
  while (link_list)
    {
    current_link = link_list->data;
    printf("link type: %d, size: %d\n", current_link->type, current_link->size);
    printf("\t");
    list2 = current_link->chain;
    while (list2)
      {
      printf("[%p] ", list2->data); 
      list2 = g_list_next(list2);
      }
    printf("\n");
    link_list = g_list_next(link_list);
    }
  printf("-------------------------\n");
  list = g_list_next(list);
  }
#endif

/* cleanup */
g_list_free(current_chain);
for (tmp_list=data->cores ; tmp_list ; tmp_list=g_slist_next(tmp_list))
  {
  core = tmp_list->data;
  core->status &= ~PRUNED;
  }

/* use cyclic/single linkages to construct ribbon */
construct_ribbon(node_list, data);
backbone = 0;

/* recalc coords (ie the ribbon coords) & draw */
coords_init(REDO_COORDS, data);
redraw_canvas(SINGLE);
}

/******************************************************/
/* general function for initializing a special object */
/******************************************************/
/* deprec by the new spatial code */
struct object_pak *create_object(gint type, struct model_pak *data)
{
struct object_pak *obj_data;

obj_data = g_malloc(sizeof(struct object_pak));
obj_data->id = 0;       /* don't bother with this? */
obj_data->type = type;
obj_data->data = NULL;

switch (type)
  {
  case RIBBON:
    break;

  default:
    g_assert_not_reached();
  }

/* is this quicker than an append??? */
data->ribbons = g_slist_reverse(data->ribbons);
data->ribbons = g_slist_prepend(data->ribbons, (gpointer *) obj_data);
data->ribbons = g_slist_reverse(data->ribbons);

return(obj_data);
}

/***************************************/
/* use node data to construct a ribbon */
/***************************************/
#define DEBUG_CONSTRUCT_RIBBON 0
void construct_ribbon(GList *node_list, struct model_pak *data)
{
gdouble vec[3], norm[4], *v1, *v2;
GList *list, *list2, *list3, *list4, *link_list;
GSList *ribbon, *rlist1, *rlist2;
struct link_pak *current_link, *link2, *group1, *group2;
struct node_pak *current_node, *node2;
struct ribbon_pak *ribbon1, *ribbon2;
struct object_pak *object_data;
struct core_pak *core, *join1, *join2;

/* TODO - scan the SINGLE links & check if they reference */
/* two cyclic groups at their end points (this should always */
/* happen for cyclic groups) & if so - submit a ribbon segment */
/* comprising the cyclic group centroids/normals for the two */
/* connected cyclic groups */
/* drawing - 2 stub's for each ribbon segment (one for each planar group) */
/* and an evaluator ribbon connecting the two stubs */

object_data = create_object(RIBBON, data);
ribbon = NULL;

list = node_list;
while (list)
  {
  current_node = list->data;

  link_list = current_node->link_list; 
  while (link_list)
    {
    current_link = link_list->data;

/* find open chains */
    if (current_link->type == OPEN)
      {
        list2 = node_list;
        group1=group2=NULL;
        join1=join2=NULL;
        while (list2)
          {
          node2 = list2->data;
          list3 = node2->link_list; 
          while (list3)
            {
            link2 = list3->data;

/* match only cyclic groups */
            if (link2->type == CLOSED)
              {
              list4 = link2->chain;
              
              while (list4)
                {
                core = list4->data;
                if (current_link->first == core)
                  {
                  group1 = link2;
                  join1 = core;
                  }
                if (current_link->last == core)
                  {
                  group2 = link2;
                  join2 = core;
                  }
                list4 = g_list_next(list4);
                }
              }
            list3 = g_list_next(list3);
            }
          list2 = g_list_next(list2);
          }

/* two found - record a ribbon segment */
/* comprises two points + normals at each point */
/* probably need one extra vector - orientation (eg rectangle) */
/* maybe some intermediate points in the joining chain? */

      if (group1 && group2)
        {
#if DEBUG_CONSTRUCT_RIBBON
printf("Two cyclic groups connected via: [%p] [%p], sizes: (%d) (%d)\n",
       join1, join2, group1->size, group2->size);
#endif

ribbon1 = g_malloc(sizeof(struct ribbon_pak));

ARR3SET(ribbon1->colour, sysenv.render.ribbon_colour);
ribbon1->colour[3] = sysenv.render.transmit;

/* append the ribbon segment */
ribbon = g_slist_append(ribbon, (gpointer *) ribbon1);

/* save link ids */
ribbon1->id1 = group1->id;
ribbon1->id2 = group2->id;

/* compute centroid of group 1 */
        VEC3SET(vec, 0.0, 0.0, 0.0);
        list2 = group1->chain; 
/* exclude the first point - it's the same as the first */
        list2 = g_list_next(list2);
        while (list2)
          {
          core = list2->data;
          ARR3ADD(vec, core->x);
          list2 = g_list_next(list2);
          }
        VEC3MUL(vec, 1.0/(group1->size-1));
        ARR3SET(ribbon1->x1, vec);

/* compute centroid of group 1 */
        VEC3SET(vec, 0.0, 0.0, 0.0);
        list2 = group2->chain; 
/* exclude the first point - it's the same as the first */
        list2 = g_list_next(list2);
        while (list2)
          {
          core = list2->data;
          ARR3ADD(vec, core->x);
          list2 = g_list_next(list2);
          }
        VEC3MUL(vec, 1.0/(group2->size-1));
        ARR3SET(ribbon1->x2, vec);

/* compute normal of each */
/* NB: result is a 4D vector Ax+By+Cz+D = 0, even tho we only want the normal */
        compute_plane(norm, group1->chain); 
        normalize(norm, 3);
        ARR3SET(ribbon1->u1, norm);

        compute_plane(norm, group2->chain); 
	normalize(norm, 3);
        ARR3SET(ribbon1->u2, norm);

/* compute orientation 1->2 */
        ARR3SET(vec, ribbon1->x2);
        ARR3SUB(vec, ribbon1->x1);
        normalize(vec, 3);
        proj_vop(ribbon1->v1, vec, ribbon1->u1);

/* compute orientation 2->1 */
        ARR3SET(vec, ribbon1->x1);
        ARR3SUB(vec, ribbon1->x2);
        normalize(vec, 3);
        proj_vop(ribbon1->v2, vec, ribbon1->u2);
        }
      }
    link_list = g_list_next(link_list);
    }
  list = g_list_next(list);
  }



/* NEW - post process - ribbon segment smoothing */
rlist1 = ribbon;
while (rlist1)
  {
  ribbon1 = rlist1->data;

#if DEBUG_CONSTRUCT_RIBBON
printf("ribbon segment: [%d] - [%d]\n", ribbon1->id1, ribbon1->id2);
#endif

/* "consistent" normal enforcement */
/* of dubious value, since ribbon lighing will be two sided */
/* also, not as easy as this (eg twisted ribbon) ... good enough though? */
  if (via(ribbon1->x1, ribbon1->u1, 3) > PI/2.0)
    {
    VEC3MUL(ribbon1->u1, -1.0);
    }
  if (via(ribbon1->x2, ribbon1->u2, 3) > PI/2.0)
    {
    VEC3MUL(ribbon1->u2, -1.0);
    }

/* search uniquely for a segment with common link id's */
  rlist2 = g_slist_next(rlist1);
  while (rlist2)
    {
    ribbon2 = rlist2->data;

/* if match - get the two orientation vectors at the ribbon meeting point */
    v1 = v2 = NULL;
    if (ribbon1->id1 == ribbon2->id1)
      {
      v1 = ribbon1->v1;
      v2 = ribbon2->v1;
      }
    if (ribbon1->id1 == ribbon2->id2)
      {
      v1 = ribbon1->v1;
      v2 = ribbon2->v2;
      }
    if (ribbon1->id2 == ribbon2->id1)
      {
      v1 = ribbon1->v2;
      v2 = ribbon2->v1;
      }
    if (ribbon1->id2 == ribbon2->id2)
      {
      v1 = ribbon1->v2;
      v2 = ribbon2->v2;
      }

/* if match - align */
    if (v1 && v2)
      {
/* subtract the normalized orientation vectors */
      normalize(v1, 3);
      normalize(v2, 3);
      ARR3SET(vec, v1);
      ARR3SUB(vec, v2);

/* assign according to subtraction order, so that v1 = -v2 */
      normalize(vec, 3);

      ARR3SET(v1, vec);
      ARR3SET(v2, vec);
      VEC3MUL(v2, -1.0);
      }
    rlist2 = g_slist_next(rlist2);
    }
  rlist1 = g_slist_next(rlist1);
  }

object_data->data = (gpointer *) ribbon;
}


