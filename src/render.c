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
#include <math.h>

#include "gdis.h"
#include "coords.h"
#include "matrix.h"
#include "opengl.h"
#include "render.h"
#include "interface.h"

/* externals */
extern struct sysenv_pak sysenv;
extern GtkWidget *window;
extern GdkPixmap *pixmap;

/************************************************************/
/* process data structures to get current model state label */
/************************************************************/
gchar *get_mode_label(struct model_pak *data)
{
gchar *label;

switch(data->mode)
  {
  case FREE:
    label = g_strdup("normal");
    break;
  case ATOM_ADD:
    label = g_strdup("add atoms");
    break;
  case BOND_DELETE:
    label = g_strdup("delete bonds");
    break;
  case BOND_SINGLE:
    label = g_strdup_printf("add single bonds (%1d/2)", data->state);
    break;
  case BOND_DOUBLE:
    label = g_strdup_printf("add double bonds (%1d/2)", data->state);
    break;
  case BOND_TRIPLE:
    label = g_strdup_printf("add triple bonds (%1d/2)", data->state);
    break;
  case BOND_INFO:
    label = g_strdup("bonds");
    break;
  case DIST_INFO:
    label = g_strdup_printf("distance selection (%1d/2)", data->state);
    break;
  case ANGLE_INFO:
    label = g_strdup_printf("angle selection (%1d/3)", data->state);
    break;
  case DIHEDRAL_INFO:
    label = g_strdup_printf("dihedral selection (%1d/4)", data->state);
    break;
  case DEFINE_RIBBON:
    label = g_strdup("define ribbon");
    break;
  case DEFINE_VECTOR:
    label = g_strdup_printf("define vector (%1d/2)", data->state);
    break;
  case DEFINE_PLANE:
    label = g_strdup_printf("define plane (%1d/3)", data->state);
    break;
  case ANIMATING:
    label = g_strdup_printf("Frame: %-4d", data->cur_frame);
    break;
  case RECORD:
    label = g_strdup("recording");
    break;
  case SELECT_FRAGMENT:
    label = g_strdup_printf("select bonded atoms (%1d/2)", data->state);
    break;

  default:
    label = g_strdup("undefined");
    break;
  }

return(label);
}

/*********************************************/
/* construct all pipe lists for bond drawing */
/*********************************************/
#define PIPE_DEPTH 4
void render_make_pipes(GSList **pipes, struct model_pak *model)
{
gint i;
gdouble radius, mp1[4], mp2[4], x[3];
gdouble colour1[4], colour2[4];
GSList *list;
struct pipe_pak *pipe;
struct bond_pak *bond;
struct core_pak *core1, *core2;

/* checks */
g_assert(model != NULL);

/*
printf("-------------------------------------------------------------\n");
dump_bonds(model);
printf("-------------------------------------------------------------\n");
*/

/* init the return list(s) */
for (i=PIPE_DEPTH ; i-- ; )
  pipes[i] = NULL;

/* bond display turned off? */
if (!model->show_bonds)
  return;

/* common pipe width */
radius = sysenv.render.stick_rad;

/* enumerate bonds to construct the pipe list */
for (list=model->bonds; list ; list=g_slist_next(list))
  {
  bond = list->data;

  if (bond->status == DELETED)
    continue;
  if (bond->status == HIDDEN)
    continue;

/* the two atoms */
  core1 = bond->atom1;
  core2 = bond->atom2;

/* NEW - bailout modes */
  if (core1->render_mode == ZONE)
    continue;
  if (core2->render_mode == ZONE)
    continue;

  if (core1->status & OFF_SCREEN)
    if (core2->status & OFF_SCREEN)
      continue;

/* NEW - cope with phonon animation moving the cores */
  ARR3SET(x, core1->offset);
  ARR3ADD(x, core2->offset);
  VEC3MUL(x, 0.5);

/* compute midpoint 1 coords */
  ARR3SET(mp1, core1->x);
  ARR3ADD(mp1, bond->offset);
  mp1[3] = 1.0;
  ARR3ADD(mp1, x);
  vec4mat(model->display_lattice, mp1);

/* compute midpoint 2 coords */
/* midpoints may be different (split by pbc) */
  ARR3SET(mp2, core2->x);
  ARR3SUB(mp2, bond->offset);
  mp2[3] = 1.0;
  ARR3ADD(mp2, x);
  vec4mat(model->display_lattice, mp2);

/* deleted/ zeol hidden */
  if (core1->status & (DELETED | ZEOL_HIDDEN))
    continue;
  if (core2->status & (DELETED | ZEOL_HIDDEN))
    continue;

/* colour setup */
  switch (bond->type)
    {
    case BOND_HBOND:
      VEC4SET(colour1, 1.0, 1.0, 0.6, 0.0);
      VEC4SET(colour2, 1.0, 1.0, 0.6, 0.0);
      break;

    default:
      ARR4SET(colour1, core1->colour);
      VEC3MUL(colour1, 1.0/65535.0);
      ARR4SET(colour2, core2->colour);
      VEC3MUL(colour2, 1.0/65535.0);
    }

/* setup half-bond (pipe) for core1 */
  if (!(core1->status & HIDDEN) && core1->render_mode != CPK)
    {
/* init pipe */
    pipe = g_malloc(sizeof(struct pipe_pak));
    ARR3SET(pipe->v1, core1->rx);
    ARR3SET(pipe->v2, mp1);
    pipe->radius = radius;
    ARR4SET(pipe->colour, colour1);

/* assign to appropriate pipe list */
    if (core1->render_mode == STICK)
      pipes[3] = g_slist_prepend(pipes[3], pipe);
    else
      {
      if (core1->ghost)
        pipes[1] = g_slist_prepend(pipes[1], pipe);
      else
        {
        if (core1->render_wire)
          pipes[2] = g_slist_prepend(pipes[2], pipe);
        else
          pipes[0] = g_slist_prepend(pipes[0], pipe);
        }
      }
    }

/* setup half-bond (pipe) for core1 */
  if (!(core2->status & HIDDEN) && core2->render_mode != CPK)
    {
/* init pipe */
    pipe = g_malloc(sizeof(struct pipe_pak));
    ARR3SET(pipe->v1, core2->rx);
    ARR3SET(pipe->v2, mp2);
    pipe->radius = radius;
    ARR4SET(pipe->colour, colour2);

/* assign to appropriate pipe list */
    if (core2->render_mode == STICK)
      pipes[3] = g_slist_prepend(pipes[3], pipe);
    else
      {
      if (core2->ghost)
        pipes[1] = g_slist_prepend(pipes[1], pipe);
      else
        {
        if (core2->render_wire)
          pipes[2] = g_slist_prepend(pipes[2], pipe);
        else
          pipes[0] = g_slist_prepend(pipes[0], pipe);
        }
      }
    }
  }
}

/*****************************/
/* pipe z-ordering primitive */ 
/*****************************/
gint render_pipe_depth_sort(struct pipe_pak *p1, struct pipe_pak *p2)
{
gdouble z1, z2;

/* use pipe z-midpoints for comparison */
z1 = 0.5 * (p1->v1[2] + p1->v2[2]);
z2 = 0.5 * (p2->v1[2] + p2->v2[2]);

if (z1 > z2)
  return(-1);
return(1);
}

GSList *render_sort_pipes(GSList *pipes)
{
return(g_slist_sort(pipes, (gpointer) render_pipe_depth_sort));
}
