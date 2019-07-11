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
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <unistd.h>
#include <gtk/gtk.h>
#include <gtk/gtkgl.h>
#include <gdk/gdk.h>
#include <gdk/gdkkeysyms.h>

#include "gdis.h"
#include "coords.h"
#include "edit.h"
#include "file.h"
#include "graph.h"
#include "task.h"
#include "morph.h"
#include "model.h"
#include "module.h"
#include "matrix.h"
#include "render.h"
#include "select.h"
#include "space.h"
#include "sginfo.h"
#include "spatial.h"
#include "opengl.h"
#include "quaternion.h"
#include "surface.h"
#include "gui_shorts.h"
#include "interface.h"
#include "dialog.h"
#include "zmatrix.h"
#include "gui_image.h"
#include "undo.h"
#include "track.h"

#include "logo_left.xpm"
#include "logo_right.xpm"

/* top level data structure */
extern struct sysenv_pak sysenv;

extern GMutex *gdk_threads_mutex;

/* bit messy - put some of these in sysenv? */
/* main window */
GtkWidget *window;
/* backing pixmap */
GdkPixmap *pixmap = NULL;
/* model tree (put in sysenv?) */
GtkWidget *tree;
/* model pane stuff */
GtkWidget *scrolled_window, *clist;
/* gdisrc? */
gint pane_width=65;
/* current viewing mode combo box */
GtkWidget *viewing_mode;

/************************/
/* EVENT - button press */
/************************/
#define DEBUG_BUTTON_PRESS_EVENT 0
gint gui_press_event(GtkWidget *w, GdkEventButton *event)
{
gint refresh=0, x, y;
gint shift=FALSE;
#ifdef UNUSED_BUT_SET
gint ctrl=FALSE;
#endif
GdkModifierType state;
struct model_pak *data;
struct bond_pak *bond;
struct core_pak *core;

/* HACK */
if (sysenv.stereo)
  return(FALSE);

/* get model */
data = sysenv.active_model;
if (!data)
  return(FALSE);

/* event info */
x = event->x;
y = event->y;
sysenv.moving = FALSE;

canvas_select(x, y);

/* analyse the current state */
state = (GdkModifierType) event->state;
if ((state & GDK_SHIFT_MASK))
  shift = TRUE;
#ifdef UNUSED_BUT_SET
if ((state & GDK_CONTROL_MASK))
  ctrl = TRUE;
#endif

/* only want button 1 (for now) */
if (event->button != 1)
  return(FALSE);

/* NEW - diffraction peak search */
if (data->graph_active)
  {/*FIXME: why only diffract?*/
  struct graph_pak *graph=(struct graph_pak *)data->graph_active;
  switch(graph->type){
  case GRAPH_REGULAR:
    diffract_select_peak(x, y, data);
    break;
  case GRAPH_IY_TYPE:
  case GRAPH_XY_TYPE:
  case GRAPH_YX_TYPE:
  case GRAPH_IX_TYPE:
  case GRAPH_XX_TYPE:
    dat_graph_select(x,y,data);
    break;
  default:
    break;
  }
  return(FALSE);
  }

/* can we assoc. with a single atom */
core = gl_seek_core(x, y, data);

/* allow shift+click to add/remove single atoms in the selection */
/* TODO - depending on mode add/remove objects in selection eg atom/mols */
if (shift)
  {
  switch(data->mode)
    {
    default:
      if (core)
        {
        select_core(core, TRUE, data);
        refresh++;
        }
      else
        {
/* otherwise start a box selection */
        update_box(x,y,data,START);
        }
      break;
    }
  }
else
  {
/* determine the type of action required */
  switch(data->mode)
    {
    case BOND_INFO:
      info_bond(w,x,y,data);
      break;

    case BOND_DELETE:
/* or by bond midpoints */
      bond = gl_seek_bond(x,y,data);
      if (bond)
        {
        connect_user_bond(bond->atom1, bond->atom2, BOND_DELETE, data);
        refresh++;
        }
      break;

    case BOND_SINGLE:
    case BOND_DOUBLE:
    case BOND_TRIPLE:
      if (core)
        {
        connect_make_bond(core, data->mode, data); 
        refresh++;
        }
      break;

    case ATOM_ADD:
      add_atom(x, y, data);
      refresh++; 
      break;

    case DEFINE_RIBBON:
      if (core)
        construct_backbone(core, data);
      break;

    case DEFINE_VECTOR:
    case DEFINE_PLANE:
      if (core)
        spatial_point_add(core, data);
      break;

/* selection stuff */
    default:
      select_clear(data);
/* don't select if a core is under the mouse */
/* instead we allow the selection box to be drawn */
      update_box(x,y,data,START);
      refresh++;
      break; 
    }
  }

/* CURRENT */
if (refresh)
  gui_active_refresh();

return(FALSE);
}

/***************************/
/* EVENT - button released */
/***************************/
gint gui_release_event(GtkWidget *w, GdkEventButton *event)
{
gint x, y;
struct model_pak *data;
#ifdef UNUSED_BUT_SET
GdkModifierType state;
#endif

/* get model */
data = sysenv.active_model;
if (!data)
  return(FALSE);

/* get event info */
x = event->x;
y = event->y;
#ifdef UNUSED_BUT_SET
state = (GdkModifierType) event->state;
#endif

/* NEW - motion flag */
sysenv.moving = FALSE;

/* first mouse button */
switch (event->button)
  {
  case 1:

/* HACK */
if (sysenv.stereo)
  return(FALSE);

/* clean up after move operations */
  switch(data->mode)
    {
    case DIST_INFO:
      info_dist(x,y,data);
      break;
    case ANGLE_INFO:
      info_angle(x,y,data);
      break;
    case DIHEDRAL_INFO:
      info_torsion(x,y,data);
      break;

    default:
/* commit the box contents to the selection */
/* if shift pressed - append selection, otherwise create new */
      if (data->box_on)
        gl_select_box(w);
/* turn the box off */
      data->box_on = FALSE;
      break;
    }
  break;

  case 2:
/* hack to only update zoom factor when button released */
/* ie continuous update causes significant slowdown */
  gui_camera_refresh();
  break;
  }

redraw_canvas(SINGLE);
return(FALSE);
}

/************************/
/* EVENT - mouse scroll */
/************************/
gint gui_scroll_event(GtkWidget *w, GdkEventScroll *event)
{
/* change zoom -- based on "zoom section" of gui_press_event() */
const gdouble scroll_factor = 10.0;
gdouble scroll, v[3];
struct model_pak *data;
struct camera_pak *camera;

/* get model */
data = sysenv.active_model;
if (!data)
  return(FALSE);

camera = data->camera;

sysenv.moving = FALSE;

switch (event->direction)
{
  case GDK_SCROLL_UP:
    scroll = -scroll_factor;
    break;
  case GDK_SCROLL_DOWN:
    scroll = scroll_factor;
    break;
  default: /* left and right are not used yet */
    return FALSE;
}
scroll *= PIX2SCALE;

if (camera->perspective)
  {
  ARR3SET(v, camera->v);
  VEC3MUL(v, -scroll*10.0);
  ARR3ADD(camera->x, v);
  }
else
  camera->zoom += scroll;

data->zoom = data->rmax;

gui_camera_refresh();

sysenv.moving = TRUE;
redraw_canvas(SINGLE);
return FALSE;
}

/************************/
/* EVENT - mouse motion */
/************************/
#define DEBUG_MOTION 0
gint gui_motion_event(GtkWidget *w, GdkEventMotion *event)
{
gint x, y, dx, dy, fx, fy, refresh;
gint shift=FALSE, ctrl=FALSE;
gdouble da, dv, zoom, v[3], mat[9];
GdkModifierType state;
static gint ox=0, oy=0;
struct model_pak *data;
struct camera_pak *camera;

/* get model */
data = sysenv.active_model;
if (!data)
  return(FALSE);

camera = data->camera;

/* get mouse data */
if (event->is_hint)
  gdk_window_get_pointer(event->window, &x, &y, &state);
else
  {
/* MS windows apparently reaches this */
  x = event->x;
  y = event->y;
  state = (GdkModifierType) event->state;
  }

/* discard discontinuous jumps (ie ox,oy are invalid on 1st call) */
if (!sysenv.moving)
  {
  ox = x;
  oy = y;
  sysenv.moving = TRUE;
  return(TRUE);
  }

/* convert relative mouse motion to an increment */
dx = x-ox;
dy = oy-y;        /* inverted y */

/* single update */
refresh=0;

/* analyse the current state */
if ((state & GDK_SHIFT_MASK))
  shift = TRUE;
if ((state & GDK_CONTROL_MASK))
  ctrl = TRUE;

/* first mouse button - mostly selection stuff */
if (state & GDK_BUTTON1_MASK)
  {
  switch(data->mode)
    {
    case FREE:
      update_box(x, y, data, UPDATE);
      refresh++;
      break;

    default:
      break;
    }
/* don't switch to low quality drawing */
  sysenv.moving = FALSE;
  if (refresh)
    redraw_canvas(SINGLE);
  ox=x;
  oy=y;
  return(TRUE);
  }

/* NEW - disallow rotations (redraws) when graph is displayed */
if (data->graph_active)
  return(TRUE);

/* second mouse button */
if (state & GDK_BUTTON2_MASK)
  {
  if (shift)
    {
/* zoom */
    zoom = PIX2SCALE * (dx+dy);
    if (camera->perspective)
      {
      ARR3SET(v, camera->v);
      VEC3MUL(v, zoom*10.0);
      ARR3ADD(camera->x, v);
      }
    else
      {
      camera->zoom -= zoom;
      data->zoom = data->rmax;
      }
    sysenv.moving = TRUE;
    refresh++;
    }
  else
    {
    if (ctrl)
      {
/* selection only translation */
      select_translate(x-ox, y-oy, data);
      sysenv.moving = TRUE;
      refresh++;
      }
    else
      {
/* horizontal camera translation */
      dv = ox-x;
      ARR3SET(v, camera->e);
      VEC3MUL(v, dv*PIX2ANG);
      ARR3ADD(camera->x, v);
/* vertical  camera translation */
      dv = y-oy;
      ARR3SET(v, camera->o);
      VEC3MUL(v, dv*PIX2ANG);
      ARR3ADD(camera->x, v);

      sysenv.moving = TRUE;
      refresh++;
      }
    }
  }

/* third mouse button clicked? */
if (state & GDK_BUTTON3_MASK)
  {
/* shift clicked? */
  if (shift)
    {
/* yaw */
    if (dx || dy)
      {
/* rotation amount */
      da = abs(dx) + abs(dy);
/* vector from center to mouse pointer (different for OpenGL window) */
/* FIXME - if ctrl is pressed, the center should be the selection */
/* centroid & not the middle of the window */
      fx = x - data->offset[0] - (sysenv.x + sysenv.width/2);
      fy = data->offset[1] + (sysenv.y + sysenv.height/2 - y);

/* rotation direction via z component of cross product (+=clock, -=anti) */
      if ((fx*dy - dx*fy) < 0)
        da *= -1.0;

#if DEBUG_MOTION
printf("(%d,%d) x (%d,%d) : %f\n",dx,dy,fx,fy, da);
#endif

/* assign and calculate */
      da *= D2R * ROTATE_SCALE;
      if (ctrl)
        {
        matrix_relative_rotation(mat, da, ROLL, data);
        rotate_select(data, mat);
        }
      else
        {
        if (camera->mode == LOCKED)
          quat_concat_euler(camera->q, ROLL, da);
        else
          {
          matrix_v_rotation(mat, camera->v, -da);
          vecmat(mat, camera->o);
          vecmat(mat, camera->e);
          }
        }
      }
    sysenv.moving = TRUE;
    refresh++;
    }
  else
    {
#if DEBUG_MOTION
printf("(%d,%d)\n",dx,dy);
#endif
/* pitch and roll */
    if (dy)
      {
      da = D2R * ROTATE_SCALE * dy;
      if (ctrl)
        {
        matrix_relative_rotation(mat, -da, PITCH, data);
        rotate_select(data, mat);
        }
      else
        {
        if (camera->mode == LOCKED)
          quat_concat_euler(camera->q, PITCH, da);
        else
          {
          matrix_v_rotation(mat, camera->e, -da);
          vecmat(mat, camera->v);
          vecmat(mat, camera->o);
          }
        }
      sysenv.moving = TRUE;
      refresh++;
      }
    if (dx)
      {
      da = D2R * ROTATE_SCALE * dx;
      if (ctrl)
        {
        matrix_relative_rotation(mat, -da, YAW, data);
        rotate_select(data, mat);
        }
      else
        {
        if (camera->mode == LOCKED)
          quat_concat_euler(camera->q, YAW, -da);
        else
          {
          matrix_v_rotation(mat, camera->o, da);
          vecmat(mat, camera->v);
          vecmat(mat, camera->e);
          }
        }
      sysenv.moving = TRUE;
      refresh++;
      }
    }
  }

/* save old values */
ox=x;
oy=y;

/* redraw? */
if (refresh)
  redraw_canvas(SINGLE);

return(TRUE);
}

/***************************/
/* expose all hidden atoms */
/***************************/
void unhide_atoms(void)
{
GSList *list;
struct model_pak *data;
struct core_pak *core;

/* deletion for the active model only */
data = sysenv.active_model;

/* unhide */
for (list=data->cores ; list ; list=g_slist_next(list))
  {
  core = list->data;
  core->status &= ~HIDDEN;
  }

/* update */
redraw_canvas(SINGLE);
}

/****************/
/* change mode */
/***************/
void gui_mode_switch(gint new_mode)
{
GSList *list;
struct model_pak *data;
struct core_pak *core1;

/* get selected model */
data = sysenv.active_model;
if (!data)
  return;

/* special case for morphology */
if (data->id == MORPH)
  {
  switch(new_mode)
    {
    case FREE:
    case RECORD:
      break;

    default:
      gui_text_show(WARNING, "Disallowed mode.\n");
      return;
    }
  }

/* clean up (if necessary) from previous mode */
switch (data->mode)
  {
  case DIST_INFO:
  case BOND_INFO:
    for (list=data->cores ; list ; list=g_slist_next(list))
      {
      core1 = list->data;
      core1->status &= ~SELECT;
      }
    break;

  case RECORD:
    if (data->num_frames > 2)
      {
/* CURRENT - I don't know why, but we seem to get a duplicate of the 1st frame */
/* revert to original orientation */
      camera_init(data);
      data->num_frames = g_slist_length(data->transform_list);
      }
    else
      {
/* if no extra frames - it's not an animation any more */
      data->num_frames = 1;
      data->animation = FALSE;
      }
    break;
  }

/* special initialization */
switch (new_mode)
  {
  case RECORD:
/* NEW - disallow transformation record mode if file is */
/* a standard animation as this will confuse read_frame() */
    if (data->frame_list || data->gulp.trj_file)
      {
      gui_text_show(ERROR, "Disallowed mode change.\n");
      return;
      }
    data->animation = TRUE;
    break;
  }

/* general initialization */
select_clear(data);
data->mode = new_mode;
data->state = 0;
redraw_canvas(SINGLE);
}

/* Need this for the menu item factory widget */
void gtk_mode_switch(GtkWidget *w, gint new_mode)
{
gui_mode_switch(new_mode);
}

/***************/
/* change view */
/***************/
void switch_view(gint new_mode)
{
}

/************/
/* GTK hook */
/************/
void gtk_switch_view(GtkWidget *w, gint new_mode)
{
switch_view(new_mode);
}

/**************************************/
/* symmetry list store free primitive */
/**************************************/
gboolean symmetry_free(GtkTreeModel *treemodel, GtkTreePath *path, GtkTreeIter *iter, gpointer data)
{
gchar *a, *b;

gtk_tree_model_get(treemodel, iter, 0, &a, 1, &b, -1);
if (a)
  g_free(a);
if (b)
  g_free(b);
return(TRUE);
}

/*****************************/
/* redraw the content widget */
/*****************************/
void gui_content_refresh(GtkWidget *box)
{
GSList *list;
GtkTreeIter iter;
GtkCellRenderer *renderer;
GtkTreeViewColumn *column;
static GtkListStore *ls=NULL;
static GtkWidget *tv=NULL;

if (box)
  {
g_assert(ls == NULL);
g_assert(tv == NULL);

/* new tree list */
  ls = gtk_list_store_new(2, G_TYPE_STRING, G_TYPE_STRING);
  tv = gtk_tree_view_new_with_model(GTK_TREE_MODEL(ls));
  gtk_box_pack_start(GTK_BOX(box), tv, TRUE, TRUE, 0);

/* setup cell renderers */
  renderer = gtk_cell_renderer_text_new();
  column = gtk_tree_view_column_new_with_attributes(" ", renderer, "text", 0, NULL);
  gtk_tree_view_append_column(GTK_TREE_VIEW(tv), column);
  renderer = gtk_cell_renderer_text_new();
  column = gtk_tree_view_column_new_with_attributes(" ", renderer, "text", 1, NULL);
  gtk_tree_view_append_column(GTK_TREE_VIEW(tv), column);
  gtk_tree_view_set_headers_visible(GTK_TREE_VIEW(tv), FALSE);

/* currently, allow nothing to be selected */
  gtk_tree_selection_set_mode(gtk_tree_view_get_selection(GTK_TREE_VIEW(tv)),
                              GTK_SELECTION_NONE);
  }
else
  {
  struct model_pak *model = sysenv.active_model;

g_assert(ls != NULL);
g_assert(tv != NULL);

  gtk_list_store_clear(ls);

  if (!model)
    return;

/* model content */
  model_content_refresh(model);

  for (list=model->property_list ; list ; list=g_slist_next(list))
    {
    if (property_rank(list->data))
      {
      gtk_list_store_append(ls, &iter);
      gtk_list_store_set(ls, &iter, 0, property_label(list->data), -1);
      gtk_list_store_set(ls, &iter, 1, property_value(list->data), -1);
      }
    }
  }
}

/*************************/
/* module widget gloabls */
/*************************/
GtkTreeStore *module_ts=NULL;
GtkWidget *module_tv=NULL;

/******************************/
/* module invocation callback */
/******************************/
void cb_module_activate(GtkTreeView *tree_view,
                        GtkTreePath *tree_path,
                        GtkTreeViewColumn *tree_view_column,
                        gpointer data)
{
gchar *symbol;
gpointer module;
GtkTreeModel *tree_model;
GtkTreeIter iter;
struct model_pak *model;

model = sysenv.active_model;

tree_model = gtk_tree_view_get_model(GTK_TREE_VIEW(tree_view));
if (gtk_tree_model_get_iter(tree_model, &iter, tree_path))
  {
  gtk_tree_model_get(tree_model, &iter, 0, &symbol, 1, &module, -1);

  if (module)
    module_invoke(symbol, model, module);
  else
    gui_text_show(ERROR, "Please select function(s) listed below the module.\n");
  }
}

/***********************/
/* refresh module list */
/***********************/
void module_widget_redraw(void)
{
GSList *list1, *list2;
GtkTreeIter root, branch;

/* checks */
g_assert(GTK_IS_TREE_STORE(module_ts));
g_assert(GTK_IS_TREE_VIEW(module_tv));

/* populate with valid modules */
for (list1=sysenv.module_list ; list1 ; list1=g_slist_next(list1))
  {
  gtk_tree_store_append(module_ts, &root, NULL);
  gtk_tree_store_set(module_ts, &root, 0, module_label(list1->data), -1);
  gtk_tree_store_set(module_ts, &root, 1, NULL, -1);
/* populate with module symbols */
  for (list2=module_symbols(list1->data) ; list2 ; list2=g_slist_next(list2))
    {
    gtk_tree_store_append(module_ts, &branch, &root);
    gtk_tree_store_set(module_ts, &branch, 0, (gchar *) list2->data, -1);
    gtk_tree_store_set(module_ts, &branch, 1, (gpointer) list1->data, -1);
    }
  }
}

/***********************************************/
/* init the available module (plug-in) listing */
/***********************************************/
void module_widget_setup(GtkWidget *box)
{
GtkWidget *label;
GtkCellRenderer *renderer;
GtkTreeViewColumn *column;

label = gtk_label_new("Loaded modules");
gtk_box_pack_start(GTK_BOX(box), label, FALSE, FALSE, 0);

/* underlying data storage */
module_ts = gtk_tree_store_new(2, G_TYPE_STRING, G_TYPE_POINTER);

/* actual tree widget */
module_tv = gtk_tree_view_new_with_model(GTK_TREE_MODEL(module_ts));
gtk_box_pack_start(GTK_BOX(box), module_tv, TRUE, TRUE, 0);

/* setup the text rendering colum */
renderer = gtk_cell_renderer_text_new();
column = gtk_tree_view_column_new_with_attributes("a", renderer, "text", 0, NULL);
gtk_tree_view_append_column(GTK_TREE_VIEW(module_tv), column);
gtk_tree_view_set_headers_visible(GTK_TREE_VIEW(module_tv), FALSE);

/* fill in with valid modules */
module_widget_redraw();

/* setup the selection handler */
g_signal_connect(G_OBJECT(module_tv), "row-activated",
                 G_CALLBACK(cb_module_activate),
                 NULL);
}

/************************************/
/* pick a (loaded) model to display */
/************************************/
#define DEBUG_PICK_MODEL 0
void gui_model_select(struct model_pak *model)
{
/* checks */
if (!model)
  return;

/* make the model active */
sysenv.active_model = model;
model->graph_active = NULL;

/* TODO - combine all these into a dialog refresh callback */
gui_relation_update(model);

/* this update could replace the auto version of the relation shortcuts */
dialog_refresh_all();

canvas_shuffle();

gui_refresh(GUI_MODEL_PROPERTIES);
gui_refresh(GUI_CANVAS);
}

/************************************************/
/* reset the scale and viewing angle to default */
/************************************************/
void gui_view_default(void)
{
struct model_pak *model;

/* do we have a loaded model? */
model = sysenv.active_model;
if (!model)
  return;

/* NEW - recentre (eg if atoms added -> centroid change) */
coords_init(CENT_COORDS, model);

/* make the changes */
VEC3SET(model->offset, 0.0, 0.0, 0.0);

/* update */
camera_init(model);

/* trigger any dialog updates */
gui_model_select(model);
}

/************************/
/* view down the x axis */
/************************/
void gui_view_x(void)
{
struct camera_pak *camera;
struct model_pak *model = sysenv.active_model;

if (model)
  {
  camera_init(model);
  camera = model->camera;
  quat_concat_euler(camera->q, YAW, 0.5*G_PI);
  gui_model_select(model);
  }
}

/************************/
/* view down the y axis */
/************************/
void gui_view_y(void)
{
struct camera_pak *camera;
struct model_pak *model = sysenv.active_model;

if (model)
  {
  camera_init(model);
  camera = model->camera;
  quat_concat_euler(camera->q, YAW, G_PI);
  gui_model_select(model);
  }
}

/************************/
/* view down the z axis */
/************************/
void gui_view_z(void)
{
struct camera_pak *camera;
struct model_pak *model = sysenv.active_model;

if (model)
  {
  camera_init(model);
  camera = model->camera;
  quat_concat_euler(camera->q, PITCH, -0.5*G_PI);
  gui_model_select(model);
  }
}

/************************/
/* view down the a axis */
/************************/
void gui_view_a(void)
{
gui_view_x();
}

/************************/
/* view down the b axis */
/************************/
void gui_view_b(void)
{
struct camera_pak *camera;
struct model_pak *model = sysenv.active_model;

if (model)
  {
  if (model->periodic > 1 || model->id == MORPH)
    {
    gui_view_x();
    camera = model->camera;
    quat_concat_euler(camera->q, YAW, model->pbc[5]);
    gui_model_select(model);
    }
  else
    gui_view_y();
  }
}

/************************/
/* view down the c axis */
/************************/
void gui_view_c(void)
{
gdouble a, c[3], v[3];
struct camera_pak *camera;
struct model_pak *model = sysenv.active_model;

if (model)
  {
  if (model->periodic > 2 || model->id == MORPH)
    {
    camera_init(model);
    camera = model->camera;

/* c axis vector */
    VEC3SET(c, 0.0, 0.0, -1.0);
    vecmat(model->latmat, c);

/* angle */
    a = via(camera->v,c,3);

/* rotation axis */
    crossprod(v, camera->v, c);
    normalize(v, 3);

/* align */
    quat_concat(camera->q, v, a);

    gui_model_select(model);
    }
  else
    gui_view_z();
  }
}

/**********************/
/* viewing widget box */
/**********************/
/* TODO - add camera properties? */
void gui_view_widget(GtkWidget *w)
{
GtkWidget *hbox;

hbox = gtk_hbox_new(TRUE, PANEL_SPACING);
gtk_box_pack_start(GTK_BOX(w), hbox, FALSE, FALSE, 0);
gui_button(" x ", gui_view_x, NULL, hbox, TT);
gui_button(" y ", gui_view_y, NULL, hbox, TT);
gui_button(" z ", gui_view_z, NULL, hbox, TT);
gui_button(" a ", gui_view_a, NULL, hbox, TT);
gui_button(" b ", gui_view_b, NULL, hbox, TT);
gui_button(" c ", gui_view_c, NULL, hbox, TT);
}

/***************************/
/* create/update text tray */
/***************************/
void gui_text_show(gint type, gchar *message)
{
static GtkWidget *view=NULL;
static GtkTextBuffer *buffer=NULL;
static GtkTextIter iter;

/* checks */
if (!message)
  return;

/* try this for improved stability of tasks */
/*
if (sysenv.running_tasks)
  return;
*/

if (sysenv.canvas)
  {
  if (!buffer)
    {
/* first time init */
    view = gtk_text_view_new();
    gtk_text_view_set_editable(GTK_TEXT_VIEW(view), FALSE);
    buffer = gtk_text_view_get_buffer(GTK_TEXT_VIEW(view));
    gtk_container_add(GTK_CONTAINER(sysenv.tpane), view);
    gtk_widget_show(view);
/* tags */
    gtk_text_buffer_create_tag(buffer, "fg_blue", "foreground", "blue", NULL);  
    gtk_text_buffer_create_tag(buffer, "fg_red", "foreground", "red", NULL);
    gtk_text_buffer_create_tag(buffer, "italic", "style", PANGO_STYLE_ITALIC, NULL);
/* position iterator */
    gtk_text_buffer_get_iter_at_line(buffer, &iter, 0);
    }

/* assign colour via message type */
  switch(type)
    {
    case ERROR:
      gtk_text_buffer_insert_with_tags_by_name
        (buffer, &iter, message, -1, "fg_red", NULL); 
      break;

    case WARNING:
      gtk_text_buffer_insert_with_tags_by_name
       (buffer, &iter, message, -1, "fg_blue", NULL); 
      break;

    case ITALIC:
      gtk_text_buffer_insert_with_tags_by_name
        (buffer, &iter, message, -1, "italic", NULL); 
      break;

    default:
      gtk_text_buffer_insert(buffer, &iter, message, -1);

    }

/* NEW - scroll to the end, so recent message is visible */
if (GTK_IS_TEXT_VIEW(view))
  {
/* FIXME - doesn't quite work as this only ensures the 1st printed line is visible */
gtk_text_view_scroll_to_iter(GTK_TEXT_VIEW(view), &iter, 0.0, FALSE, 0.0, 0.0);
  }

/* TODO - delete items from buffer when line gets too big */
  }
else
  printf("[%s]\n", message);
}

/***********************************/
/* general panel visibility toggle */
/***********************************/
void pane_toggle(GtkWidget *w, GtkWidget *frame)
{
sysenv.write_gdisrc = TRUE;

if (gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(w)))
  gtk_widget_show(frame);
else
  gtk_widget_hide(frame);
}

/*****************************************/
/* cleanup/save state before actual exit */
/*****************************************/
void gdis_exit(void)
{
/* if user has changed element colours, but doesn't want them saved */
/* this will do it anyway - but they can just click reset next time */
if (sysenv.write_gdisrc)
  {
  sysenv.x = sysenv.display_box->allocation.x;
  sysenv.y = sysenv.display_box->allocation.y;
  sysenv.width = sysenv.display_box->allocation.width;
  sysenv.height = sysenv.display_box->allocation.height;

  write_gdisrc();
  }

/* cleanup */
sys_free();
/*
unlink("gdis.ps");
*/

task_queue_free();

if (sysenv.canvas)
  {
  gtk_widget_pop_visual();
  gtk_widget_pop_colormap();
  gtk_exit(0);
  }
exit(0);
}

/************************************************/
/* exit confirmation if jobs are still runnning */
/************************************************/
gint gdis_exit_event(void)
{
gint button;
GtkWidget *dialog;

if (g_thread_pool_get_num_threads(sysenv.thread_pool))
  {
  dialog = gtk_message_dialog_new(GTK_WINDOW(window),
                                  GTK_DIALOG_DESTROY_WITH_PARENT,
                                  GTK_MESSAGE_WARNING,
                                  GTK_BUTTONS_YES_NO,
                                  "There are running tasks. Really exit?");
  button = gtk_dialog_run(GTK_DIALOG(dialog));
  gtk_widget_destroy(dialog);
  if (button != -8)
    return(TRUE);
  }
return(FALSE);
}

/************************************************/
/* exit confirmation if jobs are still runnning */
/************************************************/
void gdis_exit_test(void)
{
gint button;
GtkWidget *dialog;

if (g_thread_pool_get_num_threads(sysenv.thread_pool))
  {
  dialog = gtk_message_dialog_new(GTK_WINDOW(window),
                                  GTK_DIALOG_DESTROY_WITH_PARENT,
                                  GTK_MESSAGE_WARNING,
                                  GTK_BUTTONS_YES_NO,
                                  "There are running tasks. Really exit?");
  button = gtk_dialog_run(GTK_DIALOG(dialog));
  gtk_widget_destroy(dialog);
  if (button != -8)
    return;
  }
gdis_exit();
}

/***************************************/
/* display available modules & symbols */
/***************************************/
void mod_widget(void)
{
GSList *list, *list2;

for (list=sysenv.module_list ; list ; list=g_slist_next(list))
  {
  printf("module: %s\n", module_label(list->data));
  for (list2=module_symbols(list->data) ; list2 ; list2=g_slist_next(list2))
    {
    printf("  - %s\n", (gchar *) list2->data);
    }
  }
}

/**********************/
/* geomview importing */
/**********************/
/* TODO - put elsewhere */
void import_off(gchar *filename)
{
struct model_pak *model;

model = sysenv.active_model;
if (!model)
  return;

read_off(filename, model);

redraw_canvas(SINGLE);
}

void gui_import_geomview(void)
{
file_dialog("Import Geomview", NULL, FILE_LOAD, (gpointer) import_off, GEOMVIEW_OFF);
}

/*********************/
/* project importing */
/*********************/
void import_pcf(gchar *filename)
{
struct model_pak *model;

dialog_destroy_type(FILE_SELECT);

model = sysenv.active_model;
if (!model)
  return;

project_read(filename, model);

redraw_canvas(SINGLE);
}

void gui_import_project(void)
{
file_dialog("Import Project", NULL, FILE_LOAD, (gpointer) import_pcf, PROJECT);
}

void gui_import_graph(void)
{
file_dialog("Import Graph", NULL, FILE_LOAD, (gpointer) graph_read, TEXT);
}

/*****************/
/* mode switches */
/*****************/
void gui_mode_record(void)
{
gui_mode_switch(RECORD);
}

void gui_mode_default(void)
{
gui_mode_switch(FREE);
}

/******************/
/* MENU structure */
/******************/
static GtkItemFactoryEntry menu_items[] = 
{
  { "/_File",                 NULL, NULL, 0, "<Branch>" },
/*
  { "/File/_New",             NULL, create_new_model, 1, NULL },
  { "/File/sep1",             NULL, NULL, 0, "<Separator>" },
*/
  { "/File/_Open...",         NULL, file_load_dialog, 1, NULL },
  { "/File/_Save...",         NULL, file_save_dialog, 1, NULL },
  { "/File/_Close",           NULL, tree_select_delete, 1, NULL },
  { "/File/sep1",             NULL, NULL, 0, "<Separator>" },

  { "/File/Import",              NULL, NULL, 0, "<Branch>" },
  { "/File/Import/Geomview...",  NULL, gui_import_geomview, 1, NULL },
  { "/File/Import/Project...",   NULL, gui_import_project, 1, NULL },
  { "/File/Import/Graph...",     NULL, gui_import_graph, 1, NULL },

  { "/File/Export",                     NULL, NULL, 0, "<Branch>" },
  { "/File/Export/Canvas snapshot...",  NULL, image_export_dialog, 1, NULL },
  { "/File/Export/Graph data...",       NULL, analysis_export_dialog, 1, NULL },

  { "/File/sep1",             NULL, NULL, 0, "<Separator>" },
  { "/File/_Quit",            NULL, gdis_exit_test, 0, NULL },

  { "/_Edit",               NULL, NULL, 0, "<Branch>" },
  { "/Edit/_Copy",          NULL, select_copy, 0, NULL },
  { "/Edit/_Paste",         NULL, select_paste, 0, NULL },
  { "/Edit/sep1",           NULL, NULL, 0, "<Separator>" },

  { "/Edit/Delete",         NULL, select_delete, 0, NULL },
  { "/Edit/Undo",           NULL, undo_active, 0, NULL },

  { "/Edit/sep1",           NULL, NULL, 0, "<Separator>" },

  { "/Edit/Colour...",      NULL, select_colour, 0, NULL },
  { "/Edit/Hide",           NULL, select_hide, 0, NULL },
  { "/Edit/Unhide all",     NULL, unhide_atoms, 0, NULL },
  { "/Edit/sep1",           NULL, NULL, 0, "<Separator>" },
  { "/Edit/Select all",     NULL, select_all, 0, NULL },
  { "/Edit/Invert",         NULL, select_invert, 0, NULL },

  { "/_Tools",                                NULL, NULL, 0, "<Branch>" },
  { "/Tools/Visualization",                   NULL, NULL, 0, "<Branch>" },
  { "/Tools/Visualization/Animation...",      NULL, gui_animate_dialog, 0, NULL },
  { "/Tools/Visualization/Iso-surfaces...",   NULL, gui_isosurf_dialog, 0, NULL },
  { "/Tools/Visualization/Periodic table...", NULL, gui_gperiodic_dialog, 0, NULL },

  { "/Tools/Building",                           NULL, NULL, 0, "<Branch>" },
  { "/Tools/Building/Editing...",                NULL, gui_edit_dialog, 0, NULL },
  { "/Tools/Building/Dislocations...",           NULL, gui_defect_dialog, 0, NULL },
  { "/Tools/Building/Docking...",                NULL, gui_dock_dialog, 0, NULL },
  { "/Tools/Building/Dynamics...",               NULL, gui_mdi_dialog, 0, NULL },
  { "/Tools/Building/Surfaces...",               NULL, surface_dialog, 0, NULL },
  { "/Tools/Building/Zmatrix...",                NULL, gui_zmat_dialog, 0, NULL },

  { "/Tools/Computation",                     NULL, NULL, 0, "<Branch>" },
  { "/Tools/Computation/Diffraction...",      NULL, gui_diffract_dialog, 0, NULL },
  { "/Tools/Computation/GULP...",             NULL, gulp_dialog, 0, NULL },
  { "/Tools/Computation/GAMESS...",           NULL, gamess_dialog, 0, NULL },
  { "/Tools/Computation/Monty...",            NULL, monty_dialog, 0, NULL },
  { "/Tools/Computation/SIESTA...",           NULL, gui_siesta_dialog, 0, NULL },
  { "/Tools/Computation/VASP...",             NULL, gui_vasp_dialog, 0, NULL },
  { "/Tools/Computation/USPEX...",            NULL, gui_uspex_dialog, 0, NULL },

  { "/Tools/Analysis",                        NULL, NULL, 0, "<Branch>" },
  { "/Tools/Analysis/Dynamics...",            NULL, gui_analysis_dialog, 0, NULL },
  { "/Tools/Analysis/Measurements...",        NULL, gui_measure_dialog, 0, NULL },
  { "/Tools/Analysis/Plots...",               NULL, gui_plots_dialog, 0, NULL },

  { "/_View",                       NULL, NULL, 0, "<Branch>"},
  { "/View/Display properties...",  NULL, gui_render_dialog, 0, NULL},
  { "/View/sep1",                   NULL, NULL, 0, "<Separator>"},
  { "/View/Reset model images",     NULL, space_image_widget_reset, 0, NULL},
  { "/View/sep1",                   NULL, NULL, 0, "<Separator>"},
  { "/View/Normal mode",            NULL, gui_mode_default, 0, NULL},
  { "/View/Recording mode",         NULL, gui_mode_record, 0, NULL},
  { "/View/sep1",                   NULL, NULL, 0, "<Separator>"},
  { "/View/Task manager...",        NULL, task_dialog, 0, NULL},

#ifdef WITH_GRISU
  { "/View/Grid manager...",      NULL, gui_grid_dialog, 0, NULL},
#endif

  { "/View/Executable paths...",    NULL, gui_setup_dialog, 0, NULL},

  { "/_Help",                  NULL, NULL, 0, "<Branch>"},

/* about info -> manual acknowlegements */
/*
  { "/Help/About...",          NULL, gui_help_dialog, 0, NULL},
*/
  { "/Help/Manual...",         NULL, gui_help_dialog, 0, NULL},
};

/********************************/
/* process key presses directly */
/********************************/
gint cb_key_press(GtkWidget *w, GdkEventKey *event, gpointer dummy)
{
switch(event->keyval)
  {
/* selection delete */
  case GDK_Insert:
    undo_active();
    break;

/* selection delete */
  case GDK_Delete:
    select_delete();
    break;

/* fullscreen */
  case GDK_F1:
    if (!sysenv.stereo)
      {
      sysenv.stereo = TRUE;
      sysenv.stereo_fullscreen = TRUE;
      sysenv.render.perspective = TRUE;
      stereo_open_window();
      }
    break;

/* windowed */
/* CURRENT - hacked to work with new camera code */
  case GDK_F2:
/*
    if (!sysenv.stereo && sysenv.stereo_windowed)
*/

    if (!sysenv.stereo)
      {
      sysenv.stereo = TRUE;
      sysenv.stereo_fullscreen = FALSE;
      sysenv.render.perspective = TRUE;

      redraw_canvas(SINGLE);
      }
    break;

/* NEW - switching between windowed/fullscreen stereo while in */
/* stereo mode can screw things up - so just use esc to exit */
/* for both. ie F1/F2 do not act as stereo toggles */
  case GDK_Escape:
    sysenv.stereo = FALSE;
    sysenv.render.perspective = FALSE;
    stereo_close_window();
    redraw_canvas(SINGLE);
    break;
  }
return(FALSE);
}

/*******************************************************/
/* record current model/property pane divider position */
/*******************************************************/
gboolean cb_pane_refresh(GtkPaned *w, gpointer data)
{
sysenv.tree_divider = gtk_paned_get_position(w);
return(FALSE);
}

/***************************************/
/* update the GDIS information widgets */
/***************************************/
gint gui_widget_handler(void)
{
/* TODO - all components of the gui */
if (sysenv.refresh_properties)
  {
  gui_active_refresh();
  sysenv.refresh_properties = FALSE;
  }
return(TRUE);
}

/*****************************/
/* schedule widget update(s) */
/*****************************/
/* TODO - include all update request (including canvas) */
/* TODO - make type a mask so that multiple updates can be done */
void gui_refresh(gint type)
{
switch (type)
  {
  case GUI_CANVAS:
    redraw_canvas(SINGLE);
    break;
  case GUI_MODEL_TREE:
    sysenv.refresh_tree = TRUE;
    break;
  case GUI_MODEL_PROPERTIES:
    sysenv.refresh_properties = TRUE;
    break;
  case GUI_TEXT_BUFFER:
    sysenv.refresh_text = TRUE;
    break;
  }
}

/****************************/
/* selection model callback */
/****************************/
void gui_selection_mode_set(GtkWidget *w, gpointer dummy)
{
const gchar *line;

g_assert(GTK_IS_ENTRY(w));

line = gtk_entry_get_text(GTK_ENTRY(w));

gui_mode_switch(FREE);

if (g_strrstr(line, "Atoms"))
  {
  sysenv.select_mode = CORE;
  return;
  }
if (g_strrstr(line, "Label"))
  {
  sysenv.select_mode = ATOM_LABEL;
  return;
  }
if (g_strrstr(line, "Type"))
  {
  sysenv.select_mode = ATOM_TYPE;
  return;
  }
if (g_strrstr(line, "Elements"))
  {
  if (g_strrstr(line, "molecule"))
    sysenv.select_mode = ELEM_MOL;
  else
    sysenv.select_mode = ELEM;
  return;
  }
if (g_strrstr(line, "Molecules"))
  {
  sysenv.select_mode = MOL;
  return;
  }
if (g_strrstr(line, "fragments"))
  {
  sysenv.select_mode = FRAGMENT;
  gui_mode_switch(SELECT_FRAGMENT);
  return;
  }
if (g_strrstr(line, "Regions"))
  {
  sysenv.select_mode = REGION;
  return;
  }
}


/* CURRENT */
/* functionality for extending active model pulldown capabilities */
/* TODO - put in shortcuts */

struct active_pak
{
gchar *label;
void (*setup)(gpointer);
void (*refresh)(gpointer);
GtkWidget *vbox;
};

GList *active_list=NULL;
struct active_pak *active_data=NULL;

/***********************************/
/* refresh the active model widget */
/***********************************/
void gui_active_refresh(void)
{
GList *list;
struct active_pak *current;

for (list=active_list ; list ; list=g_list_next(list))
  {
  current = list->data;

  if (active_data == current)
    {
/* restore "proper" size */
    gtk_widget_set_size_request(current->vbox, -1, -1);
    gtk_widget_show(current->vbox);
    }
  else
    gtk_widget_hide(current->vbox);
  }

if (active_data)
  {
  if (active_data->refresh)
    active_data->refresh(NULL);
  }
}
/*********************************************/
/* toggle tracking mode (if possible) --OVHPA*/
/*********************************************/
void gui_track_output(void){
struct model_pak *model;
gchar *ptr;
/**/
model = sysenv.active_model;
if (!model) return;
/*invert current tracking**/
model->track_me=(model->track_me==FALSE);
switch(model->id){
case USPEX:
	if(model->track_me) {
		g_timeout_add_full(G_PRIORITY_DEFAULT,TRACKING_TIMEOUT,track_uspex,model,track_uspex_cleanup);
		if(model->track_me) {
			ptr=g_strdup_printf("USPEX TRACKING: START.\n");
			gui_text_show(ITALIC,ptr);
			g_free(ptr);
		}
	}else{
		track_uspex(model);/*will return FALSE**/
		model->track_me=FALSE;/*just in case we need to be very clear*/
	}
	break;
case VASP:
	if(model->track_me) {
		g_timeout_add_full(G_PRIORITY_DEFAULT,TRACKING_TIMEOUT,track_vasp,model,track_vasp_cleanup);
		if(model->track_me) {
			ptr=g_strdup_printf("VASP TRACKING: START.\n");
			gui_text_show(ITALIC,ptr);
			g_free(ptr);
		}
	}else{
		track_vasp(model);/*will return FALSE**/
		model->track_me=FALSE;/*just in case we need to be very clear*/
	}
	break;
default:
	return;/*nothing to do*/
}

}
/***************************************/
/* quick - export view into eps --OVHPA*/
/***************************************/
void gui_eps_export(void){
struct model_pak *model;
GtkWidget *file_chooser;
GtkFileFilter *filter;

/**/
model = sysenv.active_model;
if (!model) return;
filter=gtk_file_filter_new();
if(sysenv.have_eps){
	gtk_file_filter_add_pattern(filter,"*.eps");
	gtk_file_filter_set_name (filter,"eps image file");
	file_chooser = gtk_file_chooser_dialog_new("snapshot to eps",
		GTK_WINDOW(sysenv.main_window),
		GTK_FILE_CHOOSER_ACTION_SAVE,GTK_STOCK_CANCEL,
		GTK_RESPONSE_CANCEL,GTK_STOCK_SAVE,GTK_RESPONSE_ACCEPT,
		NULL);
}else{
	gtk_file_filter_add_pattern(filter,"*.png");
	gtk_file_filter_set_name (filter,"png image file");
	file_chooser = gtk_file_chooser_dialog_new("snapshot to png",
		GTK_WINDOW(sysenv.main_window),
		GTK_FILE_CHOOSER_ACTION_SAVE,GTK_STOCK_CANCEL,
		GTK_RESPONSE_CANCEL,GTK_STOCK_SAVE,GTK_RESPONSE_ACCEPT,
		NULL);
}
gtk_file_chooser_add_filter (GTK_FILE_CHOOSER(file_chooser),filter);
if(gtk_dialog_run(GTK_DIALOG(file_chooser)) == GTK_RESPONSE_ACCEPT) {
	if(model->eps_file!=NULL) g_free(model->eps_file);
	model->eps_file=gtk_file_chooser_get_filename (GTK_FILE_CHOOSER (file_chooser));
	model->snapshot_eps=TRUE;
	model->redraw = TRUE;
}else{
	model->snapshot_eps=FALSE;/*useless*/
}
gtk_widget_destroy (GTK_WIDGET(file_chooser));
redraw_canvas(ALL);
}

/****************************************/
/* active pulldown change event handler */
/****************************************/
void gui_active_handler(GtkWidget *w, gpointer dummy)
{
const gchar *tmp;
GList *list;
struct active_pak *active;

tmp = gtk_entry_get_text(GTK_ENTRY(w));

for (list=active_list ; list ; list=g_list_next(list))
  {
  active = list->data;

  if (g_ascii_strcasecmp(active->label, tmp) == 0)
    {
    active_data = active;
    gui_active_refresh();
    return;
    }
  }
}

/*********************************/
/* setup the active model widget */
/*********************************/
void gui_active_setup(GtkWidget *box)
{
gpointer entry;
GList *list, *menu;
struct active_pak *active;

/* build the list of labels */
menu = NULL;
for (list=active_list ; list ; list=g_list_next(list))
  {
  active = list->data;
  menu = g_list_append(menu, g_strdup(active->label));
  }

/* create the header pulldown */
entry = gui_pulldown_new(NULL, menu, FALSE, box);

/* create and pack the data boxes */
for (list=active_list ; list ; list=g_list_next(list))
  {
  active = list->data;

  gtk_box_pack_start(GTK_BOX(box), active->vbox, FALSE, FALSE, 0);

/* stop GTK maxing out the allocated height for all active model children */
  gtk_widget_set_size_request(active->vbox, -1, 1);

  if (active->setup)
    active->setup(active->vbox);
  }

/* top item is active by default */
if (active_list)
  active_data = active_list->data;

g_signal_connect(GTK_OBJECT(entry), "changed", GTK_SIGNAL_FUNC(gui_active_handler), NULL);
}

/*************************************************/
/* create a new entry in the active model widget */
/*************************************************/
void gui_active_new(const gchar *label, gpointer setup, gpointer refresh)
{
struct active_pak *active;

active = g_malloc(sizeof(struct active_pak));

active->label = g_strdup(label);
active->setup = setup;
active->refresh = refresh;
active->vbox = gtk_vbox_new(FALSE, 0);

active_list = g_list_append(active_list, active);
}

/*********/
/* SETUP */
/*********/
void gui_init(int argc, char *argv[])
{
gint nmenu_items = sizeof (menu_items) / sizeof (menu_items[0]);
gchar *text;
GList *list;
gpointer ptr;
GtkWidget *gdis_wid;
GtkWidget *hpaned, *vpaned;
GtkWidget *vbox, *hbox, *vbox_lower, *menu_bar, *toolbar;
GtkWidget *frame, *event_box;
GdkBitmap *mask;
GdkPixmap *gdis_pix=NULL;
GdkPixbuf *pixbuf;
GtkStyle *style;
GtkItemFactory *item;
GdkColor colour;

gtk_init(&argc, &argv);
gdk_gl_init(&argc, &argv);
gtk_gl_init(&argc, &argv);

/* enforce true colour (fixes SG problems) */
sysenv.visual = gdk_visual_get_best_with_type(GDK_VISUAL_TRUE_COLOR);
if (!sysenv.visual)
  {
  printf("Error: could not get requested visual.\n");
  exit(1);
  }
sysenv.depth = sysenv.visual->depth;
sysenv.colourmap = gdk_colormap_new(sysenv.visual, TRUE);
if (!sysenv.colourmap)
  {
  printf("Error: could not allocate colourmap.\n");
  exit(1);
  }
gtk_widget_push_colormap(sysenv.colourmap);
gtk_widget_push_visual(sysenv.visual);


/* NB: GTK graphical init needs to be done before this is called, */
/* NB: but this must be done before we build the GDIS interface */
image_table_init();


/* main window */
window = gtk_window_new(GTK_WINDOW_TOPLEVEL);
sysenv.main_window = window;/*FIXME: move to sysenv. --OVHPA*/
gtk_window_set_policy(GTK_WINDOW(window), TRUE, TRUE, FALSE);
gtk_window_set_title(GTK_WINDOW(window),"GTK Display Interface for Structures");
gtk_window_set_position(GTK_WINDOW(window), GTK_WIN_POS_CENTER);

g_signal_connect(GTK_OBJECT(window), "key_press_event",
                (GtkSignalFunc) cb_key_press, NULL);


/* vbox for the menubar */
vbox = gtk_vbox_new(FALSE, 0);
gtk_container_add(GTK_CONTAINER(window), vbox);

/* item factory menu creation */
item = gtk_item_factory_new(GTK_TYPE_MENU_BAR, "<main>", NULL);
gtk_item_factory_create_items(item, nmenu_items, menu_items, NULL);
menu_bar = gtk_item_factory_get_widget(item, "<main>");

/* FALSE,FALSE => don't expand to fill (eg on resize) */
gtk_box_pack_start(GTK_BOX(vbox), menu_bar, FALSE, FALSE, 0);

/* hbox for the toolbar */
hbox = gtk_hbox_new(FALSE, 0);
gtk_box_pack_start(GTK_BOX(vbox), hbox, FALSE, FALSE, 0);
toolbar = gtk_toolbar_new();
gtk_box_pack_start(GTK_BOX(hbox), toolbar, FALSE, FALSE, 0);

/* extract some important style info */
gtk_widget_realize(window);
style = gtk_widget_get_style(window);

sysenv.gtk_fontsize = pango_font_description_get_size(style->font_desc) 
                    / PANGO_SCALE;

/*
printf("psize = %d\n", pango_font_description_get_size(style->font_desc));
printf("scale = %d\n", PANGO_SCALE);
printf(" size = %d\n", sysenv.gtk_fontsize);
*/

/* load button */
pixbuf = image_table_lookup("image_folder");
gdis_wid = gtk_image_new_from_pixbuf(pixbuf);
gtk_toolbar_append_item(GTK_TOOLBAR(toolbar), NULL,
                        "Open new model", "Private",
                        gdis_wid, GTK_SIGNAL_FUNC(file_load_dialog), NULL);

/* save button */
pixbuf = image_table_lookup("image_disk");
gdis_wid = gtk_image_new_from_pixbuf(pixbuf);
gtk_toolbar_append_item(GTK_TOOLBAR (toolbar), NULL,
                        "Save active model", "Private", gdis_wid,
                        GTK_SIGNAL_FUNC(file_save_dialog), NULL);


/* NEW: a 'new model' button */
pixbuf = image_table_lookup("image_plus");
gdis_wid = gtk_image_new_from_pixbuf(pixbuf);
gtk_toolbar_append_item(GTK_TOOLBAR (toolbar),
                        NULL,
                        "New model",
                        "Private",
                        gdis_wid,
                        GTK_SIGNAL_FUNC(edit_model_create),
                        NULL);

/* delete button */
pixbuf = image_table_lookup("image_cross");
gdis_wid = gtk_image_new_from_pixbuf(pixbuf);
gtk_toolbar_append_item(GTK_TOOLBAR (toolbar),
                        NULL,
                        "Close active model",
                        "Private",
                        gdis_wid,
                        GTK_SIGNAL_FUNC(tree_select_delete),
                        NULL);

gtk_toolbar_append_space(GTK_TOOLBAR(toolbar));

/* animation */
pixbuf = image_table_lookup("image_animate");
gdis_wid = gtk_image_new_from_pixbuf(pixbuf);
gtk_toolbar_append_item(GTK_TOOLBAR (toolbar),
                        NULL,
                        "Animation",
                        "Private",
                        gdis_wid,
                        GTK_SIGNAL_FUNC(gui_animate_dialog),
                        NULL);

/* model editing */
pixbuf = image_table_lookup("image_tools");
gdis_wid = gtk_image_new_from_pixbuf(pixbuf);
gtk_toolbar_append_item(GTK_TOOLBAR (toolbar),
                        NULL,
                        "Model editing",
                        "Private",
                        gdis_wid,
                        GTK_SIGNAL_FUNC(gui_edit_dialog),
                        NULL);

/* iso surfaces */
pixbuf = image_table_lookup("image_isosurface");
gdis_wid = gtk_image_new_from_pixbuf(pixbuf);
gtk_toolbar_append_item(GTK_TOOLBAR (toolbar),
                        NULL,
                        "Iso-surfaces",
                        "Private",
                        gdis_wid,
                        GTK_SIGNAL_FUNC(gui_isosurf_dialog),
                        NULL);

/* gperiodic button */
pixbuf = image_table_lookup("image_element");
gdis_wid = gtk_image_new_from_pixbuf(pixbuf);
gtk_toolbar_append_item(GTK_TOOLBAR (toolbar),
                        NULL,
                        "Periodic table",
                        "Private",
                        gdis_wid,
                        GTK_SIGNAL_FUNC(gui_gperiodic_dialog),
                        NULL);

gtk_toolbar_append_space(GTK_TOOLBAR(toolbar));

/* diffraction button */
pixbuf = image_table_lookup("image_diffraction");
gdis_wid = gtk_image_new_from_pixbuf(pixbuf);
gtk_toolbar_append_item(GTK_TOOLBAR (toolbar),
                        NULL,
                        "Diffraction",
                        "Private",
                        gdis_wid,
                        GTK_SIGNAL_FUNC(gui_diffract_dialog),
                        NULL);

/* surface button */
pixbuf = image_table_lookup("image_surface");
gdis_wid = gtk_image_new_from_pixbuf(pixbuf);
gtk_toolbar_append_item(GTK_TOOLBAR (toolbar),
                        NULL,
                        "Surface construction",
                        "Private",
                        gdis_wid,
                        GTK_SIGNAL_FUNC(surface_dialog),
                        NULL);

gtk_toolbar_append_space(GTK_TOOLBAR(toolbar));

/* geometry button */
pixbuf = image_table_lookup("image_compass");
gdis_wid = gtk_image_new_from_pixbuf(pixbuf);
gtk_toolbar_append_item(GTK_TOOLBAR (toolbar),
                        NULL,
                        "Measurements",
                        "Private",
                        gdis_wid,
                        GTK_SIGNAL_FUNC(gui_measure_dialog),
                        NULL);

gtk_toolbar_append_space(GTK_TOOLBAR(toolbar));

/* display properties */
pixbuf = image_table_lookup("image_palette");
gdis_wid = gtk_image_new_from_pixbuf(pixbuf);
gtk_toolbar_append_item(GTK_TOOLBAR (toolbar),
                        NULL,
                        "Display properties",
                        "Private",
                        gdis_wid,
                        GTK_SIGNAL_FUNC(gui_render_dialog),
                        NULL);

/* model geomtry */
pixbuf = image_table_lookup("image_axes");
gdis_wid = gtk_image_new_from_pixbuf(pixbuf);
gtk_toolbar_append_item(GTK_TOOLBAR(toolbar),
                        NULL,
                        "Reset model geometry",
                        "Private",
                        gdis_wid,
                        GTK_SIGNAL_FUNC(gui_view_default),
                        NULL);

/* model images */
pixbuf = image_table_lookup("image_periodic");
gdis_wid = gtk_image_new_from_pixbuf(pixbuf);
gtk_toolbar_append_item(GTK_TOOLBAR(toolbar),
                        NULL,
                        "Reset model images",
                        "Private",
                        gdis_wid,
                        GTK_SIGNAL_FUNC(space_image_widget_reset),
                        NULL);

/* transformation record button */
pixbuf = image_table_lookup("image_camera");
gdis_wid = gtk_image_new_from_pixbuf(pixbuf);
gtk_toolbar_append_item(GTK_TOOLBAR (toolbar),
                        NULL,
                        "Record mode",
                        "Private",
                        gdis_wid,
                        GTK_SIGNAL_FUNC(gtk_mode_switch),
                        GINT_TO_POINTER(RECORD));

gtk_toolbar_append_space(GTK_TOOLBAR(toolbar));

/* single canvas */
pixbuf = image_table_lookup("image_canvas_single");
gdis_wid = gtk_image_new_from_pixbuf(pixbuf);
gtk_toolbar_append_item(GTK_TOOLBAR (toolbar),
                        NULL,
                        "Single canvas",
                        "Private",
                        gdis_wid,
                        GTK_SIGNAL_FUNC(canvas_single),
                        GINT_TO_POINTER(CANVAS_SINGLE));

/* create canvas */
pixbuf = image_table_lookup("image_canvas_create");
gdis_wid = gtk_image_new_from_pixbuf(pixbuf);
gtk_toolbar_append_item(GTK_TOOLBAR (toolbar),
                        NULL,
                        "Increase number of canvasses",
                        "Private",
                        gdis_wid,
                        GTK_SIGNAL_FUNC(canvas_create),
                        NULL);

/* delete canvas */
pixbuf = image_table_lookup("image_canvas_delete");
gdis_wid = gtk_image_new_from_pixbuf(pixbuf);
gtk_toolbar_append_item(GTK_TOOLBAR (toolbar),
                        NULL,
                        "Decrease number of canvasses",
                        "Private",
                        gdis_wid,
                        GTK_SIGNAL_FUNC(canvas_delete),
                        NULL);

gtk_toolbar_append_space(GTK_TOOLBAR(toolbar));

/* normal viewing button */
pixbuf = image_table_lookup("image_arrow");
gdis_wid = gtk_image_new_from_pixbuf(pixbuf);
gtk_toolbar_append_item(GTK_TOOLBAR (toolbar),
                        NULL,
                        "Reset model mode",
                        "Private",
                        gdis_wid,
                        GTK_SIGNAL_FUNC(gtk_mode_switch),
                        GINT_TO_POINTER(FREE));

/* plot control button */
pixbuf = image_table_lookup("image_plots");
gdis_wid = gtk_image_new_from_pixbuf(pixbuf);
gtk_toolbar_append_item(GTK_TOOLBAR (toolbar),
			NULL,
			"Plot controls",
			"Private",
			gdis_wid,
			GTK_SIGNAL_FUNC(gui_graph_controls),
			NULL);
/* eps/png export tool */
if(sysenv.have_eps){
	pixbuf = image_table_lookup("image_to_eps");
	gdis_wid = gtk_image_new_from_pixbuf(pixbuf);
	gtk_toolbar_append_item(GTK_TOOLBAR (toolbar),
			NULL,
			"Export to EPS",
			"Private",
			gdis_wid,
			GTK_SIGNAL_FUNC(gui_eps_export),
			NULL);
}else{
	pixbuf = image_table_lookup("image_to_eps");
	gdis_wid = gtk_image_new_from_pixbuf(pixbuf);
	gtk_toolbar_append_item(GTK_TOOLBAR (toolbar),
			NULL,
			"Export to PNG",
			"Private",
			gdis_wid,
			GTK_SIGNAL_FUNC(gui_eps_export),
			NULL);
}
/* tracking system */
pixbuf = image_table_lookup("image_track");
gdis_wid = gtk_image_new_from_pixbuf(pixbuf);
gtk_toolbar_append_item(GTK_TOOLBAR (toolbar),
			NULL,
			"Track output",
			"Private",
			gdis_wid,
			GTK_SIGNAL_FUNC(gui_track_output),
			NULL);
/* MAIN LEFT/RIGHT HBOX PANE */
/* paned window */
hpaned = gtk_hpaned_new();
gtk_container_add(GTK_CONTAINER(GTK_BOX(vbox)), hpaned);

/* LEFT PANE */
sysenv.mpane = gtk_vbox_new(FALSE, 0);
gtk_paned_pack1(GTK_PANED(hpaned), sysenv.mpane, TRUE, TRUE);
gtk_widget_set_size_request(sysenv.mpane, sysenv.tree_width, -1);

/* model tree box */
vbox = gtk_vbox_new(FALSE, 0);
gtk_box_pack_start(GTK_BOX(sysenv.mpane), vbox, TRUE, TRUE, 0);
tree_init(vbox);

/* model pulldown */
vbox_lower = gtk_vbox_new(FALSE, 0);
gtk_box_pack_start(GTK_BOX(sysenv.mpane), vbox_lower, FALSE, FALSE, 0);

/* general model data box */
vbox = gtk_vbox_new(FALSE, 0);
gtk_box_pack_start(GTK_BOX(vbox_lower), vbox, FALSE, FALSE, 0);

/* CURRENT */
gui_active_new("Model : Content", gui_content_refresh, gui_content_refresh);
gui_active_new("Model : Editing", gui_edit_widget, gui_refresh_selection);
gui_active_new("Model : Display", gui_display_widget, NULL);
gui_active_new("Model : Images",  space_image_widget_setup, space_image_widget_redraw);
gui_active_new("Model : Symmetry", gui_symmetry_refresh, gui_symmetry_refresh);
gui_active_new("Model : Viewing", gui_view_widget, NULL);
gui_active_setup(vbox);


/* GTK is a bit dumb when it sets the initial sizes of the two panes */
/* set the pane position and store when changed */
/*
gtk_paned_set_position(GTK_PANED(vpaned), sysenv.tree_divider);
g_signal_connect(GTK_OBJECT(vpaned), "size-request",
                 G_CALLBACK(cb_pane_refresh), NULL);
*/

/* selection pulldown */
vbox_lower = gtk_vbox_new(FALSE, 0);
gtk_box_pack_start(GTK_BOX(sysenv.mpane), vbox_lower, FALSE, FALSE, 0);
list = NULL;
list = g_list_append(list, "Select : Atoms");
list = g_list_append(list, "Select : Atom Label");
list = g_list_append(list, "Select : Atom FF Type");
list = g_list_append(list, "Select : Elements");
list = g_list_append(list, "Select : Elements in Molecule");
list = g_list_append(list, "Select : Molecules");
list = g_list_append(list, "Select : Molecule Fragments");
list = g_list_append(list, "Select : Regions");
ptr = gui_pulldown_new(NULL, list, FALSE, vbox_lower);
g_signal_connect(GTK_OBJECT(ptr), "changed", GTK_SIGNAL_FUNC(gui_selection_mode_set), NULL);


/* GDIS logo */
frame = gtk_frame_new(NULL);
gtk_box_pack_end(GTK_BOX(sysenv.mpane), frame, FALSE, TRUE, 0);

/* event box (get an x window for setting black background) */
event_box = gtk_event_box_new();
gtk_container_add(GTK_CONTAINER(frame), event_box);

/* black background */
colour.red = 0.0;
colour.green = 0.0;
colour.blue = 0.0;
gtk_widget_modify_bg(event_box, GTK_STATE_NORMAL, &colour);

hbox = gtk_hbox_new (FALSE, 10);
gtk_container_add (GTK_CONTAINER(event_box), hbox);
gtk_container_set_border_width(GTK_CONTAINER(hbox), PANEL_SPACING);


gdis_pix = gdk_pixmap_create_from_xpm_d(window->window, &mask,
                          &style->bg[GTK_STATE_NORMAL], logo_left_xpm);
gdis_wid = gtk_image_new_from_pixmap(gdis_pix, mask);
gtk_box_pack_start(GTK_BOX(hbox), gdis_wid, FALSE, FALSE, 0);


gdis_pix = gdk_pixmap_create_from_xpm_d(window->window, &mask,
                          &style->bg[GTK_STATE_NORMAL], logo_right_81_xpm);
gdis_wid = gtk_image_new_from_pixmap(gdis_pix, mask);
gtk_box_pack_end(GTK_BOX(hbox), gdis_wid, FALSE, FALSE, 0);


/* RIGHT PANE */
vbox = gtk_vbox_new(FALSE, 0);
gtk_paned_pack2(GTK_PANED(hpaned), vbox, TRUE, TRUE);

/* paned window */
vpaned = gtk_vpaned_new();
gtk_container_add(GTK_CONTAINER(vbox), vpaned);

sysenv.display_box = gtk_vbox_new(FALSE, 0);
gtk_paned_pack1(GTK_PANED(vpaned), sysenv.display_box, TRUE, TRUE);

/* create the main drawing area */
if (gl_init_visual())
  exit(1);

canvas_init(sysenv.display_box);
canvas_new(0, 0, sysenv.width, sysenv.height);

/* text status pane */
sysenv.tpane = gtk_scrolled_window_new(NULL, NULL);
gtk_scrolled_window_set_policy(GTK_SCROLLED_WINDOW(sysenv.tpane),
                               GTK_POLICY_AUTOMATIC, GTK_POLICY_AUTOMATIC);
gtk_paned_pack2(GTK_PANED(vpaned), sysenv.tpane, TRUE, TRUE);
gtk_widget_set_size_request(sysenv.tpane, -1, sysenv.tray_height);
gtk_widget_show(sysenv.tpane);

text = g_strdup_printf("This is free software, distributed under the terms of the GNU public license (GPL).\nFor more information visit http://www.gnu.org\n");
gui_text_show(WARNING, text);
g_free(text);

text = g_strdup_printf("Welcome to GDIS version %4.2f.%d, Copyright (C) %d by Sean Fleming and Andrew Rohl\n",VERSION,PATCH,YEAR); 
gui_text_show(STANDARD, text);
g_free(text);

/* confirmation callback (ie really exit?) */
g_signal_connect(GTK_OBJECT(window), "delete-event", GTK_SIGNAL_FUNC(gdis_exit_event), NULL);
/* save gdisrc callback */
g_signal_connect(GTK_OBJECT(window), "destroy", GTK_SIGNAL_FUNC(gdis_exit), NULL);

/* show all */
gtk_widget_show_all(window);

/* CURRENT */
gui_active_refresh();
}

