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
#include <gtk/gtkgl.h>
#ifdef __APPLE__
#include <OpenGL/gl.h>
#else
#include <GL/gl.h>
#endif

#include "gdis.h"
#include "opengl.h"
#include "interface.h"

extern struct sysenv_pak sysenv;

/***************************/
/* free a canvas structure */
/***************************/
void canvas_free(gpointer data)
{
struct canvas_pak *canvas = data;
g_free(canvas);
}

/****************************/
/* schedule redraw requests */ 
/****************************/
void redraw_canvas(gint action)
{
GSList *list;
struct model_pak *model;

switch (action)
  {
  case SINGLE:
/*
    data = sysenv.active_model;
    if (data)
      data->redraw = TRUE;
    break;
*/
  case ALL:
    for (list=sysenv.mal ; list ; list=g_slist_next(list))
      {
      model = list->data;
      model->redraw = TRUE;
      }
    break;
  }
sysenv.refresh_canvas = TRUE;
}

/*******************/
/* configure event */
/*******************/
#define DEBUG_GL_CONFIG_EVENT 0
gint canvas_configure(GtkWidget *w, GdkEventConfigure *event, gpointer data)
{
gint size;
GSList *list;
struct model_pak *model;

g_assert(w != NULL);

/* store new drawing area size */
if (w->allocation.width > w->allocation.height)
  size = w->allocation.height;
else
  size = w->allocation.width;
sysenv.x = 0;
sysenv.y = 0;
sysenv.width = w->allocation.width;
sysenv.height = w->allocation.height;
sysenv.size = size;

/* update canvases */
canvas_resize();

#if DEBUG_GL_CONFIG_EVENT
printf("Relative canvas origin: (%d,%d)\n",sysenv.x,sysenv.y);
printf("     Canvas dimensions:  %dx%d\n",sysenv.width,sysenv.height);
#endif

/* update coords */
for (list=sysenv.mal ; list ; list=g_slist_next(list))
  {
  model = list->data;
  model->redraw = TRUE;
  }

/* new screen size to be saved as default */
sysenv.write_gdisrc = TRUE;

return(TRUE);
}

/*****************/
/* expose event */
/****************/
#define DEBUG_GL_EXPOSE 0
gint canvas_expose(GtkWidget *w, GdkEventExpose *event, gpointer data)
{
/*
gl_clear_canvas();
for (list=sysenv.mal ; list ; list=g_slist_next(list))
  {
  model = list->data;
  model->redraw = TRUE;
  canvas = model->canvas;
  if (canvas->active)
    gl_draw(canvas, model);
  }
*/

redraw_canvas(ALL);

return(TRUE);
}

/*****************************************************/
/* create a new canvas and place in the canvas table */
/*****************************************************/
void canvas_new(gint x, gint y, gint w, gint h)
{
struct canvas_pak *canvas;

/* create an OpenGL capable drawing area */
canvas = g_malloc(sizeof(struct canvas_pak));
/*
printf("creating canvas: %p (%d,%d) [%d x %d] \n", canvas, x, y, w, h);
*/
canvas->x = x;
canvas->y = y;
canvas->width = w;
canvas->height = h;
if (w > h)
  canvas->size = h;
else
  canvas->size = w;
canvas->active = FALSE;
canvas->resize = TRUE;

/*
canvas->model = sysenv.active_model;
*/
canvas->model = NULL;

sysenv.canvas_list = g_slist_prepend(sysenv.canvas_list, canvas);
}

/**************************************/
/* initialize the OpenGL drawing area */
/**************************************/
void canvas_init(GtkWidget *box)
{
/* create the drawing area */
sysenv.glarea = gtk_drawing_area_new();
gtk_widget_set_gl_capability(sysenv.glarea, sysenv.glconfig, NULL, TRUE, GDK_GL_RGBA_TYPE);
gtk_widget_set_size_request(sysenv.glarea, sysenv.width, sysenv.height);
gtk_box_pack_start(GTK_BOX(box), sysenv.glarea, TRUE, TRUE, 0);

/* init signals */
g_signal_connect(GTK_OBJECT(sysenv.glarea), "expose_event",
                 GTK_SIGNAL_FUNC(canvas_expose), NULL);
g_signal_connect(GTK_OBJECT(sysenv.glarea), "configure_event",
                 GTK_SIGNAL_FUNC(canvas_configure), NULL);

/* TODO - what about the "realize" event??? */

g_signal_connect(GTK_OBJECT(sysenv.glarea), "motion_notify_event",
                 GTK_SIGNAL_FUNC(gui_motion_event), NULL);
g_signal_connect(GTK_OBJECT(sysenv.glarea), "button_press_event",
                 GTK_SIGNAL_FUNC(gui_press_event), NULL);
g_signal_connect(GTK_OBJECT(sysenv.glarea), "button_release_event",
                 GTK_SIGNAL_FUNC(gui_release_event), NULL);
g_signal_connect(GTK_OBJECT(sysenv.glarea), "scroll_event",
                 GTK_SIGNAL_FUNC(gui_scroll_event), NULL);

gtk_widget_set_events(GTK_WIDGET(sysenv.glarea), GDK_EXPOSURE_MASK
                                               | GDK_LEAVE_NOTIFY_MASK
                                               | GDK_BUTTON_PRESS_MASK
                                               | GDK_BUTTON_RELEASE_MASK
                                               | GDK_SCROLL_MASK
                                               | GDK_POINTER_MOTION_MASK
                                               | GDK_POINTER_MOTION_HINT_MASK);

gtk_widget_show(sysenv.glarea);
}

/**************************/
/* table resize primitive */
/**************************/
#define DEBUG_CANVAS_RESIZE 0
void canvas_resize(void)
{
gint i, j, n, rows, cols, width, height;
GSList *list;
struct canvas_pak *canvas;

n = g_slist_length(sysenv.canvas_list);

rows = cols = 1;
switch (n)
  {
  case 2:
    rows = 1;
    cols = 2;
    break;
  case 3:
  case 4:
    rows = 2;
    cols = 2;
    break;
  }

width = sysenv.width / cols;
height = sysenv.height / rows;

#if DEBUG_CANVAS_RESIZE
printf("Splitting (%d, %d) : %d x %d\n", rows, cols, width, height);
#endif

list = sysenv.canvas_list;
for (i=rows ; i-- ; )
  {
  for (j=0 ; j<cols ; j++)
    {
    if (list)
      {
      canvas = list->data;
      canvas->x = j*width;
      canvas->y = i*height;
      canvas->width = width;
      canvas->height = height;

canvas->resize = TRUE;

#if DEBUG_CANVAS_RESIZE
printf(" - canvas (%d, %d) : [%d, %d]\n", i, j, canvas->x, canvas->y);
#endif

      list = g_slist_next(list);    
      }
    }
  }
canvas_shuffle();
}

/*******************************/
/* revert to a single viewport */
/*******************************/
void canvas_single(void)
{
gint i, n;
struct canvas_pak *canvas;

n = g_slist_length(sysenv.canvas_list);

if (n > 1)
  {
  for (i=n-1 ; i-- ; )
    {
    canvas = g_slist_nth_data(sysenv.canvas_list, i);
    sysenv.canvas_list = g_slist_remove(sysenv.canvas_list, canvas);
    }
  }
canvas_resize();
redraw_canvas(SINGLE);
}

/************************************/
/* increase the number of viewports */
/************************************/
void canvas_create(void)
{
gint n;

n = g_slist_length(sysenv.canvas_list);
switch (n)
  {
  case 2:
    canvas_new(0, 0, 0, 0);
  case 1:
  case 0:
    canvas_new(0, 0, 0, 0);
    canvas_resize();
    redraw_canvas(ALL);
    break;
  }
}

/************************************/
/* decrease the number of viewports */
/************************************/
void canvas_delete(void)
{
gint n;
gpointer canvas;

n = g_slist_length(sysenv.canvas_list);

switch (n)
  {
  case 4:
    canvas = sysenv.canvas_list->data;
    sysenv.canvas_list = g_slist_remove(sysenv.canvas_list, canvas);
    canvas_free(canvas);

  case 2:
    canvas = sysenv.canvas_list->data;
    sysenv.canvas_list = g_slist_remove(sysenv.canvas_list, canvas);
    canvas_free(canvas);

/* resize & redraw */
    canvas_resize();
    redraw_canvas(ALL);
    break;
  }
}

/*********************************************/
/* select model at the given canvas position */
/*********************************************/
void canvas_select(gint x, gint y)
{
gint ry;
GSList *list;
struct canvas_pak *canvas;

/* height invert correction */
ry = sysenv.height - y - 1;

for (list=sysenv.canvas_list ; list ; list=g_slist_next(list))
  {
  canvas = list->data;

  if (x >= canvas->x && x < canvas->x+canvas->width)
    {
    if (ry >= canvas->y && ry < canvas->y+canvas->height)
      {
      if (canvas->model)
        {
/* only select if not already active -avoid's deselecting graphs */
        if (canvas->model != sysenv.active_model)
          tree_select_model(canvas->model);
        }
      }
    }
  }
}

/***********************************************/
/* get the canvas a model is drawn in (if any) */
/***********************************************/
gpointer canvas_find(struct model_pak *model)
{
GSList *list;
struct canvas_pak *canvas;

for (list=sysenv.canvas_list ; list ; list=g_slist_next(list))
  {
  canvas = list->data;
  if (canvas->model == model)
    return(canvas);
  }
return(NULL);
}
