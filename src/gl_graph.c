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
#include <string.h>
#include <stdlib.h>
#include <math.h>

#include "gdis.h"
#include "graph.h"
#include "matrix.h"
#include "numeric.h"
#include "opengl.h"
#include "file.h"
#include "parse.h"
#include "dialog.h"
#include "interface.h"
#ifdef __APPLE__
#include <OpenGL/gl.h>
#else
#include <GL/gl.h>
#endif
/* externals */
extern struct sysenv_pak sysenv;
extern gint gl_fontsize;

/******************************************/
/* extract info from graph data structure */
/******************************************/
gchar *graph_treename(gpointer data)
{
struct graph_pak *graph = data;
return(graph->treename);
}

/***************************/
/* free a particular graph */
/***************************/
void graph_free(gpointer data, struct model_pak *model)
{
struct graph_pak *graph = data;

model->graph_list = g_slist_remove(model->graph_list, graph);

if (model->graph_active == graph)
  model->graph_active = NULL;

free_slist(graph->set_list);
g_free(graph->treename);
g_free(graph); 
}

/******************************/
/* free all graphs in a model */
/******************************/
void graph_free_list(struct model_pak *model)
{
GSList *list;
struct graph_pak *graph;

for (list=model->graph_list ; list ; list=g_slist_next(list))
  {
  graph = list->data;

  free_slist(graph->set_list);
  g_free(graph->treename);
  g_free(graph); 
  }
g_slist_free(model->graph_list);

model->graph_list = NULL;
model->graph_active = NULL;
}

/************************/
/* allocate a new graph */
/************************/
/* return an effective graph id */
gpointer graph_new(const gchar *name, struct model_pak *model)
{
static gint n=0;
struct graph_pak *graph;

g_assert(model != NULL);

graph = g_malloc(sizeof(struct graph_pak));

graph->treename = g_strdup_printf("%s_%d", name, n);
graph->wavelength = 0.0;
graph->grafted = FALSE;
graph->xlabel = TRUE;
graph->ylabel = TRUE;
graph->xmin = 0.0;
graph->xmax = 0.0;
graph->ymin = 0.0;
graph->ymax = 0.0;
graph->xticks = 5;
graph->yticks = 5;
graph->size = 0;
graph->select = -1;
graph->select_label = NULL;
graph->set_list = NULL;
graph->type=GRAPH_REGULAR;

/* append to preserve intuitive graph order on the model tree */
model->graph_list = g_slist_append(model->graph_list, graph);

model->graph_active = graph;

n++;

return(graph);
}

/**************/
/* axes setup */
/**************/
void graph_init_y(gdouble *x, struct graph_pak *graph)
{
gdouble ymin, ymax;

g_assert(graph != NULL);

ymin = min(graph->size, x);
ymax = max(graph->size, x);

if (ymin < graph->ymin)
  graph->ymin = ymin;
  
if (ymax > graph->ymin)
  graph->ymax = ymax;
}

/*******************/
/* tree graft flag */
/*******************/
void graph_set_grafted(gint value, gpointer data)
{
struct graph_pak *graph = data;

g_assert(graph != NULL);

graph->grafted = value;
}

/**************************/
/* control label printing */
/**************************/
void graph_set_xticks(gint label, gint ticks, gpointer ptr_graph)
{
struct graph_pak *graph = ptr_graph;

g_assert(graph != NULL);
g_assert(ticks > 1);/* FIXME: ticks shouldn't matter if label=FALSE */

graph->xlabel = label;
graph->xticks = ticks;
}

/**************************/
/* control label printing */
/**************************/
void graph_set_yticks(gint label, gint ticks, gpointer ptr_graph)
{
struct graph_pak *graph = ptr_graph;

g_assert(graph != NULL);
g_assert(ticks > 1);/* FIXME: ticks shouldn't matter if label=FALSE */

graph->ylabel = label;
graph->yticks = ticks;
}

/**********************/
/* special graph data */
/**********************/
void graph_set_wavelength(gdouble wavelength, gpointer ptr_graph)
{
struct graph_pak *graph = ptr_graph;

g_assert(graph != NULL);

graph->wavelength = wavelength;
}

void graph_set_select(gdouble x, gchar *label, gpointer data)
{
gdouble n;
struct graph_pak *graph = data;

g_assert(graph != NULL);

/*
printf("x = %f, min, max = %f, %f\n", x, graph->xmin, graph->xmax);
*/

/* locate the value's position in the data point list */
n = (x - graph->xmin) / (graph->xmax - graph->xmin);
n *= graph->size;

graph->select = (gint) n;
g_free(graph->select_label);
if (label)
  graph->select_label = g_strdup(label);
else
  graph->select_label = NULL;

/*
printf("select -> %d : [0, %d]\n", graph->select, graph->size);
*/
}

/*****************************/
/* add dependent data set(s) */
/*****************************/
void graph_add_data(gint size, gdouble *x, gdouble x1, gdouble x2, gpointer data)
{
gdouble *ptr;
struct graph_pak *graph = data;

g_assert(graph != NULL);

/* try to prevent the user supplying different sized data */
/* TODO - sample in some fashion if different? */
if (graph->size)
  g_assert(graph->size == size);
else
  graph->size = size;

ptr = g_malloc(size*sizeof(gdouble));

memcpy(ptr, x, size*sizeof(gdouble));

graph->xmin = x1;
graph->xmax = x2;
graph_init_y(x, graph);

graph->set_list = g_slist_append(graph->set_list, ptr);
}

/*********************************/
/* add borned (x,y) data (ovhpa) */
/*********************************/
void graph_add_borned_data(gint size,gdouble *x,gdouble x_min,gdouble x_max,gdouble y_min,gdouble y_max,gint type,gpointer data)
{
gdouble *ptr;
struct graph_pak *graph = data;

g_assert(graph != NULL);

/* try to prevent the user supplying different sized data */
/* TODO - sample in some fashion if different? */
if(type==GRAPH_BANDOS){
	/*the first *two* sets of a GRAPH_BANDOS are _not_ the same size*/
	/*also it doesn't make much sense to allow anything on that one.*/
	if(graph->set_list==NULL) graph->size = size;/*only first size matter*/ 
}
if(type==GRAPH_BAND){
	/*the first set of a GRAPH_BAND is _not_ the same size*/
	if(graph->set_list!=NULL) {/*in other cases, proceed as usual*/
		if (graph->size) g_assert(graph->size == size);
		else graph->size = size;
	}
}
if((type!=GRAPH_BAND)&&(type!=GRAPH_BANDOS)){
if (graph->size)
  g_assert(graph->size == size);
else
  graph->size = size;
}

if(type!=GRAPH_USPEX){
ptr = g_malloc(size*sizeof(gdouble));

memcpy(ptr, x, size*sizeof(gdouble));
}else{
ptr = g_malloc(x[0]*sizeof(gdouble));
memcpy(ptr,x,x[0]*sizeof(gdouble));
}

graph->xmin = x_min;
graph->xmax = x_max;
/* because graph_init_y destroy supplied graph limits */
graph->ymin = y_min;
graph->ymax = y_max;
graph->type=type;

graph->set_list = g_slist_append(graph->set_list, ptr);
}
/*********************************************/
/* select a value in a GRAPH_FREQUENCY graph */
/*********************************************/
void graph_frequency_select(gint x, gint y, struct model_pak *model){
/*if selected, gives an IR "peak" value*/
        struct graph_pak *graph;
        struct canvas_pak *canvas;
        gdouble freq_sel;
        gdouble ox,dx;
        gint freq_index=-1;
        int i;
	gdouble *ptr;
	GSList *list;

        /*get the selected frequency, similar to diffraction graphs*/
g_assert(model!=NULL);
        graph=(struct graph_pak *)model->graph_active;
g_assert(graph!=NULL);
        canvas = g_slist_nth_data(sysenv.canvas_list, 0);
g_assert(canvas!=NULL);
        ox=canvas->x + 4*gl_fontsize;
        if (graph->ylabel) ox+=4*gl_fontsize;
	dx=(canvas->width-2.0*ox);
        freq_sel=((x-ox)/dx)*(graph->xmax-graph->xmin);
        freq_sel+=graph->xmin;
for (list=graph->set_list ; list ; list=g_slist_next(list))
  {
        ptr = (gdouble *) list->data;
        for(i=0;i<graph->size;i+=2)
                if((freq_sel>ptr[i]-5.0)&&(freq_sel<ptr[i]+5.0)) freq_index=i;
        if(freq_index>0){/*something was selected*/
		graph->select=freq_index;
		g_free(graph->select_label);
		graph->select_label=g_strdup_printf("[%f]",ptr[freq_index]);
		return;/*no need to continue for all data!*/
        }
  }
}
/***************************************/
/* Select a value in GRAPH_USPEX graph */
/***************************************/
void graph_uspex_select(gint x, gint y, struct model_pak *model){
/*when selected, print value*/
/*TODO: open the corresponding structure file*/
        struct graph_pak *graph;
        struct canvas_pak *canvas;
        gint struct_sel;
        gdouble ox,dx;
	gdouble oy,dy,yf;
	gint yy;
        int i,j;
        gdouble *ptr;
        GSList *list;
/*just a try*/
	FILE *vf;
	gint n_struct;
        /*get the selected structure*/
g_assert(model!=NULL);
        graph=(struct graph_pak *)model->graph_active;
g_assert(graph!=NULL);
        canvas = g_slist_nth_data(sysenv.canvas_list, 0);
g_assert(canvas!=NULL);



/*calculate real data*/
	ox=canvas->x + 4*gl_fontsize;
	if (graph->ylabel) ox+=4*gl_fontsize;
	oy=canvas->y + canvas->height - 2*gl_fontsize;
	if (graph->xlabel) oy-=2*gl_fontsize;
	dy=(canvas->height-8.0*gl_fontsize);
	dx=(canvas->width-2.0*ox);

	struct_sel=(gint)(0.25+((x-ox)/dx)*(graph->size));/*generation*/
	struct_sel--;

/*check if the data point exists*/
if((struct_sel<0)||(struct_sel>graph->size-1)) return;/*not in there*/
j=0;n_struct=0;
for (list=graph->set_list ; list ; list=g_slist_next(list))
  {
	if(j>struct_sel) return;
	if(j<struct_sel) {
		ptr = (gdouble *) list->data;
		n_struct+=(gint)ptr[0];
		j++;
		continue;
	}else{
		/*we are on the good set*/
		ptr = (gdouble *) list->data;
		for(i=1;i<(gint)ptr[0];i++){
			yf = ptr[i];
			yf -= graph->ymin;
			yf /= (graph->ymax - graph->ymin);
			yf *= dy;
			yy = (gint) yf;
			yy *= -1;
			yy += oy;
			if((y>yy-3)&&(y<yy+3)){
				/*we have a hit*/
				graph->select=struct_sel;
				graph->select_2=ptr[i];
				g_free(graph->select_label);
				graph->select_label=g_strdup_printf("[%i,%f]",struct_sel+1,ptr[i]);
				vf=fopen(model->filename, "r");
				if(!vf) return;
				model->cur_frame=n_struct+i-1;
				read_raw_frame(vf,n_struct+i-1,model);
				fclose(vf);
				tree_model_refresh(model);
				model_prep(model);
				return;/*no need to continue for all data!*/
			}
		}
		j++;
	}
  }
}



/************************************/
/* graph data extraction primitives */
/************************************/
gdouble graph_xmin(gpointer data)
{
struct graph_pak *graph = data;
return(graph->xmin);
}

gdouble graph_xmax(gpointer data)
{
struct graph_pak *graph = data;
return(graph->xmax);
}

gint graph_ylabel(gpointer data)
{
struct graph_pak *graph = data;
return(graph->ylabel);
}

gdouble graph_wavelength(gpointer data)
{
struct graph_pak *graph = data;
return(graph->wavelength);
}

gint graph_grafted(gpointer data)
{
struct graph_pak *graph = data;
return(graph->grafted);
}

/****************/
/* draw a graph */
/****************/
#define DEBUG_GRAPH_1D 0
void graph_draw_1d_regular(struct canvas_pak *canvas, struct graph_pak *graph)
{/*previously the only type of graph*/
gint i, x, y, oldx, oldy, ox, oy, sx, sy;
gint flag;
gchar *text;
gdouble *ptr;
gdouble xf, yf, dx, dy;
GSList *list;

#if DEBUG_GRAPH_1D
printf("x range: [%f - %f]\n", graph->xmin, graph->xmax);
printf("y range: [%f - %f]\n", graph->ymin, graph->ymax);
#endif

/* compute origin */
ox = canvas->x + 4*gl_fontsize;
if (graph->ylabel)
  ox += 4*gl_fontsize;

oy = canvas->y + canvas->height - 2*gl_fontsize;
if (graph->xlabel)
  oy -= 2*gl_fontsize;

/* increments for screen drawing */
dy = (canvas->height-8.0*gl_fontsize);
dx = (canvas->width-2.0*ox);

/* axes label colour */
glColor3f(sysenv.render.fg_colour[0],
          sysenv.render.fg_colour[1],
          sysenv.render.fg_colour[2]);
glLineWidth(2.0);


/* x labels */
oldx = ox;
for (i=0 ; i<graph->xticks ; i++)
  {
/* get real index */
  xf = (gdouble) i / (gdouble) (graph->xticks-1);

  x = ox + xf*dx;

  if (graph->xlabel)
    {/*only calculate real value when needed*/
    xf *= (graph->xmax-graph->xmin);
    xf += graph->xmin;
    text = g_strdup_printf("%.2f", xf);
    gl_print_window(text, x-2*gl_fontsize, oy+2*gl_fontsize, canvas);
    g_free(text);
    }

/* axis segment + tick */
  glBegin(GL_LINE_STRIP);
  gl_vertex_window(oldx, oy, canvas);
  gl_vertex_window(x, oy, canvas);
  gl_vertex_window(x, oy+5, canvas);
  glEnd();

  oldx = x;
  }

/* y labels */
oldy = oy;
for (i=0 ; i<graph->yticks ; i++)
  {
/* get screen position */
  yf = (gdouble) i / (gdouble) (graph->yticks-1);
  y = -yf*dy;
  y += oy;

/* label */
  if (graph->ylabel)
    {/*only calculate real value when needed*/
    yf *= (graph->ymax - graph->ymin);
    yf += graph->ymin;
    if (graph->ymax > 999.999999)
      text = g_strdup_printf("%.2e", yf);
    else
      text = g_strdup_printf("%7.2f", yf);
    gl_print_window(text, 0, y-1, canvas);
    g_free(text);
    }

/* axis segment + tick */
  glBegin(GL_LINE_STRIP);
  gl_vertex_window(ox, oldy, canvas);
  gl_vertex_window(ox, y-1, canvas);
  gl_vertex_window(ox-5, y-1, canvas);
  glEnd();

  oldy = y;
  }

/* data drawing colour */
glColor3f(sysenv.render.title_colour[0],
          sysenv.render.title_colour[1],
          sysenv.render.title_colour[2]);
glLineWidth(1.0);


flag = FALSE;
sx = sy = 0;

for (list=graph->set_list ; list ; list=g_slist_next(list))
  {
  ptr = (gdouble *) list->data;

  glBegin(GL_LINE_STRIP);
  for (i=0 ; i<graph->size ; i++)
    {
    xf = (gdouble) i / (gdouble) (graph->size-1);
    x = ox + xf*dx;

/* scale real value to screen coords */
    yf = ptr[i];
    yf -= graph->ymin;
    yf /= (graph->ymax - graph->ymin);
    yf *= dy;

    y = (gint) yf;
    y *= -1;
    y += oy;

if (i == graph->select)
  {
  sx = x;
  sy = y-1;
  flag = TRUE;
  }

/* lift y axis 1 pixel up so y=0 won't overwrite the x axis */
    gl_vertex_window(x, y-1, canvas);
    }
  glEnd();

/* draw peak selector (if any) */
/* TODO - turn off is click outside range? */
  if (flag)
    {
    glEnable(GL_LINE_STIPPLE);
    glLineStipple(1, 0x0303);
    glColor3f(0.9, 0.7, 0.4);
    glLineWidth(2.0);
    glBegin(GL_LINES);
    gl_vertex_window(sx, sy-10, canvas);
    gl_vertex_window(sx, 3*gl_fontsize, canvas);
    glEnd();

    if (graph->select_label)
      {
      gint xoff;

      xoff = (gint) strlen(graph->select_label);
      xoff *= gl_fontsize;
      xoff /= 4;
      gl_print_window(graph->select_label, sx-xoff, 2*gl_fontsize, canvas);
      }
    glDisable(GL_LINE_STIPPLE);
    }
  }


}

void graph_draw_1d_uspex(struct canvas_pak *canvas, struct graph_pak *graph)
{/*similar to 1d but with selectable points*/
gint i, j, x, y, oldx, oldy, ox, oy, sx, sy;
gint flag;
gchar *text;
gdouble *ptr;
gdouble xf, yf, dx, dy;
GSList *list;

/* compute origin */
ox = canvas->x + 4*gl_fontsize;
if (graph->ylabel)
  ox += 4*gl_fontsize;

oy = canvas->y + canvas->height - 2*gl_fontsize;
if (graph->xlabel)
  oy -= 2*gl_fontsize;

/* increments for screen drawing */
dy = (canvas->height-8.0*gl_fontsize);
dx = (canvas->width-2.0*ox);

/* axes label colour */
glColor3f(sysenv.render.fg_colour[0],
          sysenv.render.fg_colour[1],
          sysenv.render.fg_colour[2]);
glLineWidth(2.0);


/* x labels */
oldx = ox;
for (i=0 ; i<graph->xticks ; i++)
  {
/* get real index */
  xf = (gdouble) i / (gdouble) (graph->xticks-1);

  x = ox + xf*dx;

  if (graph->xlabel)
    {/*only calculate real value when needed*/
    xf *= (graph->xmax-graph->xmin);
    xf += graph->xmin;
    text = g_strdup_printf("%i",(gint)xf);
    gl_print_window(text, x-gl_fontsize, oy+2*gl_fontsize, canvas);
    g_free(text);
    }

/* axis segment + tick */
  glBegin(GL_LINE_STRIP);
  gl_vertex_window(oldx, oy, canvas);
  gl_vertex_window(x, oy, canvas);
  gl_vertex_window(x, oy+5, canvas);
  glEnd();

  oldx = x;
  }

/* y labels */
oldy = oy;
for (i=0 ; i<graph->yticks ; i++)
  {
/* get screen position */
  yf = (gdouble) i / (gdouble) (graph->yticks-1);
  y = -yf*dy;
  y += oy;

/* label */
  if (graph->ylabel)
    {/*only calculate real value when needed*/
    yf *= (graph->ymax - graph->ymin);
    yf += graph->ymin;
    if (graph->ymax > 999.999999)
      text = g_strdup_printf("%.2e", yf);
    else
      text = g_strdup_printf("%7.2f", yf);
    gl_print_window(text, 0, y-1, canvas);
    g_free(text);
    }

/* axis segment + tick */
  glBegin(GL_LINE_STRIP);
  gl_vertex_window(ox, oldy, canvas);
  gl_vertex_window(ox, y-1, canvas);
  gl_vertex_window(ox-5, y-1, canvas);
  glEnd();

  oldy = y;
  }

/* data drawing colour */
glColor3f(sysenv.render.title_colour[0],
          sysenv.render.title_colour[1],
          sysenv.render.title_colour[2]);
glLineWidth(1.0);


flag = FALSE;
sx = sy = 0;

i = 1;
for (list=graph->set_list ; list ; list=g_slist_next(list))
  {
  ptr = (gdouble *) list->data;
  for (j=1 ; j<(gint)ptr[0] ; j++)/*the first data is number of structures in a generation */
    {
    xf = (gdouble) (i) / (gdouble) (graph->size);
    x = ox + xf*dx;

/* scale real value to screen coords */
    yf = ptr[j];
    yf -= graph->ymin;
    yf /= (graph->ymax - graph->ymin);
    yf *= dy;

    y = (gint) yf;
    y *= -1;
    y += oy;
  /*draw a rectangle*/
  glBegin(GL_LINE_STRIP);
  gl_vertex_window(x-2, y-2, canvas);
  gl_vertex_window(x+2, y-2, canvas);
  gl_vertex_window(x+2, y+2, canvas);
  gl_vertex_window(x-2, y+2, canvas);
  gl_vertex_window(x-2, y-2, canvas);
  glEnd();


    }
  i++;
  }

if(graph->select_label!=NULL){
	xf = (gdouble) (graph->select+1) / (gdouble) (graph->size);
	sx = ox + xf*dx;
	yf = graph->select_2;
	yf -= graph->ymin; 
	yf /= (graph->ymax - graph->ymin);
	yf *= dy;
	sy = (gint) yf;
	sy *= -1;
	sy += oy;
	glEnable(GL_LINE_STIPPLE);
	glLineStipple(1, 0x0303);
	glColor3f(0.9, 0.7, 0.4);
	glLineWidth(2.0);
	/*horizontal*/
	glBegin(GL_LINES);
	gl_vertex_window(ox, sy, canvas);
	gl_vertex_window(sx, sy, canvas);
	glEnd();
	/*vertical*/
	glBegin(GL_LINES);
	gl_vertex_window(sx, sy, canvas);
	gl_vertex_window(sx, -1*(gint)(dy)+oy, canvas);
	glEnd();
	/*label*/
	gint xoff;
	xoff = (gint) strlen(graph->select_label);
	xoff *= gl_fontsize;
	xoff /= 4;
	gl_print_window(graph->select_label, sx-xoff, 2*gl_fontsize, canvas);
	glDisable(GL_LINE_STIPPLE);
}




}


void graph_draw_1d_frequency(struct canvas_pak *canvas, struct graph_pak *graph)
{
gint i, x, y, oldx, oldy, ox, oy, sx, sy;
gint flag;
gchar *text;
gdouble *ptr;
gdouble xf, yf, dx, dy;
GSList *list;

/* compute origin */
ox = canvas->x + 4*gl_fontsize;
if (graph->ylabel) ox+=4*gl_fontsize;

oy = canvas->y + canvas->height - 2*gl_fontsize;
if (graph->xlabel) oy-=2*gl_fontsize;

/* increments for screen drawing */
dy = (canvas->height-8.0*gl_fontsize);
dx = (canvas->width-2.0*ox);

/* axes label colour */
glColor3f(sysenv.render.fg_colour[0],sysenv.render.fg_colour[1],sysenv.render.fg_colour[2]);
glLineWidth(2.0);

/* x labels */
oldx = ox;
for (i=0 ; i<graph->xticks ; i++)
  {
/* get real index */
  xf = (gdouble) i / (gdouble) (graph->xticks-1);
  x = ox + xf*dx;
  if (graph->xlabel)
    {/*only calculate real value when needed*/
    xf *= (graph->xmax-graph->xmin);
    xf += graph->xmin;
    text = g_strdup_printf("%.2f", xf);
    gl_print_window(text, x-2*gl_fontsize, oy+2*gl_fontsize, canvas);
    g_free(text);
    }
/* axis segment + tick */
  glBegin(GL_LINE_STRIP);
  gl_vertex_window(oldx, oy, canvas);
  gl_vertex_window(x, oy, canvas);
  gl_vertex_window(x, oy+5, canvas);
  glEnd();
  oldx = x;
  }

/* y labels */
oldy = oy;
for (i=0 ; i<graph->yticks ; i++)
  {
/* get screen position */
  yf = (gdouble) i / (gdouble) (graph->yticks-1);
  y = -yf*dy;
  y += oy;
/* label */
  if (graph->ylabel)
    {/*only calculate real value when needed*/
    yf *= (graph->ymax - graph->ymin);
    yf += graph->ymin;
    if (graph->ymax > 999.999999)
      text = g_strdup_printf("%.2e", yf);
    else
      text = g_strdup_printf("%7.2f", yf);
    gl_print_window(text, 0, y-1, canvas);
    g_free(text);
    }
/* axis segment + tick */
  glBegin(GL_LINE_STRIP);
  gl_vertex_window(ox, oldy, canvas);
  gl_vertex_window(ox, y-1, canvas);
  gl_vertex_window(ox-5, y-1, canvas);
  glEnd();
  oldy = y;
  }

/* data drawing colour */
glColor3f(sysenv.render.title_colour[0],sysenv.render.title_colour[1],sysenv.render.title_colour[2]);
glLineWidth(1.0);

flag = FALSE;
sx = sy = 0;

/* plot 1 unit rectangle */
for (list=graph->set_list ; list ; list=g_slist_next(list))
  {
        ptr = (gdouble *) list->data;
        glBegin(GL_LINE_STRIP);
        i=0;
        while(i<graph->size){
                /*each peak is listed here*/
                xf = ptr[i]-1.0;/*1/2 size of rectangle*/
                xf -= graph->xmin;
                xf /= (graph->xmax - graph->xmin);
                x = ox + xf*dx;
                y = oy;/*base*/
if (i==graph->select){
	xf=(ptr[i]-graph->xmin)/(graph->xmax-graph->xmin);
	sx = ox+xf*dx;sy = y-1;flag = TRUE;
}
                gl_vertex_window(x, y-1, canvas);/*go to*/
                yf = ptr[i+1];/*aka intensity*/
                yf -= graph->ymin;
                yf /= (graph->ymax - graph->ymin);
                yf *= dy;
                y = (gint) yf;
                y *= -1;
                y += oy;
                gl_vertex_window(x, y-1, canvas);/*go up*/
                xf = ptr[i]+1.0;/*1/2 size of rectangle*/
                xf -= graph->xmin;
                xf /= (graph->xmax - graph->xmin);
                x = ox + xf*dx;
                gl_vertex_window(x, y-1, canvas);/*move up*/
                y = oy;/*base*/
                gl_vertex_window(x, y-1, canvas);/*go down*/
                i+=2;
        }
        glEnd();
if (flag){
	glEnable(GL_LINE_STIPPLE);
	glLineStipple(1, 0x0303);
	glColor3f(0.9, 0.7, 0.4);
	glLineWidth(2.0);
	glBegin(GL_LINES);
	gl_vertex_window(sx, sy-10, canvas);
	gl_vertex_window(sx, 3*gl_fontsize, canvas);
	glEnd();
	if (graph->select_label)
	{
		gint xoff;
		xoff = strlen(graph->select_label);
		xoff *= gl_fontsize;
		xoff /= 4;
		gl_print_window(graph->select_label, sx-xoff, 2*gl_fontsize, canvas);
	}
	glDisable(GL_LINE_STIPPLE);
}
  }
}
void graph_draw_1d_band(struct canvas_pak *canvas, struct graph_pak *graph)
{
gint i, x, y, oldx, oldy, ox, oy;
gchar *text;
gdouble *ptr;
gdouble xf, yf, dx, dy;
gint shift;
gdouble sz;
gint nbands, nkpoints, efermi, ispin;
gdouble *eval;
GSList *list;

/* compute origin */
ox = canvas->x + 4*gl_fontsize;
if (graph->ylabel) ox+=4*gl_fontsize;

oy = canvas->y + canvas->height - 2*gl_fontsize;
if (graph->xlabel) oy-=2*gl_fontsize;

/* increments for screen drawing */
dy = (canvas->height-8.0*gl_fontsize);
dx = (canvas->width-2.0*ox);

/* axes label colour */
glColor3f(sysenv.render.fg_colour[0],sysenv.render.fg_colour[1],sysenv.render.fg_colour[2]);
glLineWidth(2.0);

/* x labels */
oldx = ox;
for (i=0 ; i<graph->xticks ; i++)
  {
/* get real index */
  xf = (gdouble) i / (gdouble) (graph->xticks-1);
  x = ox + xf*dx;
  if (graph->xlabel)
    {/*only calculate real value when needed*/
    xf *= (graph->xmax-graph->xmin);
    xf += graph->xmin;
    text = g_strdup_printf("%.2f", xf);
    gl_print_window(text, x-2*gl_fontsize, oy+2*gl_fontsize, canvas);
    g_free(text);
    }
/* axis segment + tick */
  glBegin(GL_LINE_STRIP);
  gl_vertex_window(oldx, oy, canvas);
  gl_vertex_window(x, oy, canvas);
  gl_vertex_window(x, oy+5, canvas);
  glEnd();
  oldx = x;
  }

/* y labels */
/* we WANT 0 to be a tick, if nticks>2 */
if(graph->yticks>2) {
	sz=(graph->ymax-graph->ymin)/(graph->yticks);/*size of a tick*/
	shift=(gint)(graph->ymin/sz);
	oldy=oy;
	for (i=0 ; i<=graph->yticks ; i++)
	{/* get real index */
		yf=(i*sz+shift*sz);
		yf -= graph->ymin;
		yf /= (graph->ymax-graph->ymin);
		y = oy - yf*dy;
		if(y<oy-dy) continue;
		if (graph->ylabel)
		{
			yf=(i+shift)*sz;/*label name*/
			if (graph->ymax > 999.999999)
				text = g_strdup_printf("%.2e", yf);
			else
				text = g_strdup_printf("%7.2f", yf);
			gl_print_window(text, 0, y-1, canvas);
			g_free(text);
		}
		/* axis segment + tick */
		glBegin(GL_LINE_STRIP);
		gl_vertex_window(ox, oldy, canvas);
		gl_vertex_window(ox, y-1, canvas);
		gl_vertex_window(ox-5, y-1, canvas);
		glEnd();
		oldy = y;
	}
	/*we need to add a last segment*/
	y=oy-dy;
	glBegin(GL_LINE_STRIP);
	gl_vertex_window(ox, oldy, canvas);
	gl_vertex_window(ox, y-1, canvas);
	glEnd();
}else{/*do the usual ymin,ymax ticks*/
	oldy = oy;
	for (i=0 ; i<graph->yticks ; i++)
	{/* get screen position */
		yf = (gdouble) i / (gdouble) (graph->yticks-1);
		y = -yf*dy;
		y += oy;
		/* label */
		if (graph->ylabel)
		{/*only calculate real value when needed*/
			yf *= (graph->ymax - graph->ymin);
			yf += graph->ymin;
			if (graph->ymax > 999.999999)
				text = g_strdup_printf("%.2e", yf);
			else
				text = g_strdup_printf("%7.2f", yf);
			gl_print_window(text, 0, y-1, canvas);
			g_free(text);
		}
		/* axis segment + tick */
		glBegin(GL_LINE_STRIP);
		gl_vertex_window(ox, oldy, canvas);
		gl_vertex_window(ox, y-1, canvas);
		gl_vertex_window(ox-5, y-1, canvas);
		glEnd();
		oldy = y;
		}
}

/* data drawing colour */
glColor3f(sysenv.render.title_colour[0],sysenv.render.title_colour[1],sysenv.render.title_colour[2]);
glLineWidth(1.0);

/*the first set of GRAPH_BAND data have a special 4 value header*/
list=graph->set_list;
ptr = (gdouble *) list->data;
nbands=(gint)ptr[0];
nkpoints=(gint)ptr[1];
efermi=ptr[2];
ispin=ptr[3];
/*x is then stored as the next nkpoints value */
eval=g_malloc(nkpoints*sizeof(gdouble));/*TODO: eliminate the need of g_malloc/g_free*/
for(i=4;i<(graph->size+4);i++) eval[i-4]=ptr[i];
/*now back to "regular"*/
for (list=g_slist_next(list) ; list ; list=g_slist_next(list))
  {
  ptr = (gdouble *) list->data;

  glBegin(GL_LINE_STRIP);
  for (i=0 ; i<graph->size ; i++)
    {
    xf = eval[i];
	if((xf<graph->xmin)||(xf>graph->xmax)) continue;
    xf -= graph->xmin;
    xf /= (graph->xmax - graph->xmin);
    x = ox + xf*dx;

    yf = ptr[i];
	if((yf<graph->ymin)||(yf>graph->ymax)) continue;
    yf -= graph->ymin;
    yf /= (graph->ymax - graph->ymin);
    yf *= dy;

    y = (gint) yf;
    y *= -1;
    y += oy;

    gl_vertex_window(x, y-1, canvas);
    }
  glEnd();

  }
/*BAND requires a special axis at fermi level*/
  glBegin(GL_LINE_STRIP);
  glColor3f(0.9, 0.7, 0.4);
  glLineWidth(2.0);
  x = ox;
  yf = (0.0-1.0*graph->ymin)/(graph->ymax - graph->ymin);
  y = oy - dy*yf;
  gl_vertex_window(x, y-1, canvas);
  x = ox+dx;
  gl_vertex_window(x, y-1, canvas);
  glEnd();

g_free(eval);
}
void graph_draw_1d_dos(struct canvas_pak *canvas, struct graph_pak *graph)
{
gint i, x, y, oldx, oldy, ox, oy;
gchar *text;
gdouble *ptr;
gdouble xf, yf, dx, dy;
gint shift;
gdouble sz;
GSList *list;

/* compute origin */
ox = canvas->x + 4*gl_fontsize;
if (graph->ylabel) ox += 4*gl_fontsize;
oy = canvas->y + canvas->height - 2*gl_fontsize;
if (graph->xlabel) oy -= 2*gl_fontsize;

/* increments for screen drawing */
dy = (canvas->height-8.0*gl_fontsize);
dx = (canvas->width-2.0*ox);

/* axes label colour */
glColor3f(sysenv.render.fg_colour[0],sysenv.render.fg_colour[1],sysenv.render.fg_colour[2]);
glLineWidth(2.0);

/* x labels */
/* we WANT 0 to be a tick, if nticks>2 */
if(graph->xticks>2) {
	sz=(graph->xmax-graph->xmin)/(graph->xticks);/*size of a tick*/
	shift=(gint)(graph->xmin/sz);
	oldx=ox;
	for (i=0 ; i<=graph->xticks ; i++)
	  {
	/* get real index */
		xf=(i*sz+shift*sz);
		xf -= graph->xmin;
		xf /= (graph->xmax-graph->xmin);
		x = ox + xf*dx;
		if(x>ox+dx) continue;
	  if (graph->xlabel)
	    {
	    xf=(i+shift)*sz;/*label name*/
	    text = g_strdup_printf("%.2f", xf);
	    gl_print_window(text, x-2*gl_fontsize, oy+2*gl_fontsize, canvas);
	    g_free(text);
	    }
	/* axis segment + tick */
	  glBegin(GL_LINE_STRIP);
	  gl_vertex_window(oldx, oy, canvas);
	  gl_vertex_window(x, oy, canvas);
	  gl_vertex_window(x, oy+5, canvas);
	  glEnd();
	  oldx = x;
	  }
	/*we need to add a last segment*/
	x=ox+dx;
	glBegin(GL_LINE_STRIP);
	gl_vertex_window(oldx, oy, canvas);
	gl_vertex_window(x, oy, canvas);
	glEnd();
}else{/*do the usual xmin,xmax ticks*/
	oldx = ox;
	for (i=0 ; i<graph->xticks ; i++)
	  {
	/* get real index */
	  xf = (gdouble) i / (gdouble) (graph->xticks-1);
	  x = ox + xf*dx;
	  if (graph->xlabel)
	    {/*only calculate real value when needed*/
	    xf *= (graph->xmax-graph->xmin);
	    xf += graph->xmin;
	    text = g_strdup_printf("%.2f", xf);
	    gl_print_window(text, x-2*gl_fontsize, oy+2*gl_fontsize, canvas);
	    g_free(text);
	    }
	/* axis segment + tick */
	  glBegin(GL_LINE_STRIP);
	  gl_vertex_window(oldx, oy, canvas);
	  gl_vertex_window(x, oy, canvas);
	  gl_vertex_window(x, oy+5, canvas);
	  glEnd();
	  oldx = x;
	  }
}

/* y labels - useful? */
oldy = oy;
for (i=0 ; i<graph->yticks ; i++)
  {
/* get screen position */
  yf = (gdouble) i / (gdouble) (graph->yticks-1);
  y = -yf*dy;
  y += oy;
/* label */
  if (graph->ylabel)
    {/*only calculate real value when needed*/
    yf *= (graph->ymax - graph->ymin);
    yf += graph->ymin;
    if (graph->ymax > 999.999999)
      text = g_strdup_printf("%.2e", yf);
    else
      text = g_strdup_printf("%7.2f", yf);
    gl_print_window(text, 0, y-1, canvas);
    g_free(text);
    }
/* axis segment + tick */
  glBegin(GL_LINE_STRIP);
  gl_vertex_window(ox, oldy, canvas);
  gl_vertex_window(ox, y-1, canvas);
  gl_vertex_window(ox-5, y-1, canvas);
  glEnd();
  oldy = y;
  }

/* data drawing colour */
glColor3f(sysenv.render.title_colour[0],sysenv.render.title_colour[1],sysenv.render.title_colour[2]);
glLineWidth(1.0);

  for (list=graph->set_list ; list ; list=g_slist_next(list))
        {
        ptr = (gdouble *) list->data;
        glBegin(GL_LINE_STRIP);
        for (i=0 ; i<graph->size ; i+=2) {
                xf = ptr[i];
		if((xf<graph->xmin)||(xf>graph->xmax)) continue;
                xf -= graph->xmin;
                xf /= (graph->xmax - graph->xmin);
                x = ox+xf*dx;
                yf = ptr[i+1];
                yf -= graph->ymin;
                yf /= (graph->ymax - graph->ymin);
                yf *= dy;
                y = (gint) yf;
                y *= -1;
                y += oy;
                gl_vertex_window(x, y-1, canvas);
        }
        glEnd();
  }
/*DOS requires a special axis at e=0 (fermi level)*/
  glBegin(GL_LINE_STRIP);
  glColor3f(0.9, 0.7, 0.4);
  glLineWidth(2.0);
  xf = -1.0*graph->xmin/(graph->xmax - graph->xmin);
  x = ox+xf*dx;
  yf = -1.0*graph->ymin/(graph->ymax - graph->ymin);
  y = oy - dy*yf;
  gl_vertex_window(x, y-1, canvas);
  y = oy - dy;
  gl_vertex_window(x, y-1, canvas);
  glEnd();

}
void graph_draw_1d_bandos(struct canvas_pak *canvas, struct graph_pak *graph)
{
gint i, x, y, oldx, oldy, ox, oy;
gchar *text;
gdouble *ptr,*xptr;
gdouble xf, yf, dx, dy;
gint shift;
gdouble sz;
gdouble xmin=0.;
gdouble xmax=1.;
gint nbands, nkpoints, efermi, ispin;
GSList *list;
GSList *xlist;

/* compute origin */
ox = canvas->x + 4*gl_fontsize;
if (graph->ylabel) ox+=4*gl_fontsize;

oy = canvas->y + canvas->height - 2*gl_fontsize;
if (graph->xlabel) oy-=2*gl_fontsize;

/* increments for screen drawing */
dy = (canvas->height-8.0*gl_fontsize);
dx = (canvas->width-2.0*ox);

/* axes label colour */
glColor3f(sysenv.render.fg_colour[0],sysenv.render.fg_colour[1],sysenv.render.fg_colour[2]);
glLineWidth(2.0);

/* x axis */
oldx = ox;
for (i=0 ; i<graph->xticks ; i++)
  {
/* get real index */
  xf = (gdouble) i / (gdouble) (graph->xticks-1);
  x = ox + xf*dx*0.3;
  if (graph->xlabel)
    {/*only calculate real value when needed*/
    xf *= (xmax-xmin);
    xf += xmin;
    text = g_strdup_printf("%.2f", xf);
    gl_print_window(text, x-2*gl_fontsize, oy+2*gl_fontsize, canvas);
    g_free(text);
    }
/* axis segment + tick */
  glBegin(GL_LINE_STRIP);
  gl_vertex_window(oldx, oy, canvas);
  gl_vertex_window(x, oy, canvas);
  gl_vertex_window(x, oy+5, canvas);
  glEnd();
  oldx = x;
  }

/*second x axis*/
oldx = ox+dx*0.4;
for (i=0 ; i<graph->xticks ; i++)
  {
/* get real index */
  xf = (gdouble) i / (gdouble) (graph->xticks-1);
  x = ox + dx*0.4 + xf*dx*0.6;
  if (graph->xlabel)
    {/*only calculate real value when needed*/
    xf *= (xmax-xmin);
    xf += xmin;
    text = g_strdup_printf("%.2f", xf);
    gl_print_window(text, x-2*gl_fontsize, oy+2*gl_fontsize, canvas);
    g_free(text);
    }
/* axis segment + tick */
  glBegin(GL_LINE_STRIP);
  gl_vertex_window(oldx, oy, canvas);
  gl_vertex_window(x, oy, canvas);
  gl_vertex_window(x, oy+5, canvas);
  glEnd();
  oldx = x;
  }

/* y labels */
/* we WANT 0 to be a tick, if nticks>2 */
if(graph->yticks>2) {
	sz=(graph->ymax-graph->ymin)/(graph->yticks);/*size of a tick*/
	shift=(gint)(graph->ymin/sz);
	oldy=oy;
	for (i=0 ; i<=graph->yticks ; i++)
	{/* get real index */
		yf=(i*sz+shift*sz);
		yf -= graph->ymin;
		yf /= (graph->ymax-graph->ymin);
		y = oy - yf*dy;
		if(y<oy-dy) continue;
		if (graph->ylabel)
		{
			yf=(i+shift)*sz;/*label name*/
			if (graph->ymax > 999.999999)
				text = g_strdup_printf("%.2e", yf);
			else
				text = g_strdup_printf("%7.2f", yf);
			gl_print_window(text, 0, y-1, canvas);
			g_free(text);
		}
		/* axis segment + tick */
		glBegin(GL_LINE_STRIP);
		gl_vertex_window(ox, oldy, canvas);
		gl_vertex_window(ox, y-1, canvas);
		gl_vertex_window(ox-5, y-1, canvas);
		glEnd();
		oldy = y;
	}
	/*we need to add a last segment*/
	y=oy-dy;
	glBegin(GL_LINE_STRIP);
	gl_vertex_window(ox, oldy, canvas);
	gl_vertex_window(ox, y-1, canvas);
	glEnd();
}else{/*do the usual ymin,ymax ticks*/
	oldy = oy;
	for (i=0 ; i<graph->yticks ; i++)
	{/* get screen position */
		yf = (gdouble) i / (gdouble) (graph->yticks-1);
		y = -yf*dy;
		y += oy;
		/* label */
		if (graph->ylabel)
		{/*only calculate real value when needed*/
			yf *= (graph->ymax - graph->ymin);
			yf += graph->ymin;
			if (graph->ymax > 999.999999)
				text = g_strdup_printf("%.2e", yf);
			else
				text = g_strdup_printf("%7.2f", yf);
			gl_print_window(text, 0, y-1, canvas);
			g_free(text);
		}
		/* axis segment + tick */
		glBegin(GL_LINE_STRIP);
		gl_vertex_window(ox, oldy, canvas);
		gl_vertex_window(ox, y-1, canvas);
		gl_vertex_window(ox-5, y-1, canvas);
		glEnd();
		oldy = y;
		}
}

/*second y axis*/
  x=ox+0.4*dx;
  y=oy-dy;
  oldy=oy;
  glBegin(GL_LINE_STRIP);
  gl_vertex_window(x, oldy, canvas);
  gl_vertex_window(x, y-1, canvas);
  gl_vertex_window(x-5, y-1, canvas);
  glEnd();
/* second axis y ticks - NO labels */
/* we WANT 0 to be a tick, if nticks>2 */
if(graph->yticks>2) {
	sz=(graph->ymax-graph->ymin)/(graph->yticks);/*size of a tick*/
	shift=(gint)(graph->ymin/sz);
	oldy=oy;
	for (i=0 ; i<=graph->yticks ; i++)
	{/* get real index */
		yf=(i*sz+shift*sz);
		yf -= graph->ymin;
		yf /= (graph->ymax-graph->ymin);
		y = oy - yf*dy;
		if(y<oy-dy) continue;
		/* axis segment + tick */
		glBegin(GL_LINE_STRIP);
		gl_vertex_window(ox+0.4*dx, oldy, canvas);
		gl_vertex_window(ox+0.4*dx, y-1, canvas);
		gl_vertex_window(ox+0.4*dx-5, y-1, canvas);
		glEnd();
		oldy = y;
	}
	/*we need to add a last segment*/
	y=oy-dy;
	glBegin(GL_LINE_STRIP);
	gl_vertex_window(ox+0.4*dx, oldy, canvas);
	gl_vertex_window(ox+0.4*dx, y-1, canvas);
	glEnd();
}else{/*do the usual ymin,ymax ticks*/
	oldy = oy;
	for (i=0 ; i<graph->yticks ; i++)
	{/* get screen position */
		yf = (gdouble) i / (gdouble) (graph->yticks-1);
		y = -yf*dy;
		y += oy;
		/* axis segment + tick */
		glBegin(GL_LINE_STRIP);
		gl_vertex_window(ox+0.4*dx, oldy, canvas);
		gl_vertex_window(ox+0.4*dx, y-1, canvas);
		gl_vertex_window(ox+0.4*dx-5, y-1, canvas);
		glEnd();
		oldy = y;
		}
}

/* data drawing colour */
glColor3f(sysenv.render.title_colour[0],sysenv.render.title_colour[1],sysenv.render.title_colour[2]);
glLineWidth(1.0);

/*first set is the DOS, which we need to put on [0,0.3] of x, rotated by 90 degrees*/
list=graph->set_list;
ptr = (gdouble *) list->data;
glBegin(GL_LINE_STRIP);
        for (i=0 ; i<graph->size ; i+=2) {
                xf = ptr[i+1];
                xf -= graph->xmin;
                xf /= (graph->xmax - graph->xmin);
                x = ox+xf*dx*0.3;
                yf = ptr[i];
		if((yf<graph->ymin)||(yf>graph->ymax)) continue;
                yf -= graph->ymin;
                yf /= (graph->ymax - graph->ymin);
                yf *= dy;
                y = (gint) yf;
                y *= -1;
                y += oy;
                gl_vertex_window(x, y-1, canvas);
        }
  glEnd();
/*BANDOS requires two special axis at fermi level*/
  glBegin(GL_LINE_STRIP);
  glColor3f(0.9, 0.7, 0.4);
  glLineWidth(2.0);
  x = ox;
  yf = (0.0-1.0*graph->ymin)/(graph->ymax - graph->ymin);
  y = oy - dy*yf;
  gl_vertex_window(x, y-1, canvas);
  x = ox+dx*0.3;
  gl_vertex_window(x, y-1, canvas);
  glEnd();
/* data drawing colour */
glColor3f(sysenv.render.title_colour[0],sysenv.render.title_colour[1],sysenv.render.title_colour[2]);
glLineWidth(1.0);
/* next set contain the header of bands + x values */
xlist=g_slist_next(list);
xptr = (gdouble *) xlist->data;
nbands=(gint)xptr[0];
nkpoints=(gint)xptr[1];
efermi=xptr[2];
ispin=xptr[3];
xmin=xptr[4];
xmax=xptr[nkpoints+3];
/*now to BAND plotting*/
for (list=g_slist_next(xlist) ; list ; list=g_slist_next(list))
  {
  ptr = (gdouble *) list->data;

  glBegin(GL_LINE_STRIP);
  for (i=0 ; i<nkpoints ; i++)
    {
    xf = xptr[i+4];
    xf -= xmin;
    xf /= (xmax - xmin);
    x = ox + dx*0.4 + xf*dx*0.6;

    yf = ptr[i];
	if((yf<graph->ymin)||(yf>graph->ymax)) continue;
    yf -= graph->ymin;
    yf /= (graph->ymax - graph->ymin);
    yf *= dy;

    y = (gint) yf;
    y *= -1;
    y += oy;

    gl_vertex_window(x, y-1, canvas);
    }
  glEnd();

  }
/*BANDOS requires two special axis at fermi level*/
  glBegin(GL_LINE_STRIP);
  glColor3f(0.9, 0.7, 0.4);
  glLineWidth(2.0);
  x = ox+dx*0.4;
  yf = (0.0-1.0*graph->ymin)/(graph->ymax - graph->ymin);
  y = oy - dy*yf;
  gl_vertex_window(x, y-1, canvas);
  x = ox+dx;
  gl_vertex_window(x, y-1, canvas);
  glEnd();

}
void graph_draw_1d(struct canvas_pak *canvas, struct graph_pak *graph)
{/* NEW graph selector */
	switch(graph->type){
	case GRAPH_REGULAR:
		graph_draw_1d_regular(canvas,graph);
		break;
	case GRAPH_FREQUENCY:
		graph_draw_1d_frequency(canvas,graph);
		break;
	case GRAPH_BAND:
		graph_draw_1d_band(canvas,graph);
		break;
	case GRAPH_DOS:
		graph_draw_1d_dos(canvas,graph);
		break;
	case GRAPH_BANDOS:
		graph_draw_1d_bandos(canvas,graph);
		break;
	case GRAPH_USPEX:
		graph_draw_1d_uspex(canvas,graph);
		break;
	case GRAPH_UNKNOWN:
	default:
		fprintf(stderr,"ERROR: graph type unkonwn");
	}
}

/*********************************/
/* init for OpenGL graph drawing */
/*********************************/
void graph_draw(struct canvas_pak *canvas, struct model_pak *model)
{
/* checks */
g_assert(canvas != NULL);
g_assert(model != NULL);
if (!g_slist_find(model->graph_list, model->graph_active))
  return;

/* init drawing model */
glDisable(GL_LIGHTING);
glDisable(GL_LINE_STIPPLE);
glDisable(GL_DEPTH_TEST);

glEnable(GL_BLEND);
glEnable(GL_LINE_SMOOTH);
glEnable(GL_POINT_SMOOTH);

glDisable(GL_COLOR_LOGIC_OP);

/* draw the appropriate type */
graph_draw_1d(canvas, model->graph_active);
}

/*******************/
/* graph exporting */
/*******************/
void graph_write(gchar *name, gpointer ptr_graph)
{
gint i;
gchar *filename;
gdouble x, *y;
GSList *list;
FILE *fp;
struct graph_pak *graph = ptr_graph;

/* checks */
g_assert(graph != NULL);
g_assert(name != NULL);

/* init */
filename = g_build_filename(sysenv.cwd, name, NULL);
fp = fopen(filename,"wt");
if (!fp)
  return;

/* write */
for (list=graph->set_list ; list ; list=g_slist_next(list))
  {
  y = (gdouble *) list->data;

  for (i=0 ; i<graph->size ; i++)
    {
    x = (gdouble) i / (gdouble) graph->size;
    x *= (graph->xmax - graph->xmin);
    x += graph->xmin;

    fprintf(fp, "%f %f\n", x, y[i]);
    }
  }
fclose(fp);
g_free(filename);
}

/*******************/
/* graph importing */
/*******************/
#define GRAPH_SIZE_MAX 4
void graph_read(gchar *filename)
{
gint i, j, n, num_tokens;
gchar *fullpath, *line, **buff;
gdouble xstart, xstop, x[GRAPH_SIZE_MAX], *y;
GArray *garray;
gpointer graph;
struct model_pak *model;
FILE *fp;

/* TODO - close dialog */
dialog_destroy_type(FILE_SELECT);

/* read / create model / plot etc */
model = sysenv.active_model;
if (!model)
  {
  edit_model_create();
  model = sysenv.active_model;
  }

garray = g_array_new(FALSE, FALSE, sizeof(gdouble));

/* read */
fullpath = g_build_filename(sysenv.cwd, filename, NULL);
fp = fopen(fullpath, "rt");
line = file_read_line(fp);
n=0;
while (line)
  {
  buff = tokenize(line, &num_tokens);
  g_free(line);

/* get all numbers */
  j = 0;
  for (i=0 ; i<num_tokens ; i++)
    {
    if (str_is_float(*(buff+i)))
      {
      if (j < GRAPH_SIZE_MAX)
        x[j++] = str_to_float(*(buff+i));
      }
    }

/* need x & y (2 items) minimum */
  if (j > 1)
    {
    if (!n)
      xstart = x[0];
    g_array_append_val(garray, x[1]);
    n++;
    }

  g_strfreev(buff);
  line = file_read_line(fp);
  }
xstop = x[0];
fclose(fp);

y = g_malloc(n * sizeof(gdouble));
for (i=0 ; i<n ; i++)
  y[i] = g_array_index(garray, gdouble, i);

graph = graph_new("Graph", model);
graph_add_data(n, y, xstart, xstop, graph);

g_array_free(garray, TRUE);
g_free(y);
}

