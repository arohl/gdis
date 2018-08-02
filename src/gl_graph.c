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
#include "model.h"
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
/*for uspex structure -> should be in interface?*/
#include "file_uspex.h"
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

graph->treename = g_strdup_printf("%s", name);
graph->treenumber = n;
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

if((type==GRAPH_USPEX)||(type==GRAPH_USPEX_BEST)){
ptr = g_malloc(1+((gint)x[0])*sizeof(gdouble));
memcpy(ptr,x,(1+(gint)x[0])*sizeof(gdouble));
}else{
ptr = g_malloc(size*sizeof(gdouble));
memcpy(ptr, x, size*sizeof(gdouble));
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
/*when selected, print value and open the corresponding structure file*/
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

/*USPEX_BEST: simple search*/
if(graph->type==GRAPH_USPEX_BEST){
list=graph->set_list;/*there should be only one list*/
ptr = (gdouble *) list->data;
i = struct_sel + 1;
yf = ptr[i];
yf -= graph->ymin;
yf /= (graph->ymax - graph->ymin);
yf *= dy;
yy = (gint) yf;
yy *= -1;
yy += oy;
if((y>yy-3)&&(y<yy+3)){
	uspex_output_struct *uspex_calc=model->uspex;
	/*we have a hit*/
	graph->select=struct_sel;
	graph->select_2=ptr[i];
	g_free(graph->select_label);
	graph->select_label=g_strdup_printf("[%i,%f]",struct_sel+1,ptr[i]);
	vf=fopen(model->filename, "r");
	if(!vf) return;
	model->cur_frame=(*uspex_calc).best_ind[struct_sel+1];/*this is the structure number*/
	read_raw_frame(vf,model->cur_frame,model);
	fclose(vf);
	tree_model_refresh(model);
	model_prep(model);
	return;/*no need to continue for all data!*/
}
}else{
/*USPEX: more complex search*/
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
		for(i=1;i<=(gint)ptr[0];i++){
			yf = ptr[i];
			yf -= graph->ymin;
			yf /= (graph->ymax - graph->ymin);
			yf *= dy;
			yy = (gint) yf;
			yy *= -1;
			yy += oy;
			if((y>yy-3)&&(y<yy+3)){
				uspex_output_struct *uspex_calc=model->uspex;
				/*we have a hit*/
				graph->select=struct_sel;
				graph->select_2=ptr[i];
				g_free(graph->select_label);
				graph->select_label=g_strdup_printf("[%i,%f]",struct_sel+1,ptr[i]);
				vf=fopen(model->filename, "r");
				if(!vf) return;
				if((*uspex_calc).method==US_CM_META) n_struct++;
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
}
/******************************************/
/* Select a value in GRAPH_USPEX_2D graph */
/******************************************/
void graph_uspex_2d_select(gint x, gint y, struct model_pak *model){
/*when selected, print value and open the corresponding structure file*/
        struct graph_pak *graph;
        struct canvas_pak *canvas;
        gdouble ox,dx, xf;
        gdouble oy,dy, yf;
        gint xx,yy;
        int i;
	gdouble *xval;
        gdouble *yval;
	gdouble *tag;
        GSList *list;
        FILE *vf;
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

/*first set is tags*/
list=graph->set_list;
tag = (gdouble *) list->data;
list=g_slist_next(list);
/*second set is x*/
xval = (gdouble *) list->data;
list=g_slist_next(list);

for ( ; list ; list=g_slist_next(list)){
  i=0;
  yval = (gdouble *) list->data;
  while(i<graph->size){
        xf = xval[i];
	xf -= graph->xmin;
	xf /= (graph->xmax - graph->xmin);
	xx = ox + xf*dx;
	if((x>xx-3)&&(x<xx+3)){/*x is good*/
		yf = yval[i];
		yf -= graph->ymin;
		yf /= (graph->ymax - graph->ymin);
		yf *= dy;
		yy = (gint) yf;
		yy *= -1;
		yy += oy;
		if((y>yy-3)&&(y<yy+3)){/*y also good*/
			graph->select=i;
			graph->select_2=yval[i];
			g_free(graph->select_label);
			graph->select_label=g_strdup_printf("[%f,%f]",xval[i],yval[i]);
			vf=fopen(model->filename, "r");
			if(!vf) return;
			model->cur_frame=tag[i];
			read_raw_frame(vf,tag[i],model);
			fclose(vf);
			tree_model_refresh(model);
			model_prep(model);
			return;/*no need to continue for all data!*/
		}
	}
	i++;
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
/****************************************************/
/* draw a graph using the new, general graph system */
/****************************************************/
void graph_draw_new(struct canvas_pak *canvas, struct graph_pak *graph){
/*draw a graph using the new, more general graph system*/
gint i, j, x, y, oldx, oldy, ox, oy, sx, sy;
gint flag;
gchar *text;
gdouble *ptr;
gdouble xf, yf, dx, dy;
GSList *list;
/*specific*/
gint size;
gint shift;
gint nkpoints;
gdouble sz;
gdouble *xval=NULL;
gdouble xmin=0.;
gdouble xmax=0.;
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
if((graph->type==GRAPH_DOS)&&(graph->xticks>2)){
	/* we WANT 0 to be a tick, if nticks>2 */
        sz=(graph->xmax-graph->xmin)/(graph->xticks);/*size of a tick*/
        shift=(gint)(graph->xmin/sz);
        oldx=ox;
        for (i=0 ; i<=graph->xticks ; i++){
        /* get real index */
                xf=(i*sz+shift*sz);
                xf -= graph->xmin;
                xf /= (graph->xmax-graph->xmin);
                x = ox + xf*dx;
                if(x>ox+dx) continue;
          if (graph->xlabel){
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
  for (i=0 ; i<graph->xticks ; i++){
  /* get real index */
  xf = (gdouble) i / (gdouble) (graph->xticks-1);
  if(graph->type==GRAPH_BANDOS){
    x = ox + xf*dx*0.3;
  }else{
    x = ox + xf*dx;
  }
  if (graph->xlabel){
	/*only calculate real value when needed*/
	if(graph->type!=GRAPH_BANDOS){
	    xf *= (graph->xmax-graph->xmin);
	    xf += graph->xmin;
	}
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
/* another x axis */
if(graph->type==GRAPH_BANDOS){
  /*second x axis*/
  oldx = ox+dx*0.4;
  for (i=0 ; i<graph->xticks ; i++){
  /* get real index */
  xf = (gdouble) i / (gdouble) (graph->xticks-1);
  x = ox + dx*0.4 + xf*dx*0.6;
  if (graph->xlabel){
    /*only calculate real value when needed*/
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
/* y labels */
if(((graph->type==GRAPH_BANDOS)||(graph->type==GRAPH_BAND))&&(graph->yticks>2)){
  /* we WANT 0 to be a tick, if nticks>2 */
  sz=(graph->ymax-graph->ymin)/(graph->yticks);/*size of a tick*/
  shift=(gint)(graph->ymin/sz);
  oldy=oy;
  for (i=0 ; i<=graph->yticks ; i++){
    /* get real index */
    yf=(i*sz+shift*sz);
    yf -= graph->ymin;
    yf /= (graph->ymax-graph->ymin);
    y = oy - yf*dy;
    if(y<oy-dy) continue;
    if (graph->ylabel){
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
  for (i=0 ; i<graph->yticks ; i++){
    /* get screen position */
    yf = (gdouble) i / (gdouble) (graph->yticks-1);
    y = -yf*dy;
    y += oy;
    /* label */
    if (graph->ylabel){
	/*only calculate real value when needed*/
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
/* another y axis */
if(graph->type==GRAPH_BANDOS){
  x=ox+0.4*dx;
  y=oy-dy;
  oldy=oy;
  /* second axis y ticks - NO labels */
  /* we WANT 0 to be a tick, if nticks>2 */
  if(graph->yticks>2){
    sz=(graph->ymax-graph->ymin)/(graph->yticks);/*size of a tick*/
    shift=(gint)(graph->ymin/sz);
    oldy=oy;
    for (i=0 ; i<=graph->yticks ; i++){
	/* get real index */
	yf=(i*sz+shift*sz);
	yf -= graph->ymin;
	yf /= (graph->ymax-graph->ymin);
	y = oy - yf*dy;
	if(y<oy-dy) continue;
	/* axis segment + tick */
	glBegin(GL_LINE_STRIP);
	gl_vertex_window(x, oldy, canvas);
	gl_vertex_window(x, y-1, canvas);
	gl_vertex_window(x-5, y-1, canvas);
	glEnd();
	oldy = y;
    }
    /*we need to add a last segment*/
    y=oy-dy;
    glBegin(GL_LINE_STRIP);
    gl_vertex_window(x, oldy, canvas);
    gl_vertex_window(x, y-1, canvas);
    glEnd();
  }else{/*do the usual ymin,ymax ticks*/
    oldy = oy;
    glBegin(GL_LINE_STRIP);
    gl_vertex_window(x, oldy, canvas);
    gl_vertex_window(x, y-1, canvas);
    gl_vertex_window(x-5, y-1, canvas);
    glEnd();
    for (i=0 ; i<graph->yticks ; i++){
	/* get screen position */
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
}
/* data drawing colour */
glColor3f(sysenv.render.title_colour[0],sysenv.render.title_colour[1],sysenv.render.title_colour[2]);
glLineWidth(1.0);
flag = FALSE;
sx = sy = 0;
list=graph->set_list;
size=graph->size;
if(graph->type==GRAPH_BANDOS){
/*first set is the DOS, to put on [0,0.3] of x, rotated by 90 degrees*/
	ptr = (gdouble *) list->data;
	glBegin(GL_LINE_STRIP);
	for (i=0 ; i<size ; i+=2) {
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
/* next set contain nkpoints, new xmin xmax, and x values */
	list=g_slist_next(list);
	ptr = (gdouble *) list->data;
/*we do not read the full header (yet)*/
	nkpoints=(gint)ptr[1];
	xmin=ptr[4];
	xmax=ptr[nkpoints+3];
/*setup x values*/
	xval=&(ptr[4]);
	size=nkpoints;
	list=g_slist_next(list);
}
if(graph->type==GRAPH_BAND){
/*the first set of GRAPH_BAND data contain X values*/
/*we do not read the full header (yet)*/
	ptr = (gdouble *) list->data;
	xval=&(ptr[4]);
	list=g_slist_next(list);
}
if(graph->type==GRAPH_USPEX_2D){
/*first set is the tags, ignore*/
list=g_slist_next(list);
/*second set is x*/
	ptr = (gdouble *) list->data;
	xval=&(ptr[0]);
	list=g_slist_next(list);
}
j=0;
for ( ; list ; list=g_slist_next(list))
  {
  ptr = (gdouble *) list->data;
  oldx=-1;
  oldy=-1;
  i=0;
  if(graph->type==GRAPH_USPEX) size=(gint)ptr[0];
  while(i<size)
    {
/*get real values*/
switch (graph->type){
	case GRAPH_USPEX_BEST:
	xf = (gdouble) (i);
	yf = ptr[i+1];
	break;
	case GRAPH_USPEX:
	xf = (gdouble) (j);
	yf = ptr[i+1];
	break;
	case GRAPH_FREQUENCY:
	case GRAPH_DOS:
	xf = ptr[i];
	yf = ptr[i+1];
	break;
	case GRAPH_USPEX_2D:
	case GRAPH_BANDOS:
	case GRAPH_BAND:
	xf = xval[i];
	yf = ptr[i];
	break;
/*complete until full*/
	case GRAPH_REGULAR:
	default:
	xf = (gdouble) i;
	yf = ptr[i];
}
/*new: skip drawing for points outside of graph extrema (also interrupt line when line drawing)*/
if((xf<graph->xmin)||(xf>graph->xmax)||(yf<graph->ymin)||(yf>graph->ymax)) {
	/*update index*/
	switch(graph->type){
	  case GRAPH_FREQUENCY:
	  case GRAPH_DOS:
		i++;/* ie twice ++ */
	  case GRAPH_USPEX_2D:
	  case GRAPH_USPEX_BEST:
	  case GRAPH_USPEX:
	  case GRAPH_BAND:
	  case GRAPH_BANDOS:
	  case GRAPH_REGULAR:
	  default:
		i++;
	}
	oldx=-1;
	continue;
}
/*calculate screen values*/
switch (graph->type){
	case GRAPH_FREQUENCY:
	/*we can draw impulse directly*/
	glBegin(GL_LINE_STRIP);
	xf = ptr[i];/*1/2 size of rectangle*/
	xf -= graph->xmin;
	xf /= (graph->xmax - graph->xmin);
	x = ox + xf*dx;
	y = oy;/*base*/
	/*selector for frequency*/
	if (i==graph->select){
		sx=x;sy=y-1;flag = TRUE;
	}
	gl_vertex_window(x-1,y-1,canvas);/*go to*/
	yf -= graph->ymin;
	yf /= (graph->ymax - graph->ymin);
	yf *= dy;
	y = (gint) yf;
	y *= -1;
	y += oy;
	gl_vertex_window(x-1,y-1, canvas);/*go up*/
	gl_vertex_window(x+1,y-1, canvas);/*move right*/
	y = oy;/*base*/
	gl_vertex_window(x+1,y-1, canvas);/*go down*/
	glEnd();
	break;
	case GRAPH_DOS:
	case GRAPH_BAND:
	xf -= graph->xmin;
	xf /= (graph->xmax - graph->xmin);
	x = ox + xf*dx;
	yf -= graph->ymin;
	yf /= (graph->ymax - graph->ymin);
	yf *= dy;
	y = (gint) yf;
	y *= -1;
	y += oy;
	break;
	case GRAPH_BANDOS:
	xf -= xmin;
	xf /= (xmax - xmin);
	x = ox + dx*0.4 + xf*dx*0.6;
	yf -= graph->ymin;
	yf /= (graph->ymax - graph->ymin);
	yf *= dy;
	y = (gint) yf;
	y *= -1;
	y += oy;
	break;
	case GRAPH_USPEX_BEST:
	case GRAPH_USPEX:
	xf += 1.;
	xf /= (gdouble) (graph->size);/*NOT size*/
	x = ox + xf*dx;
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
	break;
	case GRAPH_USPEX_2D:
	xf -= graph->xmin;
	xf /= (graph->xmax - graph->xmin);
	x = ox + xf*dx;
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
	if((i == graph->select)&&(ptr[i] == graph->select_2)){
		sx=x;
		sy=y-1;
		flag=TRUE;
	}
	break;
	case GRAPH_REGULAR:
	default:
	xf = (gdouble) i / (gdouble) (graph->size-1);
	x = ox + xf*dx;
	yf -= graph->ymin;
	yf /= (graph->ymax - graph->ymin);
	yf *= dy;
	y = (gint) yf;
	y *= -1;
	y += oy;
	/*selector will be only within graph limits*/
	if (i == graph->select){
		sx=x;
		sy=y-1;
		flag=TRUE;
	}
}
if(oldx!=-1){
	/*plot LINE data*/
  switch(graph->type){
	case GRAPH_USPEX:
	case GRAPH_USPEX_2D:
	break;
	case GRAPH_FREQUENCY:
		/*we should never reach here*/
	break;
	case GRAPH_USPEX_BEST:
	glBegin(GL_LINE_STRIP);
	gl_vertex_window(oldx,oldy-1,canvas);
	gl_vertex_window(x,y-1,canvas);
	glEnd();
	break;
	case GRAPH_DOS:
	case GRAPH_BAND:
	case GRAPH_REGULAR:
	default:
	/*draw a line*/
	glBegin(GL_LINE_STRIP);
	/* lift y axis 1 pixel up so y=0 won't overwrite the x axis */
	gl_vertex_window(oldx, oldy-1, canvas);
	gl_vertex_window(x, y-1, canvas);
	glEnd();
  }
}
	oldx=x;
	oldy=y;
/*update index*/
	switch(graph->type){
	case GRAPH_FREQUENCY:
	case GRAPH_DOS:
		i++;/* ie twice ++ */
	case GRAPH_USPEX:
	case GRAPH_USPEX_BEST:
	case GRAPH_USPEX_2D:
	case GRAPH_BAND:
	case GRAPH_BANDOS:
	case GRAPH_REGULAR:
	default:
		i++;
	}
    }
    j++;
  }
/*outside of loop*/
switch (graph->type){
	/*plot add-axis*/
	case GRAPH_DOS:
	/* need 0 y axis at Fermi level */
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
	break;
	case GRAPH_BAND:
	/* need 0 x axis at Fermi level */
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
	break;
	case GRAPH_BANDOS:
	/*BANDOS requires two special 0 y axis at fermi level*/
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
	/* second axis */
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
	break;
	case GRAPH_USPEX_BEST:
	case GRAPH_USPEX_2D:/*TODO: add a 0.5 composition axis?*/
	case GRAPH_USPEX:
	case GRAPH_FREQUENCY:
	case GRAPH_REGULAR:
	default:
	break;
}

/* draw peak selector */
/*regular, frequency*/
  if (flag) {
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
/* uspex, uspex_best */
  if(((graph->type==GRAPH_USPEX)||(graph->type==GRAPH_USPEX_BEST))&&(graph->select_label!=NULL)){
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

/****************/
/* generic call */
/****************/
void graph_draw_1d(struct canvas_pak *canvas, struct graph_pak *graph)
{/* NEW graph selector */
	switch(graph->type){
	case GRAPH_REGULAR:
	case GRAPH_FREQUENCY:
	case GRAPH_BAND:
	case GRAPH_DOS:
	case GRAPH_BANDOS:
	case GRAPH_USPEX:
	case GRAPH_USPEX_BEST:
	case GRAPH_USPEX_2D:
		graph_draw_new(canvas,graph);
		break;
	default:
		graph_draw_new(canvas,graph);
		fprintf(stderr,"WARNING: graph type unknown!\n");
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
gdouble xstart=0.;
gdouble xstop, x[GRAPH_SIZE_MAX], *y;
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

