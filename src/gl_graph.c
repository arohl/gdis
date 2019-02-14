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
/*NEW: graph_controls*/
graph->title=NULL;
graph->sub_title=NULL;
graph->x_title=NULL;
graph->y_title=NULL;

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
if(label) g_assert(ticks > 1);

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
if(label) g_assert(ticks > 1);

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
/********************************/
/* regular selection of 1D data */
/********************************/
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

ptr = g_malloc(size*sizeof(gdouble));
memcpy(ptr, x, size*sizeof(gdouble));

graph->xmin = x_min;
graph->xmax = x_max;
/* because graph_init_y destroy supplied graph limits */
graph->ymin = y_min;
graph->ymax = y_max;
graph->type=type;

graph->set_list = g_slist_append(graph->set_list, ptr);
}
/******************************/
/* NEW - dat_graph: set title */
/******************************/
void dat_graph_set_title(const gchar *title,gpointer pgraph){
struct graph_pak *graph = pgraph;
g_assert(graph != NULL);
/*duplicate string*/
graph->title=g_strdup(title);
}
/**********************************/
/* NEW - dat_graph: set sub-title */
/**********************************/
void dat_graph_set_sub_title(const gchar *title,gpointer pgraph){
struct graph_pak *graph = pgraph;
g_assert(graph != NULL);
/*duplicate string*/
graph->sub_title=g_strdup(title);
}
/********************************/
/* NEW - dat_graph: set x_title */
/********************************/
void dat_graph_set_x_title(const gchar *x_title,gpointer pgraph){
struct graph_pak *graph = pgraph;
g_assert(graph != NULL);
/*duplicate string*/
graph->x_title=g_strdup(x_title);
}
/********************************/
/* NEW - dat_graph: set y_title */
/********************************/
void dat_graph_set_y_title(const gchar *y_title,gpointer pgraph){
struct graph_pak *graph = pgraph;
g_assert(graph != NULL);
/*duplicate string*/
graph->y_title=g_strdup(y_title);
}
/*************************/
/* NEW - dat_graph: set x*/
/*************************/
void dat_graph_set_x(g_data_x dx,gpointer pgraph){
struct graph_pak *graph = pgraph;
int i;
g_data_x *px;

g_assert(graph != NULL);

/*duplicate and register*/
px=g_malloc(sizeof(g_data_x));
g_assert(px != NULL);
px->x_size = dx.x_size;
px->x = g_malloc(dx.x_size*sizeof(gdouble));
//memcpy(px->x,dx.x,dx.x_size*sizeof(gdouble));
for(i=0;i<dx.x_size;i++) px->x[i]=dx.x[i];

graph->size=0;/*graph size hold the number of y arrays!*/
graph->set_list = g_slist_append(graph->set_list, px);
}
/**************************/
/* NEW - dat_graph: add y */
/**************************/
void dat_graph_add_y(g_data_y dy,gpointer pgraph){
struct graph_pak *graph = pgraph;
g_data_y *py;

g_assert(graph != NULL);

if(dy.y_size==0) {
	/*empty size data set are OK*/
	py=g_malloc(sizeof(g_data_y));
	py->y_size=0;
	py->y=NULL;
	py->idx=NULL;
	py->type=GRAPH_REGULAR;
	py->symbol=NULL;
	py->mixed_symbol=FALSE;
	py->line=GRAPH_LINE_NONE;
	py->color=GRAPH_COLOR_DEFAULT;
	graph->size++;/*graph size hold the number of y arrays!*/
	graph->set_list = g_slist_append(graph->set_list, py);
	return;
}
/*duplicate and register*/
py=g_malloc(sizeof(g_data_y));
g_assert(py != NULL);
py->y_size=dy.y_size;
py->y=g_malloc(dy.y_size*sizeof(gdouble));
memcpy(py->y,dy.y,dy.y_size*sizeof(gdouble));
if(dy.idx==NULL) py->idx=NULL;
else{
	py->idx=g_malloc(dy.y_size*sizeof(gint32));
	memcpy(py->idx,dy.idx,dy.y_size*sizeof(gint32));
}
if(dy.symbol==NULL) py->symbol=NULL;
else{
	py->symbol=g_malloc(dy.y_size*sizeof(graph_symbol));
	memcpy(py->symbol,dy.symbol,dy.y_size*sizeof(graph_symbol));
}
py->mixed_symbol=dy.mixed_symbol;
py->line=dy.line;
py->color=dy.color;
py->type=dy.type;

graph->size++;/*graph size hold the number of y arrays!*/
graph->set_list = g_slist_append(graph->set_list, py);
}
/*******************************/
/* NEW - dat_graph: set limits */
/*******************************/
void dat_graph_set_limits(gdouble x_min,gdouble x_max,gdouble y_min,gdouble y_max,gpointer pgraph){
struct graph_pak *graph = pgraph;

g_assert(graph != NULL);

graph->xmin = x_min;
graph->xmax = x_max;
graph->ymin = y_min;
graph->ymax = y_max;
}
/*******************************/
/* NEW - set global graph type */
/*******************************/
void dat_graph_set_type(graph_type type,gpointer pgraph){
struct graph_pak *graph = pgraph;

g_assert(graph != NULL);

graph->type=type;
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
/***********************************/
/* NEW - dat_graph: select a value */
/***********************************/
void dat_graph_select(gint x, gint y, struct model_pak *model){
        struct graph_pak *graph;
        struct canvas_pak *canvas;
        gdouble ox,dx, xf;
        gdouble oy,dy, yf;
        gint xx,yy;
        int i,j;
        GSList *list;
        FILE *vf;
	gint x_index;
	gint y_index;
        g_data_x *p_x;
	g_data_y *p_y;

g_assert(model!=NULL);
	graph=(struct graph_pak *)model->graph_active;
g_assert(graph!=NULL);
	canvas = g_slist_nth_data(sysenv.canvas_list, 0);
g_assert(canvas!=NULL);
/*calculate x data*/
	ox=canvas->x + 4*gl_fontsize;
	if (graph->ylabel) ox+=2*gl_fontsize;
	oy=canvas->y + canvas->height - 2*gl_fontsize;
	if (graph->xlabel) oy-=2*gl_fontsize;
	dy=(canvas->height-8.0*gl_fontsize);
	dx=(canvas->width-2.0*ox);

if (graph->title) dy = (canvas->height-10.0*gl_fontsize);
if (graph->x_title) oy = canvas->y + canvas->height - 5*gl_fontsize;

	list=graph->set_list;
	p_x = (g_data_x *) list->data;
	/*get the corresponding x index*/
	x_index=-1;
	for(i=0;i<p_x->x_size;i++){
		xf = p_x->x[i];
		xf -= graph->xmin;
		xf /= (graph->xmax - graph->xmin);
		xx = ox + xf*dx;
		if((x>xx-3)&&(x<xx+3)) {
			/*got it*/
			x_index=i;
			break;
		}
	}
	if(x_index<0) {
		/*not found*/
		if(graph->select_label!=NULL) g_free(graph->select_label);
		graph->select_label=NULL;
		return;
	}
	/*scan for the proper y value*/
	j=0;y_index=-1;
  if((graph->type==GRAPH_IY_TYPE)||(graph->type==GRAPH_XY_TYPE)){
	for ( ; list ; list=g_slist_next(list)){
		p_y = (g_data_y *) list->data;
		yf = p_y->y[x_index];/*only need to look here*/
		yf -= graph->ymin;
		yf /= (graph->ymax - graph->ymin);
		yf *= dy;
		yy = (gint) yf;
		yy *= -1;
		yy += oy;
		if((y>yy-3)&&(y<yy+3)){
			/*got it*/
			y_index=1;
			graph->select=x_index;
			graph->select_2=p_y->y[x_index];
			if(graph->select_label) g_free(graph->select_label);
			if(graph->type==GRAPH_IY_TYPE) graph->select_label=g_strdup_printf("[%i,%f]",(gint)p_x->x[x_index],p_y->y[x_index]);
			else graph->select_label=g_strdup_printf("[%G,%G]",p_x->x[x_index],p_y->y[x_index]);
			if(p_y->idx!=NULL) {
				if(p_y->idx[x_index]<0) {
					update_frame_uspex(p_y->idx[x_index]-1,model);
					return;/*unavailable*/
				}
				/*load structure*/
				vf=fopen(model->filename, "r");
				if(!vf) return;
				model->cur_frame=p_y->idx[x_index]-1;
				read_raw_frame(vf,p_y->idx[x_index]-1,model);
				fclose(vf);
				tree_model_refresh(model);
				model_prep(model);
			}
			break;
		}
		j++;
	}
	if(y_index<0) {
		/*not found*/
		if(graph->select_label!=NULL) g_free(graph->select_label);
		graph->select_label=NULL;
		return;
	}
  }else if((graph->type==GRAPH_IX_TYPE)||(graph->type==GRAPH_XX_TYPE)){
	/*in this case, we know the data is on the x_index set*/
	list=g_slist_nth(graph->set_list,x_index+1);
	if(!list) return;/*missing data?*/
	p_y = (g_data_y *) list->data;
	for(j=0;j<p_y->y_size;j++){
		yf = p_y->y[j];
		yf -= graph->ymin;
		yf /= (graph->ymax - graph->ymin);
		yf *= dy;
		yy = (gint) yf;
		yy *= -1;
		yy += oy;
		if((y>yy-3)&&(y<yy+3)){
			/*got it*/
			y_index=j;
			graph->select=x_index;
			graph->select_2=p_y->y[y_index];
			if(graph->select_label) g_free(graph->select_label);
			if(graph->type==GRAPH_IX_TYPE) graph->select_label=g_strdup_printf("[%i,%f]",(gint)p_x->x[x_index],p_y->y[y_index]);
			else graph->select_label=g_strdup_printf("[%G,%G]",p_x->x[x_index],p_y->y[y_index]);
			if(p_y->idx!=NULL) {
				if(p_y->idx[y_index]<0) {
					update_frame_uspex(p_y->idx[y_index],model);
					return;/*unavailable*/
				}
				/*load structure*/
				vf=fopen(model->filename, "r");
				if(!vf) return;
				model->cur_frame=p_y->idx[y_index];
				read_raw_frame(vf,p_y->idx[y_index],model);
				fclose(vf);
				tree_model_refresh(model);
				model_prep(model);
			}
			break;
		}
	}
	if(y_index<0) {
		/*not found*/
		if(graph->select_label!=NULL) g_free(graph->select_label);
		graph->select_label=NULL;
		return;
	}
  }else {
	/*unknown graph type*/
	if(graph->select_label!=NULL) g_free(graph->select_label);
	graph->select_label=NULL;
	return;
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
/*NEW - dat_graph*/
g_data_x *gx=NULL;
g_data_y *gy=NULL;
graph_type type=graph->type;
graph_line line=GRAPH_LINE_NONE;
/* compute origin */
ox = canvas->x + 4*gl_fontsize;
if (graph->ylabel) ox += 2*gl_fontsize;
oy = canvas->y + canvas->height - 2*gl_fontsize;
if (graph->xlabel) oy -= 2*gl_fontsize;
/* increments for screen drawing */
dy = (canvas->height-8.0*gl_fontsize);
dx = (canvas->width-2.0*ox);
if (graph->title) {
//  oy += 2*gl_fontsize;
  dy = (canvas->height-10.0*gl_fontsize);
}
if (graph->x_title) {
/*we have a x_title we need to offset y*/
  oy = canvas->y + canvas->height - 5*gl_fontsize;
}


/* axes label colour */
glColor3f(sysenv.render.fg_colour[0],sysenv.render.fg_colour[1],sysenv.render.fg_colour[2]);
glLineWidth(2.0);
/* x labels */
if((type==GRAPH_DOS)&&(graph->xticks>2)){
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
	    pango_print(text, x-2*gl_fontsize, oy+gl_fontsize, canvas, gl_fontsize-2,0);
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
  if(type==GRAPH_BANDOS){
    x = ox + xf*dx*0.3;
  }else{
    x = ox + xf*dx;
  }
  if (graph->xlabel){
	/*only calculate real value when needed*/
	if(type!=GRAPH_BANDOS){
	    xf *= (graph->xmax-graph->xmin);
	    xf += graph->xmin;
	}
if((type==GRAPH_IY_TYPE)||(type==GRAPH_IX_TYPE)) {
  text = g_strdup_printf("%i",(gint)xf);
  pango_print(text, x, oy+gl_fontsize, canvas, gl_fontsize-2,0);
} else {
  text = g_strdup_printf("%.2f", xf);
  pango_print(text, x-2*gl_fontsize, oy+gl_fontsize, canvas, gl_fontsize-2,0);
}
    g_free(text);
  }
  /* axis segment + tick */
  glBegin(GL_LINE_STRIP);
  gl_vertex_window(oldx, oy, canvas);
  gl_vertex_window(x, oy, canvas);
  gl_vertex_window(x, oy+5, canvas);
  glEnd();
  if((type==GRAPH_IY_TYPE)||(type==GRAPH_XY_TYPE)||(type==GRAPH_IX_TYPE)||(type==GRAPH_XX_TYPE)){
	  glBegin(GL_LINE_STRIP);
	  gl_vertex_window(oldx,oy-dy-1,canvas);
	  gl_vertex_window(x,oy-dy-1,canvas);
	  gl_vertex_window(x,oy-dy-5,canvas);
	  glEnd();
  }
  oldx = x;
  }
}
/* another x axis */
if(type==GRAPH_BANDOS){
  /*second x axis*/
  oldx = ox+dx*0.4;
  for (i=0 ; i<graph->xticks ; i++){
  /* get real index */
  xf = (gdouble) i / (gdouble) (graph->xticks-1);
  x = ox + dx*0.4 + xf*dx*0.6;
  if (graph->xlabel){
    /*only calculate real value when needed*/
    text = g_strdup_printf("%.2f", xf);
    pango_print(text, x-2*gl_fontsize, oy+gl_fontsize, canvas, gl_fontsize-2,0);
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
if(((type==GRAPH_BANDOS)||(type==GRAPH_BAND))&&(graph->yticks>2)){
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
	pango_print(text, 2*gl_fontsize, y-1, canvas,gl_fontsize-2,0);
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
	pango_print(text, 2*gl_fontsize, y-1, canvas,gl_fontsize-2,0);
	g_free(text);
    }
    /* axis segment + tick */
    glBegin(GL_LINE_STRIP);
    gl_vertex_window(ox, oldy, canvas);
    gl_vertex_window(ox, y-1, canvas);
    gl_vertex_window(ox-5, y-1, canvas);
    glEnd();
    if((type==GRAPH_IY_TYPE)||(type==GRAPH_XY_TYPE)||(type==GRAPH_IX_TYPE)||(type==GRAPH_XX_TYPE)){
	    glBegin(GL_LINE_STRIP);
		gl_vertex_window(ox+dx,oldy,canvas);
		gl_vertex_window(ox+dx,y-1,canvas);
		gl_vertex_window(ox+dx+4,y-1,canvas);
		glEnd();
    }
    oldy = y;
  }
}
/* another y axis */
if(type==GRAPH_BANDOS){
  x=ox+0.4*dx;
  y=oy-dy;
  /*FIX: 186c9c*/
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
if(type==GRAPH_BANDOS){
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
if(type==GRAPH_BAND){
/*the first set of GRAPH_BAND data contain X values*/
/*we do not read the full header (yet)*/
	ptr = (gdouble *) list->data;
	xval=&(ptr[4]);
	list=g_slist_next(list);
}
if((type==GRAPH_IY_TYPE)
    ||(type==GRAPH_XY_TYPE)
    ||(type==GRAPH_IX_TYPE)
    ||(type==GRAPH_XX_TYPE)){
/*we have dedicated x,y*/
	gx=(g_data_x *)list->data;
	xval=&(gx->x[0]);
	list=g_slist_next(list);
}
j=0;
for ( ; list ; list=g_slist_next(list))
  {
	/*FIXME here with ALL graph*/
  if((type==GRAPH_IY_TYPE)
      ||(type==GRAPH_XY_TYPE)
      ||(type==GRAPH_IX_TYPE)
      ||(type==GRAPH_XX_TYPE)){
	gy=(g_data_y *)list->data;
	ptr=&(gy->y[0]);
	size=gy->y_size;
	line=gy->line;
	if((ptr==NULL)||(size==0)){
		/*empty data set is OK*/
		j++;
		continue;
	}
	if(gy->type!=GRAPH_REGULAR) type=gy->type;/*update type within type?*/
	switch(gy->color){
		case GRAPH_COLOR_WHITE:/*1.0,  1.0,    1.0*/
			glColor3f(1.0,1.0,1.0);break;
		case GRAPH_COLOR_BLUE:/*0.0,  0.0,    1.0*/
			glColor3f(0.0,0.0,1.0);break;
		case GRAPH_COLOR_GREEN:/*0.0,  0.5,    0.0*/
			glColor3f(0.0,0.5,0.0);break;
		case GRAPH_COLOR_RED:/*1.0,  0.0,    0.0*/
			glColor3f(1.0,0.0,0.0);break;
		case GRAPH_COLOR_YELLOW:/*1.0,  1.0,    0.0*/
			glColor3f(1.0,1.0,0.0);break;
		case GRAPH_COLOR_GRAY:/*0.5,  0.5,    0.5*/
			glColor3f(0.5,0.5,0.5);break;
		case GRAPH_COLOR_NAVY:/*0.0,  0.0,    0.5*/
			glColor3f(0.0,0.0,0.5);break;
		case GRAPH_COLOR_LIME:/*0.0,  1.0,    0.0*/
			glColor3f(0.0,1.0,0.0);break;
		case GRAPH_COLOR_TEAL:/*0.0,  0.5,    0.5*/
			glColor3f(0.0,0.5,0.5);break;
		case GRAPH_COLOR_AQUA:/*0.0,  1.0,    1.0*/
			glColor3f(0.0,1.0,1.0);break;
		case GRAPH_COLOR_MAROON:/*0.5,  0.0,    0.0*/
			glColor3f(0.5,0.0,0.0);break;
		case GRAPH_COLOR_PURPLE:/*0.5,  0.0,    0.5*/
			glColor3f(0.5,0.0,0.5);break;
		case GRAPH_COLOR_OLIVE:/*0.5,  0.5,    0.0*/
			glColor3f(0.5,0.5,0.0);break;
		case GRAPH_COLOR_SILVER:/*0.75, 0.75,   0.75*/
			glColor3f(0.75,0.75,0.75);break;
		case GRAPH_COLOR_FUSHIA:/*1.0,  0.0,    1.0*/
			glColor3f(1.0,0.0,1.0);break;
		case GRAPH_COLOR_BLACK:/*0.0,0.0,0.0*/
			glColor3f(0.,0.,0.);
		case GRAPH_COLOR_DEFAULT:/*whatever title color*/
		default:
			glColor3f(sysenv.render.title_colour[0],sysenv.render.title_colour[1],sysenv.render.title_colour[2]);
	}
  } else {
	ptr = (gdouble *) list->data;
  }
  oldx=-1;
  oldy=-1;
  i=0;
  while(i<size)
    {
/*exclude non-value*/
if(isnan(ptr[i])) {
	i++;
	oldx=-1;
	continue;
}
/*get real values*/
switch (type){
	case GRAPH_IY_TYPE:
	case GRAPH_XY_TYPE:
	xf = xval[i];
	yf = ptr[i];
	break;
	case GRAPH_IX_TYPE:
	case GRAPH_XX_TYPE:
	xf = xval[j];
	yf = ptr[i];
	break;
	case GRAPH_FREQUENCY:
	case GRAPH_DOS:
	xf = ptr[i];
	yf = ptr[i+1];
	break;
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
	switch(type){
	  case GRAPH_FREQUENCY:
	  case GRAPH_DOS:
		i++;/* ie twice ++ */
	  case GRAPH_IY_TYPE:
	  case GRAPH_XY_TYPE:
	  case GRAPH_IX_TYPE:
	  case GRAPH_XX_TYPE:
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
switch (type){
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
	case GRAPH_IY_TYPE:
	case GRAPH_XY_TYPE:
	case GRAPH_IX_TYPE:
	case GRAPH_XX_TYPE:
//	xf += 1.;
	xf -= graph->xmin;
	xf /= (graph->xmax - graph->xmin);
	x = ox + xf*dx;
	yf -= graph->ymin;
	yf /= (graph->ymax - graph->ymin);
	yf *= dy;
	y = (gint) yf;
	y *= -1;
	y += oy;
if(gy->symbol){
	switch(gy->symbol[i]){
	case GRAPH_SYMB_SQUARE:
		glBegin(GL_LINE_STRIP);
		gl_vertex_window(x-2, y-2, canvas);
		gl_vertex_window(x+2, y-2, canvas);
		gl_vertex_window(x+2, y+2, canvas);
		gl_vertex_window(x-2, y+2, canvas);
		gl_vertex_window(x-2, y-2, canvas);
		glEnd();
	break;
	case GRAPH_SYMB_CROSS:
		glBegin(GL_LINE_STRIP);
		gl_vertex_window(x-2, y-2, canvas);
		gl_vertex_window(x+2, y+2, canvas);
		glEnd();
		glBegin(GL_LINE_STRIP);
		gl_vertex_window(x+2, y-2, canvas);
		gl_vertex_window(x-2, y+2, canvas);
		glEnd();
	break;
	case GRAPH_SYMB_TRI_DN:
		glBegin(GL_LINE_STRIP);
		gl_vertex_window(x, y-3, canvas);
		gl_vertex_window(x-3, y+3, canvas);
		gl_vertex_window(x+3, y+3, canvas);
		gl_vertex_window(x, y-3, canvas);
		glEnd();
	break;
	case GRAPH_SYMB_TRI_UP:
		glBegin(GL_LINE_STRIP);
		gl_vertex_window(x, y+3, canvas);
		gl_vertex_window(x-3, y-3, canvas);
		gl_vertex_window(x+3, y-3, canvas);
		gl_vertex_window(x, y+3, canvas);
		glEnd();
	break;
	case GRAPH_SYMB_DIAM:
		glBegin(GL_LINE_STRIP);
		gl_vertex_window(x, y+5, canvas);
		gl_vertex_window(x-5, y, canvas);
		gl_vertex_window(x, y-5, canvas);
		gl_vertex_window(x+5, y, canvas);
		gl_vertex_window(x, y+5, canvas);
		glEnd();
	break;
	case GRAPH_SYMB_NONE:
	default:
		/*no symbol*/
	break;
	}
}else{
	/*draw a default rectangle*/
	glBegin(GL_LINE_STRIP);
	gl_vertex_window(x-2, y-2, canvas);
	gl_vertex_window(x+2, y-2, canvas);
	gl_vertex_window(x+2, y+2, canvas);
	gl_vertex_window(x-2, y+2, canvas);
	gl_vertex_window(x-2, y-2, canvas);
	glEnd();
}
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
  switch(type){
	case GRAPH_IY_TYPE:
	case GRAPH_XY_TYPE:
	/* NEW: use line*/
		switch(line){
		case GRAPH_LINE_THICK:
			glLineWidth(2.0);
		case GRAPH_LINE_SINGLE:
			glBegin(GL_LINE_STRIP);
			gl_vertex_window(oldx, oldy-1, canvas);
			gl_vertex_window(x, y-1, canvas);
			glEnd();
			glLineWidth(1.0);
		break;
		case GRAPH_LINE_DASH:
			glEnable(GL_LINE_STIPPLE);
			glLineStipple(1,0xCCCC);
			glBegin(GL_LINE_STRIP);
			gl_vertex_window(oldx, oldy-1, canvas);
			gl_vertex_window(x, y-1, canvas);
			glEnd();
			glDisable(GL_LINE_STIPPLE);
		case GRAPH_LINE_DOT:
			glEnable(GL_LINE_STIPPLE);
			glLineStipple(1,0xAAAA);
			glBegin(GL_LINE_STRIP);
			gl_vertex_window(oldx, oldy-1, canvas);
			gl_vertex_window(x, y-1, canvas);
			glEnd();
			glDisable(GL_LINE_STIPPLE);
		break;
		case GRAPH_LINE_NONE:
		default:
			break;
		}
	break;
	case GRAPH_IX_TYPE:
	case GRAPH_XX_TYPE:
		/*not ready yet*/
	break;
	case GRAPH_FREQUENCY:
		/*we should never reach here*/
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
	switch(type){
	case GRAPH_FREQUENCY:
	case GRAPH_DOS:
		i++;/* ie twice ++ */
	case GRAPH_IY_TYPE:
	case GRAPH_XY_TYPE:
	case GRAPH_IX_TYPE:
	case GRAPH_XX_TYPE:
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
switch (type){
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
	case GRAPH_IY_TYPE:
	case GRAPH_XY_TYPE:
	case GRAPH_IX_TYPE:
	case GRAPH_XX_TYPE:
		flag=(graph->select_label!=NULL);
		break;
	case GRAPH_FREQUENCY:
	case GRAPH_REGULAR:
	default:
	break;
}
/* NEW - selector*/
	
	if((flag)||(graph->select_label!=NULL)){
		gint xoff;
		list=graph->set_list;
		glEnable(GL_LINE_STIPPLE);
		glLineStipple(1, 0x0303);
		glColor3f(0.9, 0.7, 0.4);
		glLineWidth(2.0);
		/*depend on graph->type ; TODO: manage changing types*/
		switch(graph->type){
		case GRAPH_IX_TYPE:
		case GRAPH_XX_TYPE:
			gx=(g_data_x *)list->data;
			xf=gx->x[graph->select];
			yf=graph->select_2;
			xf -= graph->xmin;
			xf /= (graph->xmax - graph->xmin);
			sx = ox + xf*dx;
			yf -= graph->ymin;
			yf /= (graph->ymax - graph->ymin);
			yf *= dy;
			sy = (gint) yf;
			sy *= -1;
			sy += oy;
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
		break;
		case GRAPH_IY_TYPE:
		case GRAPH_XY_TYPE:
			gx=(g_data_x *)list->data;
			xf=gx->x[graph->select];
			yf=graph->select_2;
			xf -= graph->xmin;
			xf /= (graph->xmax - graph->xmin);
			sx = ox + xf*dx;
			yf -= graph->ymin;
			yf /= (graph->ymax - graph->ymin);
			yf *= dy;
			sy = (gint) yf;
			sy *= -1;
			sy += oy;
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
		break;
		default:
			glBegin(GL_LINES);
			gl_vertex_window(sx, sy-10, canvas);
			gl_vertex_window(sx, 3*gl_fontsize, canvas);
			glEnd();
		}
	/*label is always same (I think)*/
	xoff = (gint) strlen(graph->select_label);
	xoff *= gl_fontsize;
	xoff /= 4;
	pango_print(graph->select_label, sx-xoff, 0, canvas,gl_fontsize-2,0);
	glDisable(GL_LINE_STIPPLE);
	}
/* NEW: set titles */
if(graph->title){
  pango_print(graph->title,ox+2,oy-(gint)dy-4*gl_fontsize,canvas,gl_fontsize,0);
}
if(graph->sub_title){
  pango_print(graph->sub_title,ox+0.5*dx+2, oy-dy-2*gl_fontsize, canvas,gl_fontsize,0);
}
if(graph->x_title){
  pango_print(graph->x_title,ox+0.33*dx+2,oy+2*gl_fontsize,canvas,gl_fontsize,0);
}
if(graph->y_title){
  pango_print(graph->y_title,0,oy-0.33*dy,canvas,gl_fontsize,90);
}
/*axis?*/
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
	case GRAPH_IY_TYPE:
	case GRAPH_XY_TYPE:
	case GRAPH_IX_TYPE:
	case GRAPH_XX_TYPE:
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
n=0;x[0]=0.;/*FIX 09bf1e*/
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
