/*
Copyright (C) 2018 by Okadome Valencia

hubert.valencia _at_ imass.nagoya-u.ac.jp

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

/* plot simple data (energy, force, pressure, volume) */
/* and complex data (bands, dos, bandos, frequency) */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "gdis.h"
#include "coords.h"
#include "model.h"
#include "file.h"
#include "graph.h"
#include "interface.h"
#include "plots.h"

extern struct elem_pak elements[];
extern struct sysenv_pak sysenv;


void plot_prepare_data(struct plot_pak *plot){
	int i,j,idx;
	gdouble x,y;
	struct model_pak *model;
	g_assert(plot != NULL);
	g_assert(plot->model != NULL);
	model=plot->model;
/* task cleanup set plots that does not depend on frames */
/* add BAND data */
	if(plot->plot_mask&PLOT_BAND){
	/*we have an array of [kpts_d,band_up,band_down] which size is [nkpoints*nbands]*/
		/*pass all parameters as the begining of an array*/
		plot->band.data[0]=model->nbands;
		plot->band.data[1]=model->nkpoints;
		plot->band.data[2]=model->efermi;/*necessary?*/
		if(model->spin_polarized) plot->band.data[3]=2;
		else plot->band.data[3]=1;
		plot->band.xmin=model->kpts_d[0];
		plot->band.xmax=model->kpts_d[0];
		plot->band.ymin=model->band_up[0]-model->efermi;
		plot->band.ymax=model->band_up[0]-model->efermi;
		/*first prepare the x*/
		for(i=0;i<(model->nkpoints);++i){
			x=model->kpts_d[i];/*kpoint distance*/
			if(x<plot->band.xmin) plot->band.xmin=x;
			if(x>plot->band.xmax) plot->band.xmax=x;
			plot->band.data[i+4]=x;/*mind the shift*/
		}
		/*very tricky: reorder all kpoints energy per band*/
		idx=4+model->nkpoints;/*origin shift*/
		for(i=0;i<model->nbands;i++)
			for(j=0;j<model->nkpoints;j++){
				y=model->band_up[i+j*model->nbands]-model->efermi;
				if(y<plot->band.ymin) plot->band.ymin=y;
				if(y>plot->band.ymax) plot->band.ymax=y;
				plot->band.data[idx]=y;/*mind the shift*/
//fprintf(stdout,"#DBG: REG: %lf %lf (%i %i %i)\n",plot->band.data[j+4],y,j,i,idx);
				idx++;
			}
		/*all done*/
	}
/* add DOS data */
        if(plot->plot_mask&PLOT_DOS){
                j=0;
		/*initial values*/
                plot->dos.xmin=model->dos_eval[0]-model->efermi;
		plot->dos.xmax=model->dos_eval[model->ndos-1]-model->efermi;
		plot->dos.ymin=0.;/*for now: all dos positive*/
		plot->dos.ymax=model->dos_spin_up[0];
		/*all values*/
                for(i=0;i<model->ndos;++i){
                        x=model->dos_eval[i];
			x-=model->efermi;/*scale to fermi level!*/
			if(x<plot->dos.xmin) plot->dos.xmin=x;
			if(x>plot->dos.xmax) plot->dos.xmax=x;
			y=model->dos_spin_up[i];
			if(model->spin_polarized) y+=model->dos_spin_down[i];
			if(y>plot->dos.ymax) plot->dos.ymax=y;
			if(y<plot->dos.ymin) y=0.;/*no negative value*/
			plot->dos.data[j]=x;
			plot->dos.data[j+1]=y;
                        j+=2;
                }

        }
	/* TODO: FREQ: add intensities */
	/* TODO: how to signal task is completed?*/
	plot->task->progress=100.0;/*job done*/
	plot->task->status=COMPLETED;
        plot->data_changed=FALSE;/*protect data until changed*/
}
void plot_load_data(struct plot_pak *plot,struct task_pak *task){
	struct model_pak *model;
        FILE *fp;
        int i;
	gint status;
	GString *err_text;
	gchar *key;
	gdouble y=0.0;
        gint nframes;
	gint cur_frame;
/* This part will populate everything that depends on frame reading (ie dynamics) */
if(plot->data_changed==FALSE) return;/*we do not need to reload data*/
	/* TODO: prevent update while reading frame? */
	g_assert(plot != NULL);
	g_assert(plot->model != NULL);
/* simply use model provided in plot */
	model=plot->model;
	plot->task=task;
	if((model->num_frames<2)||(plot->plot_mask==0)) {
		task->progress=50.0;/*job half done (that was fast)*/
		return;/*skip this part*/
	}
	/*open file*/
        fp=fopen(model->filename, "r");
        if (!fp){
                gui_text_show(ERROR,g_strdup_printf("I/O ERROR: can't open!\n"));
                return;
        }
	nframes=model->num_frames;
	cur_frame=model->cur_frame;
	status=read_raw_frame(fp,0,model);
	if(status) {
		err_text = g_string_new("");
		g_string_printf(err_text, "Error reading frame: %d\n", 0);
		gui_text_show(ERROR, err_text->str);
		g_string_free(err_text, TRUE);
		return;
	}
/* populate first data */
	if(plot->plot_mask&PLOT_ENERGY){
		/*prepare preload data for energy*/
		key=property_lookup("Energy",model);
		if(key) {/*initialize range*/
			sscanf(key,"%lf %*s",&y);
			plot->energy.ymin=y;
			plot->energy.ymax=y;
			g_free(key);
			if(plot->energy.data!=NULL) {
				g_free(plot->energy.data);
				plot->energy.data=NULL;
			}
			plot->energy.data=(gdouble *)g_malloc(nframes*sizeof(gdouble));
			plot->energy.size=nframes;
			plot->energy.xmin=0;
			plot->energy.xmax=nframes-1;
			plot->energy.data[0]=y;
		}else{/*invalidate energy*/
			plot->plot_mask^=PLOT_ENERGY;
		}
	}
	if(plot->plot_mask&PLOT_FORCE){
		/*prepare preload data for forces*/
		key=property_lookup("Force",model);
		if(key) {/*initialize range*/
			sscanf(key,"%lf %*s",&y);
			plot->force.ymin=y;
			plot->force.ymax=y;
			g_free(key);
			if(plot->force.data!=NULL) {
				g_free(plot->force.data);
				plot->force.data=NULL;
			}
			plot->force.data=(gdouble *)g_malloc(nframes*sizeof(gdouble));
			plot->force.size=nframes;
			plot->force.xmin=0;
			plot->force.xmax=nframes-1;
			plot->force.data[0]=y;
		}else{/*invalidate force*/
			plot->plot_mask^=PLOT_FORCE;
		}
	}
        if(plot->plot_mask&PLOT_VOLUME){
                /*prepare preload data for volume*/
		/* volume data always exists (but can be just 0)*/
		y=model->volume;
                if(plot->volume.data!=NULL) {
                        g_free(plot->volume.data);
                        plot->volume.data=NULL;
                }
                plot->volume.data=(gdouble *)g_malloc(nframes*sizeof(gdouble));
                plot->volume.size=nframes;
		plot->volume.xmin=0;
		plot->volume.xmax=nframes-1;
		plot->volume.ymin = y;
		plot->volume.ymax = y;
		plot->volume.data[0]=y;
        }
        if(plot->plot_mask&PLOT_PRESSURE){
                /*prepare preload data for pressure*/
		key=property_lookup("Pressure",model);
		if(key) {/*initialize range*/
			sscanf(key,"%lf %*s",&y);
			plot->pressure.ymin=y;
			plot->pressure.ymax=y;
			g_free(key);
			if(plot->pressure.data!=NULL) {
				g_free(plot->pressure.data);
				plot->pressure.data=NULL;
			}
			plot->pressure.data=(gdouble *)g_malloc(nframes*sizeof(gdouble));
			plot->pressure.size=nframes;
			plot->pressure.xmin=0;
			plot->pressure.xmax=nframes-1;
			plot->pressure.data[0]=y;
		}else{/*Pressure is also sometimes called stress*/
			g_free(key);
			key=property_lookup("Stress",model);
			if(key) {/*initialize range*/
				sscanf(key,"%lf %*s",&y);
				plot->pressure.ymin=y;
				plot->pressure.ymax=y;
				g_free(key);
				if(plot->pressure.data!=NULL) {
					g_free(plot->pressure.data);
					plot->pressure.data=NULL;
				}
				plot->pressure.data=(gdouble *)g_malloc(nframes*sizeof(gdouble));
				plot->pressure.size=nframes;
				plot->pressure.xmin=0;
				plot->pressure.xmax=nframes-1;
				plot->pressure.data[0]=y;
			}else{/*invalidate pressure*/
				plot->plot_mask^=PLOT_PRESSURE;
			}
		}
        }
/* preload necessary data */
        for (i=1;i<nframes;i++){
                status=read_raw_frame(fp,i,model);/*this will fail.. always*/
		if(status) {
			err_text = g_string_new("");
			g_string_printf(err_text, "Error reading frame: %d\n", i);
			gui_text_show(ERROR, err_text->str);
			g_string_free(err_text, TRUE);
			return;
		}
		if(plot->plot_mask&PLOT_ENERGY){/*populate energy*/
			key=property_lookup("Energy",model);
			if(key) {
				sscanf(key,"%lf %*s",&y);
				if(y<plot->energy.ymin) plot->energy.ymin=y;
				if(y>plot->energy.ymax) plot->energy.ymax=y;
				plot->energy.data[i]=y;
				g_free(key);
			}
		}
                if(plot->plot_mask&PLOT_FORCE){/*populate force*/
			key=property_lookup("Force",model);
			if(key) {
				sscanf(key,"%lf %*s",&y);
				if(y<plot->force.ymin) plot->force.ymin=y;
				if(y>plot->force.ymax) plot->force.ymax=y;
				plot->force.data[i]=y;
				g_free(key);
			}
                }
                if(plot->plot_mask&PLOT_VOLUME){/*populate volume*/
			y=model->volume;/*easy*/
                        if(y<plot->volume.ymin) plot->volume.ymin=y;
                        if(y>plot->volume.ymax) plot->volume.ymax=y;
                        plot->volume.data[i]=y;
                }
                if(plot->plot_mask&PLOT_PRESSURE){/*populate pressure*/
			key=property_lookup("Pressure",model);
			if(key) {
				sscanf(key,"%lf %*s",&y);
				if(y<plot->pressure.ymin) plot->pressure.ymin=y;
				if(y>plot->pressure.ymax) plot->pressure.ymax=y;
				plot->pressure.data[i]=y;
				g_free(key);
			}else{ /*pressure is also refered to as stress*/
				key=property_lookup("Stress",model);
				if(key) {
					sscanf(key,"%lf %*s",&y);
					if(y<plot->pressure.ymin) plot->pressure.ymin=y;
					if(y>plot->pressure.ymax) plot->pressure.ymax=y;
					plot->pressure.data[i]=y;
					g_free(key);
				}
			}
                }
        }
/* ending */
	model->cur_frame=cur_frame;/*go back to current frame <- useful?*/
	read_raw_frame(fp,cur_frame,model);/*here?*/
        fclose(fp);
	/* all done */
	task->progress=50.0;/*job half done*/
}
/***************/
/* energy plot */
/***************/
void draw_plot_energy(struct model_pak *model,struct plot_pak *plot){
g_data_x gx;
g_data_y gy;
gint idx;
gdouble d;
gdouble xmin,xmax;
gdouble ymin,ymax;
/*avoid rare case _BUG_*/
if(plot->energy.size<1) return;
/*always auto*/
xmin=plot->energy.xmin;
xmax=plot->energy.xmax;
ymin=plot->energy.ymin;
ymax=plot->energy.ymax;
if(ymin==ymax) {
	ymin=ymin-1.0;
	ymax=ymax+1.0;
}
/*add 5% to y limits for "readability"*/
d=ymax-ymin;
ymin=ymin-d*0.05;
ymax=ymax+d*0.05;
plot->graph=graph_new("ENERGY", model);
dat_graph_set_title("<big>Energy <i>vs.</i> Ionic step</big>",plot->graph);
//dat_graph_set_sub_title("<small>(From <b>GDIS</b> data)</small>",plot->graph);
dat_graph_set_x_title("Ionic steps",plot->graph);
dat_graph_set_y_title("Energy (eV)",plot->graph);/*FIXME: add proper unit text!*/
dat_graph_set_type(GRAPH_IY_TYPE,plot->graph);
/*set x*/
gx.x_size=plot->energy.size+1;
gx.x=g_malloc(gx.x_size*sizeof(gdouble));
for(idx=0;idx<gx.x_size;idx++) gx.x[idx]=(gdouble)(idx);
dat_graph_set_x(gx,plot->graph);
g_free(gx.x);
/*set y*/
gy.y_size=plot->energy.size+1;
gy.y=g_malloc(gy.y_size*sizeof(gdouble));
gy.idx=g_malloc(gy.y_size*sizeof(gint32));/*<- will set structure frame automagically*/
gy.symbol=g_malloc(gy.y_size*sizeof(graph_symbol));
gy.mixed_symbol=FALSE;
gy.sym_color=NULL;
gy.y[0]=NAN;/*not a value*/
for(idx=1;idx<gy.y_size;idx++){
	gy.y[idx]=plot->energy.data[idx-1];
	gy.idx[idx]=idx;
	gy.symbol[idx]=GRAPH_SYMB_DIAM;
}
gy.type=GRAPH_IY_TYPE;
dat_graph_set_limits(xmin,xmax,ymin,ymax,plot->graph);
gy.line=GRAPH_LINE_THICK;
gy.color=GRAPH_COLOR_DEFAULT;
dat_graph_add_y(gy,plot->graph);
g_free(gy.y);
g_free(gy.idx);
g_free(gy.symbol);
/*force TRUE on label*/
if(plot->xtics <= 1) graph_set_xticks(TRUE,2,plot->graph);
else graph_set_xticks(TRUE,plot->xtics,plot->graph);
if(plot->ytics <= 1) graph_set_yticks(TRUE,2,plot->graph);
else graph_set_yticks(TRUE,plot->ytics,plot->graph);
}
/**************/
/* force plot */
/**************/
void draw_plot_force(struct model_pak *model,struct plot_pak *plot){
g_data_x gx;
g_data_y gy;
gint idx;
gdouble d;
gdouble xmin,xmax;
gdouble ymin,ymax;
/*avoid rare case _BUG_*/
if(plot->force.size<1) return;
/*always auto*/
xmin=plot->force.xmin;
xmax=plot->force.xmax;
ymin=plot->force.ymin;
ymax=plot->force.ymax;
if(ymin==ymax) {
	ymin=ymin-1.0;
	ymax=ymax+1.0;
}
/*add 5% to y limits for "readability"*/
d=ymax-ymin;
ymin=ymin-d*0.05;
ymax=ymax+d*0.05;
plot->graph=graph_new("FORCE", model);
dat_graph_set_title("<big>Force <i>vs.</i> Ionic step</big>",plot->graph);
//dat_graph_set_sub_title("<small>(From <b>GDIS</b> data)</small>",plot->graph);
dat_graph_set_x_title("Ionic steps",plot->graph);
dat_graph_set_y_title("Force (eV/Ang)",plot->graph);/*FIXME: add proper unit text!*/
dat_graph_set_type(GRAPH_IY_TYPE,plot->graph);
/*set x*/
gx.x_size=plot->force.size+1;
gx.x=g_malloc(gx.x_size*sizeof(gdouble));
for(idx=0;idx<gx.x_size;idx++) gx.x[idx]=(gdouble)(idx);
dat_graph_set_x(gx,plot->graph);
g_free(gx.x);
/*set y*/
gy.y_size=plot->force.size+1;
gy.y=g_malloc(gy.y_size*sizeof(gdouble));
gy.idx=g_malloc(gy.y_size*sizeof(gint32));/*<- will set structure frame automagically*/
gy.symbol=g_malloc(gy.y_size*sizeof(graph_symbol));
gy.mixed_symbol=FALSE;
gy.sym_color=NULL;
gy.y[0]=NAN;/*not a value*/
for(idx=1;idx<gy.y_size;idx++){
        gy.y[idx]=plot->force.data[idx-1];
        gy.idx[idx]=idx;
        gy.symbol[idx]=GRAPH_SYMB_DIAM;
}
gy.type=GRAPH_IY_TYPE;
dat_graph_set_limits(xmin,xmax,ymin,ymax,plot->graph);
gy.line=GRAPH_LINE_THICK;
gy.color=GRAPH_COLOR_DEFAULT;
dat_graph_add_y(gy,plot->graph);
g_free(gy.y);
g_free(gy.idx);
g_free(gy.symbol);
/*force TRUE on label*/
if(plot->xtics <= 1) graph_set_xticks(TRUE,2,plot->graph);
else graph_set_xticks(TRUE,plot->xtics,plot->graph);
if(plot->ytics <= 1) graph_set_yticks(TRUE,2,plot->graph);
else graph_set_yticks(TRUE,plot->ytics,plot->graph);
}
/***************/
/* plot volume */
/***************/
void draw_plot_volume(struct model_pak *model,struct plot_pak *plot){
g_data_x gx;
g_data_y gy;
gint idx;
gdouble d;
gdouble xmin,xmax;
gdouble ymin,ymax;
/*avoid rare case _BUG_*/
if(plot->volume.size<1) return;
/*always auto*/
xmin=plot->volume.xmin;
xmax=plot->volume.xmax;
ymin=plot->volume.ymin;
ymax=plot->volume.ymax;
if(ymin==ymax) {
	ymin=ymin-1.0;
	ymax=ymax+1.0;
}
/*add 5% to y limits for "readability"*/
d=ymax-ymin;
ymin=ymin-d*0.05;
ymax=ymax+d*0.05;
plot->graph=graph_new("VOLUME", model);
dat_graph_set_title("<big>Volume <i>vs.</i> Ionic step</big>",plot->graph);
//dat_graph_set_sub_title("<small>(From <b>GDIS</b> data)</small>",plot->graph);
dat_graph_set_x_title("Ionic steps",plot->graph);
dat_graph_set_y_title("Volume (Ang<sup>3</sup>)",plot->graph);/*FIXME: add proper unit text!*/
dat_graph_set_type(GRAPH_IY_TYPE,plot->graph);
/*set x*/
gx.x_size=plot->volume.size+1;
gx.x=g_malloc(gx.x_size*sizeof(gdouble));
for(idx=0;idx<gx.x_size;idx++) gx.x[idx]=(gdouble)(idx);
dat_graph_set_x(gx,plot->graph);
g_free(gx.x);
/*set y*/
gy.y_size=plot->volume.size+1;
gy.y=g_malloc(gy.y_size*sizeof(gdouble));
gy.idx=g_malloc(gy.y_size*sizeof(gint32));/*<- will set structure frame automagically*/
gy.symbol=g_malloc(gy.y_size*sizeof(graph_symbol));
gy.mixed_symbol=FALSE;
gy.sym_color=NULL;
gy.y[0]=NAN;/*not a value*/
for(idx=1;idx<gy.y_size;idx++){
        gy.y[idx]=plot->volume.data[idx-1];
        gy.idx[idx]=idx;
        gy.symbol[idx]=GRAPH_SYMB_DIAM;
}
gy.type=GRAPH_IY_TYPE;
dat_graph_set_limits(xmin,xmax,ymin,ymax,plot->graph);
gy.line=GRAPH_LINE_THICK;
gy.color=GRAPH_COLOR_DEFAULT;
dat_graph_add_y(gy,plot->graph);
g_free(gy.y);
g_free(gy.idx);
g_free(gy.symbol);
/*force TRUE on label*/
if(plot->xtics <= 1) graph_set_xticks(TRUE,2,plot->graph);
else graph_set_xticks(TRUE,plot->xtics,plot->graph);
if(plot->ytics <= 1) graph_set_yticks(TRUE,2,plot->graph);
else graph_set_yticks(TRUE,plot->ytics,plot->graph);
}
/*****************/
/* plot pressure */
/*****************/
void draw_plot_pressure(struct model_pak *model,struct plot_pak *plot){
g_data_x gx;
g_data_y gy;
gint idx;
gdouble d;
gdouble xmin,xmax;
gdouble ymin,ymax;
/*avoid rare case _BUG_*/
if(plot->pressure.size<1) return;
/*always auto*/
xmin=plot->pressure.xmin;
xmax=plot->pressure.xmax;
ymin=plot->pressure.ymin;
ymax=plot->pressure.ymax;
if(ymin==ymax) {
	ymin=ymin-1.0;
	ymax=ymax+1.0;
}
/*add 5% to y limits for "readability"*/
d=ymax-ymin;
ymin=ymin-d*0.05;
ymax=ymax+d*0.05;
plot->graph=graph_new("PRESSURE", model);
dat_graph_set_title("<big>Pressure <i>vs.</i> Ionic step</big>",plot->graph);
//dat_graph_set_sub_title("<small>(From <b>GDIS</b> data)</small>",plot->graph);
dat_graph_set_x_title("Ionic steps",plot->graph);
dat_graph_set_y_title("Pressure (kB)",plot->graph);/*FIXME: add proper unit text!*/
dat_graph_set_type(GRAPH_IY_TYPE,plot->graph);
/*set x*/
gx.x_size=plot->pressure.size+1;
gx.x=g_malloc(gx.x_size*sizeof(gdouble));
for(idx=0;idx<gx.x_size;idx++) gx.x[idx]=(gdouble)(idx);
dat_graph_set_x(gx,plot->graph);
g_free(gx.x);
/*set y*/
gy.y_size=plot->pressure.size+1;
gy.y=g_malloc(gy.y_size*sizeof(gdouble));
gy.idx=g_malloc(gy.y_size*sizeof(gint32));/*<- will set structure frame automagically*/
gy.symbol=g_malloc(gy.y_size*sizeof(graph_symbol));
gy.mixed_symbol=FALSE;
gy.sym_color=NULL;
gy.y[0]=NAN;/*not a value*/
for(idx=1;idx<gy.y_size;idx++){
        gy.y[idx]=plot->pressure.data[idx-1];
        gy.idx[idx]=idx;
        gy.symbol[idx]=GRAPH_SYMB_DIAM;
}
gy.type=GRAPH_IY_TYPE;
dat_graph_set_limits(xmin,xmax,ymin,ymax,plot->graph);
gy.line=GRAPH_LINE_THICK;
gy.color=GRAPH_COLOR_DEFAULT;
dat_graph_add_y(gy,plot->graph);
g_free(gy.y);
g_free(gy.idx);
g_free(gy.symbol);
/*force TRUE on label*/
if(plot->xtics <= 1) graph_set_xticks(TRUE,2,plot->graph);
else graph_set_xticks(TRUE,plot->xtics,plot->graph);
if(plot->ytics <= 1) graph_set_yticks(TRUE,2,plot->graph);
else graph_set_yticks(TRUE,plot->ytics,plot->graph);
}
/************/
/* plot DOS */
/************/
void draw_plot_dos(struct model_pak *model,struct plot_pak *plot){
g_data_x gx;
g_data_y gy;
gint idx;
gdouble xmin,xmax;
gdouble ymin,ymax;
/*avoid rare case _BUG_*/
if(model->ndos<1) return;
/*always auto*/
xmin=plot->dos.xmin;
xmax=plot->dos.xmax;
ymin=plot->dos.ymin;
ymax=plot->dos.ymax;
if(ymin==ymax) {
	ymin=ymin-1.0;
	ymax=ymax+1.0;
}
/*add 5% to ymax limit for "readability"*/
ymax=ymax+(ymax-ymin)*0.05;
plot->graph=graph_new("DOS", model);
dat_graph_toggle_yaxis(plot->graph);
dat_graph_set_title("<big>Density of states</big>",plot->graph);
//dat_graph_set_sub_title("<small>(From <b>GDIS</b> data)</small>",plot->graph);
dat_graph_set_x_title("Energy (eV)",plot->graph);/*FIXME: add proper unit text!*/
dat_graph_set_y_title("Density (state/eV)",plot->graph);/*FIXME: add proper unit text!*/
dat_graph_set_type(GRAPH_XY_TYPE,plot->graph);
/*set x*/
gx.x_size=model->ndos;
gx.x=g_malloc(gx.x_size*sizeof(gdouble));
for(idx=0;idx<gx.x_size;idx++) gx.x[idx]=plot->dos.data[idx*2];
dat_graph_set_x(gx,plot->graph);
g_free(gx.x);
/*set y*/
gy.y_size=model->ndos;
gy.y=g_malloc(gy.y_size*sizeof(gdouble));
gy.idx=g_malloc(gy.y_size*sizeof(gint32));/*<- will NOT set any structure*/
gy.symbol=g_malloc(gy.y_size*sizeof(graph_symbol));
gy.mixed_symbol=FALSE;
gy.sym_color=NULL;
for(idx=0;idx<gy.y_size;idx++){
        gy.y[idx]=plot->dos.data[idx*2+1];
        gy.idx[idx]=-1;
        gy.symbol[idx]=GRAPH_SYMB_NONE;
}
gy.type=GRAPH_XY_TYPE;
dat_graph_set_limits(xmin,xmax,ymin,ymax,plot->graph);
gy.line=GRAPH_LINE_THICK;
gy.color=GRAPH_COLOR_DEFAULT;
dat_graph_add_y(gy,plot->graph);
g_free(gy.y);
g_free(gy.idx);
g_free(gy.symbol);
/*force TRUE on label*/
if(plot->xtics <= 1) graph_set_xticks(TRUE,2,plot->graph);
else graph_set_xticks(TRUE,plot->xtics,plot->graph);
if(plot->ytics <= 1) graph_set_yticks(TRUE,2,plot->graph);
else graph_set_yticks(TRUE,plot->ytics,plot->graph);
}
/***************/
/* plot BANDOS */
/***************/
void draw_plot_bandos(struct model_pak *model,struct plot_pak *plot){
g_data_x gx;
g_data_y gy;
gint idx,jdx;
gdouble scale;
gdouble xmin,xmax;
gdouble ymin,ymax;
#ifdef UNUSED_BUT_SET
/*always auto*/
xmin=plot->dos.ymin;
xmax=plot->dos.ymax;
ymin=plot->band.ymin;
ymax=plot->band.ymax;
if(ymin==ymax) {
	ymin=ymin-1.0;
	ymax=ymax+1.0;
}
#endif //UNUSED_BUT_SET
/*calculate scale*/
scale=plot->band.xmax*0.5/plot->dos.ymax;
plot->graph=graph_new("BANDOS", model);
dat_graph_toggle_xaxis(plot->graph);
dat_graph_toggle_yaxis(plot->graph);
dat_graph_set_title("<big>Density of states</big>",plot->graph);
dat_graph_set_sub_title("<big>Band structure</big>",plot->graph);
dat_graph_set_x_title("density (state/eV)\tk-points distance (Ang)",plot->graph);/*FIXME: add proper unit text!*/
dat_graph_set_y_title("Eigenvalues (eV)",plot->graph);/*FIXME: add proper unit text!*/
dat_graph_set_type(GRAPH_YX_TYPE,plot->graph);/*<-note the inversion!!*/
/*set x*/
gx.x_size=model->nkpoints+1;/*ADD DOS origin for the first set*/
gx.x=g_malloc(gx.x_size*sizeof(gdouble));
gx.x[0]=-1.0*scale*plot->dos.ymax;
for(idx=0;idx<gx.x_size-1;idx++) gx.x[idx+1]=model->kpts_d[idx];/*FIX*/
dat_graph_set_x(gx,plot->graph);
g_free(gx.x);
/*first set is DOS*/
gy.y_size=2*model->ndos;/*interleaved data*/
gy.y=g_malloc((model->ndos*2)*sizeof(gdouble));
gy.idx=g_malloc((model->ndos*2)*sizeof(gint32));/*<- will NOT set any structure*/
gy.symbol=g_malloc((model->ndos*2)*sizeof(graph_symbol));
gy.mixed_symbol=FALSE;
gy.sym_color=NULL;
idx=0;
while(idx<(model->ndos*2)){
	/*idx   <- x: plot tags READ*/
	gy.y[idx]=plot->dos.data[idx];
	gy.idx[idx]=-1;
	gy.symbol[idx]=GRAPH_SYMB_NONE;
	idx++;
	/*idx+1 <- y: plot tags IGNORED*/
	gy.y[idx]=-1.0*scale*plot->dos.data[idx];
	gy.idx[idx]=-1;
	gy.symbol[idx]=GRAPH_SYMB_NONE;
	idx++;
}
gy.type=GRAPH_YX_TYPE;/*<-note the inversion!!*/
gy.line=GRAPH_LINE_THICK;
gy.color=GRAPH_COLOR_DEFAULT;
dat_graph_add_y(gy,plot->graph);
g_free(gy.y);
g_free(gy.idx);
g_free(gy.symbol);
/*next sets are BAND*/
gy.y_size=model->nkpoints+1;/*add DOS origin*/
gy.y=g_malloc(gy.y_size*sizeof(gdouble));
gy.idx=g_malloc(gy.y_size*sizeof(gint32));/*<- will NOT set any structure*/
gy.symbol=g_malloc(gy.y_size*sizeof(graph_symbol));
gy.mixed_symbol=FALSE;
gy.sym_color=NULL;
for(jdx=0;jdx<model->nbands;jdx++){
gy.y[0]=NAN;
for(idx=0;idx<model->nkpoints;idx++){
        gy.y[idx+1]=model->band_up[jdx+idx*model->nbands]-model->efermi;
        gy.idx[idx+1]=-1;
        gy.symbol[idx+1]=GRAPH_SYMB_NONE;
        gy.type=GRAPH_XY_TYPE;/*<- regular type!*/
}
gy.line=GRAPH_LINE_THICK;
gy.color=GRAPH_COLOR_DEFAULT;
dat_graph_add_y(gy,plot->graph);
}
xmin=-1.0*scale*plot->dos.ymax;
xmax=plot->band.xmax;
ymin=plot->band.ymin;
ymax=plot->band.ymax;
dat_graph_set_limits(xmin,xmax,ymin,ymax,plot->graph);/*HARD part*/
g_free(gy.y);
g_free(gy.idx);
g_free(gy.symbol);
/*force TRUE on label*/
if(plot->xtics <= 1) graph_set_xticks(TRUE,2,plot->graph);
else graph_set_xticks(TRUE,plot->xtics,plot->graph);
if(plot->ytics <= 1) graph_set_yticks(TRUE,2,plot->graph);
else graph_set_yticks(TRUE,plot->ytics,plot->graph);
}
/*************/
/* plot BAND */
/*************/
void draw_plot_band(struct model_pak *model,struct plot_pak *plot){
g_data_x gx;
g_data_y gy;
gint idx,jdx;
gdouble d;
gdouble xmin,xmax;
gdouble ymin,ymax;
/*avoid rare case _BUG_*/
if(model->nkpoints<1) return;
if(model->nbands<1) return;
/*always auto*/
xmin=plot->band.xmin;
xmax=plot->band.xmax;
ymin=plot->band.ymin;
ymax=plot->band.ymax;
if(ymin==ymax) {
	ymin=ymin-1.0;
	ymax=ymax+1.0;
}
/*add 5% to y limits for "readability"*/
d=ymax-ymin;
ymin=ymin-d*0.05;
ymax=ymax+d*0.05;
plot->graph=graph_new("BAND", model);
dat_graph_toggle_xaxis(plot->graph);
dat_graph_set_title("<big>Band diagram</big>",plot->graph);
//dat_graph_set_sub_title("<small>(From <b>GDIS</b> data)</small>",plot->graph);
dat_graph_set_x_title("k-points distance (Ang)",plot->graph);
dat_graph_set_y_title("Band eigenvalue (eV)",plot->graph);/*FIXME: add proper unit text!*/
dat_graph_set_type(GRAPH_XY_TYPE,plot->graph);
/*set x*/
gx.x_size=model->nkpoints;
gx.x=g_malloc(gx.x_size*sizeof(gdouble));
for(idx=0;idx<gx.x_size;idx++) gx.x[idx]=model->kpts_d[idx];
dat_graph_set_x(gx,plot->graph);
g_free(gx.x);
/*set y*/
gy.y_size=model->nkpoints;
gy.y=g_malloc(gy.y_size*sizeof(gdouble));
gy.idx=g_malloc(gy.y_size*sizeof(gint32));/*<- will NOT set any structure*/
gy.symbol=g_malloc(gy.y_size*sizeof(graph_symbol));
gy.mixed_symbol=FALSE;
gy.sym_color=NULL;
for(jdx=0;jdx<model->nbands;jdx++){
for(idx=0;idx<model->nkpoints;idx++){
        gy.y[idx]=model->band_up[jdx+idx*model->nbands]-model->efermi;
        gy.idx[idx]=-1;
        gy.symbol[idx]=GRAPH_SYMB_NONE;
}
gy.type=GRAPH_XY_TYPE;
gy.line=GRAPH_LINE_THICK;
gy.color=GRAPH_COLOR_DEFAULT;
dat_graph_add_y(gy,plot->graph);
}
dat_graph_set_limits(xmin,xmax,ymin,ymax,plot->graph);
g_free(gy.y);
g_free(gy.idx);
g_free(gy.symbol);
/*force TRUE on label*/
if(plot->xtics <= 1) graph_set_xticks(TRUE,2,plot->graph);
else graph_set_xticks(TRUE,plot->xtics,plot->graph);
if(plot->ytics <= 1) graph_set_yticks(TRUE,2,plot->graph);
else graph_set_yticks(TRUE,plot->ytics,plot->graph);
}
/******************/
/* plot frequency */
/******************/
void draw_plot_frequency(struct model_pak *model,struct plot_pak *plot){
g_data_x gx;
g_data_y gy;
gint idx,jdx;
gdouble xmin,xmax;
gdouble ymin,ymax;
/*avoid rare case _BUG_*/
if(model->nfreq<1) return;
/*always auto*/
xmin=plot->frequency.xmin;
xmax=plot->frequency.xmax;
ymin=plot->frequency.ymin;
ymax=plot->frequency.ymax;
if(ymin==ymax) {
	ymin=ymin-1.0;
	ymax=ymax+1.0;
}
/*add 5% to y limits for "readability"*/
ymax=ymax+(ymax-ymin)*0.05;
plot->graph=graph_new("FREQ", model);
dat_graph_set_title("<big>Vibrational frequency</big>",plot->graph);
//dat_graph_set_sub_title("<small>(From <b>GDIS</b> data)</small>",plot->graph);
dat_graph_set_x_title("Frequency (cm<sup>-1</sup>)",plot->graph);
dat_graph_set_y_title("Intensity (arb. unit)",plot->graph);
dat_graph_set_type(GRAPH_XY_TYPE,plot->graph);
/*set x*/
gx.x_size=1+model->nfreq*3;
/*we triple the data as each freq f,
  and intensity i is represented by: 
  {(f,0),(f,i),(f,0)} points.
  added one point: (0,0) for origin.*/
gx.x=g_malloc(gx.x_size*sizeof(gdouble));
gx.x[0]=0;
idx=1;jdx=0;
while(idx<gx.x_size){
	gx.x[idx]=plot->frequency.data[jdx];
	gx.x[idx+1]=plot->frequency.data[jdx];
	gx.x[idx+2]=plot->frequency.data[jdx];
	jdx+=2;
	idx+=3;
}
dat_graph_set_x(gx,plot->graph);
g_free(gx.x);
/*set y*/
gy.y_size=1+model->nfreq*3;
gy.y=g_malloc(gy.y_size*sizeof(gdouble));
gy.idx=g_malloc(gy.y_size*sizeof(gint32));/*<- will NOT set any structure*/
gy.symbol=g_malloc(gy.y_size*sizeof(graph_symbol));
gy.mixed_symbol=TRUE;/*symbol are mixed by default*/
gy.sym_color=NULL;
gy.y[0]=0;
idx=1;jdx=1;
while(idx<gx.x_size){
        gy.y[idx]=0.;
	gy.y[idx+1]=plot->frequency.data[jdx];
	gy.y[idx+2]=0.;
        gy.idx[idx]=-1;
	gy.idx[idx+1]=-1;
	gy.idx[idx+2]=-1;
        gy.symbol[idx]=GRAPH_SYMB_NONE;
	gy.symbol[idx+1]=GRAPH_SYMB_SQUARE;/*the actual peak*/
	gy.symbol[idx+2]=GRAPH_SYMB_NONE;
	jdx+=2;
	idx+=3;
}
gy.type=GRAPH_XY_TYPE;
dat_graph_set_limits(xmin,xmax,ymin,ymax,plot->graph);
gy.line=GRAPH_LINE_THICK;
gy.color=GRAPH_COLOR_DEFAULT;
dat_graph_add_y(gy,plot->graph);
g_free(gy.y);
g_free(gy.idx);
g_free(gy.symbol);
/*force TRUE on label*/
if(plot->xtics <= 1) graph_set_xticks(TRUE,2,plot->graph);
else graph_set_xticks(TRUE,plot->xtics,plot->graph);
if(plot->ytics <= 1) graph_set_yticks(FALSE,2,plot->graph);
else graph_set_yticks(FALSE,plot->ytics,plot->graph);
}
/****************************/
/* Generic plotting routine */
/****************************/
void plot_draw_graph(struct plot_pak *plot,struct task_pak *task){
	struct model_pak *data;
	g_assert(plot != NULL);
	data = plot->model;
	g_assert(data != NULL);
        /* plot the corresponding graph */
        switch(plot->type){
                case PLOT_ENERGY:
                        draw_plot_energy(data,plot);
                        return;
                case PLOT_FORCE:
			draw_plot_force(data,plot);
                        return;
                case PLOT_VOLUME:
                        draw_plot_volume(data,plot);
                        return;
                case PLOT_PRESSURE:
			draw_plot_pressure(data,plot);
			return;
                break;
                case PLOT_BAND:
			draw_plot_band(data,plot);
			return;
		break;
                case PLOT_DOS:
			draw_plot_dos(data,plot);
			return;
                case PLOT_BANDOS:
			draw_plot_bandos(data,plot);
			return;
		break;
                case PLOT_FREQUENCY:
			draw_plot_frequency(data,plot);
			return;
		case PLOT_RAMAN:
                fprintf(stdout,"ERROR: plot type not available (yet)!\n");
                break;
                case PLOT_NONE:
                default:
                fprintf(stdout,"ERROR: no plot type selected!\n");
        }
}

void plot_show_graph(struct plot_pak *plot){
	/* task cleanup: will have to re-process some data, maybe? */
}

