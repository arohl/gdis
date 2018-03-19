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
/* TODO plot complex data: (bands, dos, bandos, frequency, ...) */

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
	int i,j;
	gdouble x,y;
	struct model_pak *model;
	g_assert(plot != NULL);
	g_assert(plot->model != NULL);
	model=plot->model;

/* task cleanup will interpolate some graph (DOS/BAND + FREQ...) */
	/* TODO: DOS/BAND */
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
	if(plot->plot_mask==0) return;/*no data to crunch*/
	g_assert(plot->model != NULL);
/* simply use model provided in plot */
	model=plot->model;
	if(model->num_frames<2) return;/*no frame -> skip this part*/
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
}
/***************/
/* energy plot */
/***************/
void draw_plot_energy(struct model_pak *model,struct plot_pak *plot){
        /*draw energy for each ionic step iterations*/
	gdouble xmin,xmax;
	gdouble ymin,ymax;
	/* get limits */
        if(plot->auto_x==FALSE){
                xmin=plot->xmin;
                xmax=plot->xmax;
        }else{  
                xmin=plot->energy.xmin;
                xmax=plot->energy.xmax;
        }
        if(plot->auto_y==FALSE){
                ymin=plot->ymin;
                ymax=plot->ymax;
        }else{  
                ymin=plot->energy.ymin;
                ymax=plot->energy.ymax;
        }
        if(ymin==ymax) {
                ymin=ymin-1.0;
                ymax=ymax+1.0;
        }
        plot->graph=graph_new("ENERGY",model);
        graph_add_borned_data(plot->energy.size,plot->energy.data,xmin,xmax,ymin,ymax,GRAPH_REGULAR,plot->graph);
        /* set tics */
        /* due to assert xtics and ytics needs to be set > 1 even if they are unused */
	/* see FIXME in gl_graph.c */
        if(plot->xtics <= 1) graph_set_xticks(FALSE,2,plot->graph);
        else graph_set_xticks(TRUE,plot->xtics,plot->graph);
        if(plot->ytics <= 1) graph_set_yticks(FALSE,2,plot->graph);
        else graph_set_yticks(TRUE,plot->ytics,plot->graph);
}
/**************/
/* force plot */
/**************/
void draw_plot_force(struct model_pak *model,struct plot_pak *plot){
        /*draw force for each ionic step iterations*/
	gdouble xmin,xmax;
	gdouble ymin,ymax;
        /* get limits */
        if(plot->auto_x==FALSE){
		xmin=plot->xmin;
		xmax=plot->xmax;
        }else{
		xmin=plot->force.xmin;
		xmax=plot->force.xmax;
	}
	if(plot->auto_y==FALSE){
		ymin=plot->ymin;
		ymax=plot->ymax;
	}else{
		ymin=plot->force.ymin;
		ymax=plot->force.ymax;
	}
        if(ymin==ymax) {
                ymin=ymin-1.0;
                ymax=ymax+1.0;
        }
        plot->graph=graph_new("FORCE",model);
        graph_add_borned_data(plot->force.size,plot->force.data,xmin,xmax,ymin,ymax,GRAPH_REGULAR,plot->graph);
        /* set tics */
        if(plot->xtics <= 1) graph_set_xticks(FALSE,2,plot->graph);
        else graph_set_xticks(TRUE,plot->xtics,plot->graph);
        if(plot->ytics <= 1) graph_set_yticks(FALSE,2,plot->graph);
        else graph_set_yticks(TRUE,plot->ytics,plot->graph);
}
/***************/
/* plot volume */
/***************/
void draw_plot_volume(struct model_pak *model,struct plot_pak *plot){
        /*draw volume for each ionic step iterations*/
        gdouble xmin,xmax;
        gdouble ymin,ymax;
        /* get limits */
        if(plot->auto_x==FALSE){
                xmin=plot->xmin;
                xmax=plot->xmax;
        }else{  
                xmin=plot->volume.xmin;
                xmax=plot->volume.xmax;
        }
        if(plot->auto_y==FALSE){
                ymin=plot->ymin;
                ymax=plot->ymax;
        }else{  
                ymin=plot->volume.ymin;
                ymax=plot->volume.ymax;
        }
	if(ymin==ymax) {
		ymin=ymin-1.0;
		ymax=ymax+1.0;
	}
        plot->graph=graph_new("VOLUME",model);
        graph_add_borned_data(plot->volume.size,plot->volume.data,xmin,xmax,ymin,ymax,GRAPH_REGULAR,plot->graph);
        /* set tics */
        if(plot->xtics <= 1) graph_set_xticks(FALSE,2,plot->graph);
        else graph_set_xticks(TRUE,plot->xtics,plot->graph);
        if(plot->ytics <= 1) graph_set_yticks(FALSE,2,plot->graph);
        else graph_set_yticks(TRUE,plot->ytics,plot->graph);
}
/*****************/
/* plot pressure */
/*****************/
void draw_plot_pressure(struct model_pak *model,struct plot_pak *plot){
        /*draw volume for each ionic step iterations*/
        gdouble xmin,xmax;
        gdouble ymin,ymax;
        /* get limits */
        if(plot->auto_x==FALSE){
                xmin=plot->xmin;
                xmax=plot->xmax;
        }else{
                xmin=plot->pressure.xmin;
                xmax=plot->pressure.xmax;
        }
        if(plot->auto_y==FALSE){
                ymin=plot->ymin;
                ymax=plot->ymax;
        }else{
                ymin=plot->pressure.ymin;
                ymax=plot->pressure.ymax;
        }
        if(ymin==ymax) {
                ymin=ymin-1.0;
                ymax=ymax+1.0;
        }
        plot->graph=graph_new("PRESSURE",model);
        graph_add_borned_data(plot->pressure.size,plot->pressure.data,xmin,xmax,ymin,ymax,GRAPH_REGULAR,plot->graph);
        /* set tics */
        if(plot->xtics <= 1) graph_set_xticks(FALSE,2,plot->graph);
        else graph_set_xticks(TRUE,plot->xtics,plot->graph);
        if(plot->ytics <= 1) graph_set_yticks(FALSE,2,plot->graph);
        else graph_set_yticks(TRUE,plot->ytics,plot->graph);
}
/************/
/* plot DOS */
/************/
void draw_plot_dos(struct model_pak *model,struct plot_pak *plot){
        /*draw density of state (dos)*/
        gdouble xmin,xmax;
        gdouble ymin,ymax;
        /* get limits */
        if(plot->auto_x==FALSE){
                xmin=plot->xmin;
                xmax=plot->xmax;
        }else{
                xmin=plot->dos.xmin;
                xmax=plot->dos.xmax;
        }
        if(plot->auto_y==FALSE){
                ymin=plot->ymin;
                ymax=plot->ymax;
        }else{
                ymin=plot->dos.ymin;
                ymax=plot->dos.ymax;
        }
        if(ymin==ymax) {
                ymin=ymin-1.0;
                ymax=ymax+1.0;
        }
        plot->graph=graph_new("DOS",model);
        graph_add_borned_data(plot->dos.size,plot->dos.data,xmin,xmax,ymin,ymax,GRAPH_DOS,plot->graph);
        /* set tics */
        if(plot->xtics <= 1) graph_set_xticks(FALSE,2,plot->graph);
        else graph_set_xticks(TRUE,plot->xtics,plot->graph);
        if(plot->ytics <= 1) graph_set_yticks(FALSE,2,plot->graph);
        else graph_set_yticks(TRUE,plot->ytics,plot->graph);
}
/******************/
/* plot frequency */
/******************/
void draw_plot_frequency(struct model_pak *model,struct plot_pak *plot){
        /*draw all frequency values*/
        gdouble xmin,xmax;
        gdouble ymin,ymax;
        /* get limits */
        if(plot->auto_x==FALSE){
                xmin=plot->xmin;
                xmax=plot->xmax;
        }else{
                xmin=plot->frequency.xmin;
                xmax=plot->frequency.xmax;
        }
        if(plot->auto_y==FALSE){
                ymin=plot->ymin;
                ymax=plot->ymax;
        }else{
                ymin=plot->frequency.ymin;
                ymax=plot->frequency.ymax;
        }
        if(ymin==ymax) {
                ymin=ymin-1.0;
                ymax=ymax+1.0;
        }
        plot->graph=graph_new("FREQUENCY",model);
        graph_add_borned_data(plot->frequency.size,plot->frequency.data,xmin,xmax,ymin,ymax,GRAPH_FREQUENCY,plot->graph);
        /* set tics */
        if(plot->xtics <= 1) graph_set_xticks(FALSE,2,plot->graph);
        else graph_set_xticks(TRUE,plot->xtics,plot->graph);
        if(plot->ytics <= 1) graph_set_yticks(FALSE,2,plot->graph);
        else graph_set_yticks(TRUE,plot->ytics,plot->graph);
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
		fprintf(stdout,"ERROR: plot type not available (yet)!\n");
		break;
                case PLOT_DOS:
			draw_plot_dos(data,plot);
			return;
                case PLOT_BANDOS:
		fprintf(stdout,"ERROR: plot type not available (yet)!\n");
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

