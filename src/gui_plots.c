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

/* Plot simple data, such as energy, force, volume, etc. */

#include <math.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
//#include <time.h>
#include <unistd.h>

#include "gdis.h"
#include "graph.h"
#include "model.h"
#include "file.h"
#include "parse.h"
#include "render.h"
#include "quaternion.h"
#include "gui_shorts.h"
#include "interface.h"
#include "dialog.h"
#include "plots.h"
#include "gui_defs.h"

/* global pak structures */
extern struct sysenv_pak sysenv;
extern struct elem_pak elements[];
struct plot_pak *plot;
struct graph_ui *graph_ui;
/* text widgets */
GtkWidget *plot_xmin;
GtkWidget *plot_xmax;
GtkWidget *plot_ymin;
GtkWidget *plot_ymax;
/* radio buttons */
/*1st page*/
GtkWidget *b_1_e;
GtkWidget *b_1_f;
GtkWidget *b_1_v;
GtkWidget *b_1_p;
/*2nd page*/
GtkWidget *b_2_d;
GtkWidget *b_2_b;
GtkWidget *b_2_bd;
/*3rd page*/
GtkWidget *b_3_f;
GtkWidget *b_3_r;

#define GRAPH_UI (*graph_ui)

/************************/
/* plot toggle checkbox */
/************************/
void set_auto_x(void){
	switch(plot->type){
	case PLOT_ENERGY:
		gtk_entry_set_text(GTK_ENTRY(plot_xmin),g_strdup_printf("%lf",plot->energy.xmin));
		gtk_entry_set_text(GTK_ENTRY(plot_xmax),g_strdup_printf("%lf",plot->energy.xmax));
		break;
	case PLOT_FORCE:
		gtk_entry_set_text(GTK_ENTRY(plot_xmin),g_strdup_printf("%lf",plot->force.xmin));
		gtk_entry_set_text(GTK_ENTRY(plot_xmax),g_strdup_printf("%lf",plot->force.xmax));
		break;
	case PLOT_VOLUME:
		gtk_entry_set_text(GTK_ENTRY(plot_xmin),g_strdup_printf("%lf",plot->volume.xmin));
		gtk_entry_set_text(GTK_ENTRY(plot_xmax),g_strdup_printf("%lf",plot->volume.xmax));
		break;
	case PLOT_PRESSURE:
		gtk_entry_set_text(GTK_ENTRY(plot_xmin),g_strdup_printf("%lf",plot->pressure.xmin));
		gtk_entry_set_text(GTK_ENTRY(plot_xmax),g_strdup_printf("%lf",plot->pressure.xmax));
		break;
	case PLOT_DOS:
		gtk_entry_set_text(GTK_ENTRY(plot_xmin),g_strdup_printf("%lf",plot->dos.xmin));
		gtk_entry_set_text(GTK_ENTRY(plot_xmax),g_strdup_printf("%lf",plot->dos.xmax));
		break;
	case PLOT_FREQUENCY:
		gtk_entry_set_text(GTK_ENTRY(plot_xmin),g_strdup_printf("%lf",plot->frequency.xmin));
		gtk_entry_set_text(GTK_ENTRY(plot_xmax),g_strdup_printf("%lf",plot->frequency.xmax));
		break;
	case PLOT_BAND:
		gtk_entry_set_text(GTK_ENTRY(plot_xmin),g_strdup_printf("%lf",plot->band.xmin));
		gtk_entry_set_text(GTK_ENTRY(plot_xmax),g_strdup_printf("%lf",plot->band.xmax));
		break;
	case PLOT_BANDOS:
		gtk_entry_set_text(GTK_ENTRY(plot_xmin),g_strdup_printf("%lf",0.0));
		gtk_entry_set_text(GTK_ENTRY(plot_xmax),g_strdup_printf("%lf",10.0));
		break;
	case PLOT_RAMAN:
	case PLOT_NONE:
	default:
		gtk_entry_set_text(GTK_ENTRY(plot_xmin),g_strdup_printf("%lf",plot->xmin));
		gtk_entry_set_text(GTK_ENTRY(plot_xmax),g_strdup_printf("%lf",plot->xmax));
	}
}
void set_auto_y(void){
        switch(plot->type){
        case PLOT_ENERGY:
                gtk_entry_set_text(GTK_ENTRY(plot_ymin),g_strdup_printf("%lf",plot->energy.ymin));
                gtk_entry_set_text(GTK_ENTRY(plot_ymax),g_strdup_printf("%lf",plot->energy.ymax));
                break;
        case PLOT_FORCE:
                gtk_entry_set_text(GTK_ENTRY(plot_ymin),g_strdup_printf("%lf",plot->force.ymin));
                gtk_entry_set_text(GTK_ENTRY(plot_ymax),g_strdup_printf("%lf",plot->force.ymax));
                break;
        case PLOT_VOLUME:
                gtk_entry_set_text(GTK_ENTRY(plot_ymin),g_strdup_printf("%lf",plot->volume.ymin));
                gtk_entry_set_text(GTK_ENTRY(plot_ymax),g_strdup_printf("%lf",plot->volume.ymax));
                break;
        case PLOT_PRESSURE:
                gtk_entry_set_text(GTK_ENTRY(plot_ymin),g_strdup_printf("%lf",plot->pressure.ymin));
                gtk_entry_set_text(GTK_ENTRY(plot_ymax),g_strdup_printf("%lf",plot->pressure.ymax));
                break;
        case PLOT_DOS:
		gtk_entry_set_text(GTK_ENTRY(plot_ymin),g_strdup_printf("%lf",plot->dos.ymin));
		gtk_entry_set_text(GTK_ENTRY(plot_ymax),g_strdup_printf("%lf",plot->dos.ymax));
		break;
        case PLOT_FREQUENCY:
		gtk_entry_set_text(GTK_ENTRY(plot_ymin),g_strdup_printf("%lf",plot->frequency.ymin));
		gtk_entry_set_text(GTK_ENTRY(plot_ymax),g_strdup_printf("%lf",plot->frequency.ymax));
		break;
	case PLOT_BAND:
		gtk_entry_set_text(GTK_ENTRY(plot_ymin),g_strdup_printf("%lf",plot->band.ymin));
		gtk_entry_set_text(GTK_ENTRY(plot_ymax),g_strdup_printf("%lf",plot->band.ymax));
		break;
	case PLOT_BANDOS:
		gtk_entry_set_text(GTK_ENTRY(plot_ymin),g_strdup_printf("%lf",plot->band.ymin));
		gtk_entry_set_text(GTK_ENTRY(plot_ymax),g_strdup_printf("%lf",plot->band.ymax));
		break;
	case PLOT_RAMAN:
        case PLOT_NONE:
        default:
                gtk_entry_set_text(GTK_ENTRY(plot_ymin),g_strdup_printf("%lf",plot->ymin));
                gtk_entry_set_text(GTK_ENTRY(plot_ymax),g_strdup_printf("%lf",plot->ymax));
        }
}
void auto_x_toggle(GtkWidget *w, GtkWidget *box){
        if(plot->auto_x==TRUE){
                /* disallow modification */
                gtk_widget_set_sensitive(plot_xmin,FALSE);
                gtk_widget_set_sensitive(plot_xmax,FALSE);
		/* fill in extreme */
		set_auto_x();
        } else {
                /* allow modification */
                gtk_widget_set_sensitive(plot_xmin,TRUE);
                gtk_widget_set_sensitive(plot_xmax,TRUE);
        }
}
void auto_y_toggle(GtkWidget *w, GtkWidget *box){
        if(plot->auto_y==TRUE){
                /* disallow modification */
                gtk_widget_set_sensitive(plot_ymin,FALSE);
                gtk_widget_set_sensitive(plot_ymax,FALSE);
		/* fill in extrema */
		set_auto_y();
        } else {
                /* allow modification */
                gtk_widget_set_sensitive(plot_ymin,TRUE);
                gtk_widget_set_sensitive(plot_ymax,TRUE);
        }
}
/*******************/
/* plot type radio */
/*******************/
void plot_energy(){
	/* switched to PLOT_ENERGY */
	plot->type=PLOT_ENERGY;
	plot->plot_sel=plot->plot_sel&240;/*reset 1st page*/
	plot->plot_sel+=PLOT_ENERGY;/*selected*/
	plot->xmin=plot->energy.xmin;
	plot->xmax=plot->energy.xmax;
	plot->ymin=plot->energy.ymin;
	plot->ymax=plot->energy.ymax;
	plot->auto_x=TRUE;
	plot->auto_y=TRUE;
	auto_x_toggle(NULL,NULL);
	auto_y_toggle(NULL,NULL);
	gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(b_1_e),TRUE);
	gui_relation_update(plot->model);
}
void plot_force(){
	plot->type=PLOT_FORCE;
        plot->plot_sel=plot->plot_sel&240;/*reset 1st page*/
        plot->plot_sel+=PLOT_FORCE;/*selected*/
        plot->xmin=plot->force.xmin;
        plot->xmax=plot->force.xmax;
        plot->ymin=plot->force.ymin;
        plot->ymax=plot->force.ymax;
        plot->auto_x=TRUE;
        plot->auto_y=TRUE;
        auto_x_toggle(NULL,NULL);
        auto_y_toggle(NULL,NULL);
	gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(b_1_f),TRUE);
	gui_relation_update(plot->model);
}
void plot_volume(){
	plot->type=PLOT_VOLUME;
        plot->plot_sel=plot->plot_sel&240;/*reset 1st page*/
        plot->plot_sel+=PLOT_VOLUME;/*selected*/
        plot->xmin=plot->volume.xmin;
        plot->xmax=plot->volume.xmax;
        plot->ymin=plot->volume.ymin;
        plot->ymax=plot->volume.ymax;
        plot->auto_x=TRUE;
        plot->auto_y=TRUE;
        auto_x_toggle(NULL,NULL);
        auto_y_toggle(NULL,NULL);
	gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(b_1_v),TRUE);
	gui_relation_update(plot->model);
}
void plot_pressure(){
	plot->type=PLOT_PRESSURE;
        plot->plot_sel=plot->plot_sel&240;/*reset 1st page*/
        plot->plot_sel+=PLOT_PRESSURE;/*selected*/
        plot->xmin=plot->pressure.xmin;
        plot->xmax=plot->pressure.xmax;
        plot->ymin=plot->pressure.ymin;
        plot->ymax=plot->pressure.ymax;
        plot->auto_x=TRUE;
        plot->auto_y=TRUE;
        auto_x_toggle(NULL,NULL);
        auto_y_toggle(NULL,NULL);
	gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(b_1_p),TRUE);
	gui_relation_update(plot->model);
}
void plot_band(){
        plot->type=PLOT_BAND;
	plot->plot_sel=plot->plot_sel&207;/*reset 2nd page*/
	plot->plot_sel+=PLOT_BAND;/*selected*/
        plot->xmin=plot->band.xmin;
        plot->xmax=plot->band.xmax;
        plot->ymin=plot->band.ymin;
        plot->ymax=plot->band.ymax;
        plot->auto_x=TRUE;
        plot->auto_y=TRUE;
        auto_x_toggle(NULL,NULL);
        auto_y_toggle(NULL,NULL);
	gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(b_2_b),TRUE);
	gui_relation_update(plot->model);
}
void plot_dos(){
        plot->type=PLOT_DOS;
	plot->plot_sel=plot->plot_sel&207;/*reset 2nd page*/
	plot->plot_sel+=PLOT_DOS;/*selected*/
        plot->xmin=plot->dos.xmin;
        plot->xmax=plot->dos.xmax;
        plot->ymin=plot->dos.ymin;
        plot->ymax=plot->dos.ymax;
        plot->auto_x=TRUE;
        plot->auto_y=TRUE;
        auto_x_toggle(NULL,NULL);
        auto_y_toggle(NULL,NULL);
	gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(b_2_d),TRUE);
	gui_relation_update(plot->model);
}
void plot_bandos(){
        plot->type=PLOT_BANDOS;
	plot->plot_sel=plot->plot_sel&207;/*reset 2nd page*/
	plot->plot_sel+=PLOT_BANDOS;/*selected*/
        plot->xmin=0.;
        plot->xmax=10.;
        plot->ymin=plot->band.ymin;
        plot->ymax=plot->band.ymax;
        plot->auto_x=TRUE;
        plot->auto_y=TRUE;
        auto_x_toggle(NULL,NULL);
        auto_y_toggle(NULL,NULL);
	gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(b_2_bd),TRUE);
	gui_relation_update(plot->model);
}
void plot_frequency(){
        plot->type=PLOT_FREQUENCY;
	plot->plot_sel=plot->plot_sel&63;/*reset 3rd page*/
	plot->plot_sel+=PLOT_FREQUENCY;/*selected*/
        plot->xmin=plot->frequency.xmin;
        plot->xmax=plot->frequency.xmax;
        plot->ymin=plot->frequency.ymin;
        plot->ymax=plot->frequency.ymax;
        plot->auto_x=TRUE;
        plot->auto_y=TRUE;
        auto_x_toggle(NULL,NULL);
        auto_y_toggle(NULL,NULL);
	gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(b_3_f),TRUE);
	gui_relation_update(plot->model);
}
void plot_raman(){
        plot->type=PLOT_RAMAN;
	plot->plot_sel=plot->plot_sel&63;/*reset 3rd page*/
	plot->plot_sel+=PLOT_FREQUENCY;/*selected*/
        plot->xmin=plot->raman.xmin;
        plot->xmax=plot->raman.xmax;
        plot->ymin=plot->raman.ymin;
        plot->ymax=plot->raman.ymax;
        plot->auto_x=TRUE;
        plot->auto_y=TRUE;
        auto_x_toggle(NULL,NULL);
        auto_y_toggle(NULL,NULL);
	gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(b_3_r),TRUE);
	gui_relation_update(plot->model);
}
void plot_none(){
	plot->type=PLOT_NONE;
	/* no plot type available */
	plot->xmin=0.;
        plot->xmax=1.;
        plot->ymin=0.;
        plot->ymax=1.;
        plot->auto_x=TRUE;
        plot->auto_y=TRUE;
        auto_x_toggle(NULL,NULL);
        auto_y_toggle(NULL,NULL);
	gui_relation_update(plot->model);
}

void plot_page_switch(GtkNotebook *notebook,GtkWidget *page,guint page_num,gpointer user_data){
//	GtkWidget *child=gtk_notebook_get_nth_page(notebook,page_num);/*in case we need/want to resize the current page*/
	switch(page_num){/*pb is we don't know who is selected in page*/
		case 0:/*Dynamics*/
			if(plot->plot_sel&PLOT_ENERGY) plot_energy();
			else if(plot->plot_sel&PLOT_FORCE) plot_force();
			else if(plot->plot_sel&PLOT_VOLUME) plot_volume();
			else if(plot->plot_sel&PLOT_PRESSURE) plot_pressure();
			else {/*no plot selected, but plot_selection possible?*/
			        if(!(plot->plot_mask&(15))) plot_none();
			        else if(plot->plot_mask&PLOT_ENERGY) plot_energy();
			        else if(plot->plot_mask&PLOT_FORCE) plot_force();
			        else if(plot->plot_mask&PLOT_VOLUME) plot_volume();
			        else if(plot->plot_mask&PLOT_PRESSURE) plot_pressure();
			}
			break;
		case 1:/*Electronic*/
			if(plot->plot_sel&PLOT_BAND) plot_band();
			else if(plot->plot_sel&PLOT_DOS) plot_dos();
			else if((plot->plot_sel&PLOT_DOS)&&(plot->plot_sel&PLOT_BAND)) plot_bandos();
			else {/*no plot selected, but plot_selection possible?*/
                                if(!(plot->plot_mask&(48))) plot_none();
				else if(plot->plot_mask&PLOT_DOS) plot_dos();
				else if(plot->plot_mask&PLOT_BAND) plot_band();
			}
			break;
		case 2:/*Frequency*/
			if(plot->plot_sel&PLOT_FREQUENCY) plot_frequency();
			else if(plot->plot_sel&PLOT_RAMAN) plot_raman();
			else {/*no plot selected, but plot_selection possible?*/
                                if(!(plot->plot_mask&(192))) plot_none();
                                else if(plot->plot_mask&PLOT_FREQUENCY) plot_frequency();
                                else if(plot->plot_mask&PLOT_RAMAN) plot_raman();
			}
			break;
		default:
			plot_none();
	}
}
/* prepare plot structure */
void plot_initialize(struct model_pak *data){
	int i,j;
	gdouble y;
	/* init the plot pak */
	if (!data) return;
	if(plot!=NULL) g_free(plot);/* FIXME: safe? */
	plot=(struct plot_pak *)g_malloc(sizeof(struct plot_pak));
	if(plot==NULL) return;
	plot->auto_x=TRUE;
	plot->auto_y=TRUE;
	plot->xmin=0.0;
	plot->xmax=1.0;
	plot->ymin=0.0;
	plot->ymax=1.0;
	plot->xtics=1.0;
	plot->ytics=1.0;
	plot->graph=data->graph_active;
	plot->data_changed=TRUE;/*forces a reload of data*/
	plot->energy.data=NULL;
	plot->force.data=NULL;
	plot->volume.data=NULL;
	plot->pressure.data=NULL;
	plot->band.data=NULL;
	plot->dos.data=NULL;
	plot->frequency.data=NULL;
	plot->ndos=0;
	plot->nbands=0;
	plot->nfreq=0;
	plot->nraman=0;
	plot->plot_mask=PLOT_NONE;
	plot->plot_sel=0;
	plot->task=NULL;
	/*dynamics*/
	if(data->num_frames>1){
		if(property_lookup("Energy",data)) plot->plot_mask+=PLOT_ENERGY;
		if(property_lookup("Force",data)) plot->plot_mask+=PLOT_FORCE;
		if(property_lookup("Pressure",data)||property_lookup("Stress",data)) plot->plot_mask+=PLOT_PRESSURE;
		if(data->periodic == 3) plot->plot_mask+=PLOT_VOLUME;
		if(plot->plot_mask==PLOT_NONE) plot->type=PLOT_NONE;
		else if((plot->plot_mask&PLOT_ENERGY)!=0) plot->type=PLOT_ENERGY;
		else if((plot->plot_mask&PLOT_FORCE)!=0) plot->type=PLOT_FORCE;
		else if((plot->plot_mask&PLOT_VOLUME)!=0) plot->type=PLOT_VOLUME;
		else if((plot->plot_mask&PLOT_PRESSURE)!=0) plot->type=PLOT_PRESSURE;
	}
	/*electronics*/
	if(data->ndos>1) {
		/*we have some dos... preload (will be scaled later on a thread)*/
		plot->plot_mask+=PLOT_DOS;
		plot->ndos=data->ndos;
		plot->dos.size=plot->ndos*2;
		if(plot->dos.data!=NULL) {
			g_free(plot->dos.data);
			plot->dos.data=NULL;
		}
		plot->dos.data=(gdouble *)g_malloc(plot->dos.size*sizeof(gdouble));
		
	}
	if((data->nbands>1)&&(data->band_up!=NULL)) {
		plot->plot_mask+=PLOT_BAND;
		plot->nbands=data->nbands;
		/*the sizes!*/
		plot->band.size=data->nkpoints;
		if(data->spin_polarized) plot->band.data=(gdouble *)g_malloc(((1+plot->nbands)*plot->band.size+4)*sizeof(gdouble));
		else plot->band.data=(gdouble *)g_malloc(((1+2*plot->nbands)*plot->band.size+4)*sizeof(gdouble));
	}
	if(!((plot->plot_mask&PLOT_BANDOS)^PLOT_BANDOS)) {
		/*bandos type is possible*/
	}
	/*frequency*/
	if(data->have_frequency==TRUE) {
		/*this don't need to be loaded from outside.. it won't change*/
		plot->plot_mask+=PLOT_FREQUENCY;
                plot->nfreq=data->nfreq;
		plot->frequency.size=plot->nfreq*2;
                if(plot->frequency.data!=NULL) {
                        g_free(plot->frequency.data);
                        plot->frequency.data=NULL;
                }
                plot->frequency.data=(gdouble *)g_malloc(plot->frequency.size*sizeof(gdouble));
                plot->frequency.xmin=0.0;/*the default is to start plotting frequency in [0,ymax+100]*/
                plot->frequency.xmax=data->freq[0];
		plot->frequency.ymin=0.;/*default is relative intensity ie [0.0,1.0]*/
		plot->frequency.ymax=1.;
		j=0;
                for(i=0;i<data->nfreq;i++){
                        y=data->freq[i];
                        if(y>plot->frequency.xmax) plot->frequency.xmax=y;
                        plot->frequency.data[j]=y;/*this data is actually x*/
			if(data->freq_intens==NULL) plot->frequency.data[j+1]=1.0;
			else {
				plot->frequency.data[j+1]=data->freq_intens[i];
				if(data->freq_intens[i]>plot->frequency.ymax) plot->frequency.ymax=data->freq_intens[i];
			}
			j=j+2;
		}
		plot->frequency.xmax=((gint)(plot->frequency.xmax/1000)+1)*1000.0;
	}
}
/************************/
/* plots dialog cleanup */
/************************/
void plot_cleanup(struct model_pak *model){
	/* TODO: no cleanup to do? */
	g_assert(model != NULL);
	/*cleanup plot mess*/
	if(plot!=NULL){
		/* plot->graph is destroyed when closing graph */
		g_free(plot->energy.data);
		g_free(plot->force.data);
		g_free(plot->volume.data);
		g_free(plot->pressure.data);
		g_free(plot->band.data);/*so what?*/
		g_free(plot->dos.data);
		g_free(plot->frequency.data);
		g_free(plot);
		plot=NULL;
	}
	model->plots=FALSE;
}
void plot_quit(GtkWidget *w, gpointer data){
	struct model_pak *model=data;
	g_assert(model != NULL);
	/*cleanup plot mess*/
	plot_cleanup(model);
	dialog_destroy(w,model);
}
/******************/
/* start plotting */
/******************/
void exec_plot(){
	FILE *fp;
	struct model_pak *data;
	gpointer camera=NULL;
	g_assert(plot != NULL);
	data = sysenv.active_model;
	g_assert(data != NULL);
	/*get limit values if not auto*/
	if(plot->type==PLOT_NONE) return;
	camera = camera_dup(data->camera);/*save camera (from gui_animate.c)*/
	if(plot->auto_x==FALSE){
		plot->xmin=str_to_float(gtk_entry_get_text(GTK_ENTRY(plot_xmin)));
		plot->xmax=str_to_float(gtk_entry_get_text(GTK_ENTRY(plot_xmax)));
		/*ensure xmin<xmax*/
		if(plot->xmin>=plot->xmax) plot->xmin=plot->xmax;
	}
	if(plot->auto_y==FALSE){
		plot->ymin=str_to_float(gtk_entry_get_text(GTK_ENTRY(plot_ymin)));
		plot->ymax=str_to_float(gtk_entry_get_text(GTK_ENTRY(plot_ymax)));
		/*ensure ymin<ymax*/
		if(plot->ymin>=plot->ymax) plot->ymin=plot->ymax;
	}
        data->locked=TRUE;/*useful with task?*/
	if(plot->data_changed){
	        plot->model=data;
		task_system_new("PLOT-INIT",&plot_load_data,plot,&plot_prepare_data,plot,data);
	}
/* draw the graph in task */
	task_system_new("PLOT-DRAW",&plot_draw_graph,plot,&plot_show_graph,plot,data);
        data->locked=FALSE;
/*finalize drawing graph*/
/* There is a _BUG_ in QE redisplay that makes the atom distances shrink for each
 * plot that is displayed after the first one.
 * The only solution I found is to reload the current frame. (OVHPA)*/
if(data->id==QE_OUT){
        fp=fopen(data->filename, "r");
        if (!fp){
                gui_text_show(ERROR,g_strdup_printf("I/O ERROR: can't open!\n"));
                return;
        }
        read_raw_frame(fp,data->cur_frame,data);
	fclose(fp);
}
	tree_model_refresh(data);
	model_prep(data);
/*restore camera*/
        data->camera = camera;
        g_free(data->camera_default);
        data->camera_default = camera;
/*update everything*/
        redraw_canvas(ALL);
	
}
/******************/
/* configure plot */
/******************/
void gui_plots_dialog(void){
/* main dialog for drawing plots */
	gchar *tmp, *title;
	gpointer dialog;
	GtkWidget *window, *frame, *vbox, *hbox;
	GtkWidget *notebook, *page, *label;
	/* special */
	GtkWidget *ionic_box, *electronic_box, *frequency_box;
	struct model_pak *data;
	gpointer camera=NULL;
	gint cur_frame;
/* checks */
	data = sysenv.active_model;
	if (!data) return;
	/*Only one active at a time*/
	if (data->plots) return;
	data->plots=TRUE;
	cur_frame=data->cur_frame;
/* initialization */
	camera = camera_dup(data->camera);/*save camera (from gui_animate.c)*/
	gui_mode_switch(FREE);
	plot=NULL;
	plot_initialize(data);
/*1- lock and try to NOT update model*/
        data->redraw=0;
        data->locked=TRUE;
	sysenv.refresh_dialog=FALSE;
	plot->model=data;	
/*2- do some work*/
	task_system_new("PLOT-INIT",&plot_load_data,plot,&plot_prepare_data,plot,data);
/*3- unlock model*/
//	data->locked=FALSE;
//	model_content_refresh(data);
/* dialog setup */
	title = g_strdup_printf("Plots: %s", data->basename);
	dialog = dialog_request(PLOTS, title, NULL, plot_cleanup, data);
	g_free(title);
	window = dialog_window(dialog);
/* notebook frame */
	frame = gtk_frame_new(NULL);
	gtk_box_pack_start(GTK_BOX(GTK_DIALOG(window)->vbox), frame, FALSE, FALSE, 0);
	gtk_container_set_border_width(GTK_CONTAINER(frame), PANEL_SPACING);
/* create notebook */
	notebook = gtk_notebook_new();
	gtk_notebook_set_tab_pos(GTK_NOTEBOOK(notebook), GTK_POS_TOP);
	gtk_container_add(GTK_CONTAINER(frame), notebook);
	gtk_notebook_set_show_border(GTK_NOTEBOOK(notebook), TRUE);
/* --- page 1 -> ENERGY/FORCE/VOLUME/PRESSURE */
	page = gtk_vbox_new(FALSE, PANEL_SPACING);
	label = gtk_label_new("Dynamics");
	gtk_notebook_append_page(GTK_NOTEBOOK(notebook),page,label);
/* cycle options */
	frame = gtk_frame_new(NULL);
	gtk_box_pack_start(GTK_BOX(page),frame,FALSE,FALSE,0);
/* create a sensitive vbox in the frame */
	ionic_box = gtk_vbox_new(TRUE, 0);
	gtk_container_add(GTK_CONTAINER(frame), ionic_box);
/* num frames */
	hbox = gtk_hbox_new(FALSE, 0);
	gtk_container_set_border_width(GTK_CONTAINER(hbox), PANEL_SPACING);
	gtk_box_pack_start(GTK_BOX(ionic_box), hbox, FALSE, FALSE, 0);
	label = gtk_label_new("Ionic steps:");
	gtk_box_pack_start(GTK_BOX(hbox), label, FALSE, FALSE, 0);
	tmp = g_strdup_printf("%9d", data->num_frames);
	label = gtk_label_new(tmp);
	g_free(tmp);
	gtk_box_pack_end(GTK_BOX(hbox), label, FALSE, FALSE, 0);
/* Create a new vbox */
	vbox = gtk_vbox_new(FALSE,0);
	gtk_box_pack_start(GTK_BOX(ionic_box), vbox, TRUE, TRUE, 0);
	gtk_container_set_border_width(GTK_CONTAINER(vbox), PANEL_SPACING);
/* radio buttons */
	new_radio_group(0, ionic_box, FF);
/*by default any data not available _now_ will not be available later*/
	b_1_e = add_radio_button("Energy / step", (gpointer) plot_energy,NULL);
	if (plot->type == PLOT_ENERGY) gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(b_1_e), TRUE);
	if (!(plot->plot_mask&PLOT_ENERGY)) gtk_widget_set_sensitive(b_1_e,FALSE);
	b_1_f = add_radio_button("Forces / step", (gpointer) plot_force,NULL);
	if (plot->type == PLOT_FORCE) gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(b_1_f), TRUE);
	if (!(plot->plot_mask&PLOT_FORCE)) gtk_widget_set_sensitive(b_1_f,FALSE);
	b_1_v = add_radio_button("Volume / step", (gpointer) plot_volume,NULL);
	if (plot->type == PLOT_VOLUME) gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(b_1_v), TRUE);
	if (!(plot->plot_mask&PLOT_VOLUME)) gtk_widget_set_sensitive(b_1_v,FALSE);
	b_1_p = add_radio_button("Pressure / step", (gpointer) plot_pressure,NULL);
	if (plot->type == PLOT_PRESSURE) gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(b_1_p), TRUE);
	if (!(plot->plot_mask&PLOT_PRESSURE)) gtk_widget_set_sensitive(b_1_p,FALSE);
/*Set availability of the whole selection TODO: base on plot->plot_mask?*/
	if(data->num_frames>1) gtk_widget_set_sensitive(ionic_box, TRUE);
	else gtk_widget_set_sensitive(ionic_box, FALSE);
/* --- page 2 -> TODO: DOS/BAND */
	page = gtk_vbox_new(FALSE, PANEL_SPACING);
	label = gtk_label_new("Electronic");
	gtk_notebook_append_page(GTK_NOTEBOOK(notebook), page, label);
/* cycle options */
	frame = gtk_frame_new(NULL);
	gtk_box_pack_start(GTK_BOX(page),frame,FALSE,FALSE,0);
/* create another sensitive vbox */
	electronic_box = gtk_vbox_new(TRUE, 0);
	gtk_container_add(GTK_CONTAINER(frame), electronic_box);
/* add an hbox with DOS information */
	hbox = gtk_hbox_new(FALSE, 0);
	gtk_container_set_border_width(GTK_CONTAINER(hbox), PANEL_SPACING);
	gtk_box_pack_start(GTK_BOX(electronic_box), hbox, FALSE, FALSE, 0);
	label = gtk_label_new("DOS present:");
	gtk_box_pack_start(GTK_BOX(hbox), label, FALSE, FALSE, 0);
	tmp = g_strdup_printf("%9d", plot->ndos);
	label = gtk_label_new(tmp);
	g_free(tmp);
	gtk_box_pack_end(GTK_BOX(hbox), label, FALSE, FALSE, 0);
/* another hbox with BAND information */
	hbox = gtk_hbox_new(FALSE, 0);
	gtk_container_set_border_width(GTK_CONTAINER(hbox), PANEL_SPACING);
	gtk_box_pack_start(GTK_BOX(electronic_box), hbox, FALSE, FALSE, 0);
	label = gtk_label_new("BAND present:");
	gtk_box_pack_start(GTK_BOX(hbox), label, FALSE, FALSE, 0);
	tmp = g_strdup_printf("%9d", plot->nbands);
	label = gtk_label_new(tmp);
	g_free(tmp);
	gtk_box_pack_end(GTK_BOX(hbox), label, FALSE, FALSE, 0);
/* Create a new vbox */
	vbox = gtk_vbox_new(FALSE,0);
	gtk_box_pack_start(GTK_BOX(electronic_box), vbox, TRUE, TRUE, 0);
	gtk_container_set_border_width(GTK_CONTAINER(vbox), PANEL_SPACING);
/* radio buttons */
	new_radio_group(0, electronic_box, FF);
	b_2_d = add_radio_button("Density Of States", (gpointer) plot_dos,NULL);
	if (plot->type == PLOT_DOS) gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(b_2_d), TRUE);
	if (!(plot->plot_mask&PLOT_DOS)) gtk_widget_set_sensitive(b_2_d,FALSE);
	b_2_b = add_radio_button("Band structure", (gpointer) plot_band,NULL);
	if (plot->type == PLOT_BAND) gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(b_2_b), TRUE);
	if (!(plot->plot_mask&PLOT_BAND)) gtk_widget_set_sensitive(b_2_b,FALSE);
	b_2_bd = add_radio_button("Band / DOS", (gpointer) plot_bandos,NULL);
	if (plot->type == PLOT_BANDOS) gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(b_2_bd), TRUE);
	if ((plot->plot_mask&PLOT_BANDOS)^PLOT_BANDOS) gtk_widget_set_sensitive(b_2_bd,FALSE);
/* Set availability */
/* TODO: DOS/BAND availability */
	if (plot->plot_mask&PLOT_BANDOS) gtk_widget_set_sensitive(electronic_box, TRUE);
	else gtk_widget_set_sensitive(electronic_box, FALSE);
/* --- Page 3 -> TODO: Frequency */
        page = gtk_vbox_new(FALSE, PANEL_SPACING);
        label = gtk_label_new("Frequency");
        gtk_notebook_append_page(GTK_NOTEBOOK(notebook), page, label);
/* cycle options */
        frame = gtk_frame_new(NULL);
        gtk_box_pack_start(GTK_BOX(page),frame,FALSE,FALSE,0);
/* create another sensitive vbox */
        frequency_box = gtk_vbox_new(TRUE, 0);
        gtk_container_add(GTK_CONTAINER(frame), frequency_box);
/* add an hbox with FREQ information */
        hbox = gtk_hbox_new(FALSE, 0);
        gtk_container_set_border_width(GTK_CONTAINER(hbox), PANEL_SPACING);
        gtk_box_pack_start(GTK_BOX(frequency_box), hbox, FALSE, FALSE, 0);
        label = gtk_label_new("FREQs present:");
        gtk_box_pack_start(GTK_BOX(hbox), label, FALSE, FALSE, 0);
        tmp = g_strdup_printf("%9d", plot->nfreq);
        label = gtk_label_new(tmp);
        g_free(tmp);
        gtk_box_pack_end(GTK_BOX(hbox), label, FALSE, FALSE, 0);
/* another hbox with RAMAN information */
        hbox = gtk_hbox_new(FALSE, 0);
        gtk_container_set_border_width(GTK_CONTAINER(hbox), PANEL_SPACING);
        gtk_box_pack_start(GTK_BOX(frequency_box), hbox, FALSE, FALSE, 0);
        label = gtk_label_new("RAMAN present:");
        gtk_box_pack_start(GTK_BOX(hbox), label, FALSE, FALSE, 0);
        tmp = g_strdup_printf("%9d", plot->nraman);
        label = gtk_label_new(tmp);
        g_free(tmp);
        gtk_box_pack_end(GTK_BOX(hbox), label, FALSE, FALSE, 0);
/* Create a new vbox */
        vbox = gtk_vbox_new(FALSE,0);
        gtk_box_pack_start(GTK_BOX(frequency_box), vbox, TRUE, TRUE, 0);
        gtk_container_set_border_width(GTK_CONTAINER(vbox), PANEL_SPACING);
/* radio buttons */
        new_radio_group(0, frequency_box, FF);
        b_3_f = add_radio_button("Vibrational", (gpointer) plot_frequency,NULL);
        if (plot->type == PLOT_FREQUENCY) gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(b_3_f), TRUE);
	if (!(plot->plot_mask&PLOT_FREQUENCY)) gtk_widget_set_sensitive(b_3_f,FALSE);
        b_3_r = add_radio_button("Raman", (gpointer) plot_raman,NULL);
        if (plot->type == PLOT_RAMAN) gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(b_3_r), TRUE);
	if (!(plot->plot_mask&PLOT_RAMAN)) gtk_widget_set_sensitive(b_3_r,FALSE);
/* Set availability */
/* TODO: FREQUENCY/RAMAN availability */
	if (plot->plot_mask&PLOT_FREQUENCY) gtk_widget_set_sensitive(frequency_box, TRUE);
	else gtk_widget_set_sensitive(frequency_box, FALSE);
/* --- Outside of notebook */
	frame = gtk_frame_new(NULL);
	gtk_box_pack_start(GTK_BOX(GTK_DIALOG(window)->vbox), frame, FALSE, FALSE, 0);
	vbox = gtk_vbox_new(FALSE, PANEL_SPACING);
	gtk_container_add(GTK_CONTAINER(frame), vbox);
/* plot -> XMIN */
	hbox = gtk_hbox_new(FALSE, 0);
	gtk_box_pack_start(GTK_BOX(vbox), hbox, FALSE, FALSE, PANEL_SPACING);
	label = gtk_label_new(g_strdup_printf("X min"));
	gtk_box_pack_start(GTK_BOX(hbox), label, FALSE, FALSE, 0);
	plot_xmin = gtk_entry_new();
	gtk_entry_set_text(GTK_ENTRY(plot_xmin),g_strdup_printf("%-9.6f",plot->xmin));
	gtk_box_pack_end(GTK_BOX(hbox),plot_xmin, FALSE, FALSE, 0);
/* plot -> XMAX */
	hbox = gtk_hbox_new(FALSE, 0);
	gtk_box_pack_start(GTK_BOX(vbox), hbox, FALSE, FALSE, PANEL_SPACING);
	label = gtk_label_new(g_strdup_printf("X max"));
	gtk_box_pack_start(GTK_BOX(hbox), label, FALSE, FALSE, 0);
	plot_xmax = gtk_entry_new();
	gtk_entry_set_text(GTK_ENTRY(plot_xmax),g_strdup_printf("%-9.6f",plot->xmax));
	gtk_box_pack_end(GTK_BOX(hbox),plot_xmax, FALSE, FALSE, 0);
/* plot -> AUTOSCALE X */
	gui_direct_check("Auto X range",&plot->auto_x,auto_x_toggle,NULL,vbox);
/* initial state */
	if(plot->auto_x){
		gtk_widget_set_sensitive(plot_xmin,FALSE);
	        gtk_widget_set_sensitive(plot_xmax,FALSE);
		set_auto_x();
	} else {
	        gtk_widget_set_sensitive(plot_xmin,TRUE);
	        gtk_widget_set_sensitive(plot_xmax,TRUE);
	}
/* plot -> YMIN */
	hbox = gtk_hbox_new(FALSE, 0);
	gtk_box_pack_start(GTK_BOX(vbox), hbox, FALSE, FALSE, PANEL_SPACING);
	label = gtk_label_new(g_strdup_printf("Y min"));
	gtk_box_pack_start(GTK_BOX(hbox), label, FALSE, FALSE, 0);
	plot_ymin = gtk_entry_new();
	gtk_entry_set_text(GTK_ENTRY(plot_ymin),g_strdup_printf("%-9.6f",plot->ymin));
	gtk_box_pack_end(GTK_BOX(hbox),plot_ymin, FALSE, FALSE, 0);
/* plot -> YMAX */
	hbox = gtk_hbox_new(FALSE, 0);
	gtk_box_pack_start(GTK_BOX(vbox), hbox, FALSE, FALSE, PANEL_SPACING);
	label = gtk_label_new(g_strdup_printf("Y max"));
	gtk_box_pack_start(GTK_BOX(hbox), label, FALSE, FALSE, 0);
	plot_ymax = gtk_entry_new();
	gtk_entry_set_text(GTK_ENTRY(plot_ymax),g_strdup_printf("%-9.6f",plot->ymax));
	gtk_box_pack_end(GTK_BOX(hbox),plot_ymax, FALSE, FALSE, 0);
/* plot -> AUTOSCALE Y */
	gui_direct_check("Auto Y range",&plot->auto_y,auto_y_toggle,NULL,vbox);
/* initial state */
	if(plot->auto_y){
	        gtk_widget_set_sensitive(plot_ymin,FALSE);
	        gtk_widget_set_sensitive(plot_ymax,FALSE);
		set_auto_y();
	} else {
	        gtk_widget_set_sensitive(plot_ymin,TRUE);
	        gtk_widget_set_sensitive(plot_ymax,TRUE);
	}
/* plot -> tics */
	gui_direct_spin("X tics",&plot->xtics,1.0,50.0,1.0,NULL,NULL,vbox);
	gui_direct_spin("Y tics",&plot->ytics,1.0,50.0,1.0,NULL,NULL,vbox);
/* Action buttons */
	gui_stock_button(GTK_STOCK_EXECUTE, exec_plot, NULL, GTK_DIALOG(window)->action_area);
	gui_stock_button(GTK_STOCK_CLOSE, plot_quit, dialog, GTK_DIALOG(window)->action_area);
/* connect to signals */
        g_signal_connect(GTK_NOTEBOOK(notebook),"switch-page",GTK_SIGNAL_FUNC(plot_page_switch),NULL);
	/*select in page 1*/
	if(plot->plot_mask==PLOT_NONE) plot->type=PLOT_NONE;
	else if(plot->plot_mask&PLOT_ENERGY) plot_energy();
	else if(plot->plot_mask&PLOT_FORCE) plot_force();
	else if(plot->plot_mask&PLOT_VOLUME) plot_volume();
	else if(plot->plot_mask&PLOT_PRESSURE) plot_pressure();
/* That's all folks */
//        sysenv.refresh_dialog=TRUE;
if((plot->task)!=NULL) {
	while(plot->task->progress<50.0) usleep(50*1000);/*sleep 50ms until half job is done*/
}

	gtk_widget_show_all(window);
	sysenv.refresh_dialog=TRUE;
	data->locked=FALSE;
	data->cur_frame=cur_frame;

        tree_model_refresh(data);
        model_prep(data); 
/*restore camera*/
        data->camera = camera;
        g_free(data->camera_default);
        data->camera_default = camera;
/*update everything*/
	redraw_canvas(ALL);
}
/**************************************/
/* <- FROM HERE GRAPH_CONTROLS GUI -> */
/**************************************/
/**********************/
/* setup the controls */
/**********************/
void sync_graph_controls(struct graph_pak *graph){
	gchar *text;
	GSList *list;
	g_data_y *p_y;
	g_assert(graph != NULL);
	/*text*/
	if(graph->title!=NULL) text=g_strdup(graph->title);
	else text=g_strdup("");
	GUI_ENTRY_TEXT(GRAPH_UI.title,text);g_free(text);
	if(graph->sub_title!=NULL) text=g_strdup(graph->sub_title);
	else text=g_strdup("");
	GUI_ENTRY_TEXT(GRAPH_UI.sub_title,text);g_free(text);
	if(graph->x_title!=NULL) text=g_strdup(graph->x_title);
	else text=g_strdup("");
	GUI_ENTRY_TEXT(GRAPH_UI.x_title,text);g_free(text);
	if(graph->y_title!=NULL) text=g_strdup(graph->y_title);
	else text=g_strdup("");
	GUI_ENTRY_TEXT(GRAPH_UI.y_title,text);g_free(text);
	/*values*/
	text=g_strdup_printf("%G",graph->xmin);
	GUI_ENTRY_TEXT(GRAPH_UI.xmin,text);g_free(text);
	text=g_strdup_printf("%G",graph->xmax);
	GUI_ENTRY_TEXT(GRAPH_UI.xmax,text);g_free(text);
	text=g_strdup_printf("%G",graph->ymin);
	GUI_ENTRY_TEXT(GRAPH_UI.ymin,text);g_free(text);
	text=g_strdup_printf("%G",graph->ymax);
	GUI_ENTRY_TEXT(GRAPH_UI.ymax,text);g_free(text);
	/*ticks & auto*/
	GRAPH_UI.auto_x=TRUE;
	GRAPH_UI.auto_y=TRUE;
	GRAPH_UI.xticks=(gdouble)graph->xticks;
	GUI_SPIN_SET(GRAPH_UI.xtics,GRAPH_UI.xticks);
	GRAPH_UI.yticks=(gdouble)graph->yticks;
	GUI_SPIN_SET(GRAPH_UI.ytics,GRAPH_UI.yticks);
	/*sets, which depends on graph->type*/
	GRAPH_UI.set_number=1.;
	switch(graph->type){
		case GRAPH_IX_TYPE:
		case GRAPH_IY_TYPE:
		case GRAPH_XX_TYPE:
		case GRAPH_XY_TYPE:
			/*we have multiple sets*/
			GUI_SPIN_RANGE(GRAPH_UI.set,1.,(gdouble)graph->size);
			list=graph->set_list;/*this is g_data_x*/
			list=g_slist_next(list);/*this is g_data_y*/
			p_y = (g_data_y *) list->data;
			if(p_y->y_size==0){
				/*this is normal for non META/MINHOP*/
				list=g_slist_next(list);
				p_y = (g_data_y *) list->data;
			}
			GUI_UNLOCK(GRAPH_UI.type);
			switch(p_y->type){
				case GRAPH_IY_TYPE:
				case GRAPH_XY_TYPE:
					GUI_COMBOBOX_SET(GRAPH_UI.type,2);
					break;
				case GRAPH_IX_TYPE:
				case GRAPH_XX_TYPE:
				default:
					GUI_COMBOBOX_SET(GRAPH_UI.type,1);
			}
			GUI_LOCK(GRAPH_UI.type);
			GUI_UNLOCK(GRAPH_UI.size);
			text=g_strdup_printf("%i",p_y->y_size);
			GUI_ENTRY_TEXT(GRAPH_UI.size,text);g_free(text);
			GUI_LOCK(GRAPH_UI.size);
			switch(p_y->line){
				case GRAPH_LINE_SINGLE:
					GUI_COMBOBOX_SET(GRAPH_UI.line,1);break;
				case GRAPH_LINE_DASH:
					GUI_COMBOBOX_SET(GRAPH_UI.line,2);break;
				case GRAPH_LINE_DOT:
					GUI_COMBOBOX_SET(GRAPH_UI.line,3);break;
				case GRAPH_LINE_THICK:
					GUI_COMBOBOX_SET(GRAPH_UI.line,4);break;
				case GRAPH_LINE_NONE:
				default:
					GUI_COMBOBOX_SET(GRAPH_UI.line,0);
			}
			switch(p_y->color){
				case GRAPH_COLOR_BLACK:
					GUI_COMBOBOX_SET(GRAPH_UI.color,0);break;
				case GRAPH_COLOR_WHITE:
					GUI_COMBOBOX_SET(GRAPH_UI.color,1);break;
				case GRAPH_COLOR_BLUE:
					GUI_COMBOBOX_SET(GRAPH_UI.color,2);break;
				case GRAPH_COLOR_GREEN:
					GUI_COMBOBOX_SET(GRAPH_UI.color,3);break;
				case GRAPH_COLOR_RED:
					GUI_COMBOBOX_SET(GRAPH_UI.color,4);break;
				case GRAPH_COLOR_YELLOW:
					GUI_COMBOBOX_SET(GRAPH_UI.color,5);break;
				case GRAPH_COLOR_GRAY:
					GUI_COMBOBOX_SET(GRAPH_UI.color,6);break;
				case GRAPH_COLOR_NAVY:
					GUI_COMBOBOX_SET(GRAPH_UI.color,7);break;
				case GRAPH_COLOR_LIME:
					GUI_COMBOBOX_SET(GRAPH_UI.color,8);break;
				case GRAPH_COLOR_TEAL:
					GUI_COMBOBOX_SET(GRAPH_UI.color,9);break;
				case GRAPH_COLOR_AQUA:
					GUI_COMBOBOX_SET(GRAPH_UI.color,10);break;
				case GRAPH_COLOR_MAROON:
					GUI_COMBOBOX_SET(GRAPH_UI.color,11);break;
				case GRAPH_COLOR_PURPLE:
					GUI_COMBOBOX_SET(GRAPH_UI.color,12);break;
				case GRAPH_COLOR_OLIVE:
					GUI_COMBOBOX_SET(GRAPH_UI.color,13);break;
				case GRAPH_COLOR_SILVER:
					GUI_COMBOBOX_SET(GRAPH_UI.color,14);break;
				case GRAPH_COLOR_FUSHIA:
					GUI_COMBOBOX_SET(GRAPH_UI.color,15);break;
				case GRAPH_COLOR_DEFAULT:
				default:
					GUI_COMBOBOX_SET(GRAPH_UI.color,16);
			}
			/*TODO: act on selection?*/
			switch(p_y->symbol[0]){
				case GRAPH_SYMB_CROSS:
					GUI_COMBOBOX_SET(GRAPH_UI.symbol,1);break;
				case GRAPH_SYMB_SQUARE:
					GUI_COMBOBOX_SET(GRAPH_UI.symbol,2);break;
				case GRAPH_SYMB_TRI_DN:
					GUI_COMBOBOX_SET(GRAPH_UI.symbol,3);break;
				case GRAPH_SYMB_TRI_UP:
					GUI_COMBOBOX_SET(GRAPH_UI.symbol,4);break;
				case GRAPH_SYMB_DIAM:
					GUI_COMBOBOX_SET(GRAPH_UI.symbol,5);break;
				case GRAPH_SYMB_NONE:
				default:
					GUI_COMBOBOX_SET(GRAPH_UI.symbol,0);
			}
			break;
		case GRAPH_FREQUENCY:
		case GRAPH_BAND:
		case GRAPH_DOS:
		case GRAPH_BANDOS:
			/*these types will be soon replace by above GRAPH_??_TYPE*/
			GUI_SPIN_RANGE(GRAPH_UI.set,1.,1.);
			GUI_UNLOCK(GRAPH_UI.type);
			GUI_COMBOBOX_SET(GRAPH_UI.type,0);
			GUI_LOCK(GRAPH_UI.type);
			GUI_UNLOCK(GRAPH_UI.size);
			text=g_strdup("1");
			GUI_ENTRY_TEXT(GRAPH_UI.size,text);g_free(text);
			GUI_LOCK(GRAPH_UI.size);
			GUI_COMBOBOX_SET(GRAPH_UI.symbol,0);
			GUI_COMBOBOX_SET(GRAPH_UI.line,0);
			GUI_COMBOBOX_SET(GRAPH_UI.color,16);
			break;
		case GRAPH_REGULAR:
		default:
			/*assume 1 set x[] only*/
			GUI_SPIN_RANGE(GRAPH_UI.set,1.,1.);
			GUI_UNLOCK(GRAPH_UI.type);
			GUI_COMBOBOX_SET(GRAPH_UI.type,0);
			GUI_LOCK(GRAPH_UI.type);
			GUI_UNLOCK(GRAPH_UI.size);
			text=g_strdup("1");
			GUI_ENTRY_TEXT(GRAPH_UI.size,text);g_free(text);
			GUI_LOCK(GRAPH_UI.size);
			GUI_COMBOBOX_SET(GRAPH_UI.symbol,0);
			GUI_COMBOBOX_SET(GRAPH_UI.line,0);
			GUI_COMBOBOX_SET(GRAPH_UI.color,16);
	}
	GUI_SPIN_SET(GRAPH_UI.set,GRAPH_UI.set_number);
}
/****************/
/* update xtics */
/****************/
void spin_update_xtics(void){



}
/****************/
/* update ytics */
/****************/
void spin_update_ytics(void){



}
/*****************/
/* toggle auto_x */
/*****************/
void toggle_auto_x(){
	struct model_pak *data;
	struct graph_pak *graph=(struct graph_pak *)GRAPH_UI.graph_active;
	gint idx;
	gchar *text;
	GSList *list;
	g_data_x *p_x;
	/**/
	g_assert(graph != NULL);
	if(GRAPH_UI.auto_x){
		data = sysenv.active_model;
		if(data){
			if(data->graph_active){
				if(graph!=data->graph_active){
					/*graph data changed*/
					sync_graph_controls((struct graph_pak *)data->graph_active);
					GRAPH_UI.graph_active=data->graph_active;
					GUI_LOCK(GRAPH_UI.xmin);
					GUI_LOCK(GRAPH_UI.xmax);
					return;
				}
			}
		}
		/*recalculate xmin,xmax*/
		switch(graph->type){
			case GRAPH_IY_TYPE:
			case GRAPH_XY_TYPE:
			case GRAPH_IX_TYPE:
			case GRAPH_XX_TYPE:
				list=graph->set_list;
				if(!list) return;/*BAD: no x data!*/
				p_x=(g_data_x *) list->data;
				graph->xmin=p_x->x[0];
				graph->xmax=p_x->x[0];
				for(idx=0;idx<p_x->x_size;idx++){
					if(p_x->x[idx]<graph->xmin) graph->xmin=p_x->x[idx];
					if(p_x->x[idx]>graph->xmax) graph->xmax=p_x->x[idx];
				}
				text=g_strdup_printf("%G",graph->xmin);
				GUI_ENTRY_TEXT(GRAPH_UI.xmin,text);g_free(text);
				text=g_strdup_printf("%G",graph->xmax);
				GUI_ENTRY_TEXT(GRAPH_UI.xmax,text);g_free(text);
				break;
			case GRAPH_FREQUENCY:
				if(data->have_frequency==TRUE){
					graph->xmin=0.;
					graph->xmax=data->freq[0];
					for(idx=0;idx<data->nfreq;idx++) {
						if(data->freq[idx]>graph->xmax) graph->xmax=data->freq[idx];
					}
					idx=((gint)(graph->xmax/1000)+1)*1000.0;/*magic formula*/
					graph->xmax=(gdouble)idx;
					text=g_strdup_printf("%G",graph->xmin);
					GUI_ENTRY_TEXT(GRAPH_UI.xmin,text);g_free(text);
					text=g_strdup_printf("%G",graph->xmax);
					GUI_ENTRY_TEXT(GRAPH_UI.xmax,text);g_free(text);
				}/*else... there is a problem*/
				break;
			case GRAPH_BAND:
				if(data->kpts_d){
					graph->xmin=data->kpts_d[0];
					graph->xmax=data->kpts_d[0];
					for(idx=0;idx<(data->nkpoints);idx++){
						if(data->kpts_d[idx]<graph->xmin) graph->xmin=data->kpts_d[idx];
						if(data->kpts_d[idx]>graph->xmax) graph->xmax=data->kpts_d[idx];
					}
					text=g_strdup_printf("%G",graph->xmin);
					GUI_ENTRY_TEXT(GRAPH_UI.xmin,text);g_free(text);
					text=g_strdup_printf("%G",graph->xmax);
					GUI_ENTRY_TEXT(GRAPH_UI.xmax,text);g_free(text);
				}/*else... there is a problem*/
				break;
			case GRAPH_DOS:
				if(data->dos_eval){
					gdouble x;
					/*everything is scaled to Fermi lavel*/
					graph->xmin=data->dos_eval[0]-data->efermi;
					graph->xmax=data->dos_eval[0]-data->efermi;
					for(idx=0;idx<data->ndos;idx++){
						x=data->dos_eval[idx]-data->efermi;
						if(x<graph->xmin) graph->xmin=x;
						if(x>graph->xmax) graph->xmax=x;
					}
					text=g_strdup_printf("%G",graph->xmin);
					GUI_ENTRY_TEXT(GRAPH_UI.xmin,text);g_free(text);
					text=g_strdup_printf("%G",graph->xmax);
					GUI_ENTRY_TEXT(GRAPH_UI.xmax,text);g_free(text);
				}/*else... there is a problem*/
				break;
			case GRAPH_BANDOS:
				if(data->dos_spin_up){
					gdouble x;
					x=data->dos_spin_up[0];
					graph->xmin=x;
					graph->xmax=x;
					for(idx=0;idx<data->ndos;idx++){
						x=data->dos_spin_up[idx];
						if(data->spin_polarized) x+=data->dos_spin_down[idx];
						if(x<graph->xmin) graph->xmin=x;
						if(x>graph->xmax) graph->xmax=x;
					}
					text=g_strdup_printf("%G",graph->xmin);
					GUI_ENTRY_TEXT(GRAPH_UI.xmin,text);g_free(text);
					text=g_strdup_printf("%G",graph->xmax);
					GUI_ENTRY_TEXT(GRAPH_UI.xmax,text);g_free(text);
				}/*else... there is a problem*/
				break;
			case GRAPH_REGULAR:
			default:
				/*there is one y set, xmin and xmax are determined only by graph->size*/
				graph->xmin=0.;
				graph->xmax=(gdouble)graph->size;
		}
		GUI_LOCK(GRAPH_UI.xmin);
		GUI_LOCK(GRAPH_UI.xmax);
	}else{
		/*don't recalculate xmin,xmax*/
		GUI_UNLOCK(GRAPH_UI.xmin);
		GUI_UNLOCK(GRAPH_UI.xmax);
	}
}
/*****************/
/* toggle auto_y */
/*****************/
void toggle_auto_y(){
	struct model_pak *data;
	struct graph_pak *graph=(struct graph_pak *)GRAPH_UI.graph_active;
	gint idx;
	gchar *text;
	GSList *list;
	g_data_y *p_y;
	/**/
	g_assert(graph != NULL);
	if(GRAPH_UI.auto_y){
		data = sysenv.active_model;
		if(data){
			if(data->graph_active){
				if(graph!=data->graph_active){
					/*graph data changed*/
					sync_graph_controls((struct graph_pak *)data->graph_active);
					GRAPH_UI.graph_active=data->graph_active;
					GUI_LOCK(GRAPH_UI.xmin);
					GUI_LOCK(GRAPH_UI.xmax);
					return;
				}
			}
		}
		switch(graph->type){
			case GRAPH_IY_TYPE:
			case GRAPH_XY_TYPE:
			case GRAPH_IX_TYPE:
			case GRAPH_XX_TYPE:
				list=graph->set_list;/*this is g_data_x*/
				list=g_slist_next(list);/*this is g_data_y*/
				if(!list) return;/*no Y data?*/
				p_y = (g_data_y *) list->data;
				if(p_y->y_size==0) {
					/*this is normal for non META/MINHOP*/
					list=g_slist_next(list);
					p_y = (g_data_y *) list->data;
				}
				if(p_y->y[0]==-999.9){
					graph->ymin=p_y->y[1];
					graph->ymax=p_y->y[1];
				}else{
					graph->ymin=p_y->y[0];
					graph->ymax=p_y->y[0];
				}
				/*go through all y sets*/
				for(;list;list=g_slist_next(list)){
					p_y = (g_data_y *) list->data;
					if(p_y->y[0]==-999.9){
						for(idx=1;idx<p_y->y_size;idx++){
							if(p_y->y[idx]<graph->ymin) graph->ymin=p_y->y[idx];
							if(p_y->y[idx]>graph->ymax) graph->ymax=p_y->y[idx];
						}
					}else{
						for(idx=0;idx<p_y->y_size;idx++){
							if(p_y->y[idx]<graph->ymin) graph->ymin=p_y->y[idx];
							if(p_y->y[idx]>graph->ymax) graph->ymax=p_y->y[idx];
						}
					}
				}
				/*add 5%*/
				graph->ymin=graph->ymin-(graph->ymax-graph->ymin)*0.05;
				graph->ymax=graph->ymax+(graph->ymax-graph->ymin)*0.05;
				text=g_strdup_printf("%G",graph->ymin);
				GUI_ENTRY_TEXT(GRAPH_UI.ymin,text);g_free(text);
				text=g_strdup_printf("%G",graph->ymax);
				GUI_ENTRY_TEXT(GRAPH_UI.ymax,text);g_free(text);
				break;
			case GRAPH_FREQUENCY:
				/*easy*/
				graph->ymin=0.;
				graph->ymax=1.;
				text=g_strdup_printf("%G",graph->ymin);
				GUI_ENTRY_TEXT(GRAPH_UI.ymin,text);g_free(text);
				text=g_strdup_printf("%G",graph->ymax);
				GUI_ENTRY_TEXT(GRAPH_UI.ymax,text);g_free(text);
				break;
			case GRAPH_BANDOS:
				/* PASS THROUGH */
			case GRAPH_BAND:
				if(data->band_up){
					gdouble y;
					y=data->band_up[0]-data->efermi;
					graph->ymin=y;
					graph->ymax=y;
					for(idx=1;idx<(data->nbands*data->nkpoints);idx++){
						y=data->band_up[idx]-data->efermi;
						if(y<graph->ymin) graph->ymin=y;
						if(y>graph->ymax) graph->ymax=y;
					}
					if(data->band_down){
						for(idx=0;idx<(data->nbands*data->nkpoints);idx++){
							y=data->band_down[idx]-data->efermi;
							if(y<graph->ymin) graph->ymin=y;
							if(y>graph->ymax) graph->ymax=y;
						}
					}
					text=g_strdup_printf("%G",graph->ymin);
					GUI_ENTRY_TEXT(GRAPH_UI.ymin,text);g_free(text);
					text=g_strdup_printf("%G",graph->ymax);
					GUI_ENTRY_TEXT(GRAPH_UI.ymax,text);g_free(text);
				}/*else... there is a problem*/
				break;
			case GRAPH_DOS:
				if(data->dos_spin_up){
					gdouble y;
					y=data->dos_spin_up[0];
					graph->ymin=y;
					graph->ymax=y;
					for(idx=0;idx<data->ndos;idx++){
						y=data->dos_spin_up[idx];
						if(data->spin_polarized) y+=data->dos_spin_down[idx];
						if(y<graph->ymin) graph->ymin=y;
						if(y>graph->ymax) graph->ymax=y;
					}
					text=g_strdup_printf("%G",graph->ymin);
					GUI_ENTRY_TEXT(GRAPH_UI.ymin,text);g_free(text);
					text=g_strdup_printf("%G",graph->ymax);
					GUI_ENTRY_TEXT(GRAPH_UI.ymax,text);g_free(text);
				}/*else... there is a problem*/
				break;
			case GRAPH_REGULAR:
			default:
				graph->ymin=0.;
		}
		/*recalculate ymin,ymax*/
		GUI_LOCK(GRAPH_UI.ymin);
		GUI_LOCK(GRAPH_UI.ymax);
	}else{
		GUI_UNLOCK(GRAPH_UI.ymin);
		GUI_UNLOCK(GRAPH_UI.ymax);
	}
}
/*********************/
/* update set number */
/*********************/
void spin_update_set(void){
	struct model_pak *data;
	struct graph_pak *graph=(struct graph_pak *)GRAPH_UI.graph_active;
	gchar *text;
	GSList *list;
	g_data_y *p_y;
	graph_symbol symb;
	gint idx=(gint)GRAPH_UI.set_number;
	/**/
/* --- detect if graph_active has change and react accordingly*/
	g_assert(graph != NULL);
	data = sysenv.active_model;
	if(data){
		if(data->graph_active){
			if(graph!=data->graph_active){
				/*graph data changed*/
				sync_graph_controls((struct graph_pak *)data->graph_active);
				GRAPH_UI.graph_active=data->graph_active;
				return;/*does not make sense to modify set number anymore*/
			}
		}/*else user has moved, but graph should still be valid*/
	}/*no active model doesn't mean that graph data is invalid*/
	if(idx>graph->size){
		/*we do have a problem with sets...*/
		if((graph->type==GRAPH_IX_TYPE)||(graph->type==GRAPH_IY_TYPE)||(graph->type==GRAPH_XX_TYPE)||(graph->type==GRAPH_XY_TYPE)){
			GUI_SPIN_RANGE(GRAPH_UI.set,1.,(gdouble)graph->size);
			GRAPH_UI.set_number=(gdouble)graph->size;
		}else{
			GUI_SPIN_RANGE(GRAPH_UI.set,1.,1.);
			GRAPH_UI.set_number=1.;
		}
		GUI_SPIN_SET(GRAPH_UI.set,GRAPH_UI.set_number);
		return;
	}
        switch(graph->type){
                case GRAPH_IX_TYPE:
                case GRAPH_IY_TYPE:
                case GRAPH_XX_TYPE:
                case GRAPH_XY_TYPE:
			/*select i-th set*/
			list=graph->set_list;/*this is g_data_x*/
			idx=0;
			while((list)&&(idx<(gint)GRAPH_UI.set_number)){
				list=g_slist_next(list);
				idx++;
			}
			if(!list) {
				/*Bad... less data than predicted?*/
				GUI_SPIN_RANGE(GRAPH_UI.set,1.,(gdouble)graph->size);
				GRAPH_UI.set_number=1.;/*return to a safe one*/
				GUI_SPIN_SET(GRAPH_UI.set,GRAPH_UI.set_number);
				return;
			}
			p_y = (g_data_y *) list->data;
			/*update type, this is important because graph type can change*/
			GUI_UNLOCK(GRAPH_UI.type);
			switch(p_y->type){
				case GRAPH_IY_TYPE:
				case GRAPH_XY_TYPE:
					GUI_COMBOBOX_SET(GRAPH_UI.type,2);
					break;
				case GRAPH_IX_TYPE:
				case GRAPH_XX_TYPE:
				default:
					GUI_COMBOBOX_SET(GRAPH_UI.type,1);
			}
			GUI_LOCK(GRAPH_UI.type);
			/*the size of each set can also differ (especially for IX,XX types)*/
			GUI_UNLOCK(GRAPH_UI.size);
			text=g_strdup_printf("%i",p_y->y_size);
			GUI_ENTRY_TEXT(GRAPH_UI.size,text);g_free(text);
			GUI_LOCK(GRAPH_UI.size);
			/*other are different from one set to another*/
			switch(p_y->line){
				case GRAPH_LINE_SINGLE:
					GUI_COMBOBOX_SET(GRAPH_UI.line,1);break;
				case GRAPH_LINE_DASH:
					GUI_COMBOBOX_SET(GRAPH_UI.line,2);break;
				case GRAPH_LINE_DOT:
					GUI_COMBOBOX_SET(GRAPH_UI.line,3);break;
				case GRAPH_LINE_THICK:
					GUI_COMBOBOX_SET(GRAPH_UI.line,4);break;
				case GRAPH_LINE_NONE:
				default:
					GUI_COMBOBOX_SET(GRAPH_UI.line,0);
			}
			switch(p_y->color){
				case GRAPH_COLOR_BLACK:
					GUI_COMBOBOX_SET(GRAPH_UI.color,0);break;
				case GRAPH_COLOR_WHITE:
					GUI_COMBOBOX_SET(GRAPH_UI.color,1);break;
				case GRAPH_COLOR_BLUE:
					GUI_COMBOBOX_SET(GRAPH_UI.color,2);break;
				case GRAPH_COLOR_GREEN:
					GUI_COMBOBOX_SET(GRAPH_UI.color,3);break;
				case GRAPH_COLOR_RED:
					GUI_COMBOBOX_SET(GRAPH_UI.color,4);break;
				case GRAPH_COLOR_YELLOW:    
					GUI_COMBOBOX_SET(GRAPH_UI.color,5);break;
				case GRAPH_COLOR_GRAY:
					GUI_COMBOBOX_SET(GRAPH_UI.color,6);break;
				case GRAPH_COLOR_NAVY:
					GUI_COMBOBOX_SET(GRAPH_UI.color,7);break;
				case GRAPH_COLOR_LIME:
					GUI_COMBOBOX_SET(GRAPH_UI.color,8);break;
				case GRAPH_COLOR_TEAL:
					GUI_COMBOBOX_SET(GRAPH_UI.color,9);break;
				case GRAPH_COLOR_AQUA:
					GUI_COMBOBOX_SET(GRAPH_UI.color,10);break;
				case GRAPH_COLOR_MAROON:
					GUI_COMBOBOX_SET(GRAPH_UI.color,11);break;
				case GRAPH_COLOR_PURPLE:
					GUI_COMBOBOX_SET(GRAPH_UI.color,12);break;
				case GRAPH_COLOR_OLIVE:
					GUI_COMBOBOX_SET(GRAPH_UI.color,13);break;
				case GRAPH_COLOR_SILVER:
					GUI_COMBOBOX_SET(GRAPH_UI.color,14);break;
				case GRAPH_COLOR_FUSHIA:
					GUI_COMBOBOX_SET(GRAPH_UI.color,15);break;
				case GRAPH_COLOR_DEFAULT:
				default:
					GUI_COMBOBOX_SET(GRAPH_UI.color,16);
			}
			/*TODO: act on selection?*/
			if(p_y->symbol!=NULL) symb=p_y->symbol[0];
			else symb=GRAPH_SYMB_NONE;/*actually a set can be empty AND not the last one*/
			switch(symb){
				case GRAPH_SYMB_CROSS:
					GUI_COMBOBOX_SET(GRAPH_UI.symbol,1);break;
				case GRAPH_SYMB_SQUARE:
					GUI_COMBOBOX_SET(GRAPH_UI.symbol,2);break;
				case GRAPH_SYMB_TRI_DN:
					GUI_COMBOBOX_SET(GRAPH_UI.symbol,3);break;
				case GRAPH_SYMB_TRI_UP:
					GUI_COMBOBOX_SET(GRAPH_UI.symbol,4);break;
				case GRAPH_SYMB_DIAM:
					GUI_COMBOBOX_SET(GRAPH_UI.symbol,5);break;
				case GRAPH_SYMB_NONE:
				default:
					GUI_COMBOBOX_SET(GRAPH_UI.symbol,0);
			}
			break;
                case GRAPH_FREQUENCY:
                case GRAPH_BAND:
                case GRAPH_DOS:
                case GRAPH_BANDOS:
                        /*these types will be soon replace by above GRAPH_??_TYPE*/
			GUI_UNLOCK(GRAPH_UI.type);
			GUI_COMBOBOX_SET(GRAPH_UI.type,0);
			GUI_LOCK(GRAPH_UI.type);
			GUI_UNLOCK(GRAPH_UI.size);
			text=g_strdup("1");
			GUI_ENTRY_TEXT(GRAPH_UI.size,text);g_free(text);
			GUI_COMBOBOX_SET(GRAPH_UI.symbol,0);
			GUI_COMBOBOX_SET(GRAPH_UI.line,0);
			GUI_COMBOBOX_SET(GRAPH_UI.color,16);
                        break;
                case GRAPH_REGULAR:
                default:
			GUI_UNLOCK(GRAPH_UI.type);
			GUI_COMBOBOX_SET(GRAPH_UI.type,0);
			GUI_LOCK(GRAPH_UI.type);
			GUI_UNLOCK(GRAPH_UI.size);
			text=g_strdup("1");
			GUI_ENTRY_TEXT(GRAPH_UI.size,text);g_free(text);
			GUI_COMBOBOX_SET(GRAPH_UI.symbol,0);
			GUI_COMBOBOX_SET(GRAPH_UI.line,0);
			GUI_COMBOBOX_SET(GRAPH_UI.color,16);
			break;
	}

}
/**************************/
/* apply graph dimensions */
/**************************/
void graph_control_apply_dim(void){
        struct model_pak *data;
        struct graph_pak *graph=(struct graph_pak *)GRAPH_UI.graph_active;
        GSList *list;
	gint idx;
        g_data_y *p_y;
	gboolean have_changed=FALSE;
	/**/
/* --- detect if graph_active has changed*/
	g_assert(graph != NULL);
	data = sysenv.active_model;
	if(data){
		if(data->graph_active)
			if(graph!=data->graph_active) have_changed=TRUE;/*graph has changed*/
	}
/*Apply global, even if graph has changed*/
	if(graph->title!=NULL) g_free(graph->title);
	if(GUI_ENTRY_LENGTH(GRAPH_UI.title)>0) GUI_ENTRY_GET_TEXT(GRAPH_UI.title,graph->title);
	else graph->title=NULL;
	if(graph->sub_title!=NULL) g_free(graph->sub_title);
	if(GUI_ENTRY_LENGTH(GRAPH_UI.sub_title)>0) GUI_ENTRY_GET_TEXT(GRAPH_UI.sub_title,graph->sub_title);
	else graph->sub_title=NULL;
	if(graph->x_title!=NULL) g_free(graph->x_title);
	if(GUI_ENTRY_LENGTH(GRAPH_UI.x_title)>0) GUI_ENTRY_GET_TEXT(GRAPH_UI.x_title,graph->x_title);
	else graph->x_title=NULL;
	if(graph->y_title!=NULL) g_free(graph->y_title);
	if(GUI_ENTRY_LENGTH(GRAPH_UI.y_title)>0) GUI_ENTRY_GET_TEXT(GRAPH_UI.y_title,graph->y_title);
	else graph->y_title=NULL;
	GUI_REG_VAL(GRAPH_UI.xmin,graph->xmin,"%lf");
	GUI_REG_VAL(GRAPH_UI.xmax,graph->xmax,"%lf");
	GUI_REG_VAL(GRAPH_UI.ymin,graph->ymin,"%lf");
	GUI_REG_VAL(GRAPH_UI.ymax,graph->ymax,"%lf");
	graph->xticks=(gint)GRAPH_UI.xticks;
	graph->yticks=(gint)GRAPH_UI.yticks;
/*Apply Set, only if have_changed==FALSE*/
if(have_changed==FALSE){
	if((graph->type==GRAPH_IX_TYPE)||(graph->type==GRAPH_IY_TYPE)||(graph->type==GRAPH_XX_TYPE)||(graph->type==GRAPH_XY_TYPE)){
		/*get set number first*/
		list=graph->set_list;/*this is g_data_x*/
		idx=0;
		while((list)&&(idx<(gint)GRAPH_UI.set_number)){
			list=g_slist_next(list);
			idx++;
		}
		if(!list){
			/*BAD, data set are missing...*/
			return;
		}
		p_y = (g_data_y *) list->data;
		/*we can't change type or size*/
		GUI_COMBOBOX_GET(GRAPH_UI.line,idx);
		switch(idx){
			case 1:
				p_y->line=GRAPH_LINE_SINGLE;break;
			case 2:
				p_y->line=GRAPH_LINE_DASH;break;
			case 3:
				p_y->line=GRAPH_LINE_DOT;break;
			case 4:
				p_y->line=GRAPH_LINE_THICK;break;
			case 0:
			default:
				p_y->line=GRAPH_LINE_NONE;
		}
		GUI_COMBOBOX_GET(GRAPH_UI.color,idx);
		switch(idx){
			case 0:
				p_y->color=GRAPH_COLOR_BLACK;break;
			case 1:
				p_y->color=GRAPH_COLOR_WHITE;break;
			case 2:
				p_y->color=GRAPH_COLOR_BLUE;break;
			case 3:
				p_y->color=GRAPH_COLOR_GREEN;break;
			case 4:
				p_y->color=GRAPH_COLOR_RED;break;
			case 5:
				p_y->color=GRAPH_COLOR_YELLOW;break;
			case 6:
				p_y->color=GRAPH_COLOR_GRAY;break;
			case 7:
				p_y->color=GRAPH_COLOR_NAVY;break;
			case 8:
				p_y->color=GRAPH_COLOR_LIME;break;
			case 9:
				p_y->color=GRAPH_COLOR_TEAL;break;
			case 10:
				p_y->color=GRAPH_COLOR_AQUA;break;
			case 11:
				p_y->color=GRAPH_COLOR_MAROON;break;
			case 12:
				p_y->color=GRAPH_COLOR_PURPLE;break;
			case 13:
				p_y->color=GRAPH_COLOR_OLIVE;break;
			case 14:
				p_y->color=GRAPH_COLOR_SILVER;break;
			case 15:
				p_y->color=GRAPH_COLOR_FUSHIA;break;
			case 16:
			default:
				p_y->color=GRAPH_COLOR_DEFAULT;
		}
		if(p_y->symbol!=NULL){
			gint i;
			GUI_COMBOBOX_GET(GRAPH_UI.symbol,idx);
			switch(idx){
				case 1:
					for(i=0;i<p_y->y_size;i++) p_y->symbol[i]=GRAPH_SYMB_CROSS;
					break;
				case 2:
					for(i=0;i<p_y->y_size;i++) p_y->symbol[i]=GRAPH_SYMB_SQUARE;
					break;
				case 3:
					for(i=0;i<p_y->y_size;i++) p_y->symbol[i]=GRAPH_SYMB_TRI_DN;
					break;
				case 4:
					for(i=0;i<p_y->y_size;i++) p_y->symbol[i]=GRAPH_SYMB_TRI_UP;
					break;
				case 5:
					for(i=0;i<p_y->y_size;i++) p_y->symbol[i]=GRAPH_SYMB_DIAM;
					break;
				case 0:
				default:
					for(i=0;i<p_y->y_size;i++) p_y->symbol[i]=GRAPH_SYMB_NONE;
			}
		}


		
	}
}
/*no need to re-sync graph data*/
	tree_model_refresh(data);
	redraw_canvas(ALL);
}
/**************************/
/* reset graph dimensions */
/**************************/
void graph_control_reset_dim(void){

}
/*************************/
/* plot controls cleanup */
/*************************/
void graph_controls_cleanup(struct model_pak *model){
	g_assert(model != NULL);
	if(graph_ui!=NULL) g_free(graph_ui);
	graph_ui=NULL;
	model->graph_ui_active=NULL;
}
void quit_graph_control_gui(GUI_OBJ *w, gpointer data){
        struct model_pak *model=data;
        graph_controls_cleanup(data);
        dialog_destroy(w,model);
}
/**********************/
/* graph controls gui */
/**********************/
void gui_graph_controls(void){
        gchar *title;
        gpointer dialog;
        GUI_OBJ *frame, *vbox, *table;
        GUI_OBJ *button;
	GUI_OBJ *window;
        /* special */
        struct model_pak *data;
	struct graph_pak *graph;
	/* checks */
	data = sysenv.active_model;
	if (!data) return;
	if(!data->graph_active) return;
	if(data->graph_ui_active!=NULL) return;/*graph controls already opened TODO: bring it to front!*/
	/* init */
	graph=(struct graph_pak *)data->graph_active;
	graph_ui=g_malloc(sizeof(struct graph_ui));
	g_assert(graph_ui!=NULL);
	data->graph_ui_active=(gpointer)graph_ui;
	graph_ui->graph_active=(gpointer)graph;
	title = g_strdup("GRAPH CONTROLS");
	dialog = dialog_request(PLOTS, title, NULL, graph_controls_cleanup, data);
	g_free(title);
	window = dialog_window(dialog);
#define GRAPH_UI (*graph_ui)
	/*frame*/
	GUI_FRAME_WINDOW(window,frame);
	GUI_VBOX_FRAME(frame,vbox);
	GUI_TABLE_NOTE(vbox,table,12,3);
/* --- Title & Axis */
	GUI_LABEL_TABLE(table,"Title & Axis",0,3,0,1);
	/* line 1 */
	GUI_TEXT_TABLE(table,GRAPH_UI.title,graph->title,"MAIN TITLE:",0,3,1,2);
GUI_TOOLTIP(GRAPH_UI.title,"Title of the graph.");
	/* line 2 */
	GUI_TEXT_TABLE(table,GRAPH_UI.sub_title,graph->sub_title,"(Opt.) SUB-TITLE:",0,3,2,3);
GUI_TOOLTIP(GRAPH_UI.sub_title,"Optional sub-title of the graph.");
	/* line 3 */
	GUI_TEXT_TABLE(table,GRAPH_UI.x_title,graph->x_title,"X axis:",0,3,3,4);
GUI_TOOLTIP(GRAPH_UI.x_title,"Title of the X axis.");
	/* line 4 */
	GUI_TEXT_TABLE(table,GRAPH_UI.y_title,graph->y_title,"Y axis:",0,3,4,5);
GUI_TOOLTIP(GRAPH_UI.y_title,"Title of the Y axis.");
/* --- Dimensions */
	GUI_LABEL_TABLE(table,"Dimensions",0,3,5,6);
	/* line 6 */
	GUI_ENTRY_TABLE(table,GRAPH_UI.xmin,graph->xmin,"%.4f","X_MIN:",0,1,6,7);
GUI_TOOLTIP(GRAPH_UI.xmin,"Minimum X value of the graph.");
	GUI_ENTRY_TABLE(table,GRAPH_UI.xmax,graph->xmax,"%.4f","X_MAX:",1,2,6,7);
GUI_TOOLTIP(GRAPH_UI.xmax,"Maximum X value of the graph.");
	GUI_SPIN_TABLE(table,GRAPH_UI.xtics,GRAPH_UI.xticks,spin_update_xtics,"XTICS",2,3,6,7);
GUI_TOOLTIP(GRAPH_UI.xtics,"Number of tics on X axis.");
	/* line 7 */
	GUI_ENTRY_TABLE(table,GRAPH_UI.ymin,graph->ymin,"%.4f","Y_MIN:",0,1,7,8);
GUI_TOOLTIP(GRAPH_UI.ymin,"Minimum Y value of the graph.");
	GUI_ENTRY_TABLE(table,GRAPH_UI.ymax,graph->ymax,"%.4f","Y_MAX:",1,2,7,8);
GUI_TOOLTIP(GRAPH_UI.ymax,"Maximum Y value of the graph.");
	GUI_SPIN_TABLE(table,GRAPH_UI.ytics,GRAPH_UI.yticks,spin_update_ytics,"YTICS",2,3,7,8);
GUI_TOOLTIP(GRAPH_UI.ytics,"Number of tics on Y axis.");
	/* line 8 */
	GUI_CHECK_TABLE(table,button,GRAPH_UI.auto_x,toggle_auto_x,"Auto_X",0,1,8,9);/*not calling anything*/
GUI_TOOLTIP(button,"Set X axis limits automatically.");
GUI_TOGGLE_ON(button);
	GUI_CHECK_TABLE(table,button,GRAPH_UI.auto_y,toggle_auto_y,"Auto_Y",1,2,8,9);/*not calling anything*/
GUI_TOOLTIP(button,"Set Y axis limits automatically.");
GUI_TOGGLE_ON(button);
	GUI_2BUTTONS_TABLE(table,graph_control_apply_dim,graph_control_reset_dim,2,3,8,9);
/* --- special */
	GUI_LABEL_TABLE(table,"Special",0,3,9,10);
	/* line 10 */
	GUI_SPIN_TABLE(table,GRAPH_UI.set,GRAPH_UI.set_number,spin_update_set,"SET",0,1,10,11);
GUI_TOOLTIP(GRAPH_UI.set,"Select the graph SET number.");
	GUI_COMBOBOX_TABLE(table,GRAPH_UI.type,"TYPE:",1,2,10,11);
	GUI_COMBOBOX_ADD(GRAPH_UI.type,"NORMAL");
	GUI_COMBOBOX_ADD(GRAPH_UI.type,"X_LINE");
	GUI_COMBOBOX_ADD(GRAPH_UI.type,"Y_LINE");
GUI_TOOLTIP(GRAPH_UI.type,"Type of graph, depending on DATA:\nNORMAL -> only one y[i] set, {i} being x[i]=i;\nX_LINE -> One {Y} set correspond to all X_i values;\nY_LINE -> One {Y} set correspond to one X_i value.");
	GUI_ENTRY_TABLE(table,GRAPH_UI.size,0,"%4i","SIZE:",2,3,10,11);
GUI_TOOLTIP(GRAPH_UI.size,"Size of the current set.");
	/* line 11 */
	GUI_COMBOBOX_TABLE(table,GRAPH_UI.symbol,"SYMBOL:",0,1,11,12);
	GUI_COMBOBOX_ADD(GRAPH_UI.symbol,"NONE");
	GUI_COMBOBOX_ADD(GRAPH_UI.symbol,"CROSS");
	GUI_COMBOBOX_ADD(GRAPH_UI.symbol,"SQUARE");
	GUI_COMBOBOX_ADD(GRAPH_UI.symbol,"TRIANGLE (UP)");
	GUI_COMBOBOX_ADD(GRAPH_UI.symbol,"TRIANGLE (DN)");
	GUI_COMBOBOX_ADD(GRAPH_UI.symbol,"DIAMOND");
GUI_TOOLTIP(GRAPH_UI.symbol,"Default symbol for the graph.");
	GUI_COMBOBOX_TABLE(table,GRAPH_UI.line,"LINE:",1,2,11,12);
	GUI_COMBOBOX_ADD(GRAPH_UI.line,"NONE");
	GUI_COMBOBOX_ADD(GRAPH_UI.line,"SINGLE");
	GUI_COMBOBOX_ADD(GRAPH_UI.line,"DASH");
	GUI_COMBOBOX_ADD(GRAPH_UI.line,"DOT");
	GUI_COMBOBOX_ADD(GRAPH_UI.line,"THICK");
GUI_TOOLTIP(GRAPH_UI.line,"Default line type for the graph.");
	GUI_COMBOBOX_TABLE(table,GRAPH_UI.color,"COLOR:",2,3,11,12);
	GUI_COMBOBOX_ADD(GRAPH_UI.color,"BLACK");
	GUI_COMBOBOX_ADD(GRAPH_UI.color,"WHITE");
	GUI_COMBOBOX_ADD(GRAPH_UI.color,"BLUE");
	GUI_COMBOBOX_ADD(GRAPH_UI.color,"GREEN");
	GUI_COMBOBOX_ADD(GRAPH_UI.color,"RED");
	GUI_COMBOBOX_ADD(GRAPH_UI.color,"YELLOW");
	GUI_COMBOBOX_ADD(GRAPH_UI.color,"GRAY");
	GUI_COMBOBOX_ADD(GRAPH_UI.color,"NAVY");
	GUI_COMBOBOX_ADD(GRAPH_UI.color,"LIME");
	GUI_COMBOBOX_ADD(GRAPH_UI.color,"TEAL");
	GUI_COMBOBOX_ADD(GRAPH_UI.color,"AQUA");
	GUI_COMBOBOX_ADD(GRAPH_UI.color,"MAROON");
	GUI_COMBOBOX_ADD(GRAPH_UI.color,"PURPLE");
	GUI_COMBOBOX_ADD(GRAPH_UI.color,"OLIVE");
	GUI_COMBOBOX_ADD(GRAPH_UI.color,"SILVER");
	GUI_COMBOBOX_ADD(GRAPH_UI.color,"FUSHIA");
	GUI_COMBOBOX_ADD(GRAPH_UI.color,"DEFAULT");
GUI_TOOLTIP(GRAPH_UI.color,"Default color for the graph.");
/* initialize everything */
	sync_graph_controls(graph);
//	GUI_COMBOBOX_SETUP(xxxxx,0,xxxxx_function);
	toggle_auto_x();
	toggle_auto_y();
	GUI_LOCK(GRAPH_UI.type);
	GUI_LOCK(GRAPH_UI.size);
/* --- Outside of page */
        GUI_FRAME_WINDOW(window,frame);
        GUI_VBOX_FRAME(frame,vbox);
/* Action buttons */
        GUI_CLOSE_ACTION(window,button,quit_graph_control_gui,dialog);
/* all done */
        GUI_SHOW(window);/*display*/
        sysenv.refresh_dialog=TRUE;/*necessary?*/
}
