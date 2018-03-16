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

/* global pak structures */
extern struct sysenv_pak sysenv;
extern struct elem_pak elements[];
struct plot_pak *plot;
/* add widgets */
GtkWidget *plot_xmin;
GtkWidget *plot_xmax;
GtkWidget *plot_ymin;
GtkWidget *plot_ymax;

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
	case PLOT_BAND:
	case PLOT_DOS:
	case PLOT_BANDOS:
	case PLOT_FREQUENCY:
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
        case PLOT_BAND:
        case PLOT_DOS:
        case PLOT_BANDOS:
        case PLOT_FREQUENCY:
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
}
void plot_bandos(){
        plot->type=PLOT_BANDOS;
	plot->plot_sel=plot->plot_sel&207;/*reset 2nd page*/
	plot->plot_sel+=PLOT_BANDOS;/*selected*/
        plot->xmin=plot->band.xmin;
        plot->xmax=plot->band.xmax;
        plot->ymin=plot->dos.ymin;
        plot->ymax=plot->dos.ymax;
        plot->auto_x=TRUE;
        plot->auto_y=TRUE;
        auto_x_toggle(NULL,NULL);
        auto_y_toggle(NULL,NULL);
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
			        if(!(plot->plot_mask&(15))) plot->type=PLOT_NONE;
			        else if(plot->plot_mask&PLOT_ENERGY) plot_energy();
			        else if(plot->plot_mask&PLOT_FORCE) plot_force();
			        else if(plot->plot_mask&PLOT_VOLUME) plot_volume();
			        else if(plot->plot_mask&PLOT_PRESSURE) plot_pressure();
			}
			break;
		case 1:/*Electronic*/
			if(plot->plot_sel&PLOT_BAND) plot_band();
			else if(plot->plot_sel&PLOT_DOS) plot_dos();
			else if(plot->plot_sel&PLOT_BANDOS) plot_bandos();
			else {/*no plot selected, but plot_selection possible?*/
                                if(!(plot->plot_mask&(48))) plot->type=PLOT_NONE;
                                else if(plot->plot_mask&PLOT_BAND) plot_band();
                                else if(plot->plot_mask&PLOT_DOS) plot_dos();
                                else if(plot->plot_mask&PLOT_BANDOS) plot_bandos();
			}
			break;
		case 2:/*Frequency*/
			if(plot->plot_sel&PLOT_FREQUENCY) plot_frequency();
			else if(plot->plot_sel&PLOT_RAMAN) plot_raman();
			else {/*no plot selected, but plot_selection possible?*/
                                if(!(plot->plot_mask&(192))) plot->type=PLOT_NONE;
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
	plot->nband=0;
	plot->nfreq=0;
	plot->nraman=0;
	plot->plot_mask=PLOT_NONE;
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
			j=j+2;
			}
			plot->frequency.xmax=((gint)(plot->frequency.xmax/1000)+1)*1000.0;
		}
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
		g_free(plot->band.data);
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
	struct model_pak *data;
	g_assert(plot != NULL);
	data = sysenv.active_model;
	g_assert(data != NULL);
	/*get limit values if not auto*/
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
		task_new("PLOT-INIT",&plot_load_data,plot,&plot_prepare_data,plot,data);
	}
/* draw the graph in task */
	task_new("PLOT-DRAW",&plot_draw_graph,plot,&plot_show_graph,plot,data);
        data->locked=FALSE;
/*finalize drawing graph*/
	model_prep(data);
	model_content_refresh(data);
	gui_relation_update(data);
        data->picture_active = NULL;
        tree_model_refresh(data);
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
	GtkWidget *notebook, *page, *button, *label;
	/* special */
	GtkWidget *ionic_box, *electronic_box, *frequency_box;
	struct model_pak *data;
/* checks */
	data = sysenv.active_model;
	if (!data) return;
	/*Only one active at a time*/
	if (data->plots) return;
	data->plots=TRUE;
/* initialization */
	gui_mode_switch(FREE);
	plot=NULL;
	plot_initialize(data);
/*1- lock and try to NOT update model*/
        data->redraw=0;
        data->locked=TRUE;
	sysenv.refresh_dialog=FALSE;
	plot->model=data;	
/*2- do some work*/
	task_new("PLOT-INIT",&plot_load_data,plot,&plot_prepare_data,plot,data);
/*3- unlock model*/
	data->locked=FALSE;
	model_content_refresh(data);
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
	button = add_radio_button("Energy / step", (gpointer) plot_energy,NULL);
	if (plot->type == PLOT_ENERGY) gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(button), TRUE);
	if (!(plot->plot_mask&PLOT_ENERGY)) gtk_widget_set_sensitive(button,FALSE);
	button = add_radio_button("Forces / step", (gpointer) plot_force,NULL);
	if (plot->type == PLOT_FORCE) gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(button), TRUE);
	if (!(plot->plot_mask&PLOT_FORCE)) gtk_widget_set_sensitive(button,FALSE);
	button = add_radio_button("Volume / step", (gpointer) plot_volume,NULL);
	if (plot->type == PLOT_VOLUME) gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(button), TRUE);
	if (!(plot->plot_mask&PLOT_VOLUME)) gtk_widget_set_sensitive(button,FALSE);
	button = add_radio_button("Pressure / step", (gpointer) plot_pressure,NULL);
	if (plot->type == PLOT_PRESSURE) gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(button), TRUE);
	if (!(plot->plot_mask&PLOT_PRESSURE)) gtk_widget_set_sensitive(button,FALSE);
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
	tmp = g_strdup_printf("%9d", plot->nband);
	label = gtk_label_new(tmp);
	g_free(tmp);
	gtk_box_pack_end(GTK_BOX(hbox), label, FALSE, FALSE, 0);
/* Create a new vbox */
	vbox = gtk_vbox_new(FALSE,0);
	gtk_box_pack_start(GTK_BOX(electronic_box), vbox, TRUE, TRUE, 0);
	gtk_container_set_border_width(GTK_CONTAINER(vbox), PANEL_SPACING);
/* radio buttons */
	new_radio_group(0, electronic_box, FF);
	button = add_radio_button("Density Of States", (gpointer) plot_dos,NULL);
	if (plot->type == PLOT_DOS) gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(button), TRUE);
	if (!(plot->plot_mask&PLOT_DOS)) gtk_widget_set_sensitive(button,FALSE);
	button = add_radio_button("Band structure", (gpointer) plot_band,NULL);
	if (plot->type == PLOT_BAND) gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(button), TRUE);
	if (!(plot->plot_mask&PLOT_BAND)) gtk_widget_set_sensitive(button,FALSE);
	button = add_radio_button("Band / DOS", (gpointer) plot_bandos,NULL);
	if (plot->type == PLOT_BANDOS) gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(button), TRUE);
	if (!(plot->plot_mask&PLOT_BANDOS)) gtk_widget_set_sensitive(button,FALSE);
/* Set availability */
/* TODO: DOS/BAND availability */
	gtk_widget_set_sensitive(electronic_box, FALSE);
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
        button = add_radio_button("Vibrational", (gpointer) plot_frequency,NULL);
        if (plot->type == PLOT_FREQUENCY) gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(button), TRUE);
	if (!(plot->plot_mask&PLOT_FREQUENCY)) gtk_widget_set_sensitive(button,FALSE);
        button = add_radio_button("Raman", (gpointer) plot_raman,NULL);
        if (plot->type == PLOT_RAMAN) gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(button), TRUE);
	if (!(plot->plot_mask&PLOT_RAMAN)) gtk_widget_set_sensitive(button,FALSE);
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
	gtk_widget_show_all(window);
	sysenv.refresh_dialog=TRUE;
	model_prep(data);
	gui_relation_update(data);
	data->redraw=1;/*we need a redraw*/
	tree_model_refresh(data);
	redraw_canvas(ALL);
}
