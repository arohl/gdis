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

/* simple VASP launcher */
#define BREAKPOINT __asm__("int $3")

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <math.h>
#include <sys/stat.h>
#include <sys/types.h>


#include "gdis.h"
#include "coords.h"
#include "model.h"
#include "file.h"
#include "matrix.h"
#include "module.h"
#include "numeric.h"
#include "parse.h"
#include "project.h"
#include "render.h"
#include "spatial.h"
#include "quaternion.h"
#include "task.h"
#include "interface.h"
#include "dialog.h"
#include "gui_shorts.h"
#include "gui_vasp.h"

extern struct sysenv_pak sysenv;

/*globals variables*/
struct vasp_calc_gui vasp_gui;
vasp_calc_struct vasp_calc;

/*************************/
/* initialize gui values */
/*************************/
void gui_vasp_init(){
	if(vasp_gui.result_model!=NULL) g_free(vasp_gui.result_model);
	vasp_gui.result_model=g_malloc(sizeof(struct model_pak));
	model_init(vasp_gui.result_model);
	vasp_gui.have_result=FALSE;
	if(!vasp_gui.have_xml){
		vasp_gui.rions=TRUE;
		vasp_gui.rshape=FALSE;
		vasp_gui.rvolume=FALSE;
		vasp_gui.poscar_dirty=TRUE;
	}
}
/*********************************************************/
/* Load calculation setting from an included vasprun.xml */
/*********************************************************/
void load_vasprun_dialog(GtkButton *button, gpointer data){
  GtkWidget *file_chooser;
  GtkFileFilter *filter;
/*set filter*/
  filter = gtk_file_filter_new();
  gtk_file_filter_add_pattern(filter, "*vasprun.xml");
  gtk_file_filter_set_name (filter,"vasprun.xml files");
  file_chooser = gtk_file_chooser_dialog_new("Select a vasprun.xml File",GTK_WINDOW(vasp_gui.window),GTK_FILE_CHOOSER_ACTION_OPEN,GTK_STOCK_CANCEL,GTK_RESPONSE_CANCEL,GTK_STOCK_OPEN,GTK_RESPONSE_ACCEPT,NULL);
gtk_file_chooser_add_filter (GTK_FILE_CHOOSER(file_chooser),filter);/*apply filter*/
  if (gtk_dialog_run (GTK_DIALOG (file_chooser)) == GTK_RESPONSE_ACCEPT)
  {
    char *filename;
    filename = gtk_file_chooser_get_filename (GTK_FILE_CHOOSER (file_chooser));
    if(filename) {
	gtk_entry_set_text(GTK_ENTRY(vasp_gui.file_entry),g_strdup_printf("%s",filename));
	vasprun_update(filename,&vasp_calc);/*<- update everything according to the new vasprun.xml*/
	g_free (filename);
    }
  }
  gtk_widget_destroy (GTK_WIDGET(file_chooser));
}
/**************************/
/* load mpirun executable */
/**************************/
void load_mpirun_dialog(GtkButton *button, gpointer data){
  GtkWidget *file_chooser;
  GtkFileFilter *filter;
/*set filter*/
  filter = gtk_file_filter_new();
  gtk_file_filter_add_pattern(filter, "mpirun");
  gtk_file_filter_set_name (filter,"mpirun exec");
  file_chooser = gtk_file_chooser_dialog_new("Select mpirun Executable",GTK_WINDOW(vasp_gui.window),GTK_FILE_CHOOSER_ACTION_OPEN,GTK_STOCK_CANCEL,GTK_RESPONSE_CANCEL,GTK_STOCK_OPEN,GTK_RESPONSE_ACCEPT,NULL);
gtk_file_chooser_add_filter (GTK_FILE_CHOOSER(file_chooser),filter);/*apply filter*/
  if (gtk_dialog_run (GTK_DIALOG (file_chooser)) == GTK_RESPONSE_ACCEPT)
  {
    char *filename;
    filename = gtk_file_chooser_get_filename (GTK_FILE_CHOOSER (file_chooser));
    if(filename) {
        gtk_entry_set_text(GTK_ENTRY(vasp_gui.job_mpirun),g_strdup_printf("%s",filename));
	/* sync job_mpirun */
	VASP_REG_TEXT(job_mpirun);
        g_free (filename);
    }
  }
  gtk_widget_destroy (GTK_WIDGET(file_chooser));
}
/****************************/
/* load the vasp executable */
/****************************/
void load_vasp_exe_dialog(GtkButton *button, gpointer data){
  GtkWidget *file_chooser;
  file_chooser = gtk_file_chooser_dialog_new("Select VASP Executable",GTK_WINDOW(vasp_gui.window),GTK_FILE_CHOOSER_ACTION_OPEN,GTK_STOCK_CANCEL,GTK_RESPONSE_CANCEL,GTK_STOCK_OPEN,GTK_RESPONSE_ACCEPT,NULL);
  if (gtk_dialog_run (GTK_DIALOG (file_chooser)) == GTK_RESPONSE_ACCEPT)
  {
    char *filename;
    filename = gtk_file_chooser_get_filename (GTK_FILE_CHOOSER (file_chooser));
    if(filename) {
        gtk_entry_set_text(GTK_ENTRY(vasp_gui.job_vasp_exe),g_strdup_printf("%s",filename));
	/*sync job_vasp_exe*/
	VASP_REG_TEXT(job_vasp_exe);
        g_free (filename);
    }
  }
  gtk_widget_destroy (GTK_WIDGET(file_chooser));
}
/*****************************************/
/* point to the path to run the VASP job */
/*****************************************/
void load_path_dialog(GtkButton *button, gpointer data){
  GtkWidget *file_chooser;
file_chooser=gtk_file_chooser_dialog_new("Select working folder",GTK_WINDOW(vasp_gui.window),GTK_FILE_CHOOSER_ACTION_OPEN,GTK_STOCK_CANCEL,GTK_RESPONSE_CANCEL,GTK_STOCK_OPEN,GTK_RESPONSE_ACCEPT,NULL);
  gtk_file_chooser_set_action (GTK_FILE_CHOOSER(file_chooser),GTK_FILE_CHOOSER_ACTION_SELECT_FOLDER);
  if (gtk_dialog_run (GTK_DIALOG (file_chooser)) == GTK_RESPONSE_ACCEPT)
  {
    char *foldername = gtk_file_chooser_get_filename (GTK_FILE_CHOOSER (file_chooser));
    if(foldername) {
        gtk_entry_set_text(GTK_ENTRY(vasp_gui.job_path),g_strdup_printf("%s",foldername));
	/*sync job_path*/
	VASP_REG_TEXT(job_path);
        g_free (foldername);
    }
  }
  gtk_widget_destroy (GTK_WIDGET(file_chooser));
}
/******************/
/* selecting PREC */
/******************/
void vasp_prec_selected(GtkWidget *w, struct model_pak *model){
/*how much more simple it could be?*/
	const gchar *label = gtk_entry_get_text(GTK_ENTRY(w));
if (g_ascii_strcasecmp(label,"Normal") == 0) vasp_calc.prec=VP_NORM;
else if (g_ascii_strcasecmp(label,"Single") == 0) vasp_calc.prec=VP_SINGLE;
else if (g_ascii_strcasecmp(label,"Accurate") == 0)vasp_calc.prec=VP_ACCURATE;
else if (g_ascii_strcasecmp(label,"High") == 0)vasp_calc.prec=VP_HIGH;
else if (g_ascii_strcasecmp(label,"Medium") == 0)vasp_calc.prec=VP_MED;
else if (g_ascii_strcasecmp(label,"Low") == 0) vasp_calc.prec=VP_LOW;
}
/*****************************************/
/* Toggle use of automatic prec settings */
/*****************************************/
void prec_toggle(GtkWidget *w, GtkWidget *box){
	if(!(vasp_calc.use_prec)){
		/*switch to manual values*/
		gtk_widget_set_sensitive(vasp_gui.encut,TRUE);
		gtk_widget_set_sensitive(vasp_gui.enaug,TRUE);
		gtk_entry_set_text(GTK_ENTRY(vasp_gui.encut),g_strdup_printf("%.2lf",vasp_calc.encut));
		gtk_entry_set_text(GTK_ENTRY(vasp_gui.enaug),g_strdup_printf("%.2lf",vasp_calc.enaug));
	}else{
		/*switch to automatic values*/
		gtk_widget_set_sensitive(vasp_gui.encut,FALSE);
		gtk_widget_set_sensitive(vasp_gui.enaug,FALSE);
	}


}
/******************/
/* selecting ALGO */
/******************/
void vasp_algo_selected(GtkWidget *w, struct model_pak *model){
	const gchar *label = gtk_entry_get_text(GTK_ENTRY(w));
gtk_widget_set_sensitive(vasp_gui.ialgo,FALSE);
if (g_ascii_strcasecmp(label,"NORMAL") == 0) vasp_calc.algo=VA_NORM;
else if (g_ascii_strcasecmp(label,"USE_IALGO") == 0) vasp_calc.algo=VA_IALGO;
else if (g_ascii_strcasecmp(label,"VERYFAST") == 0) vasp_calc.algo=VA_VERYFAST;
else if (g_ascii_strcasecmp(label,"FAST") == 0) vasp_calc.algo=VA_FAST;
else if (g_ascii_strcasecmp(label,"CONJUGATE") == 0) vasp_calc.algo=VA_CONJ;
else if (g_ascii_strcasecmp(label,"ALL") == 0) vasp_calc.algo=VA_ALL;
else if (g_ascii_strcasecmp(label,"DAMPED") == 0) vasp_calc.algo=VA_DAMPED;
else if (g_ascii_strcasecmp(label,"SUBROT") == 0) vasp_calc.algo=VA_SUBROT;
else if (g_ascii_strcasecmp(label,"EIGENVAL") == 0) vasp_calc.algo=VA_EIGEN;
else if (g_ascii_strcasecmp(label,"NONE") == 0) vasp_calc.algo=VA_NONE;
else if (g_ascii_strcasecmp(label,"NOTHING") == 0) vasp_calc.algo=VA_NOTHING;
else if (g_ascii_strcasecmp(label,"EXACT") == 0) vasp_calc.algo=VA_EXACT;
else if (g_ascii_strcasecmp(label,"DIAG") == 0) vasp_calc.algo=VA_DIAG;
if(vasp_calc.algo==VA_IALGO) gtk_widget_set_sensitive(vasp_gui.ialgo,TRUE);
}
/*******************/
/* selecting IALGO */
/*******************/
void vasp_ialgo_selected(GtkWidget *w, struct model_pak *model){
	const gchar *label = gtk_entry_get_text(GTK_ENTRY(w));
if (g_ascii_strcasecmp(label,"2:FIXED ORB/1E") == 0) vasp_calc.ialgo=VIA_OE_FIXED;
else if (g_ascii_strcasecmp(label,"3:FIXED ORB") == 0) vasp_calc.ialgo=VIA_O_FIXED;
else if (g_ascii_strcasecmp(label,"4:SUBROT ONLY") == 0) vasp_calc.ialgo=VIA_SUBROT;
else if (g_ascii_strcasecmp(label,"5:STEEPEST") == 0) vasp_calc.ialgo=VIA_STEEP;
else if (g_ascii_strcasecmp(label,"6:CONJUGUATE GRADIENTS") == 0) vasp_calc.ialgo=VIA_CG;
else if (g_ascii_strcasecmp(label,"7:PRECOND. STEEPEST") == 0) vasp_calc.ialgo=VIA_PSTEEP;
else if (g_ascii_strcasecmp(label,"8:PRECOND. CG") == 0) vasp_calc.ialgo=VIA_PCG;
else if (g_ascii_strcasecmp(label,"38:KOSUGI") == 0) vasp_calc.ialgo=VIA_KOSUGI;
else if (g_ascii_strcasecmp(label,"44:RMM STEEPEST") == 0) vasp_calc.ialgo=VIA_ESTEEP;
else if (g_ascii_strcasecmp(label,"46:RMM PRECOND.") == 0) vasp_calc.ialgo=VIA_RMMP;
else if (g_ascii_strcasecmp(label,"48:PRECOND. RMM") == 0) vasp_calc.ialgo=VIA_PRMM;
else if (g_ascii_strcasecmp(label,"53:VAR DAMPED MD") == 0) vasp_calc.ialgo=VIA_VAR_DAMP;
else if (g_ascii_strcasecmp(label,"54:VAR DAMPED QUENCH MD") == 0) vasp_calc.ialgo=VIA_VAR_QUENCH;
else if (g_ascii_strcasecmp(label,"58:VAR PRECOND. CG") == 0) vasp_calc.ialgo=VIA_VAR_PCG;
else if (g_ascii_strcasecmp(label,"90:EXACT") == 0) vasp_calc.ialgo=VIA_EXACT;
}
/***************/
/* elec toggle */
/***************/
void elec_toggle(GtkWidget *w, GtkWidget *box){
	if(!(vasp_calc.auto_elec)){
		/*switch to manual values*/
		gtk_widget_set_sensitive(vasp_gui.nsim,TRUE);
		gtk_widget_set_sensitive(vasp_gui.vtime,TRUE);
		gtk_widget_set_sensitive(vasp_gui.nbands,TRUE);
		gtk_widget_set_sensitive(vasp_gui.nelect,TRUE);
		gtk_widget_set_sensitive(vasp_gui.iwavpr,TRUE);
		gtk_entry_set_text(GTK_ENTRY(vasp_gui.nsim),g_strdup_printf("%i",vasp_calc.nsim));
		gtk_entry_set_text(GTK_ENTRY(vasp_gui.vtime),g_strdup_printf("%.4lf",vasp_calc.vtime));
		gtk_entry_set_text(GTK_ENTRY(vasp_gui.nbands),g_strdup_printf("%i",vasp_calc.nbands));
		gtk_entry_set_text(GTK_ENTRY(vasp_gui.nelect),g_strdup_printf("%.4lf",vasp_calc.nelect));
		gtk_entry_set_text(GTK_ENTRY(vasp_gui.iwavpr),g_strdup_printf("%i",vasp_calc.iwavpr));
	}else{
		/*switch to automatic values*/
		gtk_widget_set_sensitive(vasp_gui.nsim,FALSE);
		gtk_widget_set_sensitive(vasp_gui.vtime,FALSE);
		gtk_widget_set_sensitive(vasp_gui.nbands,FALSE);
		gtk_widget_set_sensitive(vasp_gui.nelect,FALSE);
		gtk_widget_set_sensitive(vasp_gui.iwavpr,FALSE);
	}
}
/*******************/
/* selecting LREAL */
/*******************/
void vasp_lreal_selected(GtkWidget *w, struct model_pak *model){
	const gchar *label = gtk_entry_get_text(GTK_ENTRY(w));
if (g_ascii_strcasecmp(label,"Auto") == 0) vasp_calc.lreal=VLR_AUTO;
else if (g_ascii_strcasecmp(label,"On") == 0) vasp_calc.lreal=VLR_ON;
else if (g_ascii_strcasecmp(label,"TRUE") == 0) vasp_calc.lreal=VLR_TRUE;
else if (g_ascii_strcasecmp(label,"FALSE") == 0) vasp_calc.lreal=VLR_FALSE;
}
/***************/
/* grid toggle */
/***************/
void grid_toggle(GtkWidget *w, GtkWidget *box){
        if(!(vasp_calc.auto_grid)){
                /*switch to manual values*/
                gtk_widget_set_sensitive(vasp_gui.ngx,TRUE);
                gtk_widget_set_sensitive(vasp_gui.ngy,TRUE);
                gtk_widget_set_sensitive(vasp_gui.ngz,TRUE);
                gtk_widget_set_sensitive(vasp_gui.ngxf,TRUE);
                gtk_widget_set_sensitive(vasp_gui.ngyf,TRUE);
                gtk_widget_set_sensitive(vasp_gui.ngzf,TRUE);
                gtk_entry_set_text(GTK_ENTRY(vasp_gui.ngx),g_strdup_printf("%i",vasp_calc.ngx));
		gtk_entry_set_text(GTK_ENTRY(vasp_gui.ngy),g_strdup_printf("%i",vasp_calc.ngy));
		gtk_entry_set_text(GTK_ENTRY(vasp_gui.ngz),g_strdup_printf("%i",vasp_calc.ngz));
		gtk_entry_set_text(GTK_ENTRY(vasp_gui.ngxf),g_strdup_printf("%i",vasp_calc.ngxf));
		gtk_entry_set_text(GTK_ENTRY(vasp_gui.ngyf),g_strdup_printf("%i",vasp_calc.ngyf));
		gtk_entry_set_text(GTK_ENTRY(vasp_gui.ngzf),g_strdup_printf("%i",vasp_calc.ngzf));
        }else{
                /*switch to automatic values*/
                gtk_widget_set_sensitive(vasp_gui.ngx,FALSE);
                gtk_widget_set_sensitive(vasp_gui.ngy,FALSE);
                gtk_widget_set_sensitive(vasp_gui.ngz,FALSE);
                gtk_widget_set_sensitive(vasp_gui.ngxf,FALSE);
                gtk_widget_set_sensitive(vasp_gui.ngyf,FALSE);
                gtk_widget_set_sensitive(vasp_gui.ngzf,FALSE);
        }
}
/***************************/
/* ismear value is changed */
/***************************/
void ismear_changed(GtkWidget *w, struct model_pak *model){
        const gchar *label = gtk_entry_get_text(GTK_ENTRY(w));
        sscanf(label,"%i",&(vasp_calc.ismear));
	if(vasp_calc.ismear==-5) gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(vasp_gui.have_tetra),TRUE);
	else gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(vasp_gui.have_tetra),FALSE);
}
/*****************************/
/* kspacing value is changed */
/*****************************/
void kspacing_changed(GtkWidget *w, struct model_pak *model){
	const gchar *label = gtk_entry_get_text(GTK_ENTRY(w));
        sscanf(label,"%lf",&(vasp_calc.kspacing));
}
/************************/
/* toggle LNONCOLLINEAR */
/************************/
void lnoncoll_toggle(GtkWidget *w, GtkWidget *box){
        /* unselecting LNONCOLLINEAR automatically unselect LSORBIT */
        if(!vasp_calc.non_collinear) gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(vasp_gui.lsorbit),FALSE);
}
/******************/
/* toggle LSORBIT */
/******************/
void lsorbit_toggle(GtkWidget *w, GtkWidget *box){
	/* LSORBIT automatically select LNONCOLLINEAR */
	if(vasp_calc.lsorbit) gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(vasp_gui.lnoncoll),TRUE);
}
/*****************/
/* selecting GGA */
/*****************/
void vasp_gga_selected(GtkWidget *w, struct model_pak *model){
        const gchar *label = gtk_entry_get_text(GTK_ENTRY(w));
if (g_ascii_strcasecmp(label,"91:PW91") == 0) vasp_calc.gga=VG_91;
else if (g_ascii_strcasecmp(label,"PE:PBE") == 0) vasp_calc.gga=VG_PE;
else if (g_ascii_strcasecmp(label,"RP:rPBE") == 0) vasp_calc.gga=VG_RP;
else if (g_ascii_strcasecmp(label,"AM:AM05") == 0) vasp_calc.gga=VG_AM;
else if (g_ascii_strcasecmp(label,"PS:PBEsol") == 0) vasp_calc.gga=VG_PS;
}
void vasp_metagga_selected(GtkWidget *w, struct model_pak *model){
        const gchar *label = gtk_entry_get_text(GTK_ENTRY(w));
gtk_widget_set_sensitive(vasp_gui.cmbj,FALSE);
gtk_widget_set_sensitive(vasp_gui.cmbja,FALSE);
gtk_widget_set_sensitive(vasp_gui.cmbjb,FALSE);
if (g_ascii_strcasecmp(label,"TPSS") == 0) vasp_calc.mgga=VMG_TPSS;
else if (g_ascii_strcasecmp(label,"rTPSS") == 0) vasp_calc.mgga=VMG_RTPSS;
else if (g_ascii_strcasecmp(label,"M06L") == 0) vasp_calc.mgga=VMG_M06L;
else if (g_ascii_strcasecmp(label,"MBJ") == 0) vasp_calc.mgga=VMG_MBJ;
if(vasp_calc.mgga==VMG_MBJ){
	gtk_widget_set_sensitive(vasp_gui.cmbj,TRUE);
	gtk_widget_set_sensitive(vasp_gui.cmbja,TRUE);
	gtk_widget_set_sensitive(vasp_gui.cmbjb,TRUE);
}
}
/*******************/
/* toggle LMETAGGA */
/*******************/
void metagga_toggle(GtkWidget *w, GtkWidget *box){
        if(!vasp_calc.lmetagga){
                /*switch to manual values*/
                gtk_widget_set_sensitive(vasp_gui.metagga,FALSE);
		gtk_widget_set_sensitive(vasp_gui.cmbj,FALSE);
		gtk_widget_set_sensitive(vasp_gui.cmbja,FALSE);
		gtk_widget_set_sensitive(vasp_gui.cmbjb,FALSE);
        }else{
                /*switch to automatic values*/
                gtk_widget_set_sensitive(vasp_gui.metagga,TRUE);
		gtk_entry_set_text(GTK_ENTRY(GTK_COMBO(vasp_gui.metagga)->entry),"MBJ");/*actually, there is no default*/
		gtk_widget_set_sensitive(vasp_gui.cmbj,TRUE);
		gtk_widget_set_sensitive(vasp_gui.cmbja,TRUE);
		gtk_widget_set_sensitive(vasp_gui.cmbjb,TRUE);
        }
}
/*******************/
/* toggle L(S)DA+U */
/*******************/
void ldau_toggle(GtkWidget *w, GtkWidget *box){
        if(!vasp_calc.ldau){
                /*switch to manual values*/
                gtk_widget_set_sensitive(vasp_gui.ldau,FALSE);
		gtk_widget_set_sensitive(vasp_gui.ldau_print,FALSE);
        }else{
                /*switch to automatic values*/
                gtk_widget_set_sensitive(vasp_gui.ldau,TRUE);
		gtk_widget_set_sensitive(vasp_gui.ldau_print,TRUE);
		gtk_entry_set_text(GTK_ENTRY(GTK_COMBO(vasp_gui.ldau)->entry),"2:LSDA+U Dudarev");
		gtk_entry_set_text(GTK_ENTRY(GTK_COMBO(vasp_gui.ldau_print)->entry),"0:Silent");
        }
}
/***************************/
/* selecting L(S)DA+U type */
/***************************/
void vasp_ldau_selected(GtkWidget *w, struct model_pak *model){
	const gchar *label = gtk_entry_get_text(GTK_ENTRY(w));
if (g_ascii_strcasecmp(label,"1:LSDA+U Liechtenstein") == 0) vasp_calc.ldau_type=VU_1;
else if (g_ascii_strcasecmp(label,"2:LSDA+U Dudarev") == 0) vasp_calc.ldau_type=VU_2;
else if (g_ascii_strcasecmp(label,"4:LDA+U Liechtenstein") == 0) vasp_calc.ldau_type=VU_4;
}
/***************************************/
/* selecting L(S)DA+U output verbosity */
/***************************************/
void vasp_ldau_print_selected(GtkWidget *w, struct model_pak *model){
	const gchar *label = gtk_entry_get_text(GTK_ENTRY(w));
if (g_ascii_strcasecmp(label,"0:Silent") == 0) vasp_calc.ldau_output=VUO_0;
else if (g_ascii_strcasecmp(label,"1:Occupancy matrix") == 0) vasp_calc.ldau_output=VUO_1;
else if (g_ascii_strcasecmp(label,"2:Full output") == 0) vasp_calc.ldau_output=VUO_2;
}
/************************************/
/* toggle manual selection of mixer */
/************************************/
void mixer_toggle(GtkWidget *w, GtkWidget *box){
        if((vasp_calc.auto_mixer)){
                /*switch to manual values*/
                gtk_widget_set_sensitive(vasp_gui.mixer,FALSE);
		gtk_widget_set_sensitive(vasp_gui.amix,FALSE);
		gtk_widget_set_sensitive(vasp_gui.bmix,FALSE);
		gtk_widget_set_sensitive(vasp_gui.amin,FALSE);
		gtk_widget_set_sensitive(vasp_gui.amix_mag,FALSE);
		gtk_widget_set_sensitive(vasp_gui.bmix_mag,FALSE);
		gtk_widget_set_sensitive(vasp_gui.maxmix,FALSE);
		gtk_widget_set_sensitive(vasp_gui.wc,FALSE);
		gtk_widget_set_sensitive(vasp_gui.inimix,FALSE);
		gtk_widget_set_sensitive(vasp_gui.mixpre,FALSE);
        }else{
                /*switch to automatic values*/
                gtk_widget_set_sensitive(vasp_gui.mixer,TRUE);
                gtk_widget_set_sensitive(vasp_gui.amix,TRUE);
                gtk_widget_set_sensitive(vasp_gui.bmix,TRUE);
                gtk_widget_set_sensitive(vasp_gui.amin,TRUE);
                gtk_widget_set_sensitive(vasp_gui.amix_mag,TRUE);
                gtk_widget_set_sensitive(vasp_gui.bmix_mag,TRUE);
		gtk_widget_set_sensitive(vasp_gui.maxmix,TRUE);
                gtk_entry_set_text(GTK_ENTRY(GTK_COMBO(vasp_gui.mixer)->entry),"4:BROYDEN");
		gtk_entry_set_text(GTK_ENTRY(vasp_gui.amix),g_strdup_printf("%.4lf",vasp_calc.amix));
		gtk_entry_set_text(GTK_ENTRY(vasp_gui.bmix),g_strdup_printf("%.4lf",vasp_calc.bmix));
		gtk_entry_set_text(GTK_ENTRY(vasp_gui.amin),g_strdup_printf("%.4lf",vasp_calc.amin));
		gtk_entry_set_text(GTK_ENTRY(vasp_gui.amix_mag),g_strdup_printf("%.4lf",vasp_calc.amix_mag));
		gtk_entry_set_text(GTK_ENTRY(vasp_gui.bmix_mag),g_strdup_printf("%.4lf",vasp_calc.bmix_mag));
		gtk_entry_set_text(GTK_ENTRY(vasp_gui.maxmix),g_strdup_printf("%i",vasp_calc.maxmix));
		if(vasp_calc.imix==VM_4){
			gtk_widget_set_sensitive(vasp_gui.wc,TRUE);
			gtk_widget_set_sensitive(vasp_gui.inimix,TRUE);
			gtk_widget_set_sensitive(vasp_gui.mixpre,TRUE);
		}
		
        }
}
/*******************/
/* selecting mixer */
/*******************/
void vasp_mixer_selected(GtkWidget *w, struct model_pak *model){
        const gchar *label = gtk_entry_get_text(GTK_ENTRY(w));
gtk_widget_set_sensitive(vasp_gui.wc,FALSE);
gtk_widget_set_sensitive(vasp_gui.inimix,FALSE);
gtk_widget_set_sensitive(vasp_gui.mixpre,FALSE);
if (g_ascii_strcasecmp(label,"0:NO MIXING") == 0) vasp_calc.imix=VM_0;
else if (g_ascii_strcasecmp(label,"1:KERKER") == 0) vasp_calc.imix=VM_1;
else if (g_ascii_strcasecmp(label,"2:TCHEBYCHEV") == 0) vasp_calc.imix=VM_2;
else if (g_ascii_strcasecmp(label,"4:BROYDEN") == 0) vasp_calc.imix=VM_4;
if(vasp_calc.imix==VM_4) {
	gtk_widget_set_sensitive(vasp_gui.wc,TRUE);
	gtk_widget_set_sensitive(vasp_gui.inimix,TRUE);
	gtk_widget_set_sensitive(vasp_gui.mixpre,TRUE);
}
}
/***************************/
/* selecting initial mixer */
/***************************/
void vasp_inimix_selected(GtkWidget *w, struct model_pak *model){
	const gchar *label = gtk_entry_get_text(GTK_ENTRY(w));
if (g_ascii_strcasecmp(label,"0:LINEAR") == 0) vasp_calc.inimix=VIM_0;
else if (g_ascii_strcasecmp(label,"1:KERKER") == 0) vasp_calc.inimix=VIM_1;
}
/************************/
/* selecting pre-mixing */
/************************/
void vasp_mixpre_selected(GtkWidget *w, struct model_pak *model){
	const gchar *label = gtk_entry_get_text(GTK_ENTRY(w));
if (g_ascii_strcasecmp(label,"0:NONE") == 0) vasp_calc.mixpre=VMP_0;
else if (g_ascii_strcasecmp(label,"1:F=20 INVERSE KERKER") == 0) vasp_calc.mixpre=VMP_1;
else if (g_ascii_strcasecmp(label,"2:F=200 INVERSE KERKER") == 0) vasp_calc.mixpre=VMP_2;
}
/********************/
/* selecting idipol */
/********************/
void vasp_idipol_selected(GtkWidget *w, struct model_pak *model){
	const gchar *label = gtk_entry_get_text(GTK_ENTRY(w));
gtk_widget_set_sensitive(vasp_gui.ldipol,FALSE);
gtk_widget_set_sensitive(vasp_gui.lmono,FALSE);
gtk_widget_set_sensitive(vasp_gui.dipol,FALSE);
gtk_widget_set_sensitive(vasp_gui.epsilon,FALSE);
gtk_widget_set_sensitive(vasp_gui.efield,FALSE);
if (g_ascii_strcasecmp(label,"0:no calcul") == 0) vasp_calc.idipol=VID_0;
else if (g_ascii_strcasecmp(label,"1:u axis") == 0) vasp_calc.idipol=VID_1;
else if (g_ascii_strcasecmp(label,"2:v axis") == 0) vasp_calc.idipol=VID_2;
else if (g_ascii_strcasecmp(label,"3:w axis") == 0) vasp_calc.idipol=VID_3;
else if (g_ascii_strcasecmp(label,"4:all axis") == 0) vasp_calc.idipol=VID_4;
/*setting other than 0 should trigger dipol*/
if(vasp_calc.idipol!=VID_0) {
	gtk_widget_set_sensitive(vasp_gui.ldipol,TRUE);
	gtk_widget_set_sensitive(vasp_gui.lmono,TRUE);
	gtk_widget_set_sensitive(vasp_gui.dipol,TRUE);
	gtk_widget_set_sensitive(vasp_gui.epsilon,TRUE);
	if(vasp_calc.idipol!=VID_4) gtk_widget_set_sensitive(vasp_gui.efield,TRUE);
}
}
/********************/
/* selecting ibrion */
/********************/
void vasp_ibrion_selected(GtkWidget *w, struct model_pak *model){
	const gchar *label = gtk_entry_get_text(GTK_ENTRY(w));
gtk_widget_set_sensitive(vasp_gui.tebeg,FALSE);
gtk_widget_set_sensitive(vasp_gui.teend,FALSE);
gtk_widget_set_sensitive(vasp_gui.smass,FALSE);
gtk_widget_set_sensitive(vasp_gui.nblock,FALSE);
gtk_widget_set_sensitive(vasp_gui.kblock,FALSE);
gtk_widget_set_sensitive(vasp_gui.npaco,FALSE);
gtk_widget_set_sensitive(vasp_gui.apaco,FALSE);
if (g_ascii_strcasecmp(label,"-1:Nothing") == 0) vasp_calc.ibrion=-1;
else if (g_ascii_strcasecmp(label,"0:Molecular Dynamics") == 0) vasp_calc.ibrion=0;
else if (g_ascii_strcasecmp(label,"1:Quasi-Newton") == 0) vasp_calc.ibrion=1;
else if (g_ascii_strcasecmp(label,"2:Conjugate Gradients") == 0) vasp_calc.ibrion=2;
else if (g_ascii_strcasecmp(label,"3:Damped") == 0) vasp_calc.ibrion=3;
else if (g_ascii_strcasecmp(label,"5:Finite Differences") == 0) vasp_calc.ibrion=5;
else if (g_ascii_strcasecmp(label,"6:Finite Diff. with SYM") == 0) vasp_calc.ibrion=6;
else if (g_ascii_strcasecmp(label,"7:Perturbation Theory") == 0) vasp_calc.ibrion=7;
else if (g_ascii_strcasecmp(label,"8:Pert. Theory with SYM") == 0) vasp_calc.ibrion=8;
else if (g_ascii_strcasecmp(label,"44:Transition State") == 0) vasp_calc.ibrion=44;
if (vasp_calc.ibrion==0) {/*selecting MD enable molecular dynamics settings*/
gtk_widget_set_sensitive(vasp_gui.tebeg,TRUE);
gtk_widget_set_sensitive(vasp_gui.teend,TRUE);
gtk_widget_set_sensitive(vasp_gui.smass,TRUE);
gtk_widget_set_sensitive(vasp_gui.nblock,TRUE);
gtk_widget_set_sensitive(vasp_gui.kblock,TRUE);
gtk_widget_set_sensitive(vasp_gui.npaco,TRUE);
gtk_widget_set_sensitive(vasp_gui.apaco,TRUE);
}
}
/******************/
/* selecting ISIF */
/******************/
void vasp_isif_selected(GtkWidget *w, struct model_pak *model){
/*is it even going to work?*/
	const gchar *label = gtk_entry_get_text(GTK_ENTRY(w));
if (g_ascii_strcasecmp(label,"0:F_0_I_0_0") == 0) vasp_calc.isif=0;
else if (g_ascii_strcasecmp(label,"1:F_P_I_0_0") == 0) vasp_calc.isif=1;
else if (g_ascii_strcasecmp(label,"2:F_S_I_0_0") == 0) vasp_calc.isif=2;
else if (g_ascii_strcasecmp(label,"3:F_S_I_S_V") == 0) vasp_calc.isif=3;
else if (g_ascii_strcasecmp(label,"4:F_S_I_S_0") == 0) vasp_calc.isif=4;
else if (g_ascii_strcasecmp(label,"5:F_S_0_S_0") == 0) vasp_calc.isif=5;
else if (g_ascii_strcasecmp(label,"6:F_S_0_S_V") == 0) vasp_calc.isif=6;
else if (g_ascii_strcasecmp(label,"7:F_S_0_0_V") == 0) vasp_calc.isif=7;
switch (vasp_calc.isif){
	case 0:
	case 1:
	case 2:
		gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(vasp_gui.relax_ions),TRUE);
		gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(vasp_gui.relax_shape),FALSE);
		gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(vasp_gui.relax_volume),FALSE);
		break;
	case 3:
		gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(vasp_gui.relax_ions),TRUE);
		gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(vasp_gui.relax_shape),TRUE);
		gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(vasp_gui.relax_volume),TRUE);
		break;
	case 4:
		gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(vasp_gui.relax_ions),TRUE);
		gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(vasp_gui.relax_shape),TRUE);
		gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(vasp_gui.relax_volume),FALSE);
		break;
	case 5:
		gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(vasp_gui.relax_ions),FALSE);
		gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(vasp_gui.relax_shape),TRUE);
		gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(vasp_gui.relax_volume),FALSE);
		break;
	case 6:
		gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(vasp_gui.relax_ions),FALSE);
		gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(vasp_gui.relax_shape),TRUE);
		gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(vasp_gui.relax_volume),TRUE);
		break;
	case 7:
		gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(vasp_gui.relax_ions),FALSE);
		gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(vasp_gui.relax_shape),FALSE);
		gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(vasp_gui.relax_volume),TRUE);
		break;
	default:
		gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(vasp_gui.relax_ions),TRUE);
		gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(vasp_gui.relax_shape),FALSE);
		gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(vasp_gui.relax_volume),FALSE);
		vasp_calc.isif=2;/*this should never happen*/
}
}
/****************/
/* isif toggles */
/****************/
void isif_toggle(GtkWidget *w, GtkWidget *box){
/*recalculate isif function of the user switches*/
if((vasp_gui.rions)&&(!vasp_gui.rshape)&&(!vasp_gui.rvolume)) {if(vasp_calc.isif>2) vasp_calc.isif=2;}
else if((vasp_gui.rions)&&(vasp_gui.rshape)&&(vasp_gui.rvolume)) vasp_calc.isif=3;
else if ((vasp_gui.rions)&&(vasp_gui.rshape)&&(!vasp_gui.rvolume)) vasp_calc.isif=4;
else if ((!vasp_gui.rions)&&(vasp_gui.rshape)&&(!vasp_gui.rvolume)) vasp_calc.isif=5;
else if ((!vasp_gui.rions)&&(vasp_gui.rshape)&&(vasp_gui.rvolume)) vasp_calc.isif=6;
else if ((!vasp_gui.rions)&&(!vasp_gui.rshape)&&(vasp_gui.rvolume)) vasp_calc.isif=7;
/*deal with forbiden combination*/
if((vasp_gui.rions)&&(!vasp_gui.rshape)&&(vasp_gui.rvolume)) {
	gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(vasp_gui.relax_ions),FALSE);
	vasp_calc.isif=7;
}
if((!vasp_gui.rions)&&(!vasp_gui.rshape)&&(!vasp_gui.rvolume)) {
	gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(vasp_gui.relax_ions),TRUE);
	vasp_calc.isif=2;
}
switch (vasp_calc.isif){
	case 0:
		gtk_entry_set_text(GTK_ENTRY(GTK_COMBO(vasp_gui.isif)->entry),"0:F_0_I_0_0");
		break;
	case 1:
		gtk_entry_set_text(GTK_ENTRY(GTK_COMBO(vasp_gui.isif)->entry),"1:F_P_I_0_0");
		break;
	case 2:
		gtk_entry_set_text(GTK_ENTRY(GTK_COMBO(vasp_gui.isif)->entry),"2:F_S_I_0_0");
		break;
	case 3:
		gtk_entry_set_text(GTK_ENTRY(GTK_COMBO(vasp_gui.isif)->entry),"3:F_S_I_S_V");
		break;
	case 4:
		gtk_entry_set_text(GTK_ENTRY(GTK_COMBO(vasp_gui.isif)->entry),"4:F_S_I_S_0");
		break;
	case 5:
		gtk_entry_set_text(GTK_ENTRY(GTK_COMBO(vasp_gui.isif)->entry),"5:F_S_0_S_0");
		break;
	case 6:
		gtk_entry_set_text(GTK_ENTRY(GTK_COMBO(vasp_gui.isif)->entry),"6:F_S_0_S_V");
		break;
	case 7:
		gtk_entry_set_text(GTK_ENTRY(GTK_COMBO(vasp_gui.isif)->entry),"7:F_S_0_0_V");
		break;
	default:
		vasp_calc.isif=2;//this should never happen
		gtk_entry_set_text(GTK_ENTRY(GTK_COMBO(vasp_gui.isif)->entry),"2:F_S_I_0_0");
}
}
/*****************************/
/* selective dynamics toggle */
/*****************************/
void toggle_poscar_sd(){
	if((vasp_calc.poscar_free==VPF_MAN)&&(vasp_calc.poscar_sd)){
		/*we need some tags*/
		gtk_widget_set_sensitive(vasp_gui.poscar_tx,TRUE);
		gtk_widget_set_sensitive(vasp_gui.poscar_ty,TRUE);
		gtk_widget_set_sensitive(vasp_gui.poscar_tz,TRUE);
	}else{
		/*tags are no longer important*/
		gtk_widget_set_sensitive(vasp_gui.poscar_tx,FALSE);
		gtk_widget_set_sensitive(vasp_gui.poscar_ty,FALSE);
		gtk_widget_set_sensitive(vasp_gui.poscar_tz,FALSE);
	}
	gtk_widget_set_sensitive(vasp_gui.poscar_free,vasp_calc.poscar_sd);
}
/*************************/
/* selecting poscar_free */
/*************************/
void vasp_poscar_free_selected(GtkWidget *w, struct model_pak *model){
	gint index=gtk_combo_box_get_active(GTK_COMBO_BOX(vasp_gui.poscar_atoms));/*for the push-pull*/
	const gchar *label = gtk_entry_get_text(GTK_ENTRY(GTK_COMBO(vasp_gui.poscar_free)->entry));
if (g_ascii_strcasecmp(label,"All atom FIXED") == 0) vasp_calc.poscar_free=VPF_FIXED;
else if (g_ascii_strcasecmp(label,"All atom FREE") == 0) vasp_calc.poscar_free=VPF_FREE;
else if (g_ascii_strcasecmp(label,"Manual selection") == 0) vasp_calc.poscar_free=VPF_MAN;
        /* set tags sensitivity */
        gtk_widget_set_sensitive(vasp_gui.poscar_tx,(vasp_calc.poscar_free==VPF_MAN));
        gtk_widget_set_sensitive(vasp_gui.poscar_ty,(vasp_calc.poscar_free==VPF_MAN));
        gtk_widget_set_sensitive(vasp_gui.poscar_tz,(vasp_calc.poscar_free==VPF_MAN));
	if(vasp_calc.poscar_free==VPF_FREE){/*set all tags to TRUE*/
		gchar *text;
		gdouble x,y,z;
		gchar symbol[3];
		gint idx=0;
		gtk_combo_box_set_active(GTK_COMBO_BOX(vasp_gui.poscar_atoms),idx);
		text=gtk_combo_box_text_get_active_text(GTK_COMBO_BOX_TEXT(vasp_gui.poscar_atoms));
		
		while(g_ascii_strcasecmp(text,"ADD atom")!=0){
			/*modify each line tags*/
			sscanf(text,"%lf %lf %lf %*c %*c %*c ! atom: %*i (%[^)])",&x,&y,&z,&(symbol[0]));
                	gtk_combo_box_insert_text(GTK_COMBO_BOX(vasp_gui.poscar_atoms),idx,
                        	g_strdup_printf("%.8lf %.8lf %.8lf  T   T   T ! atom: %i (%s)",
                        	x,y,z,idx,symbol));
                	gtk_combo_box_text_remove(GTK_COMBO_BOX_TEXT(vasp_gui.poscar_atoms),idx+1);
			idx++;
			gtk_combo_box_set_active(GTK_COMBO_BOX(vasp_gui.poscar_atoms),idx);
			text=gtk_combo_box_text_get_active_text(GTK_COMBO_BOX_TEXT(vasp_gui.poscar_atoms));
		}
	}
	if(vasp_calc.poscar_free==VPF_FIXED){/*set all tags to FALSE*/
		gchar *text;
		gdouble x,y,z;
		gchar symbol[3];
		gint idx=0;
		gtk_combo_box_set_active(GTK_COMBO_BOX(vasp_gui.poscar_atoms),idx);
		text=gtk_combo_box_text_get_active_text(GTK_COMBO_BOX_TEXT(vasp_gui.poscar_atoms));
		while(g_ascii_strcasecmp(text,"ADD atom")!=0){
			/*modify each line tags*/
			sscanf(text,"%lf %lf %lf %*c %*c %*c ! atom: %*i (%[^)])",&x,&y,&z,&(symbol[0]));
			gtk_combo_box_insert_text(GTK_COMBO_BOX(vasp_gui.poscar_atoms),idx,
				g_strdup_printf("%.8lf %.8lf %.8lf  F   F   F ! atom: %i (%s)",
				x,y,z,idx,symbol));
			gtk_combo_box_text_remove(GTK_COMBO_BOX_TEXT(vasp_gui.poscar_atoms),idx+1);
			idx++;
			gtk_combo_box_set_active(GTK_COMBO_BOX(vasp_gui.poscar_atoms),idx);
			text=gtk_combo_box_text_get_active_text(GTK_COMBO_BOX_TEXT(vasp_gui.poscar_atoms));
		}
	}
	/*return to selected item*/
	gtk_combo_box_set_active(GTK_COMBO_BOX(vasp_gui.poscar_atoms),index);
}
/************************/
/* poscar_direct toggle */
/************************/
void toggle_poscar_direct(){
	/* TODO: recalculate all positions on a single click */
}
/*************************/
/* selecting poscar atom */
/*************************/
void vasp_poscar_atoms_selected(GtkWidget *w, struct model_pak *model){
	/* using combobox*/
	gint index=gtk_combo_box_get_active(GTK_COMBO_BOX(w));
	gchar *text=gtk_combo_box_text_get_active_text(GTK_COMBO_BOX_TEXT(w));
	gdouble x,y,z;
	gchar tx,ty,tz;
	gchar symbol[3];
	gtk_entry_set_text(GTK_ENTRY(vasp_gui.poscar_index),g_strdup_printf("%i",index));
	if (g_ascii_strcasecmp(text,"ADD atom") == 0) {
		/*allow symbol*/
		gtk_widget_set_sensitive(vasp_gui.poscar_symbol,TRUE);
	}else{
		/*populate atom properties*/
		sscanf(text,"%lf %lf %lf %c %c %c ! atom: %*i (%[^)])",&x,&y,&z,&tx,&ty,&tz,&(symbol[0]));
		gtk_entry_set_text(GTK_ENTRY(vasp_gui.poscar_symbol),g_strdup_printf("%s",symbol));
		gtk_entry_set_text(GTK_ENTRY(vasp_gui.poscar_x),g_strdup_printf("%.8lf",x));
		gtk_entry_set_text(GTK_ENTRY(vasp_gui.poscar_y),g_strdup_printf("%.8lf",y));
		gtk_entry_set_text(GTK_ENTRY(vasp_gui.poscar_z),g_strdup_printf("%.8lf",z));
		/* register all that */
		vasp_calc.poscar_index=index;
		vasp_calc.poscar_x=x;
		vasp_calc.poscar_y=y;
		vasp_calc.poscar_z=z;
		vasp_calc.poscar_tx=(tx=='T');
		vasp_calc.poscar_ty=(ty=='T');
		vasp_calc.poscar_tz=(tz=='T');
		/* disallow symbol */
		gtk_widget_set_sensitive(vasp_gui.poscar_symbol,FALSE);
	}
	/* set tags sensitivity */
	gtk_widget_set_sensitive(vasp_gui.poscar_tx,(vasp_calc.poscar_free==VPF_MAN));
	gtk_widget_set_sensitive(vasp_gui.poscar_ty,(vasp_calc.poscar_free==VPF_MAN));
	gtk_widget_set_sensitive(vasp_gui.poscar_tz,(vasp_calc.poscar_free==VPF_MAN));
}
/**************************/
/* modify/add poscar atom */
/**************************/
void vasp_atom_modified(){
	gint index=gtk_combo_box_get_active(GTK_COMBO_BOX(vasp_gui.poscar_atoms));
	gchar *text=gtk_combo_box_text_get_active_text(GTK_COMBO_BOX_TEXT(vasp_gui.poscar_atoms));
	gchar *tamp;
	gchar tx,ty,tz;

	/*mini-sync*/
	VASP_REG_VAL(poscar_index,"%i");
	VASP_REG_VAL(poscar_x,"%lf");
	VASP_REG_VAL(poscar_y,"%lf");
	VASP_REG_VAL(poscar_z,"%lf");
	sprintf(vasp_calc.poscar_symbol,"%s",gtk_entry_get_text(GTK_ENTRY(vasp_gui.poscar_symbol)));
	switch (vasp_calc.poscar_free){
	case VPF_MAN:
		if(vasp_calc.poscar_tx) tx='T';
		else tx='F';
		if(vasp_calc.poscar_ty) ty='T';
		else ty='F';
		if(vasp_calc.poscar_tz) tz='T';
		else tz='F';
		break;
	case VPF_FREE:
		tx='T';ty='T';tz='T';
		break;
	case VPF_FIXED:
	default:
		tx='F';ty='F';tz='F';
	}
	/* add/modify */
	if(g_ascii_strcasecmp(text,"ADD atom") == 0){
		/*we are going to insert some data*/
		gtk_combo_box_insert_text(GTK_COMBO_BOX(vasp_gui.poscar_atoms),index,
			g_strdup_printf("%.8lf %.8lf %.8lf  %c   %c   %c ! atom: %i (%s)",
			vasp_calc.poscar_x,vasp_calc.poscar_y,vasp_calc.poscar_z,tx,ty,tz,index,vasp_calc.poscar_symbol));
	}else{
		/*we are goign to modify some data*/
		gtk_combo_box_insert_text(GTK_COMBO_BOX(vasp_gui.poscar_atoms),index,
			g_strdup_printf("%.8lf %.8lf %.8lf  %c   %c   %c ! atom: %i (%s)",
			vasp_calc.poscar_x,vasp_calc.poscar_y,vasp_calc.poscar_z,tx,ty,tz,index,vasp_calc.poscar_symbol));
		/* Note: GTK 3.0 has gtk_combo_box_text_insert */
		gtk_combo_box_text_remove(GTK_COMBO_BOX_TEXT(vasp_gui.poscar_atoms),index+1);
	}
	/*re-select current entry (to refresh content)*/
	gtk_combo_box_set_active(GTK_COMBO_BOX(vasp_gui.poscar_atoms),index);
	vasp_gui.poscar_dirty=TRUE;
}
/**********************/
/* delete poscar atom */
/**********************/
void vasp_atom_delete(){
	gint index=gtk_combo_box_get_active(GTK_COMBO_BOX(vasp_gui.poscar_atoms));
	gchar *text=gtk_combo_box_text_get_active_text(GTK_COMBO_BOX_TEXT(vasp_gui.poscar_atoms));
	if(g_ascii_strcasecmp(text,"ADD atom") == 0) return;/*can't delete this one*/
	if(index==-1) return;/*nothing selected anyway*/
	gtk_combo_box_text_remove(GTK_COMBO_BOX_TEXT(vasp_gui.poscar_atoms),index);
	gtk_combo_box_set_active(GTK_COMBO_BOX(vasp_gui.poscar_atoms),index-1);
	vasp_gui.poscar_dirty=TRUE;
}
/**************************/
/* selecting kpoints_mode */
/**************************/
void vasp_kpoints_mode_selected(GtkWidget *w, struct model_pak *model){
	const gchar *label = gtk_combo_box_text_get_active_text(GTK_COMBO_BOX_TEXT(w));
if (g_ascii_strcasecmp(label,"Manual Entry") == 0) vasp_calc.kpoints_mode=VKP_MAN;
else if (g_ascii_strcasecmp(label,"Line Mode") == 0) vasp_calc.kpoints_mode=VKP_LINE;
else if (g_ascii_strcasecmp(label,"Automatic (M&P)") == 0) vasp_calc.kpoints_mode=VKP_AUTO;
else if (g_ascii_strcasecmp(label,"Gamma (M&P)") == 0) vasp_calc.kpoints_mode=VKP_GAMMA;
else if (g_ascii_strcasecmp(label,"Monkhorst-Pack classic") == 0) vasp_calc.kpoints_mode=VKP_MP;
else if (g_ascii_strcasecmp(label,"Basis set definition") == 0) vasp_calc.kpoints_mode=VKP_BASIS;
gtk_widget_set_sensitive(vasp_gui.kpoints_w,FALSE);
gtk_widget_set_sensitive(vasp_gui.kpoints_cart,FALSE);
gtk_widget_set_sensitive(vasp_gui.kpoints_tetra,FALSE);
gtk_widget_set_sensitive(vasp_gui.have_tetra,FALSE);
gtk_widget_set_sensitive(vasp_gui.kpoints_coord,FALSE);
gtk_widget_set_sensitive(vasp_gui.kpoints_kx,FALSE);
gtk_widget_set_sensitive(vasp_gui.kpoints_ky,FALSE);
gtk_widget_set_sensitive(vasp_gui.kpoints_kz,FALSE);
gtk_widget_set_sensitive(vasp_gui.kpoints_sx,FALSE);
gtk_widget_set_sensitive(vasp_gui.kpoints_sy,FALSE);
gtk_widget_set_sensitive(vasp_gui.kpoints_sz,FALSE);
switch (vasp_calc.kpoints_mode){
	case VKP_MAN:
	gtk_widget_set_sensitive(vasp_gui.kpoints_w,TRUE);//only VKP_MAN require weights
	if(vasp_calc.ismear==-5) gtk_widget_set_sensitive(vasp_gui.kpoints_tetra,TRUE);//only VKP_MAN require tetra
	if(vasp_calc.ismear==-5) gtk_widget_set_sensitive(vasp_gui.have_tetra,TRUE);//and only if ISMEAR=-5
	case VKP_LINE:
	case VKP_BASIS:
	gtk_widget_set_sensitive(vasp_gui.kpoints_coord,TRUE);//need coord
	gtk_widget_set_sensitive(vasp_gui.kpoints_nkpts,TRUE);//need nkpts
	gtk_widget_set_sensitive(vasp_gui.kpoints_cart,TRUE);//can be cart
	break;
	case VKP_GAMMA:
	case VKP_MP:
	gtk_widget_set_sensitive(vasp_gui.kpoints_ky,TRUE);//require ky
	gtk_widget_set_sensitive(vasp_gui.kpoints_kz,TRUE);//require kz
	gtk_widget_set_sensitive(vasp_gui.kpoints_sx,TRUE);//require sx
	gtk_widget_set_sensitive(vasp_gui.kpoints_sy,TRUE);//require sy
	gtk_widget_set_sensitive(vasp_gui.kpoints_sz,TRUE);//require sz
	case VKP_AUTO:
	gtk_widget_set_sensitive(vasp_gui.kpoints_kx,TRUE);//AUTO only require kx
	vasp_calc.kpoints_nkpts=0;
	gtk_entry_set_text(GTK_ENTRY(vasp_gui.kpoints_nkpts),g_strdup_printf("%i",vasp_calc.kpoints_nkpts));
	gtk_widget_set_sensitive(vasp_gui.kpoints_nkpts,FALSE);//NKPTS is 0 and disable
	break;
	default:
		fprintf(stderr,"NO KPOINTS GENERATION MODE SELECTED!");//should never happen
	}
}
/**********************/
/* selecting a kpoint */
/**********************/
void vasp_kpoints_kpts_selected(GtkWidget *w, struct model_pak *model){
        /* using combobox*/
        gint index=gtk_combo_box_get_active(GTK_COMBO_BOX(w));
        gchar *text=gtk_combo_box_text_get_active_text(GTK_COMBO_BOX_TEXT(w));
        gdouble x,y,z,wt;
	gtk_entry_set_text(GTK_ENTRY(vasp_gui.kpoints_i),g_strdup_printf("%i",index));
	if (g_ascii_strcasecmp(text,"ADD kpoint") == 0) {
                /*no need to update anything but index*/
	} else {
		/*update index,x,y,z,w*/
		wt=0.0;
		if(vasp_calc.kpoints_mode==VKP_LINE) sscanf(text,"%lf %lf %lf ! kpoint: %*i",&x,&y,&z);
		else sscanf(text,"%lf %lf %lf %lf ! kpoint: %*i",&x,&y,&z,&wt);
                gtk_entry_set_text(GTK_ENTRY(vasp_gui.kpoints_x),g_strdup_printf("%.8lf",x));
                gtk_entry_set_text(GTK_ENTRY(vasp_gui.kpoints_y),g_strdup_printf("%.8lf",y));
                gtk_entry_set_text(GTK_ENTRY(vasp_gui.kpoints_z),g_strdup_printf("%.8lf",z));
		gtk_entry_set_text(GTK_ENTRY(vasp_gui.kpoints_w),g_strdup_printf("%.8lf",wt));
	}
}
/*********************/
/* Add/modify kpoint */
/*********************/
void vasp_kpoint_modified(){
        gint index=gtk_combo_box_get_active(GTK_COMBO_BOX(vasp_gui.kpoints_kpts));
        gchar *text=gtk_combo_box_text_get_active_text(GTK_COMBO_BOX_TEXT(vasp_gui.kpoints_kpts));
        gchar *tamp;
        /*mini-sync*/
        VASP_REG_VAL(kpoints_i,"%i");
        VASP_REG_VAL(kpoints_x,"%lf");
        VASP_REG_VAL(kpoints_y,"%lf");
        VASP_REG_VAL(kpoints_z,"%lf");
	VASP_REG_VAL(kpoints_w,"%lf");
        /* add/modify */
        if(g_ascii_strcasecmp(text,"ADD kpoint") == 0){
                /*we are going to insert some data*/
		if(vasp_calc.kpoints_mode==VKP_LINE)
			gtk_combo_box_insert_text(GTK_COMBO_BOX(vasp_gui.kpoints_kpts),index,
                        	g_strdup_printf("%.8lf %.8lf %.8lf ! kpoint: %i",
                        	vasp_calc.kpoints_x,vasp_calc.kpoints_y,vasp_calc.kpoints_z,index));
		else
			gtk_combo_box_insert_text(GTK_COMBO_BOX(vasp_gui.kpoints_kpts),index,
				g_strdup_printf("%.8lf %.8lf %.8lf %.8lf ! kpoint: %i",
				vasp_calc.kpoints_x,vasp_calc.kpoints_y,vasp_calc.kpoints_z,vasp_calc.kpoints_w,index));
        }else{
                /*we are goign to modify some data*/
		if(vasp_calc.kpoints_mode==VKP_LINE)
			gtk_combo_box_insert_text(GTK_COMBO_BOX(vasp_gui.kpoints_kpts),index,
				g_strdup_printf("%.8lf %.8lf %.8lf ! kpoint: %i",
				vasp_calc.kpoints_x,vasp_calc.kpoints_y,vasp_calc.kpoints_z,index));
		else
			gtk_combo_box_insert_text(GTK_COMBO_BOX(vasp_gui.kpoints_kpts),index,
				g_strdup_printf("%.8lf %.8lf %.8lf %.8lf ! kpoint: %i",
				vasp_calc.kpoints_x,vasp_calc.kpoints_y,vasp_calc.kpoints_z,vasp_calc.kpoints_w,index));
                gtk_combo_box_text_remove(GTK_COMBO_BOX_TEXT(vasp_gui.kpoints_kpts),index+1);
        }
        /*re-select current entry (to refresh content)*/
        gtk_combo_box_set_active(GTK_COMBO_BOX(vasp_gui.kpoints_kpts),index);
}
/******************/
/* deleted kpoint */
/******************/
void vasp_kpoint_delete(){
        gint index=gtk_combo_box_get_active(GTK_COMBO_BOX(vasp_gui.kpoints_kpts));
        gchar *text=gtk_combo_box_text_get_active_text(GTK_COMBO_BOX_TEXT(vasp_gui.kpoints_kpts));
        if(g_ascii_strcasecmp(text,"ADD kpoint") == 0) return;/*can't delete this one*/
        if(index==-1) return;/*nothing selected anyway*/
        gtk_combo_box_text_remove(GTK_COMBO_BOX_TEXT(vasp_gui.kpoints_kpts),index);
        gtk_combo_box_set_active(GTK_COMBO_BOX(vasp_gui.kpoints_kpts),index-1);
}
/****************/
/* tetra toggle */
/****************/
void toogle_tetra(){
	/*enable/disable tetrahedron*/
	if(vasp_calc.ismear!=-5) vasp_calc.kpoints_tetra=FALSE;/*ISMEAR=-5 required*/
	if(vasp_calc.kpoints_mode!=VKP_MAN) vasp_calc.kpoints_tetra=FALSE;/*manual generation must be selected*/
	if(vasp_calc.kpoints_tetra==FALSE) gtk_widget_set_sensitive(vasp_gui.kpoints_tetra,FALSE);
	else gtk_widget_set_sensitive(vasp_gui.kpoints_tetra,TRUE);
}
/*************************/
/* selecting tetrahedron */
/*************************/
void vasp_tetra_selected(GtkWidget *w, struct model_pak *model){
	/* using combobox*/
	gint index=gtk_combo_box_get_active(GTK_COMBO_BOX(w));
	gchar *text=gtk_combo_box_text_get_active_text(GTK_COMBO_BOX_TEXT(w));
	gdouble wt;
	gint a,b,c,d;
	if (g_ascii_strcasecmp(text,"ADD tetrahedron") == 0) {
		/*no need to update anything*/
	} else {
		sscanf(text,"%lf %i %i %i %i ! tetrahedron: %*i",&wt,&a,&b,&c,&d);
		gtk_entry_set_text(GTK_ENTRY(vasp_gui.tetra_i),g_strdup_printf("%i",index));
		gtk_entry_set_text(GTK_ENTRY(vasp_gui.tetra_w),g_strdup_printf("%.4lf",wt));
		gtk_entry_set_text(GTK_ENTRY(vasp_gui.tetra_a),g_strdup_printf("%i",a));
		gtk_entry_set_text(GTK_ENTRY(vasp_gui.tetra_b),g_strdup_printf("%i",b));
		gtk_entry_set_text(GTK_ENTRY(vasp_gui.tetra_c),g_strdup_printf("%i",c));
		gtk_entry_set_text(GTK_ENTRY(vasp_gui.tetra_d),g_strdup_printf("%i",d));
	}
}
/**************************/
/* Add/modify tetrahedron */
/**************************/
void vasp_tetra_modified(){
        gint index=gtk_combo_box_get_active(GTK_COMBO_BOX(vasp_gui.tetra));
        gchar *text=gtk_combo_box_text_get_active_text(GTK_COMBO_BOX_TEXT(vasp_gui.tetra));
        gchar *tamp;

        /*mini-sync*/
        VASP_REG_VAL(tetra_i,"%i");
        VASP_REG_VAL(tetra_w,"%lf");
	VASP_REG_VAL(tetra_a,"%i");
        VASP_REG_VAL(tetra_b,"%i");
        VASP_REG_VAL(tetra_c,"%i");
        VASP_REG_VAL(tetra_d,"%i");
        /* add/modify */
        if(g_ascii_strcasecmp(text,"ADD tetrahedron") == 0){
                /*we are going to insert some data*/
		gtk_combo_box_insert_text(GTK_COMBO_BOX(vasp_gui.tetra),index,
			g_strdup_printf("%lf %i %i %i %i ! tetrahedron: %i",
			vasp_calc.tetra_w,vasp_calc.tetra_a,vasp_calc.tetra_b,vasp_calc.tetra_c,vasp_calc.tetra_d,index));
        }else{  
                /*we are goign to modify some data*/
		gtk_combo_box_insert_text(GTK_COMBO_BOX(vasp_gui.tetra),index,
			g_strdup_printf("%lf %i %i %i %i ! tetrahedron: %i",
			vasp_calc.tetra_w,vasp_calc.tetra_a,vasp_calc.tetra_b,vasp_calc.tetra_c,vasp_calc.tetra_d,index));
		gtk_combo_box_text_remove(GTK_COMBO_BOX_TEXT(vasp_gui.tetra),index+1);
        }
        /*re-select current entry (to refresh content)*/
        gtk_combo_box_set_active(GTK_COMBO_BOX(vasp_gui.tetra),index);
}
/**********************/
/* delete tetrahedron */
/**********************/
void vasp_tetra_delete(){
        gint index=gtk_combo_box_get_active(GTK_COMBO_BOX(vasp_gui.tetra));
        gchar *text=gtk_combo_box_text_get_active_text(GTK_COMBO_BOX_TEXT(vasp_gui.tetra));
        if(g_ascii_strcasecmp(text,"ADD tetrahedron") == 0) return;/*can't delete this one*/
        if(index==-1) return;/*nothing selected anyway*/
        gtk_combo_box_text_remove(GTK_COMBO_BOX_TEXT(vasp_gui.tetra),index);
        gtk_combo_box_set_active(GTK_COMBO_BOX(vasp_gui.tetra),index-1);
}
/****************************************************/
/* get species information (ie. symbol) from POTCAR */
/****************************************************/
void potcar_get_species(){
	FILE *fp;
	gchar *tamp;
	gchar *line;
	gchar sym[3];
	/*the resulting species information should match that of POSCAR*/
	fp=fopen(vasp_calc.potcar_file,"r");
        if (!fp) {
		if(vasp_calc.potcar_species!=NULL) g_free(vasp_calc.potcar_species);
		vasp_calc.potcar_species=NULL;
                return;
        }
	/*file is opened*/
	sym[2]='\0';
	if(vasp_calc.potcar_species!=NULL) g_free(vasp_calc.potcar_species);
	vasp_calc.potcar_species=g_strdup_printf(" ");
	line = file_read_line(fp);
	while(line){
		if (strstr(line,"TITEL") != NULL) {
			sscanf(line," TITEL = %*s %c%c%*s",&(sym[0]),&(sym[1]));
			tamp=g_strdup_printf("%s %s",vasp_calc.potcar_species,sym);
			g_free(vasp_calc.potcar_species);
			vasp_calc.potcar_species=tamp;
		}
		g_free(line);
		line = file_read_line(fp);
	}
	fclose(fp);
	/*update vasp_gui.potcar_species*/
	gtk_entry_set_text(GTK_ENTRY(vasp_gui.potcar_species),g_strdup_printf("%s",vasp_calc.potcar_species));
}
/************************************************************************/
/* point to a file where POTCAR information is stored for each elements */
/************************************************************************/
void load_potcar_file_dialog(GtkButton *button, gpointer data){
  GtkWidget *file_chooser;
  GtkFileFilter *filter;
/*set filter*/
  filter = gtk_file_filter_new();
  gtk_file_filter_add_pattern(filter, "*POTCAR");
  gtk_file_filter_set_name (filter,"POTCAR files");
  file_chooser = gtk_file_chooser_dialog_new("Select a POTCAR File",GTK_WINDOW(vasp_gui.window),GTK_FILE_CHOOSER_ACTION_OPEN,GTK_STOCK_CANCEL,GTK_RESPONSE_CANCEL,GTK_STOCK_OPEN,GTK_RESPONSE_ACCEPT,NULL);
gtk_file_chooser_add_filter (GTK_FILE_CHOOSER(file_chooser),filter);/*apply filter*/
  if (gtk_dialog_run (GTK_DIALOG (file_chooser)) == GTK_RESPONSE_ACCEPT)
  {
    char *filename;
    filename = gtk_file_chooser_get_filename (GTK_FILE_CHOOSER (file_chooser));
    if(filename) {
        gtk_entry_set_text(GTK_ENTRY(vasp_gui.potcar_file),g_strdup_printf("%s",filename));
/*dont forget to sync!!*/
        VASP_REG_TEXT(potcar_file);
        g_free (filename);
        potcar_get_species();
    }
  }
  gtk_widget_destroy (GTK_WIDGET(file_chooser));
}
/***************************************************************/
/* point to a folder where POTCAR are stored for each elements */
/***************************************************************/
void load_potcar_folder_dialog(GtkButton *button, gpointer data){
  GtkWidget *file_chooser;
file_chooser=gtk_file_chooser_dialog_new("Select the POTCAR folder",GTK_WINDOW(vasp_gui.window),GTK_FILE_CHOOSER_ACTION_OPEN,GTK_STOCK_CANCEL,GTK_RESPONSE_CANCEL,GTK_STOCK_OPEN,GTK_RESPONSE_ACCEPT,NULL);
  gtk_file_chooser_set_action (GTK_FILE_CHOOSER(file_chooser),GTK_FILE_CHOOSER_ACTION_SELECT_FOLDER);
  if (gtk_dialog_run (GTK_DIALOG (file_chooser)) == GTK_RESPONSE_ACCEPT)
  {
    char *foldername = gtk_file_chooser_get_filename (GTK_FILE_CHOOSER (file_chooser));
    if(foldername) {
        gtk_entry_set_text(GTK_ENTRY(vasp_gui.potcar_folder),g_strdup_printf("%s",foldername));
        g_free (foldername);
    }
  }
  gtk_widget_destroy (GTK_WIDGET(file_chooser));
}
void parallel_eq(GtkWidget *w, GtkWidget *box){
	gchar *tamp;
	/*parallel equilibration*/
	gint np;//,ncore,kpar;
	/* mini sync */
	VASP_REG_VAL(ncore,"%lf");
	VASP_REG_VAL(kpar,"%lf");
	VASP_REG_VAL(job_nproc,"%lf");
	np=(gint)vasp_calc.job_nproc;
//	ncore=(gint)vasp_calc.ncore;
//	kpar=(gint)vasp_calc.kpar;
	/*we have np CPU*/
	if(np==1){
		//NCORE=1 ; KPAR=1
		gtk_spin_button_set_value(GTK_SPIN_BUTTON(vasp_gui.ncore),1.0);
		gtk_spin_button_set_value(GTK_SPIN_BUTTON(vasp_gui.kpar),1.0);
		gtk_spin_button_set_range(GTK_SPIN_BUTTON(vasp_gui.ncore),1.0,1.0);
		gtk_spin_button_set_range(GTK_SPIN_BUTTON(vasp_gui.kpar),1.0,1.0);
		return;
	}
	/*(re)set default barrier*/
	gtk_spin_button_set_range(GTK_SPIN_BUTTON(vasp_gui.ncore),1.0,(gdouble)np);
	gtk_spin_button_set_range(GTK_SPIN_BUTTON(vasp_gui.kpar),1.0,(gdouble)np);
	/* TODO: SMART re-definition */
}
/***************************/
/* sync poscar information */
/***************************/
void vasp_poscar_sync(){
        gchar *tamp;
        gchar *text;
        gchar *text2;
        gchar *tamp2;
        gchar symbol[3];
        gint idx,jdx;
	if(vasp_gui.poscar_dirty==FALSE) return;/*ne need to update*/
/* pass_1: GROUP atom list by atom types */
        idx=0;
        gtk_combo_box_set_active(GTK_COMBO_BOX(vasp_gui.poscar_atoms),idx);
        text=gtk_combo_box_text_get_active_text(GTK_COMBO_BOX_TEXT(vasp_gui.poscar_atoms));
        while(g_ascii_strcasecmp(text,"ADD atom")!=0){
                sscanf(text,"%*f %*f %*f %*c %*c %*c ! atom: %*i (%[^)])",&(vasp_calc.poscar_symbol[0]));
                jdx=0;
                gtk_combo_box_set_active(GTK_COMBO_BOX(vasp_gui.poscar_atoms),jdx);
                text2=gtk_combo_box_text_get_active_text(GTK_COMBO_BOX_TEXT(vasp_gui.poscar_atoms));
                while(g_ascii_strcasecmp(text2,"ADD atom")!=0){
                        /**/
                        sscanf(text2,"%*f %*f %*f %*c %*c %*c ! atom: %*i (%[^)])",&(symbol[0]));
                        if(g_ascii_strcasecmp(vasp_calc.poscar_symbol,symbol) == 0){
                                gtk_combo_box_text_remove(GTK_COMBO_BOX_TEXT(vasp_gui.poscar_atoms),jdx);
                                gtk_combo_box_text_prepend_text(GTK_COMBO_BOX_TEXT(vasp_gui.poscar_atoms),text2);
                                idx++;
                        }
                        jdx++;
                        gtk_combo_box_set_active(GTK_COMBO_BOX(vasp_gui.poscar_atoms),jdx);
                        text2=gtk_combo_box_text_get_active_text(GTK_COMBO_BOX_TEXT(vasp_gui.poscar_atoms));
                }
                idx++;
                gtk_combo_box_set_active(GTK_COMBO_BOX(vasp_gui.poscar_atoms),idx);
                text=gtk_combo_box_text_get_active_text(GTK_COMBO_BOX_TEXT(vasp_gui.poscar_atoms));
        }
/* pass_2: REORDER atom list <- useful for 'consistent' ordering*/
        idx=0;
        gtk_combo_box_set_active(GTK_COMBO_BOX(vasp_gui.poscar_atoms),idx);
        text=gtk_combo_box_text_get_active_text(GTK_COMBO_BOX_TEXT(vasp_gui.poscar_atoms));
        while(g_ascii_strcasecmp(text,"ADD atom")!=0){
                gtk_combo_box_text_remove(GTK_COMBO_BOX_TEXT(vasp_gui.poscar_atoms),idx);
                gtk_combo_box_text_prepend_text(GTK_COMBO_BOX_TEXT(vasp_gui.poscar_atoms),text);
                idx++;
                gtk_combo_box_set_active(GTK_COMBO_BOX(vasp_gui.poscar_atoms),idx);
                text=gtk_combo_box_text_get_active_text(GTK_COMBO_BOX_TEXT(vasp_gui.poscar_atoms));
        }
/* pass_3: SETUP unexposed properties */
        idx=0;jdx=0;
        vasp_calc.species_total=1;
        gtk_combo_box_set_active(GTK_COMBO_BOX(vasp_gui.poscar_atoms),idx);
        text=gtk_combo_box_text_get_active_text(GTK_COMBO_BOX_TEXT(vasp_gui.poscar_atoms));
        sscanf(text,"%*f %*f %*f %*c %*c %*c ! atom: %*i (%[^)])",&(vasp_calc.poscar_symbol[0]));
        tamp=g_strdup_printf("  %s",vasp_calc.poscar_symbol);
        tamp2=g_strdup_printf(" ");
        /*1st pass: get number of species*/
        while(g_ascii_strcasecmp(text,"ADD atom")!=0){
                /*get each line*/
                sscanf(text,"%*f %*f %*f %*c %*c %*c ! atom: %*i (%[^)])",&(symbol[0]));
                if(g_ascii_strcasecmp(vasp_calc.poscar_symbol,symbol) != 0) {
                        /*new species*/
                        text=g_strdup_printf("%s %s",tamp,symbol);
                        g_free(tamp);
                        tamp=text;
                        text=g_strdup_printf("%s %i",tamp2,jdx);
                        g_free(tamp2);
                        tamp2=text;
                        vasp_calc.species_total++;
                        sprintf(vasp_calc.poscar_symbol,"%s",symbol);
                        jdx=0;
                }
                idx++;
                jdx++;
                gtk_combo_box_set_active(GTK_COMBO_BOX(vasp_gui.poscar_atoms),idx);
                text=gtk_combo_box_text_get_active_text(GTK_COMBO_BOX_TEXT(vasp_gui.poscar_atoms));
        }
        vasp_calc.atoms_total=idx;/*so we have the total number of atoms*/
        /*we need to get the last jdx*/
        text=g_strdup_printf("%s %i",tamp2,jdx);
        g_free(tamp2);
        tamp2=text;
        if(vasp_calc.species_symbols!=NULL) g_free(vasp_calc.species_symbols);
        if(vasp_calc.species_numbers!=NULL) g_free(vasp_calc.species_numbers);
        vasp_calc.species_symbols=g_strdup(tamp);g_free(tamp);
        vasp_calc.species_numbers=g_strdup(tamp2);g_free(tamp2);
	vasp_gui.poscar_dirty=FALSE;
}

/************************/
/* Switch notebook page */
/************************/
void vasp_calc_page_switch(GtkNotebook *notebook,GtkWidget *page,guint page_num,gpointer user_data){
        /* some (few) values need to be updated from a page to another*/
	if(page_num==VASP_PAGE_ELECTRONIC){
		gtk_entry_set_text(GTK_ENTRY(vasp_gui.ismear),g_strdup_printf("%i",vasp_calc.ismear));
		gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(vasp_gui.kgamma),vasp_calc.kgamma);
		gtk_entry_set_text(GTK_ENTRY(vasp_gui.kspacing),g_strdup_printf("%.4lf",vasp_calc.kspacing));
	} else if (page_num==VASP_PAGE_IONIC){
		vasp_poscar_sync();/*arrange poscar atom values*/
	} else if (page_num==VASP_PAGE_KPOINTS){
		gtk_entry_set_text(GTK_ENTRY(vasp_gui.ismear_3),g_strdup_printf("%i",vasp_calc.ismear));
		gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(vasp_gui.kgamma_3),vasp_calc.kgamma);
		gtk_entry_set_text(GTK_ENTRY(vasp_gui.kspacing_3),g_strdup_printf("%.4lf",vasp_calc.kspacing));
		toogle_tetra();/*in case we need*/
	} else if (page_num==VASP_PAGE_POTCAR){
		vasp_poscar_sync();/*we need to know species_symbols*/
		gtk_entry_set_text(GTK_ENTRY(vasp_gui.species),g_strdup_printf("%s",vasp_calc.species_symbols));
	}
	/*to be continued...*/
}
/**********************************************************/
/* convert current vasp_gui into current vasp_calc_struct */
/**********************************************************/
void vasp_gui_sync(){
	gchar *tamp;
	/* sync start here */
	VASP_REG_TEXT(name);
	//prec already sync
	VASP_REG_VAL(encut,"%lf");
	VASP_REG_VAL(enaug,"%lf");
	VASP_REG_VAL(ediff,"%lE");
	//algo already sync
	//ialgo already sync
	//ldiag already sync
	if(!vasp_calc.auto_elec){
		VASP_REG_VAL(nbands,"%i");
		VASP_REG_VAL(nelect,"%lf");
		VASP_REG_VAL(iwavpr,"%i");
		VASP_REG_VAL(nsim,"%i");
		VASP_REG_VAL(vtime,"%lf");
	}
//	VASP_REG_VAL(ismear,"%i"); /*since we follow the changing of ismear*/
	VASP_REG_VAL(sigma,"%lf");
	VASP_REG_TEXT(fermwe);
	VASP_REG_TEXT(fermdo);
//	VASP_REG_VAL(kspacing,"%lf"); /*since we follow the changing of kspacing*/
	//kgamma already sync
	//lreal already sync
	VASP_REG_TEXT(ropt);
	//addgrid already sync
	VASP_REG_VAL(lmaxmix,"%i");
	VASP_REG_VAL(lmaxpaw,"%i");
	VASP_REG_VAL(istart,"%i");
	VASP_REG_VAL(icharg,"%i");
        //iniwav already sync
	//ispin already sync
	//non_collinear already sync
	VASP_REG_TEXT(magmom);
	VASP_REG_VAL(nupdown,"%lf");
	//lsorbit already sync
	VASP_REG_TEXT(saxis);
	//gga_compat already sync
	//gga already sync
	//voskown already sync
	//lasph already sync
	VASP_REG_VAL(lmaxtau,"%i");
	//lmixtau already sync
	//lmetagga already sync
	//mgga already sync
	VASP_REG_TEXT(cmbj);
	VASP_REG_VAL(cmbja,"%lf");
	VASP_REG_VAL(cmbjb,"%lf");
	//ldau already sync
	//ldau_type already sync
	//ldau_output already sync
	VASP_REG_TEXT(ldaul);
	VASP_REG_TEXT(ldauu);
	VASP_REG_TEXT(ldauj);
	VASP_REG_VAL(nelm,"%i");
	VASP_REG_VAL(nelmdl,"%i");
	VASP_REG_VAL(nelmin,"%i");
	VASP_REG_VAL(amix,"%lf");
	VASP_REG_VAL(bmix,"%lf");
	VASP_REG_VAL(amin,"%lf");
	VASP_REG_VAL(amix_mag,"%lf");
	VASP_REG_VAL(bmix_mag,"%lf");
	//imix already sync
	VASP_REG_VAL(maxmix,"%i");
	VASP_REG_VAL(wc,"%lf");
	//inimix already sync
	//mixpre already sync
	//ldipol already sync
	//lmono already sync
	//idipol already sync
	VASP_REG_VAL(epsilon,"%lf");
	VASP_REG_TEXT(dipol);
	VASP_REG_VAL(efield,"%lf");
	//auto_grid already sync
	VASP_REG_VAL(ngx,"%i");
	VASP_REG_VAL(ngy,"%i");
	VASP_REG_VAL(ngz,"%i");
	VASP_REG_VAL(ngxf,"%i");
	VASP_REG_VAL(ngyf,"%i");
	VASP_REG_VAL(ngzf,"%i");
	VASP_REG_VAL(nsw,"%i");
	//ibrion already sync
	//isif already sync
	VASP_REG_VAL(pstress,"%lf");
	VASP_REG_VAL(ediffg,"%lf");
	VASP_REG_VAL(nfree,"%i");
	VASP_REG_VAL(potim,"%lf");
	VASP_REG_VAL(tebeg,"%lf");
	VASP_REG_VAL(teend,"%lf");
	VASP_REG_VAL(smass,"%lf");
	VASP_REG_VAL(nblock,"%i");
	VASP_REG_VAL(kblock,"%i");
	VASP_REG_VAL(npaco,"%i");
	VASP_REG_VAL(apaco,"%lf");
	VASP_REG_VAL(isym,"%i");
	VASP_REG_VAL(sym_prec,"%lE");
	//poscar_sd already sync
	//poscar_free already sync
	//poscar_direct already sync
	VASP_REG_VAL(poscar_a0,"%lf");
	VASP_REG_VAL(poscar_ux,"%lf");
	VASP_REG_VAL(poscar_uy,"%lf");
	VASP_REG_VAL(poscar_uz,"%lf");
	VASP_REG_VAL(poscar_vx,"%lf");
	VASP_REG_VAL(poscar_vy,"%lf");
	VASP_REG_VAL(poscar_vz,"%lf");
	VASP_REG_VAL(poscar_wx,"%lf");
	VASP_REG_VAL(poscar_wy,"%lf");
	VASP_REG_VAL(poscar_wz,"%lf");
	VASP_REG_VAL(poscar_index,"%i");
	/*PROCESS POSCAR*/
	vasp_poscar_sync();
	/*no need to sync individual poscar entry*/
	/*PROCESS KPOINTS*/
        VASP_REG_VAL(kpoints_nkpts,"%i");
	VASP_REG_VAL(kpoints_kx,"%lf");
	VASP_REG_VAL(kpoints_ky,"%lf");
	VASP_REG_VAL(kpoints_kz,"%lf");
	VASP_REG_VAL(kpoints_sx,"%lf");
	VASP_REG_VAL(kpoints_sy,"%lf");
	VASP_REG_VAL(kpoints_sz,"%lf");
        VASP_REG_VAL(kpoints_i,"%i");
        VASP_REG_VAL(kpoints_x,"%lf");
        VASP_REG_VAL(kpoints_y,"%lf");
        VASP_REG_VAL(kpoints_z,"%lf");
        VASP_REG_VAL(kpoints_w,"%lf");
	VASP_REG_VAL(tetra_total,"%i");
	VASP_REG_VAL(tetra_volume,"%lf");
	/*no need to sync individual tetrahedron entry*/
	/*PERF*/
	VASP_REG_VAL(ncore,"%lf");
	VASP_REG_VAL(kpar,"%lf");
	VASP_REG_VAL(job_nproc,"%lf");
}
/***************************************************/
/* convert vasp structure to POSCAR (not gui_free) */
/***************************************************/
void vasp_calc_to_poscar(FILE *output,vasp_calc_struct calc){
	gint index=gtk_combo_box_get_active(GTK_COMBO_BOX(vasp_gui.poscar_atoms));
	gchar *text;
	gdouble x,y,z;
	gchar tx,ty,tz;
	gint idx;
	/* generate a VASP5 POSCAR file */
        if(calc.name!=NULL) fprintf(output,"%s ",calc.name);
	else fprintf(output,"UNNAMED MODEL ");
	fprintf(output,"! GENERATED BY GDIS %4.2f.%d (C) %d\n",VERSION,PATCH,YEAR);
	fprintf(output,"%lf\n",calc.poscar_a0);
	fprintf(output,"   %12.8lf   %12.8lf   %12.8lf\n",calc.poscar_ux,calc.poscar_uy,calc.poscar_uz);
	fprintf(output,"   %12.8lf   %12.8lf   %12.8lf\n",calc.poscar_vx,calc.poscar_vy,calc.poscar_vz);
	fprintf(output,"   %12.8lf   %12.8lf   %12.8lf\n",calc.poscar_wx,calc.poscar_wy,calc.poscar_wz);
	/*SPECIES TAGS for VASP5 POSCAR FORMAT*/
	fprintf(output,"%s\n",vasp_calc.species_symbols);
	/*NUM ATOM PER SPECIES*/
	fprintf(output,"%s\n",vasp_calc.species_numbers);
	if(calc.poscar_sd) fprintf(output,"Selective dynamics\n");
	if(calc.poscar_direct) fprintf(output,"Direct\n");
	else fprintf(output,"Cartesian\n");
	idx=0;
        gtk_combo_box_set_active(GTK_COMBO_BOX(vasp_gui.poscar_atoms),idx);
        text=gtk_combo_box_text_get_active_text(GTK_COMBO_BOX_TEXT(vasp_gui.poscar_atoms));
	sscanf(text,"%*f %*f %*f %*c %*c %*c ! atom: %*i (%[^)])",&(vasp_calc.poscar_symbol[0]));
        while(g_ascii_strcasecmp(text,"ADD atom")!=0){
		/*get each line*/
		sscanf(text,"%lf %lf %lf %c %c %c %*s",&x,&y,&z,&tx,&ty,&tz);
		fprintf(output,"   % .8lf   % .8lf   % .8lf  %c   %c   %c\n",x,y,z,tx,ty,tz);
		idx++;
		gtk_combo_box_set_active(GTK_COMBO_BOX(vasp_gui.poscar_atoms),idx);
		text=gtk_combo_box_text_get_active_text(GTK_COMBO_BOX_TEXT(vasp_gui.poscar_atoms));
                }
	gtk_combo_box_set_active(GTK_COMBO_BOX(vasp_gui.poscar_atoms),index);
}
/****************************************************/
/* convert vasp structure to KPOINTS (not gui_free) */
/****************************************************/
void vasp_calc_to_kpoints(FILE *output,vasp_calc_struct calc){
	gint index=gtk_combo_box_get_active(GTK_COMBO_BOX(vasp_gui.kpoints_kpts));
	gchar *text;
	gdouble x,y,z,wt;
	gint a,b,c,d;
	gint idx;
	/* generate a VASP KPOINTS file */
	if(calc.name!=NULL) fprintf(output,"%s ",calc.name);
	else fprintf(output,"UNNAMED MODEL ");
	fprintf(output,"! GENERATED BY GDIS %4.2f.%d (C) %d\n",VERSION,PATCH,YEAR);
	/*nb of kpoints*/
	if((vasp_calc.kpoints_mode==VKP_MAN)||(vasp_calc.kpoints_mode==VKP_LINE)) fprintf(output,"%i\n",calc.kpoints_nkpts);
	else fprintf(output,"0\n");
	
	switch (vasp_calc.kpoints_mode){
	case VKP_LINE:
		fprintf(output,"Line-mode\n");
	case VKP_BASIS:
        case VKP_MAN:
		if(vasp_calc.kpoints_cart) fprintf(output,"Cartesian\n");
		else fprintf(output,"Reciprocal\n");
		break;
        case VKP_GAMMA:
		fprintf(output,"Gamma\n");
		fprintf(output,"%i %i %i\n",(gint)vasp_calc.kpoints_kx,(gint)vasp_calc.kpoints_ky,(gint)vasp_calc.kpoints_kz);
		if((vasp_calc.kpoints_sx!=0.0)&&(vasp_calc.kpoints_sy!=0.0)&&(vasp_calc.kpoints_sz!=0.0))
			fprintf(output,"%lf %lf %lf\n",vasp_calc.kpoints_sx,vasp_calc.kpoints_sy,vasp_calc.kpoints_sz);
		return;
        case VKP_MP:
		fprintf(output,"Monkhorst-Pack\n");
		fprintf(output,"%i %i %i\n",(gint)vasp_calc.kpoints_kx,(gint)vasp_calc.kpoints_ky,(gint)vasp_calc.kpoints_kz);
		if((vasp_calc.kpoints_sx!=0.0)&&(vasp_calc.kpoints_sy!=0.0)&&(vasp_calc.kpoints_sz!=0.0))
			fprintf(output,"%lf %lf %lf\n",vasp_calc.kpoints_sx,vasp_calc.kpoints_sy,vasp_calc.kpoints_sz);
		return;
        case VKP_AUTO:
		fprintf(output,"Auto\n");
		fprintf(output,"%i",(gint)vasp_calc.kpoints_kx);
		return;
        }
	/*print all coordinates*/
	idx=0;
	gtk_combo_box_set_active(GTK_COMBO_BOX(vasp_gui.kpoints_kpts),idx);
	text=gtk_combo_box_text_get_active_text(GTK_COMBO_BOX_TEXT(vasp_gui.kpoints_kpts));
        while(g_ascii_strcasecmp(text,"ADD kpoint")!=0){
		if(vasp_calc.kpoints_mode==VKP_MAN) {
			sscanf(text,"%lf %lf %lf %lf ! kpoint: %*i",&x,&y,&z,&wt);
			fprintf(output,"%lf %lf %lf %lf\n",x,y,z,wt);
		} else {
			sscanf(text,"%lf %lf %lf ! kpoint: %*i",&x,&y,&z);
			fprintf(output,"%lf %lf %lf\n",x,y,z);
		}
		idx++;
		gtk_combo_box_set_active(GTK_COMBO_BOX(vasp_gui.kpoints_kpts),idx);
		text=gtk_combo_box_text_get_active_text(GTK_COMBO_BOX_TEXT(vasp_gui.kpoints_kpts));
	}
	/*return to previous postition*/
        gtk_combo_box_set_active(GTK_COMBO_BOX(vasp_gui.kpoints_kpts),index);
	/*print all tetrahedron (if any)*/
	if(!vasp_calc.kpoints_tetra) return;
	fprintf(output,"Tetrahedron\n");
	fprintf(output,"%i %lf\n",vasp_calc.tetra_total,vasp_calc.tetra_volume);
	index=gtk_combo_box_get_active(GTK_COMBO_BOX(vasp_gui.tetra));
	idx=0;
	gtk_combo_box_set_active(GTK_COMBO_BOX(vasp_gui.tetra),idx);
	text=gtk_combo_box_text_get_active_text(GTK_COMBO_BOX_TEXT(vasp_gui.tetra));
	while(g_ascii_strcasecmp(text,"ADD tetrahedron")!=0){
		sscanf(text,"%lf %i %i %i %i ! tetrahedron: %*i",&wt,&a,&b,&c,&d);
		fprintf(output,"%lf %i %i %i %i\n",wt,a,b,c,d);
		idx++;
		gtk_combo_box_set_active(GTK_COMBO_BOX(vasp_gui.tetra),idx);
		text=gtk_combo_box_text_get_active_text(GTK_COMBO_BOX_TEXT(vasp_gui.tetra));
	}
	/*return to previous postition*/
	gtk_combo_box_set_active(GTK_COMBO_BOX(vasp_gui.tetra),index);
}
/**/
/**/
/**/
void vasp_calc_to_potcar(FILE *output,vasp_calc_struct calc){
	FILE *org_potcar;
	gchar *line;
	/*saving of the POTCAR is for now a mere copy of provided POTCAR*/
	org_potcar=fopen(calc.potcar_file,"r");
	if(!org_potcar) {
		fprintf(stderr,"#ERR: can't open file %s for READ!\n",calc.potcar_file);
		return;
	}
	line = file_read_line(org_potcar);
	while(line){
		fprintf(output,"%s",line);
		g_free(line);
		line = file_read_line(org_potcar);
	}
	fclose(org_potcar);
	/*that's all folks*/
}
/***************************/
/* save current parameters */
/***************************/
void save_vasp_calc(){
#define DEBUG_VASP_SAVE 0
	FILE *save_fp;
	gchar *filename;
	/*1-synchronize*/
	vasp_gui_sync();
	vasp_poscar_sync();
	/*2-output*/
	//INCAR
	filename=g_strdup_printf("%s/INCAR",vasp_calc.job_path);
	save_fp=fopen(filename,"w");
	if(!save_fp){
		fprintf(stderr,"#ERR: can't open file INCAR for WRITE!\n");
		return;
	}
	vasp_calc_to_incar(save_fp,vasp_calc);
	fclose(save_fp);
	g_free(filename);
	//POSCAR
	filename=g_strdup_printf("%s/POSCAR",vasp_calc.job_path);
	save_fp=fopen(filename,"w");
	if(!save_fp){
		fprintf(stderr,"#ERR: can't open file POSCAR for WRITE!\n");
		return;
	}
	vasp_calc_to_poscar(save_fp,vasp_calc);
	fclose(save_fp);
	g_free(filename);
	//KPOINTS
	filename=g_strdup_printf("%s/KPOINTS",vasp_calc.job_path);
	save_fp=fopen(filename,"w");
	if(!save_fp){
		fprintf(stderr,"#ERR: can't open file KPOINTS for WRITE!\n");
		return;
	}
	vasp_calc_to_kpoints(save_fp,vasp_calc);
	fclose(save_fp);
	g_free(filename);
	//POTCAR
	filename=g_strdup_printf("%s/POTCAR",vasp_calc.job_path);
	save_fp=fopen(filename,"w");
	if(!save_fp){
		fprintf(stderr,"#ERR: can't open file KPOINTS for WRITE!\n");
		return;
	}
	vasp_calc_to_potcar(save_fp,vasp_calc);
	fclose(save_fp);
	g_free(filename);
	/*debug print*/
#if DEBUG_VASP_SAVE
	vasp_calc_to_incar(stdout,vasp_calc);
	vasp_calc_to_poscar(stdout,vasp_calc);
	vasp_calc_to_kpoints(stdout,vasp_calc);

#endif
}
/*****************************************/
/* Execute or enqueue a vasp calculation */
/*****************************************/
void run_vasp_exec(vasp_calc_struct *calc,struct task_pak *task){
	/* Execute a vasp task TODO distant job */
	gchar *cmd;
	gchar *cwd;/*for push/pull*/

#if __WIN32
/*	At present running vasp in __WIN32 environment is impossible.
 * 	However,  it should be possible to launch a (remote) job on a 
 * 	distant server. See TODO */
	fprintf(stderr,"VASP calculation can't be done in this environment.\n");return;
#else 
	/*1CPU*/
	if((*calc).job_nproc<2) cmd = g_strdup_printf("%s > vasp.log",(*calc).job_vasp_exe);
	else cmd = g_strdup_printf("%s -np %i %s > vasp.log",(*calc).job_mpirun,(gint)(*calc).job_nproc,(*calc).job_vasp_exe);
#endif
	cwd=sysenv.cwd;/*push*/
	sysenv.cwd=g_strdup_printf("%s",(*calc).job_path);
	task_sync(cmd);
	g_free(sysenv.cwd);
	sysenv.cwd=cwd;/*pull*/
	g_free(cmd);
}
/********************************************************/
/* hide spinner event (will load result on have_result) */
/********************************************************/
void vasp_spin_hide(){
if(!vasp_gui.have_result) return;
else {
	gchar *filename=g_strdup_printf("%s/vasprun.xml",gtk_entry_get_text(GTK_ENTRY(vasp_gui.job_path)));

/* FIXME: Why I can't see the resulting result model until I remove original one?*/

	file_load(filename,vasp_gui.result_model);
	model_prep(vasp_gui.result_model);
	tree_model_add(vasp_gui.result_model);
	tree_select_model(vasp_gui.result_model);
	tree_select_active();
	tree_model_refresh(vasp_gui.result_model);
//	gui_model_select(vasp_gui.result_model);
//	gui_active_refresh();
}
vasp_gui.have_result=FALSE;
}
/*****************************************************/
/* cleanup task: trigger spinner (which load result) */
/*****************************************************/
void cleanup_vasp_exec(struct vasp_calc_gui *gui_vasp){
	/*VASP process has exit, try to load the result?*/
	(*gui_vasp).have_result=TRUE;
	gtk_spinner_stop(GTK_SPINNER((*gui_vasp).spinner));
	gtk_widget_hide((*gui_vasp).spinner);
}
/******************************/
/* Enqueue a vasp calculation */
/******************************/
void exec_vasp_calc(){
	struct model_pak *data = sysenv.active_model;
	/*this will sync then enqueue a VASP calculation*/
	if(vasp_calc.potcar_file==NULL) return;/*not ready*/
	save_vasp_calc();/*sync and save all file*/
	/*launch vasp in a task?*/
	gtk_widget_show(vasp_gui.spinner);
	gtk_spinner_start(GTK_SPINNER(vasp_gui.spinner));
	task_new("VASP", &run_vasp_exec,&vasp_calc,&cleanup_vasp_exec,&vasp_gui,data);
}
/*********************************/
/* Cleanup calculation structure */
/*********************************/
void vasp_cleanup(struct model_pak *model){

}
/********************************/
/* Quit vasp calculation dialog */
/********************************/
void quit_vasp_calc(GtkWidget *w, gpointer data){
	struct model_pak *model=data;
	g_assert(model != NULL);
	/* first cleanup everything */
	vasp_cleanup(model);
	/* Then bailout */
	dialog_destroy(w,model);
}
void gui_vasp_dialog(void){
/*launch the main interface*/
/*TODO: use distant connections*/
	gchar *title;
	gpointer dialog;
	GtkWidget *frame, *vbox, *hbox;
	GtkWidget *notebook, *page, *button, *label;
	GtkWidget *combo, *separator, *spinner;
	/* special */
	GList *list;
	GSList *list2;
	struct model_pak *data;
	struct core_pak *core;
	gint idx;

/* checks */
	data = sysenv.active_model;
	if (!data) return;
	/*TODO Only one active at a time*/
/* do we have a vasprun.xml model?  */
	if (data->id == VASP) {
		vasp_gui.have_xml=TRUE;
		gui_vasp_init();
		vasprun_update(data->filename,&vasp_calc);/*initialize according to vasprun.xml*/
	} else {
		vasp_gui.have_xml=FALSE;/*initialize with default values*/
		gui_vasp_init();
		vasprun_update(NULL,&vasp_calc);
	}
/* dialog setup */
	title = g_strdup_printf("VASP: %s", data->basename);
	dialog = dialog_request(CVASP, title, NULL, vasp_cleanup, data);
	g_free(title);
	vasp_gui.window = dialog_window(dialog);
/* --- Outside of notebook */
	frame = gtk_frame_new(NULL);
	gtk_box_pack_start(GTK_BOX(GTK_DIALOG(vasp_gui.window)->vbox), frame, FALSE, FALSE, 0);
	vbox = gtk_vbox_new(FALSE, _SPACE);
	gtk_container_add(GTK_CONTAINER(frame), vbox);
/* --- MODEL NAME */
	hbox = gtk_hbox_new(FALSE, 0);
	gtk_box_pack_start(GTK_BOX(vbox), hbox, FALSE, FALSE, _SPACE);
	label = gtk_label_new(g_strdup_printf("MODEL NAME"));
	gtk_box_pack_start(GTK_BOX(hbox), label, FALSE, FALSE, 0);
	VASP_TEXT_ENTRY(vasp_gui.name,g_strdup_printf("v_%s",data->basename),TRUE);
/* --- CONNECTED XML */
	hbox = gtk_hbox_new(FALSE, 0);
        gtk_box_pack_start(GTK_BOX(vbox), hbox, FALSE, FALSE, _SPACE);
        label = gtk_label_new(g_strdup_printf("CONNECTED XML"));
        gtk_box_pack_start(GTK_BOX(hbox), label, FALSE, FALSE, 0);
	VASP_TEXT_ENTRY(vasp_gui.file_entry,data->filename,TRUE);
	if(data->id!=VASP) gtk_entry_set_text(GTK_ENTRY(vasp_gui.file_entry),"");
	button=gtk_button_new_from_stock(GTK_STOCK_OPEN);
	gtk_box_pack_end(GTK_BOX(hbox), button, FALSE, FALSE, 0);
	g_signal_connect(GTK_OBJECT(button),"clicked",GTK_SIGNAL_FUNC(load_vasprun_dialog),NULL);
/* notebook frame */
	frame = gtk_frame_new(NULL);
	gtk_box_pack_start(GTK_BOX(GTK_DIALOG(vasp_gui.window)->vbox), frame, FALSE, FALSE, 0);
	gtk_container_set_border_width(GTK_CONTAINER(frame), _SPACE);
/* create notebook */
	notebook = gtk_notebook_new();
	gtk_notebook_set_tab_pos(GTK_NOTEBOOK(notebook), GTK_POS_TOP);
	gtk_container_add(GTK_CONTAINER(frame), notebook);
	gtk_notebook_set_show_border(GTK_NOTEBOOK(notebook), TRUE);

/*----------------------*/
/* page 0 -> SIMPLIFIED */
/*----------------------*/
        page = gtk_vbox_new(FALSE, _SPACE);
        label = gtk_label_new("SIMPLIFIED");
        gtk_notebook_append_page(GTK_NOTEBOOK(notebook),page,label);
/* --- DISCLAMER */
        frame = gtk_frame_new("DISCLAMER");
        gtk_box_pack_start(GTK_BOX(page),frame,FALSE,FALSE,0);
/* create a vbox in the frame */
        vbox = gtk_vbox_new(TRUE, 0);
        gtk_container_add(GTK_CONTAINER(frame), vbox);
/* text */
        VASP_NEW_LINE();
	VASP_NEW_SEPARATOR();
	spinner=gtk_spinner_new();
	gtk_box_pack_start(GTK_BOX(hbox),spinner, FALSE, FALSE, 0);
	VASP_NEW_LABEL("UNDER CONSTRUCTION");
	VASP_NEW_SEPARATOR();
	VASP_NEW_LINE();
	VASP_NEW_SEPARATOR();
	VASP_NEW_LINE();
	VASP_NEW_LABEL("Please report bugs at https://github.com/ovhpa/gdis/tree/gui_vasp");
	VASP_NEW_SEPARATOR();
	VASP_NEW_LINE();
	VASP_NEW_SEPARATOR();
/* --- end frame */
        hbox = gtk_hbox_new(FALSE, 0);
        gtk_box_pack_end(GTK_BOX(page),hbox,FALSE,FALSE,0);

/*-----------------------*/
/* page 1 -> CONVERGENCE */
/*-----------------------*/
	page = gtk_vbox_new(FALSE, _SPACE);
	label = gtk_label_new("CONVERGENCE");
	gtk_notebook_append_page(GTK_NOTEBOOK(notebook),page,label);
/* --- general */
	frame = gtk_frame_new("General");
	gtk_box_pack_start(GTK_BOX(page),frame,FALSE,FALSE,0);
/* create a vbox in the frame */
	vbox = gtk_vbox_new(TRUE, 0);
	gtk_container_add(GTK_CONTAINER(frame), vbox);
/* 1st line */
	VASP_NEW_LINE();
	VASP_NEW_LABEL("PREC=");
	/*TODO: folding */
	list = NULL;
	list = g_list_append(list,"Normal");
	list = g_list_append(list,"Single");
	list = g_list_append(list,"Accurate");
	list = g_list_append(list,"High");
	list = g_list_append(list,"Medium");
	list = g_list_append(list,"Low");
	VASP_REG_COMBO(combo,"Normal",vasp_prec_selected);
	VASP_NEW_LABEL("ENCUT=");
	VASP_ENTRY(vasp_gui.encut,vasp_calc.encut,"%.2lf",6);
	VASP_NEW_LABEL("eV ");
	VASP_NEW_LABEL("ENAUG=");
	VASP_ENTRY(vasp_gui.enaug,vasp_calc.enaug,"%.2lf",6);
	VASP_NEW_LABEL("eV ");
	VASP_NEW_SEPARATOR();
	VASP_NEW_LABEL("EDIFF=");
	VASP_ENTRY(vasp_gui.ediff,vasp_calc.ediff,"%.2lE",8);
	VASP_NEW_LABEL("eV");
/* 2nd line */
	VASP_NEW_LINE();
	gui_direct_check("USE PREC",&(vasp_calc.use_prec),prec_toggle,NULL,hbox);
	prec_toggle(NULL,NULL);/*initialize TODO: move to end*/
	VASP_NEW_LABEL("ALGO=");
        list = NULL;
        list = g_list_append(list,"NORMAL");
	list = g_list_append(list,"USE_IALGO");
        list = g_list_append(list,"VERYFAST");
        list = g_list_append(list,"FAST");
        list = g_list_append(list,"CONJUGATE");
	list = g_list_append(list,"ALL");
	list = g_list_append(list,"DAMPED");
	list = g_list_append(list,"SUBROT");
	list = g_list_append(list,"EIGENVAL");
	list = g_list_append(list,"NONE");
	list = g_list_append(list,"NOTHING");
	list = g_list_append(list,"EXACT");
        list = g_list_append(list,"DIAG");
	VASP_REG_COMBO(combo,"NORMAL",vasp_algo_selected);/*no need to register widget*/
	VASP_NEW_LABEL("IALGO=");
	list = NULL;
	list = g_list_append(list,"2:FIXED ORB/1E");
	list = g_list_append(list,"3:FIXED ORB");
	list = g_list_append(list,"4:SUBROT ONLY");
	list = g_list_append(list,"4:SUBROT ONLY");
	list = g_list_append(list,"5:STEEPEST");
	list = g_list_append(list,"6:CONJUGUATE GRADIENTS");
	list = g_list_append(list,"7:PRECOND. STEEPEST");
	list = g_list_append(list,"8:PRECOND. CG");
	list = g_list_append(list,"38:KOSUGI");
	list = g_list_append(list,"44:RMM STEEPEST");
	list = g_list_append(list,"46:RMM PRECOND.");
	list = g_list_append(list,"48:PRECOND. RMM");
	list = g_list_append(list,"53:VAR DAMPED MD");
	list = g_list_append(list,"54:VAR DAMPED QUENCH MD");
	list = g_list_append(list,"58:VAR PRECOND. CG");
	list = g_list_append(list,"90:EXACT");
	VASP_REG_COMBO(vasp_gui.ialgo,"38:KOSUGI",vasp_ialgo_selected);
	if(vasp_calc.algo!=VA_IALGO) gtk_widget_set_sensitive(vasp_gui.ialgo,FALSE);
	else gtk_widget_set_sensitive(vasp_gui.ialgo,TRUE);
	gui_direct_check("LDIAG",&(vasp_calc.ldiag),NULL,NULL,hbox);/*not calling anything*/
/* 3rd line */
	VASP_NEW_LINE();
	gui_direct_check("AUTO",&(vasp_calc.auto_elec),elec_toggle,NULL,hbox);
	VASP_NEW_LABEL("NSIM=");
	VASP_ENTRY(vasp_gui.nsim,vasp_calc.nsim,"%i",3);
	VASP_NEW_LABEL("TIME=");
	VASP_ENTRY(vasp_gui.vtime,vasp_calc.vtime,"%.4lf",6);
	VASP_NEW_LABEL("NBANDS=");
	VASP_ENTRY(vasp_gui.nbands,vasp_calc.nbands,"%i",6);
	VASP_NEW_LABEL("NELECT=");
	VASP_ENTRY(vasp_gui.nelect,vasp_calc.nelect,"%.4lf",6);
	VASP_NEW_LABEL("INIWAVPR=");
	VASP_ENTRY(vasp_gui.iwavpr,vasp_calc.iwavpr,"%i",3);/*<- should be a combo?*/
	elec_toggle(NULL,NULL);/*initialize TODO: move to end*/
/* --- end frame */
	hbox = gtk_hbox_new(FALSE, 0);
	gtk_box_pack_end(GTK_BOX(page),hbox,FALSE,FALSE,0);
/* --- unnamed frame */
        frame = gtk_frame_new(NULL);
        gtk_box_pack_start(GTK_BOX(page),frame,FALSE,FALSE,0);
/* create a vbox in the frame */
        vbox = gtk_vbox_new(TRUE, 0);
        gtk_container_add(GTK_CONTAINER(frame), vbox);
/* 1st line */
        /* TODO: this whole thing should move to another page */
        VASP_NEW_LINE();
        VASP_NEW_LABEL("ISTART=");
        VASP_ENTRY(vasp_gui.istart,vasp_calc.istart,"%i",3);/*<- should be a combo?*/
        VASP_NEW_LABEL("ICHARG=");
        VASP_ENTRY(vasp_gui.icharg,vasp_calc.icharg,"%i",3);/*<- should be a combo?*/
        gui_direct_check("INIWAV",&(vasp_calc.iniwav),NULL,NULL,hbox);/*not calling anything*/
        VASP_NEW_SEPARATOR();
        label = gtk_label_new(" ");/*TODO trick to be removed*/
        gtk_box_pack_end(GTK_BOX(hbox), label, FALSE, FALSE, 0);
/* --- end frame */
        hbox = gtk_hbox_new(FALSE, 0);
        gtk_box_pack_end(GTK_BOX(page),hbox,FALSE,FALSE,0);
/* --- mixing */
        frame = gtk_frame_new("Mixing");
        gtk_box_pack_start(GTK_BOX(page),frame,FALSE,FALSE,0);
/* create a vbox in the frame */
        vbox = gtk_vbox_new(TRUE, 0);
        gtk_container_add(GTK_CONTAINER(frame), vbox);
/* 1st line */
	VASP_NEW_LINE();
	VASP_NEW_LABEL("NELM=");
	VASP_ENTRY(vasp_gui.nelm,vasp_calc.nelm,"%i",3);
	VASP_NEW_LABEL("NELMDL=");
	VASP_ENTRY(vasp_gui.nelmdl,vasp_calc.nelmdl,"%i",3);
	VASP_NEW_LABEL("NELMIN=");
	VASP_ENTRY(vasp_gui.nelmin,vasp_calc.nelmin,"%i",3);
	gui_direct_check("AUTO MIXER",&(vasp_calc.auto_mixer),mixer_toggle,NULL,hbox);
	VASP_NEW_LABEL("IMIX=");
        list = NULL;
        list = g_list_append(list,"0:NO MIXING");
        list = g_list_append(list,"1:KERKER");
        list = g_list_append(list,"2:TCHEBYCHEV");
        list = g_list_append(list,"4:BROYDEN");
	VASP_REG_COMBO(vasp_gui.mixer,"4:BROYDEN",vasp_mixer_selected);
/* 2nd line */
	VASP_NEW_LINE();
	VASP_NEW_LABEL("AMIX=");
	VASP_ENTRY(vasp_gui.amix,vasp_calc.amix,"%.4lf",6);
	VASP_NEW_LABEL("BMIX=");
	VASP_ENTRY(vasp_gui.bmix,vasp_calc.bmix,"%.4lf",6);
	VASP_NEW_LABEL("AMIN=");
	VASP_ENTRY(vasp_gui.amin,vasp_calc.amin,"%.4lf",6);
        VASP_NEW_LABEL("AMIX_MAG=");
	VASP_ENTRY(vasp_gui.amix_mag,vasp_calc.amix_mag,"%.4lf",6);
        VASP_NEW_LABEL("BMIX_MAG=");
	VASP_ENTRY(vasp_gui.bmix_mag,vasp_calc.bmix_mag,"%.4lf",6);
/* 3rd line */
	VASP_NEW_LINE();
	VASP_NEW_LABEL("MAXMIX=");
	VASP_ENTRY(vasp_gui.maxmix,vasp_calc.maxmix,"%i",3);
	VASP_NEW_LABEL("WC=");
	VASP_ENTRY(vasp_gui.wc,vasp_calc.wc,"%.1lf",8);
	VASP_NEW_SEPARATOR();
	VASP_NEW_LABEL("INIMIX=");
        list = NULL;
        list = g_list_append(list,"0:LINEAR");
        list = g_list_append(list,"1:KERKER");
	VASP_REG_COMBO(vasp_gui.inimix,"1:KERKER",vasp_inimix_selected);
	VASP_NEW_LABEL("MIXPRE=");
        list = NULL;
        list = g_list_append(list,"0:NONE");
        list = g_list_append(list,"1:F=20 INVERSE KERKER");
        list = g_list_append(list,"2:F=200 INVERSE KERKER");
	VASP_REG_COMBO(vasp_gui.mixpre,"1:F=20 INVERSE KERKER",vasp_mixpre_selected);
	mixer_toggle(NULL,NULL);/*initialize TODO: move to end*/
/* --- end frame */
        hbox = gtk_hbox_new(FALSE, 0);
        gtk_box_pack_end(GTK_BOX(page),hbox,FALSE,FALSE,0);
/* --- projector */
        frame = gtk_frame_new("Projector");
        gtk_box_pack_start(GTK_BOX(page),frame,FALSE,FALSE,0);
/* create a vbox in the frame */
        vbox = gtk_vbox_new(TRUE, 0);
        gtk_container_add(GTK_CONTAINER(frame), vbox);
/* 1st line */
        VASP_NEW_LINE();
        VASP_NEW_LABEL("LREAL=");
        list = NULL;
        list = g_list_append(list,"Auto");
        list = g_list_append(list,"On");
        list = g_list_append(list,"TRUE");
        list = g_list_append(list,"FALSE");
        VASP_REG_COMBO(combo,"FALSE",vasp_lreal_selected);
        VASP_NEW_LABEL("ROPT=");
        VASP_TEXT_ENTRY(vasp_gui.ropt,vasp_calc.ropt,TRUE);
/* 2nd line */
        VASP_NEW_LINE();
        VASP_NEW_SEPARATOR();
        gui_direct_check("ADDGRID",&(vasp_calc.addgrid),NULL,NULL,hbox);/*not calling anything*/
        VASP_NEW_LABEL("LMAXMIX=");
        VASP_ENTRY(vasp_gui.lmaxmix,vasp_calc.lmaxmix,"%i",3);
        VASP_NEW_LABEL("LMAXPAW=");
        VASP_ENTRY(vasp_gui.lmaxpaw,vasp_calc.lmaxpaw,"%i",4);
/* --- end frame */
        hbox = gtk_hbox_new(FALSE, 0);
        gtk_box_pack_end(GTK_BOX(page),hbox,FALSE,FALSE,0);
/* --- grids */
        frame = gtk_frame_new("Grids");
        gtk_box_pack_start(GTK_BOX(page),frame,FALSE,FALSE,0);
/* create a vbox in the frame */
        vbox = gtk_vbox_new(TRUE, 0);
        gtk_container_add(GTK_CONTAINER(frame), vbox);
/* 1st line */
        VASP_NEW_LINE();
	gui_direct_check("AUTO",&(vasp_calc.auto_grid),grid_toggle,NULL,hbox);
	VASP_NEW_SEPARATOR();
	VASP_NEW_LABEL("NGX=");
	VASP_ENTRY(vasp_gui.ngx,vasp_calc.ngx,"%i",4);
	VASP_NEW_LABEL("NGY=");
	VASP_ENTRY(vasp_gui.ngy,vasp_calc.ngy,"%i",4);
	VASP_NEW_LABEL("NGZ=");
	VASP_ENTRY(vasp_gui.ngz,vasp_calc.ngz,"%i",4);
/* 2nd line */
        VASP_NEW_LINE();
	VASP_NEW_SEPARATOR();
        VASP_NEW_LABEL("NGXF=");
        VASP_ENTRY(vasp_gui.ngxf,vasp_calc.ngxf,"%i",4);
        VASP_NEW_LABEL("NGYF=");
        VASP_ENTRY(vasp_gui.ngyf,vasp_calc.ngyf,"%i",4);
        VASP_NEW_LABEL("NGZF=");
        VASP_ENTRY(vasp_gui.ngzf,vasp_calc.ngzf,"%i",4);
	grid_toggle(NULL,NULL);/*initialize TODO: move to end*/
/* --- end frame */
        hbox = gtk_hbox_new(FALSE, 0);
        gtk_box_pack_end(GTK_BOX(page),hbox,FALSE,FALSE,0);
/*----------------------*/
/* page 2 -> ELECTRONIC */
/*----------------------*/
        page = gtk_vbox_new(FALSE, _SPACE);
        label = gtk_label_new("ELECTRONIC");
        gtk_notebook_append_page(GTK_NOTEBOOK(notebook),page,label);
/* --- smearing */
        frame = gtk_frame_new("Smearing");
        gtk_box_pack_start(GTK_BOX(page),frame,FALSE,FALSE,0);
/* create a vbox in the frame */
        vbox = gtk_vbox_new(TRUE, 0);
        gtk_container_add(GTK_CONTAINER(frame), vbox);
/* 1st line */
        VASP_NEW_LINE();
        VASP_NEW_LABEL("ISMEAR=");
        VASP_ENTRY(vasp_gui.ismear,vasp_calc.ismear,"%i",3);/*<- should be a combo?*/
	g_signal_connect(GTK_OBJECT(GTK_ENTRY(vasp_gui.ismear)),"changed",GTK_SIGNAL_FUNC(ismear_changed),data);
        VASP_NEW_LABEL("SIGMA=");
        VASP_ENTRY(vasp_gui.sigma,vasp_calc.sigma,"%.4lf",6);
        VASP_NEW_SEPARATOR();
        vasp_gui.kgamma=gui_direct_check("KGAMMA",&(vasp_calc.kgamma),NULL,NULL,hbox);/*not calling anything*/
        VASP_NEW_LABEL("KSPACING=");
        VASP_ENTRY(vasp_gui.kspacing,vasp_calc.kspacing,"%.4lf",6);
	g_signal_connect(GTK_OBJECT(GTK_ENTRY(vasp_gui.kspacing)),"changed",GTK_SIGNAL_FUNC(kspacing_changed),data);
/* 2nd line */
        VASP_NEW_LINE();
        VASP_NEW_LABEL("FERMWE=");
        VASP_TEXT_ENTRY(vasp_gui.fermwe,vasp_calc.fermwe,TRUE);
        VASP_NEW_LABEL("FERMDO=");
        VASP_TEXT_ENTRY(vasp_gui.fermdo,vasp_calc.fermdo,TRUE);
/* --- end frame */
        hbox = gtk_hbox_new(FALSE, 0);
        gtk_box_pack_end(GTK_BOX(page),hbox,FALSE,FALSE,0);
/* --- spin */
        frame = gtk_frame_new("Spin");
        gtk_box_pack_start(GTK_BOX(page),frame,FALSE,FALSE,0);
/* create a vbox in the frame */
        vbox = gtk_vbox_new(TRUE, 0);
        gtk_container_add(GTK_CONTAINER(frame), vbox);
/* 1st line */
        VASP_NEW_LINE();
        gui_direct_check("SPIN POLARIZED",&(vasp_calc.ispin),NULL,NULL,hbox);/*not calling anything*/
        vasp_gui.lnoncoll=gui_direct_check("NON COLLINEAR",&(vasp_calc.non_collinear),lnoncoll_toggle,NULL,hbox);
        VASP_NEW_LABEL("NUPDOWN=");
        VASP_ENTRY(vasp_gui.nupdown,vasp_calc.nupdown,"%.1lf",4);
        VASP_NEW_LABEL("MAGMOM=");
        VASP_TEXT_ENTRY(vasp_gui.magmom,vasp_calc.magmom,TRUE);
/* 2nd line */
        VASP_NEW_LINE();
        vasp_gui.lsorbit=gui_direct_check("LSORBIT",&(vasp_calc.lsorbit),lsorbit_toggle,NULL,hbox);
        VASP_NEW_LABEL("SAXIS=");
        VASP_TEXT_ENTRY(vasp_gui.saxis,vasp_calc.saxis,TRUE);
        gui_direct_check("GGA_COMPAT",&(vasp_calc.gga_compat),NULL,NULL,hbox);/*not calling anything*/
/* --- end frame */
        hbox = gtk_hbox_new(FALSE, 0);
        gtk_box_pack_end(GTK_BOX(page),hbox,FALSE,FALSE,0);
/* --- exchange correlation */
        frame = gtk_frame_new("Exchange Correlation");
        gtk_box_pack_start(GTK_BOX(page),frame,FALSE,FALSE,0);
/* create a vbox in the frame */
        vbox = gtk_vbox_new(TRUE, 0);
        gtk_container_add(GTK_CONTAINER(frame), vbox);
/* 1st line */
        VASP_NEW_LINE();
        VASP_NEW_LABEL("GGA=");
        list = NULL;
        list = g_list_append(list,"91:PW91");
        list = g_list_append(list,"PE:PBE");
        list = g_list_append(list,"RP:rPBE");
        list = g_list_append(list,"AM:AM05");
        list = g_list_append(list,"PS:PBEsol");
        VASP_REG_COMBO(combo,"PE:PBE",vasp_gga_selected);
        gui_direct_check("VOSKOWN",&(vasp_calc.voskown),NULL,NULL,hbox);/*not calling anything*/
        gui_direct_check("LASPH",&(vasp_calc.lasph),NULL,NULL,hbox);/*not calling anything*/
        VASP_NEW_LABEL("LMAXTAU=");
        VASP_ENTRY(vasp_gui.lmaxtau,vasp_calc.lmaxtau,"%i",3);
/* 2nd line */
        VASP_NEW_LINE();
        VASP_NEW_LABEL("META GGA=");
        list = NULL;
        list = g_list_append(list,"TPSS");
        list = g_list_append(list,"rTPSS");
        list = g_list_append(list,"M06L");
        list = g_list_append(list,"MBJ");
        VASP_REG_COMBO(vasp_gui.metagga,"MBJ",vasp_metagga_selected);
if(vasp_calc.lmetagga) gtk_widget_set_sensitive(vasp_gui.metagga,TRUE);
else gtk_widget_set_sensitive(vasp_gui.metagga,FALSE);
        gui_direct_check("LMETAGGA",&(vasp_calc.lmetagga),metagga_toggle,NULL,hbox);
        VASP_NEW_SEPARATOR();
        gui_direct_check("LMIXTAU",&(vasp_calc.lmixtau),NULL,NULL,hbox);/*not calling anything*/
        VASP_NEW_LINE();
        VASP_NEW_LABEL("CMBJ=");
        VASP_TEXT_ENTRY(vasp_gui.cmbj,vasp_calc.cmbj,TRUE);
        VASP_NEW_LABEL("CMBJA=");
        VASP_ENTRY(vasp_gui.cmbja,vasp_calc.cmbja,"%.4lf",6);
        VASP_NEW_LABEL("CMBJB=");
        VASP_ENTRY(vasp_gui.cmbjb,vasp_calc.cmbjb,"%.4lf",6);
        metagga_toggle(NULL,NULL);/*initialize TODO: move to end*/
/* 4th line */
        VASP_NEW_LINE();
        gui_direct_check("LDA+U",&(vasp_calc.ldau),ldau_toggle,NULL,hbox);
        VASP_NEW_LABEL("TYPE=");
        list = NULL;
        list = g_list_append(list,"1:LSDA+U Liechtenstein");
        list = g_list_append(list,"2:LSDA+U Dudarev");
        list = g_list_append(list,"4:LDA+U Liechtenstein");
        VASP_REG_COMBO(vasp_gui.ldau,"2:LSDA+U Dudarev",vasp_ldau_selected);
if(vasp_calc.ldau) gtk_widget_set_sensitive(vasp_gui.ldau,TRUE);
else gtk_widget_set_sensitive(vasp_gui.ldau,FALSE);
        VASP_NEW_LABEL("Print=");
        list = NULL;
        list = g_list_append(list,"0:Silent");
        list = g_list_append(list,"1:Occupancy matrix");
        list = g_list_append(list,"2:Full output");
        VASP_REG_COMBO(vasp_gui.ldau_print,"0:Silent",vasp_ldau_print_selected);
if(vasp_calc.ldau) gtk_widget_set_sensitive(vasp_gui.ldau_print,TRUE);
else gtk_widget_set_sensitive(vasp_gui.ldau_print,FALSE);
/* 5th line */
        VASP_NEW_LINE();
        VASP_NEW_LABEL("LDAU-L=");
        VASP_TEXT_ENTRY(vasp_gui.ldaul,vasp_calc.ldaul,FALSE);
        VASP_NEW_LABEL("LDAU-U=");
        VASP_TEXT_ENTRY(vasp_gui.ldauu,vasp_calc.ldauu,FALSE);
        VASP_NEW_LABEL("LDAU-J=");
        VASP_TEXT_ENTRY(vasp_gui.ldauj,vasp_calc.ldauj,FALSE);
/* --- end frame */
        hbox = gtk_hbox_new(FALSE, 0);
        gtk_box_pack_end(GTK_BOX(page),hbox,FALSE,FALSE,0);
/* --- dipole */
        frame = gtk_frame_new("Dipole");
        gtk_box_pack_start(GTK_BOX(page),frame,FALSE,FALSE,0);
/* create a vbox in the frame */
        vbox = gtk_vbox_new(TRUE, 0);
        gtk_container_add(GTK_CONTAINER(frame), vbox);
/* 1st line */
        VASP_NEW_LINE();
	VASP_NEW_LABEL("IDIPOL=");
	list = NULL;
	list = g_list_append(list,"0:no calcul");
	list = g_list_append(list,"1:u axis");
	list = g_list_append(list,"2:v axis");
	list = g_list_append(list,"3:w axis");
	list = g_list_append(list,"4:all axis");
	VASP_REG_COMBO(vasp_gui.idipol,"0:no calcul",vasp_idipol_selected);
	VASP_NEW_SEPARATOR();
	vasp_gui.ldipol=gui_direct_check("LDIPOL",&(vasp_calc.ldipol),NULL,NULL,hbox);/*not calling anything*/
	vasp_gui.lmono=gui_direct_check("LMONO",&(vasp_calc.lmono),NULL,NULL,hbox);/*not calling anything*/
/* 2nd line */
        VASP_NEW_LINE();
        VASP_NEW_LABEL("EPSILON=");
        VASP_ENTRY(vasp_gui.epsilon,vasp_calc.epsilon,"%.4lf",6);
        VASP_NEW_LABEL("EFIELD=");
        VASP_ENTRY(vasp_gui.efield,vasp_calc.efield,"%.4lf",6);
	VASP_NEW_SEPARATOR();
	VASP_NEW_LABEL("DIPOL=");
	VASP_TEXT_ENTRY(vasp_gui.dipol,vasp_calc.dipol,TRUE);
	/*setup sensitive widget*/
if(vasp_calc.idipol!=VID_0) {
        gtk_widget_set_sensitive(vasp_gui.ldipol,TRUE);
        gtk_widget_set_sensitive(vasp_gui.lmono,TRUE);
        gtk_widget_set_sensitive(vasp_gui.dipol,TRUE);
        gtk_widget_set_sensitive(vasp_gui.epsilon,TRUE);
        if(vasp_calc.idipol!=VID_4) gtk_widget_set_sensitive(vasp_gui.efield,TRUE);
}else{
        gtk_widget_set_sensitive(vasp_gui.ldipol,FALSE);
        gtk_widget_set_sensitive(vasp_gui.lmono,FALSE);
        gtk_widget_set_sensitive(vasp_gui.dipol,FALSE);
        gtk_widget_set_sensitive(vasp_gui.epsilon,FALSE);
        gtk_widget_set_sensitive(vasp_gui.efield,FALSE);
}
/* --- end frame */
        hbox = gtk_hbox_new(FALSE, 0);
        gtk_box_pack_end(GTK_BOX(page),hbox,FALSE,FALSE,0);

/*-----------------*/
/* page 3 -> IONIC */
/*-----------------*/
        page = gtk_vbox_new(FALSE, _SPACE);
        label = gtk_label_new("IONIC");
        gtk_notebook_append_page(GTK_NOTEBOOK(notebook),page,label);
/* --- general */
        frame = gtk_frame_new("General");
        gtk_box_pack_start(GTK_BOX(page),frame,FALSE,FALSE,0);
/* create a vbox in the frame */
        vbox = gtk_vbox_new(TRUE, 0);
        gtk_container_add(GTK_CONTAINER(frame), vbox);
/* 1st line */
        VASP_NEW_LINE();
        VASP_NEW_LABEL("NSW=");
	VASP_ENTRY(vasp_gui.nsw,vasp_calc.nsw,"%i",4);
	VASP_NEW_LABEL("IBRION=");
        list = NULL;
        list = g_list_append(list,"-1:Nothing");
	list = g_list_append(list,"0:Molecular Dynamics");
	list = g_list_append(list,"1:Quasi-Newton");
	list = g_list_append(list,"2:Conjugate Gradients");
	list = g_list_append(list,"3:Damped");
	list = g_list_append(list,"5:Finite Differences");
	list = g_list_append(list,"6:Finite Diff. with SYM");
	list = g_list_append(list,"7:Perturbation Theory");
	list = g_list_append(list,"8:Pert. Theory with SYM");
	list = g_list_append(list,"44:Transition State");
	VASP_REG_COMBO(vasp_gui.ibrion,"-1:Nothing",vasp_ibrion_selected);
	VASP_NEW_SEPARATOR();
        VASP_NEW_LABEL("PSTRESS=");
        VASP_ENTRY(vasp_gui.pstress,vasp_calc.pstress,"%.4lf",6);
	VASP_NEW_LABEL("EDIFFG=");
	VASP_ENTRY(vasp_gui.ediffg,vasp_calc.ediffg,"%.2lE",8);
/* 2nd line */
	VASP_NEW_LINE();
	VASP_NEW_LABEL("ISIF=");
	list = NULL;
	list = g_list_append(list,"0:F_0_I_0_0");
	list = g_list_append(list,"1:F_P_I_0_0");
	list = g_list_append(list,"2:F_S_I_0_0");
	list = g_list_append(list,"3:F_S_I_S_V");
	list = g_list_append(list,"4:F_S_I_S_0");
	list = g_list_append(list,"5:F_S_0_S_0");
	list = g_list_append(list,"6:F_S_0_S_V");
	list = g_list_append(list,"7:F_S_0_0_V");
	VASP_REG_COMBO(vasp_gui.isif,"-1:Nothing",vasp_isif_selected);
	VASP_NEW_SEPARATOR();
	VASP_NEW_LABEL("RELAX:");
	vasp_gui.relax_ions=gui_direct_check("ions",&(vasp_gui.rions),isif_toggle,NULL,hbox);
	vasp_gui.relax_shape=gui_direct_check("cell shape",&(vasp_gui.rshape),isif_toggle,NULL,hbox);
	vasp_gui.relax_volume=gui_direct_check("cell volume",&(vasp_gui.rvolume),isif_toggle,NULL,hbox);
/* 3rd line */
	VASP_NEW_LINE();
	VASP_NEW_LABEL("POTIM=");
	VASP_ENTRY(vasp_gui.potim,vasp_calc.potim,"%.4lf",6);
	VASP_NEW_LABEL("NFREE=");
	VASP_ENTRY(vasp_gui.nfree,vasp_calc.nfree,"%i",3);
	VASP_NEW_SEPARATOR();
	VASP_NEW_LABEL("ISYM=");
	VASP_ENTRY(vasp_gui.isym,vasp_calc.isym,"%i",3);
	VASP_NEW_LABEL("SYM_PREC=");
	VASP_ENTRY(vasp_gui.sym_prec,vasp_calc.sym_prec,"%.2lE",8);
/* --- end frame */
        hbox = gtk_hbox_new(FALSE, 0);
        gtk_box_pack_end(GTK_BOX(page),hbox,FALSE,FALSE,0);
/* --- MD */
        frame = gtk_frame_new("Molecular Dynamics");
        gtk_box_pack_start(GTK_BOX(page),frame,FALSE,FALSE,0);
/* create a vbox in the frame */
        vbox = gtk_vbox_new(TRUE, 0);
        gtk_container_add(GTK_CONTAINER(frame), vbox);
/* 1st line */
        VASP_NEW_LINE();
        VASP_NEW_LABEL("TEBEG=");
	VASP_ENTRY(vasp_gui.tebeg,vasp_calc.tebeg,"%.1lf",6);
	VASP_NEW_LABEL("TEEND=");
	VASP_ENTRY(vasp_gui.teend,vasp_calc.teend,"%.1lf",6);
	VASP_NEW_LABEL("SMASS=");
	VASP_ENTRY(vasp_gui.smass,vasp_calc.smass,"%.1f",4);
	VASP_NEW_SEPARATOR();
	VASP_NEW_LABEL("NBLOCK=");
	VASP_ENTRY(vasp_gui.nblock,vasp_calc.nblock,"%i",3);
	VASP_NEW_LABEL("KBLOCK=");
	VASP_ENTRY(vasp_gui.kblock,vasp_calc.kblock,"%i",3);
/* 2nd line */
	VASP_NEW_LINE();
	VASP_NEW_LABEL("NPACO=");
	VASP_ENTRY(vasp_gui.npaco,vasp_calc.npaco,"%i",4);
	VASP_NEW_LABEL("APACO=");
	VASP_ENTRY(vasp_gui.apaco,vasp_calc.apaco,"%.4f",6);
	VASP_NEW_SEPARATOR();
	VASP_NEW_LABEL("SET IBRION=0 FOR MD");
	VASP_NEW_SEPARATOR();
	/*setup sensitive widget*/	
	if (vasp_calc.ibrion==0) {/*selecting MD enable molecular dynamics settings*/
		gtk_widget_set_sensitive(vasp_gui.tebeg,TRUE);
		gtk_widget_set_sensitive(vasp_gui.teend,TRUE);
		gtk_widget_set_sensitive(vasp_gui.smass,TRUE);
		gtk_widget_set_sensitive(vasp_gui.nblock,TRUE);
		gtk_widget_set_sensitive(vasp_gui.kblock,TRUE);
		gtk_widget_set_sensitive(vasp_gui.npaco,TRUE);
		gtk_widget_set_sensitive(vasp_gui.apaco,TRUE);
	}else{
                gtk_widget_set_sensitive(vasp_gui.tebeg,FALSE);
                gtk_widget_set_sensitive(vasp_gui.teend,FALSE);
                gtk_widget_set_sensitive(vasp_gui.smass,FALSE);
                gtk_widget_set_sensitive(vasp_gui.nblock,FALSE);
                gtk_widget_set_sensitive(vasp_gui.kblock,FALSE);
                gtk_widget_set_sensitive(vasp_gui.npaco,FALSE);
                gtk_widget_set_sensitive(vasp_gui.apaco,FALSE);
	}
/* --- end frame */
        hbox = gtk_hbox_new(FALSE, 0);
        gtk_box_pack_end(GTK_BOX(page),hbox,FALSE,FALSE,0);
/* --- POSCAR */
/*initialize some data*/
if(!vasp_gui.have_xml){
	/* FIXME: correctly recalculate lattice for non 3D-periodic!
	 * FIXME: initialize at the begining and not here.*/
	if(data->periodic<1) {
		/*not periodic: create a big cubic box*/
		data->latmat[0]=15.;data->latmat[1]=0.0;data->latmat[2]=0.0;
		data->latmat[3]=0.0;data->latmat[4]=15.;data->latmat[5]=0.0;
		data->latmat[6]=0.0;data->latmat[7]=0.0;data->latmat[8]=15.;
		vasp_calc.poscar_direct=FALSE;
	}else if(data->periodic<2){
		/* 1D-periodic? only lattice element 0 is defined: add huge y,z dimensions*/
		data->latmat[1]=0.0;data->latmat[2]=0.0;
		data->latmat[3]=0.0;data->latmat[4]=15.;data->latmat[5]=0.0;
		data->latmat[6]=0.0;data->latmat[7]=0.0;data->latmat[8]=15.;
		vasp_calc.poscar_direct=FALSE;
	}else if(data->periodic<3){
		/* 2D-periodic? only lattice element 0,3 and 1,4 are defined: add a huge z dimension*/
		/* This only make sense is surface was created by GDIS though..*/
		data->latmat[2]=0.0;
		data->latmat[5]=0.0;
		data->latmat[6]=0.0;data->latmat[7]=0.0;data->latmat[8]=(data->surface.depth)+15.;
		vasp_calc.poscar_direct=FALSE;
	}else {
		/* 3D-periodic! So far so good */
		if(data->fractional) vasp_calc.poscar_direct=TRUE;
	}
}
        frame = gtk_frame_new("POSCAR");
        gtk_box_pack_start(GTK_BOX(page),frame,FALSE,FALSE,0);
/* create a vbox in the frame */
        vbox = gtk_vbox_new(TRUE, 0);
        gtk_container_add(GTK_CONTAINER(frame), vbox);
/* 1st line */
        VASP_NEW_LINE();
	gui_direct_check("Selective Dynamics",&(vasp_calc.poscar_sd),toggle_poscar_sd,NULL,hbox);
        list = NULL;
        list = g_list_append(list,"All atom FIXED");
        list = g_list_append(list,"All atom FREE");
        list = g_list_append(list,"Manual selection");
	VASP_REG_COMBO(vasp_gui.poscar_free,"All atom FREE",vasp_poscar_free_selected);	
	VASP_NEW_SEPARATOR();
	vasp_gui.poscar_direct=gui_direct_check("Direct Coord.",&(vasp_calc.poscar_direct),toggle_poscar_direct,NULL,hbox);
/* 2nd line */
	VASP_NEW_LINE();
	VASP_NEW_LABEL("A0=");
	VASP_ENTRY(vasp_gui.poscar_a0,vasp_calc.poscar_a0,"%.4lf",6);
	VASP_NEW_SEPARATOR();
	VASP_NEW_LABEL("UX=");
	VASP_ENTRY(vasp_gui.poscar_ux,vasp_calc.poscar_ux,"%.4lf",6);
	VASP_NEW_LABEL("UY=");
	VASP_ENTRY(vasp_gui.poscar_uy,vasp_calc.poscar_uy,"%.4lf",6);
	VASP_NEW_LABEL("UZ=");
	VASP_ENTRY(vasp_gui.poscar_uz,vasp_calc.poscar_uz,"%.4lf",6);
/* 2nd line */
        VASP_NEW_LINE();
	VASP_NEW_SEPARATOR();
        VASP_NEW_LABEL("VX=");
        VASP_ENTRY(vasp_gui.poscar_vx,vasp_calc.poscar_vx,"%.4lf",6);
        VASP_NEW_LABEL("VY=");
        VASP_ENTRY(vasp_gui.poscar_vy,vasp_calc.poscar_vy,"%.4lf",6);
        VASP_NEW_LABEL("VZ=");
        VASP_ENTRY(vasp_gui.poscar_vz,vasp_calc.poscar_vz,"%.4lf",6);
/* 3rd line */
        VASP_NEW_LINE();
        VASP_NEW_SEPARATOR();
        VASP_NEW_LABEL("WX=");
        VASP_ENTRY(vasp_gui.poscar_wx,vasp_calc.poscar_wx,"%.4lf",6);
        VASP_NEW_LABEL("WY=");
        VASP_ENTRY(vasp_gui.poscar_wy,vasp_calc.poscar_wy,"%.4lf",6);
        VASP_NEW_LABEL("WZ=");
        VASP_ENTRY(vasp_gui.poscar_wz,vasp_calc.poscar_wz,"%.4lf",6);
/* 4th line */
	VASP_NEW_LINE();
	VASP_NEW_LABEL("ATOM:");
	/*this combo is special and should be explicitely defined*/
	vasp_gui.poscar_atoms = gtk_combo_box_text_new_with_entry();/*gtkcombo replacement*/
	gtk_entry_set_editable(GTK_ENTRY(gtk_bin_get_child(GTK_BIN(vasp_gui.poscar_atoms))), FALSE);/*supposed to be simple?*/
	idx=0;
	for (list2=data->cores ; list2 ; list2=g_slist_next(list2)){
		core=list2->data;
		if(vasp_calc.poscar_free==VPF_FIXED)
			gtk_combo_box_text_append_text (GTK_COMBO_BOX_TEXT (vasp_gui.poscar_atoms),g_strdup_printf("%.8lf %.8lf %.8lf  F   F   F ! atom: %i (%c%c)",
                                        core->x[0],core->x[1],core->x[2],idx,core->atom_label[0],core->atom_label[1]));
		else
			gtk_combo_box_text_append_text (GTK_COMBO_BOX_TEXT (vasp_gui.poscar_atoms),g_strdup_printf("%.8lf %.8lf %.8lf  T   T   T ! atom: %i (%c%c)",
					core->x[0],core->x[1],core->x[2],idx,core->atom_label[0],core->atom_label[1]));
		idx++;
	}
vasp_calc.atoms_total=idx;
	gtk_combo_box_text_append_text (GTK_COMBO_BOX_TEXT (vasp_gui.poscar_atoms),"ADD atom");
        gtk_box_pack_start(GTK_BOX(hbox), vasp_gui.poscar_atoms, TRUE, TRUE, 0);
	g_signal_connect(GTK_COMBO_BOX_TEXT(vasp_gui.poscar_atoms),"changed",GTK_SIGNAL_FUNC(vasp_poscar_atoms_selected),data);
/* 5th line */
	VASP_NEW_LINE();
	VASP_NEW_LABEL("@");
	VASP_ENTRY(vasp_gui.poscar_index,vasp_calc.poscar_index,"%i",4);
	gtk_widget_set_sensitive(vasp_gui.poscar_index,FALSE);/*<- is changed only by selecting poscar_atoms*/
	VASP_NEW_LABEL("SYMBOL:");
	VASP_ENTRY(vasp_gui.poscar_symbol,vasp_calc.poscar_symbol,"%s",2);
	VASP_NEW_LABEL("X=");
	VASP_ENTRY(vasp_gui.poscar_x,vasp_calc.poscar_x,"%.8lf",10);
	VASP_NEW_LABEL("Y=");
	VASP_ENTRY(vasp_gui.poscar_y,vasp_calc.poscar_y,"%.8lf",10);
	VASP_NEW_LABEL("Z=");
	VASP_ENTRY(vasp_gui.poscar_z,vasp_calc.poscar_z,"%.8lf",10);
	vasp_gui.poscar_tx=gui_direct_check("TX",&(vasp_calc.poscar_tx),NULL,NULL,hbox);/*not calling anything*/
	vasp_gui.poscar_ty=gui_direct_check("TY",&(vasp_calc.poscar_ty),NULL,NULL,hbox);/*not calling anything*/
	vasp_gui.poscar_tz=gui_direct_check("TZ",&(vasp_calc.poscar_tz),NULL,NULL,hbox);/*not calling anything*/
/*TODO: move this if outta here*/
if((vasp_calc.poscar_free==VPF_FIXED)||(vasp_calc.poscar_free==VPF_FREE)){
	gtk_widget_set_sensitive(vasp_gui.poscar_tx,FALSE);
	gtk_widget_set_sensitive(vasp_gui.poscar_ty,FALSE);
	gtk_widget_set_sensitive(vasp_gui.poscar_tz,FALSE);
	if(vasp_calc.poscar_free==VPF_FIXED){
		vasp_calc.poscar_tx=FALSE;
		vasp_calc.poscar_ty=FALSE;
		vasp_calc.poscar_tz=FALSE;
	}else{
		vasp_calc.poscar_tx=TRUE;
		vasp_calc.poscar_tx=TRUE;
		vasp_calc.poscar_tx=TRUE;
	}
}
	gui_stock_button(GTK_STOCK_APPLY, vasp_atom_modified, NULL, hbox);
	gui_stock_button(GTK_STOCK_DELETE, vasp_atom_delete, NULL, hbox);
	/* (re)initialize all */
	toggle_poscar_sd(NULL,NULL);
	gtk_combo_box_set_active(GTK_COMBO_BOX(vasp_gui.poscar_atoms),idx);/*sets "ADD atom" as active*/
	vasp_poscar_free_selected(vasp_gui.poscar_free,NULL);
/* --- end frame */
        hbox = gtk_hbox_new(FALSE, 0);
        gtk_box_pack_end(GTK_BOX(page),hbox,FALSE,FALSE,0);

/*-------------------*/
/* page 4 -> KPOINTS */
/*-------------------*/
        page = gtk_vbox_new(FALSE, _SPACE);
        label = gtk_label_new("KPOINTS");
        gtk_notebook_append_page(GTK_NOTEBOOK(notebook),page,label);
/* --- in INCAR */
        frame = gtk_frame_new("in INCAR");
        gtk_box_pack_start(GTK_BOX(page),frame,FALSE,FALSE,0);
/* create a vbox in the frame */
        vbox = gtk_vbox_new(TRUE, 0);
        gtk_container_add(GTK_CONTAINER(frame), vbox);
/* 1st line */
        VASP_NEW_LINE();
        VASP_NEW_LABEL("ISMEAR=");
        VASP_ENTRY(vasp_gui.ismear_3,vasp_calc.ismear,"%i",3);/*<- should be a combo?*/
	g_signal_connect(GTK_OBJECT(GTK_ENTRY(vasp_gui.ismear_3)),"changed",GTK_SIGNAL_FUNC(ismear_changed),data);
	VASP_NEW_SEPARATOR();
        vasp_gui.kgamma_3=gui_direct_check("KGAMMA",&(vasp_calc.kgamma),NULL,NULL,hbox);/*not calling anything*/
        VASP_NEW_LABEL("KSPACING=");
        VASP_ENTRY(vasp_gui.kspacing_3,vasp_calc.kspacing,"%.4lf",6);
	g_signal_connect(GTK_OBJECT(GTK_ENTRY(vasp_gui.kspacing_3)),"changed",GTK_SIGNAL_FUNC(kspacing_changed),data);
/* --- end frame */
        hbox = gtk_hbox_new(FALSE, 0);
        gtk_box_pack_end(GTK_BOX(page),hbox,FALSE,FALSE,0);
/* --- General */
        frame = gtk_frame_new("General");
        gtk_box_pack_start(GTK_BOX(page),frame,FALSE,FALSE,0);
/* create a vbox in the frame */
        vbox = gtk_vbox_new(TRUE, 0);
        gtk_container_add(GTK_CONTAINER(frame), vbox);
/* 1st line */
        VASP_NEW_LINE();
	VASP_NEW_LABEL("NKPTS=");
	VASP_ENTRY(vasp_gui.kpoints_nkpts,vasp_calc.kpoints_nkpts,"%i",4);
	/*using comboboxtext*/
	vasp_gui.kpoints_mode=gtk_combo_box_text_new_with_entry();
	gtk_entry_set_editable(GTK_ENTRY(gtk_bin_get_child(GTK_BIN(vasp_gui.kpoints_mode))), FALSE);
	gtk_combo_box_text_append_text (GTK_COMBO_BOX_TEXT (vasp_gui.kpoints_mode),"Manual Entry");
	gtk_combo_box_text_append_text (GTK_COMBO_BOX_TEXT (vasp_gui.kpoints_mode),"Line Mode");
	gtk_combo_box_text_append_text (GTK_COMBO_BOX_TEXT (vasp_gui.kpoints_mode),"Automatic (M&P)");
	gtk_combo_box_text_append_text (GTK_COMBO_BOX_TEXT (vasp_gui.kpoints_mode),"Gamma (M&P)");
	gtk_combo_box_text_append_text (GTK_COMBO_BOX_TEXT (vasp_gui.kpoints_mode),"Monkhorst-Pack classic");
	gtk_combo_box_text_append_text (GTK_COMBO_BOX_TEXT (vasp_gui.kpoints_mode),"Basis set definition");
	gtk_box_pack_start(GTK_BOX(hbox), vasp_gui.kpoints_mode, FALSE, FALSE, 0);
	/*end comboboxtext*/
	vasp_gui.kpoints_cart=gui_direct_check("Cartesian",&(vasp_calc.kpoints_cart),NULL,NULL,hbox);/*not calling anything*/
	vasp_gui.have_tetra=gui_direct_check("Tetrahedron",&(vasp_calc.kpoints_tetra),toogle_tetra,NULL,hbox);
/* 2nd line */
        VASP_NEW_LINE();
	vasp_gui.kpoints_kx=gui_direct_spin("KX",&(vasp_calc.kpoints_kx),1.0,100.0,1.0,NULL,NULL,hbox);/*not calling anything*/
	vasp_gui.kpoints_ky=gui_direct_spin("KY",&(vasp_calc.kpoints_ky),1.0,100.0,1.0,NULL,NULL,hbox);/*not calling anything*/
	vasp_gui.kpoints_kz=gui_direct_spin("KZ",&(vasp_calc.kpoints_kz),1.0,100.0,1.0,NULL,NULL,hbox);/*not calling anything*/
	VASP_NEW_LABEL("SHIFT:");
	VASP_NEW_LABEL("X=");
	VASP_ENTRY(vasp_gui.kpoints_sx,vasp_calc.kpoints_sx,"%.4lf",6);
	VASP_NEW_LABEL("Y=");
	VASP_ENTRY(vasp_gui.kpoints_sy,vasp_calc.kpoints_sy,"%.4lf",6);
	VASP_NEW_LABEL("Z=");
	VASP_ENTRY(vasp_gui.kpoints_sz,vasp_calc.kpoints_sz,"%.4lf",6);
/* --- end frame */
        hbox = gtk_hbox_new(FALSE, 0);
        gtk_box_pack_end(GTK_BOX(page),hbox,FALSE,FALSE,0);
/* --- Coordinate */
        vasp_gui.kpoints_coord = gtk_frame_new("Coordinate");
        gtk_box_pack_start(GTK_BOX(page),vasp_gui.kpoints_coord,FALSE,FALSE,0);
/* create a vbox in the frame */
        vbox = gtk_vbox_new(TRUE, 0);
        gtk_container_add(GTK_CONTAINER(vasp_gui.kpoints_coord), vbox);
/* 1st line */
	VASP_NEW_LINE();
	VASP_NEW_LABEL("KPOINT:");
        /*this combo is special and should be explicitely defined*/
        vasp_gui.kpoints_kpts = gtk_combo_box_text_new_with_entry();/*gtkcombo replacement*/
        gtk_entry_set_editable(GTK_ENTRY(gtk_bin_get_child(GTK_BIN(vasp_gui.kpoints_kpts))), FALSE);/*supposed to be simple?*/
	/* start with an empty list, containing only "ADD kpoint" */
	gtk_combo_box_text_append_text (GTK_COMBO_BOX_TEXT (vasp_gui.kpoints_kpts),"ADD kpoint");
	gtk_box_pack_start(GTK_BOX(hbox), vasp_gui.kpoints_kpts, TRUE, TRUE, 0);
	gtk_combo_box_set_active(GTK_COMBO_BOX(vasp_gui.kpoints_kpts),0);
	g_signal_connect(GTK_COMBO_BOX_TEXT(vasp_gui.kpoints_kpts),"changed",GTK_SIGNAL_FUNC(vasp_kpoints_kpts_selected),data);
/* 2nd line */
        VASP_NEW_LINE();
        VASP_NEW_LABEL("KPT_i=");
	VASP_ENTRY(vasp_gui.kpoints_i,vasp_calc.kpoints_i,"%i",3);
	VASP_NEW_LABEL("x=");
	VASP_ENTRY(vasp_gui.kpoints_x,vasp_calc.kpoints_x,"%.8lf",10);
	VASP_NEW_LABEL("y=");
	VASP_ENTRY(vasp_gui.kpoints_y,vasp_calc.kpoints_y,"%.8lf",10);
	VASP_NEW_LABEL("z=");
	VASP_ENTRY(vasp_gui.kpoints_z,vasp_calc.kpoints_z,"%.8lf",10);
	VASP_NEW_LABEL("w=");
	VASP_ENTRY(vasp_gui.kpoints_w,vasp_calc.kpoints_w,"%.8lf",10);
	VASP_NEW_SEPARATOR();
        gui_stock_button(GTK_STOCK_APPLY, vasp_kpoint_modified, NULL, hbox);
        gui_stock_button(GTK_STOCK_DELETE, vasp_kpoint_delete, NULL, hbox);
/* --- end frame */
        hbox = gtk_hbox_new(FALSE, 0);
        gtk_box_pack_end(GTK_BOX(page),hbox,FALSE,FALSE,0);
/* --- Tetrahedron */
        vasp_gui.kpoints_tetra = gtk_frame_new("Tetrahedron");
        gtk_box_pack_start(GTK_BOX(page),vasp_gui.kpoints_tetra,FALSE,FALSE,0);
/* create a vbox in the frame */
        vbox = gtk_vbox_new(TRUE, 0);
        gtk_container_add(GTK_CONTAINER(vasp_gui.kpoints_tetra), vbox);
/* 1st line */
        VASP_NEW_LINE();
        VASP_NEW_LABEL("TOTAL=");
	VASP_ENTRY(vasp_gui.tetra_total,vasp_calc.tetra_total,"%i",3);
	VASP_NEW_LABEL("EXACT VOLUME=");
	VASP_ENTRY(vasp_gui.tetra_volume,vasp_calc.tetra_volume,"%.8lf",10);
	VASP_NEW_LABEL("(ABCD)");
	VASP_NEW_SEPARATOR();
/* 2nd line */
	VASP_NEW_LINE();
	VASP_NEW_LABEL("Index=");
	VASP_ENTRY(vasp_gui.tetra_i,vasp_calc.tetra_i,"%i",3);
	VASP_NEW_LABEL("Degen. Weight=");
	VASP_ENTRY(vasp_gui.tetra_w,vasp_calc.tetra_w,"%.4lf",6);
	VASP_NEW_LABEL("Pt_A=");
	VASP_ENTRY(vasp_gui.tetra_a,vasp_calc.tetra_a,"%i",3);
	VASP_NEW_LABEL("Pt_B=");
	VASP_ENTRY(vasp_gui.tetra_b,vasp_calc.tetra_b,"%i",3);
	VASP_NEW_LABEL("Pt_C=");
	VASP_ENTRY(vasp_gui.tetra_c,vasp_calc.tetra_c,"%i",3);
	VASP_NEW_LABEL("Pt_D=");
	VASP_ENTRY(vasp_gui.tetra_d,vasp_calc.tetra_d,"%i",3);
	VASP_NEW_SEPARATOR();
	gui_stock_button(GTK_STOCK_APPLY, vasp_tetra_modified, NULL, hbox);
	gui_stock_button(GTK_STOCK_DELETE, vasp_tetra_delete, NULL, hbox);
/* 3rd line */
        VASP_NEW_LINE();
        VASP_NEW_LABEL("TETRAHEDRON:");
        /*this combo is special and should be explicitely defined*/
        vasp_gui.tetra = gtk_combo_box_text_new_with_entry();/*gtkcombo replacement*/
        gtk_entry_set_editable(GTK_ENTRY(gtk_bin_get_child(GTK_BIN(vasp_gui.tetra))), FALSE);/*supposed to be simple?*/
        /* start with an empty list, containing only "ADD kpoint" */
        gtk_combo_box_text_append_text (GTK_COMBO_BOX_TEXT (vasp_gui.tetra),"ADD tetrahedron");
        gtk_box_pack_start(GTK_BOX(hbox), vasp_gui.tetra, TRUE, TRUE, 0);
        gtk_combo_box_set_active(GTK_COMBO_BOX(vasp_gui.tetra),0);
        g_signal_connect(GTK_COMBO_BOX_TEXT(vasp_gui.tetra),"changed",GTK_SIGNAL_FUNC(vasp_tetra_selected),data);
/* initialize */
	g_signal_connect(GTK_COMBO_BOX_TEXT(vasp_gui.kpoints_mode),"changed",GTK_SIGNAL_FUNC(vasp_kpoints_mode_selected),data);
	gtk_combo_box_set_active(GTK_COMBO_BOX(vasp_gui.kpoints_mode),2);
	toogle_tetra();	
/* --- end frame */
        hbox = gtk_hbox_new(FALSE, 0);
        gtk_box_pack_end(GTK_BOX(page),hbox,FALSE,FALSE,0);

/*------------------*/
/* page 5 -> POTCAR */
/*------------------*/
        page = gtk_vbox_new(FALSE, _SPACE);
        label = gtk_label_new("POTCAR");
        gtk_notebook_append_page(GTK_NOTEBOOK(notebook),page,label);
/* --- general */
        frame = gtk_frame_new("General");
        gtk_box_pack_start(GTK_BOX(page),frame,FALSE,FALSE,0);
/* create a vbox in the frame */
        vbox = gtk_vbox_new(TRUE, 0);
        gtk_container_add(GTK_CONTAINER(frame), vbox);
/* 1st line */
        VASP_NEW_LINE();
        VASP_NEW_LABEL("SPECIES:");
	VASP_TEXT_ENTRY(vasp_gui.species,vasp_calc.species_symbols,TRUE);
	gtk_widget_set_sensitive(vasp_gui.species,FALSE);
/* 2nd line */
	VASP_NEW_LINE();
	VASP_NEW_SEPARATOR();
	VASP_NEW_LABEL("(each species will require a separate POTCAR information)");
	VASP_NEW_SEPARATOR();
/* --- end frame */
        hbox = gtk_hbox_new(FALSE, 0);
        gtk_box_pack_end(GTK_BOX(page),hbox,FALSE,FALSE,0);
/* --- POTCAR folder */
        frame = gtk_frame_new("POTCAR Folder");
        gtk_box_pack_start(GTK_BOX(page),frame,FALSE,FALSE,0);
/* create a vbox in the frame */
        vbox = gtk_vbox_new(TRUE, 0);
        gtk_container_add(GTK_CONTAINER(frame), vbox);
/* 1st line */
        VASP_NEW_LINE();
        VASP_NEW_LABEL("SOURCE KPOINTS FOLDER:");
        VASP_TEXT_ENTRY(vasp_gui.potcar_folder,vasp_calc.potcar_folder,TRUE);
        gtk_entry_set_text(GTK_ENTRY(vasp_gui.potcar_folder),"");
        button=gtk_button_new_from_stock(GTK_STOCK_OPEN);
        gtk_box_pack_end(GTK_BOX(hbox), button, FALSE, FALSE, 0);
        g_signal_connect(GTK_OBJECT(button),"clicked",GTK_SIGNAL_FUNC(load_potcar_folder_dialog),NULL);
/* 2nd line */
	VASP_NEW_LINE();
	VASP_NEW_LABEL("DETECTED SPECIES:");
	VASP_NEW_SEPARATOR();
	VASP_NEW_LABEL("(This method is not supported yet)");
	/*TODO: support the method*/
	VASP_NEW_SEPARATOR();
/* --- end frame */
        hbox = gtk_hbox_new(FALSE, 0);
        gtk_box_pack_end(GTK_BOX(page),hbox,FALSE,FALSE,0);
/* --- POTCAR file */
        frame = gtk_frame_new("POTCAR File");
        gtk_box_pack_start(GTK_BOX(page),frame,FALSE,FALSE,0);
/* create a vbox in the frame */
        vbox = gtk_vbox_new(TRUE, 0);
        gtk_container_add(GTK_CONTAINER(frame), vbox);
/* 1st line */
        VASP_NEW_LINE();
        VASP_NEW_LABEL("SOURCE POTCAR File:");
        VASP_TEXT_ENTRY(vasp_gui.potcar_file,vasp_calc.potcar_file,TRUE);
        gtk_entry_set_text(GTK_ENTRY(vasp_gui.potcar_file),"");
        button=gtk_button_new_from_stock(GTK_STOCK_OPEN);
        gtk_box_pack_end(GTK_BOX(hbox), button, FALSE, FALSE, 0);
        g_signal_connect(GTK_OBJECT(button),"clicked",GTK_SIGNAL_FUNC(load_potcar_file_dialog),NULL);
/* 2nd line */
        VASP_NEW_LINE();
        VASP_NEW_LABEL("DETECTED SPECIES:");
	VASP_TEXT_ENTRY(vasp_gui.potcar_species,vasp_calc.potcar_species,TRUE);
	gtk_widget_set_sensitive(vasp_gui.potcar_species,FALSE);
/* --- end frame */
        hbox = gtk_hbox_new(FALSE, 0);
        gtk_box_pack_end(GTK_BOX(page),hbox,FALSE,FALSE,0);

/*----------------*/
/* page 6 -> EXEC */
/*----------------*/
        page = gtk_vbox_new(FALSE, _SPACE);
        label = gtk_label_new("EXEC");
        gtk_notebook_append_page(GTK_NOTEBOOK(notebook),page,label);
/* --- general */
        frame = gtk_frame_new("General");
        gtk_box_pack_start(GTK_BOX(page),frame,FALSE,FALSE,0);
/* create a vbox in the frame */
        vbox = gtk_vbox_new(TRUE, 0);
        gtk_container_add(GTK_CONTAINER(frame), vbox);
/* 1st line */
        VASP_NEW_LINE();
        VASP_NEW_LABEL("VASP_EXEC:");
	VASP_TEXT_ENTRY(vasp_gui.job_vasp_exe,vasp_calc.job_vasp_exe,TRUE);
        button=gtk_button_new_from_stock(GTK_STOCK_OPEN);
        gtk_box_pack_end(GTK_BOX(hbox), button, FALSE, FALSE, 0);
        g_signal_connect(GTK_OBJECT(button),"clicked",GTK_SIGNAL_FUNC(load_vasp_exe_dialog),NULL);
/* 2nd line */
	VASP_NEW_LINE();
	VASP_NEW_LABEL("PATH:");
	VASP_TEXT_ENTRY(vasp_gui.job_path,vasp_calc.job_path,TRUE);
	button=gtk_button_new_from_stock(GTK_STOCK_OPEN);
	gtk_box_pack_end(GTK_BOX(hbox), button, FALSE, FALSE, 0);
	g_signal_connect(GTK_OBJECT(button),"clicked",GTK_SIGNAL_FUNC(load_path_dialog),NULL);
/* 3rd line */
        VASP_NEW_LINE();
	gui_direct_check("LWAVE",&(vasp_calc.lwave),NULL,NULL,hbox);
	gui_direct_check("LCHARG",&(vasp_calc.lcharg),NULL,NULL,hbox);
	gui_direct_check("LVTOT",&(vasp_calc.lvtot),NULL,NULL,hbox);
	gui_direct_check("LVHAR",&(vasp_calc.lvhar),NULL,NULL,hbox);
	gui_direct_check("LELF",&(vasp_calc.lelf),NULL,NULL,hbox);
/* --- end frame */
        hbox = gtk_hbox_new(FALSE, 0);
        gtk_box_pack_end(GTK_BOX(page),hbox,FALSE,FALSE,0);
/* --- PARALLEL */
        frame = gtk_frame_new("Parallel Optimization");
        gtk_box_pack_start(GTK_BOX(page),frame,FALSE,FALSE,0);
/* create a vbox in the frame */
        vbox = gtk_vbox_new(TRUE, 0);
        gtk_container_add(GTK_CONTAINER(frame), vbox);
/* 1st line */
        VASP_NEW_LINE();
	VASP_NEW_LABEL("MPIRUN:");
	VASP_TEXT_ENTRY(vasp_gui.job_mpirun,vasp_calc.job_mpirun,TRUE);
	button=gtk_button_new_from_stock(GTK_STOCK_OPEN);
	gtk_box_pack_end(GTK_BOX(hbox), button, FALSE, FALSE, 0);
	g_signal_connect(GTK_OBJECT(button),"clicked",GTK_SIGNAL_FUNC(load_mpirun_dialog),NULL);
/* 2nd line */
	VASP_NEW_LINE();
	vasp_gui.job_nproc=gui_direct_spin("NP=",&(vasp_calc.job_nproc),1.0,100.0,1.0,parallel_eq,NULL,hbox);
	VASP_NEW_LABEL("NCORE=");
	vasp_gui.ncore=gui_direct_spin("NCORE=",&(vasp_calc.ncore),1.0,1.0,1.0,parallel_eq,NULL,hbox);
	VASP_NEW_LABEL("KPAR=");
	vasp_gui.kpar=gui_direct_spin("KPAR=",&(vasp_calc.kpar),1.0,1.0,1.0,parallel_eq,NULL,hbox);
	gui_direct_check("LPLANE",&(vasp_calc.lplane),NULL,NULL,hbox);
	gui_direct_check("LSCALU",&(vasp_calc.lscalu),NULL,NULL,hbox);
	gui_direct_check("LSCALAPACK",&(vasp_calc.lscalapack),NULL,NULL,hbox);
/* --- DISTANT */
        frame = gtk_frame_new("Distant calculation");
        gtk_box_pack_start(GTK_BOX(page),frame,FALSE,FALSE,0);
/* create a vbox in the frame */
        vbox = gtk_vbox_new(TRUE, 0);
        gtk_container_add(GTK_CONTAINER(frame), vbox);
/* 1st line */
        VASP_NEW_LINE();
        VASP_NEW_SEPARATOR();
        spinner=gtk_spinner_new();
        gtk_box_pack_start(GTK_BOX(hbox), spinner, FALSE, FALSE, 0);
        VASP_NEW_LABEL("UNDER CONSTRUCTION");
        VASP_NEW_SEPARATOR();
        VASP_NEW_LINE();
        VASP_NEW_SEPARATOR();


/* --- Outside of notebook */
	frame = gtk_frame_new(NULL);
	gtk_box_pack_start(GTK_BOX(GTK_DIALOG(vasp_gui.window)->vbox), frame, FALSE, FALSE, 0);
	vbox = gtk_vbox_new(FALSE, _SPACE);
	gtk_container_add(GTK_CONTAINER(frame), vbox);
/* Action buttons */
	vasp_gui.spinner=gtk_spinner_new();
	gtk_box_pack_start(GTK_BOX(GTK_DIALOG(vasp_gui.window)->action_area), vasp_gui.spinner, FALSE, FALSE, 0);
	gui_stock_button(GTK_STOCK_SAVE, save_vasp_calc, NULL, GTK_DIALOG(vasp_gui.window)->action_area);
	gui_stock_button(GTK_STOCK_EXECUTE, exec_vasp_calc, NULL, GTK_DIALOG(vasp_gui.window)->action_area);
	gui_stock_button(GTK_STOCK_CLOSE, quit_vasp_calc, dialog, GTK_DIALOG(vasp_gui.window)->action_area);
/* connect to signals */
        g_signal_connect(GTK_NOTEBOOK(notebook),"switch-page",GTK_SIGNAL_FUNC(vasp_calc_page_switch),NULL);
	g_signal_connect(GTK_WIDGET(vasp_gui.spinner),"hide",GTK_SIGNAL_FUNC(vasp_spin_hide),NULL);	

/* all done */
        gtk_widget_show_all(vasp_gui.window);
	gtk_widget_hide(vasp_gui.spinner);
        sysenv.refresh_dialog=TRUE;
}

