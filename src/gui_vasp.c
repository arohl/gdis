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
	vasp_calc.have_gui=TRUE;
	vasp_calc.have_result=TRUE;
	vasp_gui.cur_page=VASP_PAGE_SIMPLIFIED;
	if(!vasp_gui.have_xml){
		vasp_gui.rions=TRUE;
		vasp_gui.rshape=FALSE;
		vasp_gui.rvolume=FALSE;
		vasp_gui.poscar_dirty=TRUE;
		vasp_gui.have_paw=FALSE;
	}
	/*default simplified interface*/
	vasp_gui.simple_rgeom=FALSE;
	vasp_gui.dimension=3.;
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
	vasp_gui_refresh();/*refresh GUI with new vasp_calc*/
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
/*********************************************/
/* refresh the GUI with value from vasp_calc */
/*********************************************/
void vasp_gui_refresh(){
/*do not refresh the simple interface, though*/
        switch(vasp_calc.prec){
        case VP_SINGLE:
		gtk_combo_box_set_active(GTK_COMBO_BOX(vasp_gui.prec),1);break;
        case VP_ACCURATE:
		gtk_combo_box_set_active(GTK_COMBO_BOX(vasp_gui.prec),2);break;
        case VP_HIGH:
		gtk_combo_box_set_active(GTK_COMBO_BOX(vasp_gui.prec),3);break;
        case VP_MED:
		gtk_combo_box_set_active(GTK_COMBO_BOX(vasp_gui.prec),4);break;
        case VP_LOW:
		gtk_combo_box_set_active(GTK_COMBO_BOX(vasp_gui.prec),5);break;
        case VP_NORM:
        default:
		gtk_combo_box_set_active(GTK_COMBO_BOX(vasp_gui.prec),0);
        }
	//use_prec already sync
	gtk_entry_set_text(GTK_ENTRY(vasp_gui.encut),g_strdup_printf("%.2lf",vasp_calc.encut));
	gtk_entry_set_text(GTK_ENTRY(vasp_gui.enaug),g_strdup_printf("%.2lf",vasp_calc.enaug));
	gtk_entry_set_text(GTK_ENTRY(vasp_gui.ediff),g_strdup_printf("%.2lE",vasp_calc.ediff));
	switch(vasp_calc.algo){
	case VA_IALGO:
		gtk_combo_box_set_active(GTK_COMBO_BOX(vasp_gui.algo),1);
		switch(vasp_calc.ialgo){
		case VIA_OE_FIXED:
			gtk_combo_box_set_active(GTK_COMBO_BOX(vasp_gui.ialgo),0);break;
		case VIA_O_FIXED:
			gtk_combo_box_set_active(GTK_COMBO_BOX(vasp_gui.ialgo),1);break;
		case VIA_SUBROT:
			gtk_combo_box_set_active(GTK_COMBO_BOX(vasp_gui.ialgo),2);break;
		case VIA_STEEP:
			gtk_combo_box_set_active(GTK_COMBO_BOX(vasp_gui.ialgo),3);break;
		case VIA_CG:
			gtk_combo_box_set_active(GTK_COMBO_BOX(vasp_gui.ialgo),4);break;
		case VIA_PSTEEP:
			gtk_combo_box_set_active(GTK_COMBO_BOX(vasp_gui.ialgo),5);break;
		case VIA_PCG:
			gtk_combo_box_set_active(GTK_COMBO_BOX(vasp_gui.ialgo),6);break;
		case VIA_ESTEEP:
			gtk_combo_box_set_active(GTK_COMBO_BOX(vasp_gui.ialgo),8);break;
		case VIA_RMMP:
			gtk_combo_box_set_active(GTK_COMBO_BOX(vasp_gui.ialgo),9);break;
		case VIA_PRMM:
			gtk_combo_box_set_active(GTK_COMBO_BOX(vasp_gui.ialgo),10);break;
		case VIA_VAR_DAMP:
			gtk_combo_box_set_active(GTK_COMBO_BOX(vasp_gui.ialgo),11);break;
		case VIA_VAR_QUENCH:
			gtk_combo_box_set_active(GTK_COMBO_BOX(vasp_gui.ialgo),12);break;
		case VIA_VAR_PCG:
			gtk_combo_box_set_active(GTK_COMBO_BOX(vasp_gui.ialgo),13);break;
		case VIA_EXACT:
			gtk_combo_box_set_active(GTK_COMBO_BOX(vasp_gui.ialgo),14);break;
		case VIA_KOSUGI:
		default:
			gtk_combo_box_set_active(GTK_COMBO_BOX(vasp_gui.ialgo),7);
		}
		break;
	case VA_VERYFAST:
		gtk_combo_box_set_active(GTK_COMBO_BOX(vasp_gui.algo),2);break;
	case VA_FAST:
		gtk_combo_box_set_active(GTK_COMBO_BOX(vasp_gui.algo),3);break;
	case VA_CONJ:
		gtk_combo_box_set_active(GTK_COMBO_BOX(vasp_gui.algo),4);break;
	case VA_ALL:
		gtk_combo_box_set_active(GTK_COMBO_BOX(vasp_gui.algo),5);break;
	case VA_DAMPED:
		gtk_combo_box_set_active(GTK_COMBO_BOX(vasp_gui.algo),6);break;
	case VA_SUBROT:
		gtk_combo_box_set_active(GTK_COMBO_BOX(vasp_gui.algo),7);break;
	case VA_EIGEN:
		gtk_combo_box_set_active(GTK_COMBO_BOX(vasp_gui.algo),8);break;
	case VA_NONE:
		gtk_combo_box_set_active(GTK_COMBO_BOX(vasp_gui.algo),9);break;
	case VA_NOTHING:
		gtk_combo_box_set_active(GTK_COMBO_BOX(vasp_gui.algo),10);break;
	case VA_EXACT:
		gtk_combo_box_set_active(GTK_COMBO_BOX(vasp_gui.algo),11);break;
	case VA_DIAG:
		gtk_combo_box_set_active(GTK_COMBO_BOX(vasp_gui.algo),12);break;
	case VA_NORM:
	default:
		gtk_combo_box_set_active(GTK_COMBO_BOX(vasp_gui.algo),0);break;
	}
	//ldiag already sync
	gtk_entry_set_text(GTK_ENTRY(vasp_gui.nsim),g_strdup_printf("%i",vasp_calc.nsim));
	gtk_entry_set_text(GTK_ENTRY(vasp_gui.vtime),g_strdup_printf("%.4lf",vasp_calc.vtime));
	gtk_entry_set_text(GTK_ENTRY(vasp_gui.iwavpr),g_strdup_printf("%i",vasp_calc.iwavpr));
	//auto already sync
	gtk_entry_set_text(GTK_ENTRY(vasp_gui.nbands),g_strdup_printf("%i",vasp_calc.nbands));
	gtk_entry_set_text(GTK_ENTRY(vasp_gui.nelect),g_strdup_printf("%.4lf",vasp_calc.nelect));
	//iniwav already sync
	gtk_entry_set_text(GTK_ENTRY(vasp_gui.istart),g_strdup_printf("%i",vasp_calc.istart));
	gtk_entry_set_text(GTK_ENTRY(vasp_gui.icharg),g_strdup_printf("%i",vasp_calc.icharg));
	/*mixer*/
	gtk_entry_set_text(GTK_ENTRY(vasp_gui.nelm),g_strdup_printf("%i",vasp_calc.nelm));
	gtk_entry_set_text(GTK_ENTRY(vasp_gui.nelmdl),g_strdup_printf("%i",vasp_calc.nelmdl));
	gtk_entry_set_text(GTK_ENTRY(vasp_gui.nelmin),g_strdup_printf("%i",vasp_calc.nelmin));
	if(!vasp_calc.auto_mixer){
		switch(vasp_calc.imix){
		case VM_0:
			gtk_combo_box_set_active(GTK_COMBO_BOX(vasp_gui.mixer),0);break;
		case VM_1:
			gtk_combo_box_set_active(GTK_COMBO_BOX(vasp_gui.mixer),1);break;
		case VM_2:
			gtk_combo_box_set_active(GTK_COMBO_BOX(vasp_gui.mixer),2);break;
		case VM_4:
		default:
			gtk_combo_box_set_active(GTK_COMBO_BOX(vasp_gui.mixer),3);
		}
		switch(vasp_calc.inimix){
		case VIM_0:
			gtk_combo_box_set_active(GTK_COMBO_BOX(vasp_gui.inimix),0);break;
		case VIM_1:
		default:
			gtk_combo_box_set_active(GTK_COMBO_BOX(vasp_gui.inimix),1);
		}
		switch(vasp_calc.mixpre){
		case VMP_0:
			gtk_combo_box_set_active(GTK_COMBO_BOX(vasp_gui.mixpre),0);break;
		case VMP_2:
			gtk_combo_box_set_active(GTK_COMBO_BOX(vasp_gui.mixpre),2);break;
		case VMP_1:
		default:
			gtk_combo_box_set_active(GTK_COMBO_BOX(vasp_gui.mixpre),1);
		}
        }
	switch(vasp_calc.lreal){
	case VLR_AUTO:
		gtk_combo_box_set_active(GTK_COMBO_BOX(vasp_gui.lreal),0);break;
	case VLR_ON:
		gtk_combo_box_set_active(GTK_COMBO_BOX(vasp_gui.lreal),1);break;
	case VLR_TRUE:
		gtk_combo_box_set_active(GTK_COMBO_BOX(vasp_gui.lreal),2);break;
	case VLR_FALSE:
	default:
		gtk_combo_box_set_active(GTK_COMBO_BOX(vasp_gui.lreal),3);break;
	}
	if(vasp_calc.ropt!=NULL) gtk_entry_set_text(GTK_ENTRY(vasp_gui.ropt),g_strdup_printf("%s",vasp_calc.ropt));
	//addgrid already sync
	gtk_entry_set_text(GTK_ENTRY(vasp_gui.lmaxmix),g_strdup_printf("%i",vasp_calc.lmaxmix));
	gtk_entry_set_text(GTK_ENTRY(vasp_gui.lmaxpaw),g_strdup_printf("%i",vasp_calc.lmaxpaw));
	if(!vasp_calc.auto_grid){
		gtk_entry_set_text(GTK_ENTRY(vasp_gui.ngx),g_strdup_printf("%i",vasp_calc.ngx));
		gtk_entry_set_text(GTK_ENTRY(vasp_gui.ngy),g_strdup_printf("%i",vasp_calc.ngy));
		gtk_entry_set_text(GTK_ENTRY(vasp_gui.ngz),g_strdup_printf("%i",vasp_calc.ngz));
		gtk_entry_set_text(GTK_ENTRY(vasp_gui.ngxf),g_strdup_printf("%i",vasp_calc.ngxf));
		gtk_entry_set_text(GTK_ENTRY(vasp_gui.ngyf),g_strdup_printf("%i",vasp_calc.ngyf));
		gtk_entry_set_text(GTK_ENTRY(vasp_gui.ngzf),g_strdup_printf("%i",vasp_calc.ngzf));
	}
	gtk_entry_set_text(GTK_ENTRY(vasp_gui.ismear),g_strdup_printf("%i",vasp_calc.ismear));
	gtk_entry_set_text(GTK_ENTRY(vasp_gui.sigma),g_strdup_printf("%.4lf",vasp_calc.sigma));
	gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(vasp_gui.kgamma),vasp_calc.kgamma);
	gtk_entry_set_text(GTK_ENTRY(vasp_gui.kspacing),g_strdup_printf("%.4lf",vasp_calc.kspacing));
	if(vasp_calc.fermwe!=NULL) gtk_entry_set_text(GTK_ENTRY(vasp_gui.fermwe),g_strdup_printf("%s",vasp_calc.fermwe));
	if(vasp_calc.fermdo!=NULL) gtk_entry_set_text(GTK_ENTRY(vasp_gui.fermdo),g_strdup_printf("%s",vasp_calc.fermdo));
	//ispin already sync
	//non_collinear already sync
	//lsorbit already sync
	//gga_compat already sync
	gtk_entry_set_text(GTK_ENTRY(vasp_gui.nupdown),g_strdup_printf("%.1lf",vasp_calc.nupdown));
	if(vasp_calc.magmom!=NULL) gtk_entry_set_text(GTK_ENTRY(vasp_gui.magmom),g_strdup_printf("%s",vasp_calc.magmom));
	if(vasp_calc.saxis!=NULL) gtk_entry_set_text(GTK_ENTRY(vasp_gui.saxis),g_strdup_printf("%s",vasp_calc.saxis));
	switch(vasp_calc.gga){
	case VG_91:
		gtk_combo_box_set_active(GTK_COMBO_BOX(vasp_gui.gga),0);break;
	case VG_RP:
		gtk_combo_box_set_active(GTK_COMBO_BOX(vasp_gui.gga),2);break;
	case VG_AM:
		gtk_combo_box_set_active(GTK_COMBO_BOX(vasp_gui.gga),3);break;
	case VG_PS:
		gtk_combo_box_set_active(GTK_COMBO_BOX(vasp_gui.gga),4);break;
	case VG_PE:
	default:
		gtk_combo_box_set_active(GTK_COMBO_BOX(vasp_gui.gga),1);
	}
	//voskown already sync
	switch (vasp_calc.mgga){
	case VMG_TPSS:
		gtk_combo_box_set_active(GTK_COMBO_BOX(vasp_gui.metagga),0);break;
	case VMG_RTPSS:
		gtk_combo_box_set_active(GTK_COMBO_BOX(vasp_gui.metagga),1);break;
	case VMG_M06L:
		gtk_combo_box_set_active(GTK_COMBO_BOX(vasp_gui.metagga),2);break;
	case VMG_MBJ:
	default:
		gtk_combo_box_set_active(GTK_COMBO_BOX(vasp_gui.metagga),3);
	}
	//lmetagga already sync
	//lasph already sync
	//lmixtau already sync
	gtk_entry_set_text(GTK_ENTRY(vasp_gui.lmaxtau),g_strdup_printf("%i",vasp_calc.lmaxtau));
	if(vasp_calc.cmbj!=NULL) gtk_entry_set_text(GTK_ENTRY(vasp_gui.cmbj),g_strdup_printf("%s",vasp_calc.cmbj));
	gtk_entry_set_text(GTK_ENTRY(vasp_gui.cmbja),g_strdup_printf("%.4lf",vasp_calc.cmbja));
	gtk_entry_set_text(GTK_ENTRY(vasp_gui.cmbjb),g_strdup_printf("%.4lf",vasp_calc.cmbjb));
	if(vasp_calc.ldau){
		switch(vasp_calc.ldau_type){
		case VU_1:
			gtk_combo_box_set_active(GTK_COMBO_BOX(vasp_gui.ldau),0);break;
		case VU_4:
			gtk_combo_box_set_active(GTK_COMBO_BOX(vasp_gui.ldau),2);break;
		case VU_2:
		default:
			gtk_combo_box_set_active(GTK_COMBO_BOX(vasp_gui.ldau),1);
		}
		switch(vasp_calc.ldau_output){
		case VUO_1:
			gtk_combo_box_set_active(GTK_COMBO_BOX(vasp_gui.ldau_print),1);break;
		case VUO_2:
			gtk_combo_box_set_active(GTK_COMBO_BOX(vasp_gui.ldau_print),2);break;
		case VUO_0:
		default:
			gtk_combo_box_set_active(GTK_COMBO_BOX(vasp_gui.ldau_print),0);
		}
		if(vasp_calc.ldaul!=NULL) gtk_entry_set_text(GTK_ENTRY(vasp_gui.ldaul),g_strdup_printf("%s",vasp_calc.ldaul));
		if(vasp_calc.ldauu!=NULL) gtk_entry_set_text(GTK_ENTRY(vasp_gui.ldauu),g_strdup_printf("%s",vasp_calc.ldauu));
		if(vasp_calc.ldauj!=NULL) gtk_entry_set_text(GTK_ENTRY(vasp_gui.ldauj),g_strdup_printf("%s",vasp_calc.ldauj));
	}
	switch(vasp_calc.idipol){
	case VID_1:
		gtk_combo_box_set_active(GTK_COMBO_BOX(vasp_gui.idipol),1);break;
	case VID_2:
		gtk_combo_box_set_active(GTK_COMBO_BOX(vasp_gui.idipol),2);break;
	case VID_3:
		gtk_combo_box_set_active(GTK_COMBO_BOX(vasp_gui.idipol),3);break;
	case VID_4:
		gtk_combo_box_set_active(GTK_COMBO_BOX(vasp_gui.idipol),4);break;
	case VID_0:
	default:
		gtk_combo_box_set_active(GTK_COMBO_BOX(vasp_gui.idipol),0);
	}
	gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(vasp_gui.ldipol),vasp_calc.ldipol);
	gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(vasp_gui.lmono),vasp_calc.lmono);
	if(vasp_calc.dipol!=NULL) gtk_entry_set_text(GTK_ENTRY(vasp_gui.dipol),g_strdup_printf("%s",vasp_calc.dipol));
	gtk_entry_set_text(GTK_ENTRY(vasp_gui.epsilon),g_strdup_printf("%.4lf",vasp_calc.epsilon));
	gtk_entry_set_text(GTK_ENTRY(vasp_gui.efield),g_strdup_printf("%.4lf",vasp_calc.efield));
	switch(vasp_calc.lorbit){
	case 1:
		gtk_combo_box_set_active(GTK_COMBO_BOX(vasp_gui.lorbit),1);break;
	case 2:
		gtk_combo_box_set_active(GTK_COMBO_BOX(vasp_gui.lorbit),2);break;
	case 5:
		gtk_combo_box_set_active(GTK_COMBO_BOX(vasp_gui.lorbit),3);break;
	case 10:
		gtk_combo_box_set_active(GTK_COMBO_BOX(vasp_gui.lorbit),4);break;
	case 11:
		gtk_combo_box_set_active(GTK_COMBO_BOX(vasp_gui.lorbit),5);break;
	case 12:
		gtk_combo_box_set_active(GTK_COMBO_BOX(vasp_gui.lorbit),6);break;
	case 0:
	default:
		gtk_combo_box_set_active(GTK_COMBO_BOX(vasp_gui.lorbit),0);
	}
	gtk_entry_set_text(GTK_ENTRY(vasp_gui.nedos),g_strdup_printf("%i",vasp_calc.nedos));
	gtk_entry_set_text(GTK_ENTRY(vasp_gui.emin),g_strdup_printf("%.2lf",vasp_calc.emin));
	gtk_entry_set_text(GTK_ENTRY(vasp_gui.emax),g_strdup_printf("%.2lf",vasp_calc.emax));
	gtk_entry_set_text(GTK_ENTRY(vasp_gui.efermi),g_strdup_printf("%.2lf",vasp_calc.efermi));
	gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(vasp_gui.have_paw),vasp_calc.have_paw);
	if(vasp_calc.rwigs!=NULL) gtk_entry_set_text(GTK_ENTRY(vasp_gui.rwigs),g_strdup_printf("%s",vasp_calc.rwigs));
	gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(vasp_gui.loptics),vasp_calc.loptics);
	gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(vasp_gui.lepsilon),vasp_calc.lepsilon);
	gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(vasp_gui.lrpa),vasp_calc.lrpa);
	gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(vasp_gui.lnabla),vasp_calc.lnabla);
	gtk_entry_set_text(GTK_ENTRY(vasp_gui.cshift),g_strdup_printf("%.2lf",vasp_calc.cshift));
	switch(vasp_calc.ibrion){
	case 0:
		gtk_combo_box_set_active(GTK_COMBO_BOX(vasp_gui.ibrion),1);break;
	case 1:
		gtk_combo_box_set_active(GTK_COMBO_BOX(vasp_gui.ibrion),2);break;
	case 2:
		gtk_combo_box_set_active(GTK_COMBO_BOX(vasp_gui.ibrion),3);break;
	case 3:
		gtk_combo_box_set_active(GTK_COMBO_BOX(vasp_gui.ibrion),4);break;
	case 5:
		gtk_combo_box_set_active(GTK_COMBO_BOX(vasp_gui.ibrion),5);break;
	case 6:
		gtk_combo_box_set_active(GTK_COMBO_BOX(vasp_gui.ibrion),6);break;
	case 7:
		gtk_combo_box_set_active(GTK_COMBO_BOX(vasp_gui.ibrion),7);break;
	case 8:
		gtk_combo_box_set_active(GTK_COMBO_BOX(vasp_gui.ibrion),8);break;
	case 44:
		gtk_combo_box_set_active(GTK_COMBO_BOX(vasp_gui.ibrion),9);break;
	case -1:
	default:
		gtk_combo_box_set_active(GTK_COMBO_BOX(vasp_gui.ibrion),0);break;
	}
	gtk_combo_box_set_active(GTK_COMBO_BOX(vasp_gui.isif),vasp_calc.isif);
	//relax_ions already sync
	//relax_shape already sync
	//relax_volume already sync
	gtk_entry_set_text(GTK_ENTRY(vasp_gui.nsw),g_strdup_printf("%i",vasp_calc.nsw));
	gtk_entry_set_text(GTK_ENTRY(vasp_gui.ediffg),g_strdup_printf("%.2lE",vasp_calc.ediffg));
	gtk_entry_set_text(GTK_ENTRY(vasp_gui.potim),g_strdup_printf("%.4lf",vasp_calc.potim));
	gtk_entry_set_text(GTK_ENTRY(vasp_gui.pstress),g_strdup_printf("%.4lf",vasp_calc.pstress));
	gtk_entry_set_text(GTK_ENTRY(vasp_gui.nfree),g_strdup_printf("%i",vasp_calc.nfree));
	if(vasp_calc.ibrion==0){
		gtk_entry_set_text(GTK_ENTRY(vasp_gui.tebeg),g_strdup_printf("%.1lf",vasp_calc.tebeg));
		gtk_entry_set_text(GTK_ENTRY(vasp_gui.teend),g_strdup_printf("%.1lf",vasp_calc.teend));
		gtk_entry_set_text(GTK_ENTRY(vasp_gui.smass),g_strdup_printf("%.1lf",vasp_calc.smass));
		gtk_entry_set_text(GTK_ENTRY(vasp_gui.nblock),g_strdup_printf("%i",vasp_calc.nblock));
		gtk_entry_set_text(GTK_ENTRY(vasp_gui.kblock),g_strdup_printf("%i",vasp_calc.kblock));
		gtk_entry_set_text(GTK_ENTRY(vasp_gui.npaco),g_strdup_printf("%i",vasp_calc.npaco));
		gtk_entry_set_text(GTK_ENTRY(vasp_gui.apaco),g_strdup_printf("%.4lf",vasp_calc.apaco));
	}
	/*skip POSCAR*/
	gtk_entry_set_text(GTK_ENTRY(vasp_gui.ismear_3),g_strdup_printf("%i",vasp_calc.ismear));
	gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(vasp_gui.kgamma_3),vasp_calc.kgamma);
	gtk_entry_set_text(GTK_ENTRY(vasp_gui.kspacing_3),g_strdup_printf("%.4lf",vasp_calc.kspacing));
	gtk_entry_set_text(GTK_ENTRY(vasp_gui.kpoints_nkpts),g_strdup_printf("%i",vasp_calc.kpoints_nkpts));
	switch(vasp_calc.kpoints_mode){
	case VKP_MAN:
		gtk_combo_box_set_active(GTK_COMBO_BOX(vasp_gui.kpoints_mode),0);
		/*we do not update the kpoint list (here)*/
		break;
	case VKP_LINE:
		gtk_combo_box_set_active(GTK_COMBO_BOX(vasp_gui.kpoints_mode),1);
		/*we do not update the kpoint list (here)*/
		break;
	case VKP_GAMMA:
		gtk_combo_box_set_active(GTK_COMBO_BOX(vasp_gui.kpoints_mode),3);
		gtk_spin_button_set_value(GTK_SPIN_BUTTON(vasp_gui.kpoints_kx),vasp_calc.kpoints_kx);
		gtk_spin_button_set_value(GTK_SPIN_BUTTON(vasp_gui.kpoints_kx),vasp_calc.kpoints_ky);
		gtk_spin_button_set_value(GTK_SPIN_BUTTON(vasp_gui.kpoints_kx),vasp_calc.kpoints_kz);
		gtk_entry_set_text(GTK_ENTRY(vasp_gui.kpoints_sx),g_strdup_printf("%.4lf",vasp_calc.kpoints_sx));
		gtk_entry_set_text(GTK_ENTRY(vasp_gui.kpoints_sy),g_strdup_printf("%.4lf",vasp_calc.kpoints_sy));
		gtk_entry_set_text(GTK_ENTRY(vasp_gui.kpoints_sz),g_strdup_printf("%.4lf",vasp_calc.kpoints_sz));
		break;
	case VKP_MP:
		gtk_combo_box_set_active(GTK_COMBO_BOX(vasp_gui.kpoints_mode),4);
		gtk_spin_button_set_value(GTK_SPIN_BUTTON(vasp_gui.kpoints_kx),vasp_calc.kpoints_kx);
		gtk_spin_button_set_value(GTK_SPIN_BUTTON(vasp_gui.kpoints_kx),vasp_calc.kpoints_ky);
		gtk_spin_button_set_value(GTK_SPIN_BUTTON(vasp_gui.kpoints_kx),vasp_calc.kpoints_kz);
		gtk_entry_set_text(GTK_ENTRY(vasp_gui.kpoints_sx),g_strdup_printf("%.4lf",vasp_calc.kpoints_sx));
		gtk_entry_set_text(GTK_ENTRY(vasp_gui.kpoints_sy),g_strdup_printf("%.4lf",vasp_calc.kpoints_sy));
		gtk_entry_set_text(GTK_ENTRY(vasp_gui.kpoints_sz),g_strdup_printf("%.4lf",vasp_calc.kpoints_sz));
		break;
	case VKP_BASIS:
		gtk_combo_box_set_active(GTK_COMBO_BOX(vasp_gui.kpoints_mode),5);
		/*we do not update the basis definition list (here)*/
		break;
	case VKP_AUTO:
	default:
		gtk_combo_box_set_active(GTK_COMBO_BOX(vasp_gui.kpoints_mode),2);break;
		gtk_spin_button_set_value(GTK_SPIN_BUTTON(vasp_gui.kpoints_kx),vasp_calc.kpoints_kx);
	}
	/*skip POTCAR*/
	/*skip EXEC*/

/*
        gtk_entry_set_text(GTK_ENTRY(vasp_gui.ismear),g_strdup_printf("%i",vasp_calc.ismear));
        gtk_combo_box_set_active(GTK_COMBO_BOX(combobox),default_value);
        gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(vasp_gui.have_tetra),TRUE);
        gtk_spin_button_set_value(GTK_SPIN_BUTTON(vasp_gui.ncore),1.0);
*/

}
/************************************************************/
/* Preload with default parameter from the simple interface */
/************************************************************/
void vasp_apply_simple(GtkButton *button, gpointer data){
	/*simple interface*/
#define _OUT(a) gtk_text_buffer_insert_at_cursor(vasp_gui.simple_message_buff,a,-1)
	gint calcul=gtk_combo_box_get_active(GTK_COMBO_BOX(vasp_gui.simple_calcul));
	gint system=gtk_combo_box_get_active(GTK_COMBO_BOX(vasp_gui.simple_system));
	gint kgrid=gtk_combo_box_get_active(GTK_COMBO_BOX(vasp_gui.simple_kgrid));
	gint dim=(gint)vasp_gui.dimension;
	/*1st take care of the obvious*/
	if((vasp_gui.simple_rgeom)&&(calcul<3)){
		_OUT("FAIL: ROUGHT GEOMETRY: OPTIMIZE GEOMETRY FIRST!\n");
		return;
	}
	if(dim==0){
		_OUT("ATOM/MOLECULE in a box, setting only gamma point\n");
		vasp_calc.kpoints_mode=VKP_GAMMA;
		vasp_calc.kpoints_kx=1;
		vasp_calc.kpoints_ky=1;
		vasp_calc.kpoints_kz=1;
	}else{
		vasp_calc.kpoints_mode=VKP_AUTO;
		vasp_calc.kpoints_kx=dim+dim*(1+2*(2*kgrid+1));/*not a serious formula*/
	}
	if(calcul<3){
		vasp_calc.prec=VP_ACCURATE;
		vasp_calc.algo=VA_NORM;
		vasp_calc.ldiag=TRUE;
	}
	vasp_calc.use_prec=TRUE;
	/*set LREAL*/
	if(vasp_calc.atoms_total>16) vasp_calc.lreal=VLR_AUTO;
	else vasp_calc.lreal=VLR_FALSE;
	/*set smearing*/
	if(system>0) {
		vasp_calc.ismear=0;
		vasp_calc.sigma=0.05;
	} else {
		vasp_calc.ismear=1;
		vasp_calc.sigma=0.2;
	}
	switch(calcul){
	case 0://Energy (single point)
		/*use ismear=-5 if we have enough kpoints (ie more than 3)*/
		if((kgrid>0)&&(dim>0)) vasp_calc.ismear=-5;
		vasp_calc.nsw=1;/*or 0*/
		break;
	case 1://DOS/BANDS (single point)
		if(dim==0) _OUT("PB: COARSE KGRID BUT BAND/DOS REQUIRED!\n");
		if(vasp_calc.have_paw) vasp_calc.lorbit=12;
		else vasp_calc.lorbit=2;
		vasp_calc.nsw=1;/*or 0*/
		vasp_calc.nedos=2001;
		vasp_calc.emin=-10.;
		vasp_calc.emax=10.;
		break;
	case 2://Lattice Dynamics (opt.)
		vasp_calc.addgrid=TRUE;
		/*according to vasp manual, Sec. 9.9*/
		//PREC=Accurate already set
		vasp_calc.lreal=VLR_FALSE;/*force LREAL=FALSE*/
		//smearing default 
		vasp_calc.ibrion=6;
		_OUT("Default to finite Difference (IBRION=6).\nSET IBRION=8 for Linear Perturbation!\n");
		break;
	case 3://Geometry (opt.)
		if(vasp_gui.simple_rgeom){
			vasp_calc.prec=VP_NORM;
			vasp_calc.use_prec=TRUE;
			vasp_calc.algo=VA_FAST;
			/*according to vasp manual, Sec. 6.2.4*/
			vasp_calc.nelmin=5;
			vasp_calc.ediff=1E-2;
			vasp_calc.ediffg=-0.3;
			vasp_calc.nsw=10;
			vasp_calc.ibrion=2;
		}else{
			vasp_calc.prec=VP_ACCURATE;
			vasp_calc.use_prec=TRUE;
			vasp_calc.algo=VA_NORM;
			vasp_calc.ldiag=TRUE;
			vasp_calc.addgrid=TRUE;
			/*according to manual, Sec. 6.2.5*/
			vasp_calc.nelmin=8;
			vasp_calc.ediff=1E-5;
			vasp_calc.ediffg=-0.01;
			vasp_calc.nsw=20;
			vasp_calc.maxmix=80;
			vasp_calc.ibrion=1;
			vasp_calc.nfree=10;
		}
		break;
	case 4://Molecular Dynamics (opt.)
		if(vasp_gui.simple_rgeom) vasp_calc.prec=VP_LOW;
		else vasp_calc.prec=VP_NORM;
		vasp_calc.use_prec=TRUE;/*<- this is maybe too high for PREC=NORMAL*/
		/*according to vasp manual, Sec. 9.7*/
		vasp_calc.ediff=1E-5;
		_OUT("Please check smearing manually!\nDefault: sigma=0.086 Fermi smearing.\n");
		vasp_calc.ismear=-1;
		vasp_calc.sigma=0.086;
		/*choose smearing wisely*/
		vasp_calc.algo=VA_VERYFAST;
		vasp_calc.maxmix=40;
		vasp_calc.isym=0;
		vasp_calc.nelmin=4;
		vasp_calc.ibrion=0;
		vasp_calc.nsw=100;
		vasp_calc.nwrite=0;/*TODO: for (adv.)*/
		vasp_calc.lcharg=FALSE;
		vasp_calc.lwave=FALSE;
		_OUT("Please set TEBEG and TEEND manually!\nDefault: 1000 for both.\n");
		vasp_calc.tebeg=1000;
		vasp_calc.teend=1000;
		vasp_calc.smass=3;
		vasp_calc.nblock=50;
		vasp_calc.potim=1.5;
		break;
	default:
		_OUT("FAIL: UNKNOWN CALCUL SETTING!\n");
		return;
	}
	vasp_gui_refresh();
#undef _OUT
}
/******************/
/* selecting PREC */
/******************/
void vasp_prec_selected(GtkWidget *w, struct model_pak *model){
	const gint index=gtk_combo_box_get_active(GTK_COMBO_BOX(w));
	switch (index){
	case 1://Single
		vasp_calc.prec=VP_SINGLE;break;
	case 2://Accurate
		vasp_calc.prec=VP_ACCURATE;break;
	case 3://High
		vasp_calc.prec=VP_HIGH;break;
	case 4://Medium
		vasp_calc.prec=VP_MED;break;
	case 5://Low
		vasp_calc.prec=VP_LOW;break;
	case 0://Normal
	default:
		vasp_calc.prec=VP_NORM;
	}
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
	const gint index=gtk_combo_box_get_active(GTK_COMBO_BOX(w));
	gtk_widget_set_sensitive(vasp_gui.ialgo,FALSE);
	switch(index){
	case 1://USE_IALGO
		vasp_calc.algo=VA_IALGO;
		gtk_widget_set_sensitive(vasp_gui.ialgo,TRUE);
		break;
	case 2://VERYFAST
		vasp_calc.algo=VA_VERYFAST;break;
	case 3://FAST
		vasp_calc.algo=VA_FAST;break;
	case 4://CONJUGATE
		vasp_calc.algo=VA_CONJ;break;
	case 5://ALL
		vasp_calc.algo=VA_ALL;break;
	case 6://DAMPED
		vasp_calc.algo=VA_DAMPED;break;
	case 7://SUBROT
		vasp_calc.algo=VA_SUBROT;break;
	case 8://EIGENVAL
		vasp_calc.algo=VA_EIGEN;break;
	case 9://NONE
		vasp_calc.algo=VA_NONE;break;
	case 10://NOTHING
		vasp_calc.algo=VA_NOTHING;break;
	case 11://EXACT
		vasp_calc.algo=VA_EXACT;break;
	case 12://DIAG
		vasp_calc.algo=VA_DIAG;break;
	case 0://NORMAL
	default:
		vasp_calc.algo=VA_NORM;
	}
}
/*******************/
/* selecting IALGO */
/*******************/
void vasp_ialgo_selected(GtkWidget *w, struct model_pak *model){
	const gint index=gtk_combo_box_get_active(GTK_COMBO_BOX(w));
	switch(index){
	case 0://2:FIXED ORB/1E
		vasp_calc.ialgo=VIA_OE_FIXED;break;
	case 1://3:FIXED ORB
		vasp_calc.ialgo=VIA_O_FIXED;break;
	case 2://4:SUBROT ONLY
		vasp_calc.ialgo=VIA_SUBROT;break;
	case 3://5:STEEPEST
		vasp_calc.ialgo=VIA_STEEP;break;
	case 4://6:CONJUGUATE GRADIENTS
		vasp_calc.ialgo=VIA_CG;break;
	case 5://7:PRECOND. STEEPEST
		vasp_calc.ialgo=VIA_PSTEEP;break;
	case 6://8:PRECOND. CG
		vasp_calc.ialgo=VIA_PCG;break;
	case 8://44:RMM STEEPEST
		vasp_calc.ialgo=VIA_ESTEEP;break;
	case 9://46:RMM PRECOND.
		vasp_calc.ialgo=VIA_RMMP;break;
	case 10://48:PRECOND. RMM
		vasp_calc.ialgo=VIA_PRMM;break;
	case 11://53:VAR DAMPED MD
		vasp_calc.ialgo=VIA_VAR_DAMP;break;
	case 12://54:VAR DAMPED QUENCH MD
		vasp_calc.ialgo=VIA_VAR_QUENCH;break;
	case 13://58:VAR PRECOND. CG
		vasp_calc.ialgo=VIA_VAR_PCG;break;
	case 14://90:EXACT
		vasp_calc.ialgo=VIA_EXACT;break;
	case 7://38:KOSUGI
	default:
		vasp_calc.ialgo=VIA_KOSUGI;
}
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
	const gint index=gtk_combo_box_get_active(GTK_COMBO_BOX(w));
	switch(index){
	case 0://Auto
		vasp_calc.lreal=VLR_AUTO;break;
	case 1://On
		vasp_calc.lreal=VLR_ON;break;
	case 2://TRUE
		vasp_calc.lreal=VLR_TRUE;break;
	case 3://FALSE
	default:
		vasp_calc.lreal=VLR_FALSE;
	}
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
	if(vasp_calc.ismear==-5) {
		gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(vasp_gui.have_tetra),TRUE);
		gtk_widget_set_sensitive(vasp_gui.have_tetra,TRUE);
	}
	else {
		gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(vasp_gui.have_tetra),FALSE);
		gtk_widget_set_sensitive(vasp_gui.have_tetra,FALSE);
	}
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
	const gint index=gtk_combo_box_get_active(GTK_COMBO_BOX(w));
	switch(index){
	case 0://91:PW91
		vasp_calc.gga=VG_91;break;
	case 2://RP:rPBE
		vasp_calc.gga=VG_RP;break;
	case 3://AM:AM05
		vasp_calc.gga=VG_AM;break;
	case 4://PS:PBEsol
		vasp_calc.gga=VG_PS;break;
	case 1://PE:PBE
	default:
		vasp_calc.gga=VG_PE;
	}
}
/*********************/
/* selecting METAGGA */
/*********************/
void vasp_metagga_selected(GtkWidget *w, struct model_pak *model){
	const gint index=gtk_combo_box_get_active(GTK_COMBO_BOX(w));
	gtk_widget_set_sensitive(vasp_gui.cmbj,FALSE);
	gtk_widget_set_sensitive(vasp_gui.cmbja,FALSE);
	gtk_widget_set_sensitive(vasp_gui.cmbjb,FALSE);
	switch(index){
	case 0://TPSS
		vasp_calc.mgga=VMG_TPSS;break;
	case 1://rTPSS
		vasp_calc.mgga=VMG_RTPSS;break;
	case 2://M06L
		vasp_calc.mgga=VMG_M06L;break;
	case 3://MBJ
	default:
		vasp_calc.mgga=VMG_MBJ;
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
		gtk_combo_box_set_active(GTK_COMBO_BOX(vasp_gui.metagga),3);//MBJ (actually there is no default)
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
		gtk_widget_set_sensitive(vasp_gui.ldaul,FALSE);
		gtk_widget_set_sensitive(vasp_gui.ldauu,FALSE);
		gtk_widget_set_sensitive(vasp_gui.ldauj,FALSE);
        }else{
                /*switch to automatic values*/
                gtk_widget_set_sensitive(vasp_gui.ldau,TRUE);
		gtk_widget_set_sensitive(vasp_gui.ldau_print,TRUE);
		gtk_widget_set_sensitive(vasp_gui.ldaul,TRUE);
		gtk_widget_set_sensitive(vasp_gui.ldauu,TRUE);
		gtk_widget_set_sensitive(vasp_gui.ldauj,TRUE);
		gtk_combo_box_set_active(GTK_COMBO_BOX(vasp_gui.ldau),1);//2:LSDA+U Dudarev
		gtk_combo_box_set_active(GTK_COMBO_BOX(vasp_gui.ldau_print),0);//0:Silent
		/*TODO: update ldaul, ldauu, ldauj*/
        }
}
/***************************/
/* selecting L(S)DA+U type */
/***************************/
void vasp_ldau_selected(GtkWidget *w, struct model_pak *model){
	const gint index=gtk_combo_box_get_active(GTK_COMBO_BOX(w));
	switch(index){
	case 0://1:LSDA+U Liechtenstein
		vasp_calc.ldau_type=VU_1;break;
	case 2://2:LSDA+U Dudarev
		vasp_calc.ldau_type=VU_4;break;
	case 1://4:LDA+U Liechtenstein
	default:
		vasp_calc.ldau_type=VU_2;
	}
}
/***************************************/
/* selecting L(S)DA+U output verbosity */
/***************************************/
void vasp_ldau_print_selected(GtkWidget *w, struct model_pak *model){
	const gint index=gtk_combo_box_get_active(GTK_COMBO_BOX(w));
	switch(index){
	case 1://1:Occupancy matrix
		vasp_calc.ldau_output=VUO_1;break;
	case 2://2:Full output
		vasp_calc.ldau_output=VUO_2;break;
	case 0://0:Silent
	default:
		vasp_calc.ldau_output=VUO_0;
	}
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
		gtk_combo_box_set_active(GTK_COMBO_BOX(vasp_gui.mixer),3);//4:BROYDEN
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
	const gint index=gtk_combo_box_get_active(GTK_COMBO_BOX(w));
	gtk_widget_set_sensitive(vasp_gui.wc,FALSE);
	gtk_widget_set_sensitive(vasp_gui.inimix,FALSE);
	gtk_widget_set_sensitive(vasp_gui.mixpre,FALSE);
	switch(index){
	case 0://0:NO MIXING
		vasp_calc.imix=VM_0;break;
	case 1://1:KERKER
		vasp_calc.imix=VM_1;break;
	case 2://2:TCHEBYCHEV
		vasp_calc.imix=VM_2;break;
	case 3://4:BROYDEN
	default:
		vasp_calc.imix=VM_4;
		gtk_widget_set_sensitive(vasp_gui.wc,TRUE);
		gtk_widget_set_sensitive(vasp_gui.inimix,TRUE);
		gtk_widget_set_sensitive(vasp_gui.mixpre,TRUE);
	}
}
/***************************/
/* selecting initial mixer */
/***************************/
void vasp_inimix_selected(GtkWidget *w, struct model_pak *model){
	const gint index=gtk_combo_box_get_active(GTK_COMBO_BOX(w));
	switch(index){
	case 0://0:LINEAR
		vasp_calc.inimix=VIM_0;break;
	case 1://1:KERKER
	default:
		vasp_calc.inimix=VIM_1;
	}
}
/************************/
/* selecting pre-mixing */
/************************/
void vasp_mixpre_selected(GtkWidget *w, struct model_pak *model){
	const gint index=gtk_combo_box_get_active(GTK_COMBO_BOX(w));
	switch(index){
	case 0://0:NONE
		vasp_calc.mixpre=VMP_0;break;
	case 2:
		vasp_calc.mixpre=VMP_2;break;
	case 1:
	default:
		vasp_calc.mixpre=VMP_1;
	}
}
/********************/
/* selecting idipol */
/********************/
void vasp_idipol_selected(GtkWidget *w, struct model_pak *model){
	const gint index=gtk_combo_box_get_active(GTK_COMBO_BOX(w));
	gtk_widget_set_sensitive(vasp_gui.ldipol,TRUE);
	gtk_widget_set_sensitive(vasp_gui.lmono,TRUE);
	gtk_widget_set_sensitive(vasp_gui.dipol,TRUE);
	gtk_widget_set_sensitive(vasp_gui.epsilon,TRUE);
	gtk_widget_set_sensitive(vasp_gui.efield,TRUE);
	switch(index){
	case 1://1:u axis
		vasp_calc.idipol=VID_1;break;
	case 2://2:v axis
		vasp_calc.idipol=VID_2;break;
	case 3://3:w axis
		vasp_calc.idipol=VID_3;break;
	case 4://4:all axis
		vasp_calc.idipol=VID_4;
		gtk_widget_set_sensitive(vasp_gui.efield,FALSE);/*not implemented*/
		break;
	case 0://0:no calcul
	default:
		vasp_calc.idipol=VID_0;
		gtk_widget_set_sensitive(vasp_gui.ldipol,FALSE);
		gtk_widget_set_sensitive(vasp_gui.lmono,FALSE);
		gtk_widget_set_sensitive(vasp_gui.dipol,FALSE);
		gtk_widget_set_sensitive(vasp_gui.epsilon,FALSE);
		gtk_widget_set_sensitive(vasp_gui.efield,FALSE);
	}
}
/********************/
/* selecting lorbit */
/********************/
void vasp_lorbit_selected(GtkWidget *w, struct model_pak *model){
	const gint index=gtk_combo_box_get_active(GTK_COMBO_BOX(w));
	switch(index){
	case 1://1:lm-PROCAR
		vasp_calc.lorbit=1;break;
	case 2://2:phase+lm-PROCAR
		vasp_calc.lorbit=2;break;
	case 3://5:PROOUT
		vasp_calc.lorbit=5;break;
	case 4://10:PROCAR
		vasp_calc.lorbit=10;break;
	case 5://11:lm-PROCAR
		vasp_calc.lorbit=11;break;
	case 6://12:phase+lm-PROCAR
		vasp_calc.lorbit=12;break;
	case 0://0:PROCAR
	default:
		vasp_calc.lorbit=0;
	}
}
/********************/
/* selecting ibrion */
/********************/
void vasp_ibrion_selected(GtkWidget *w, struct model_pak *model){
	const gint index=gtk_combo_box_get_active(GTK_COMBO_BOX(w));
	gtk_widget_set_sensitive(vasp_gui.tebeg,FALSE);
	gtk_widget_set_sensitive(vasp_gui.teend,FALSE);
	gtk_widget_set_sensitive(vasp_gui.smass,FALSE);
	gtk_widget_set_sensitive(vasp_gui.nblock,FALSE);
	gtk_widget_set_sensitive(vasp_gui.kblock,FALSE);
	gtk_widget_set_sensitive(vasp_gui.npaco,FALSE);
	gtk_widget_set_sensitive(vasp_gui.apaco,FALSE);
	switch(index){
	case 1://0:Molecular Dynamics
		vasp_calc.ibrion=0;
		/*molecular dynamics part*/
		gtk_widget_set_sensitive(vasp_gui.tebeg,TRUE);
		gtk_widget_set_sensitive(vasp_gui.teend,TRUE);
		gtk_widget_set_sensitive(vasp_gui.smass,TRUE);
		gtk_widget_set_sensitive(vasp_gui.nblock,TRUE);
		gtk_widget_set_sensitive(vasp_gui.kblock,TRUE);
		gtk_widget_set_sensitive(vasp_gui.npaco,TRUE);
		gtk_widget_set_sensitive(vasp_gui.apaco,TRUE);
		break;
	case 2://1:Quasi-Newton
		vasp_calc.ibrion=1;break;
	case 3://2:Conjugate Gradients
		vasp_calc.ibrion=2;break;
	case 4://3:Damped
		vasp_calc.ibrion=3;break;
	case 5://5:Finite Differences
		vasp_calc.ibrion=5;break;
	case 6://6:Finite Diff. with SYM
		vasp_calc.ibrion=6;break;
	case 7://7:Perturbation Theory
		vasp_calc.ibrion=7;break;
	case 8://8:Pert. Theory with SYM
		vasp_calc.ibrion=8;break;
	case 9://44:Transition State
		vasp_calc.ibrion=44;break;
	case 0://-1:Nothing
	default:
		vasp_calc.ibrion=-1;
	}
}
/******************/
/* selecting ISIF */
/******************/
void vasp_isif_selected(GtkWidget *w, struct model_pak *model){
	const gint index=gtk_combo_box_get_active(GTK_COMBO_BOX(w));
	vasp_calc.isif=index;
	switch(index){
	case 3://3:F_S_I_S_V
		gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(vasp_gui.relax_ions),TRUE);
		gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(vasp_gui.relax_shape),TRUE);
		gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(vasp_gui.relax_volume),TRUE);
		break;
	case 4://4:F_S_I_S_0
		gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(vasp_gui.relax_ions),TRUE);
		gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(vasp_gui.relax_shape),TRUE);
		gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(vasp_gui.relax_volume),FALSE);
		break;
	case 5://5:F_S_0_S_0
		gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(vasp_gui.relax_ions),FALSE);
		gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(vasp_gui.relax_shape),TRUE);
		gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(vasp_gui.relax_volume),FALSE);
		break;	
	case 6://6:F_S_0_S_V
		gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(vasp_gui.relax_ions),FALSE);
		gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(vasp_gui.relax_shape),TRUE);
		gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(vasp_gui.relax_volume),TRUE);
		break;
	case 7://7:F_S_0_0_V
		gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(vasp_gui.relax_ions),FALSE);
		gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(vasp_gui.relax_shape),FALSE);
		gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(vasp_gui.relax_volume),TRUE);
		break;
	case 0://0:F_0_I_0_0
	case 1://1:F_P_I_0_0
	case 2://2:F_S_I_0_0
	default:
		gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(vasp_gui.relax_ions),TRUE);
		gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(vasp_gui.relax_shape),FALSE);
		gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(vasp_gui.relax_volume),FALSE);
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
gtk_combo_box_set_active(GTK_COMBO_BOX(vasp_gui.isif),vasp_calc.isif);/*more simple than before*/
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
	const gint index=gtk_combo_box_get_active(GTK_COMBO_BOX(w));
	gint ix=gtk_combo_box_get_active(GTK_COMBO_BOX(vasp_gui.poscar_atoms));/*PUSH*/
	gchar *text;
	gdouble x,y,z;
	gchar symbol[3];
	gint idx=0;
	gchar tag;
        gtk_widget_set_sensitive(vasp_gui.poscar_tx,(index==2));
        gtk_widget_set_sensitive(vasp_gui.poscar_ty,(index==2));
        gtk_widget_set_sensitive(vasp_gui.poscar_tz,(index==2));
	switch(index){
	case 0://All atom FIXED
		vasp_calc.poscar_free=VPF_FIXED;
		tag='F';
		break;
	case 2://Manual selection
		vasp_calc.poscar_free=VPF_MAN;
		/*no need to modify atom list*/
		return;
	case 1://All atom FREE
	default:
		vasp_calc.poscar_free=VPF_FREE;
		tag='T';
	}
	gtk_combo_box_set_active(GTK_COMBO_BOX(vasp_gui.poscar_atoms),idx);
	text=gtk_combo_box_text_get_active_text(GTK_COMBO_BOX_TEXT(vasp_gui.poscar_atoms));
	while(g_ascii_strcasecmp(text,"ADD atom")!=0){
		/*modify each line tags*/
		sscanf(text,"%lf %lf %lf %*c %*c %*c ! atom: %*i (%[^)])",&x,&y,&z,&(symbol[0]));
			gtk_combo_box_insert_text(GTK_COMBO_BOX(vasp_gui.poscar_atoms),idx,
			g_strdup_printf("%.8lf %.8lf %.8lf  %c   %c   %c ! atom: %i (%s)",
			x,y,z,tag,tag,tag,idx,symbol));
		gtk_combo_box_text_remove(GTK_COMBO_BOX_TEXT(vasp_gui.poscar_atoms),idx+1);
		idx++;
		gtk_combo_box_set_active(GTK_COMBO_BOX(vasp_gui.poscar_atoms),idx);
		text=gtk_combo_box_text_get_active_text(GTK_COMBO_BOX_TEXT(vasp_gui.poscar_atoms));
	}
	/*return to selected atom*/
	gtk_combo_box_set_active(GTK_COMBO_BOX(vasp_gui.poscar_atoms),ix);/*PULL*/
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
		gtk_entry_set_text(GTK_ENTRY(vasp_gui.poscar_x),g_strdup_printf("%.6lf",x));
		gtk_entry_set_text(GTK_ENTRY(vasp_gui.poscar_y),g_strdup_printf("%.6lf",y));
		gtk_entry_set_text(GTK_ENTRY(vasp_gui.poscar_z),g_strdup_printf("%.6lf",z));
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
	/* (re)set tags value! */
	if(vasp_calc.poscar_tx) gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(vasp_gui.poscar_tx),TRUE);
	else gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(vasp_gui.poscar_tx),FALSE);
	if(vasp_calc.poscar_ty) gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(vasp_gui.poscar_ty),TRUE);
	else gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(vasp_gui.poscar_ty),FALSE);
	if(vasp_calc.poscar_tz) gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(vasp_gui.poscar_tz),TRUE);
	else gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(vasp_gui.poscar_tz),FALSE);
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
	tx='F';
	ty='F';
	tz='F';
	switch (vasp_calc.poscar_free){
	case VPF_MAN:
		if(vasp_calc.poscar_tx) tx='T';
		if(vasp_calc.poscar_ty) ty='T';
		if(vasp_calc.poscar_tz) tz='T';
		break;
	case VPF_FREE:
		tx='T';ty='T';tz='T';
		break;
	case VPF_FIXED:
	default:
		break;
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
	const gint index=gtk_combo_box_get_active(GTK_COMBO_BOX(w));
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
	switch(index){
	case 0://Manual Entry
		vasp_calc.kpoints_mode=VKP_MAN;
		gtk_widget_set_sensitive(vasp_gui.kpoints_w,TRUE);
		gtk_widget_set_sensitive(vasp_gui.kpoints_coord,TRUE);//need coord
		gtk_widget_set_sensitive(vasp_gui.kpoints_nkpts,TRUE);//need nkpts
		gtk_widget_set_sensitive(vasp_gui.kpoints_cart,TRUE);//can be cart
		if(vasp_calc.ismear==-5) gtk_widget_set_sensitive(vasp_gui.kpoints_tetra,TRUE);//only VKP_MAN require tetra
		if(vasp_calc.ismear==-5) gtk_widget_set_sensitive(vasp_gui.have_tetra,TRUE);//and only if ISMEAR=-5
		break;
	case 1://Line Mode
		vasp_calc.kpoints_mode=VKP_LINE;
		gtk_widget_set_sensitive(vasp_gui.kpoints_coord,TRUE);//need coord
		gtk_widget_set_sensitive(vasp_gui.kpoints_nkpts,TRUE);//need nkpts
		gtk_widget_set_sensitive(vasp_gui.kpoints_cart,TRUE);//can be cart
		break;
	case 3://Gamma (M&P)
		vasp_calc.kpoints_mode=VKP_GAMMA;
		gtk_widget_set_sensitive(vasp_gui.kpoints_kx,TRUE);//require kx
		gtk_widget_set_sensitive(vasp_gui.kpoints_ky,TRUE);//require ky
		gtk_widget_set_sensitive(vasp_gui.kpoints_kz,TRUE);//require kz
		gtk_widget_set_sensitive(vasp_gui.kpoints_sx,TRUE);//require sx
		gtk_widget_set_sensitive(vasp_gui.kpoints_sy,TRUE);//require sy
		gtk_widget_set_sensitive(vasp_gui.kpoints_sz,TRUE);//require sz
		vasp_calc.kpoints_nkpts=0;
		gtk_entry_set_text(GTK_ENTRY(vasp_gui.kpoints_nkpts),g_strdup_printf("%i",vasp_calc.kpoints_nkpts));
		gtk_widget_set_sensitive(vasp_gui.kpoints_nkpts,FALSE);//NKPTS is 0 and disable
		break;
	case 4://Monkhorst-Pack classic
		vasp_calc.kpoints_mode=VKP_MP;
		gtk_widget_set_sensitive(vasp_gui.kpoints_kx,TRUE);//require kx
		gtk_widget_set_sensitive(vasp_gui.kpoints_ky,TRUE);//require ky
		gtk_widget_set_sensitive(vasp_gui.kpoints_kz,TRUE);//require kz
		gtk_widget_set_sensitive(vasp_gui.kpoints_sx,TRUE);//require sx
		gtk_widget_set_sensitive(vasp_gui.kpoints_sy,TRUE);//require sy
		gtk_widget_set_sensitive(vasp_gui.kpoints_sz,TRUE);//require sz
		vasp_calc.kpoints_nkpts=0;
		gtk_entry_set_text(GTK_ENTRY(vasp_gui.kpoints_nkpts),g_strdup_printf("%i",vasp_calc.kpoints_nkpts));
		gtk_widget_set_sensitive(vasp_gui.kpoints_nkpts,FALSE);//NKPTS is 0 and disable
		break;
	case 5://Basis set definition
		vasp_calc.kpoints_mode=VKP_BASIS;
		gtk_widget_set_sensitive(vasp_gui.kpoints_coord,TRUE);//need coord
		gtk_widget_set_sensitive(vasp_gui.kpoints_nkpts,TRUE);//need nkpts
		gtk_widget_set_sensitive(vasp_gui.kpoints_cart,TRUE);//can be cart
		break;
	case 2://Automatic (M&P)
	default:
		vasp_calc.kpoints_mode=VKP_AUTO;
		gtk_widget_set_sensitive(vasp_gui.kpoints_kx,TRUE);//AUTO only require kx
		vasp_calc.kpoints_nkpts=0;
		gtk_entry_set_text(GTK_ENTRY(vasp_gui.kpoints_nkpts),g_strdup_printf("%i",vasp_calc.kpoints_nkpts));
		gtk_widget_set_sensitive(vasp_gui.kpoints_nkpts,FALSE);//NKPTS is 0 and disable
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
	sscanf(line," %c%c%c%*s %*s",&(sym[0]),&(sym[1]),&(sym[2]));
	if((sym[0]=='P')&&(sym[1]=='A')&&(sym[2]=='W')) {
		vasp_calc.have_paw=TRUE;
		gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(vasp_gui.have_paw),TRUE);
	}
	sym[2]='\0';
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
	/*update vasp_gui.simple_species and vasp_gui.potcar_species*/
	gtk_entry_set_text(GTK_ENTRY(vasp_gui.simple_species),g_strdup_printf("%s",vasp_calc.potcar_species));
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
/*************************************/
/* Select POTCAR information by file */
/*************************************/
void toggle_potcar_file(GtkWidget *w, GtkWidget *box){
	gtk_widget_set_sensitive(vasp_gui.potcar_folder,FALSE);
	gtk_widget_set_sensitive(vasp_gui.potcar_folder_button,FALSE);
	gtk_widget_set_sensitive(vasp_gui.species,FALSE);
	gtk_widget_set_sensitive(vasp_gui.species_flavor,FALSE);
	gtk_widget_set_sensitive(vasp_gui.species_button,FALSE);
	gtk_widget_set_sensitive(vasp_gui.potcar_file,TRUE);
	gtk_widget_set_sensitive(vasp_gui.potcar_file_button,TRUE);
}
/*************************************/
/* Select POTCAR information by path */
/*************************************/
void toggle_potcar_folder(GtkWidget *w, GtkWidget *box){
	gtk_widget_set_sensitive(vasp_gui.potcar_file,FALSE);
	gtk_widget_set_sensitive(vasp_gui.potcar_file_button,FALSE);
	gtk_widget_set_sensitive(vasp_gui.potcar_folder,TRUE);
	gtk_widget_set_sensitive(vasp_gui.potcar_folder_button,TRUE);
	gtk_widget_set_sensitive(vasp_gui.species,TRUE);
	gtk_widget_set_sensitive(vasp_gui.species_flavor,TRUE);
	gtk_widget_set_sensitive(vasp_gui.species_button,TRUE);
}
/************************************************************************/
/* tentative to equilibrate the NPAR,NCORE,KPAR with the number of CPUs */
/************************************************************************/
void parallel_eq(GtkWidget *w, GtkWidget *box){
	/*parallel equilibration*/
	gint np,ncore,kpar;
	/* no need to sync */
	np=(gint)vasp_calc.job_nproc;
	ncore=(gint)vasp_calc.ncore;
	kpar=(gint)vasp_calc.kpar;
	/*we have np CPU*/
	if(np==1){
		//NCORE=1 ; KPAR=1
		gtk_spin_button_set_value(GTK_SPIN_BUTTON(vasp_gui.ncore),1.0);
		gtk_spin_button_set_value(GTK_SPIN_BUTTON(vasp_gui.kpar),1.0);
		gtk_spin_button_set_range(GTK_SPIN_BUTTON(vasp_gui.ncore),1.0,1.0);
		gtk_spin_button_set_range(GTK_SPIN_BUTTON(vasp_gui.kpar),1.0,1.0);
		/*same for simplified interface*/
		gtk_spin_button_set_value(GTK_SPIN_BUTTON(vasp_gui.simple_ncore),1.0);
		gtk_spin_button_set_value(GTK_SPIN_BUTTON(vasp_gui.simple_kpar),1.0);
		gtk_spin_button_set_range(GTK_SPIN_BUTTON(vasp_gui.simple_ncore),1.0,1.0);
		gtk_spin_button_set_range(GTK_SPIN_BUTTON(vasp_gui.simple_kpar),1.0,1.0);
		return;
	}
	/* ncore should be [1,ncpu/kpar] */
	if(ncore>((gint)np/kpar)) gtk_spin_button_set_value(GTK_SPIN_BUTTON(vasp_gui.ncore),(gdouble)((gint)np/kpar));
	gtk_spin_button_set_range(GTK_SPIN_BUTTON(vasp_gui.ncore),1.0,(gdouble)((gint)np/kpar));
	/* kpar can be [1,ncpu] */
	gtk_spin_button_set_range(GTK_SPIN_BUTTON(vasp_gui.kpar),1.0,(gdouble)np);
	/*same for the simplified interface*/
	if(ncore>((gint)np/kpar)) gtk_spin_button_set_value(GTK_SPIN_BUTTON(vasp_gui.simple_ncore),(gdouble)((gint)np/kpar));
	gtk_spin_button_set_range(GTK_SPIN_BUTTON(vasp_gui.simple_ncore),1.0,(gdouble)((gint)np/kpar));
	gtk_spin_button_set_range(GTK_SPIN_BUTTON(vasp_gui.simple_kpar),1.0,(gdouble)np);
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
	vasp_gui.cur_page=page_num;/*TODO: for (adv.)*/
	if(page_num==VASP_PAGE_SIMPLIFIED){
		gtk_spin_button_set_value(GTK_SPIN_BUTTON(vasp_gui.simple_np),vasp_calc.job_nproc);
		gtk_spin_button_set_value(GTK_SPIN_BUTTON(vasp_gui.simple_ncore),vasp_calc.ncore);
		gtk_spin_button_set_value(GTK_SPIN_BUTTON(vasp_gui.simple_kpar),vasp_calc.kpar);
	} else if(page_num==VASP_PAGE_ELECTRONIC){
		gtk_entry_set_text(GTK_ENTRY(vasp_gui.ismear),g_strdup_printf("%i",vasp_calc.ismear));
		gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(vasp_gui.kgamma),vasp_calc.kgamma);
		gtk_entry_set_text(GTK_ENTRY(vasp_gui.kspacing),g_strdup_printf("%.4lf",vasp_calc.kspacing));
	} else if (page_num==VASP_PAGE_IONIC){
		vasp_poscar_sync();/*arrange poscar atom values*/
	} else if (page_num==VASP_PAGE_KPOINTS){
		gtk_entry_set_text(GTK_ENTRY(vasp_gui.ismear_3),g_strdup_printf("%i",vasp_calc.ismear));
		gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(vasp_gui.kgamma_3),vasp_calc.kgamma);
		gtk_entry_set_text(GTK_ENTRY(vasp_gui.kspacing_3),g_strdup_printf("%.4lf",vasp_calc.kspacing));
		gtk_combo_box_set_active(GTK_COMBO_BOX(vasp_gui.kpoints_mode),vasp_calc.kpoints_mode);
		toogle_tetra();/*in case we need*/
	} else if (page_num==VASP_PAGE_POTCAR){
		vasp_poscar_sync();/*we need to know species_symbols*/
		gtk_entry_set_text(GTK_ENTRY(vasp_gui.poscar_species),g_strdup_printf("%s",vasp_calc.species_symbols));
	} else if (page_num==VASP_PAGE_EXEC){
		gtk_spin_button_set_value(GTK_SPIN_BUTTON(vasp_gui.job_nproc),vasp_calc.job_nproc);
		gtk_spin_button_set_value(GTK_SPIN_BUTTON(vasp_gui.ncore),vasp_calc.ncore);
		gtk_spin_button_set_value(GTK_SPIN_BUTTON(vasp_gui.kpar),vasp_calc.kpar);
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
	//lorbit already sync
	VASP_REG_VAL(nedos,"%i");
	VASP_REG_VAL(emin,"%lf");
	VASP_REG_VAL(emax,"%lf");
	VASP_REG_VAL(efermi,"%lf");
	VASP_REG_TEXT(rwigs);
	//loptics already sync
	//lepsilon already sync
	//lrpa already sync
	//lnabla already sync
	VASP_REG_VAL(cshift,"%lf");
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
	//ncore already sync
	//kpar already sync
	//job_nproc already sync
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
		fprintf(output,"%i\n",(gint)vasp_calc.kpoints_kx);
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
/***********************************/
/* copy or concatenation of POTCAR */
/***********************************/
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
/*********************************/
/* Cleanup calculation structure */
/*********************************/
void vasp_cleanup(){
        /*nothing much to do but free all*/
	if(!vasp_calc.have_result) return;/*we have ongoing calculation*/
	vasprun_free(&vasp_calc);
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
	/*direct launch*/
	if((*calc).job_nproc<2) cmd = g_strdup_printf("%s > vasp.log",(*calc).job_vasp_exe);
	else cmd = g_strdup_printf("%s -np %i %s > vasp.log",(*calc).job_mpirun,(gint)(*calc).job_nproc,(*calc).job_vasp_exe);
#endif
	cwd=sysenv.cwd;/*push*/
	sysenv.cwd=g_strdup_printf("%s",(*calc).job_path);
	task_sync(cmd);
	g_free(sysenv.cwd);
	sysenv.cwd=cwd;/*pull*/
	g_free(cmd);
	(*calc).have_result=TRUE;
}
/*****************************/
/* cleanup task: load result */
/*****************************/
void cleanup_vasp_exec(vasp_calc_struct *calc){
	/*VASP process has exit, try to load the result?*/
	gchar *filename;
	/*sync_ wait for result?*/
	while(!(*calc).have_result) usleep(50*1000);/*sleep 50ms until job is done*/
	if(vasp_calc.have_gui){/*gui still active*/
		gtk_spinner_stop(GTK_SPINNER(vasp_gui.spinner));
		gtk_widget_hide(vasp_gui.spinner);
		gtk_widget_set_sensitive(vasp_gui.button_save,TRUE);
		gtk_widget_set_sensitive(vasp_gui.button_exec,TRUE);
	}
	/*init a model*/
        if((*calc).result_model!=NULL) g_free((*calc).result_model);
	(*calc).result_model=model_new();
        model_init((*calc).result_model);
	/*put result into it*/
	filename=g_strdup_printf("%s/vasprun.xml",(*calc).job_path);
	file_load(filename,(*calc).result_model);
	/*display model: we don't make result model active in case several calculations are running*/
	model_prep((*calc).result_model);
	tree_model_add((*calc).result_model);
	tree_model_refresh((*calc).result_model);
	canvas_shuffle();
	redraw_canvas(ALL);/*ALL: necessary?*/
	vasp_calc.have_result=TRUE;
	if(!vasp_calc.have_gui){/*gui was closed*/
		vasp_cleanup();/*remove calculation structure*/
	}
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
	vasp_calc.have_result=FALSE;
	gtk_widget_show(vasp_gui.spinner);
	gtk_widget_set_sensitive(vasp_gui.button_save,FALSE);
	gtk_widget_set_sensitive(vasp_gui.button_exec,FALSE);
	gtk_spinner_start(GTK_SPINNER(vasp_gui.spinner));
	task_new("VASP", &run_vasp_exec,&vasp_calc,&cleanup_vasp_exec,&vasp_calc,data);
}
/********************************/
/* Quit vasp calculation dialog */
/********************************/
void quit_vasp_calc(GtkWidget *w, gpointer data){
	struct model_pak *model=data;
	/* first cleanup everything */
//	vasp_cleanup();/*do not clean vasp_exec due to "smart" realloc...*/
	vasp_calc.have_gui=FALSE;
	/* Then bailout */
	dialog_destroy(w,model);
}
/********************************/
/* VASP calculation main dialog */
/********************************/
void gui_vasp_dialog(void){
/*launch the main interface*/
/*TODO: use distant connections*/
	gchar *title;
	gpointer dialog;
	GtkWidget *frame, *vbox, *hbox, *table;
	GtkWidget *notebook, *page, *button, *label;
	GtkWidget *separator;
	/* special */
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
	dialog = dialog_request(CVASP, title, NULL, vasp_cleanup,data);
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
/* --- general */
        frame = gtk_frame_new("General");
        gtk_box_pack_start(GTK_BOX(page),frame,FALSE,FALSE,0);
/* create a table in the frame*/
        table = gtk_table_new(4, 4, FALSE);
        gtk_container_add(GTK_CONTAINER(frame), table);
/* 1st line */
        VASP_COMBOBOX_TABLE(vasp_gui.simple_calcul,"CALCUL:",0,1,0,1);
	VASP_COMBOBOX_ADD(vasp_gui.simple_calcul,"Energy (single point)");
	VASP_COMBOBOX_ADD(vasp_gui.simple_calcul,"DOS/BANDS (single point)");
	VASP_COMBOBOX_ADD(vasp_gui.simple_calcul,"Lattice Dynamics (opt.)");
	VASP_COMBOBOX_ADD(vasp_gui.simple_calcul,"Geometry (opt.)");
	VASP_COMBOBOX_ADD(vasp_gui.simple_calcul,"Molecular Dynamics (opt.)");
	VASP_CHECK_TABLE(button,vasp_gui.simple_rgeom,NULL,"ROUGH GEOM.",2,3,0,1);/*not calling anything*/
	VASP_SPIN_TABLE(vasp_gui.simple_dim,vasp_gui.dimension,NULL,"DIM:",3,4,0,1);/*not calling anything*/
	gtk_spin_button_set_range(GTK_SPIN_BUTTON(vasp_gui.simple_dim),0.0,3.0);/*0D~3D*/
/* 2nd line */
        VASP_COMBOBOX_TABLE(vasp_gui.simple_system,"SYSTEM:",0,1,1,2);
        VASP_COMBOBOX_ADD(vasp_gui.simple_system,"Metal");
        VASP_COMBOBOX_ADD(vasp_gui.simple_system,"Semi-conducting");
        VASP_COMBOBOX_ADD(vasp_gui.simple_system,"Insulator");
        VASP_COMBOBOX_ADD(vasp_gui.simple_system,"Unknown");
        VASP_TEXT_TABLE(vasp_gui.simple_species,vasp_calc.potcar_species,"SPECIES:",1,4,1,2);
        gtk_widget_set_sensitive(vasp_gui.simple_species,FALSE);
/* 3rd line */
        VASP_COMBOBOX_TABLE(vasp_gui.simple_kgrid,"KPT GRID:",0,1,2,3);
        VASP_COMBOBOX_ADD(vasp_gui.simple_kgrid,"COARSE");
        VASP_COMBOBOX_ADD(vasp_gui.simple_kgrid,"MEDIUM");
        VASP_COMBOBOX_ADD(vasp_gui.simple_kgrid,"FINE");
        VASP_COMBOBOX_ADD(vasp_gui.simple_kgrid,"ULTRA-FINE");
	VASP_TEXT_TABLE(vasp_gui.simple_potcar,vasp_calc.potcar_file,"POTCAR:",1,3,2,3);
	VASP_BUTTON_TABLE(vasp_gui.simple_potcar_button,GTK_STOCK_OPEN,load_potcar_file_dialog,3,4,2,3);
/* 4th line */
	VASP_LABEL_TABLE("EXEC:",0,1,3,4);
        VASP_SPIN_TABLE(vasp_gui.simple_np,vasp_calc.job_nproc,parallel_eq,"NP=",1,2,3,4);
        VASP_SPIN_TABLE(vasp_gui.simple_ncore,vasp_calc.ncore,parallel_eq,"NCORE=",2,3,3,4);
        VASP_SPIN_TABLE(vasp_gui.simple_kpar,vasp_calc.kpar,parallel_eq,"KPAR=",3,4,3,4);
/* initialize */
	gtk_combo_box_set_active(GTK_COMBO_BOX(vasp_gui.simple_calcul),0);/*not calling anything*/
	gtk_combo_box_set_active(GTK_COMBO_BOX(vasp_gui.simple_system),0);/*not calling anything*/
	gtk_combo_box_set_active(GTK_COMBO_BOX(vasp_gui.simple_kgrid),0);/*not calling anything*/
/* --- end frame */
/* --- DISCLAMER */
        frame = gtk_frame_new("DISCLAMER");
        gtk_box_pack_start(GTK_BOX(page),frame,FALSE,FALSE,0);
/* create a table in the frame*/
        table = gtk_table_new(4, 3,FALSE);
        gtk_container_add(GTK_CONTAINER(frame), table);
//      gtk_container_set_border_width(GTK_CONTAINER(vasp_gui.kpoints_coord), PANEL_SPACING);/*useful?*/
/* 1st line */
	label = gtk_label_new("This simplified interface is aimed at pre-loading parameters with default values");
	gtk_table_attach_defaults(GTK_TABLE(table),label,1,2,0,1);
/* 2nd line */
	label = gtk_label_new("for a given calculation type. User should then check said defaults extremely carefuly...");
	gtk_table_attach_defaults(GTK_TABLE(table),label,1,2,1,2);
/* 3th line */
	label = gtk_label_new("GENERATE ==>");
	gtk_table_attach_defaults(GTK_TABLE(table),label,0,1,3,4);
	VASP_BUTTON_TABLE(vasp_gui.simple_apply,GTK_STOCK_APPLY,vasp_apply_simple,1,2,3,4);
	label = gtk_label_new("<== (and check...)");
	gtk_table_attach_defaults(GTK_TABLE(table),label,2,3,3,4);
/* --- end frame */
/* --- MESSAGE */
        frame = gtk_frame_new("MESSAGE");
        gtk_box_pack_start(GTK_BOX(page),frame,FALSE,FALSE,0);
/*create a vbox in frame*/
	vbox=gtk_vbox_new(TRUE, 0);
        gtk_container_add(GTK_CONTAINER(frame),vbox);
/* num frames */
	vasp_gui.simple_message=gtk_text_view_new();
	vasp_gui.simple_message_buff=gtk_text_view_get_buffer(GTK_TEXT_VIEW(vasp_gui.simple_message));
	gtk_text_view_set_editable(GTK_TEXT_VIEW(vasp_gui.simple_message),FALSE);
	gtk_widget_set_sensitive(vasp_gui.simple_message,FALSE);
	button=gtk_scrolled_window_new(NULL, NULL);/*GTK-ABSURD*/
	gtk_scrolled_window_set_policy(GTK_SCROLLED_WINDOW(button),GTK_POLICY_AUTOMATIC,GTK_POLICY_AUTOMATIC);
	gtk_container_add(GTK_CONTAINER(button),vasp_gui.simple_message);
	gtk_container_set_border_width(GTK_CONTAINER(button),1);
	gtk_box_pack_start(GTK_BOX(vbox),button,TRUE,TRUE,0);
/* --- end frame */



/*-----------------------*/
/* page 1 -> CONVERGENCE */
/*-----------------------*/
	page = gtk_vbox_new(FALSE, _SPACE);
	label = gtk_label_new("CONVERGENCE");
	gtk_notebook_append_page(GTK_NOTEBOOK(notebook),page,label);
/* --- general */
	frame = gtk_frame_new("General");
	gtk_box_pack_start(GTK_BOX(page),frame,FALSE,FALSE,0);
/* create a table in the frame*/
	table = gtk_table_new(5, 4, FALSE);
	gtk_container_add(GTK_CONTAINER(frame), table);
/* 1st line */
        VASP_COMBOBOX_TABLE(vasp_gui.prec,"PREC=",0,1,0,1);
        VASP_COMBOBOX_ADD(vasp_gui.prec,"Normal");
        VASP_COMBOBOX_ADD(vasp_gui.prec,"Single");
        VASP_COMBOBOX_ADD(vasp_gui.prec,"Accurate");
        VASP_COMBOBOX_ADD(vasp_gui.prec,"High");
        VASP_COMBOBOX_ADD(vasp_gui.prec,"Medium");
	VASP_COMBOBOX_ADD(vasp_gui.prec,"Low");
	VASP_CHECK_TABLE(button,vasp_calc.use_prec,prec_toggle,"USE PREC",1,2,0,1);
	VASP_ENTRY_TABLE(vasp_gui.encut,vasp_calc.encut,"%.2lf","ENCUT=",2,3,0,1);
	VASP_ENTRY_TABLE(vasp_gui.enaug,vasp_calc.enaug,"%.2lf","ENAUG=",3,4,0,1);
	VASP_ENTRY_TABLE(vasp_gui.ediff,vasp_calc.ediff,"%.2lE","EDIFF=",4,5,0,1);
/* 2nd line */
	VASP_COMBOBOX_TABLE(vasp_gui.algo,"ALGO=",0,1,1,2);
	VASP_COMBOBOX_ADD(vasp_gui.algo,"NORMAL");
	VASP_COMBOBOX_ADD(vasp_gui.algo,"USE_IALGO");
	VASP_COMBOBOX_ADD(vasp_gui.algo,"VERYFAST");
	VASP_COMBOBOX_ADD(vasp_gui.algo,"FAST");
	VASP_COMBOBOX_ADD(vasp_gui.algo,"CONJUGATE");
	VASP_COMBOBOX_ADD(vasp_gui.algo,"ALL");
	VASP_COMBOBOX_ADD(vasp_gui.algo,"DAMPED");
	VASP_COMBOBOX_ADD(vasp_gui.algo,"SUBROT");
	VASP_COMBOBOX_ADD(vasp_gui.algo,"EIGENVAL");
	VASP_COMBOBOX_ADD(vasp_gui.algo,"NONE");
	VASP_COMBOBOX_ADD(vasp_gui.algo,"NOTHING");
	VASP_COMBOBOX_ADD(vasp_gui.algo,"EXACT");
	VASP_COMBOBOX_ADD(vasp_gui.algo,"DIAG");
	VASP_CHECK_TABLE(button,vasp_calc.ldiag,NULL,"LDIAG",1,2,1,2);/*not calling anything*/
	VASP_ENTRY_TABLE(vasp_gui.nsim,vasp_calc.nsim,"%i","NSIM=",2,3,1,2);
	VASP_ENTRY_TABLE(vasp_gui.vtime,vasp_calc.vtime,"%.4lf","TIME=",3,4,1,2);
	VASP_ENTRY_TABLE(vasp_gui.iwavpr,vasp_calc.iwavpr,"%i","INIWAVPR=",4,5,1,2);
/* 3rd line */
	VASP_COMBOBOX_TABLE(vasp_gui.ialgo,"IALGO=",0,1,2,3);
	VASP_COMBOBOX_ADD(vasp_gui.ialgo,"2:FIXED ORB/1E");
	VASP_COMBOBOX_ADD(vasp_gui.ialgo,"3:FIXED ORB");
	VASP_COMBOBOX_ADD(vasp_gui.ialgo,"4:SUBROT ONLY");
	VASP_COMBOBOX_ADD(vasp_gui.ialgo,"5:STEEPEST");
	VASP_COMBOBOX_ADD(vasp_gui.ialgo,"6:CONJUGUATE GRADIENTS");
	VASP_COMBOBOX_ADD(vasp_gui.ialgo,"7:PRECOND. STEEPEST");
	VASP_COMBOBOX_ADD(vasp_gui.ialgo,"8:PRECOND. CG");
	VASP_COMBOBOX_ADD(vasp_gui.ialgo,"38:KOSUGI");
	VASP_COMBOBOX_ADD(vasp_gui.ialgo,"44:RMM STEEPEST");
	VASP_COMBOBOX_ADD(vasp_gui.ialgo,"46:RMM PRECOND.");
	VASP_COMBOBOX_ADD(vasp_gui.ialgo,"48:PRECOND. RMM");
	VASP_COMBOBOX_ADD(vasp_gui.ialgo,"53:VAR DAMPED MD");
	VASP_COMBOBOX_ADD(vasp_gui.ialgo,"54:VAR DAMPED QUENCH MD");
	VASP_COMBOBOX_ADD(vasp_gui.ialgo,"58:VAR PRECOND. CG");
	VASP_COMBOBOX_ADD(vasp_gui.ialgo,"90:EXACT");
	VASP_CHECK_TABLE(button,vasp_calc.auto_elec,elec_toggle,"AUTO",1,2,2,3);
	VASP_ENTRY_TABLE(vasp_gui.nbands,vasp_calc.nbands,"%i","NBANDS=",2,3,2,3);
	VASP_ENTRY_TABLE(vasp_gui.nelect,vasp_calc.nelect,"%.2lf","NELECT=",3,4,2,3);
	/*empty: col4*/
/* 4th line */
	/*empty: col0*/
	VASP_CHECK_TABLE(button,vasp_calc.iniwav,NULL,"INIWAV",1,2,3,4);/*not calling anything*/
	VASP_ENTRY_TABLE(vasp_gui.istart,vasp_calc.istart,"%i","ISTART=",2,3,3,4);
	VASP_ENTRY_TABLE(vasp_gui.icharg,vasp_calc.icharg,"%i","ICHARG=",3,4,3,4);
	/*empty: col4*/
/* initialize */
	VASP_COMBOBOX_SETUP(vasp_gui.prec,0,vasp_prec_selected);
	VASP_COMBOBOX_SETUP(vasp_gui.algo,0,vasp_algo_selected);
	VASP_COMBOBOX_SETUP(vasp_gui.ialgo,7,vasp_ialgo_selected);
	gtk_widget_set_sensitive(vasp_gui.ialgo,FALSE);
	prec_toggle(NULL,NULL);/*initialize TODO: move to end*/
	elec_toggle(NULL,NULL);/*initialize TODO: move to end*/
/* --- end frame */
/* --- mixing */
        frame = gtk_frame_new("Mixing");
        gtk_box_pack_start(GTK_BOX(page),frame,FALSE,FALSE,0);
/* create a table in the frame*/
        table = gtk_table_new(5, 3, FALSE);
        gtk_container_add(GTK_CONTAINER(frame), table);
/* 1st line */
	VASP_COMBOBOX_TABLE(vasp_gui.mixer,"IMIX=",0,1,0,1);
	VASP_COMBOBOX_ADD(vasp_gui.mixer,"0:NONE");
	VASP_COMBOBOX_ADD(vasp_gui.mixer,"1:KERKER");
	VASP_COMBOBOX_ADD(vasp_gui.mixer,"2:TCHEBYCHEV");
	VASP_COMBOBOX_ADD(vasp_gui.mixer,"4:BROYDEN");
	VASP_CHECK_TABLE(button,vasp_calc.auto_mixer,mixer_toggle,"AUTO MIXER",1,2,0,1);
	VASP_ENTRY_TABLE(vasp_gui.nelm,vasp_calc.nelm,"%i","NELM=",2,3,0,1);
	VASP_ENTRY_TABLE(vasp_gui.nelmdl,vasp_calc.nelmdl,"%i","NELMDL=",3,4,0,1);
	VASP_ENTRY_TABLE(vasp_gui.nelmin,vasp_calc.nelmin,"%i","NELMIN=",4,5,0,1);
/* 2nd line */
	VASP_COMBOBOX_TABLE(vasp_gui.mixpre,"MIXPRE=",0,1,1,2);
	VASP_COMBOBOX_ADD(vasp_gui.mixpre,"0:NONE");
	VASP_COMBOBOX_ADD(vasp_gui.mixpre,"1:F=20 INVERSE KERKER");
	VASP_COMBOBOX_ADD(vasp_gui.mixpre,"2:F=200 INVERSE KERKER");
	/*empty: col1*/
        VASP_ENTRY_TABLE(vasp_gui.amix,vasp_calc.amix,"%.4lf","AMIX=",2,3,1,2);
        VASP_ENTRY_TABLE(vasp_gui.bmix,vasp_calc.bmix,"%.4lf","BMIX=",3,4,1,2);
        VASP_ENTRY_TABLE(vasp_gui.amin,vasp_calc.amin,"%.4lf","AMIN=",4,5,1,2);
/* 3rd line */
	VASP_COMBOBOX_TABLE(vasp_gui.inimix,"INIMIX=",0,1,2,3);
	VASP_COMBOBOX_ADD(vasp_gui.inimix,"0:LINEAR");
	VASP_COMBOBOX_ADD(vasp_gui.inimix,"1:KERKER");
	VASP_ENTRY_TABLE(vasp_gui.maxmix,vasp_calc.maxmix,"%i","MAXMIX=",1,2,2,3);
	VASP_ENTRY_TABLE(vasp_gui.amix_mag,vasp_calc.amix_mag,"%.4lf","AMIX_MAG=",2,3,2,3);
	VASP_ENTRY_TABLE(vasp_gui.bmix_mag,vasp_calc.bmix_mag,"%.4lf","BMIX_MAG=",3,4,2,3);
	VASP_ENTRY_TABLE(vasp_gui.wc,vasp_calc.wc,"%.1lf","WC=",4,5,2,3);
/* initialize */
	VASP_COMBOBOX_SETUP(vasp_gui.mixer,3,vasp_mixer_selected);
	VASP_COMBOBOX_SETUP(vasp_gui.mixpre,1,vasp_mixpre_selected);
	VASP_COMBOBOX_SETUP(vasp_gui.inimix,1,vasp_inimix_selected);
	mixer_toggle(NULL,NULL);/*initialize TODO: move to end*/
/* --- end frame */
/* --- projector */
        frame = gtk_frame_new("Projector");
        gtk_box_pack_start(GTK_BOX(page),frame,FALSE,FALSE,0);
/* create a table in the frame*/
        table = gtk_table_new(5, 4, FALSE);
        gtk_container_add(GTK_CONTAINER(frame), table);
/* 1st line */
	VASP_COMBOBOX_TABLE(vasp_gui.lreal,"LREAL=",0,1,0,1);
	VASP_COMBOBOX_ADD(vasp_gui.lreal,"Auto");
	VASP_COMBOBOX_ADD(vasp_gui.lreal,"On");
	VASP_COMBOBOX_ADD(vasp_gui.lreal,"TRUE");
	VASP_COMBOBOX_ADD(vasp_gui.lreal,"FALSE");
	/*empty: col2*/
	VASP_TEXT_TABLE(vasp_gui.ropt,vasp_calc.ropt,"ROPT=",2,5,0,1);
/* 2nd line */
	/*empty: col1*/
	VASP_CHECK_TABLE(button,vasp_calc.addgrid,NULL,"ADDGRID",1,2,1,2);/*not calling anything*/
        VASP_ENTRY_TABLE(vasp_gui.lmaxmix,vasp_calc.lmaxmix,"%i","LMAXMIX=",2,3,1,2);
	/*empty: col4*/
        VASP_ENTRY_TABLE(vasp_gui.lmaxpaw,vasp_calc.lmaxpaw,"%i","LMAXPAW=",4,5,1,2);
/* 3rd line */
	/*empty: col1*/
	VASP_CHECK_TABLE(button,vasp_calc.auto_grid,grid_toggle,"AUTO GRID",1,2,2,3);
	VASP_ENTRY_TABLE(vasp_gui.ngx,vasp_calc.ngx,"%i","NGX=",2,3,2,3);
	VASP_ENTRY_TABLE(vasp_gui.ngy,vasp_calc.ngy,"%i","NGY=",3,4,2,3);
	VASP_ENTRY_TABLE(vasp_gui.ngz,vasp_calc.ngz,"%i","NGZ=",4,5,2,3);
/* 4th line */
	/*empty: col1*/
	/*empty: col2*/
	VASP_ENTRY_TABLE(vasp_gui.ngxf,vasp_calc.ngxf,"%i","NGXF=",2,3,3,4);
	VASP_ENTRY_TABLE(vasp_gui.ngyf,vasp_calc.ngyf,"%i","NGYF=",3,4,3,4);
	VASP_ENTRY_TABLE(vasp_gui.ngzf,vasp_calc.ngzf,"%i","NGZF=",4,5,3,4);
/* initialize */
	VASP_COMBOBOX_SETUP(vasp_gui.lreal,3,vasp_lreal_selected);
	grid_toggle(NULL,NULL);/*initialize TODO: move to end*/
/* --- end frame */

/*-------------------*/
/* page 2 -> ELECT-I */
/*-------------------*/
        page = gtk_vbox_new(FALSE, _SPACE);
        label = gtk_label_new("ELECT-I");
        gtk_notebook_append_page(GTK_NOTEBOOK(notebook),page,label);
/* --- smearing */
        frame = gtk_frame_new("Smearing");
        gtk_box_pack_start(GTK_BOX(page),frame,FALSE,FALSE,0);
/* create a table in the frame*/
        table = gtk_table_new(4, 2, FALSE);
        gtk_container_add(GTK_CONTAINER(frame), table);
//      gtk_container_set_border_width(GTK_CONTAINER(table), PANEL_SPACING);/*useful?*/
/* 1st line */
	VASP_ENTRY_TABLE(vasp_gui.ismear,vasp_calc.ismear,"%i","ISMEAR=",0,1,0,1);
	VASP_ENTRY_TABLE(vasp_gui.sigma,vasp_calc.sigma,"%.4lf","SIGMA=",1,2,0,1);
	VASP_CHECK_TABLE(vasp_gui.kgamma,vasp_calc.kgamma,NULL,"KGAMMA",2,3,0,1);/*not calling anything*/
	VASP_ENTRY_TABLE(vasp_gui.kspacing,vasp_calc.kspacing,"%.4lf","KSPACING=",3,4,0,1);
/* 2nd line */
	VASP_TEXT_TABLE(vasp_gui.fermwe,vasp_calc.fermwe,"FERMWE=",0,2,1,2);
	VASP_TEXT_TABLE(vasp_gui.fermdo,vasp_calc.fermdo,"FERMDO=",2,4,1,2);
/* --- end frame */
/* --- spin */
        frame = gtk_frame_new("Spin");
        gtk_box_pack_start(GTK_BOX(page),frame,FALSE,FALSE,0);
/* create a table in the frame*/
	table = gtk_table_new(5, 2,FALSE);
	gtk_container_add(GTK_CONTAINER(frame), table);
//      gtk_container_set_border_width(GTK_CONTAINER(table), PANEL_SPACING);/*useful?*/
/* 1st line */
	VASP_CHECK_TABLE(button,vasp_calc.ispin,NULL,"SPIN POLARIZED",0,1,0,1);/*not calling anything*/
	VASP_CHECK_TABLE(vasp_gui.lnoncoll,vasp_calc.non_collinear,lnoncoll_toggle,"NON COLLINEAR",1,2,0,1);
	VASP_CHECK_TABLE(vasp_gui.lsorbit,vasp_calc.lsorbit,lsorbit_toggle,"LSORBIT",2,3,0,1);
	VASP_CHECK_TABLE(button,vasp_calc.gga_compat,NULL,"GGA_COMPAT",3,4,0,1);/*not calling anything*/
	VASP_ENTRY_TABLE(vasp_gui.nupdown,vasp_calc.nupdown,"%.1f","NUPDOWN=",4,5,0,1);
/* 2nd line */
	VASP_TEXT_TABLE(vasp_gui.magmom,vasp_calc.magmom,"MAGMOM=",0,3,1,2);
	VASP_TEXT_TABLE(vasp_gui.saxis,vasp_calc.saxis,"SAXIS=",3,5,1,2);
/* --- end frame */
/* --- exchange correlation */
        frame = gtk_frame_new("Exchange Correlation");
        gtk_box_pack_start(GTK_BOX(page),frame,FALSE,FALSE,0);
/* create a table in the frame*/
	table = gtk_table_new(5, 5,FALSE);
	gtk_container_add(GTK_CONTAINER(frame), table);
//      gtk_container_set_border_width(GTK_CONTAINER(table), PANEL_SPACING);/*useful?*/
/* 1st line */
        VASP_COMBOBOX_TABLE(vasp_gui.gga,"GGA=",0,1,0,1);
        VASP_COMBOBOX_ADD(vasp_gui.gga,"91:PW91");
	VASP_COMBOBOX_ADD(vasp_gui.gga,"PE:PBE");
	VASP_COMBOBOX_ADD(vasp_gui.gga,"RP:rPBE");
	VASP_COMBOBOX_ADD(vasp_gui.gga,"AM:AM05");
	VASP_COMBOBOX_ADD(vasp_gui.gga,"PS:PBEsol");
	VASP_CHECK_TABLE(button,vasp_calc.voskown,NULL,"VOSKOWN",1,2,0,1);/*not calling anything*/
	/*empty: col2*/
	/*empty: col3*/
	/*empty: col4*/
/* 2nd line */
	VASP_COMBOBOX_TABLE(vasp_gui.metagga,"METAGGA=",0,1,1,2);
	VASP_COMBOBOX_ADD(vasp_gui.metagga,"TPSS");
	VASP_COMBOBOX_ADD(vasp_gui.metagga,"rTPSS");
	VASP_COMBOBOX_ADD(vasp_gui.metagga,"M06L");
	VASP_COMBOBOX_ADD(vasp_gui.metagga,"MBJ");
	VASP_CHECK_TABLE(button,vasp_calc.lmetagga,metagga_toggle,"LMETAGGA",1,2,1,2);
	VASP_CHECK_TABLE(button,vasp_calc.lasph,NULL,"LASPH",2,3,1,2);/*not calling anything*/
	VASP_CHECK_TABLE(button,vasp_calc.lmixtau,NULL,"LMIXTAU",3,4,1,2);/*not calling anything*/
	VASP_ENTRY_TABLE(vasp_gui.lmaxtau,vasp_calc.lmaxtau,"%i","LMAXTAU=",4,5,1,2);
/* 3rd line */
	VASP_TEXT_TABLE(vasp_gui.cmbj,vasp_calc.cmbj,"CMBJ=",0,3,2,3);
	VASP_ENTRY_TABLE(vasp_gui.cmbja,vasp_calc.cmbja,"%.4lf","CMBJA=",3,4,2,3);
	VASP_ENTRY_TABLE(vasp_gui.cmbjb,vasp_calc.cmbjb,"%.4lf","CMBJB=",4,5,2,3);
/* 4th line */
	VASP_COMBOBOX_TABLE(vasp_gui.ldau,"LDAUTYPE=",0,1,3,4);
	VASP_COMBOBOX_ADD(vasp_gui.ldau,"1:LSDA+U Liechtenstein");
	VASP_COMBOBOX_ADD(vasp_gui.ldau,"2:LSDA+U Dudarev");
	VASP_COMBOBOX_ADD(vasp_gui.ldau,"4:LDA+U Liechtenstein");
	VASP_CHECK_TABLE(button,vasp_calc.ldau,ldau_toggle,"LDAU",1,2,3,4);
	VASP_TEXT_TABLE(vasp_gui.ldaul,vasp_calc.ldaul,"LDAUL=",2,5,3,4);
/* 5th line */
	VASP_COMBOBOX_TABLE(vasp_gui.ldau_print,"LDAUPRINT=",0,1,4,5);
	VASP_COMBOBOX_ADD(vasp_gui.ldau_print,"0:Silent");
	VASP_COMBOBOX_ADD(vasp_gui.ldau_print,"1:Occupancy matrix");
	VASP_COMBOBOX_ADD(vasp_gui.ldau_print,"2:Full output");
	VASP_TEXT_TABLE(vasp_gui.ldauu,vasp_calc.ldauu,"LDAUU=",1,3,4,5);
	VASP_TEXT_TABLE(vasp_gui.ldauj,vasp_calc.ldauj,"LDAUJ=",3,5,4,5);
/* initialize */
	VASP_COMBOBOX_SETUP(vasp_gui.gga,1,vasp_gga_selected);
	VASP_COMBOBOX_SETUP(vasp_gui.metagga,3,vasp_metagga_selected);
	VASP_COMBOBOX_SETUP(vasp_gui.ldau,1,vasp_ldau_selected);
	VASP_COMBOBOX_SETUP(vasp_gui.ldau_print,0,vasp_ldau_print_selected);
	metagga_toggle(NULL,NULL);/*initialize TODO: move to end*/
	ldau_toggle(NULL,NULL);/*initialize TODO: move to end*/
/* --- end frame */
/* --- dipole */
        frame = gtk_frame_new("Dipole");
        gtk_box_pack_start(GTK_BOX(page),frame,FALSE,FALSE,0);
/* create a table in the frame*/
        table = gtk_table_new(2, 5,FALSE);
        gtk_container_add(GTK_CONTAINER(frame), table);
//      gtk_container_set_border_width(GTK_CONTAINER(table), PANEL_SPACING);/*useful?*/
/* 1st line */
	VASP_COMBOBOX_TABLE(vasp_gui.idipol,"IDIPOL=",0,1,0,1);
	VASP_COMBOBOX_ADD(vasp_gui.idipol,"0:no calcul");
	VASP_COMBOBOX_ADD(vasp_gui.idipol,"1:u axis");
	VASP_COMBOBOX_ADD(vasp_gui.idipol,"2:v axis");
	VASP_COMBOBOX_ADD(vasp_gui.idipol,"3:w axis");
	VASP_COMBOBOX_ADD(vasp_gui.idipol,"4:all axis");
	VASP_CHECK_TABLE(vasp_gui.ldipol,vasp_calc.ldipol,NULL,"LDIPOL",1,2,0,1);/*not calling anything*/
	VASP_CHECK_TABLE(vasp_gui.lmono,vasp_calc.lmono,NULL,"LMONO",2,3,0,1);/*not calling anything*/
	VASP_TEXT_TABLE(vasp_gui.dipol,vasp_calc.dipol,"DIPOL=",3,5,0,1);
/* 2nd line */
	/*empty: col0*/
	/*empty: col1*/
	/*empty: col2*/
	VASP_ENTRY_TABLE(vasp_gui.epsilon,vasp_calc.epsilon,"%.4lf","EPSILON=",3,4,2,3);
	VASP_ENTRY_TABLE(vasp_gui.efield,vasp_calc.efield,"%.4lf","EFIELD=",4,5,2,3);
/*initialize sensitive widget*/
	VASP_COMBOBOX_SETUP(vasp_gui.idipol,0,vasp_idipol_selected);
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

/*---------------------*/
/* page 2b -> ELECT-II */
/*---------------------*/
        page = gtk_vbox_new(FALSE, _SPACE);
        label = gtk_label_new("ELECT-II");
        gtk_notebook_append_page(GTK_NOTEBOOK(notebook),page,label);
/* --- DOS */
        frame = gtk_frame_new("Density of States");
        gtk_box_pack_start(GTK_BOX(page),frame,FALSE,FALSE,0);
/* create a table in the frame*/
        table = gtk_table_new(2, 5, FALSE);
        gtk_container_add(GTK_CONTAINER(frame), table);
//      gtk_container_set_border_width(GTK_CONTAINER(table), PANEL_SPACING);/*useful?*/
/* 1st line */
	VASP_COMBOBOX_TABLE(vasp_gui.lorbit,"LORBIT=",0,1,0,1);
	VASP_COMBOBOX_ADD(vasp_gui.lorbit,"0:PROCAR");
	VASP_COMBOBOX_ADD(vasp_gui.lorbit,"1:lm-PROCAR");
	VASP_COMBOBOX_ADD(vasp_gui.lorbit,"2:phase+lm-PROCAR");
	VASP_COMBOBOX_ADD(vasp_gui.lorbit,"5:PROOUT");
	VASP_COMBOBOX_ADD(vasp_gui.lorbit,"10:PROCAR");
	VASP_COMBOBOX_ADD(vasp_gui.lorbit,"11:lm-PROCAR");
	VASP_COMBOBOX_ADD(vasp_gui.lorbit,"12:phase+lm-PROCAR");
	VASP_ENTRY_TABLE(vasp_gui.nedos,vasp_calc.nedos,"%i","NEDOS=",1,2,0,1);
	VASP_ENTRY_TABLE(vasp_gui.emin,vasp_calc.emin,"%.2lf","EMIN=",2,3,0,1);
	VASP_ENTRY_TABLE(vasp_gui.emax,vasp_calc.emax,"%.2lf","EMAX=",3,4,0,1);
	VASP_ENTRY_TABLE(vasp_gui.efermi,vasp_calc.efermi,"%.2lf","EFERMI=",4,5,0,1);
/* 2nd line */
	VASP_LABEL_TABLE("(LORBIT>5) => Req. PAW",0,1,1,2);
	VASP_CHECK_TABLE(vasp_gui.have_paw,vasp_calc.have_paw,NULL,"HAVE PAW",1,2,1,2);/*not calling anything*/
	gtk_widget_set_sensitive(vasp_gui.have_paw,FALSE);/*just an indication*/
	VASP_TEXT_TABLE(vasp_gui.rwigs,vasp_calc.rwigs,"RWIGS=",2,5,1,2);
/* initialize */
if(vasp_calc.have_paw) VASP_COMBOBOX_SETUP(vasp_gui.lorbit,4,vasp_lorbit_selected);
else VASP_COMBOBOX_SETUP(vasp_gui.lorbit,0,vasp_lorbit_selected);
/* --- end frame */
/* --- Linear Response */
        frame = gtk_frame_new("Linear Response");
        gtk_box_pack_start(GTK_BOX(page),frame,FALSE,FALSE,0);
/* create a table in the frame*/
        table = gtk_table_new(1, 5,FALSE);
        gtk_container_add(GTK_CONTAINER(frame), table);
/* 1st line */
	VASP_CHECK_TABLE(vasp_gui.loptics,vasp_calc.loptics,NULL,"LOPTICS",0,1,0,1);/*not calling anything*/
	VASP_CHECK_TABLE(vasp_gui.lepsilon,vasp_calc.lepsilon,NULL,"LEPSILON",1,2,0,1);/*not calling anything*/
	VASP_CHECK_TABLE(vasp_gui.lrpa,vasp_calc.lrpa,NULL,"LRPA",2,3,0,1);/*not calling anything*/
	VASP_CHECK_TABLE(vasp_gui.lnabla,vasp_calc.lnabla,NULL,"LNABLA",3,4,0,1);/*not calling anything*/
	VASP_ENTRY_TABLE(vasp_gui.cshift,vasp_calc.cshift,"%.2lf","CSHIFT=",4,5,0,1);
/* --- end frame */



/*-----------------*/
/* page 3 -> IONIC */
/*-----------------*/
        page = gtk_vbox_new(FALSE, _SPACE);
        label = gtk_label_new("IONIC");
        gtk_notebook_append_page(GTK_NOTEBOOK(notebook),page,label);
/* --- general */
        frame = gtk_frame_new("General");
        gtk_box_pack_start(GTK_BOX(page),frame,FALSE,FALSE,0);
/* create a table in the frame*/
        table = gtk_table_new(2, 5,FALSE);
        gtk_container_add(GTK_CONTAINER(frame), table);
//      gtk_container_set_border_width(GTK_CONTAINER(table), PANEL_SPACING);/*useful?*/
/* 1st line */
	VASP_COMBOBOX_TABLE(vasp_gui.ibrion,"IBRION=",0,1,0,1);
	VASP_COMBOBOX_ADD(vasp_gui.ibrion,"-1:Nothing");
	VASP_COMBOBOX_ADD(vasp_gui.ibrion,"0:Molecular Dynamics");
	VASP_COMBOBOX_ADD(vasp_gui.ibrion,"1:Quasi-Newton");
	VASP_COMBOBOX_ADD(vasp_gui.ibrion,"2:Conjugate Gradients");
	VASP_COMBOBOX_ADD(vasp_gui.ibrion,"3:Damped");
	VASP_COMBOBOX_ADD(vasp_gui.ibrion,"5:Finite Differences");
	VASP_COMBOBOX_ADD(vasp_gui.ibrion,"6:Finite Diff. with SYM");
	VASP_COMBOBOX_ADD(vasp_gui.ibrion,"7:Perturbation Theory");
	VASP_COMBOBOX_ADD(vasp_gui.ibrion,"8:Pert. Theory with SYM");
	VASP_COMBOBOX_ADD(vasp_gui.ibrion,"44:Transition State");
	VASP_ENTRY_TABLE(vasp_gui.nsw,vasp_calc.nsw,"%i","NSW=",1,2,0,1);
	VASP_ENTRY_TABLE(vasp_gui.ediffg,vasp_calc.ediffg,"%.2E","EDIFFG=",2,3,0,1);
	VASP_ENTRY_TABLE(vasp_gui.potim,vasp_calc.potim,"%.4lf","POTIM=",3,4,0,1);
	VASP_ENTRY_TABLE(vasp_gui.pstress,vasp_calc.pstress,"%.4lf","PSTRESS=",4,5,0,1);
/* 2nd line */
	VASP_COMBOBOX_TABLE(vasp_gui.isif,"ISIF=",0,1,1,2);
	VASP_COMBOBOX_ADD(vasp_gui.isif,"0:F_0_I_0_0");
	VASP_COMBOBOX_ADD(vasp_gui.isif,"1:F_P_I_0_0");
	VASP_COMBOBOX_ADD(vasp_gui.isif,"2:F_S_I_0_0");
	VASP_COMBOBOX_ADD(vasp_gui.isif,"3:F_S_I_S_V");
	VASP_COMBOBOX_ADD(vasp_gui.isif,"4:F_S_I_S_0");
	VASP_COMBOBOX_ADD(vasp_gui.isif,"5:F_S_0_S_0");
	VASP_COMBOBOX_ADD(vasp_gui.isif,"6:F_S_0_S_V");
	VASP_COMBOBOX_ADD(vasp_gui.isif,"7:F_S_0_0_V");
	VASP_CHECK_TABLE(vasp_gui.relax_ions,vasp_gui.rions,isif_toggle,"atom positions",1,2,1,2);
	VASP_CHECK_TABLE(vasp_gui.relax_shape,vasp_gui.rshape,isif_toggle,"cell shape",2,3,1,2);
	VASP_CHECK_TABLE(vasp_gui.relax_volume,vasp_gui.rvolume,isif_toggle,"cell volume",3,4,1,2);
	VASP_ENTRY_TABLE(vasp_gui.nfree,vasp_calc.nfree,"%i","NFREE=",4,5,1,2);
/* initialize */
	VASP_COMBOBOX_SETUP(vasp_gui.ibrion,0,vasp_ibrion_selected);
	VASP_COMBOBOX_SETUP(vasp_gui.isif,0,vasp_isif_selected);
/* --- end frame */
/* --- MD */
        frame = gtk_frame_new("Molecular Dynamics");
        gtk_box_pack_start(GTK_BOX(page),frame,FALSE,FALSE,0);
/* create a table in the frame*/
        table = gtk_table_new(2, 5,FALSE);
        gtk_container_add(GTK_CONTAINER(frame), table);
/* 1st line */
	//multicolumn label
	label = gtk_label_new("SET IBRION=0 FOR MD");
	gtk_table_attach_defaults(GTK_TABLE(table),label,0,1,0,2);
	VASP_ENTRY_TABLE(vasp_gui.tebeg,vasp_calc.tebeg,"%.1lf","TEBEG=",1,2,0,1);
	VASP_ENTRY_TABLE(vasp_gui.teend,vasp_calc.teend,"%.1lf","TEEND=",2,3,0,1);
	VASP_ENTRY_TABLE(vasp_gui.smass,vasp_calc.smass,"%.1lf","SMASS=",3,4,0,1);
	/*empty: col4*/
/* 2nd line */
	/*occ: col0*/
	VASP_ENTRY_TABLE(vasp_gui.nblock,vasp_calc.nblock,"%i","NBLOCK=",1,2,1,2);
	VASP_ENTRY_TABLE(vasp_gui.kblock,vasp_calc.kblock,"%i","KBLOCK=",2,3,1,2);
	VASP_ENTRY_TABLE(vasp_gui.npaco,vasp_calc.npaco,"%i","NPACO=",3,4,1,2);
	VASP_ENTRY_TABLE(vasp_gui.apaco,vasp_calc.apaco,"%.4lf","APACO=",4,5,1,2);
/* initialize sensitive widget*/	
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
        frame = gtk_frame_new("POSCAR & SYMMTERY");
        gtk_box_pack_start(GTK_BOX(page),frame,FALSE,FALSE,0);
/* create a table in the frame*/
        table = gtk_table_new(7, 5,FALSE);
        gtk_container_add(GTK_CONTAINER(frame), table);
//      gtk_container_set_border_width(GTK_CONTAINER(table), PANEL_SPACING);/*useful?*/
/* 1st line */
	VASP_COMBOBOX_TABLE(vasp_gui.poscar_free,"DYN:",0,1,0,1);
	VASP_COMBOBOX_ADD(vasp_gui.poscar_free,"All atom FIXED");
	VASP_COMBOBOX_ADD(vasp_gui.poscar_free,"All atom FREE");
	VASP_COMBOBOX_ADD(vasp_gui.poscar_free,"Manual selection");
	VASP_CHECK_TABLE(button,vasp_calc.poscar_sd,toggle_poscar_sd,"Sel. Dynamics",1,2,0,1);
	VASP_ENTRY_TABLE(vasp_gui.isym,vasp_calc.isym,"%i","ISYM=",2,3,0,1);
	VASP_ENTRY_TABLE(vasp_gui.sym_prec,vasp_calc.sym_prec,"%.2lE","_PREC=",3,4,0,1);
	VASP_ENTRY_TABLE(vasp_gui.poscar_a0,vasp_calc.poscar_a0,"%.4lf","A0=",4,5,0,1);
/* 2nd line */
	/*empty: col0 (TODO: SUPERCELL)*/
	VASP_CHECK_TABLE(vasp_gui.poscar_direct,vasp_calc.poscar_direct,toggle_poscar_direct,"Direct Coord.",1,2,1,2);
	VASP_ENTRY_TABLE(vasp_gui.poscar_ux,vasp_calc.poscar_ux,"%.4lf","UX=",2,3,1,2);
	VASP_ENTRY_TABLE(vasp_gui.poscar_uy,vasp_calc.poscar_uy,"%.4lf","UY=",3,4,1,2);
	VASP_ENTRY_TABLE(vasp_gui.poscar_uz,vasp_calc.poscar_uz,"%.4lf","UZ=",4,5,1,2);
/* 3rd line */
	/*empty: col0 (TODO: -X)*/
	/*empty: col1 (TODO: +X)*/
	VASP_ENTRY_TABLE(vasp_gui.poscar_vx,vasp_calc.poscar_vx,"%.4lf","VX=",2,3,2,3);
	VASP_ENTRY_TABLE(vasp_gui.poscar_vy,vasp_calc.poscar_vy,"%.4lf","VY=",3,4,2,3);
	VASP_ENTRY_TABLE(vasp_gui.poscar_vz,vasp_calc.poscar_vz,"%.4lf","VZ=",4,5,2,3);
/* 4th line */
	/*empty: col0 (TODO: -Y)*/
	/*empty: col1 (TODO: +Y)*/
	VASP_ENTRY_TABLE(vasp_gui.poscar_wx,vasp_calc.poscar_wx,"%.4lf","WX=",2,3,3,4);
	VASP_ENTRY_TABLE(vasp_gui.poscar_wy,vasp_calc.poscar_wy,"%.4lf","WY=",3,4,3,4);
	VASP_ENTRY_TABLE(vasp_gui.poscar_wz,vasp_calc.poscar_wz,"%.4lf","WZ=",4,5,3,4);
/* 5th line */
	/*empty: col0 (TODO: -Z)*/
	/*empty: col1 (TODO: +Z)*/
	VASP_CHECK_TABLE(vasp_gui.poscar_tx,vasp_calc.poscar_tx,NULL,"TX",2,3,4,5);/*not calling anything*/
	VASP_CHECK_TABLE(vasp_gui.poscar_ty,vasp_calc.poscar_ty,NULL,"TY",3,4,4,5);/*not calling anything*/
	VASP_CHECK_TABLE(vasp_gui.poscar_tz,vasp_calc.poscar_tz,NULL,"TZ",4,5,4,5);/*not calling anything*/
/* 6th line */
        /*this combo is special and should be explicitely defined*/
	hbox = gtk_hbox_new(FALSE, 0);
	label = gtk_label_new("ATOMS:");
	gtk_box_pack_start(GTK_BOX(hbox),label,FALSE,FALSE,0);
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
        gtk_combo_box_text_append_text(GTK_COMBO_BOX_TEXT (vasp_gui.poscar_atoms),"ADD atom");
	gtk_box_pack_start(GTK_BOX(hbox),vasp_gui.poscar_atoms,TRUE,TRUE,0);
        gtk_table_attach_defaults(GTK_TABLE(table),hbox,0,4,5,6);
        g_signal_connect(GTK_COMBO_BOX_TEXT(vasp_gui.poscar_atoms),"changed",GTK_SIGNAL_FUNC(vasp_poscar_atoms_selected),data);
	/*2 buttons*/
	hbox = gtk_hbox_new(FALSE, 0);
	VASP_NEW_SEPARATOR();
        gui_stock_button(GTK_STOCK_APPLY,vasp_atom_modified,NULL,hbox);
        gui_stock_button(GTK_STOCK_DELETE,vasp_atom_delete,NULL,hbox);
	gtk_table_attach_defaults(GTK_TABLE(table),hbox,4,5,5,6);
/* 7th line */
	VASP_ENTRY_TABLE(vasp_gui.poscar_index,vasp_calc.poscar_index,"%i","INDEX:",0,1,6,7);
	gtk_widget_set_sensitive(vasp_gui.poscar_index,FALSE);/*<- is changed only by selecting poscar_atoms*/
	VASP_ENTRY_TABLE(vasp_gui.poscar_symbol,vasp_calc.poscar_symbol,"%s","SYMBOL:",1,2,6,7);
	VASP_ENTRY_TABLE(vasp_gui.poscar_x,vasp_calc.poscar_x,"%.6lf","X=",2,3,6,7);
	VASP_ENTRY_TABLE(vasp_gui.poscar_y,vasp_calc.poscar_y,"%.6lf","Y=",3,4,6,7);
	VASP_ENTRY_TABLE(vasp_gui.poscar_z,vasp_calc.poscar_z,"%.6lf","Z=",4,5,6,7);
/* initialize sensitive widget*/
	VASP_COMBOBOX_SETUP(vasp_gui.poscar_free,1,vasp_poscar_free_selected);
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
/* (re)initialize */
	toggle_poscar_sd(NULL,NULL);
	gtk_combo_box_set_active(GTK_COMBO_BOX(vasp_gui.poscar_atoms),idx);/*sets "ADD atom" as active*/
	vasp_poscar_free_selected(vasp_gui.poscar_free,NULL);
/* --- end frame */

/*-------------------*/
/* page 4 -> KPOINTS */
/*-------------------*/
        page = gtk_vbox_new(FALSE, _SPACE);
        label = gtk_label_new("KPOINTS");
        gtk_notebook_append_page(GTK_NOTEBOOK(notebook),page,label);
/* --- from INCAR */
        frame = gtk_frame_new("from INCAR");
        gtk_box_pack_start(GTK_BOX(page),frame,FALSE,FALSE,0);
/* create a table in the frame*/
        table = gtk_table_new(1, 5,FALSE);
        gtk_container_add(GTK_CONTAINER(frame), table);
//      gtk_container_set_border_width(GTK_CONTAINER(table), PANEL_SPACING);/*useful?*/
/* 1st line */
	VASP_ENTRY_TABLE(vasp_gui.ismear_3,vasp_calc.ismear,"%i","ISMEAR=",0,1,0,1);
	VASP_CHECK_TABLE(vasp_gui.kgamma_3,vasp_calc.kgamma,NULL,"KGAMMA",1,2,0,1);/*not calling anything*/
	VASP_ENTRY_TABLE(vasp_gui.kspacing_3,vasp_calc.kspacing,"%.4f","KSPACING=",2,3,0,1);
	/*empty: col3*/
	/*empty: col4*/
/* --- end frame */
/* --- General */
        frame = gtk_frame_new("General");
        gtk_box_pack_start(GTK_BOX(page),frame,FALSE,FALSE,0);
/* create a table in the frame*/
        table = gtk_table_new(2, 5,FALSE);
        gtk_container_add(GTK_CONTAINER(frame), table);
//      gtk_container_set_border_width(GTK_CONTAINER(table), PANEL_SPACING);/*useful?*/
/* 1st line */
	VASP_COMBOBOX_TABLE(vasp_gui.kpoints_mode,"GEN:",0,1,0,1);
	VASP_COMBOBOX_ADD(vasp_gui.kpoints_mode,"Manual Entry");
	VASP_COMBOBOX_ADD(vasp_gui.kpoints_mode,"Line Mode");
	VASP_COMBOBOX_ADD(vasp_gui.kpoints_mode,"Automatic (M&P)");
	VASP_COMBOBOX_ADD(vasp_gui.kpoints_mode,"Gamma (M&P)");
	VASP_COMBOBOX_ADD(vasp_gui.kpoints_mode,"Monkhorst-Pack classic");
	VASP_COMBOBOX_ADD(vasp_gui.kpoints_mode,"Basis set definition");
	VASP_CHECK_TABLE(vasp_gui.kpoints_cart,vasp_calc.kpoints_cart,NULL,"Cartesian",1,2,0,1);/*not calling anything*/
	VASP_SPIN_TABLE(vasp_gui.kpoints_kx,vasp_calc.kpoints_kx,NULL,"KX",2,3,0,1);
	VASP_SPIN_TABLE(vasp_gui.kpoints_ky,vasp_calc.kpoints_ky,NULL,"KY",3,4,0,1);
	VASP_SPIN_TABLE(vasp_gui.kpoints_kz,vasp_calc.kpoints_kz,NULL,"KZ",4,5,0,1);
/* 2nd line */
	VASP_ENTRY_TABLE(vasp_gui.kpoints_nkpts,vasp_calc.kpoints_nkpts,"%i","NKPTS=",0,1,1,2);
	VASP_CHECK_TABLE(vasp_gui.have_tetra,vasp_calc.kpoints_tetra,toogle_tetra,"Tetrahedron",1,2,1,2);
	VASP_ENTRY_TABLE(vasp_gui.kpoints_sx,vasp_calc.kpoints_sx,"%.4lf","SX=",2,3,1,2);
	VASP_ENTRY_TABLE(vasp_gui.kpoints_sy,vasp_calc.kpoints_sy,"%.4lf","SY=",3,4,1,2);
	VASP_ENTRY_TABLE(vasp_gui.kpoints_sz,vasp_calc.kpoints_sz,"%.4lf","SZ=",4,5,1,2);
/* --- end frame */
/* --- Coordinate */
        vasp_gui.kpoints_coord = gtk_frame_new("Coordinate");
        gtk_box_pack_start(GTK_BOX(page),vasp_gui.kpoints_coord,FALSE,FALSE,0);
/* create a table in the frame*/
        table = gtk_table_new(2, 5,FALSE);
        gtk_container_add(GTK_CONTAINER(vasp_gui.kpoints_coord), table);
//      gtk_container_set_border_width(GTK_CONTAINER(vasp_gui.kpoints_coord), PANEL_SPACING);/*useful?*/
/* 1st line */
	VASP_COMBOBOX_TABLE(vasp_gui.kpoints_kpts,"KPOINT:",0,4,0,1);
	VASP_COMBOBOX_ADD(vasp_gui.kpoints_kpts,"ADD kpoint");
        /*2 buttons*/
        hbox = gtk_hbox_new(FALSE, 0);
	VASP_NEW_SEPARATOR();
        gui_stock_button(GTK_STOCK_APPLY,vasp_kpoint_modified,NULL,hbox);
        gui_stock_button(GTK_STOCK_DELETE,vasp_kpoint_delete,NULL,hbox);
        gtk_table_attach_defaults(GTK_TABLE(table),hbox,4,5,0,1);
/* 2nd line */
        VASP_ENTRY_TABLE(vasp_gui.kpoints_i,vasp_calc.kpoints_i,"%i","INDEX:",0,1,1,2);
        gtk_widget_set_sensitive(vasp_gui.kpoints_i,FALSE);/*<- is changed only by selecting kpoints_kpts*/
	VASP_ENTRY_TABLE(vasp_gui.kpoints_x,vasp_calc.kpoints_x,"%.4lf","X=",1,2,1,2);
	VASP_ENTRY_TABLE(vasp_gui.kpoints_y,vasp_calc.kpoints_y,"%.4lf","Y=",2,3,1,2);
	VASP_ENTRY_TABLE(vasp_gui.kpoints_z,vasp_calc.kpoints_z,"%.4lf","Z=",3,4,1,2);
	VASP_ENTRY_TABLE(vasp_gui.kpoints_w,vasp_calc.kpoints_w,"%.4lf","W=",4,5,1,2);
/* --- end frame */
/* --- Tetrahedron */
        vasp_gui.kpoints_tetra = gtk_frame_new("Tetrahedron");
        gtk_box_pack_start(GTK_BOX(page),vasp_gui.kpoints_tetra,FALSE,FALSE,0);
/* create a table in the frame*/
        table = gtk_table_new(3, 6,FALSE);
        gtk_container_add(GTK_CONTAINER(vasp_gui.kpoints_tetra), table);
//      gtk_container_set_border_width(GTK_CONTAINER(vasp_gui.kpoints_coord), PANEL_SPACING);/*useful?*/
/* 1st line */
	VASP_ENTRY_TABLE(vasp_gui.tetra_total,vasp_calc.tetra_total,"%i","TOTAL=",0,1,0,1);
	label = gtk_label_new("EXACT VOLUME:");
	gtk_table_attach_defaults(GTK_TABLE(table),label,1,2,0,1);
	VASP_ENTRY_TABLE(vasp_gui.tetra_volume,vasp_calc.tetra_volume,"%.6lf","V=",2,3,0,1);
	label = gtk_label_new("for (A,B,C,D) tetrahedron");
	gtk_table_attach_defaults(GTK_TABLE(table),label,3,6,0,1);
/* 2nd line */
	VASP_COMBOBOX_TABLE(vasp_gui.tetra,"TETRAHEDRON:",0,5,1,2);
	VASP_COMBOBOX_ADD(vasp_gui.tetra,"ADD tetrahedron");
	/*2 buttons*/
	hbox = gtk_hbox_new(FALSE, 0);
	VASP_NEW_SEPARATOR();
	gui_stock_button(GTK_STOCK_APPLY,vasp_tetra_modified,NULL,hbox);
	gui_stock_button(GTK_STOCK_DELETE,vasp_tetra_delete,NULL,hbox);
	gtk_table_attach_defaults(GTK_TABLE(table),hbox,5,6,1,2);
/* 3rd line */
	VASP_ENTRY_TABLE(vasp_gui.tetra_i,vasp_calc.tetra_i,"%i","INDEX:",0,1,2,3);
	gtk_widget_set_sensitive(vasp_gui.tetra_i,FALSE);/*<- is changed only by selecting tetra*/
	VASP_ENTRY_TABLE(vasp_gui.tetra_w,vasp_calc.tetra_w,"%.4lf","Degen. W=",1,2,2,3);
	VASP_ENTRY_TABLE(vasp_gui.tetra_a,vasp_calc.tetra_a,"%i","PtA:",2,3,2,3);
	VASP_ENTRY_TABLE(vasp_gui.tetra_b,vasp_calc.tetra_b,"%i","PtB:",3,4,2,3);
	VASP_ENTRY_TABLE(vasp_gui.tetra_c,vasp_calc.tetra_c,"%i","PtC:",4,5,2,3);
	VASP_ENTRY_TABLE(vasp_gui.tetra_d,vasp_calc.tetra_d,"%i","PtD:",5,6,2,3);
/* initialize */
        g_signal_connect(GTK_OBJECT(GTK_ENTRY(vasp_gui.ismear_3)),"changed",GTK_SIGNAL_FUNC(ismear_changed),data);
        g_signal_connect(GTK_OBJECT(GTK_ENTRY(vasp_gui.kspacing_3)),"changed",GTK_SIGNAL_FUNC(kspacing_changed),data);
	VASP_COMBOBOX_SETUP(vasp_gui.kpoints_mode,2,vasp_kpoints_mode_selected);
	VASP_COMBOBOX_SETUP(vasp_gui.kpoints_kpts,0,vasp_kpoints_kpts_selected);
	VASP_COMBOBOX_SETUP(vasp_gui.tetra,0,vasp_tetra_selected);
	toogle_tetra();
/* --- end frame */
/* TODO: import KPOINTS/IBZKPT file*/

/*------------------*/
/* page 5 -> POTCAR */
/*------------------*/
        page = gtk_vbox_new(FALSE, _SPACE);
        label = gtk_label_new("POTCAR");
        gtk_notebook_append_page(GTK_NOTEBOOK(notebook),page,label);
/* --- general */
        frame = gtk_frame_new("General");
        gtk_box_pack_start(GTK_BOX(page),frame,FALSE,FALSE,0);
/* create a table in the frame*/
        table = gtk_table_new(2, 1,FALSE);
        gtk_container_add(GTK_CONTAINER(frame), table);
//      gtk_container_set_border_width(GTK_CONTAINER(vasp_gui.kpoints_coord), PANEL_SPACING);/*useful?*/
/* 1st line */
	VASP_TEXT_TABLE(vasp_gui.poscar_species,vasp_calc.species_symbols,"SPECIES:",0,1,0,1);
	gtk_widget_set_sensitive(vasp_gui.poscar_species,FALSE);
/* 2nd line */
	label = gtk_label_new("(each species will require a separate POTCAR information)");
        gtk_table_attach_defaults(GTK_TABLE(table),label,0,1,1,2);
/* --- end frame */
/* --- get POTCAR information */
        frame = gtk_frame_new("get POTCAR information");
        gtk_box_pack_start(GTK_BOX(page),frame,FALSE,FALSE,0);
/* create a table in the frame*/
        table = gtk_table_new(3, 4,FALSE);
        gtk_container_add(GTK_CONTAINER(frame), table);
//      gtk_container_set_border_width(GTK_CONTAINER(vasp_gui.kpoints_coord), PANEL_SPACING);/*useful?*/
/* 1st line */
	/* 2 lines radiobutton */
	vbox = gtk_vbox_new(FALSE,0);
	new_radio_group(0,vbox,FF);
	vasp_gui.potcar_select_file = add_radio_button("Use POTCAR file",(gpointer)toggle_potcar_file,NULL);
	vasp_gui.potcar_select_folder = add_radio_button("Use POTCAR path",(gpointer)toggle_potcar_folder,NULL);
	gtk_table_attach_defaults(GTK_TABLE(table),vbox,0,1,0,2);
	VASP_TEXT_TABLE(vasp_gui.potcar_file,vasp_calc.potcar_file,"FILE:",1,3,0,1);
	VASP_BUTTON_TABLE(vasp_gui.potcar_file_button,GTK_STOCK_OPEN,load_potcar_file_dialog,3,4,0,1);
/* 2nd line */
	VASP_TEXT_TABLE(vasp_gui.potcar_folder,vasp_calc.potcar_folder,"PATH:",1,3,1,2);
	VASP_BUTTON_TABLE(vasp_gui.potcar_folder_button,GTK_STOCK_OPEN,load_potcar_folder_dialog,3,4,1,2);
/* 3rd line */
        VASP_COMBOBOX_TABLE(vasp_gui.species,"SPECIES:",0,2,2,3);
        VASP_COMBOBOX_ADD(vasp_gui.species,"ADD SPECIES");
	VASP_COMBOBOX_TABLE(vasp_gui.species_flavor,"FLAVOR:",2,3,2,3);
	VASP_BUTTON_TABLE(vasp_gui.species_button,GTK_STOCK_APPLY,vasp_gui.species_flavor,3,4,2,3);
/* initialize */
	toggle_potcar_file(NULL,NULL);
	gtk_widget_set_sensitive(vasp_gui.potcar_select_folder,FALSE);/*TODO: not ready yet!*/
/* --- end frame */
/* --- POTCAR results */
        frame = gtk_frame_new("POTCAR results");
        gtk_box_pack_start(GTK_BOX(page),frame,FALSE,FALSE,0);
/* create a table in the frame*/
        table = gtk_table_new(2, 2,FALSE);
        gtk_container_add(GTK_CONTAINER(frame), table);
//      gtk_container_set_border_width(GTK_CONTAINER(vasp_gui.kpoints_coord), PANEL_SPACING);/*useful?*/
/* 1st line */
	//multicolumn label
        label = gtk_label_new("DETECTED");
        gtk_table_attach_defaults(GTK_TABLE(table),label,0,1,0,2);
	VASP_TEXT_TABLE(vasp_gui.potcar_species,vasp_calc.potcar_species,"SPECIES:",1,2,0,1);
	gtk_widget_set_sensitive(vasp_gui.potcar_species,FALSE);
/* 2nd line */
	/*occ: col0*/
	VASP_TEXT_TABLE(vasp_gui.potcar_species_flavor,vasp_calc.potcar_species_flavor,"FLAVOR:",1,2,1,2);
	gtk_widget_set_sensitive(vasp_gui.potcar_species_flavor,FALSE);

/*----------------*/
/* page 6 -> EXEC */
/*----------------*/
        page = gtk_vbox_new(FALSE, _SPACE);
        label = gtk_label_new("EXEC");
        gtk_notebook_append_page(GTK_NOTEBOOK(notebook),page,label);
/* --- general */
        frame = gtk_frame_new("General");
        gtk_box_pack_start(GTK_BOX(page),frame,FALSE,FALSE,0);
/* create a table in the frame*/
        table = gtk_table_new(3, 6,FALSE);
        gtk_container_add(GTK_CONTAINER(frame), table);
//      gtk_container_set_border_width(GTK_CONTAINER(vasp_gui.kpoints_coord), PANEL_SPACING);/*useful?*/
/* 1st line */
	VASP_LABEL_TABLE("VASP EXEC:",0,1,0,1);
	VASP_TEXT_TABLE(vasp_gui.job_vasp_exe,vasp_calc.job_vasp_exe,"FILE=",1,5,0,1);
	VASP_BUTTON_TABLE(button,GTK_STOCK_OPEN,load_vasp_exe_dialog,5,6,0,1);
/* 2nd line */
	VASP_LABEL_TABLE("CALCUL PATH:",0,1,1,2);
	VASP_TEXT_TABLE(vasp_gui.job_path,vasp_calc.job_path,"PATH=",1,5,1,2);
	VASP_BUTTON_TABLE(button,GTK_STOCK_OPEN,load_path_dialog,5,6,1,2);
/* 3rd line */
	VASP_LABEL_TABLE("OUTPUT:",0,1,2,3);
	VASP_CHECK_TABLE(button,vasp_calc.lwave,NULL,"LWAVE",1,2,2,3);/*not calling anything*/
	VASP_CHECK_TABLE(button,vasp_calc.lcharg,NULL,"LCHARG",2,3,2,3);/*not calling anything*/
	VASP_CHECK_TABLE(button,vasp_calc.lvtot,NULL,"LVTOT",3,4,2,3);/*not calling anything*/
	VASP_CHECK_TABLE(button,vasp_calc.lvhar,NULL,"LVHAR",4,5,2,3);/*not calling anything*/
	VASP_CHECK_TABLE(button,vasp_calc.lelf,NULL,"LELF",5,6,2,3);/*not calling anything*/
/* --- end frame */
/* --- PARALLEL */
        frame = gtk_frame_new("Parallel Optimization");
        gtk_box_pack_start(GTK_BOX(page),frame,FALSE,FALSE,0);
/* create a table in the frame*/
	table = gtk_table_new(2, 6,FALSE);
	gtk_container_add(GTK_CONTAINER(frame), table);
//      gtk_container_set_border_width(GTK_CONTAINER(vasp_gui.kpoints_coord), PANEL_SPACING);/*useful?*/
/* 1st line */
	VASP_LABEL_TABLE("MPIRUN:",0,1,0,1);
	VASP_TEXT_TABLE(vasp_gui.job_mpirun,vasp_calc.job_mpirun,"FILE=",1,5,0,1);
	VASP_BUTTON_TABLE(button,GTK_STOCK_OPEN,load_mpirun_dialog,5,6,0,1);
/* 2nd line */
	VASP_SPIN_TABLE(vasp_gui.job_nproc,vasp_calc.job_nproc,parallel_eq,"NP=",0,1,1,2);
	VASP_SPIN_TABLE(vasp_gui.ncore,vasp_calc.ncore,parallel_eq,"NCORE=",1,2,1,2);
	VASP_SPIN_TABLE(vasp_gui.kpar,vasp_calc.kpar,parallel_eq,"KPAR=",2,3,1,2);
	VASP_CHECK_TABLE(button,vasp_calc.lplane,NULL,"LPLANE",3,4,1,2);/*not calling anything*/
	VASP_CHECK_TABLE(button,vasp_calc.lscalu,NULL,"LSCALU",4,5,1,2);/*not calling anything*/
	VASP_CHECK_TABLE(button,vasp_calc.lscalapack,NULL,"LSCALAPACK",5,6,1,2);/*not calling anything*/
/* --- end frame */
/* --- DISTANT */
        frame = gtk_frame_new("Distant calculation");
        gtk_box_pack_start(GTK_BOX(page),frame,FALSE,FALSE,0);
/* create a table in the frame*/
        table = gtk_table_new(3, 3,FALSE);
        gtk_container_add(GTK_CONTAINER(frame), table);
//      gtk_container_set_border_width(GTK_CONTAINER(vasp_gui.kpoints_coord), PANEL_SPACING);/*useful?*/
/* 1st line */
        VASP_LABEL_TABLE("REMOTE HOSTS:",0,1,0,1);
        label = gtk_label_new("UNDER CONSTRUCTION");
        gtk_table_attach_defaults(GTK_TABLE(table),label,1,3,0,1);
/* 2nd line */
	VASP_LABEL_TABLE("REMOTE PATH:",0,1,1,2);
	label = gtk_label_new("UNDER CONSTRUCTION");
	gtk_table_attach_defaults(GTK_TABLE(table),label,1,3,1,2);
/* 3rd line */
	VASP_LABEL_TABLE("SUBMIT SCRIPT:",0,1,2,3);
	label = gtk_label_new("UNDER CONSTRUCTION");
	gtk_table_attach_defaults(GTK_TABLE(table),label,1,3,2,3);

/* --- Outside of notebook */
	frame = gtk_frame_new(NULL);
	gtk_box_pack_start(GTK_BOX(GTK_DIALOG(vasp_gui.window)->vbox), frame, FALSE, FALSE, 0);
	vbox = gtk_vbox_new(FALSE, _SPACE);
	gtk_container_add(GTK_CONTAINER(frame), vbox);
/* Action buttons */
	vasp_gui.spinner=gtk_spinner_new();
	gtk_box_pack_start(GTK_BOX(GTK_DIALOG(vasp_gui.window)->action_area), vasp_gui.spinner, FALSE, FALSE, 0);
	vasp_gui.button_save=gui_stock_button(GTK_STOCK_SAVE, save_vasp_calc, NULL, GTK_DIALOG(vasp_gui.window)->action_area);
	vasp_gui.button_exec=gui_stock_button(GTK_STOCK_EXECUTE, exec_vasp_calc, NULL, GTK_DIALOG(vasp_gui.window)->action_area);
	gui_stock_button(GTK_STOCK_CLOSE, quit_vasp_calc, dialog, GTK_DIALOG(vasp_gui.window)->action_area);
/* connect to signals */
        g_signal_connect(GTK_NOTEBOOK(notebook),"switch-page",GTK_SIGNAL_FUNC(vasp_calc_page_switch),NULL);

/* all done */
        gtk_widget_show_all(vasp_gui.window);
	gtk_widget_hide(vasp_gui.spinner);
        sysenv.refresh_dialog=TRUE;
}

