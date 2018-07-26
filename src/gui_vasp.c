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
#include "file_vasp.h"
#include "gui_defs.h"
#include "gui_vasp.h"

extern struct sysenv_pak sysenv;

/*globals variables*/
struct vasp_calc_gui vasp_gui;

/*************************/
/* initialize gui values */
/*************************/
void gui_vasp_init(){
	vasp_gui.cur_page=VASP_PAGE_SIMPLIFIED;
	vasp_gui.have_potcar_folder=FALSE;
	vasp_gui.poscar_dirty=TRUE;
	if(!vasp_gui.have_xml){
		vasp_gui.rions=TRUE;
		vasp_gui.rshape=FALSE;
		vasp_gui.rvolume=FALSE;
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
	vasprun_update(filename,&vasp_gui.calc);/*<- update everything according to the new vasprun.xml*/
	g_free (filename);
	vasp_gui_refresh();/*refresh GUI with new vasp_gui.calc*/
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
/* refresh the GUI with value from vasp_gui.calc */
/*********************************************/
void vasp_gui_refresh(){
/*do not refresh the simple interface, though*/
        switch(vasp_gui.calc.prec){
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
	gtk_entry_set_text(GTK_ENTRY(vasp_gui.encut),g_strdup_printf("%.2lf",vasp_gui.calc.encut));
	gtk_entry_set_text(GTK_ENTRY(vasp_gui.enaug),g_strdup_printf("%.2lf",vasp_gui.calc.enaug));
	gtk_entry_set_text(GTK_ENTRY(vasp_gui.ediff),g_strdup_printf("%.2lE",vasp_gui.calc.ediff));
	switch(vasp_gui.calc.algo){
	case VA_IALGO:
		gtk_combo_box_set_active(GTK_COMBO_BOX(vasp_gui.algo),1);
		switch(vasp_gui.calc.ialgo){
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
	gtk_entry_set_text(GTK_ENTRY(vasp_gui.nsim),g_strdup_printf("%i",vasp_gui.calc.nsim));
	gtk_entry_set_text(GTK_ENTRY(vasp_gui.vtime),g_strdup_printf("%.4lf",vasp_gui.calc.vtime));
	gtk_entry_set_text(GTK_ENTRY(vasp_gui.iwavpr),g_strdup_printf("%i",vasp_gui.calc.iwavpr));
	//auto already sync
	gtk_entry_set_text(GTK_ENTRY(vasp_gui.nbands),g_strdup_printf("%i",vasp_gui.calc.nbands));
	gtk_entry_set_text(GTK_ENTRY(vasp_gui.nelect),g_strdup_printf("%.4lf",vasp_gui.calc.nelect));
	//iniwav already sync
	gtk_entry_set_text(GTK_ENTRY(vasp_gui.istart),g_strdup_printf("%i",vasp_gui.calc.istart));
	gtk_entry_set_text(GTK_ENTRY(vasp_gui.icharg),g_strdup_printf("%i",vasp_gui.calc.icharg));
	/*mixer*/
	gtk_entry_set_text(GTK_ENTRY(vasp_gui.nelm),g_strdup_printf("%i",vasp_gui.calc.nelm));
	gtk_entry_set_text(GTK_ENTRY(vasp_gui.nelmdl),g_strdup_printf("%i",vasp_gui.calc.nelmdl));
	gtk_entry_set_text(GTK_ENTRY(vasp_gui.nelmin),g_strdup_printf("%i",vasp_gui.calc.nelmin));
	if(!vasp_gui.calc.auto_mixer){
		switch(vasp_gui.calc.imix){
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
		switch(vasp_gui.calc.inimix){
		case VIM_0:
			gtk_combo_box_set_active(GTK_COMBO_BOX(vasp_gui.inimix),0);break;
		case VIM_1:
		default:
			gtk_combo_box_set_active(GTK_COMBO_BOX(vasp_gui.inimix),1);
		}
		switch(vasp_gui.calc.mixpre){
		case VMP_0:
			gtk_combo_box_set_active(GTK_COMBO_BOX(vasp_gui.mixpre),0);break;
		case VMP_2:
			gtk_combo_box_set_active(GTK_COMBO_BOX(vasp_gui.mixpre),2);break;
		case VMP_1:
		default:
			gtk_combo_box_set_active(GTK_COMBO_BOX(vasp_gui.mixpre),1);
		}
        }
	switch(vasp_gui.calc.lreal){
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
	if(vasp_gui.calc.ropt!=NULL) gtk_entry_set_text(GTK_ENTRY(vasp_gui.ropt),g_strdup_printf("%s",vasp_gui.calc.ropt));
	//addgrid already sync
	gtk_entry_set_text(GTK_ENTRY(vasp_gui.lmaxmix),g_strdup_printf("%i",vasp_gui.calc.lmaxmix));
	gtk_entry_set_text(GTK_ENTRY(vasp_gui.lmaxpaw),g_strdup_printf("%i",vasp_gui.calc.lmaxpaw));
	if(!vasp_gui.calc.auto_grid){
		gtk_entry_set_text(GTK_ENTRY(vasp_gui.ngx),g_strdup_printf("%i",vasp_gui.calc.ngx));
		gtk_entry_set_text(GTK_ENTRY(vasp_gui.ngy),g_strdup_printf("%i",vasp_gui.calc.ngy));
		gtk_entry_set_text(GTK_ENTRY(vasp_gui.ngz),g_strdup_printf("%i",vasp_gui.calc.ngz));
		gtk_entry_set_text(GTK_ENTRY(vasp_gui.ngxf),g_strdup_printf("%i",vasp_gui.calc.ngxf));
		gtk_entry_set_text(GTK_ENTRY(vasp_gui.ngyf),g_strdup_printf("%i",vasp_gui.calc.ngyf));
		gtk_entry_set_text(GTK_ENTRY(vasp_gui.ngzf),g_strdup_printf("%i",vasp_gui.calc.ngzf));
	}
	gtk_entry_set_text(GTK_ENTRY(vasp_gui.ismear),g_strdup_printf("%i",vasp_gui.calc.ismear));
	gtk_entry_set_text(GTK_ENTRY(vasp_gui.sigma),g_strdup_printf("%.4lf",vasp_gui.calc.sigma));
	gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(vasp_gui.kgamma),vasp_gui.calc.kgamma);
	gtk_entry_set_text(GTK_ENTRY(vasp_gui.kspacing),g_strdup_printf("%.4lf",vasp_gui.calc.kspacing));
	if(vasp_gui.calc.fermwe!=NULL) gtk_entry_set_text(GTK_ENTRY(vasp_gui.fermwe),g_strdup_printf("%s",vasp_gui.calc.fermwe));
	if(vasp_gui.calc.fermdo!=NULL) gtk_entry_set_text(GTK_ENTRY(vasp_gui.fermdo),g_strdup_printf("%s",vasp_gui.calc.fermdo));
	//ispin already sync
	//non_collinear already sync
	//lsorbit already sync
	//gga_compat already sync
	gtk_entry_set_text(GTK_ENTRY(vasp_gui.nupdown),g_strdup_printf("%.1lf",vasp_gui.calc.nupdown));
	if(vasp_gui.calc.magmom!=NULL) gtk_entry_set_text(GTK_ENTRY(vasp_gui.magmom),g_strdup_printf("%s",vasp_gui.calc.magmom));
	if(vasp_gui.calc.saxis!=NULL) gtk_entry_set_text(GTK_ENTRY(vasp_gui.saxis),g_strdup_printf("%s",vasp_gui.calc.saxis));
	switch(vasp_gui.calc.gga){
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
	switch (vasp_gui.calc.mgga){
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
	gtk_entry_set_text(GTK_ENTRY(vasp_gui.lmaxtau),g_strdup_printf("%i",vasp_gui.calc.lmaxtau));
	if(vasp_gui.calc.cmbj!=NULL) gtk_entry_set_text(GTK_ENTRY(vasp_gui.cmbj),g_strdup_printf("%s",vasp_gui.calc.cmbj));
	gtk_entry_set_text(GTK_ENTRY(vasp_gui.cmbja),g_strdup_printf("%.4lf",vasp_gui.calc.cmbja));
	gtk_entry_set_text(GTK_ENTRY(vasp_gui.cmbjb),g_strdup_printf("%.4lf",vasp_gui.calc.cmbjb));
	if(vasp_gui.calc.ldau){
		switch(vasp_gui.calc.ldau_type){
		case VU_1:
			gtk_combo_box_set_active(GTK_COMBO_BOX(vasp_gui.ldau),0);break;
		case VU_4:
			gtk_combo_box_set_active(GTK_COMBO_BOX(vasp_gui.ldau),2);break;
		case VU_2:
		default:
			gtk_combo_box_set_active(GTK_COMBO_BOX(vasp_gui.ldau),1);
		}
		switch(vasp_gui.calc.ldau_output){
		case VUO_1:
			gtk_combo_box_set_active(GTK_COMBO_BOX(vasp_gui.ldau_print),1);break;
		case VUO_2:
			gtk_combo_box_set_active(GTK_COMBO_BOX(vasp_gui.ldau_print),2);break;
		case VUO_0:
		default:
			gtk_combo_box_set_active(GTK_COMBO_BOX(vasp_gui.ldau_print),0);
		}
		if(vasp_gui.calc.ldaul!=NULL) gtk_entry_set_text(GTK_ENTRY(vasp_gui.ldaul),g_strdup_printf("%s",vasp_gui.calc.ldaul));
		if(vasp_gui.calc.ldauu!=NULL) gtk_entry_set_text(GTK_ENTRY(vasp_gui.ldauu),g_strdup_printf("%s",vasp_gui.calc.ldauu));
		if(vasp_gui.calc.ldauj!=NULL) gtk_entry_set_text(GTK_ENTRY(vasp_gui.ldauj),g_strdup_printf("%s",vasp_gui.calc.ldauj));
	}
	switch(vasp_gui.calc.idipol){
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
	gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(vasp_gui.ldipol),vasp_gui.calc.ldipol);
	gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(vasp_gui.lmono),vasp_gui.calc.lmono);
	if(vasp_gui.calc.dipol!=NULL) gtk_entry_set_text(GTK_ENTRY(vasp_gui.dipol),g_strdup_printf("%s",vasp_gui.calc.dipol));
	gtk_entry_set_text(GTK_ENTRY(vasp_gui.epsilon),g_strdup_printf("%.4lf",vasp_gui.calc.epsilon));
	gtk_entry_set_text(GTK_ENTRY(vasp_gui.efield),g_strdup_printf("%.4lf",vasp_gui.calc.efield));
	switch(vasp_gui.calc.lorbit){
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
	gtk_entry_set_text(GTK_ENTRY(vasp_gui.nedos),g_strdup_printf("%i",vasp_gui.calc.nedos));
	gtk_entry_set_text(GTK_ENTRY(vasp_gui.emin),g_strdup_printf("%.2lf",vasp_gui.calc.emin));
	gtk_entry_set_text(GTK_ENTRY(vasp_gui.emax),g_strdup_printf("%.2lf",vasp_gui.calc.emax));
	gtk_entry_set_text(GTK_ENTRY(vasp_gui.efermi),g_strdup_printf("%.2lf",vasp_gui.calc.efermi));
	gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(vasp_gui.have_paw),vasp_gui.calc.have_paw);
	if(vasp_gui.calc.rwigs!=NULL) gtk_entry_set_text(GTK_ENTRY(vasp_gui.rwigs),g_strdup_printf("%s",vasp_gui.calc.rwigs));
	gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(vasp_gui.loptics),vasp_gui.calc.loptics);
	gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(vasp_gui.lepsilon),vasp_gui.calc.lepsilon);
	gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(vasp_gui.lrpa),vasp_gui.calc.lrpa);
	gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(vasp_gui.lnabla),vasp_gui.calc.lnabla);
	gtk_entry_set_text(GTK_ENTRY(vasp_gui.cshift),g_strdup_printf("%.2lf",vasp_gui.calc.cshift));
	switch(vasp_gui.calc.ibrion){
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
	gtk_combo_box_set_active(GTK_COMBO_BOX(vasp_gui.isif),vasp_gui.calc.isif);
	//relax_ions already sync
	//relax_shape already sync
	//relax_volume already sync
	gtk_entry_set_text(GTK_ENTRY(vasp_gui.nsw),g_strdup_printf("%i",vasp_gui.calc.nsw));
	gtk_entry_set_text(GTK_ENTRY(vasp_gui.ediffg),g_strdup_printf("%.2lE",vasp_gui.calc.ediffg));
	gtk_entry_set_text(GTK_ENTRY(vasp_gui.potim),g_strdup_printf("%.4lf",vasp_gui.calc.potim));
	gtk_entry_set_text(GTK_ENTRY(vasp_gui.pstress),g_strdup_printf("%.4lf",vasp_gui.calc.pstress));
	gtk_entry_set_text(GTK_ENTRY(vasp_gui.nfree),g_strdup_printf("%i",vasp_gui.calc.nfree));
	if(vasp_gui.calc.ibrion==0){
		gtk_entry_set_text(GTK_ENTRY(vasp_gui.tebeg),g_strdup_printf("%.1lf",vasp_gui.calc.tebeg));
		gtk_entry_set_text(GTK_ENTRY(vasp_gui.teend),g_strdup_printf("%.1lf",vasp_gui.calc.teend));
		gtk_entry_set_text(GTK_ENTRY(vasp_gui.smass),g_strdup_printf("%.1lf",vasp_gui.calc.smass));
		gtk_entry_set_text(GTK_ENTRY(vasp_gui.nblock),g_strdup_printf("%i",vasp_gui.calc.nblock));
		gtk_entry_set_text(GTK_ENTRY(vasp_gui.kblock),g_strdup_printf("%i",vasp_gui.calc.kblock));
		gtk_entry_set_text(GTK_ENTRY(vasp_gui.npaco),g_strdup_printf("%i",vasp_gui.calc.npaco));
		gtk_entry_set_text(GTK_ENTRY(vasp_gui.apaco),g_strdup_printf("%.4lf",vasp_gui.calc.apaco));
	}
	/*POSCAR*/
	gtk_entry_set_text(GTK_ENTRY(vasp_gui.poscar_species),g_strdup_printf("%s",vasp_gui.calc.species_symbols));
	gtk_entry_set_text(GTK_ENTRY(vasp_gui.simple_poscar),g_strdup_printf("%s",vasp_gui.calc.species_symbols));
	/*KPOINTS*/
	gtk_entry_set_text(GTK_ENTRY(vasp_gui.ismear_3),g_strdup_printf("%i",vasp_gui.calc.ismear));
	gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(vasp_gui.kgamma_3),vasp_gui.calc.kgamma);
	gtk_entry_set_text(GTK_ENTRY(vasp_gui.kspacing_3),g_strdup_printf("%.4lf",vasp_gui.calc.kspacing));
	gtk_entry_set_text(GTK_ENTRY(vasp_gui.kpoints_nkpts),g_strdup_printf("%i",vasp_gui.calc.kpoints_nkpts));
	gtk_combo_box_set_active(GTK_COMBO_BOX(vasp_gui.kpoints_mode),vasp_gui.calc.kpoints_mode);
	switch(vasp_gui.calc.kpoints_mode){
	case VKP_MAN:
		/*we do not update the kpoint list (here)*/
		break;
	case VKP_LINE:
		/*we do not update the kpoint list (here)*/
		break;
	case VKP_GAMMA:
		gtk_spin_button_set_value(GTK_SPIN_BUTTON(vasp_gui.kpoints_kx),vasp_gui.calc.kpoints_kx);
		gtk_spin_button_set_value(GTK_SPIN_BUTTON(vasp_gui.kpoints_ky),vasp_gui.calc.kpoints_ky);
		gtk_spin_button_set_value(GTK_SPIN_BUTTON(vasp_gui.kpoints_kz),vasp_gui.calc.kpoints_kz);
		gtk_entry_set_text(GTK_ENTRY(vasp_gui.kpoints_sx),g_strdup_printf("%.4lf",vasp_gui.calc.kpoints_sx));
		gtk_entry_set_text(GTK_ENTRY(vasp_gui.kpoints_sy),g_strdup_printf("%.4lf",vasp_gui.calc.kpoints_sy));
		gtk_entry_set_text(GTK_ENTRY(vasp_gui.kpoints_sz),g_strdup_printf("%.4lf",vasp_gui.calc.kpoints_sz));
		break;
	case VKP_MP:
		gtk_spin_button_set_value(GTK_SPIN_BUTTON(vasp_gui.kpoints_kx),vasp_gui.calc.kpoints_kx);
		gtk_spin_button_set_value(GTK_SPIN_BUTTON(vasp_gui.kpoints_ky),vasp_gui.calc.kpoints_ky);
		gtk_spin_button_set_value(GTK_SPIN_BUTTON(vasp_gui.kpoints_kz),vasp_gui.calc.kpoints_kz);
		gtk_entry_set_text(GTK_ENTRY(vasp_gui.kpoints_sx),g_strdup_printf("%.4lf",vasp_gui.calc.kpoints_sx));
		gtk_entry_set_text(GTK_ENTRY(vasp_gui.kpoints_sy),g_strdup_printf("%.4lf",vasp_gui.calc.kpoints_sy));
		gtk_entry_set_text(GTK_ENTRY(vasp_gui.kpoints_sz),g_strdup_printf("%.4lf",vasp_gui.calc.kpoints_sz));
		break;
	case VKP_BASIS:
		/*we do not update the basis definition list (here)*/
		break;
	case VKP_AUTO:
	default:
		gtk_spin_button_set_value(GTK_SPIN_BUTTON(vasp_gui.kpoints_kx),vasp_gui.calc.kpoints_kx);
	}
	/*skip POTCAR*/
	/*skip EXEC*/

/*
        gtk_entry_set_text(GTK_ENTRY(vasp_gui.ismear),g_strdup_printf("%i",vasp_gui.calc.ismear));
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
	gtk_text_buffer_set_text(vasp_gui.simple_message_buff,"Simple interface started.\n",-1);
	/*1st take care of the obvious*/
	if((vasp_gui.simple_rgeom)&&(calcul<3)){
		_OUT("FAIL: ROUGHT GEOMETRY: OPTIMIZE GEOMETRY FIRST!\n");
		return;
	}
	if(dim==0){
		_OUT("ATOM/MOLECULE in a box, setting only gamma point\n");
		vasp_gui.calc.kpoints_mode=VKP_GAMMA;
		vasp_gui.calc.kpoints_kx=1.;
		vasp_gui.calc.kpoints_ky=1.;
		vasp_gui.calc.kpoints_kz=1.;
	}else{
		vasp_gui.calc.kpoints_mode=VKP_AUTO;
		vasp_gui.calc.kpoints_kx = (gdouble)(dim+dim*(1+2*(2*kgrid+1)));/*not a serious formula*/
		gtk_text_buffer_insert_at_cursor(vasp_gui.simple_message_buff,
			g_strdup_printf("SET AUTO gamma-centered grid w/ AUTO = %i\n",(gint)vasp_gui.calc.kpoints_kx),-1);
	}
	if(calcul<3){
		vasp_gui.calc.prec=VP_ACCURATE;
		_OUT("SET PREC=ACCURATE\n");
		vasp_gui.calc.algo=VA_NORM;
		_OUT("SET ALGO=NORMAL\n");
		vasp_gui.calc.ldiag=TRUE;
		_OUT("SET LDIAG=.TRUE.\n");
	}
	vasp_gui.calc.use_prec=TRUE;
	/*set LREAL*/
	if(vasp_gui.calc.atoms_total>16) {
		vasp_gui.calc.lreal=VLR_AUTO;
		_OUT("SET LREAL=Auto\n");
	} else {
		vasp_gui.calc.lreal=VLR_FALSE;
		_OUT("SET LREAL=.FALSE.\n");
	}
	/*set smearing*/
	if(system>0) {
		vasp_gui.calc.ismear=0;
		_OUT("SET ISMEAR=0\n");
		vasp_gui.calc.sigma=0.05;
		_OUT("SET SIGMA=0.05\n");
	} else {
		vasp_gui.calc.ismear=1;
		_OUT("SET ISMEAR=1\n");
		vasp_gui.calc.sigma=0.2;
		_OUT("SET SIGMA=0.2\n");
	}
	/*set dipol correction*/
	_OUT("Please manually SET all dipole related settings in ELECT-I!\n");
	switch(dim){
		case 0://atom in a box
			vasp_gui.calc.idipol=VID_4;
			_OUT("SET IDIPOL=4\n");
			break;
		case 1:
			//linear system is characterized by a major axis
			vasp_gui.calc.idipol=VID_3;
			_OUT("SET IDIPOL=3\n");
			_OUT("By convention, the major direction for a 1D material is the Z axis.\n");
			_OUT("If this is not the case, RESET IDIPOL accordingly!\n");
			break;
		case 2:
			//surface are characterized by a normal axis
			vasp_gui.calc.idipol=VID_3;
			_OUT("SET IDIPOL=3\n");
			_OUT("By convention, the Z AXIS is normal to the surface of 2D material.\n");
			_OUT("If this is not the case, RESET IDIPOL accordingly!\n");
			break;
		case 3:
		default:
			vasp_gui.calc.ldipol=FALSE;
			_OUT("SET LDIPOL=.FALSE.\n");
	}
	/* spin? */
	if((gint)(vasp_gui.calc.electron_total/2)-(gdouble)(vasp_gui.calc.electron_total/2.0)!=0.0){
		vasp_gui.calc.ispin=TRUE;
		_OUT("SET ISPIN=2\n");
	}else{
		vasp_gui.calc.ispin=FALSE;
		_OUT("SET ISPIN=1\n");
	}
	/*set calculation specific*/
	switch(calcul){
	case 0://Energy (single point)
		/*use ismear=-5 if we have enough kpoints (ie more than 3)*/
		if((kgrid>1)&&(dim>0)) {
			vasp_gui.calc.ismear=-5;
			_OUT("SET ISMEAR=-5\n");
		}
		break;
	case 1://DOS/BANDS (single point)
		if(dim==0) _OUT("PB: COARSE KGRID BUT BAND/DOS REQUIRED!\n");
		if(vasp_gui.calc.have_paw) {
			vasp_gui.calc.lorbit=12;
			_OUT("SET LORBIT=12\n");
		} else {
			vasp_gui.calc.lorbit=2;
			_OUT("SET LORBIT=2\n");
		}
		vasp_gui.calc.nedos=2001;
		_OUT("SET NEDOS=2001\n");
		vasp_gui.calc.emin=-10.;
		_OUT("SET EMIN=-10.0\n");
		vasp_gui.calc.emax=10.;
		_OUT("SET EMAX=10.0\n");
		break;
	case 2://Lattice Dynamics (opt.)
		vasp_gui.calc.addgrid=TRUE;
		_OUT("SET ADDGRID=.TRUE.\n");
		/*according to vasp manual, Sec. 9.9*/
		//PREC=Accurate already set
		vasp_gui.calc.lreal=VLR_FALSE;/*force LREAL=FALSE*/
		_OUT("SET LREAL=.FALSE.\n");
		//smearing default 
		vasp_gui.calc.ibrion=6;
		_OUT("Default to finite Difference (IBRION=6)... SET IBRION=8 for Linear Perturbation!\n");
		/*according to vasp manual, Sec. 6.22.6*/
		vasp_gui.calc.potim=0.015;
		_OUT("SET POTIM=0.015\n");
		if(vasp_gui.calc.ncore>1){
			_OUT("Lattice Dynamics does not support NCORE>1... RESETING NCORE\n");
			vasp_gui.calc.ncore=1;
			_OUT("(note that kpar can be set > 1\n");
		}
		break;
	case 3://Geometry (opt.)
		if(vasp_gui.simple_rgeom){
			vasp_gui.calc.prec=VP_NORM;
			_OUT("SET PREC=NORMAL\n");
			vasp_gui.calc.use_prec=TRUE;
			vasp_gui.calc.algo=VA_FAST;
			_OUT("SET ALGO=FAST\n");
			/*according to vasp manual, Sec. 6.2.4*/
			vasp_gui.calc.nelmin=5;
			_OUT("SET NELMIN=5\n");
			vasp_gui.calc.ediff=1E-2;
			_OUT("SET EDIFF=1E-2\n");
			vasp_gui.calc.ediffg=-0.3;
			_OUT("SET EDIFFG=-0.3\n");
			vasp_gui.calc.nsw=10;
			_OUT("SET NSW=10\n");
			vasp_gui.calc.ibrion=2;
			_OUT("SET IBRION=2\n");
		}else{
			vasp_gui.calc.prec=VP_ACCURATE;
			_OUT("SET PREC=ACCURATE\n");
			vasp_gui.calc.use_prec=TRUE;
			vasp_gui.calc.algo=VA_NORM;
			_OUT("SET ALGO=NORMAL\n");
			vasp_gui.calc.ldiag=TRUE;
			_OUT("SET LDIAG=.TRUE.\n");
			vasp_gui.calc.addgrid=TRUE;
			_OUT("SET ADDGRID=.TRUE.\n");
			/*according to manual, Sec. 6.2.5*/
			vasp_gui.calc.nelmin=8;
			_OUT("SET NELMIN=8\n");
			vasp_gui.calc.ediff=1E-5;
			_OUT("SET EDIFF=1E-5\n");
			vasp_gui.calc.ediffg=-0.01;
			_OUT("SET EDIFFG=-0.01\n");
			vasp_gui.calc.nsw=20;
			_OUT("SET NSW=20\n");
			vasp_gui.calc.maxmix=80;
			_OUT("SET MAXMIX=80\n");
			vasp_gui.calc.ibrion=1;
			_OUT("SET IBRION=1\n");
			vasp_gui.calc.nfree=10;
			_OUT("SET NFREE=10\n");
		}
		break;
	case 4://Molecular Dynamics (opt.)
		if(vasp_gui.simple_rgeom) {
			vasp_gui.calc.prec=VP_LOW;
			_OUT("SET PREC=LOW\n");
		} else {
			vasp_gui.calc.prec=VP_NORM;
			_OUT("SET PREC=NROMAL\n");
		}
		vasp_gui.calc.use_prec=TRUE;/*<- this is maybe too high for PREC=NORMAL*/
		/*according to vasp manual, Sec. 9.7*/
		vasp_gui.calc.ediff=1E-5;
		_OUT("SET EDIFF=1E-5\n");
		_OUT("Please set smearing manually!\n");
		vasp_gui.calc.ismear=-1;
		_OUT("SET ISMEAR=-1\n");
		vasp_gui.calc.sigma=0.086;
		_OUT("SET SIGMA=0.086\n");
		/*choose smearing wisely*/
		vasp_gui.calc.algo=VA_VERYFAST;
		_OUT("SET ALGO=VERYFAST");
		vasp_gui.calc.maxmix=40;
		_OUT("SET MAXMIX=40\n");
		vasp_gui.calc.isym=0;
		_OUT("SET ISYM=0\n");
		vasp_gui.calc.nelmin=4;
		_OUT("SET NELMIN=4\n");
		vasp_gui.calc.ibrion=0;
		_OUT("SET IBRION=0\n");
		vasp_gui.calc.nsw=100;
		_OUT("SET NSW=100\n");
		vasp_gui.calc.nwrite=0;/*TODO: for (adv.)*/
		_OUT("SET NWRITE=0\n");
		vasp_gui.calc.lcharg=FALSE;
		_OUT("SET LCHARG=.FALSE.\n");
		vasp_gui.calc.lwave=FALSE;
		_OUT("SET LWAVE=.FALSE.\n");
		_OUT("Please set TEBEG and TEEND manually!\n");
		vasp_gui.calc.tebeg=1000;
		_OUT("SET TEBEG=1000\n");
		vasp_gui.calc.teend=1000;
		_OUT("SET TEEND=1000\n");
		vasp_gui.calc.smass=3;
		_OUT("SET SMASS=3\n");
		vasp_gui.calc.nblock=50;
		_OUT("SET NBLOCK=50\n");
		vasp_gui.calc.potim=1.5;
		_OUT("SET POTIM=1.5\n");
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
		vasp_gui.calc.prec=VP_SINGLE;break;
	case 2://Accurate
		vasp_gui.calc.prec=VP_ACCURATE;break;
	case 3://High
		vasp_gui.calc.prec=VP_HIGH;break;
	case 4://Medium
		vasp_gui.calc.prec=VP_MED;break;
	case 5://Low
		vasp_gui.calc.prec=VP_LOW;break;
	case 0://Normal
	default:
		vasp_gui.calc.prec=VP_NORM;
	}
}
/*****************************************/
/* Toggle use of automatic prec settings */
/*****************************************/
void prec_toggle(GtkWidget *w, GtkWidget *box){
	if(!(vasp_gui.calc.use_prec)){
		/*switch to manual values*/
		GUI_UNLOCK(vasp_gui.encut);
		GUI_UNLOCK(vasp_gui.enaug);
		GUI_ENTRY_TEXT(vasp_gui.encut,g_strdup_printf("%.2lf",vasp_gui.calc.encut));
		GUI_ENTRY_TEXT(vasp_gui.enaug,g_strdup_printf("%.2lf",vasp_gui.calc.enaug));
	}else{
		/*switch to automatic values*/
		GUI_LOCK(vasp_gui.encut);
		GUI_LOCK(vasp_gui.enaug);
	}


}
/******************/
/* selecting ALGO */
/******************/
void vasp_algo_selected(GtkWidget *w, struct model_pak *model){
	const gint index=gtk_combo_box_get_active(GTK_COMBO_BOX(w));
	GUI_LOCK(vasp_gui.ialgo);
	switch(index){
	case 1://USE_IALGO
		vasp_gui.calc.algo=VA_IALGO;
		GUI_UNLOCK(vasp_gui.ialgo);
		break;
	case 2://VERYFAST
		vasp_gui.calc.algo=VA_VERYFAST;break;
	case 3://FAST
		vasp_gui.calc.algo=VA_FAST;break;
	case 4://CONJUGATE
		vasp_gui.calc.algo=VA_CONJ;break;
	case 5://ALL
		vasp_gui.calc.algo=VA_ALL;break;
	case 6://DAMPED
		vasp_gui.calc.algo=VA_DAMPED;break;
	case 7://SUBROT
		vasp_gui.calc.algo=VA_SUBROT;break;
	case 8://EIGENVAL
		vasp_gui.calc.algo=VA_EIGEN;break;
	case 9://NONE
		vasp_gui.calc.algo=VA_NONE;break;
	case 10://NOTHING
		vasp_gui.calc.algo=VA_NOTHING;break;
	case 11://EXACT
		vasp_gui.calc.algo=VA_EXACT;break;
	case 12://DIAG
		vasp_gui.calc.algo=VA_DIAG;break;
	case 0://NORMAL
	default:
		vasp_gui.calc.algo=VA_NORM;
	}
}
/*******************/
/* selecting IALGO */
/*******************/
void vasp_ialgo_selected(GtkWidget *w, struct model_pak *model){
	const gint index=gtk_combo_box_get_active(GTK_COMBO_BOX(w));
	switch(index){
	case 0://2:FIXED ORB/1E
		vasp_gui.calc.ialgo=VIA_OE_FIXED;break;
	case 1://3:FIXED ORB
		vasp_gui.calc.ialgo=VIA_O_FIXED;break;
	case 2://4:SUBROT ONLY
		vasp_gui.calc.ialgo=VIA_SUBROT;break;
	case 3://5:STEEPEST
		vasp_gui.calc.ialgo=VIA_STEEP;break;
	case 4://6:CONJUGUATE GRADIENTS
		vasp_gui.calc.ialgo=VIA_CG;break;
	case 5://7:PRECOND. STEEPEST
		vasp_gui.calc.ialgo=VIA_PSTEEP;break;
	case 6://8:PRECOND. CG
		vasp_gui.calc.ialgo=VIA_PCG;break;
	case 8://44:RMM STEEPEST
		vasp_gui.calc.ialgo=VIA_ESTEEP;break;
	case 9://46:RMM PRECOND.
		vasp_gui.calc.ialgo=VIA_RMMP;break;
	case 10://48:PRECOND. RMM
		vasp_gui.calc.ialgo=VIA_PRMM;break;
	case 11://53:VAR DAMPED MD
		vasp_gui.calc.ialgo=VIA_VAR_DAMP;break;
	case 12://54:VAR DAMPED QUENCH MD
		vasp_gui.calc.ialgo=VIA_VAR_QUENCH;break;
	case 13://58:VAR PRECOND. CG
		vasp_gui.calc.ialgo=VIA_VAR_PCG;break;
	case 14://90:EXACT
		vasp_gui.calc.ialgo=VIA_EXACT;break;
	case 7://38:KOSUGI
	default:
		vasp_gui.calc.ialgo=VIA_KOSUGI;
}
}
/***************/
/* elec toggle */
/***************/
void elec_toggle(GtkWidget *w, GtkWidget *box){
	if(!(vasp_gui.calc.auto_elec)){
		/*switch to manual values*/
		GUI_UNLOCK(vasp_gui.nsim);
		GUI_UNLOCK(vasp_gui.vtime);
		GUI_UNLOCK(vasp_gui.nbands);
		GUI_UNLOCK(vasp_gui.nelect);
		GUI_UNLOCK(vasp_gui.iwavpr);
		GUI_ENTRY_TEXT(vasp_gui.nsim,g_strdup_printf("%i",vasp_gui.calc.nsim));
		GUI_ENTRY_TEXT(vasp_gui.vtime,g_strdup_printf("%.4lf",vasp_gui.calc.vtime));
		GUI_ENTRY_TEXT(vasp_gui.nbands,g_strdup_printf("%i",vasp_gui.calc.nbands));
		GUI_ENTRY_TEXT(vasp_gui.nelect,g_strdup_printf("%.4lf",vasp_gui.calc.nelect));
		GUI_ENTRY_TEXT(vasp_gui.iwavpr,g_strdup_printf("%i",vasp_gui.calc.iwavpr));
	}else{
		/*switch to automatic values*/
		GUI_LOCK(vasp_gui.nsim);
		GUI_LOCK(vasp_gui.vtime);
		GUI_LOCK(vasp_gui.nbands);
		GUI_LOCK(vasp_gui.nelect);
		GUI_LOCK(vasp_gui.iwavpr);
	}
}
/*******************/
/* selecting LREAL */
/*******************/
void vasp_lreal_selected(GtkWidget *w, struct model_pak *model){
	const gint index=gtk_combo_box_get_active(GTK_COMBO_BOX(w));
	switch(index){
	case 0://Auto
		vasp_gui.calc.lreal=VLR_AUTO;break;
	case 1://On
		vasp_gui.calc.lreal=VLR_ON;break;
	case 2://TRUE
		vasp_gui.calc.lreal=VLR_TRUE;break;
	case 3://FALSE
	default:
		vasp_gui.calc.lreal=VLR_FALSE;
	}
}
/***************/
/* grid toggle */
/***************/
void grid_toggle(GtkWidget *w, GtkWidget *box){
        if(!(vasp_gui.calc.auto_grid)){
                /*switch to manual values*/
		GUI_UNLOCK(vasp_gui.ngx);
		GUI_UNLOCK(vasp_gui.ngy);
		GUI_UNLOCK(vasp_gui.ngz);
		GUI_UNLOCK(vasp_gui.ngxf);
		GUI_UNLOCK(vasp_gui.ngyf);
		GUI_UNLOCK(vasp_gui.ngzf);
		GUI_ENTRY_TEXT(vasp_gui.ngx,g_strdup_printf("%i",vasp_gui.calc.ngx));
		GUI_ENTRY_TEXT(vasp_gui.ngy,g_strdup_printf("%i",vasp_gui.calc.ngy));
		GUI_ENTRY_TEXT(vasp_gui.ngz,g_strdup_printf("%i",vasp_gui.calc.ngz));
		GUI_ENTRY_TEXT(vasp_gui.ngxf,g_strdup_printf("%i",vasp_gui.calc.ngxf));
		GUI_ENTRY_TEXT(vasp_gui.ngyf,g_strdup_printf("%i",vasp_gui.calc.ngyf));
		GUI_ENTRY_TEXT(vasp_gui.ngzf,g_strdup_printf("%i",vasp_gui.calc.ngzf));
        }else{
                /*switch to automatic values*/
		GUI_LOCK(vasp_gui.ngx);
		GUI_LOCK(vasp_gui.ngy);
		GUI_LOCK(vasp_gui.ngz);
		GUI_LOCK(vasp_gui.ngxf);
		GUI_LOCK(vasp_gui.ngyf);
		GUI_LOCK(vasp_gui.ngzf);
        }
}
/***************************/
/* ismear value is changed */
/***************************/
void ismear_changed(GtkWidget *w, struct model_pak *model){
        const gchar *label = gtk_entry_get_text(GTK_ENTRY(w));
        sscanf(label,"%i",&(vasp_gui.calc.ismear));
	if(vasp_gui.calc.ismear==-5) {
		GUI_TOGGLE_ON(vasp_gui.have_tetra);
		GUI_UNLOCK(vasp_gui.have_tetra);
	}
	else {
		GUI_TOGGLE_OFF(vasp_gui.have_tetra);
		GUI_LOCK(vasp_gui.have_tetra);
	}
}
/*****************************/
/* kspacing value is changed */
/*****************************/
void kspacing_changed(GtkWidget *w, struct model_pak *model){
	const gchar *label = gtk_entry_get_text(GTK_ENTRY(w));
        sscanf(label,"%lf",&(vasp_gui.calc.kspacing));
}
/************************/
/* toggle LNONCOLLINEAR */
/************************/
void lnoncoll_toggle(GtkWidget *w, GtkWidget *box){
        /* unselecting LNONCOLLINEAR automatically unselect LSORBIT */
        if(!vasp_gui.calc.non_collinear) gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(vasp_gui.lsorbit),FALSE);
}
/******************/
/* toggle LSORBIT */
/******************/
void lsorbit_toggle(GtkWidget *w, GtkWidget *box){
	/* LSORBIT automatically select LNONCOLLINEAR */
	if(vasp_gui.calc.lsorbit) gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(vasp_gui.lnoncoll),TRUE);
}
/*****************/
/* selecting GGA */
/*****************/
void vasp_gga_selected(GtkWidget *w, struct model_pak *model){
	const gint index=gtk_combo_box_get_active(GTK_COMBO_BOX(w));
	switch(index){
	case 0://91:PW91
		vasp_gui.calc.gga=VG_91;break;
	case 2://RP:rPBE
		vasp_gui.calc.gga=VG_RP;break;
	case 3://AM:AM05
		vasp_gui.calc.gga=VG_AM;break;
	case 4://PS:PBEsol
		vasp_gui.calc.gga=VG_PS;break;
	case 1://PE:PBE
	default:
		vasp_gui.calc.gga=VG_PE;
	}
}
/*********************/
/* selecting METAGGA */
/*********************/
void vasp_metagga_selected(GtkWidget *w, struct model_pak *model){
	const gint index=gtk_combo_box_get_active(GTK_COMBO_BOX(w));
	GUI_LOCK(vasp_gui.cmbj);
	GUI_LOCK(vasp_gui.cmbja);
	GUI_LOCK(vasp_gui.cmbjb);
	switch(index){
	case 0://TPSS
		vasp_gui.calc.mgga=VMG_TPSS;break;
	case 1://rTPSS
		vasp_gui.calc.mgga=VMG_RTPSS;break;
	case 2://M06L
		vasp_gui.calc.mgga=VMG_M06L;break;
	case 3://MBJ
	default:
		vasp_gui.calc.mgga=VMG_MBJ;
		GUI_UNLOCK(vasp_gui.cmbj);
		GUI_UNLOCK(vasp_gui.cmbja);
		GUI_UNLOCK(vasp_gui.cmbjb);
	}
}
/*******************/
/* toggle LMETAGGA */
/*******************/
void metagga_toggle(GtkWidget *w, GtkWidget *box){
        if(!vasp_gui.calc.lmetagga){
                /*switch to manual values*/
		GUI_LOCK(vasp_gui.metagga);
		GUI_LOCK(vasp_gui.cmbj);
		GUI_LOCK(vasp_gui.cmbja);
		GUI_LOCK(vasp_gui.cmbjb);
        }else{
                /*switch to automatic values*/
		GUI_UNLOCK(vasp_gui.metagga);
		GUI_COMBOBOX_SET(vasp_gui.metagga,3);//MBJ (actually there is no default)
		GUI_UNLOCK(vasp_gui.cmbj);
		GUI_UNLOCK(vasp_gui.cmbja);
		GUI_UNLOCK(vasp_gui.cmbjb);
        }
}
/*******************/
/* toggle L(S)DA+U */
/*******************/
void ldau_toggle(GtkWidget *w, GtkWidget *box){
        if(!vasp_gui.calc.ldau){
                /*switch to manual values*/
		GUI_LOCK(vasp_gui.ldau);
		GUI_LOCK(vasp_gui.ldau_print);
		GUI_LOCK(vasp_gui.ldaul);
		GUI_LOCK(vasp_gui.ldauu);
		GUI_LOCK(vasp_gui.ldauj);
        }else{
                /*switch to automatic values*/
		GUI_UNLOCK(vasp_gui.ldau);
		GUI_UNLOCK(vasp_gui.ldau_print);
		GUI_UNLOCK(vasp_gui.ldaul);
		GUI_UNLOCK(vasp_gui.ldauu);
		GUI_UNLOCK(vasp_gui.ldauj);
		GUI_COMBOBOX_SET(vasp_gui.ldau,1);//2:LSDA+U Dudarev
		GUI_COMBOBOX_SET(vasp_gui.ldau_print,0);//0:Silent
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
		vasp_gui.calc.ldau_type=VU_1;break;
	case 2://2:LSDA+U Dudarev
		vasp_gui.calc.ldau_type=VU_4;break;
	case 1://4:LDA+U Liechtenstein
	default:
		vasp_gui.calc.ldau_type=VU_2;
	}
}
/***************************************/
/* selecting L(S)DA+U output verbosity */
/***************************************/
void vasp_ldau_print_selected(GtkWidget *w, struct model_pak *model){
	const gint index=gtk_combo_box_get_active(GTK_COMBO_BOX(w));
	switch(index){
	case 1://1:Occupancy matrix
		vasp_gui.calc.ldau_output=VUO_1;break;
	case 2://2:Full output
		vasp_gui.calc.ldau_output=VUO_2;break;
	case 0://0:Silent
	default:
		vasp_gui.calc.ldau_output=VUO_0;
	}
}
/************************************/
/* toggle manual selection of mixer */
/************************************/
void mixer_toggle(GtkWidget *w, GtkWidget *box){
        if((vasp_gui.calc.auto_mixer)){
		/*switch to manual values*/
		GUI_LOCK(vasp_gui.mixer);
		GUI_LOCK(vasp_gui.amix);
		GUI_LOCK(vasp_gui.bmix);
		GUI_LOCK(vasp_gui.amin);
		GUI_LOCK(vasp_gui.amix_mag);
		GUI_LOCK(vasp_gui.bmix_mag);
		GUI_LOCK(vasp_gui.maxmix);
		GUI_LOCK(vasp_gui.wc);
		GUI_LOCK(vasp_gui.inimix);
		GUI_LOCK(vasp_gui.mixpre);
        }else{
                /*switch to automatic values*/
		GUI_UNLOCK(vasp_gui.mixer);
		GUI_UNLOCK(vasp_gui.amix);
		GUI_UNLOCK(vasp_gui.bmix);
		GUI_UNLOCK(vasp_gui.amin);
		GUI_UNLOCK(vasp_gui.amix_mag);
		GUI_UNLOCK(vasp_gui.bmix_mag);
		GUI_UNLOCK(vasp_gui.maxmix);
		GUI_COMBOBOX_SET(vasp_gui.mixer,3);//4:BROYDEN
		GUI_ENTRY_TEXT(vasp_gui.amix,g_strdup_printf("%.4lf",vasp_gui.calc.amix));
		GUI_ENTRY_TEXT(vasp_gui.bmix,g_strdup_printf("%.4lf",vasp_gui.calc.bmix));
		GUI_ENTRY_TEXT(vasp_gui.amin,g_strdup_printf("%.4lf",vasp_gui.calc.amin));
		GUI_ENTRY_TEXT(vasp_gui.amix_mag,g_strdup_printf("%.4lf",vasp_gui.calc.amix_mag));
		GUI_ENTRY_TEXT(vasp_gui.bmix_mag,g_strdup_printf("%.4lf",vasp_gui.calc.bmix_mag));
		GUI_ENTRY_TEXT(vasp_gui.maxmix,g_strdup_printf("%i",vasp_gui.calc.maxmix));
		if(vasp_gui.calc.imix==VM_4){
			GUI_UNLOCK(vasp_gui.wc);
			GUI_UNLOCK(vasp_gui.inimix);
			GUI_UNLOCK(vasp_gui.mixpre);
		}
        }
}
/*******************/
/* selecting mixer */
/*******************/
void vasp_mixer_selected(GtkWidget *w, struct model_pak *model){
	const gint index=gtk_combo_box_get_active(GTK_COMBO_BOX(w));
	GUI_LOCK(vasp_gui.wc);
	GUI_LOCK(vasp_gui.inimix);
	GUI_LOCK(vasp_gui.mixpre);
	switch(index){
	case 0://0:NO MIXING
		vasp_gui.calc.imix=VM_0;break;
	case 1://1:KERKER
		vasp_gui.calc.imix=VM_1;break;
	case 2://2:TCHEBYCHEV
		vasp_gui.calc.imix=VM_2;break;
	case 3://4:BROYDEN
	default:
		vasp_gui.calc.imix=VM_4;
		GUI_UNLOCK(vasp_gui.wc);
		GUI_UNLOCK(vasp_gui.inimix);
		GUI_UNLOCK(vasp_gui.mixpre);
	}
}
/***************************/
/* selecting initial mixer */
/***************************/
void vasp_inimix_selected(GtkWidget *w, struct model_pak *model){
	const gint index=gtk_combo_box_get_active(GTK_COMBO_BOX(w));
	switch(index){
	case 0://0:LINEAR
		vasp_gui.calc.inimix=VIM_0;break;
	case 1://1:KERKER
	default:
		vasp_gui.calc.inimix=VIM_1;
	}
}
/************************/
/* selecting pre-mixing */
/************************/
void vasp_mixpre_selected(GtkWidget *w, struct model_pak *model){
	const gint index=gtk_combo_box_get_active(GTK_COMBO_BOX(w));
	switch(index){
	case 0://0:NONE
		vasp_gui.calc.mixpre=VMP_0;break;
	case 2:
		vasp_gui.calc.mixpre=VMP_2;break;
	case 1:
	default:
		vasp_gui.calc.mixpre=VMP_1;
	}
}
/********************/
/* selecting idipol */
/********************/
void vasp_idipol_selected(GtkWidget *w, struct model_pak *model){
	const gint index=gtk_combo_box_get_active(GTK_COMBO_BOX(w));
	GUI_UNLOCK(vasp_gui.ldipol);
	GUI_UNLOCK(vasp_gui.lmono);
	GUI_UNLOCK(vasp_gui.dipol);
	GUI_UNLOCK(vasp_gui.epsilon);
	GUI_UNLOCK(vasp_gui.efield);
	switch(index){
	case 1://1:u axis
		vasp_gui.calc.idipol=VID_1;break;
	case 2://2:v axis
		vasp_gui.calc.idipol=VID_2;break;
	case 3://3:w axis
		vasp_gui.calc.idipol=VID_3;break;
	case 4://4:all axis
		vasp_gui.calc.idipol=VID_4;
		GUI_LOCK(vasp_gui.efield);
		break;
	case 0://0:no calcul
	default:
		vasp_gui.calc.idipol=VID_0;
		GUI_LOCK(vasp_gui.ldipol);
		GUI_LOCK(vasp_gui.lmono);
		GUI_LOCK(vasp_gui.dipol);
		GUI_LOCK(vasp_gui.epsilon);
		GUI_LOCK(vasp_gui.efield);
	}
}
/********************/
/* selecting lorbit */
/********************/
void vasp_lorbit_selected(GtkWidget *w, struct model_pak *model){
	const gint index=gtk_combo_box_get_active(GTK_COMBO_BOX(w));
	switch(index){
	case 1://1:lm-PROCAR
		vasp_gui.calc.lorbit=1;break;
	case 2://2:phase+lm-PROCAR
		vasp_gui.calc.lorbit=2;break;
	case 3://5:PROOUT
		vasp_gui.calc.lorbit=5;break;
	case 4://10:PROCAR
		vasp_gui.calc.lorbit=10;break;
	case 5://11:lm-PROCAR
		vasp_gui.calc.lorbit=11;break;
	case 6://12:phase+lm-PROCAR
		vasp_gui.calc.lorbit=12;break;
	case 0://0:PROCAR
	default:
		vasp_gui.calc.lorbit=0;
	}
}
/********************/
/* selecting ibrion */
/********************/
void vasp_ibrion_selected(GtkWidget *w, struct model_pak *model){
	const gint index=gtk_combo_box_get_active(GTK_COMBO_BOX(w));
GUI_LOCK(vasp_gui.tebeg);
GUI_LOCK(vasp_gui.teend);
GUI_LOCK(vasp_gui.smass);
GUI_LOCK(vasp_gui.nblock);
GUI_LOCK(vasp_gui.kblock);
GUI_LOCK(vasp_gui.npaco);
GUI_LOCK(vasp_gui.apaco);
	switch(index){
	case 1://0:Molecular Dynamics
		vasp_gui.calc.ibrion=0;
		/*molecular dynamics part*/
GUI_UNLOCK(vasp_gui.tebeg);
GUI_UNLOCK(vasp_gui.teend);
GUI_UNLOCK(vasp_gui.smass);
GUI_UNLOCK(vasp_gui.nblock);
GUI_UNLOCK(vasp_gui.kblock);
GUI_UNLOCK(vasp_gui.npaco);
GUI_UNLOCK(vasp_gui.apaco);
		break;
	case 2://1:Quasi-Newton
		vasp_gui.calc.ibrion=1;break;
	case 3://2:Conjugate Gradients
		vasp_gui.calc.ibrion=2;break;
	case 4://3:Damped
		vasp_gui.calc.ibrion=3;break;
	case 5://5:Finite Differences
		vasp_gui.calc.ibrion=5;break;
	case 6://6:Finite Diff. with SYM
		vasp_gui.calc.ibrion=6;break;
	case 7://7:Perturbation Theory
		vasp_gui.calc.ibrion=7;break;
	case 8://8:Pert. Theory with SYM
		vasp_gui.calc.ibrion=8;break;
	case 9://44:Transition State
		vasp_gui.calc.ibrion=44;break;
	case 0://-1:Nothing
	default:
		vasp_gui.calc.ibrion=-1;
	}
}
/******************/
/* selecting ISIF */
/******************/
void vasp_isif_selected(GtkWidget *w, struct model_pak *model){
	const gint index=gtk_combo_box_get_active(GTK_COMBO_BOX(w));
	vasp_gui.calc.isif=index;
	switch(index){
	case 3://3:F_S_I_S_V
		GUI_TOGGLE_ON(vasp_gui.relax_ions);
		GUI_TOGGLE_ON(vasp_gui.relax_shape);
		GUI_TOGGLE_ON(vasp_gui.relax_volume);
		break;
	case 4://4:F_S_I_S_0
		GUI_TOGGLE_ON(vasp_gui.relax_ions);
		GUI_TOGGLE_ON(vasp_gui.relax_shape);
		GUI_TOGGLE_OFF(vasp_gui.relax_volume);
		break;
	case 5://5:F_S_0_S_0
		GUI_TOGGLE_OFF(vasp_gui.relax_ions);
		GUI_TOGGLE_ON(vasp_gui.relax_shape);
		GUI_TOGGLE_OFF(vasp_gui.relax_volume);
		break;	
	case 6://6:F_S_0_S_V
		GUI_TOGGLE_OFF(vasp_gui.relax_ions);
		GUI_TOGGLE_ON(vasp_gui.relax_shape);
		GUI_TOGGLE_ON(vasp_gui.relax_volume);
		break;
	case 7://7:F_S_0_0_V
		GUI_TOGGLE_OFF(vasp_gui.relax_ions);
		GUI_TOGGLE_OFF(vasp_gui.relax_shape);
		GUI_TOGGLE_ON(vasp_gui.relax_volume);
		break;
	case 0://0:F_0_I_0_0
	case 1://1:F_P_I_0_0
	case 2://2:F_S_I_0_0
	default:
		GUI_TOGGLE_ON(vasp_gui.relax_ions);
		GUI_TOGGLE_OFF(vasp_gui.relax_shape);
		GUI_TOGGLE_OFF(vasp_gui.relax_volume);
	}
}
/****************/
/* isif toggles */
/****************/
void isif_toggle(GtkWidget *w, GtkWidget *box){
/*recalculate isif function of the user switches*/
if((vasp_gui.rions)&&(!vasp_gui.rshape)&&(!vasp_gui.rvolume)) {if(vasp_gui.calc.isif>2) vasp_gui.calc.isif=2;}
else if((vasp_gui.rions)&&(vasp_gui.rshape)&&(vasp_gui.rvolume)) vasp_gui.calc.isif=3;
else if ((vasp_gui.rions)&&(vasp_gui.rshape)&&(!vasp_gui.rvolume)) vasp_gui.calc.isif=4;
else if ((!vasp_gui.rions)&&(vasp_gui.rshape)&&(!vasp_gui.rvolume)) vasp_gui.calc.isif=5;
else if ((!vasp_gui.rions)&&(vasp_gui.rshape)&&(vasp_gui.rvolume)) vasp_gui.calc.isif=6;
else if ((!vasp_gui.rions)&&(!vasp_gui.rshape)&&(vasp_gui.rvolume)) vasp_gui.calc.isif=7;
/*deal with forbiden combination*/
if((vasp_gui.rions)&&(!vasp_gui.rshape)&&(vasp_gui.rvolume)) {
	GUI_TOGGLE_OFF(vasp_gui.relax_ions);
	vasp_gui.calc.isif=7;
}
if((!vasp_gui.rions)&&(!vasp_gui.rshape)&&(!vasp_gui.rvolume)) {
	GUI_TOGGLE_ON(vasp_gui.relax_ions);
	vasp_gui.calc.isif=2;
}
GUI_COMBOBOX_SET(vasp_gui.isif,vasp_gui.calc.isif);
}
/*****************************/
/* selective dynamics toggle */
/*****************************/
void toggle_poscar_sd(){
	if((vasp_gui.calc.poscar_free==VPF_MAN)&&(vasp_gui.calc.poscar_sd)){
		/*we need some tags*/
		GUI_UNLOCK(vasp_gui.poscar_tx);
		GUI_UNLOCK(vasp_gui.poscar_ty);
		GUI_UNLOCK(vasp_gui.poscar_tz);
	}else{
		/*tags are no longer important*/
		GUI_LOCK(vasp_gui.poscar_tx);
		GUI_LOCK(vasp_gui.poscar_ty);
		GUI_LOCK(vasp_gui.poscar_tz);
	}
	if(vasp_gui.calc.poscar_sd) GUI_UNLOCK(vasp_gui.poscar_free);
	else GUI_LOCK(vasp_gui.poscar_free);
}
/*************************/
/* selecting poscar_free */
/*************************/
void vasp_poscar_free_selected(GtkWidget *w, struct model_pak *model){
	gint index;
	gint ix;
	gchar *text;
	gdouble x,y,z;
	gchar symbol[3];
	gint idx=0;
	gchar tag;
	GUI_COMBOBOX_GET(w,index);
	GUI_COMBOBOX_GET(vasp_gui.poscar_atoms,ix);/*PUSH*/
	if(index==2){
		GUI_UNLOCK(vasp_gui.poscar_tx);
		GUI_UNLOCK(vasp_gui.poscar_ty);
		GUI_UNLOCK(vasp_gui.poscar_tz);
	}else{
		GUI_LOCK(vasp_gui.poscar_tx);
		GUI_LOCK(vasp_gui.poscar_ty);
		GUI_LOCK(vasp_gui.poscar_tz);
	}
	switch(index){
	case 0://All atom FIXED
		vasp_gui.calc.poscar_free=VPF_FIXED;
		tag='F';
		break;
	case 2://Manual selection
		vasp_gui.calc.poscar_free=VPF_MAN;
		/*no need to modify atom list*/
		return;
	case 1://All atom FREE
	default:
		vasp_gui.calc.poscar_free=VPF_FREE;
		tag='T';
	}
	GUI_COMBOBOX_SET(vasp_gui.poscar_atoms,idx);
	GUI_COMBOBOX_GET_TEXT(vasp_gui.poscar_atoms,text);
	while(g_ascii_strcasecmp(text,"ADD atom")!=0){
		/*modify each line tags*/
		sscanf(text,"%lf %lf %lf %*c %*c %*c ! atom: %*i (%[^)])",&x,&y,&z,&(symbol[0]));
		g_free(text);/*combobox returned text should be freed*/
		text=g_strdup_printf("%.8lf %.8lf %.8lf  %c   %c   %c ! atom: %i (%s)",x,y,z,tag,tag,tag,idx,symbol);
		GUI_COMBOBOX_ADD_TEXT(vasp_gui.poscar_atoms,idx,text);
		GUI_COMBOBOX_DEL(vasp_gui.poscar_atoms,idx+1);
		idx++;
		GUI_COMBOBOX_SET(vasp_gui.poscar_atoms,idx);
		g_free(text);/*so is g_strdup text*/
		GUI_COMBOBOX_GET_TEXT(vasp_gui.poscar_atoms,text);
	}
	/*return to selected atom*/
	g_free(text);/*free "ADD atom" string*/
	GUI_COMBOBOX_SET(vasp_gui.poscar_atoms,ix);/*PULL*/
}
/************************/
/* poscar_direct toggle */
/************************/
void toggle_poscar_direct(){
	gint index;//=gtk_combo_box_get_active(GTK_COMBO_BOX(vasp_gui.poscar_atoms));/*PUSH*/
	gint idx=0;
	gchar *text;
	gchar *tamp;
	gdouble x,y,z;
	gdouble ux,uy,uz;
	gdouble vx,vy,vz;
	gdouble wx,wy,wz;
	/* TODO: recalculate all positions on a single click */
	GUI_COMBOBOX_GET(vasp_gui.poscar_atoms,index);/*PUSH*/
	/*mini sync*/
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
	ux=vasp_gui.calc.poscar_ux*vasp_gui.calc.poscar_a0;
	uy=vasp_gui.calc.poscar_uy*vasp_gui.calc.poscar_a0;
	uz=vasp_gui.calc.poscar_uz*vasp_gui.calc.poscar_a0;
	vx=vasp_gui.calc.poscar_vx*vasp_gui.calc.poscar_a0;
	vy=vasp_gui.calc.poscar_vy*vasp_gui.calc.poscar_a0;
	vz=vasp_gui.calc.poscar_vz*vasp_gui.calc.poscar_a0;
	wx=vasp_gui.calc.poscar_wx*vasp_gui.calc.poscar_a0;
	wy=vasp_gui.calc.poscar_wy*vasp_gui.calc.poscar_a0;
	wz=vasp_gui.calc.poscar_wz*vasp_gui.calc.poscar_a0;
	/**/
	GUI_COMBOBOX_SET(vasp_gui.poscar_atoms,idx);
	GUI_COMBOBOX_GET_TEXT(vasp_gui.poscar_atoms,text);
        while(g_ascii_strcasecmp(text,"ADD atom")!=0){
                /*get each line*/
		tamp=g_strdup(text);/*to get enough space*/
		sscanf(text,"%lf %lf %lf %[^\n]",&x,&y,&z,tamp);
		g_free(text);
		/*convert*/
			if(vasp_gui.calc.poscar_direct){
				/*was cartesian, now changed to direct*/
				x=x/(ux+vx+wx);
				y=y/(uy+vy+wy);
				z=z/(uz+vz+wz);
			} else {
				/*was direct, now changed to cartesian*/
				x=x*(ux+vx+wx);
				y=y*(uy+vy+wy);
				z=z*(uz+vz+wz);
			}
		/*rewrite data*/
		text=g_strdup_printf("%.8lf %.8lf %.8lf  %s",x,y,z,tamp);
		GUI_COMBOBOX_ADD_TEXT(vasp_gui.poscar_atoms,idx,text);
		GUI_COMBOBOX_DEL(vasp_gui.poscar_atoms,idx+1);
		g_free(tamp);
		g_free(text);
		idx++;
		/*get new line*/
		GUI_COMBOBOX_SET(vasp_gui.poscar_atoms,idx);
		GUI_COMBOBOX_GET_TEXT(vasp_gui.poscar_atoms,text);
	}
	g_free(text);/*free "ADD atom" string*/
	GUI_COMBOBOX_SET(vasp_gui.poscar_atoms,index);/*PULL*/
}
/*************************/
/* selecting poscar atom */
/*************************/
void vasp_poscar_atoms_selected(GtkWidget *w, struct model_pak *model){
	gint index;
	gchar *text;
	gdouble x,y,z;
	gchar tx,ty,tz;
	gchar symbol[3];
	GUI_COMBOBOX_GET(w,index);
	text=g_strdup_printf("%i",index);
	GUI_ENTRY_TEXT(vasp_gui.poscar_index,text);
	g_free(text);
	GUI_COMBOBOX_GET_TEXT(w,text);
	if (g_ascii_strcasecmp(text,"ADD atom") == 0) {
		/*allow symbol*/
		GUI_UNLOCK(vasp_gui.poscar_symbol);
		}else{
		/*populate atom properties*/
		sscanf(text,"%lf %lf %lf %c %c %c ! atom: %*i (%[^)])",&x,&y,&z,&tx,&ty,&tz,&(symbol[0]));
		g_free(text);text=g_strdup_printf("%s",symbol);
		GUI_ENTRY_TEXT(vasp_gui.poscar_symbol,text);
		g_free(text);text=g_strdup_printf("%.6lf",x);
		GUI_ENTRY_TEXT(vasp_gui.poscar_x,text);
		g_free(text);text=g_strdup_printf("%.6lf",y);
		GUI_ENTRY_TEXT(vasp_gui.poscar_y,text);
		g_free(text);text=g_strdup_printf("%.6lf",z);
		GUI_ENTRY_TEXT(vasp_gui.poscar_z,text);
		/* register all that */
		vasp_gui.calc.poscar_index=index;
		vasp_gui.calc.poscar_x=x;
		vasp_gui.calc.poscar_y=y;
		vasp_gui.calc.poscar_z=z;
		vasp_gui.calc.poscar_tx=(tx=='T');
		vasp_gui.calc.poscar_ty=(ty=='T');
		vasp_gui.calc.poscar_tz=(tz=='T');
		/* disallow symbol */
		GUI_LOCK(vasp_gui.poscar_symbol);
	}
	g_free(text);
	/* set tags sensitivity */
	if(vasp_gui.calc.poscar_free==VPF_MAN){
		GUI_UNLOCK(vasp_gui.poscar_tx);
		GUI_UNLOCK(vasp_gui.poscar_ty);
		GUI_UNLOCK(vasp_gui.poscar_tz);
	}else{
		GUI_LOCK(vasp_gui.poscar_tx);
		GUI_LOCK(vasp_gui.poscar_ty);
		GUI_LOCK(vasp_gui.poscar_tz);
	}
	/* (re)set tags value! */
	if(vasp_gui.calc.poscar_tx) GUI_TOGGLE_ON(vasp_gui.poscar_tx);
	else GUI_TOGGLE_OFF(vasp_gui.poscar_tx);
	if(vasp_gui.calc.poscar_ty) GUI_TOGGLE_ON(vasp_gui.poscar_ty);
	else GUI_TOGGLE_OFF(vasp_gui.poscar_ty);
	if(vasp_gui.calc.poscar_tz) GUI_TOGGLE_ON(vasp_gui.poscar_tz);
	else GUI_TOGGLE_OFF(vasp_gui.poscar_tz);
}
/**************************/
/* modify/add poscar atom */
/**************************/
void vasp_atom_modified(){
	gint index;
	gchar *text;
	gchar *tamp;
	gchar tx,ty,tz;
	/*mini-sync*/
	GUI_COMBOBOX_GET(vasp_gui.poscar_atoms,index);/*PUSH*/
	VASP_REG_VAL(poscar_index,"%i");
	VASP_REG_VAL(poscar_x,"%lf");
	VASP_REG_VAL(poscar_y,"%lf");
	VASP_REG_VAL(poscar_z,"%lf");
	GUI_ENTRY_GET_TEXT(vasp_gui.poscar_symbol,text);
	sprintf(vasp_gui.calc.poscar_symbol,"%s",text);
	g_free(text);
	tx='F';
	ty='F';
	tz='F';
	switch (vasp_gui.calc.poscar_free){
	case VPF_MAN:
		if(vasp_gui.calc.poscar_tx) tx='T';
		if(vasp_gui.calc.poscar_ty) ty='T';
		if(vasp_gui.calc.poscar_tz) tz='T';
		break;
	case VPF_FREE:
		tx='T';ty='T';tz='T';
		break;
	case VPF_FIXED:
	default:
		break;
	}
	/* add/modify */
	GUI_COMBOBOX_GET_TEXT(vasp_gui.poscar_atoms,text);
	if(g_ascii_strcasecmp(text,"ADD atom") == 0){
		/*we are going to insert some data*/
		g_free(text);
		text=g_strdup_printf("%.8lf %.8lf %.8lf  %c   %c   %c ! atom: %i (%s)",
			vasp_gui.calc.poscar_x,vasp_gui.calc.poscar_y,vasp_gui.calc.poscar_z,tx,ty,tz,index,vasp_gui.calc.poscar_symbol);
		GUI_COMBOBOX_ADD_TEXT(vasp_gui.poscar_atoms,index,text);
	}else{
		/*we are goign to modify some data*/
		g_free(text);
		text=g_strdup_printf("%.8lf %.8lf %.8lf  %c   %c   %c ! atom: %i (%s)",
			vasp_gui.calc.poscar_x,vasp_gui.calc.poscar_y,vasp_gui.calc.poscar_z,tx,ty,tz,index,vasp_gui.calc.poscar_symbol);
		GUI_COMBOBOX_ADD_TEXT(vasp_gui.poscar_atoms,index,text);
		GUI_COMBOBOX_DEL(vasp_gui.poscar_atoms,index+1);
	}
	g_free(text);
	/*re-select current entry (to refresh content)*/
	GUI_COMBOBOX_SET(vasp_gui.poscar_atoms,index);/*PULL*/
	vasp_gui.poscar_dirty=TRUE;
}
/**********************/
/* delete poscar atom */
/**********************/
void vasp_atom_delete(){
	gint index;
	gchar *text;
	GUI_COMBOBOX_GET(vasp_gui.poscar_atoms,index);
	if(index==-1) return;/*nothing selected anyway*/
	GUI_COMBOBOX_GET_TEXT(vasp_gui.poscar_atoms,text);
	if(g_ascii_strcasecmp(text,"ADD atom") == 0) {
		g_free(text);
		return;/*can't delete this one*/
	}
	GUI_COMBOBOX_DEL(vasp_gui.poscar_atoms,index);
	GUI_COMBOBOX_SET(vasp_gui.poscar_atoms,index-1);
	g_free(text);
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
		vasp_gui.calc.kpoints_mode=VKP_MAN;
		gtk_widget_set_sensitive(vasp_gui.kpoints_w,TRUE);
		gtk_widget_set_sensitive(vasp_gui.kpoints_coord,TRUE);//need coord
		gtk_widget_set_sensitive(vasp_gui.kpoints_nkpts,TRUE);//need nkpts
		gtk_widget_set_sensitive(vasp_gui.kpoints_cart,TRUE);//can be cart
		if(vasp_gui.calc.ismear==-5) gtk_widget_set_sensitive(vasp_gui.kpoints_tetra,TRUE);//only VKP_MAN require tetra
		if(vasp_gui.calc.ismear==-5) gtk_widget_set_sensitive(vasp_gui.have_tetra,TRUE);//and only if ISMEAR=-5
		break;
	case 1://Line Mode
		vasp_gui.calc.kpoints_mode=VKP_LINE;
		gtk_widget_set_sensitive(vasp_gui.kpoints_coord,TRUE);//need coord
		gtk_widget_set_sensitive(vasp_gui.kpoints_nkpts,TRUE);//need nkpts
		gtk_widget_set_sensitive(vasp_gui.kpoints_cart,TRUE);//can be cart
		break;
	case 3://Gamma (M&P)
		vasp_gui.calc.kpoints_mode=VKP_GAMMA;
		gtk_widget_set_sensitive(vasp_gui.kpoints_kx,TRUE);//require kx
		gtk_widget_set_sensitive(vasp_gui.kpoints_ky,TRUE);//require ky
		gtk_widget_set_sensitive(vasp_gui.kpoints_kz,TRUE);//require kz
		gtk_widget_set_sensitive(vasp_gui.kpoints_sx,TRUE);//require sx
		gtk_widget_set_sensitive(vasp_gui.kpoints_sy,TRUE);//require sy
		gtk_widget_set_sensitive(vasp_gui.kpoints_sz,TRUE);//require sz
		vasp_gui.calc.kpoints_nkpts=0;
		gtk_entry_set_text(GTK_ENTRY(vasp_gui.kpoints_nkpts),g_strdup_printf("%i",vasp_gui.calc.kpoints_nkpts));
		gtk_widget_set_sensitive(vasp_gui.kpoints_nkpts,FALSE);//NKPTS is 0 and disable
		break;
	case 4://Monkhorst-Pack classic
		vasp_gui.calc.kpoints_mode=VKP_MP;
		gtk_widget_set_sensitive(vasp_gui.kpoints_kx,TRUE);//require kx
		gtk_widget_set_sensitive(vasp_gui.kpoints_ky,TRUE);//require ky
		gtk_widget_set_sensitive(vasp_gui.kpoints_kz,TRUE);//require kz
		gtk_widget_set_sensitive(vasp_gui.kpoints_sx,TRUE);//require sx
		gtk_widget_set_sensitive(vasp_gui.kpoints_sy,TRUE);//require sy
		gtk_widget_set_sensitive(vasp_gui.kpoints_sz,TRUE);//require sz
		vasp_gui.calc.kpoints_nkpts=0;
		gtk_entry_set_text(GTK_ENTRY(vasp_gui.kpoints_nkpts),g_strdup_printf("%i",vasp_gui.calc.kpoints_nkpts));
		gtk_widget_set_sensitive(vasp_gui.kpoints_nkpts,FALSE);//NKPTS is 0 and disable
		break;
	case 5://Basis set definition
		vasp_gui.calc.kpoints_mode=VKP_BASIS;
		gtk_widget_set_sensitive(vasp_gui.kpoints_coord,TRUE);//need coord
		gtk_widget_set_sensitive(vasp_gui.kpoints_nkpts,TRUE);//need nkpts
		gtk_widget_set_sensitive(vasp_gui.kpoints_cart,TRUE);//can be cart
		break;
	case 2://Automatic (M&P)
	default:
		vasp_gui.calc.kpoints_mode=VKP_AUTO;
		gtk_widget_set_sensitive(vasp_gui.kpoints_kx,TRUE);//AUTO only require kx
		vasp_gui.calc.kpoints_nkpts=0;
		gtk_entry_set_text(GTK_ENTRY(vasp_gui.kpoints_nkpts),g_strdup_printf("%i",vasp_gui.calc.kpoints_nkpts));
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
	gtk_entry_set_text(GTK_ENTRY(vasp_gui.kpoints_i),g_strdup_printf("%i",index+1));
	if (g_ascii_strcasecmp(text,"ADD kpoint") == 0) {
                /*no need to update anything but index*/
	} else {
		/*update index,x,y,z,w*/
		wt=0.0;
		if(vasp_gui.calc.kpoints_mode==VKP_LINE) sscanf(text,"%lf %lf %lf ! kpoint: %*i",&x,&y,&z);
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
		if(vasp_gui.calc.kpoints_mode==VKP_LINE)
			gtk_combo_box_insert_text(GTK_COMBO_BOX(vasp_gui.kpoints_kpts),index,
                        	g_strdup_printf("%.8lf %.8lf %.8lf ! kpoint: %i",
                        	vasp_gui.calc.kpoints_x,vasp_gui.calc.kpoints_y,vasp_gui.calc.kpoints_z,index+1));
		else
			gtk_combo_box_insert_text(GTK_COMBO_BOX(vasp_gui.kpoints_kpts),index,
				g_strdup_printf("%.8lf %.8lf %.8lf %.8lf ! kpoint: %i",
				vasp_gui.calc.kpoints_x,vasp_gui.calc.kpoints_y,vasp_gui.calc.kpoints_z,vasp_gui.calc.kpoints_w,index+1));
        }else{
                /*we are goign to modify some data*/
		if(vasp_gui.calc.kpoints_mode==VKP_LINE)
			gtk_combo_box_insert_text(GTK_COMBO_BOX(vasp_gui.kpoints_kpts),index,
				g_strdup_printf("%.8lf %.8lf %.8lf ! kpoint: %i",
				vasp_gui.calc.kpoints_x,vasp_gui.calc.kpoints_y,vasp_gui.calc.kpoints_z,index+1));
		else
			gtk_combo_box_insert_text(GTK_COMBO_BOX(vasp_gui.kpoints_kpts),index,
				g_strdup_printf("%.8lf %.8lf %.8lf %.8lf ! kpoint: %i",
				vasp_gui.calc.kpoints_x,vasp_gui.calc.kpoints_y,vasp_gui.calc.kpoints_z,vasp_gui.calc.kpoints_w,index+1));
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
	if(vasp_gui.calc.ismear!=-5) vasp_gui.calc.kpoints_tetra=FALSE;/*ISMEAR=-5 required*/
	if(vasp_gui.calc.kpoints_mode!=VKP_MAN) vasp_gui.calc.kpoints_tetra=FALSE;/*manual generation must be selected*/
	if(vasp_gui.calc.kpoints_tetra==FALSE) gtk_widget_set_sensitive(vasp_gui.kpoints_tetra,FALSE);
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
			vasp_gui.calc.tetra_w,vasp_gui.calc.tetra_a,vasp_gui.calc.tetra_b,vasp_gui.calc.tetra_c,vasp_gui.calc.tetra_d,index));
        }else{  
                /*we are goign to modify some data*/
		gtk_combo_box_insert_text(GTK_COMBO_BOX(vasp_gui.tetra),index,
			g_strdup_printf("%lf %i %i %i %i ! tetrahedron: %i",
			vasp_gui.calc.tetra_w,vasp_gui.calc.tetra_a,vasp_gui.calc.tetra_b,vasp_gui.calc.tetra_c,vasp_gui.calc.tetra_d,index));
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
	fp=fopen(vasp_gui.calc.potcar_file,"r");
        if (!fp) {
		if(vasp_gui.calc.potcar_species!=NULL) g_free(vasp_gui.calc.potcar_species);
		vasp_gui.calc.potcar_species=NULL;
                return;
        }
	/*file is opened*/
	sym[2]='\0';
	if(vasp_gui.calc.potcar_species!=NULL) g_free(vasp_gui.calc.potcar_species);
	vasp_gui.calc.potcar_species=g_strdup_printf(" ");
	line = file_read_line(fp);
	sscanf(line," %c%c%c%*s %*s",&(sym[0]),&(sym[1]),&(sym[2]));
	if((sym[0]=='P')&&(sym[1]=='A')&&(sym[2]=='W')) {
		vasp_gui.calc.have_paw=TRUE;
		gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(vasp_gui.have_paw),TRUE);
	}
	sym[2]='\0';
	while(line){
		if (strstr(line,"TITEL") != NULL) {
			sscanf(line," TITEL = %*s %c%c%*s",&(sym[0]),&(sym[1]));
			tamp=g_strdup_printf("%s %s",vasp_gui.calc.potcar_species,sym);
			g_free(vasp_gui.calc.potcar_species);
			vasp_gui.calc.potcar_species=tamp;
		}
		g_free(line);
		line = file_read_line(fp);
	}
	fclose(fp);
	/*update vasp_gui.simple_species and vasp_gui.potcar_species*/
	gtk_entry_set_text(GTK_ENTRY(vasp_gui.simple_species),g_strdup_printf("%s",vasp_gui.calc.potcar_species));
	gtk_entry_set_text(GTK_ENTRY(vasp_gui.potcar_species),g_strdup_printf("%s",vasp_gui.calc.potcar_species));
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
	gtk_entry_set_text(GTK_ENTRY(vasp_gui.simple_potcar),g_strdup_printf("%s",filename));
        gtk_entry_set_text(GTK_ENTRY(vasp_gui.potcar_file),g_strdup_printf("%s",filename));
/*dont forget to sync!!*/
        VASP_REG_TEXT(potcar_file);
        g_free (filename);
        potcar_get_species();
    }
  }
  gtk_widget_destroy (GTK_WIDGET(file_chooser));
}
/*****************************************************/
/* register a specific pseudopotential for a species */
/*****************************************************/
void apply_species_flavor(GtkButton *button, gpointer data){
	gchar *text=gtk_combo_box_text_get_active_text(GTK_COMBO_BOX_TEXT(vasp_gui.species_flavor));
	gchar *symbol=g_strdup(text);
	gchar *flavor=g_strdup(text);
	gchar sym[3];
	gint idx;
	gchar *result=NULL;
	if(text==NULL) return;/*nothing selected*/
	if(text[0]=='\0') return;/*nothing selected*/
	sscanf(text,"%[^:]: %s",symbol,flavor);
	text=g_strdup_printf("%s",vasp_gui.calc.potcar_species_flavor);
	g_free(vasp_gui.calc.potcar_species_flavor);
	sym[2]='\0';
	if(g_ascii_islower(symbol[1])){/*2 letters*/
		sscanf(symbol,"%c%c",&(sym[0]),&(sym[1]));
		idx=0;
		while((g_ascii_strncasecmp(&(text[idx]),sym,2)!=0)&&(text[idx]!='\0')) idx++;
		if(text[idx]=='\0'){
			fprintf(stderr,"ERROR in setting POTCAR flavor.\n");
			return;/*this should never happen*/
		}
		if(idx!=0) {
			text[idx]='\0';
			result=g_strdup_printf("%s%s",text,flavor);
		}else{
			result=g_strdup_printf("%s",flavor);
		}
		idx++;
		while((text[idx]!=' ')&&(text[idx]!='\0')) idx++;
		if(text[idx]=='\0'){
			vasp_gui.calc.potcar_species_flavor=g_strdup_printf("%s",result);
		}else{
			vasp_gui.calc.potcar_species_flavor=g_strdup_printf("%s%s",result,&(text[idx]));
		}
	}else{/*if symbol is only 1 letter long, we need to cheat a little*/
		sscanf(symbol,"%c",&(sym[0]));sym[1]='\0';
		idx=0;
		while(text[idx]!='\0'){
			if((text[idx]==sym[0])&&((text[idx+1]==' ')||(text[idx+1]=='_')||text[idx+1]=='\0')) break;
			/*H is special*/
			if((text[0]=='H')&&((text[1]=='\0')||(text[1]=='.')||(text[1]=='1')||(text[1]=='_'))) break;
			idx++;
		}
		if(text[idx]=='\0'){
			fprintf(stderr,"ERROR in setting POTCAR flavor.\n");
			return;/*this should never happen*/
		}
		if(idx!=0) {
			text[idx]='\0';
			result=g_strdup_printf("%s%s",text,flavor);
		}else{
			result=g_strdup_printf("%s",flavor);
		}
		idx++;
		while((text[idx]!=' ')&&(text[idx]!='\0')) idx++;
		if(text[idx]=='\0'){
			vasp_gui.calc.potcar_species_flavor=g_strdup_printf("%s",result);
		}else{
			vasp_gui.calc.potcar_species_flavor=g_strdup_printf("%s%s",result,&(text[idx]));
		}
	}
	g_free(result);
	gtk_entry_set_text(GTK_ENTRY(vasp_gui.potcar_species_flavor),g_strdup_printf("%s",vasp_gui.calc.potcar_species_flavor));
}
/******************************************************/
/* register a species and a flavor from POTCAR FOLDER */
/******************************************************/
#define DEBUG_POTCAR_FOLDER 0
void potcar_folder_register(gchar *item){
	gchar elem[3];
	elem[0]=item[0];
	if(g_ascii_islower (item[1])) elem[1]=item[1];
	else elem[1]='\0';
	elem[2]='\0';
#if DEBUG_POTCAR_FOLDER
	fprintf(stdout,"#DBG: REG_ %s as %s \n",item,elem);
#endif
	VASP_COMBOBOX_ADD(vasp_gui.species_flavor,g_strdup_printf("%s: %s",elem,item));
}
/**************************************/
/* get information from POTCAR folder */
/**************************************/
void potcar_folder_get_info(){
	/*using POTCAR files, given they are not stored in a ZCOMPRESS format*/
	const gchar *path=gtk_entry_get_text(GTK_ENTRY(vasp_gui.potcar_folder));
	GSList *folders=file_dir_list(path,FALSE);
	gchar *item;
	gchar sym[3];
	vasp_poscar_sync();
#if DEBUG_POTCAR_FOLDER
	fprintf(stdout,"#DBG: looking for %s\n",vasp_gui.calc.species_symbols);
#endif
	/*wipe all entry in vasp_gui.species_flavor*/
	gtk_list_store_clear(GTK_LIST_STORE(gtk_combo_box_get_model(GTK_COMBO_BOX(vasp_gui.species_flavor))));
	sym[2]='\0';
	for ( ; folders ; folders=g_slist_next(folders)){
		item=g_strdup_printf("%s",(char *)folders->data);
		/*check if this item is of interest*/
		if(g_ascii_islower (item[1])){
			sym[0]=item[0];
			sym[1]=item[1];
			if(g_strrstr(vasp_gui.calc.species_symbols,sym)!=NULL){
				/* we have a hit, check validity */
				if((item[2]=='\0')||(item[2]=='_')) potcar_folder_register(item);
			}
		}else{
			gint idx=0;
			while(vasp_gui.calc.species_symbols[idx]!='\0'){
				if(item[0]==vasp_gui.calc.species_symbols[idx]){/*one match*/
					if(!g_ascii_islower (vasp_gui.calc.species_symbols[idx+1])){/*is not 2 letter wide*/
						if((item[0]=='H')&&((item[1]=='.')||(item[1]=='1')))
							potcar_folder_register(item);/*H is special*/
						if((item[1]=='\0')||(item[1]=='_')) potcar_folder_register(item);
					}
				}
				idx++;
			}
		}
		g_free(item);
	}
	vasp_gui.calc.potcar_species=g_strdup(vasp_gui.calc.species_symbols);/*just copy is enough for now*/
	vasp_gui.calc.potcar_species_flavor=g_strdup(vasp_gui.calc.species_symbols);/*same here*/
        gtk_entry_set_text(GTK_ENTRY(vasp_gui.simple_species),g_strdup_printf("%s",vasp_gui.calc.potcar_species));
        gtk_entry_set_text(GTK_ENTRY(vasp_gui.potcar_species),g_strdup_printf("%s",vasp_gui.calc.potcar_species));
	gtk_entry_set_text(GTK_ENTRY(vasp_gui.potcar_species_flavor),g_strdup_printf("%s",vasp_gui.calc.potcar_species_flavor));
	gtk_combo_box_set_active(GTK_COMBO_BOX(vasp_gui.species_flavor),0);
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
	if(vasp_gui.calc.potcar_folder!=NULL) g_free(vasp_gui.calc.potcar_folder);
	vasp_gui.calc.potcar_folder=g_strdup_printf("%s",foldername);
        g_free (foldername);
	potcar_folder_get_info();
    }
  }
  gtk_widget_destroy (GTK_WIDGET(file_chooser));
}
/*************************************/
/* Select POTCAR information by file */
/*************************************/
void toggle_potcar_file(GtkWidget *w, GtkWidget *box){
	vasp_gui.have_potcar_folder=FALSE;
	gtk_widget_set_sensitive(vasp_gui.potcar_folder,FALSE);
	gtk_widget_set_sensitive(vasp_gui.potcar_folder_button,FALSE);
	gtk_widget_set_sensitive(vasp_gui.species_flavor,FALSE);
	gtk_widget_set_sensitive(vasp_gui.species_button,FALSE);
	gtk_widget_set_sensitive(vasp_gui.potcar_file,TRUE);
	gtk_widget_set_sensitive(vasp_gui.potcar_file_button,TRUE);
}
/*************************************/
/* Select POTCAR information by path */
/*************************************/
void toggle_potcar_folder(GtkWidget *w, GtkWidget *box){
	vasp_gui.have_potcar_folder=TRUE;
	gtk_widget_set_sensitive(vasp_gui.potcar_file,FALSE);
	gtk_widget_set_sensitive(vasp_gui.potcar_file_button,FALSE);
	gtk_widget_set_sensitive(vasp_gui.potcar_folder,TRUE);
	gtk_widget_set_sensitive(vasp_gui.potcar_folder_button,TRUE);
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
	np=(gint)vasp_gui.calc.job_nproc;
	ncore=(gint)vasp_gui.calc.ncore;
	kpar=(gint)vasp_gui.calc.kpar;
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
                sscanf(text,"%*f %*f %*f %*c %*c %*c ! atom: %*i (%[^)])",&(vasp_gui.calc.poscar_symbol[0]));
                jdx=0;
                gtk_combo_box_set_active(GTK_COMBO_BOX(vasp_gui.poscar_atoms),jdx);
                text2=gtk_combo_box_text_get_active_text(GTK_COMBO_BOX_TEXT(vasp_gui.poscar_atoms));
                while(g_ascii_strcasecmp(text2,"ADD atom")!=0){
                        /**/
                        sscanf(text2,"%*f %*f %*f %*c %*c %*c ! atom: %*i (%[^)])",&(symbol[0]));
                        if(g_ascii_strcasecmp(vasp_gui.calc.poscar_symbol,symbol) == 0){
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
/* pass_3: SETUP unexposed properties (fix atom index)*/
        idx=0;jdx=0;vasp_gui.calc.electron_total=0;
        vasp_gui.calc.species_total=1;
        gtk_combo_box_set_active(GTK_COMBO_BOX(vasp_gui.poscar_atoms),idx);
        text=gtk_combo_box_text_get_active_text(GTK_COMBO_BOX_TEXT(vasp_gui.poscar_atoms));
        sscanf(text,"%*f %*f %*f %*c %*c %*c ! atom: %*i (%[^)])",&(vasp_gui.calc.poscar_symbol[0]));
        tamp=g_strdup_printf("  %s",vasp_gui.calc.poscar_symbol);
        tamp2=g_strdup_printf(" ");
        while(g_ascii_strcasecmp(text,"ADD atom")!=0){
                /*get each line*/
                sscanf(text,"%*f %*f %*f %*c %*c %*c ! atom: %*i (%[^)])",&(symbol[0]));
		text2=g_strdup(text);/*to get enough space*/
		sscanf(text,"%[^!]%*s",text2);
                gtk_combo_box_insert_text(GTK_COMBO_BOX(vasp_gui.poscar_atoms),idx,g_strdup_printf("%s! atom: %i (%s)",text2,idx,symbol));
                gtk_combo_box_text_remove(GTK_COMBO_BOX_TEXT(vasp_gui.poscar_atoms),idx+1);
		g_free(text2);
		vasp_gui.calc.electron_total+=elem_symbol_test(symbol);
                if(g_ascii_strcasecmp(vasp_gui.calc.poscar_symbol,symbol) != 0) {
                        /*new species*/
                        text=g_strdup_printf("%s %s",tamp,symbol);
                        g_free(tamp);
                        tamp=text;
                        text=g_strdup_printf("%s %i",tamp2,jdx);
                        g_free(tamp2);
                        tamp2=text;
                        vasp_gui.calc.species_total++;
                        sprintf(vasp_gui.calc.poscar_symbol,"%s",symbol);
                        jdx=0;
                }
                idx++;
                jdx++;
                gtk_combo_box_set_active(GTK_COMBO_BOX(vasp_gui.poscar_atoms),idx);
                text=gtk_combo_box_text_get_active_text(GTK_COMBO_BOX_TEXT(vasp_gui.poscar_atoms));
        }
        vasp_gui.calc.atoms_total=idx;/*so we have the total number of atoms*/
        /*we need to get the last jdx*/
        text=g_strdup_printf("%s %i",tamp2,jdx);
        g_free(tamp2);
        tamp2=text;
        if(vasp_gui.calc.species_symbols!=NULL) g_free(vasp_gui.calc.species_symbols);
        if(vasp_gui.calc.species_numbers!=NULL) g_free(vasp_gui.calc.species_numbers);
        vasp_gui.calc.species_symbols=g_strdup(tamp);g_free(tamp);
        vasp_gui.calc.species_numbers=g_strdup(tamp2);g_free(tamp2);
	vasp_gui.poscar_dirty=FALSE;
}

/************************/
/* Switch notebook page */
/************************/
void vasp_gui_page_switch(GtkNotebook *notebook,GtkWidget *page,guint page_num,gpointer user_data){
        /* some (few) values need to be updated from a page to another*/
	vasp_gui.cur_page=page_num;/*TODO: for (adv.)*/
	if(page_num==VASP_PAGE_SIMPLIFIED){
		vasp_poscar_sync();/*we need to know species_symbols*/
		gtk_entry_set_text(GTK_ENTRY(vasp_gui.simple_poscar),g_strdup_printf("%s",vasp_gui.calc.species_symbols));
		gtk_spin_button_set_value(GTK_SPIN_BUTTON(vasp_gui.simple_np),vasp_gui.calc.job_nproc);
		gtk_spin_button_set_value(GTK_SPIN_BUTTON(vasp_gui.simple_ncore),vasp_gui.calc.ncore);
		gtk_spin_button_set_value(GTK_SPIN_BUTTON(vasp_gui.simple_kpar),vasp_gui.calc.kpar);
	} else if(page_num==VASP_PAGE_ELECTRONIC){
		gtk_entry_set_text(GTK_ENTRY(vasp_gui.ismear),g_strdup_printf("%i",vasp_gui.calc.ismear));
		gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(vasp_gui.kgamma),vasp_gui.calc.kgamma);
		gtk_entry_set_text(GTK_ENTRY(vasp_gui.kspacing),g_strdup_printf("%.4lf",vasp_gui.calc.kspacing));
	} else if (page_num==VASP_PAGE_IONIC){
		vasp_poscar_sync();/*arrange poscar atom values*/
	} else if (page_num==VASP_PAGE_KPOINTS){
		gtk_entry_set_text(GTK_ENTRY(vasp_gui.ismear_3),g_strdup_printf("%i",vasp_gui.calc.ismear));
		gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(vasp_gui.kgamma_3),vasp_gui.calc.kgamma);
		gtk_entry_set_text(GTK_ENTRY(vasp_gui.kspacing_3),g_strdup_printf("%.4lf",vasp_gui.calc.kspacing));
		gtk_combo_box_set_active(GTK_COMBO_BOX(vasp_gui.kpoints_mode),vasp_gui.calc.kpoints_mode);
                gtk_spin_button_set_value(GTK_SPIN_BUTTON(vasp_gui.kpoints_kx),vasp_gui.calc.kpoints_kx);
                gtk_spin_button_set_value(GTK_SPIN_BUTTON(vasp_gui.kpoints_ky),vasp_gui.calc.kpoints_ky);
                gtk_spin_button_set_value(GTK_SPIN_BUTTON(vasp_gui.kpoints_kz),vasp_gui.calc.kpoints_kz);
		toogle_tetra();/*in case we need*/
	} else if (page_num==VASP_PAGE_POTCAR){
		vasp_poscar_sync();/*we need to know species_symbols*/
		gtk_entry_set_text(GTK_ENTRY(vasp_gui.poscar_species),g_strdup_printf("%s",vasp_gui.calc.species_symbols));
	} else if (page_num==VASP_PAGE_EXEC){
		gtk_spin_button_set_value(GTK_SPIN_BUTTON(vasp_gui.job_nproc),vasp_gui.calc.job_nproc);
		gtk_spin_button_set_value(GTK_SPIN_BUTTON(vasp_gui.ncore),vasp_gui.calc.ncore);
		gtk_spin_button_set_value(GTK_SPIN_BUTTON(vasp_gui.kpar),vasp_gui.calc.kpar);
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
	if(!vasp_gui.calc.auto_elec){
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
void calc_to_poscar(FILE *output,vasp_calc_struct calc){
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
	fprintf(output,"%s\n",vasp_gui.calc.species_symbols);
	/*NUM ATOM PER SPECIES*/
	fprintf(output,"%s\n",vasp_gui.calc.species_numbers);
	if(calc.poscar_sd) fprintf(output,"Selective dynamics\n");
	if(calc.poscar_direct) fprintf(output,"Direct\n");
	else fprintf(output,"Cartesian\n");
	idx=0;
        gtk_combo_box_set_active(GTK_COMBO_BOX(vasp_gui.poscar_atoms),idx);
        text=gtk_combo_box_text_get_active_text(GTK_COMBO_BOX_TEXT(vasp_gui.poscar_atoms));
	sscanf(text,"%*f %*f %*f %*c %*c %*c ! atom: %*i (%[^)])",&(vasp_gui.calc.poscar_symbol[0]));
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
void calc_to_kpoints(FILE *output,vasp_calc_struct calc){
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
	if((vasp_gui.calc.kpoints_mode==VKP_MAN)||(vasp_gui.calc.kpoints_mode==VKP_LINE)) fprintf(output,"%i\n",calc.kpoints_nkpts);
	else fprintf(output,"0\n");
	
	switch (vasp_gui.calc.kpoints_mode){
	case VKP_LINE:
		fprintf(output,"Line-mode\n");
	case VKP_BASIS:
        case VKP_MAN:
		if(vasp_gui.calc.kpoints_cart) fprintf(output,"Cartesian\n");
		else fprintf(output,"Reciprocal\n");
		break;
        case VKP_GAMMA:
		fprintf(output,"Gamma\n");
		fprintf(output,"%i %i %i\n",(gint)vasp_gui.calc.kpoints_kx,(gint)vasp_gui.calc.kpoints_ky,(gint)vasp_gui.calc.kpoints_kz);
		if((vasp_gui.calc.kpoints_sx!=0.0)&&(vasp_gui.calc.kpoints_sy!=0.0)&&(vasp_gui.calc.kpoints_sz!=0.0))
			fprintf(output,"%lf %lf %lf\n",vasp_gui.calc.kpoints_sx,vasp_gui.calc.kpoints_sy,vasp_gui.calc.kpoints_sz);
		return;
        case VKP_MP:
		fprintf(output,"Monkhorst-Pack\n");
		fprintf(output,"%i %i %i\n",(gint)vasp_gui.calc.kpoints_kx,(gint)vasp_gui.calc.kpoints_ky,(gint)vasp_gui.calc.kpoints_kz);
		if((vasp_gui.calc.kpoints_sx!=0.0)&&(vasp_gui.calc.kpoints_sy!=0.0)&&(vasp_gui.calc.kpoints_sz!=0.0))
			fprintf(output,"%lf %lf %lf\n",vasp_gui.calc.kpoints_sx,vasp_gui.calc.kpoints_sy,vasp_gui.calc.kpoints_sz);
		return;
        case VKP_AUTO:
		fprintf(output,"Auto\n");
		fprintf(output,"%i\n",(gint)vasp_gui.calc.kpoints_kx);
		return;
        }
	/*print all coordinates*/
	idx=0;
	gtk_combo_box_set_active(GTK_COMBO_BOX(vasp_gui.kpoints_kpts),idx);
	text=gtk_combo_box_text_get_active_text(GTK_COMBO_BOX_TEXT(vasp_gui.kpoints_kpts));
        while(g_ascii_strcasecmp(text,"ADD kpoint")!=0){
		if(vasp_gui.calc.kpoints_mode==VKP_MAN) {
			sscanf(text,"%lf %lf %lf %lf ! kpoint: %*i",&x,&y,&z,&wt);
			fprintf(output,"%lf %lf %lf %lf ! kpoint: %i\n",x,y,z,wt,idx+1);
		} else {
			sscanf(text,"%lf %lf %lf ! kpoint: %*i",&x,&y,&z);
			fprintf(output,"%lf %lf %lf ! kpoint: %i\n",x,y,z,idx+1);
		}
		idx++;
		gtk_combo_box_set_active(GTK_COMBO_BOX(vasp_gui.kpoints_kpts),idx);
		text=gtk_combo_box_text_get_active_text(GTK_COMBO_BOX_TEXT(vasp_gui.kpoints_kpts));
	}
	/*return to previous postition*/
        gtk_combo_box_set_active(GTK_COMBO_BOX(vasp_gui.kpoints_kpts),index);
	/*print all tetrahedron (if any)*/
	if(!vasp_gui.calc.kpoints_tetra) return;
	fprintf(output,"Tetrahedron\n");
	fprintf(output,"%i %lf\n",vasp_gui.calc.tetra_total,vasp_gui.calc.tetra_volume);
	index=gtk_combo_box_get_active(GTK_COMBO_BOX(vasp_gui.tetra));
	idx=0;
	gtk_combo_box_set_active(GTK_COMBO_BOX(vasp_gui.tetra),idx);
	text=gtk_combo_box_text_get_active_text(GTK_COMBO_BOX_TEXT(vasp_gui.tetra));
	while(g_ascii_strcasecmp(text,"ADD tetrahedron")!=0){
		sscanf(text,"%lf %i %i %i %i ! tetrahedron: %*i",&wt,&a,&b,&c,&d);
		fprintf(output,"%lf %i %i %i %i ! tetrahedron: %i\n",wt,a,b,c,d,idx+1);
		idx++;
		gtk_combo_box_set_active(GTK_COMBO_BOX(vasp_gui.tetra),idx);
		text=gtk_combo_box_text_get_active_text(GTK_COMBO_BOX_TEXT(vasp_gui.tetra));
	}
	/*return to previous postition*/
	gtk_combo_box_set_active(GTK_COMBO_BOX(vasp_gui.tetra),index);
}
/*************************************/
/* append one POTCAR file to another */
/*************************************/
void append_potcar(FILE *src,FILE *dest){
	gchar *line;
	line = file_read_line(src);
	while(line){
		fprintf(dest,"%s",line);
		g_free(line);
		line = file_read_line(src);
	}
}
/***********************************/
/* copy or concatenation of POTCAR */
/***********************************/
void calc_to_potcar(FILE *output,vasp_calc_struct calc){
	FILE *org_potcar;
	gchar *line;
if(!vasp_gui.have_potcar_folder){
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
}else{/*use potcar_folder interface*/
	gboolean hasend;
	gchar *text=g_strdup(vasp_gui.calc.potcar_species_flavor);
	gchar *ptr;
	gchar *filename;
	gint idx;
	if(vasp_gui.calc.potcar_folder==NULL) return;
	if(vasp_gui.calc.potcar_species_flavor==NULL) return;
	hasend=FALSE;
	idx=0;
	while((text[idx]==' ')) {/*skip trailing space*/
		if(text[idx]=='\0') {
			fprintf(stderr,"ERROR while reading POTCAR flavors!\n");
			g_free(text);
			return;
		}
		idx++;
	}
	ptr=&(text[idx]);
	while(!hasend){
		while((text[idx]!=' ')){
			if(text[idx]=='\0') {
				hasend=TRUE;
				break;	
			}
			idx++;
		}
		text[idx]='\0';
		idx++;
		filename=g_strdup_printf("%s/%s/POTCAR",vasp_gui.calc.potcar_folder,ptr);
		org_potcar=fopen(filename,"r");
		if(!org_potcar) {
			fprintf(stderr,"#ERR: can't open file %s for READ!\n",calc.potcar_file);
			return;
		}
		append_potcar(org_potcar,output);
		fclose(org_potcar);
		ptr=&(text[idx]);
	}
	g_free(text);
	vasp_gui.calc.potcar_file=g_strdup_printf("%s/POTCAR",vasp_gui.calc.job_path);
}
	/*that's all folks*/
}
/***************************/
/* save current parameters */
/***************************/
gint save_vasp_calc(){
#define DEBUG_VASP_SAVE 0
	FILE *save_fp;
	gchar *filename;
	/*1-synchronize*/
	vasp_gui_sync();
	vasp_poscar_sync();
	/*2-output*/
	//INCAR
	filename=g_strdup_printf("%s/INCAR",vasp_gui.calc.job_path);
	save_fp=fopen(filename,"w");
	if(!save_fp){
		fprintf(stderr,"#ERR: can't open file INCAR for WRITE!\n");
		return -1;
	}
	vasp_calc_to_incar(save_fp,vasp_gui.calc);
	fclose(save_fp);
	g_free(filename);
	//POSCAR
	filename=g_strdup_printf("%s/POSCAR",vasp_gui.calc.job_path);
	save_fp=fopen(filename,"w");
	if(!save_fp){
		fprintf(stderr,"#ERR: can't open file POSCAR for WRITE!\n");
		return -1;
	}
	calc_to_poscar(save_fp,vasp_gui.calc);
	fclose(save_fp);
	g_free(filename);
	//KPOINTS
	filename=g_strdup_printf("%s/KPOINTS",vasp_gui.calc.job_path);
	save_fp=fopen(filename,"w");
	if(!save_fp){
		fprintf(stderr,"#ERR: can't open file KPOINTS for WRITE!\n");
		return -1;
	}
	calc_to_kpoints(save_fp,vasp_gui.calc);
	fclose(save_fp);
	g_free(filename);
	//POTCAR
	filename=g_strdup_printf("%s/POTCAR",vasp_gui.calc.job_path);
	save_fp=fopen(filename,"w");
	if(!save_fp){
		fprintf(stderr,"#ERR: can't open file KPOINTS for WRITE!\n");
		return -1;
	}
	calc_to_potcar(save_fp,vasp_gui.calc);
	fclose(save_fp);
	g_free(filename);
	/*debug print*/
#if DEBUG_VASP_SAVE
	vasp_gui.calc_to_incar(stdout,vasp_gui.calc);
	calc_to_poscar(stdout,vasp_gui.calc);
	calc_to_kpoints(stdout,vasp_gui.calc);

#endif
	return 0;
}
/*********************************/
/* Cleanup calculation structure */
/*********************************/
void vasp_cleanup(){
	/*we don't have to free anything.. everything will be gone after the dialog is closed*/
}
/*****************************************/
/* Execute or enqueue a vasp calculation */
/*****************************************/
void run_vasp_exec(vasp_exec_struct *vasp_exec,struct task_pak *task){
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
	if((*vasp_exec).job_nproc<2) cmd = g_strdup_printf("%s > vasp.log",(*vasp_exec).job_vasp_exe);
	else cmd = g_strdup_printf("%s -np %i %s > vasp.log",(*vasp_exec).job_mpirun,(gint)(*vasp_exec).job_nproc,(*vasp_exec).job_vasp_exe);
#endif
	cwd=sysenv.cwd;/*push*/
	sysenv.cwd=g_strdup_printf("%s",(*vasp_exec).job_path);
	task_sync(cmd);
	g_free(sysenv.cwd);
	sysenv.cwd=cwd;/*pull*/
	g_free(cmd);
	(*vasp_exec).have_result=TRUE;
}
/*****************************/
/* cleanup task: load result */
/*****************************/
void cleanup_vasp_exec(vasp_exec_struct *vasp_exec){
	/*VASP process has exit, try to load the result*/
	struct model_pak *result_model;
	gchar *filename;
	/*sync_ wait for result?*/
	while(!(*vasp_exec).have_result) usleep(500*1000);/*sleep 500ms until job is done*/
	/*init a model*/
	result_model=model_new();
        model_init(result_model);
	/*put result into it*/
	filename=g_strdup_printf("%s/vasprun.xml",(*vasp_exec).job_path);
	/*TODO: detect if a calculation have failed*/
	file_load(filename,result_model);/*TODO: load result without annoying tree_select_active*/
	model_prep(result_model);
	tree_model_add(result_model);
	tree_model_refresh(result_model);
	canvas_shuffle();
	redraw_canvas(ALL);/*ALL: necessary?*/
/*just wipe the structure*/
	sysenv.vasp_calc_list=g_slist_remove(sysenv.vasp_calc_list,vasp_exec);/*does not free, does it?*/
	g_free(vasp_exec);
}
/******************************/
/* Enqueue a vasp calculation */
/******************************/
void exec_calc(){
	vasp_exec_struct *vasp_exec;
	/*this will sync then enqueue a VASP calculation*/
	if(save_vasp_calc()) return;/*sync and save all file*/
/*copy structure to the list*/
	vasp_exec=g_malloc(sizeof(vasp_exec_struct));
	vasp_exec->job_id=g_slist_length(sysenv.vasp_calc_list);
	vasp_exec->have_result=FALSE;
	vasp_exec->have_gui=TRUE;
	vasp_exec->job_vasp_exe=g_strdup_printf("%s",vasp_gui.calc.job_vasp_exe);
	vasp_exec->job_mpirun=g_strdup_printf("%s",vasp_gui.calc.job_mpirun);
	vasp_exec->job_path=g_strdup_printf("%s",vasp_gui.calc.job_path);
	vasp_exec->job_nproc=vasp_gui.calc.job_nproc;
	/*prepend to calc list*/
	sysenv.vasp_calc_list = g_slist_prepend (sysenv.vasp_calc_list,vasp_exec);
	/*launch vasp in a task*/
	gtk_widget_set_sensitive(vasp_gui.button_save,FALSE);
	gtk_widget_set_sensitive(vasp_gui.button_exec,FALSE);
	task_new("VASP", &run_vasp_exec,vasp_exec,&cleanup_vasp_exec,vasp_exec,sysenv.active_model);
	/*when task is launched, close dialog*/
	vasp_cleanup();
	gtk_widget_destroy(vasp_gui.window);
}
/********************************/
/* Quit vasp calculation dialog */
/********************************/
void quit_vasp_gui(GtkWidget *w, gpointer data){
	struct model_pak *model=data;
	vasp_cleanup();
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
/* do we have a vasprun.xml model?  */
	if (data->id == VASP) {
		vasp_gui.have_xml=TRUE;
		gui_vasp_init();
		vasprun_update(data->filename,&vasp_gui.calc);/*initialize according to vasprun.xml*/
	} else {
		vasp_gui.have_xml=FALSE;/*initialize with default values*/
		gui_vasp_init();
		vasprun_update(NULL,&vasp_gui.calc);/*initialize with default values*/
	}
/* dialog setup */
	title = g_strdup_printf("VASP: %s", data->basename);
	dialog = dialog_request(CVASP, title, NULL, vasp_cleanup,data);
	g_free(title);
	vasp_gui.window = dialog_window(dialog);
/* --- Outside of notebook */
	GUI_FRAME_WINDOW(vasp_gui.window,frame);
	GUI_VBOX_FRAME(frame,vbox);
/* --- MODEL NAME */
	GUI_LINE_BOX(vbox,hbox);
	GUI_LABEL_BOX(hbox,label,"MODEL NAME");
	GUI_TEXT_ENTRY(hbox,vasp_gui.name,data->basename,TRUE);
GUI_TOOLTIP(vasp_gui.name,"The model name will be used in VASP files\nas well as GDIS display.");
/* --- CONNECTED XML */
	GUI_LINE_BOX(vbox,hbox);
	GUI_LABEL_BOX(hbox,label,"CONNECTED XML");
	GUI_TEXT_ENTRY(hbox,vasp_gui.file_entry,data->filename,TRUE);
	if(data->id!=VASP) GUI_ENTRY_TEXT(vasp_gui.file_entry,"");
GUI_TOOLTIP(vasp_gui.file_entry,"Use the result of a calculation (vasprun.xml)\nas a model to fill up settings of the current one.\nNOTE: important settings (such as NELM) will be wrong\nand all settings should be checked again. Carefully.");
	GUI_OPEN_BUTTON_BOX(hbox,button,load_vasprun_dialog);
/* create frame in box */
	GUI_FRAME_WINDOW(vasp_gui.window,frame);
/* create notebook in frame */
	GUI_NOTE_FRAME(frame,notebook);

/*------------------*/
/* page 0 -> PRESET */
/*------------------*/
	GUI_PAGE_NOTE(notebook,page,"PRESET");
/* --- general */
	GUI_FRAME_NOTE(page,frame,"General");
/* create a table in the frame*/
	GUI_TABLE_FRAME(frame,table,5,4);
/* 1st line */
	GUI_COMBOBOX_TABLE(table,vasp_gui.simple_calcul,"CALCUL:",0,1,0,1);
	GUI_COMBOBOX_ADD(vasp_gui.simple_calcul,"Energy (single point)");
	GUI_COMBOBOX_ADD(vasp_gui.simple_calcul,"DOS/BANDS (single point)");
	GUI_COMBOBOX_ADD(vasp_gui.simple_calcul,"Lattice Dynamics (opt.)");
	GUI_COMBOBOX_ADD(vasp_gui.simple_calcul,"Geometry (opt.)");
	GUI_COMBOBOX_ADD(vasp_gui.simple_calcul,"Molecular Dynamics (opt.)");
GUI_TOOLTIP(vasp_gui.simple_calcul,"What is the intended calculation type?");
	GUI_CHECK_TABLE(table,button,vasp_gui.simple_rgeom,NULL,"ROUGH GEOM.",2,3,0,1);/*not calling anything*/
GUI_TOOLTIP(button,"Is the current geometry:\nFAR from groundstate geometry\nor NOT-SO-FAR from it?");
	GUI_SPIN_TABLE(table,vasp_gui.simple_dim,vasp_gui.dimension,NULL,"DIM:",3,4,0,1);/*not calling anything*/
GUI_TOOLTIP(vasp_gui.simple_dim,"What is the dimension of system?\nThe DIM (0D, 1D, 2D, or 3D) will be used\nto determine some parameters in INCAR and KPOINTS\nNO CHANGE will be made to the current geometry!");
	GUI_SPIN_RANGE(vasp_gui.simple_dim,0.0,3.0);/*0D~3D*/
/* 2nd line */
	GUI_COMBOBOX_TABLE(table,vasp_gui.simple_system,"SYSTEM:",0,1,1,2);
	GUI_COMBOBOX_ADD(vasp_gui.simple_system,"Metal");
	GUI_COMBOBOX_ADD(vasp_gui.simple_system,"Semi-conducting");
	GUI_COMBOBOX_ADD(vasp_gui.simple_system,"Insulator");
	GUI_COMBOBOX_ADD(vasp_gui.simple_system,"Unknown");
GUI_TOOLTIP(vasp_gui.simple_system,"What is the electronic properties of the intended material?");
	GUI_TEXT_TABLE(table,vasp_gui.simple_poscar,vasp_gui.calc.species_symbols,"POSCAR:",1,4,1,2);
GUI_TOOLTIP(vasp_gui.simple_poscar,"Atomic symbols of the species\ndetected from the model geometry.");
	GUI_LOCK(vasp_gui.simple_poscar);
/* 3rd line */
	GUI_COMBOBOX_TABLE(table,vasp_gui.simple_kgrid,"KPT GRID:",0,1,2,3);
	GUI_COMBOBOX_ADD(vasp_gui.simple_kgrid,"COARSE");
	GUI_COMBOBOX_ADD(vasp_gui.simple_kgrid,"MEDIUM");
	GUI_COMBOBOX_ADD(vasp_gui.simple_kgrid,"FINE");
	GUI_COMBOBOX_ADD(vasp_gui.simple_kgrid,"ULTRA-FINE");
GUI_TOOLTIP(vasp_gui.simple_kgrid,"What is the intended k-point grid density?");
	GUI_TEXT_TABLE(table,vasp_gui.simple_species,vasp_gui.calc.potcar_species,"POTCAR:",1,4,2,3);
GUI_TOOLTIP(vasp_gui.simple_species,"From POTCAR file, the following species are detected.\nThis should match the POSCAR definition above.");
	GUI_LOCK(vasp_gui.simple_species);
/* 4th line */
	GUI_TEXT_TABLE(table,vasp_gui.simple_potcar,vasp_gui.calc.potcar_file,"POTCAR FILE:",0,3,3,4);
GUI_TOOLTIP(vasp_gui.simple_potcar,"Where is the POTCAR file for this calculation?");
	GUI_OPEN_BUTTON_TABLE(table,vasp_gui.simple_potcar_button,load_potcar_file_dialog,3,4,3,4);
/* 5th line */
	GUI_LEFT_LABEL_TABLE(table,"EXEC:",0,1,4,5);
	GUI_SPIN_TABLE(table,vasp_gui.simple_np,vasp_gui.calc.job_nproc,parallel_eq,"NP=",1,2,4,5);
GUI_TOOLTIP(vasp_gui.simple_np,"How many CPU(s) should be used for calculation?");
	GUI_SPIN_TABLE(table,vasp_gui.simple_ncore,vasp_gui.calc.ncore,parallel_eq,"NCORE=",2,3,4,5);
GUI_TOOLTIP(vasp_gui.simple_ncore,"How many CPU(s) should work on one band?");
	GUI_SPIN_TABLE(table,vasp_gui.simple_kpar,vasp_gui.calc.kpar,parallel_eq,"KPAR=",3,4,4,5);
GUI_TOOLTIP(vasp_gui.simple_kpar,"How many k-points should be worked on at a time?\nNote: each k-point is worked on by NCORE CPU(s).");
/* initialize */
	GUI_COMBOBOX_SETUP(vasp_gui.simple_calcul,0,NULL);/*not calling anything*/
	GUI_COMBOBOX_SETUP(vasp_gui.simple_system,0,NULL);/*not calling anything*/
	GUI_COMBOBOX_SETUP(vasp_gui.simple_kgrid,0,NULL);/*not calling anything*/
/* --- end frame */
/* --- DISCLAMER */
	GUI_FRAME_NOTE(page,frame,"DISCLAMER");
/* create a table in the frame*/
	GUI_TABLE_FRAME(frame,table,4,3);
/* 1st line */
	GUI_LABEL_TABLE(table,"This simplified interface is aimed at pre-loading parameters with default values",1,2,0,1);
/* 2nd line */
	GUI_LABEL_TABLE(table,"for a given calculation type. User should then check said defaults extremely carefuly...",1,2,1,2);
/* 3th line */
	GUI_LABEL_TABLE(table,"GENERATE ==>",0,1,3,4);
	GUI_APPLY_BUTTON_TABLE(table,vasp_gui.simple_apply,vasp_apply_simple,1,2,3,4);
GUI_TOOLTIP(vasp_gui.simple_apply,"Use this to obtain a preset of setting.\nIf calculation is already set using this will\noverwrite settings resulting in an unpredictable state.");
	GUI_LABEL_TABLE(table,"<== (and check...)",2,3,3,4);
/* --- end frame */
/* --- MESSAGE */
	GUI_FRAME_NOTE(page,frame,"MESSAGE");
/* create a vbox in frame*/
	GUI_VBOX_FRAME(frame,vbox);
/* create a textview in a frame */
	GUI_TEXTVIEW_BOX(vbox,vasp_gui.simple_message,vasp_gui.simple_message_buff);
/* --- end frame */


/*-----------------------*/
/* page 1 -> CONVERGENCE */
/*-----------------------*/
	GUI_PAGE_NOTE(notebook,page,"CONVERGENCE");
/* --- general */
	GUI_FRAME_NOTE(page,frame,"General");
/* create a table in the frame*/
	GUI_TABLE_FRAME(frame,table,5,4);
/* 1st line */
	GUI_COMBOBOX_TABLE(table,vasp_gui.prec,"PREC=",0,1,0,1);
	GUI_COMBOBOX_ADD(vasp_gui.prec,"Normal");
	GUI_COMBOBOX_ADD(vasp_gui.prec,"Single");
	GUI_COMBOBOX_ADD(vasp_gui.prec,"Accurate");
	GUI_COMBOBOX_ADD(vasp_gui.prec,"High");
	GUI_COMBOBOX_ADD(vasp_gui.prec,"Medium");
	GUI_COMBOBOX_ADD(vasp_gui.prec,"Low");
GUI_TOOLTIP(vasp_gui.prec,"PREC: Ch. 6.11 DEFAULT: Normal\nSets numerical accuracy by changing\nENCUT, NG{X,Y,Z}, NG{X,Y,Z}F, and ROPT.\n(can be override manually)");
	GUI_CHECK_TABLE(table,button,vasp_gui.calc.use_prec,prec_toggle,"USE PREC",1,2,0,1);
GUI_TOOLTIP(button,"Let PREC decides on ENCUT, NG{X,Y,Z}, NG{X,Y,Z}F, and ROPT.");
	GUI_ENTRY_TABLE(table,vasp_gui.encut,vasp_gui.calc.encut,"%.2lf","ENCUT=",2,3,0,1);
GUI_TOOLTIP(vasp_gui.encut,"ENCUT: Ch. 6.9 DEFAULT: MAX(ENMAX)\nPlane-wave basis cutoff energy (eV).\nDefault MAX(ENMAX) is taken from POTCAR file.");
	GUI_ENTRY_TABLE(table,vasp_gui.enaug,vasp_gui.calc.enaug,"%.2lf","ENAUG=",3,4,0,1);
GUI_TOOLTIP(vasp_gui.enaug,"ENAUG: Ch. 6.10 DEFAULT: EAUG\nKinetic cutoff for augmentation charges.\nDefault EAUG is taken from POTCAR file.");
	GUI_ENTRY_TABLE(table,vasp_gui.ediff,vasp_gui.calc.ediff,"%.2lE","EDIFF=",4,5,0,1);
GUI_TOOLTIP(vasp_gui.ediff,"EDIFF: Ch. 6.18 DEFAULT: 10^{-4}\nSet the criterion (eV) for electronic convergence.");
/* 2nd line */
	GUI_COMBOBOX_TABLE(table,vasp_gui.algo,"ALGO=",0,1,1,2);
	GUI_COMBOBOX_ADD(vasp_gui.algo,"NORMAL");
	GUI_COMBOBOX_ADD(vasp_gui.algo,"USE_IALGO");
	GUI_COMBOBOX_ADD(vasp_gui.algo,"VERYFAST");
	GUI_COMBOBOX_ADD(vasp_gui.algo,"FAST");
	GUI_COMBOBOX_ADD(vasp_gui.algo,"CONJUGATE");
	GUI_COMBOBOX_ADD(vasp_gui.algo,"ALL");
	GUI_COMBOBOX_ADD(vasp_gui.algo,"DAMPED");
	GUI_COMBOBOX_ADD(vasp_gui.algo,"SUBROT");
	GUI_COMBOBOX_ADD(vasp_gui.algo,"EIGENVAL");
	GUI_COMBOBOX_ADD(vasp_gui.algo,"NONE");
	GUI_COMBOBOX_ADD(vasp_gui.algo,"NOTHING");
	GUI_COMBOBOX_ADD(vasp_gui.algo,"EXACT");
	GUI_COMBOBOX_ADD(vasp_gui.algo,"DIAG");
GUI_TOOLTIP(vasp_gui.algo,"ALGO: Ch. 6.46 DEFAULT: Normal\nSelect the electronic minimization algorithm.\nIt is overrided by IALGO when specified.");
	GUI_CHECK_TABLE(table,button,vasp_gui.calc.ldiag,NULL,"LDIAG",1,2,1,2);/*not calling anything*/
GUI_TOOLTIP(button,"LDIAG: Ch. 6.47 DEFAULT: TRUE\nSwitch on the sub-space rotation.");
	GUI_ENTRY_TABLE(table,vasp_gui.nsim,vasp_gui.calc.nsim,"%i","NSIM=",2,3,1,2);
GUI_TOOLTIP(vasp_gui.nsim,"NSIM: Ch. 6.48 DEFAULT: 4\nNumber of bands optimized at a time\nin blocked minimization algorithm.");
	GUI_ENTRY_TABLE(table,vasp_gui.vtime,vasp_gui.calc.vtime,"%.4lf","TIME=",3,4,1,2);
GUI_TOOLTIP(vasp_gui.vtime,"TIME: Ch. 6.51 DEFAULT: 0.4\nSelect the trial time or value step\nfor minimization algorithm IALGO=5X.");
	GUI_ENTRY_TABLE(table,vasp_gui.iwavpr,vasp_gui.calc.iwavpr,"%i","IWAVPR=",4,5,1,2);
GUI_TOOLTIP(vasp_gui.iwavpr,"IWAVPR: Ch. 6.26 DEFAULT: 2(IBRION=0-2)\nDetermine how charge density/orbitals are predicted\nfrom one ionic configuration to the next.\nDefault is 0 if IBRION is different from 0-2.");
/* 3rd line */
	GUI_COMBOBOX_TABLE(table,vasp_gui.ialgo,"IALGO=",0,1,2,3);
	GUI_COMBOBOX_ADD(vasp_gui.ialgo,"2:FIXED ORB/1E");
	GUI_COMBOBOX_ADD(vasp_gui.ialgo,"3:FIXED ORB");
	GUI_COMBOBOX_ADD(vasp_gui.ialgo,"4:SUBROT ONLY");
	GUI_COMBOBOX_ADD(vasp_gui.ialgo,"5:STEEPEST");
	GUI_COMBOBOX_ADD(vasp_gui.ialgo,"6:CONJUGUATE GRADIENTS");
	GUI_COMBOBOX_ADD(vasp_gui.ialgo,"7:PRECOND. STEEPEST");
	GUI_COMBOBOX_ADD(vasp_gui.ialgo,"8:PRECOND. CG");
	GUI_COMBOBOX_ADD(vasp_gui.ialgo,"38:KOSUGI");
	GUI_COMBOBOX_ADD(vasp_gui.ialgo,"44:RMM STEEPEST");
	GUI_COMBOBOX_ADD(vasp_gui.ialgo,"46:RMM PRECOND.");
	GUI_COMBOBOX_ADD(vasp_gui.ialgo,"48:PRECOND. RMM");
	GUI_COMBOBOX_ADD(vasp_gui.ialgo,"53:VAR DAMPED MD");
	GUI_COMBOBOX_ADD(vasp_gui.ialgo,"54:VAR DAMPED QUENCH MD");
	GUI_COMBOBOX_ADD(vasp_gui.ialgo,"58:VAR PRECOND. CG");
	GUI_COMBOBOX_ADD(vasp_gui.ialgo,"90:EXACT");
GUI_TOOLTIP(vasp_gui.ialgo,"IALGO: Ch. 6.47 DEFAULT: 38\nSets the minimization algorithm number.\nIt is advised to use ALGO instead of IALGO\ndue to instabilities in extra IALGO algorithms.");
	GUI_CHECK_TABLE(table,button,vasp_gui.calc.auto_elec,elec_toggle,"AUTO",1,2,2,3);
GUI_TOOLTIP(button,"Automatically sets NSIM, TIME, IWAVPR, NBANDS, and NELECT");
	GUI_ENTRY_TABLE(table,vasp_gui.nbands,vasp_gui.calc.nbands,"%i","NBANDS=",2,3,2,3);
GUI_TOOLTIP(vasp_gui.nbands,"NBANDS: Ch. 6.5 DEFAULT: NELEC/2+NIONS/2\nNumber of bands in the calculation.\nDefault for spin-polarized calcualtions is 0.6*NELECT+NMAG");
	GUI_ENTRY_TABLE(table,vasp_gui.nelect,vasp_gui.calc.nelect,"%.2lf","NELECT=",3,4,2,3);
GUI_TOOLTIP(vasp_gui.nelect,"NELECT: Ch. 6.35 DEFAULT: Ne valence\nSets the number of electrons\nCan be use add an extra charge background\nas an homogeneous background-charge.");
	/*empty: col4*/
/* 4th line */
	/*empty: col0*/
	GUI_CHECK_TABLE(table,button,vasp_gui.calc.iniwav,NULL,"INIWAV",1,2,3,4);/*not calling anything*/
GUI_TOOLTIP(button,"INIWAV: Ch. 6.16 DEFAULT: 1\nDetermine how initial wavefunctions are filled\nThe only option is the default random filling.\nOther option (INIWAV=0) is usually not available.");
	GUI_ENTRY_TABLE(table,vasp_gui.istart,vasp_gui.calc.istart,"%i","ISTART=",2,3,3,4);
GUI_TOOLTIP(vasp_gui.istart,"ISTART: Ch. 6.14 DEFAULT 1\nDetermine whether WAVECAR should be read.\nDefault is 0 when no WAVECAR file is found.");
	GUI_ENTRY_TABLE(table,vasp_gui.icharg,vasp_gui.calc.icharg,"%i","ICHARG=",3,4,3,4);
GUI_TOOLTIP(vasp_gui.icharg,"ICHARG: Ch. 6.15 DEFAULT 2\nDetermine how the initial charge density is calculated\nDefault is 0 if ISTART is not 0.");
	/*empty: col4*/
/* initialize */
	GUI_COMBOBOX_SETUP(vasp_gui.prec,0,vasp_prec_selected);
	GUI_COMBOBOX_SETUP(vasp_gui.algo,0,vasp_algo_selected);
	GUI_COMBOBOX_SETUP(vasp_gui.ialgo,7,vasp_ialgo_selected);
	GUI_LOCK(vasp_gui.ialgo);
	prec_toggle(NULL,NULL);
	elec_toggle(NULL,NULL);
/* --- end frame */
/* --- mixing */
	GUI_FRAME_NOTE(page,frame,"Mixing");
/* create a table in the frame*/
	GUI_TABLE_FRAME(frame,table,5,3);
/* 1st line */
	GUI_COMBOBOX_TABLE(table,vasp_gui.mixer,"IMIX=",0,1,0,1);
	GUI_COMBOBOX_ADD(vasp_gui.mixer,"0:NONE");
	GUI_COMBOBOX_ADD(vasp_gui.mixer,"1:KERKER");
	GUI_COMBOBOX_ADD(vasp_gui.mixer,"2:TCHEBYCHEV");
	GUI_COMBOBOX_ADD(vasp_gui.mixer,"4:BROYDEN");
GUI_TOOLTIP(vasp_gui.mixer,"IMIX: Ch. 6.49 DEFAULT: 4\nDetermine the type of mixing.");
	GUI_CHECK_TABLE(table,button,vasp_gui.calc.auto_mixer,mixer_toggle,"AUTO MIXER",1,2,0,1);
GUI_TOOLTIP(button,"Use automatic mixing settings.");
	GUI_ENTRY_TABLE(table,vasp_gui.nelm,vasp_gui.calc.nelm,"%i","NELM=",2,3,0,1);
GUI_TOOLTIP(vasp_gui.nelm,"NELM: Ch. 6.17 DEFAULT: 60\nMaximum number of electronic steps.");
	GUI_ENTRY_TABLE(table,vasp_gui.nelmdl,vasp_gui.calc.nelmdl,"%i","NELMDL=",3,4,0,1);
GUI_TOOLTIP(vasp_gui.nelmdl,"NELMDL: Ch. 6.17 DEFAULT -5(if IALGO=x8)\nSets the number of initial non-selfconsistent steps\ndefault is 0 if ISTART is not 0.");
	GUI_ENTRY_TABLE(table,vasp_gui.nelmin,vasp_gui.calc.nelmin,"%i","NELMIN=",4,5,0,1);
GUI_TOOLTIP(vasp_gui.nelmin,"NELMIN: Ch. 6.17 DEFAULT: 2\nMinimum number of electronic steps.");
/* 2nd line */
	GUI_COMBOBOX_TABLE(table,vasp_gui.mixpre,"MIXPRE=",0,1,1,2);
	GUI_COMBOBOX_ADD(vasp_gui.mixpre,"0:NONE");
	GUI_COMBOBOX_ADD(vasp_gui.mixpre,"1:F=20 INVERSE KERKER");
	GUI_COMBOBOX_ADD(vasp_gui.mixpre,"2:F=200 INVERSE KERKER");
GUI_TOOLTIP(vasp_gui.mixpre,"MIXPRE: Ch. 6.49 DEFAULT: 1\nSelect the preconditioning for Broyden mixing.");
	/*empty: col1*/
	GUI_ENTRY_TABLE(table,vasp_gui.amix,vasp_gui.calc.amix,"%.4lf","AMIX=",2,3,1,2);
GUI_TOOLTIP(vasp_gui.amix,"AMIX: Ch. 6.49 DEFAULT: 0.4\nlinear mixing parameter\nDefault is 0.8 for ultrasoft pseudopotentials.");
	GUI_ENTRY_TABLE(table,vasp_gui.bmix,vasp_gui.calc.bmix,"%.4lf","BMIX=",3,4,1,2);
GUI_TOOLTIP(vasp_gui.bmix,"BMIX: Ch. 6.49 DEFAULT: 1.0\ncutoff wavevector for Kerker mixing.");
	GUI_ENTRY_TABLE(table,vasp_gui.amin,vasp_gui.calc.amin,"%.4lf","AMIN=",4,5,1,2);
GUI_TOOLTIP(vasp_gui.amin,"AMIN: Ch. 6.49 DEFAULT: 0.1\nminimal mixing parameter.");
/* 3rd line */
	GUI_COMBOBOX_TABLE(table,vasp_gui.inimix,"INIMIX=",0,1,2,3);
	GUI_COMBOBOX_ADD(vasp_gui.inimix,"0:LINEAR");
	GUI_COMBOBOX_ADD(vasp_gui.inimix,"1:KERKER");
GUI_TOOLTIP(vasp_gui.inimix,"INIMIX: Ch. 6.49 DEFAULT: 1\ntype of initial mixing in Broyden mixing.");
	GUI_ENTRY_TABLE(table,vasp_gui.maxmix,vasp_gui.calc.maxmix,"%i","MAXMIX=",1,2,2,3);
GUI_TOOLTIP(vasp_gui.maxmix,"MAXMIX: Ch. 6.49 DEFAULT: -45\nMaximum number of steps stored in Broyden mixer.");
	GUI_ENTRY_TABLE(table,vasp_gui.amix_mag,vasp_gui.calc.amix_mag,"%.4lf","AMIX_MAG=",2,3,2,3);
GUI_TOOLTIP(vasp_gui.amix_mag,"AMIX_MAG: Ch. 6.49 DEFAULT: 1.6\nlinear mixing parameter for magnetization.");
	GUI_ENTRY_TABLE(table,vasp_gui.bmix_mag,vasp_gui.calc.bmix_mag,"%.4lf","BMIX_MAG=",3,4,2,3);
GUI_TOOLTIP(vasp_gui.bmix_mag,"BMIX_MAG: Ch. 6.49 DEFAULT: 1.0\ncutoff wavevector for Kerker mixing with magnetization.");
	GUI_ENTRY_TABLE(table,vasp_gui.wc,vasp_gui.calc.wc,"%.1lf","WC=",4,5,2,3);
GUI_TOOLTIP(vasp_gui.wc,"WC: Ch. 6.49 DEFAULT: 1000\nstep weight factor for Broyden mixing.");
/* initialize */
	GUI_COMBOBOX_SETUP(vasp_gui.mixer,3,vasp_mixer_selected);
	GUI_COMBOBOX_SETUP(vasp_gui.mixpre,1,vasp_mixpre_selected);
	GUI_COMBOBOX_SETUP(vasp_gui.inimix,1,vasp_inimix_selected);
	mixer_toggle(NULL,NULL);
/* --- end frame */
/* --- projector */
	GUI_FRAME_NOTE(page,frame,"Projector");
/* create a table in the frame*/
	GUI_TABLE_FRAME(frame,table,5,4);
/* 1st line */
	GUI_COMBOBOX_TABLE(table,vasp_gui.lreal,"LREAL=",0,1,0,1);
	GUI_COMBOBOX_ADD(vasp_gui.lreal,"Auto");
	GUI_COMBOBOX_ADD(vasp_gui.lreal,"On");
	GUI_COMBOBOX_ADD(vasp_gui.lreal,"TRUE");
	GUI_COMBOBOX_ADD(vasp_gui.lreal,"FALSE");
GUI_TOOLTIP(vasp_gui.lreal,"LREAL: Ch. 6.39 DEFAULT: FALSE\nDetermine whether projection operators are evaluated in\nreal-space or reciprocal-space.");
	/*empty: col2*/
	GUI_TEXT_TABLE(table,vasp_gui.ropt,vasp_gui.calc.ropt,"ROPT=",2,5,0,1);
GUI_TOOLTIP(vasp_gui.ropt,"ROPT: Ch. 6.39 DEFAULT: -\nROPT determine the precision of the real-space operator\nfor LREAL=AUTO and LREAL=On.");
/* 2nd line */
	/*empty: col1*/
	GUI_CHECK_TABLE(table,button,vasp_gui.calc.addgrid,NULL,"ADDGRID",1,2,1,2);/*not calling anything*/
GUI_TOOLTIP(button,"ADDGRID: Ch. 6.63 DEFAULT: FALSE\nSelect the third fine grid for augmentation charge.");
	GUI_ENTRY_TABLE(table,vasp_gui.lmaxmix,vasp_gui.calc.lmaxmix,"%i","LMAXMIX=",2,3,1,2);
GUI_TOOLTIP(vasp_gui.lmaxmix,"LMAXMIX: Ch. 6.63 DEFAULT: 2\nMaximum l-quantum number passed to charge density mixer.");
	/*empty: col4*/
	GUI_ENTRY_TABLE(table,vasp_gui.lmaxpaw,vasp_gui.calc.lmaxpaw,"%i","LMAXPAW=",4,5,1,2);
GUI_TOOLTIP(vasp_gui.lmaxpaw,"LMAXPAW: Ch. 6.63 DEFAULT: 2*lmax\nMaximum l-quantum number for the evaluation of the PAW\non-site terms on the radial grid.");
/* 3rd line */
	/*empty: col1*/
	GUI_CHECK_TABLE(table,button,vasp_gui.calc.auto_grid,grid_toggle,"AUTO GRID",1,2,2,3);
GUI_TOOLTIP(button,"Automatically sets NG{X,Y,Z}, and NG{X,Y,Z}F.");
	GUI_ENTRY_TABLE(table,vasp_gui.ngx,vasp_gui.calc.ngx,"%i","NGX=",2,3,2,3);
GUI_TOOLTIP(vasp_gui.ngx,"NGX: Ch. 6.3 DEFAULT: -\nnumber of FFT-mesh grid-points along the x direction.");
	GUI_ENTRY_TABLE(table,vasp_gui.ngy,vasp_gui.calc.ngy,"%i","NGY=",3,4,2,3);
GUI_TOOLTIP(vasp_gui.ngy,"NGY: Ch. 6.3 DEFAULT: -\nnumber of FFT-mesh grid-points along the y direction.");
	GUI_ENTRY_TABLE(table,vasp_gui.ngz,vasp_gui.calc.ngz,"%i","NGZ=",4,5,2,3);
GUI_TOOLTIP(vasp_gui.ngz,"NGZ: Ch. 6.3 DEFAULT: -\nnumber of FFT-mesh grid-points along the z direction.");
/* 4th line */
	/*empty: col1*/
	/*empty: col2*/
	GUI_ENTRY_TABLE(table,vasp_gui.ngxf,vasp_gui.calc.ngxf,"%i","NGXF=",2,3,3,4);
GUI_TOOLTIP(vasp_gui.ngxf,"NGXF: Ch. 6.3 DEFAULT: -\nnumber of fine FFT-mesh grid-points along the x direction.");
	GUI_ENTRY_TABLE(table,vasp_gui.ngyf,vasp_gui.calc.ngyf,"%i","NGYF=",3,4,3,4);
GUI_TOOLTIP(vasp_gui.ngyf,"NGYF: Ch. 6.3 DEFAULT: -\nnumber of fine FFT-mesh grid-points along the y direction.");
	GUI_ENTRY_TABLE(table,vasp_gui.ngzf,vasp_gui.calc.ngzf,"%i","NGZF=",4,5,3,4);
GUI_TOOLTIP(vasp_gui.ngzf,"NGZF: Ch. 6.3 DEFAULT: -\nnumber of fine FFT-mesh grid-points along the z direction.");
/* initialize */
	GUI_COMBOBOX_SETUP(vasp_gui.lreal,3,vasp_lreal_selected);
	grid_toggle(NULL,NULL);
/* --- end frame */

/*-------------------*/
/* page 2 -> ELECT-I */
/*-------------------*/
	GUI_PAGE_NOTE(notebook,page,"ELECT-I");
/* --- smearing */
	GUI_FRAME_NOTE(page,frame,"Smearing");
/* create a table in the frame*/
	GUI_TABLE_FRAME(frame,table,4,2);
/* 1st line */
	GUI_ENTRY_TABLE(table,vasp_gui.ismear,vasp_gui.calc.ismear,"%i","ISMEAR=",0,1,0,1);
GUI_TOOLTIP(vasp_gui.ismear,"ISMEAR: Ch. 6.38 DEFAULT: 1\nDetermine how the partial occupations are set.\nIt is advised to change it to ISMEAR=0\nfor semiconducting or insulating materials.");
	GUI_ENTRY_TABLE(table,vasp_gui.sigma,vasp_gui.calc.sigma,"%.4lf","SIGMA=",1,2,0,1);
GUI_TOOLTIP(vasp_gui.sigma,"SIGMA: Ch. 6.38 DEFAULT: 0.2\nDetermine the width of the smearing (eV)\nFor ISMEAR=0 (semiconductors or insulators)\na small SIGMA=0.05 can be chosen.");
	GUI_CHECK_TABLE(table,vasp_gui.kgamma,vasp_gui.calc.kgamma,NULL,"KGAMMA",2,3,0,1);/*not calling anything*/
GUI_TOOLTIP(vasp_gui.kgamma,"KGAMMA: Ch. 6.4 DEFAULT: TRUE\nWhen KPOINTS file is not present KGAMMA\nindicate if the k-grid is centered around gamma point.");
	GUI_ENTRY_TABLE(table,vasp_gui.kspacing,vasp_gui.calc.kspacing,"%.4lf","KSPACING=",3,4,0,1);
GUI_TOOLTIP(vasp_gui.kspacing,"KSPACING: Ch. 6.4 DEFAULT: 0.5\nWhen KPOINTS file is not present KSPACING\nDetermine the smallest spacing between k-points (Ang^-1).");
/* 2nd line */
	GUI_TEXT_TABLE(table,vasp_gui.fermwe,vasp_gui.calc.fermwe,"FERMWE=",0,2,1,2);
GUI_TOOLTIP(vasp_gui.fermwe,"FERMWE: Ch. 6.38 DEFAULT -\nFor ISMEAR=-2 explicitly sets\noccupancy for each band.");
	GUI_TEXT_TABLE(table,vasp_gui.fermdo,vasp_gui.calc.fermdo,"FERMDO=",2,4,1,2);
GUI_TOOLTIP(vasp_gui.fermdo,"FERDO: Ch. 6.38 DEFAULT -\nFor ISMEAR=-2 and ISPIN=2 sets\noccupancy for down spin of each band.");
/* --- end frame */
/* --- spin */
	GUI_FRAME_NOTE(page,frame,"Spin");
/* create a table in the frame*/
	GUI_TABLE_FRAME(frame,table,5,2);
/* 1st line */
	GUI_CHECK_TABLE(table,button,vasp_gui.calc.ispin,NULL,"SPIN POLARIZED",0,1,0,1);/*not calling anything*/
GUI_TOOLTIP(button,"ISPIN: Ch. 6.12 DEFAULT 1\nSpin polarized calculation sets ISPIN=2.");
	GUI_CHECK_TABLE(table,vasp_gui.lnoncoll,vasp_gui.calc.non_collinear,lnoncoll_toggle,"NON COLLINEAR",1,2,0,1);
GUI_TOOLTIP(vasp_gui.lnoncoll,"LNONCOLLINEAR: Ch. 6.68.1 DEFAULT FALSE\nAllow to perform fully, non-collinear\nmagnetic structure calculations.");
	GUI_CHECK_TABLE(table,vasp_gui.lsorbit,vasp_gui.calc.lsorbit,lsorbit_toggle,"LSORBIT",2,3,0,1);
 GUI_TOOLTIP(vasp_gui.lsorbit,"LSORBIT: Ch. 6.68.2 DEFAULT: FALSE\nSwitch on spin-orbit coupling (for PAW)\nautomatically sets LNONCOLLINEAR when TRUE.");
	GUI_CHECK_TABLE(table,button,vasp_gui.calc.gga_compat,NULL,"GGA_COMPAT",3,4,0,1);/*not calling anything*/
GUI_TOOLTIP(button,"GGA_COMPAT: Ch. 6.42 DEFAULT: TRUE\nRestores the lattice symmetry\nfor gradient corrected functionals.");
	GUI_ENTRY_TABLE(table,vasp_gui.nupdown,vasp_gui.calc.nupdown,"%.1f","NUPDOWN=",4,5,0,1);
GUI_TOOLTIP(vasp_gui.nupdown,"NUPDOWN: Ch. 6.36 DEFAULT -1\nSets the difference between number\nof electrons in up and down spin component.");
/* 2nd line */
	GUI_TEXT_TABLE(table,vasp_gui.magmom,vasp_gui.calc.magmom,"MAGMOM=",0,3,1,2);
GUI_TOOLTIP(vasp_gui.magmom,"MAGMOM: Ch. 6.13 DEFAULT: NIONS\nSets the initial magnetic moment for each atom\nuseful only for calculation with no WAVECAR or CHGCAR\nor magnetic calculations from a non-magnetic start.\nFor non-collinear calculations default value is 3*NIONS.");
	GUI_TEXT_TABLE(table,vasp_gui.saxis,vasp_gui.calc.saxis,"SAXIS=",3,5,1,2);
GUI_TOOLTIP(vasp_gui.saxis,"SAXIS: Ch. 6.68.2 DEFAULT -\nSets the quantisation axis for spin.\nThe default is SAXIS = eps 0 1\nwhere eps is a small quantity.");
/* --- end frame */
/* --- exchange correlation */
	GUI_FRAME_NOTE(page,frame,"Exchange Correlation");
/* create a table in the frame*/
	GUI_TABLE_FRAME(frame,table,5,5);
/* 1st line */
	GUI_COMBOBOX_TABLE(table,vasp_gui.gga,"GGA=",0,1,0,1);
	GUI_COMBOBOX_ADD(vasp_gui.gga,"91:PW91");
	GUI_COMBOBOX_ADD(vasp_gui.gga,"PE:PBE");
	GUI_COMBOBOX_ADD(vasp_gui.gga,"RP:rPBE");
	GUI_COMBOBOX_ADD(vasp_gui.gga,"AM:AM05");
	GUI_COMBOBOX_ADD(vasp_gui.gga,"PS:PBEsol");
GUI_TOOLTIP(vasp_gui.gga,"GGA: Ch. 6.40 DEFAULT -\nSets the gradient corrected functional.\nThe default is selected according to the POTCAR file.");
	GUI_CHECK_TABLE(table,button,vasp_gui.calc.voskown,NULL,"VOSKOWN",1,2,0,1);/*not calling anything*/
GUI_TOOLTIP(button,"VOSKOWN: Ch. 6.41 DEFAULT 0\nSets the Vosko, Wilk and Nusair method\nfor interpolation of correlation functional).\nUseful only if GGA=91 (PW91).");
	/*empty: col2*/
	/*empty: col3*/
	/*empty: col4*/
/* 2nd line */
	GUI_COMBOBOX_TABLE(table,vasp_gui.metagga,"METAGGA=",0,1,1,2);
	GUI_COMBOBOX_ADD(vasp_gui.metagga,"TPSS");
	GUI_COMBOBOX_ADD(vasp_gui.metagga,"rTPSS");
	GUI_COMBOBOX_ADD(vasp_gui.metagga,"M06L");
	GUI_COMBOBOX_ADD(vasp_gui.metagga,"MBJ");
GUI_TOOLTIP(vasp_gui.metagga,"METAGGA: Ch. 6.43 DEFAULT none\nSets the meta-GGA functional.");
	GUI_CHECK_TABLE(table,button,vasp_gui.calc.lmetagga,metagga_toggle,"LMETAGGA",1,2,1,2);
GUI_TOOLTIP(button,"LMETAGGA (secret) DEFAULT FALSE\nIgnored (but recognized and read) by VASP.\nUsed (here) to enable a METAGGA calculation.");
	GUI_CHECK_TABLE(table,button,vasp_gui.calc.lasph,NULL,"LASPH",2,3,1,2);/*not calling anything*/
GUI_TOOLTIP(button,"LASPH: Ch. 6.44 DEFAULT FALSE\nSets calculation of non-spherical contributions\nfrom the gradient corrections in PAW spheres.");
	GUI_CHECK_TABLE(table,button,vasp_gui.calc.lmixtau,NULL,"LMIXTAU",3,4,1,2);/*not calling anything*/
GUI_TOOLTIP(button,"LMIXTAU: Ch. 6.43.2 DEFAULT FALSE\nSets inclusion of kinetic energy density\nto the mixer (useful to converge METAGGA).");
	GUI_ENTRY_TABLE(table,vasp_gui.lmaxtau,vasp_gui.calc.lmaxtau,"%i","LMAXTAU=",4,5,1,2);
GUI_TOOLTIP(vasp_gui.lmaxtau,"LMAXTAU: Ch. 6.43.1 DEFAULT 0\nSets maximum l-quantum number\nto calculate PAW one-center expansion\nof the kinetic energy density.\nDefault is 6 if LASPH is set to TRUE.");
/* 3rd line */
	GUI_TEXT_TABLE(table,vasp_gui.cmbj,vasp_gui.calc.cmbj,"CMBJ=",0,3,2,3);
GUI_TOOLTIP(vasp_gui.cmbj,"CMBJ: Ch 6.43 DEFAULT: -\nFor MBJ meta-GGA calculations CMBJ\nsets the mixing coefficient for Becke-Roussel potential\nfor each species or constant if CMBJ has only one value.");
	GUI_ENTRY_TABLE(table,vasp_gui.cmbja,vasp_gui.calc.cmbja,"%.4lf","CMBJA=",3,4,2,3);
GUI_TOOLTIP(vasp_gui.cmbja,"CMBJA: Ch. 6.43 DEFAULT -0.012\nAlpha parameter to calculate CMBJ\nSelfconsistently at each electronic step.\nUnused if CMBJ is specified.");
	GUI_ENTRY_TABLE(table,vasp_gui.cmbjb,vasp_gui.calc.cmbjb,"%.4lf","CMBJB=",4,5,2,3);
GUI_TOOLTIP(vasp_gui.cmbjb,"CMBJB: Ch. 6.43 DEFAULT 1.023\nBeta parameter to calculate CMBJ\nSelfconsistently at each electronic step.\nUnused if CMBJ is specified.");
/* 4th line */
	GUI_COMBOBOX_TABLE(table,vasp_gui.ldau,"LDAUTYPE=",0,1,3,4);
	GUI_COMBOBOX_ADD(vasp_gui.ldau,"1:LSDA+U Liechtenstein");
	GUI_COMBOBOX_ADD(vasp_gui.ldau,"2:LSDA+U Dudarev");
	GUI_COMBOBOX_ADD(vasp_gui.ldau,"4:LDA+U Liechtenstein");
GUI_TOOLTIP(vasp_gui.ldau,"LDAUTYPE: Ch. 6.70 DEFAULT: 2\nSets the type of L(S)DA+U calculation.");
	GUI_CHECK_TABLE(table,button,vasp_gui.calc.ldau,ldau_toggle,"LDAU",1,2,3,4);
GUI_TOOLTIP(button,"LDAU: Ch. 6.70 DEFAULT: FALSE\nSwitches on L(S)DA+U calculation.");
	GUI_TEXT_TABLE(table,vasp_gui.ldaul,vasp_gui.calc.ldaul,"LDAUL=",2,5,3,4);
GUI_TOOLTIP(vasp_gui.ldaul,"LDAUL: Ch. 6.70 DEFAULT: 2\nSets the l-quantum number (for each species)\nfor which the on-site interaction is added.");
/* 5th line */
	GUI_COMBOBOX_TABLE(table,vasp_gui.ldau_print,"LDAUPRINT=",0,1,4,5);
	GUI_COMBOBOX_ADD(vasp_gui.ldau_print,"0:Silent");
	GUI_COMBOBOX_ADD(vasp_gui.ldau_print,"1:Occupancy matrix");
	GUI_COMBOBOX_ADD(vasp_gui.ldau_print,"2:Full output");
GUI_TOOLTIP(vasp_gui.ldau_print,"LDAUPRINT: Ch. 6.70 DEFAULT 0\nSets the level of verbosity of L(S)DA+U.");
	GUI_TEXT_TABLE(table,vasp_gui.ldauu,vasp_gui.calc.ldauu,"LDAUU=",1,3,4,5);
GUI_TOOLTIP(vasp_gui.ldauu,"LDAUU: Ch. 6.70 DEFAULT: -\nSets effective on-site Coulomb interaction (for each species).");
	GUI_TEXT_TABLE(table,vasp_gui.ldauj,vasp_gui.calc.ldauj,"LDAUJ=",3,5,4,5);
GUI_TOOLTIP(vasp_gui.ldauj,"LDAUJ: Ch. 6.70 DEFAULT: -\nSets effective on-site Exchange interaction (for each species).");
/* initialize */
	GUI_COMBOBOX_SETUP(vasp_gui.gga,1,vasp_gga_selected);
	GUI_COMBOBOX_SETUP(vasp_gui.metagga,3,vasp_metagga_selected);
	GUI_COMBOBOX_SETUP(vasp_gui.ldau,1,vasp_ldau_selected);
	GUI_COMBOBOX_SETUP(vasp_gui.ldau_print,0,vasp_ldau_print_selected);
	metagga_toggle(NULL,NULL);/*initialize TODO: move to end*/
	ldau_toggle(NULL,NULL);/*initialize TODO: move to end*/
/* --- end frame */
/* --- dipole */
	GUI_FRAME_NOTE(page,frame,"Dipole");
/* create a table in the frame*/
	GUI_TABLE_FRAME(frame,table,2,5);
/* 1st line */
	GUI_COMBOBOX_TABLE(table,vasp_gui.idipol,"IDIPOL=",0,1,0,1);
	GUI_COMBOBOX_ADD(vasp_gui.idipol,"0:no calcul");
	GUI_COMBOBOX_ADD(vasp_gui.idipol,"1:u axis");
	GUI_COMBOBOX_ADD(vasp_gui.idipol,"2:v axis");
	GUI_COMBOBOX_ADD(vasp_gui.idipol,"3:w axis");
	GUI_COMBOBOX_ADD(vasp_gui.idipol,"4:all axis");
GUI_TOOLTIP(vasp_gui.idipol,"IDIPOL: Ch. 6.64 DEFAULT: -\nSets the direction of the calculated dipole moment\nparallel to the 1st (1), 2nd (2) or 3rd (3) lattice vector.\nIDIPOL=4 trigger full dipole calculation\n0 switches off dipole calculation (remove IDIPOL).");
	GUI_CHECK_TABLE(table,vasp_gui.ldipol,vasp_gui.calc.ldipol,NULL,"LDIPOL",1,2,0,1);/*not calling anything*/
GUI_TOOLTIP(vasp_gui.ldipol,"LDIPOL: Ch. 6.64 DEFAULT: FALSE\nSets potential dipole correction of the\nslab/molecule periodic-boundary induced errors.");
	GUI_CHECK_TABLE(table,vasp_gui.lmono,vasp_gui.calc.lmono,NULL,"LMONO",2,3,0,1);/*not calling anything*/
GUI_TOOLTIP(vasp_gui.lmono,"LMONO: Ch. 6.64 DEFAULT: FALSE\nSets potential monopole correction of the\nslab/molecule periodic-boundary induced errors.");
	GUI_TEXT_TABLE(table,vasp_gui.dipol,vasp_gui.calc.dipol,"DIPOL=",3,5,0,1);
GUI_TOOLTIP(vasp_gui.dipol,"DIPOL: Ch. 6.64 DEFAULT: -\nSets the center of the net charge distribution.\nIf left blank a guess value is sets at runtime.");
/* 2nd line */
	/*empty: col0*/
	/*empty: col1*/
	/*empty: col2*/
	GUI_ENTRY_TABLE(table,vasp_gui.epsilon,vasp_gui.calc.epsilon,"%.4lf","EPSILON=",3,4,2,3);
GUI_TOOLTIP(vasp_gui.epsilon,"EPSILON: Ch. 6.64 DEFAULT: 1\nSets the dielectric constant of the medium.");
	GUI_ENTRY_TABLE(table,vasp_gui.efield,vasp_gui.calc.efield,"%.4lf","EFIELD=",4,5,2,3);
GUI_TOOLTIP(vasp_gui.efield,"EFIELD: Ch. 6.64 DEFAULT:  0\nApply an external electrostatic field\nin a slab, or molecule calculation.\nOnly a single component can be given\nparallel to IDIPOL direction (1-3).");
/*initialize sensitive widget*/
	GUI_COMBOBOX_SETUP(vasp_gui.idipol,0,vasp_idipol_selected);
if(vasp_gui.calc.idipol!=VID_0) {
	GUI_UNLOCK(vasp_gui.ldipol);
	GUI_UNLOCK(vasp_gui.lmono);
	GUI_UNLOCK(vasp_gui.dipol);
	GUI_UNLOCK(vasp_gui.epsilon);
	if(vasp_gui.calc.idipol!=VID_4) GUI_UNLOCK(vasp_gui.efield);
}else{
	GUI_LOCK(vasp_gui.ldipol);
	GUI_LOCK(vasp_gui.lmono);
	GUI_LOCK(vasp_gui.dipol);
	GUI_LOCK(vasp_gui.epsilon);
	GUI_LOCK(vasp_gui.efield);
}
/* --- end frame */

/*---------------------*/
/* page 2b -> ELECT-II */
/*---------------------*/
	GUI_PAGE_NOTE(notebook,page,"ELECT-II");
/* --- DOS */
	GUI_FRAME_NOTE(page,frame,"Density of States");
/* create a table in the frame*/
	GUI_TABLE_FRAME(frame,table,2,5);
/* 1st line */
	GUI_COMBOBOX_TABLE(table,vasp_gui.lorbit,"LORBIT=",0,1,0,1);
	GUI_COMBOBOX_ADD(vasp_gui.lorbit,"0:PROCAR");
	GUI_COMBOBOX_ADD(vasp_gui.lorbit,"1:lm-PROCAR");
	GUI_COMBOBOX_ADD(vasp_gui.lorbit,"2:phase+lm-PROCAR");
	GUI_COMBOBOX_ADD(vasp_gui.lorbit,"5:PROOUT");
	GUI_COMBOBOX_ADD(vasp_gui.lorbit,"10:PROCAR");
	GUI_COMBOBOX_ADD(vasp_gui.lorbit,"11:lm-PROCAR");
	GUI_COMBOBOX_ADD(vasp_gui.lorbit,"12:phase+lm-PROCAR");
GUI_TOOLTIP(vasp_gui.lorbit,"LORBIT: Ch. 6.34 DEFAULT: 0\nSets the spd- and site projected\nwavefunction character of each band\nand the local partial DOS.");
	GUI_ENTRY_TABLE(table,vasp_gui.nedos,vasp_gui.calc.nedos,"%i","NEDOS=",1,2,0,1);
GUI_TOOLTIP(vasp_gui.nedos,"NEDOS: Ch. 6.37 DEFAULT: 301\nNumber of points for which DOS is calculated.");
	GUI_ENTRY_TABLE(table,vasp_gui.emin,vasp_gui.calc.emin,"%.2lf","EMIN=",2,3,0,1);
GUI_TOOLTIP(vasp_gui.emin,"EMIN: Ch. 6.37 DEFAULT: low eig-eps\nSets the smallest energy for which DOS is calculated.\nDefault is lowest Kohn-Sham eigenvalue minus a small quantity.");
	GUI_ENTRY_TABLE(table,vasp_gui.emax,vasp_gui.calc.emax,"%.2lf","EMAX=",3,4,0,1);
GUI_TOOLTIP(vasp_gui.emax,"EMAX: Ch. 6.37 DEFAULT: high eig+eps\nSets the highest energy for which DOS is calculated.\nDefault is highest Kohn-Sham eigenvalue plus a small quantity.");
	GUI_ENTRY_TABLE(table,vasp_gui.efermi,vasp_gui.calc.efermi,"%.2lf","EFERMI=",4,5,0,1);
GUI_TOOLTIP(vasp_gui.efermi,"EFERMI: (secret) DEFAULT: EREF\nSets initial Fermi energy before it is calculated.\nUsually there is no use to change this setting.\nEREF is actually another \"secret\" parameter.");
/* 2nd line */
	GUI_LEFT_LABEL_TABLE(table,"(LORBIT>5) => Req. PAW",0,1,1,2);
	GUI_CHECK_TABLE(table,vasp_gui.have_paw,vasp_gui.calc.have_paw,NULL,"HAVE PAW",1,2,1,2);/*not calling anything*/
GUI_TOOLTIP(vasp_gui.have_paw,"Is set if a PAW POTCAR file is detected.");
	GUI_LOCK(vasp_gui.have_paw);/*just an indication*/
	GUI_TEXT_TABLE(table,vasp_gui.rwigs,vasp_gui.calc.rwigs,"RWIGS=",2,5,1,2);
GUI_TOOLTIP(vasp_gui.rwigs,"RWIGS: Ch. 6.33 DEFAULT: -\nSets the Wigner Seitz radius for each species.\nDefault value is read from POTCAR.");
/* initialize */
if(vasp_gui.calc.have_paw) GUI_COMBOBOX_SETUP(vasp_gui.lorbit,4,vasp_lorbit_selected);
else GUI_COMBOBOX_SETUP(vasp_gui.lorbit,0,vasp_lorbit_selected);
/* --- end frame */
/* --- Linear Response */
	GUI_FRAME_NOTE(page,frame,"Linear Response");
/* create a table in the frame*/
	GUI_TABLE_FRAME(frame,table,1,5);
/* 1st line */
	GUI_CHECK_TABLE(table,vasp_gui.loptics,vasp_gui.calc.loptics,NULL,"LOPTICS",0,1,0,1);/*not calling anything*/
GUI_TOOLTIP(vasp_gui.loptics,"LOPTICS: Ch. 6.72.1 DEFAULT: FALSE\nSwitch frequency dependent dielectric matrix calculation\nover a grid determined by NEDOS.\nIt is advised to increase NEDOS & NBANDS.");
	GUI_CHECK_TABLE(table,vasp_gui.lepsilon,vasp_gui.calc.lepsilon,NULL,"LEPSILON",1,2,0,1);/*not calling anything*/
GUI_TOOLTIP(vasp_gui.lepsilon,"LEPSILON: Ch. 6.72.4 DEFAULT: FALSE\nSwitch dielectric matrix calculation using\ndensity functional perturbation theory.\nConcurrent to LOPTICS, it is not recommended\nto run both at the same time.");
	GUI_CHECK_TABLE(table,vasp_gui.lrpa,vasp_gui.calc.lrpa,NULL,"LRPA",2,3,0,1);/*not calling anything*/
GUI_TOOLTIP(vasp_gui.lrpa,"LRPA: Ch. 6.72.5 DEFAULT: FALSE\nSwitch local field effects calculation\non the Hartree level only (no XC).");
	GUI_CHECK_TABLE(table,vasp_gui.lnabla,vasp_gui.calc.lnabla,NULL,"LNABLA",3,4,0,1);/*not calling anything*/
GUI_TOOLTIP(vasp_gui.lnabla,"LNABLA: Ch. 6.72.3 DEFAULT: FALSE\nSwitch to the simple transversal expressions\nof the frequency dependent dielectric matrix.\nUsually there is no use to change this setting.");
	GUI_ENTRY_TABLE(table,vasp_gui.cshift,vasp_gui.calc.cshift,"%.2lf","CSHIFT=",4,5,0,1);
GUI_TOOLTIP(vasp_gui.cshift,"CSHIFT: Ch. 6.72.2 DEFAULT: 0.1\nSets the complex shift \%nu of the Kramers-Kronig\ntransformation of the dielectric function.\nIf CSHIFT is decreased, NEDOS should be increased.");
/* --- end frame */


/*-----------------*/
/* page 3 -> IONIC */
/*-----------------*/
	GUI_PAGE_NOTE(notebook,page,"IONIC");
/* --- general */
	GUI_FRAME_NOTE(page,frame,"General");
/* create a table in the frame*/
	GUI_TABLE_FRAME(frame,table,2,5);
/* 1st line */
	GUI_COMBOBOX_TABLE(table,vasp_gui.ibrion,"IBRION=",0,1,0,1);
	GUI_COMBOBOX_ADD(vasp_gui.ibrion,"-1:Nothing");
	GUI_COMBOBOX_ADD(vasp_gui.ibrion,"0:Molecular Dynamics");
	GUI_COMBOBOX_ADD(vasp_gui.ibrion,"1:Quasi-Newton");
	GUI_COMBOBOX_ADD(vasp_gui.ibrion,"2:Conjugate Gradients");
	GUI_COMBOBOX_ADD(vasp_gui.ibrion,"3:Damped");
	GUI_COMBOBOX_ADD(vasp_gui.ibrion,"5:Finite Differences");
	GUI_COMBOBOX_ADD(vasp_gui.ibrion,"6:Finite Diff. with SYM");
	GUI_COMBOBOX_ADD(vasp_gui.ibrion,"7:Perturbation Theory");
	GUI_COMBOBOX_ADD(vasp_gui.ibrion,"8:Pert. Theory with SYM");
	GUI_COMBOBOX_ADD(vasp_gui.ibrion,"44:Transition State");
GUI_TOOLTIP(vasp_gui.ibrion,"IBRION: Ch. 6.22 DEFAULT -1\nSets the algorithm for the ions motion.\nDefault is IBRION=0 if NSW>1.");
	GUI_ENTRY_TABLE(table,vasp_gui.nsw,vasp_gui.calc.nsw,"%i","NSW=",1,2,0,1);
GUI_TOOLTIP(vasp_gui.nsw,"NSW: Ch. 6.20 DEFAULT 0\nSet the maximum number of ionic steps.");
	GUI_ENTRY_TABLE(table,vasp_gui.ediffg,vasp_gui.calc.ediffg,"%.2E","EDIFFG=",2,3,0,1);
GUI_TOOLTIP(vasp_gui.ediffg,"EDIFFG: Ch. 6.19 DEFAULT EDIFF*10\nSets the energy criterion for ending ionic relaxation\n>0 value are for energy difference (eV)\n<0 values are for absolute forces value (eV/Ang).");
	GUI_ENTRY_TABLE(table,vasp_gui.potim,vasp_gui.calc.potim,"%.4lf","POTIM=",3,4,0,1);
GUI_TOOLTIP(vasp_gui.potim,"POTIM: Ch. 6.23 DEFAULT 0.5\nFor IBRION=0 (MD) sets time step (fs)\nfor IBRION=1-3 sets force scaling constant.");
	GUI_ENTRY_TABLE(table,vasp_gui.pstress,vasp_gui.calc.pstress,"%.4lf","PSTRESS=",4,5,0,1);
GUI_TOOLTIP(vasp_gui.pstress,"PSTRESS: Ch. 6.25 DEFAULT 0\nSets the Pulay stress correction (kBar)\nalso used to set an external pressure.");
/* 2nd line */
	GUI_COMBOBOX_TABLE(table,vasp_gui.isif,"ISIF=",0,1,1,2);
	GUI_COMBOBOX_ADD(vasp_gui.isif,"0:F_0_I_0_0");
	GUI_COMBOBOX_ADD(vasp_gui.isif,"1:F_P_I_0_0");
	GUI_COMBOBOX_ADD(vasp_gui.isif,"2:F_S_I_0_0");
	GUI_COMBOBOX_ADD(vasp_gui.isif,"3:F_S_I_S_V");
	GUI_COMBOBOX_ADD(vasp_gui.isif,"4:F_S_I_S_0");
	GUI_COMBOBOX_ADD(vasp_gui.isif,"5:F_S_0_S_0");
	GUI_COMBOBOX_ADD(vasp_gui.isif,"6:F_S_0_S_V");
	GUI_COMBOBOX_ADD(vasp_gui.isif,"7:F_S_0_0_V");
GUI_TOOLTIP(vasp_gui.isif,"ISIF: Ch. 6.24 DEFAULT 2\nSwitch to calculate F(force) and S(stress)\nand to relax I(ion) S(cell shape) and V(cell voume)\n0 = no calcul/no relaxation P = only total Pressure\nDefault is 0(F_0_I_0_0) if IBRION=0.");
	GUI_CHECK_TABLE(table,vasp_gui.relax_ions,vasp_gui.rions,isif_toggle,"atom positions",1,2,1,2);
GUI_TOOLTIP(vasp_gui.relax_ions,"Switches atomic relaxation.\nUpdates ISIF accordingly (but not IBRION).");
	GUI_CHECK_TABLE(table,vasp_gui.relax_shape,vasp_gui.rshape,isif_toggle,"cell shape",2,3,1,2);
GUI_TOOLTIP(vasp_gui.relax_ions,"Switches cell shape relaxation.\nUpdates ISIF accordingly (but not IBRION).");
	GUI_CHECK_TABLE(table,vasp_gui.relax_volume,vasp_gui.rvolume,isif_toggle,"cell volume",3,4,1,2);
GUI_TOOLTIP(vasp_gui.relax_ions,"Switches cell volume relaxation.\nUpdates ISIF accordingly (but not IBRION).");
	GUI_ENTRY_TABLE(table,vasp_gui.nfree,vasp_gui.calc.nfree,"%i","NFREE=",4,5,1,2);
GUI_TOOLTIP(vasp_gui.nfree,"NFREE: Ch. 6.22 DEFAULT: -\nSets the number of ionic steps kept in history\nfor IBRION=1 Quasi-Newton algorithm.\nFor finite difference IBRION=5-6 NFREE\nsets the number of displacement/direction/atom.");
/* initialize */
	GUI_COMBOBOX_SETUP(vasp_gui.ibrion,0,vasp_ibrion_selected);
	GUI_COMBOBOX_SETUP(vasp_gui.isif,0,vasp_isif_selected);
/* --- end frame */
/* --- MD */
	GUI_FRAME_NOTE(page,frame,"Molecular Dynamics");
/* create a table in the frame*/
	GUI_TABLE_FRAME(frame,table,2,5);
/* 1st line */
	//multicolumn label
	GUI_LABEL_TABLE(table,"SET IBRION=0 FOR MD",0,1,0,2);
	GUI_ENTRY_TABLE(table,vasp_gui.tebeg,vasp_gui.calc.tebeg,"%.1lf","TEBEG=",1,2,0,1);
GUI_TOOLTIP(vasp_gui.tebeg,"TEBEG: Ch. 6.29 DEFAULT: 0\nSets the start temperature of MD run.\nTo be corrected by (Nion-1)/Nion.");
	GUI_ENTRY_TABLE(table,vasp_gui.teend,vasp_gui.calc.teend,"%.1lf","TEEND=",2,3,0,1);
GUI_TOOLTIP(vasp_gui.teend,"TEEND: Ch. 6.29 DEFAULT: TEBEG\nSets the end temperature of MD run.\nTo be corrected by (Nion-1)/Nion.");
	GUI_ENTRY_TABLE(table,vasp_gui.smass,vasp_gui.calc.smass,"%.1lf","SMASS=",3,4,0,1);
GUI_TOOLTIP(vasp_gui.smass,"SMASS: Ch. 6.30 DEFAULT 0\nControls the velocities during MD.\nFor IBRION=3 SMASS is used as a damping factor.");
	/*empty: col4*/
/* 2nd line */
	/*occ: col0*/
	GUI_ENTRY_TABLE(table,vasp_gui.nblock,vasp_gui.calc.nblock,"%i","NBLOCK=",1,2,1,2);
GUI_TOOLTIP(vasp_gui.nblock,"NBLOCK: Ch. 6.21 DEFAULT: 1\nSets number of steps after which DOS/pair correlation\nfunction are calculated, and XDATCAR updated.\nFor SMASS=-1 kinetic energy is scaled every NBLOCK.");
	GUI_ENTRY_TABLE(table,vasp_gui.kblock,vasp_gui.calc.kblock,"%i","KBLOCK=",2,3,1,2);
GUI_TOOLTIP(vasp_gui.kblock,"KBLOCK: Ch. 6.21 DEFAULT: NSW\nAfter KBLOCK*NBLOCK iterations DOS and\naverage pair correlation function is written\nto DOSCAR and PCDAT, respectively.");
	GUI_ENTRY_TABLE(table,vasp_gui.npaco,vasp_gui.calc.npaco,"%i","NPACO=",3,4,1,2);
GUI_TOOLTIP(vasp_gui.npaco,"NPACO Ch. 6.31 DEFAULT: 256\nNumber of slots of pair correlation function.");
	GUI_ENTRY_TABLE(table,vasp_gui.apaco,vasp_gui.calc.apaco,"%.4lf","APACO=",4,5,1,2);
GUI_TOOLTIP(vasp_gui.apaco,"APACO: Ch. 6.31 DEFAULT: 16\nMax distance (Ang) for calculation\nof pair correlation function.");
/* initialize sensitive widget*/	
	if (vasp_gui.calc.ibrion==0) {/*selecting MD enable molecular dynamics settings*/
		GUI_UNLOCK(vasp_gui.tebeg);
		GUI_UNLOCK(vasp_gui.teend);
		GUI_UNLOCK(vasp_gui.smass);
		GUI_UNLOCK(vasp_gui.nblock);
		GUI_UNLOCK(vasp_gui.kblock);
		GUI_UNLOCK(vasp_gui.npaco);
		GUI_UNLOCK(vasp_gui.apaco);
	}else{
		GUI_LOCK(vasp_gui.tebeg);
		GUI_LOCK(vasp_gui.teend);
		GUI_LOCK(vasp_gui.smass);
		GUI_LOCK(vasp_gui.nblock);
		GUI_LOCK(vasp_gui.kblock);
		GUI_LOCK(vasp_gui.npaco);
		GUI_LOCK(vasp_gui.apaco);
	}
/* --- end frame */
/* --- POSCAR */
	GUI_FRAME_NOTE(page,frame,"POSCAR & SYMMTERY");
/* create a table in the frame*/
	GUI_TABLE_FRAME(frame,table,7,5);
/* 1st line */
	GUI_COMBOBOX_TABLE(table,vasp_gui.poscar_free,"DYN:",0,1,0,1);
	GUI_COMBOBOX_ADD(vasp_gui.poscar_free,"All atom FIXED");
	GUI_COMBOBOX_ADD(vasp_gui.poscar_free,"All atom FREE");
	GUI_COMBOBOX_ADD(vasp_gui.poscar_free,"Manual selection");
GUI_TOOLTIP(vasp_gui.poscar_free,"Set position tags to \"F F F\" (FIXED)\n\"T T T\" (FREE) or user modifications.");
	GUI_CHECK_TABLE(table,button,vasp_gui.calc.poscar_sd,toggle_poscar_sd,"Sel. Dynamics",1,2,0,1);
GUI_TOOLTIP(button,"Allow use of tags for ion motion (Ch. 5.7).");
	GUI_ENTRY_TABLE(table,vasp_gui.isym,vasp_gui.calc.isym,"%i","ISYM=",2,3,0,1);
GUI_TOOLTIP(vasp_gui.isym,"ISYM: Ch. 6.27 DEFAULT 2\nSwitch symmetry method (0=OFF)\nDefault is 1 for ultrasoft, 2 for PAW PP.");
	GUI_ENTRY_TABLE(table,vasp_gui.sym_prec,vasp_gui.calc.sym_prec,"%.2lE","_PREC=",3,4,0,1);
GUI_TOOLTIP(vasp_gui.sym_prec,"SYMPREC: Ch. 6.27 DEFAULT 10^-5\nSets the accuracy of ion positions\nin determining the system symmetry.");
	GUI_ENTRY_TABLE(table,vasp_gui.poscar_a0,vasp_gui.calc.poscar_a0,"%.4lf","A0=",4,5,0,1);
GUI_TOOLTIP(vasp_gui.poscar_a0,"Lattice parameter (Ch. 5.7)\nSets the lattice constant (Ang.) parameter\nevery lattice vector will be scaled with.");
/* 2nd line */
	/*empty: col0 (TODO: SUPERCELL)*/
	GUI_CHECK_TABLE(table,vasp_gui.poscar_direct,vasp_gui.calc.poscar_direct,toggle_poscar_direct,"Direct Coord.",1,2,1,2);
GUI_TOOLTIP(vasp_gui.poscar_direct,"Sets whether ions positions are given\nin direct lattice or cartesian coordinate.");
	GUI_ENTRY_TABLE(table,vasp_gui.poscar_ux,vasp_gui.calc.poscar_ux,"%.4lf","UX=",2,3,1,2);
GUI_TOOLTIP(vasp_gui.poscar_ux,"First lattice vector x component.");
	GUI_ENTRY_TABLE(table,vasp_gui.poscar_uy,vasp_gui.calc.poscar_uy,"%.4lf","UY=",3,4,1,2);
GUI_TOOLTIP(vasp_gui.poscar_uy,"First lattice vector y component.");
	GUI_ENTRY_TABLE(table,vasp_gui.poscar_uz,vasp_gui.calc.poscar_uz,"%.4lf","UZ=",4,5,1,2);
GUI_TOOLTIP(vasp_gui.poscar_uz,"First lattice vector z component.");
/* 3rd line */
	/*empty: col0 (TODO: -X)*/
	/*empty: col1 (TODO: +X)*/
	GUI_ENTRY_TABLE(table,vasp_gui.poscar_vx,vasp_gui.calc.poscar_vx,"%.4lf","VX=",2,3,2,3);
GUI_TOOLTIP(vasp_gui.poscar_vx,"Second lattice vector x component.");
	GUI_ENTRY_TABLE(table,vasp_gui.poscar_vy,vasp_gui.calc.poscar_vy,"%.4lf","VY=",3,4,2,3);
GUI_TOOLTIP(vasp_gui.poscar_vy,"Second lattice vector y component.");
	GUI_ENTRY_TABLE(table,vasp_gui.poscar_vz,vasp_gui.calc.poscar_vz,"%.4lf","VZ=",4,5,2,3);
GUI_TOOLTIP(vasp_gui.poscar_vz,"Second lattice vector z component.");
/* 4th line */
	/*empty: col0 (TODO: -Y)*/
	/*empty: col1 (TODO: +Y)*/
	GUI_ENTRY_TABLE(table,vasp_gui.poscar_wx,vasp_gui.calc.poscar_wx,"%.4lf","WX=",2,3,3,4);
GUI_TOOLTIP(vasp_gui.poscar_wx,"Third lattice vector x component.");
	GUI_ENTRY_TABLE(table,vasp_gui.poscar_wy,vasp_gui.calc.poscar_wy,"%.4lf","WY=",3,4,3,4);
GUI_TOOLTIP(vasp_gui.poscar_wy,"Third lattice vector y component.");
	GUI_ENTRY_TABLE(table,vasp_gui.poscar_wz,vasp_gui.calc.poscar_wz,"%.4lf","WZ=",4,5,3,4);
GUI_TOOLTIP(vasp_gui.poscar_wz,"Third lattice vector z component.");
/* 5th line */
	/*empty: col0 (TODO: -Z)*/
	/*empty: col1 (TODO: +Z)*/
	GUI_CHECK_TABLE(table,vasp_gui.poscar_tx,vasp_gui.calc.poscar_tx,NULL,"TX",2,3,4,5);/*not calling anything*/
GUI_TOOLTIP(vasp_gui.poscar_tx,"Tag for the motion of current ion\nin first lattice coordinate.");
	GUI_CHECK_TABLE(table,vasp_gui.poscar_ty,vasp_gui.calc.poscar_ty,NULL,"TY",3,4,4,5);/*not calling anything*/
GUI_TOOLTIP(vasp_gui.poscar_tx,"Tag for the motion of current ion\nin second lattice coordinate.");
	GUI_CHECK_TABLE(table,vasp_gui.poscar_tz,vasp_gui.calc.poscar_tz,NULL,"TZ",4,5,4,5);/*not calling anything*/
GUI_TOOLTIP(vasp_gui.poscar_tx,"Tag for the motion of current ion\nin third lattice coordinate.");
/* 6th line */
	/*NEW: rewrite of ATOM list combobox*/
	GUI_COMBOBOX_TABLE(table,vasp_gui.poscar_atoms,"ATOMS:",0,4,5,6);
	idx=0;
	for (list2=data->cores ; list2 ; list2=g_slist_next(list2)){
		gchar sym[3];
		core=list2->data;
		sym[0]=core->atom_label[0];
		if(g_ascii_islower(core->atom_label[1])) sym[1]=core->atom_label[1];
		else sym[1]='\0';
		sym[2]='\0';
		if(vasp_gui.calc.poscar_free==VPF_FIXED){
			GUI_COMBOBOX_ADD(vasp_gui.poscar_atoms,
				g_strdup_printf("%.8lf %.8lf %.8lf  F   F   F ! atom: %i (%s)",core->x[0],core->x[1],core->x[2],idx,sym));
		}else{
			GUI_COMBOBOX_ADD(vasp_gui.poscar_atoms,
				g_strdup_printf("%.8lf %.8lf %.8lf  T   T   T ! atom: %i (%s)",core->x[0],core->x[1],core->x[2],idx,sym));
		}
		idx++;
	}
	vasp_gui.calc.atoms_total=idx;
	GUI_COMBOBOX_ADD(vasp_gui.poscar_atoms,"ADD atom");
	GUI_COMBOBOX_SETUP(vasp_gui.poscar_atoms,idx,vasp_poscar_atoms_selected);
        /*this combo is special and should be explicitely defined*/
GUI_TOOLTIP(vasp_gui.poscar_atoms,"Atoms positions in POSCAR format.");
	/*2 buttons*/
	GUI_2BUTTONS_TABLE(table,vasp_atom_modified,vasp_atom_delete,4,5,5,6);
/* 7th line */
	GUI_ENTRY_TABLE(table,vasp_gui.poscar_index,vasp_gui.calc.poscar_index,"%i","INDEX:",0,1,6,7);
	GUI_LOCK(vasp_gui.poscar_index);/*<- is changed only by selecting poscar_atoms*/
GUI_TOOLTIP(vasp_gui.poscar_index,"Atom number (in POSCAR ORDER).");
	GUI_ENTRY_TABLE(table,vasp_gui.poscar_symbol,vasp_gui.calc.poscar_symbol,"%s","SYMBOL:",1,2,6,7);
GUI_TOOLTIP(vasp_gui.poscar_index,"Atom symbol (will correspond to POSCAR).\nChanging symbol will change index accordingly.");
	GUI_ENTRY_TABLE(table,vasp_gui.poscar_x,vasp_gui.calc.poscar_x,"%.6lf","X=",2,3,6,7);
GUI_TOOLTIP(vasp_gui.poscar_x,"atom component X: x in cartesian mode\nand first lattice coordinate in direct mode.");
	GUI_ENTRY_TABLE(table,vasp_gui.poscar_y,vasp_gui.calc.poscar_y,"%.6lf","Y=",3,4,6,7);
GUI_TOOLTIP(vasp_gui.poscar_y,"atom component Y: y in cartesian mode\nand second lattice coordinate in direct mode.");
	GUI_ENTRY_TABLE(table,vasp_gui.poscar_z,vasp_gui.calc.poscar_z,"%.6lf","Z=",4,5,6,7);
GUI_TOOLTIP(vasp_gui.poscar_z,"atom component Z: z in cartesian mode\nand third lattice coordinate in direct mode.");
/* initialize sensitive widget*/
	GUI_COMBOBOX_SETUP(vasp_gui.poscar_free,1,vasp_poscar_free_selected);
if((vasp_gui.calc.poscar_free==VPF_FIXED)||(vasp_gui.calc.poscar_free==VPF_FREE)){
	GUI_LOCK(vasp_gui.poscar_tx);
	GUI_LOCK(vasp_gui.poscar_ty);
	GUI_LOCK(vasp_gui.poscar_tz);
	if(vasp_gui.calc.poscar_free==VPF_FIXED){
		vasp_gui.calc.poscar_tx=FALSE;
		vasp_gui.calc.poscar_ty=FALSE;
		vasp_gui.calc.poscar_tz=FALSE;
	}else{
		vasp_gui.calc.poscar_tx=TRUE;
		vasp_gui.calc.poscar_tx=TRUE;
		vasp_gui.calc.poscar_tx=TRUE;
	}
}
/* (re)initialize */
	toggle_poscar_sd(NULL,NULL);
	vasp_poscar_free_selected(vasp_gui.poscar_free,NULL);
/* --- end frame */

/*-------------------*/
/* page 4 -> KPOINTS */
/*-------------------*/
	GUI_PAGE_NOTE(notebook,page,"KPOINTS");
/* --- from INCAR */
	GUI_FRAME_NOTE(page,frame,"from INCAR");
/* create a table in the frame*/
	GUI_TABLE_FRAME(frame,table,1,5);
/* 1st line */
	GUI_ENTRY_TABLE(table,vasp_gui.ismear_3,vasp_gui.calc.ismear,"%i","ISMEAR=",0,1,0,1);
GUI_TOOLTIP(vasp_gui.ismear_3,"ISMEAR: Ch. 6.38 DEFAULT: 1\nDetermine how the partial occupations are set.\nIt is advised to change it to ISMEAR=0\nfor semiconducting or insulating materials.");
	GUI_CHECK_TABLE(table,vasp_gui.kgamma_3,vasp_gui.calc.kgamma,NULL,"KGAMMA",1,2,0,1);/*not calling anything*/
GUI_TOOLTIP(vasp_gui.kgamma_3,"KGAMMA: Ch. 6.4 DEFAULT: TRUE\nWhen KPOINTS file is not present KGAMMA\nindicate if the k-grid is centered around gamma point.");
	GUI_ENTRY_TABLE(table,vasp_gui.kspacing_3,vasp_gui.calc.kspacing,"%.4f","KSPACING=",2,3,0,1);
GUI_TOOLTIP(vasp_gui.kspacing_3,"KSPACING: Ch. 6.4 DEFAULT: 0.5\nWhen KPOINTS file is not present KSPACING\nDetermine the smallest spacing between k-points (Ang^-1).");
	/*empty: col3*/
	/*empty: col4*/
/* --- end frame */
/* --- General */
	GUI_FRAME_NOTE(page,frame,"General");
/* create a table in the frame*/
	GUI_TABLE_FRAME(frame,table,2,5);
/* 1st line */
	GUI_COMBOBOX_TABLE(table,vasp_gui.kpoints_mode,"GEN:",0,1,0,1);
	GUI_COMBOBOX_ADD(vasp_gui.kpoints_mode,"Manual Entry");
	GUI_COMBOBOX_ADD(vasp_gui.kpoints_mode,"Line Mode");
	GUI_COMBOBOX_ADD(vasp_gui.kpoints_mode,"Automatic (M&P)");
	GUI_COMBOBOX_ADD(vasp_gui.kpoints_mode,"Gamma (M&P)");
	GUI_COMBOBOX_ADD(vasp_gui.kpoints_mode,"Monkhorst-Pack classic");
	GUI_COMBOBOX_ADD(vasp_gui.kpoints_mode,"Basis set definition");
GUI_TOOLTIP(vasp_gui.kpoints_mode,"Sets how the k-points are generated/entered as defined in Ch. 5.5.");
	GUI_CHECK_TABLE(table,vasp_gui.kpoints_cart,vasp_gui.calc.kpoints_cart,NULL,"Cartesian",1,2,0,1);/*not calling anything*/
GUI_TOOLTIP(vasp_gui.kpoints_cart,"Switch cartesian/reciprocal lattice coordinate mode.");
	GUI_SPIN_TABLE(table,vasp_gui.kpoints_kx,vasp_gui.calc.kpoints_kx,NULL,"KX",2,3,0,1);
GUI_TOOLTIP(vasp_gui.kpoints_kx,"Grid size for automated generation\n1st component only for M&P.");
	GUI_SPIN_TABLE(table,vasp_gui.kpoints_ky,vasp_gui.calc.kpoints_ky,NULL,"KY",3,4,0,1);
GUI_TOOLTIP(vasp_gui.kpoints_ky,"2nd Component of the M&P grid size.");
	GUI_SPIN_TABLE(table,vasp_gui.kpoints_kz,vasp_gui.calc.kpoints_kz,NULL,"KZ",4,5,0,1);
GUI_TOOLTIP(vasp_gui.kpoints_kz,"3rd Component of the M&P grid size.");
/* 2nd line */
	GUI_ENTRY_TABLE(table,vasp_gui.kpoints_nkpts,vasp_gui.calc.kpoints_nkpts,"%i","NKPTS=",0,1,1,2);
GUI_TOOLTIP(vasp_gui.kpoints_nkpts,"Total number of k-points (0 for automated methods).");
	GUI_CHECK_TABLE(table,vasp_gui.have_tetra,vasp_gui.calc.kpoints_tetra,toogle_tetra,"Tetrahedron",1,2,1,2);
GUI_TOOLTIP(vasp_gui.have_tetra,"Switch on the tetrahedron method for entering k-points.\nISMEAR must be set to -5 for this method!");
	GUI_ENTRY_TABLE(table,vasp_gui.kpoints_sx,vasp_gui.calc.kpoints_sx,"%.4lf","SX=",2,3,1,2);
GUI_TOOLTIP(vasp_gui.kpoints_sx,"In case a grid is used to generate k-points\nThis will set a shift from grid origin\nin 1st reciprocal lattice component.");
	GUI_ENTRY_TABLE(table,vasp_gui.kpoints_sy,vasp_gui.calc.kpoints_sy,"%.4lf","SY=",3,4,1,2);
GUI_TOOLTIP(vasp_gui.kpoints_sy,"In case a grid is used to generate k-points\nThis will set a shift from grid origin\nin 2nd reciprocal lattice component.");
	GUI_ENTRY_TABLE(table,vasp_gui.kpoints_sz,vasp_gui.calc.kpoints_sz,"%.4lf","SZ=",4,5,1,2);
GUI_TOOLTIP(vasp_gui.kpoints_sz,"In case a grid is used to generate k-points\nThis will set a shift from grid origin\nin 3rd reciprocal lattice component.");
/* --- end frame */
/* --- Coordinate */
	GUI_FRAME_NOTE(page,vasp_gui.kpoints_coord,"Coordinate");
/* create a table in the frame*/
	GUI_TABLE_FRAME(vasp_gui.kpoints_coord,table,2,5);
/* 1st line */
	GUI_COMBOBOX_TABLE(table,vasp_gui.kpoints_kpts,"KPOINT:",0,4,0,1);
	GUI_COMBOBOX_ADD(vasp_gui.kpoints_kpts,"ADD kpoint");
GUI_TOOLTIP(vasp_gui.kpoints_kpts,"List of the k-points coordinates.\nThe index after comment sign \"! kpoint:\" can be used\nfor entering tetrahedron edge definition.");
        /*2 buttons*/
	GUI_2BUTTONS_TABLE(table,vasp_kpoint_modified,vasp_kpoint_delete,4,5,0,1);
/* 2nd line */
	GUI_ENTRY_TABLE(table,vasp_gui.kpoints_i,vasp_gui.calc.kpoints_i,"%i","INDEX:",0,1,1,2);
	GUI_LOCK(vasp_gui.kpoints_i);/*<- is changed only by selecting kpoints_kpts*/
GUI_TOOLTIP(vasp_gui.kpoints_i,"The index of current k-point.\n(can be changed by deletion/addition only)");
	GUI_ENTRY_TABLE(table,vasp_gui.kpoints_x,vasp_gui.calc.kpoints_x,"%.4lf","X=",1,2,1,2);
GUI_TOOLTIP(vasp_gui.kpoints_x,"Coordinate of k-point in 1st component of selected mode.\n(cartesian or reciprocal)");
	GUI_ENTRY_TABLE(table,vasp_gui.kpoints_y,vasp_gui.calc.kpoints_y,"%.4lf","Y=",2,3,1,2);
GUI_TOOLTIP(vasp_gui.kpoints_y,"Coordinate of k-point in 2nd component of selected mode.\n(cartesian or reciprocal)");
	GUI_ENTRY_TABLE(table,vasp_gui.kpoints_z,vasp_gui.calc.kpoints_z,"%.4lf","Z=",3,4,1,2);
GUI_TOOLTIP(vasp_gui.kpoints_z,"Coordinate of k-point in 3rd component of selected mode.\n(cartesian or reciprocal)");
	GUI_ENTRY_TABLE(table,vasp_gui.kpoints_w,vasp_gui.calc.kpoints_w,"%.4lf","W=",4,5,1,2);
GUI_TOOLTIP(vasp_gui.kpoints_w,"Symmetry degeneration weight of k-point.\nSum of all weight does not have to be 1\nbut relative ratio should be respected.");
/* --- end frame */
/* --- Tetrahedron */
	GUI_FRAME_NOTE(page,vasp_gui.kpoints_tetra,"Tetrahedron");
/* create a table in the frame*/
	GUI_TABLE_FRAME(vasp_gui.kpoints_tetra,table,3,6);
/* 1st line */
	GUI_ENTRY_TABLE(table,vasp_gui.tetra_total,vasp_gui.calc.tetra_total,"%i","TOTAL=",0,1,0,1);
GUI_TOOLTIP(vasp_gui.tetra_total,"Total number of tetrahedrons.");
	GUI_LABEL_TABLE(table,"EXACT VOLUME:",1,2,0,1);
	GUI_ENTRY_TABLE(table,vasp_gui.tetra_volume,vasp_gui.calc.tetra_volume,"%.6lf","V=",2,3,0,1);
GUI_TOOLTIP(vasp_gui.tetra_volume,"Volume weight for a single tetrahedron.\nAll have a volume defined by the ratio of\ntetrahedron volume / Brillouin zone Volume.");
	GUI_LABEL_TABLE(table,"for (A,B,C,D) tetrahedron",3,6,0,1);
/* 2nd line */
	GUI_COMBOBOX_TABLE(table,vasp_gui.tetra,"TETRAHEDRON:",0,5,1,2);
	GUI_COMBOBOX_ADD(vasp_gui.tetra,"ADD tetrahedron");
GUI_TOOLTIP(vasp_gui.tetra,"Tetrahedron list, consisting of\nsymmetry degeneration weight followed by\nfour corner points of each tetrahedron\nin k-point indices given above.");
	/*2 buttons*/
	GUI_2BUTTONS_TABLE(table,vasp_tetra_modified,vasp_tetra_delete,5,6,1,2);
/* 3rd line */
	GUI_ENTRY_TABLE(table,vasp_gui.tetra_i,vasp_gui.calc.tetra_i,"%i","INDEX:",0,1,2,3);
	GUI_LOCK(vasp_gui.tetra_i);/*<- is changed only by selecting tetra*/
GUI_TOOLTIP(vasp_gui.tetra_i,"The index of current tetrahedron.\n(can be changed by deletion/addition only)");
	GUI_ENTRY_TABLE(table,vasp_gui.tetra_w,vasp_gui.calc.tetra_w,"%.4lf","Degen. W=",1,2,2,3);
GUI_TOOLTIP(vasp_gui.tetra_w,"Symmetry degeneration weight.\nUnlike k-points, weight must be normalized\nie. total Brillouin zone Volume should be recovered\nby summing all degeneration * V");
	GUI_ENTRY_TABLE(table,vasp_gui.tetra_a,vasp_gui.calc.tetra_a,"%i","PtA:",2,3,2,3);
GUI_TOOLTIP(vasp_gui.tetra_a,"1st corner of the tetrahedron (A) in k-point index.");
	GUI_ENTRY_TABLE(table,vasp_gui.tetra_b,vasp_gui.calc.tetra_b,"%i","PtB:",3,4,2,3);
GUI_TOOLTIP(vasp_gui.tetra_b,"2nd corner of the tetrahedron (B) in k-point index.");
	GUI_ENTRY_TABLE(table,vasp_gui.tetra_c,vasp_gui.calc.tetra_c,"%i","PtC:",4,5,2,3);
GUI_TOOLTIP(vasp_gui.tetra_c,"3rd corner of the tetrahedron (C) in k-point index.");
	GUI_ENTRY_TABLE(table,vasp_gui.tetra_d,vasp_gui.calc.tetra_d,"%i","PtD:",5,6,2,3);
GUI_TOOLTIP(vasp_gui.tetra_d,"4th corner of the tetrahedron (D) in k-point index.");
/* initialize */
	GUI_ENTRY_CHANGE(vasp_gui.ismear_3,ismear_changed,data);
	GUI_ENTRY_CHANGE(vasp_gui.kspacing_3,kspacing_changed,data);
	GUI_COMBOBOX_SETUP(vasp_gui.kpoints_mode,vasp_gui.calc.kpoints_mode,vasp_kpoints_mode_selected);
	GUI_COMBOBOX_SETUP(vasp_gui.kpoints_kpts,0,vasp_kpoints_kpts_selected);
	GUI_COMBOBOX_SETUP(vasp_gui.tetra,0,vasp_tetra_selected);
	toogle_tetra();
	ismear_changed(vasp_gui.ismear_3,data);/*necessary*/
/* --- end frame */
/* TODO: import KPOINTS/IBZKPT file*/

/*------------------*/
/* page 5 -> POTCAR */
/*------------------*/
	GUI_PAGE_NOTE(notebook,page,"POTCAR");
/* --- general */
	GUI_FRAME_NOTE(page,frame,"General");
/* create a table in the frame*/
	GUI_TABLE_FRAME(frame,table,2,1);
/* 1st line */
	GUI_TEXT_TABLE(table,vasp_gui.poscar_species,vasp_gui.calc.species_symbols,"SPECIES:",0,1,0,1);
	GUI_LOCK(vasp_gui.poscar_species);
GUI_TOOLTIP(vasp_gui.poscar_species,"List of atomic symbol of the species\ndetected in the POSCAR ionic information.");
/* 2nd line */
	GUI_LABEL_TABLE(table,"(each species will require a separate POTCAR information)",0,1,1,2);
/* --- end frame */
/* --- get POTCAR information */
	GUI_FRAME_NOTE(page,frame,"get POTCAR information");
/* create a table in the frame*/
	GUI_TABLE_FRAME(frame,table,3,4);
/* 1st line */
	/* 2 lines radiobutton */
	GUI_2RADIO_TABLE(table,
		vasp_gui.potcar_select_file,"Use POTCAR file",toggle_potcar_file,
		vasp_gui.potcar_select_folder,"Use POTCAR path",toggle_potcar_folder,0,1,0,2);
GUI_TOOLTIP(vasp_gui.potcar_select_file,"Load a previously prepared POTCAR file.\nThe species (and species order) should match\nEXACTLY those of the POSCAR ion information.");
GUI_TOOLTIP(vasp_gui.potcar_select_folder,"Use the directory where all POTCAR information are available.\nUsually a directory containing atomic symbols sub-directories.");
	GUI_TEXT_TABLE(table,vasp_gui.potcar_file,vasp_gui.calc.potcar_file,"FILE:",1,3,0,1);
	GUI_OPEN_BUTTON_TABLE(table,vasp_gui.potcar_file_button,load_potcar_file_dialog,3,4,0,1);
GUI_TOOLTIP(vasp_gui.potcar_file,"POTCAR file, including its full path.");
/* 2nd line */
	GUI_TEXT_TABLE(table,vasp_gui.potcar_folder,vasp_gui.calc.potcar_folder,"PATH:",1,3,1,2);
	GUI_OPEN_BUTTON_TABLE(table,vasp_gui.potcar_folder_button,load_potcar_folder_dialog,3,4,1,2);
GUI_TOOLTIP(vasp_gui.potcar_folder,"Full PATH to the POTCAR sub-directories.");
/* 3rd line */
	GUI_COMBOBOX_TABLE(table,vasp_gui.species_flavor,"FLAVOR:",2,3,2,3);
	GUI_OPEN_BUTTON_TABLE(table,vasp_gui.species_button,apply_species_flavor,3,4,2,3);
GUI_TOOLTIP(vasp_gui.species_flavor,"When Using PATH method, for each species a choice\nbetween different pseudopotential setting (flavors) exists\nEx: _h hard, _s soft _pv include pseudo-core valence, etc.");
/* initialize */
	toggle_potcar_file(NULL,NULL);
/* --- end frame */
/* --- POTCAR results */
	GUI_FRAME_NOTE(page,frame,"POTCAR results");
/* create a table in the frame*/
	GUI_TABLE_FRAME(frame,table,2,2);
/* 1st line */
	//multicolumn label
	GUI_LABEL_TABLE(table,"DETECTED",0,1,0,2);
	GUI_TEXT_TABLE(table,vasp_gui.potcar_species,vasp_gui.calc.potcar_species,"SPECIES:",1,2,0,1);
GUI_TOOLTIP(vasp_gui.potcar_species,"List of atomic symbol of the species\ndetected from the POTCAR setting.");
	GUI_LOCK(vasp_gui.potcar_species);
/* 2nd line */
	/*occ: col0*/
	GUI_TEXT_TABLE(table,vasp_gui.potcar_species_flavor,vasp_gui.calc.potcar_species_flavor,"FLAVOR:",1,2,1,2);
GUI_TOOLTIP(vasp_gui.potcar_species_flavor,"List of the pseudopotential flavors (PATH method)\ndetected from the POTCAR setting.");
	GUI_LOCK(vasp_gui.potcar_species_flavor);
/* --- end frame */

/*----------------*/
/* page 6 -> EXEC */
/*----------------*/
	GUI_PAGE_NOTE(notebook,page,"EXEC");
/* --- general */
	GUI_FRAME_NOTE(page,frame,"General");
/* create a table in the frame*/
	GUI_TABLE_FRAME(frame,table,3,6);
/* 1st line */
	GUI_LABEL_TABLE(table,"VASP EXEC:",0,1,0,1);
	GUI_TEXT_TABLE(table,vasp_gui.job_vasp_exe,vasp_gui.calc.job_vasp_exe,"FILE=",1,5,0,1);
GUI_TOOLTIP(vasp_gui.job_vasp_exe,"Location (full path) of the vasp executable.");
	GUI_OPEN_BUTTON_TABLE(table,button,load_vasp_exe_dialog,5,6,0,1);
/* 2nd line */
	GUI_LABEL_TABLE(table,"CALCUL PATH:",0,1,1,2);
	GUI_TEXT_TABLE(table,vasp_gui.job_path,vasp_gui.calc.job_path,"PATH=",1,5,1,2);
GUI_TOOLTIP(vasp_gui.job_path,"Location (full path) of the calculation directory.\nit is IMPORTANT to check this parameter so that\n* no other calculation is overwritten\n* calculation will not start in an unexpected directory.");
	GUI_OPEN_BUTTON_TABLE(table,button,load_path_dialog,5,6,1,2);
/* 3rd line */
	GUI_LABEL_TABLE(table,"OUTPUT:",0,1,2,3);
	GUI_CHECK_TABLE(table,button,vasp_gui.calc.lwave,NULL,"LWAVE",1,2,2,3);/*not calling anything*/
GUI_TOOLTIP(button,"LWAVE: Ch. 6.52 DEFAULT: TRUE\nIf set the WAVECAR (wavefunctions) file is written.");
	GUI_CHECK_TABLE(table,button,vasp_gui.calc.lcharg,NULL,"LCHARG",2,3,2,3);/*not calling anything*/
GUI_TOOLTIP(button,"LCHARG: Ch. 6.52 DEFAULT: TRUE\nIf set the CHGCAR (charge density) file is written.");
	GUI_CHECK_TABLE(table,button,vasp_gui.calc.lvtot,NULL,"LVTOT",3,4,2,3);/*not calling anything*/
GUI_TOOLTIP(button,"LVTOT: Ch. 6.53 DEFAULT: FALSE\nIf set the LOCPOT (local potential) file is written.");
	GUI_CHECK_TABLE(table,button,vasp_gui.calc.lvhar,NULL,"LVHAR",4,5,2,3);/*not calling anything*/
GUI_TOOLTIP(button,"LVHAR: Ch. 6.54 DEFAULT: FALSE\nIf set the full local potential is written to LOCPOT\n(ie. ionic+Hartree+XC) otherwise only electrostatic\n(ie. ionic+Hartree) is written.");
	GUI_CHECK_TABLE(table,button,vasp_gui.calc.lelf,NULL,"LELF",5,6,2,3);/*not calling anything*/
GUI_TOOLTIP(button,"LELF: Ch 6.55 DEFAULT: FALSE\nIf set the ELFCAR (electron localization function)\nfile is written.");
/* --- end frame */
/* --- PARALLEL */
	GUI_FRAME_NOTE(page,frame,"Parallel Optimization");
/* create a table in the frame*/
	GUI_TABLE_FRAME(frame,table,2,6);
/* 1st line */
	GUI_LABEL_TABLE(table,"MPIRUN:",0,1,0,1);
	GUI_TEXT_TABLE(table,vasp_gui.job_mpirun,vasp_gui.calc.job_mpirun,"FILE=",1,5,0,1);
GUI_TOOLTIP(vasp_gui.job_mpirun,"Location (full path) of mpirun executable.\nRequired for parallel calculation only.");
	GUI_OPEN_BUTTON_TABLE(table,button,load_mpirun_dialog,5,6,0,1);
/* 2nd line */
	GUI_SPIN_TABLE(table,vasp_gui.job_nproc,vasp_gui.calc.job_nproc,parallel_eq,"NP=",0,1,1,2);
GUI_TOOLTIP(vasp_gui.job_nproc,"Total number of CPU used for calculation.\nwill be used as N in \"mpirun -np N vasp\"\nto start vasp calculation.");
	GUI_SPIN_TABLE(table,vasp_gui.ncore,vasp_gui.calc.ncore,parallel_eq,"NCORE=",1,2,1,2);
GUI_TOOLTIP(vasp_gui.ncore,"NCORE: Ch. 6.57 DEFAULT: 1\nSets how many cores works on one orbital.\nNote that NCORE is limited by KPAR.");
	GUI_SPIN_TABLE(table,vasp_gui.kpar,vasp_gui.calc.kpar,parallel_eq,"KPAR=",2,3,1,2);
GUI_TOOLTIP(vasp_gui.kpar,"KPAR: Ch. 6.57 DEFAULT: 1\nSets the number of k-points treated in parallel.\nEach k-point is worked on by NCORE CPUs.");
	GUI_CHECK_TABLE(table,button,vasp_gui.calc.lplane,NULL,"LPLANE",3,4,1,2);/*not calling anything*/
GUI_TOOLTIP(button,"LPLANE: Ch. 6.57 DEFAULT: TRUE\nIf set the real-space data distribution is done plane wise.");
	GUI_CHECK_TABLE(table,button,vasp_gui.calc.lscalu,NULL,"LSCALU",4,5,1,2);/*not calling anything*/
GUI_TOOLTIP(button,"LSCALU: Ch. 6.59 DEFAULT: FALSE\nIf set parallel LU decomposition will be used.");
	GUI_CHECK_TABLE(table,button,vasp_gui.calc.lscalapack,NULL,"LSCALAPACK",5,6,1,2);/*not calling anything*/
GUI_TOOLTIP(button,"LSCALAPACK: Ch. 6.59 DEFAULT: FALSE\nIf set scaLAPACK routines will be used.");
/* --- end frame */
/* --- DISTANT */
	GUI_FRAME_NOTE(page,frame,"Distant calculation");
/*create a vbox in frame*/
	GUI_VBOX_FRAME(frame,vbox);
/* 1st line */
	GUI_LINE_BOX(vbox,hbox);
/* 2nd line */
	GUI_LINE_BOX(vbox,hbox);
	GUI_LABEL_BOX(hbox,label,"Remote execution of VASP is currently UNDER CONSTRUCTION.");
/* 3rd line */
	GUI_LINE_BOX(vbox,hbox);
	GUI_LABEL_BOX(hbox,label,"For now, please save files and manually transfer them to the distant server, then submit and retreive files as usual.");
/* 4th line */
	GUI_LINE_BOX(vbox,hbox);
	GUI_LABEL_BOX(hbox,label,"The future GDIS job interface will hopefully automate the whole process...");
/* 5th line */
	GUI_LINE_BOX(vbox,hbox);
	GUI_LABEL_BOX(hbox,label,"...and will be made available in the VASP calculation interface at that time.");
	GUI_NEW_SEPARATOR(hbox,separator);
	GUI_LABEL_BOX(hbox,label,"--OVHPA");
/* 6th line */
	GUI_LINE_BOX(vbox,hbox);
/* --- end frame */

/* --- Outside of notebook */
	GUI_FRAME_WINDOW(vasp_gui.window,frame);
	GUI_VBOX_FRAME(frame,vbox);
/* Action buttons */
	GUI_SAVE_ACTION(vasp_gui.window,vasp_gui.button_save,save_vasp_calc,NULL);
	GUI_EXEC_ACTION(vasp_gui.window,vasp_gui.button_exec,exec_calc,NULL);
	GUI_CLOSE_ACTION(vasp_gui.window,button,quit_vasp_gui,dialog);
/* connect to signals */
	GUI_PAGE_CHANGE(notebook,vasp_gui_page_switch,NULL);
/* all done */
	vasp_gui_refresh();/*refresh once more*/
	GUI_SHOW(vasp_gui.window);/*display*/
        sysenv.refresh_dialog=TRUE;
}

