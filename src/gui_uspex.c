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

/* simple USPEX launcher */

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
#include "file_uspex.h"
#include "gui_defs.h"
#include "gui_uspex.h"

extern struct sysenv_pak sysenv;

/*globals variables*/
struct uspex_calc_gui uspex_gui;

/*************************/
/* initialize gui values */
/*************************/
void gui_uspex_init(){
	uspex_gui.cur_page=USPEX_PAGE_SET_I;
	uspex_gui.is_dirty=TRUE;
}
/************************************************************/
/* Load calculation setting from an included Parameters.txt */
/************************************************************/
void load_parameters_dialog(void){
	GUI_OBJ *file_chooser;
	gint have_answer;
	gchar *filename;
	gchar *text;
	/**/
	GUI_PREPARE_OPEN_DIALOG(uspex_gui.window,file_chooser,"Select a Parameters.txt File","*.txt","Parameters.txt files");
	GUI_OPEN_DIALOG_RUN(file_chooser,have_answer,filename);
	if(have_answer){
		/*there should be no case where have_answer==TRUE AND filename==NULL, but just in case.*/
		if(filename) {
			uspex_calc_struct *test_calc;
			text=g_strdup_printf("%s",filename);
			GUI_ENTRY_TEXT(uspex_gui.file_entry,text);
			g_free(text);
			test_calc = read_uspex_parameters(filename);
			if(test_calc != NULL) {
				/*unable to open Parameters.txt*/
			}else{
				copy_uspex_parameters(test_calc,&(uspex_gui.calc));
				uspex_gui_refresh();/*refresh GUI with new uspex_gui.calc*/
			}
			g_free (filename);
			free_uspex_parameters(test_calc);
			g_free(test_calc);
		}
	}
	GUI_KILL_OPEN_DIALOG(file_chooser);
}
/****************************/
/* load the uspex executable */
/****************************/
void load_uspex_exe_dialog(void){
	GUI_OBJ *file_chooser;
	gint have_answer;
	gchar *filename;
	gchar *text;
	/**/
	GUI_PREPARE_OPEN_DIALOG(uspex_gui.window,file_chooser,"Select USPEX Executable","uspex","uspex_exec");
	GUI_OPEN_DIALOG_RUN(file_chooser,have_answer,filename);
	if(have_answer){
		/*there should be no case where have_answer==TRUE AND filename==NULL, but just in case.*/
		if(filename) {
			text=g_strdup_printf("%s",filename);
			GUI_ENTRY_TEXT(uspex_gui.job_uspex_exe,text);
			g_free(text);
			USPEX_REG_TEXT(job_uspex_exe);
			g_free (filename);
		}
	}
	GUI_KILL_OPEN_DIALOG(file_chooser);
}
/*****************************************/
/* point to the path to run the USPEX job */
/*****************************************/
void uspex_path_dialog(void){
	GUI_OBJ *file_chooser;
	gint have_answer;
	gchar *filename;
	gchar *text;
	/**/
	GUI_PREPARE_OPEN_FOLDER(uspex_gui.window,file_chooser,"Select working folder");
	GUI_OPEN_DIALOG_RUN(file_chooser,have_answer,filename);
	if(have_answer){
		/*there should be no case where have_answer==TRUE AND filename==NULL, but just in case.*/
		if(filename) {
			text=g_strdup_printf("%s",filename);
			GUI_ENTRY_TEXT(uspex_gui.job_path,text);
			g_free(text);
			USPEX_REG_TEXT(job_path);
			g_free (filename);
		}
	}
	GUI_KILL_OPEN_DIALOG(file_chooser);
}
/**************************************************/
/* refresh the GUI with value from uspex_gui.calc */
/**************************************************/
void uspex_gui_refresh(){
//	gchar *text;
	/*refresh the whole GUI based on uspex_gui.calc information*/
}
/*******************************/
/* selecting calculationMethod */
/*******************************/
void uspex_method_selected(GUI_OBJ *w){
	gint index;
	GUI_COMBOBOX_GET(w,index);
	switch (index){
	case 0://USPEX
		uspex_gui.calc.calculationMethod = US_CM_USPEX;break;
	case 1://META
		uspex_gui.calc.calculationMethod = US_CM_META;break;
	case 2://VCNEB
		uspex_gui.calc.calculationMethod = US_CM_VCNEB;break;
	case 3://PSO
		uspex_gui.calc.calculationMethod = US_CM_PSO;break;
	case 4://UNKNOWN
	default:
		uspex_gui.calc.calculationMethod = US_CM_UNKNOWN;
	}
}
/******************************************************/
/* Transfer gdouble spin value into proper gint: _dim */
/******************************************************/
void spin_update_dim(){
	uspex_gui.calc._calctype_dim=(gint)uspex_gui._dim;
}
/*****************************************/
/* toggle mol -> update calculation type */
/*****************************************/
void mol_toggle(void){
	gchar *text;
	if(uspex_gui.calc._calctype_mol){
		if(uspex_gui.calc._calctype_var) uspex_gui.calc.calculationType = uspex_gui.calc._calctype_dim*100+10+1;
		else uspex_gui.calc.calculationType = uspex_gui.calc._calctype_dim*100+10;
	}else{
		if(uspex_gui.calc._calctype_var) uspex_gui.calc.calculationType = uspex_gui.calc._calctype_dim*100+1;
		else uspex_gui.calc.calculationType = uspex_gui.calc._calctype_dim*100;
	}
	text=g_strdup_printf("%i",uspex_gui.calc.calculationType);
	GUI_ENTRY_TEXT(uspex_gui.calculationType,text);
	g_free(text);
}
/*****************************************/
/* toggle var -> update calculation type */
/*****************************************/
void var_toggle(void){
	gchar *text;
	if(uspex_gui.calc._calctype_var){
		if(uspex_gui.calc._calctype_mol) uspex_gui.calc.calculationType = uspex_gui.calc._calctype_dim*100+10+1;
		else uspex_gui.calc.calculationType = uspex_gui.calc._calctype_dim*100+1;
	}else{
		if(uspex_gui.calc._calctype_mol) uspex_gui.calc.calculationType = uspex_gui.calc._calctype_dim*100+10;
		else uspex_gui.calc.calculationType = uspex_gui.calc._calctype_dim*100;
	}
	text=g_strdup_printf("%i",uspex_gui.calc.calculationType);
	GUI_ENTRY_TEXT(uspex_gui.calculationType,text);
	g_free(text);
}
/***********************************/
/* selecting optimization property */
/***********************************/
void uspex_optimization_selected(GUI_OBJ *w){
	gint index;
	GUI_COMBOBOX_GET(w,index);
	switch (index){
	case 0://1
		uspex_gui.calc.optType=US_OT_ENTHALPY;break;
	case 1://2
		uspex_gui.calc.optType=US_OT_VOLUME;break;
	case 2://3
		uspex_gui.calc.optType=US_OT_HARDNESS;break;
	case 3://4
		uspex_gui.calc.optType=US_OT_ORDER;break;
	case 4://5
		uspex_gui.calc.optType=US_OT_DISTANCE;break;
	case 5://6
		uspex_gui.calc.optType=US_OT_DIELEC_S;break;
	case 6://7
		uspex_gui.calc.optType=US_OT_GAP;break;
	case 7://8
		uspex_gui.calc.optType=US_OT_DIELEC_GAP;break;
	case 8://9
		uspex_gui.calc.optType=US_OT_MAG;break;
	case 9://10
		uspex_gui.calc.optType=US_OT_QE;break;
	case 10://1101
		uspex_gui.calc.optType=US_OT_BULK_M;break;
	case 11://1102
		uspex_gui.calc.optType=US_OT_SHEAR_M;break;
	case 12://1103
		uspex_gui.calc.optType=US_OT_YOUNG_M;break;
	case 13://1104
		uspex_gui.calc.optType=US_OT_POISSON;break;
	case 14://1105
		uspex_gui.calc.optType=US_OT_PUGH_R;break;
	case 15://1106
		uspex_gui.calc.optType=US_OT_VICKERS_H;break;
	case 16://1107
		uspex_gui.calc.optType=US_OT_FRACTURE;break;
	case 17://1108
		uspex_gui.calc.optType=US_OT_DEBYE_T;break;
	case 18://1109
		uspex_gui.calc.optType=US_OT_SOUND_V;break;
	case 19://1110
		uspex_gui.calc.optType=US_OT_SWAVE_V;break;
	case 20://1111
		uspex_gui.calc.optType=US_OT_PWAVE_V;break;
	default:
		uspex_gui.calc.optType=US_OT_UNKNOWN;
	}
}

/************************/
/* Switch notebook page */
/************************/
void uspex_gui_page_switch(GUI_NOTE *notebook,GUI_OBJ *page,guint page_num){
	gint from_page;
	gchar *text;
	/* some (few) values need to be updated from a page to another*/
	GUI_NOTE_PAGE_NUMBER(notebook,from_page);
	if(from_page==(gint)page_num) return;/*not moved*/
	uspex_gui.cur_page=page_num;
	GUI_LOCK(page);
	if(page_num==USPEX_PAGE_SET_I){
		/**/
	} else if(page_num==USPEX_PAGE_SET_II){
		/**/
	} else if (page_num==USPEX_PAGE_SPECIFICS){
		/**/
	} else if (page_num==USPEX_PAGE_ABINITIO){
		/**/
	} else if (page_num==USPEX_PAGE_EXEC){
		/**/
	}
	GUI_UNLOCK(page);
	/*to be continued...*/
}
/****************************************************/
/* convert current uspex_gui into uspex_calc_struct */
/****************************************************/
void uspex_gui_sync(){
	/* sync start here */
}
/***************************/
/* save current parameters */
/***************************/
gint save_uspex_calc(){
#define DEBUG_USPEX_SAVE 0
	FILE *save_fp;
	gchar *filename;
	/*1-synchronize*/
	uspex_gui_sync();
	/*2-output*/
	//SAVE calculation parameters to INPUT.txt
	filename=g_strdup_printf("%s/INPUT.txt",uspex_gui.calc.job_path);
	if(dump_uspex_parameters(filename,&(uspex_gui.calc))<0){
		fprintf(stderr,"#ERR: can't open file INPUT.txt for WRITE!\n");
		return -1;
	}
	g_free(filename);
	/*debug print*/
#if DEBUG_USPEX_SAVE
	dump_uspex_parameters(stdout,&(uspex_gui.calc));
#endif
	return 0;
}
/*********************************/
/* Cleanup calculation structure */
/*********************************/
void uspex_cleanup(){
	/*we don't have to free anything.. everything will be hopefully gone after the dialog is closed*/
}
/*****************************************/
/* Execute or enqueue a uspex calculation */
/*****************************************/
void run_uspex_exec(uspex_exec_struct *uspex_exec){
	/* Execute a uspex task TODO distant job */
	gchar *cmd;
	gchar *cwd;/*for push/pull*/

#if __WIN32
/*	At present running uspex in __WIN32 environment is impossible.
 *	However,  it should be possible to launch a (remote) job on a 
 *	distant server. See TODO */
	fprintf(stderr,"USPEX calculation can't be done in this environment.\n");return;
#else 
	/*direct launch*/
	cmd = g_strdup_printf("%s > uspex.log",(*uspex_exec).job_uspex_exe);
#endif
	cwd=sysenv.cwd;/*push*/
	sysenv.cwd=g_strdup_printf("%s",(*uspex_exec).job_path);
	task_sync(cmd);
	g_free(sysenv.cwd);
	sysenv.cwd=cwd;/*pull*/
	g_free(cmd);
	(*uspex_exec).have_result=TRUE;
}
/*****************************/
/* cleanup task: load result */
/*****************************/
void cleanup_uspex_exec(uspex_exec_struct *uspex_exec){
	/*USPEX process has exit, try to load the result*/
	struct model_pak *result_model;
	gchar *filename;
	int index;
	FILE *vf;
	/*sync_ wait for result?*/
	while(!(*uspex_exec).have_result) usleep(500*1000);/*sleep 500ms until job is done*/
	/*init a model*/
	result_model=model_new();
	model_init(result_model);
	/*put result into it*/
	/*1- find the last result folder*/
	index=1;
	filename=g_strdup_printf("%s/results1/OUTPUT.txt",(*uspex_exec).job_path);
	vf=fopen(filename,"r");
	while(vf!=NULL){
		g_free(filename);
		fclose(vf);
		index++;
		filename=g_strdup_printf("%s/results%i/OUTPUT.txt",(*uspex_exec).job_path,index);
		vf=fopen(filename,"r");
	}
	index--;
	g_free(filename);
	/*this is the last openable folder/OUTPUT.txt*/
	filename=g_strdup_printf("%s/results%i/OUTPUT.txt",(*uspex_exec).job_path,index);
	(*uspex_exec).index=index;
	/*TODO: detect if a calculation have failed*/
	file_load(filename,result_model);/*TODO: load result without annoying tree_select_active*/
	model_prep(result_model);
	tree_model_add(result_model);
	tree_model_refresh(result_model);
	canvas_shuffle();
	redraw_canvas(ALL);/*ALL: necessary?*/
/*just wipe the structure*/
	sysenv.uspex_calc_list=g_slist_remove(sysenv.uspex_calc_list,uspex_exec);/*does not free, does it?*/
	g_free(uspex_exec);
}
/******************************/
/* Enqueue a uspex calculation */
/******************************/
void uspex_exec_calc(){
	uspex_exec_struct *uspex_exec;
	/*this will sync then enqueue a USPEX calculation*/
	if(save_uspex_calc()) return;/*sync and save all file*/
/*copy structure to the list*/
	uspex_exec=g_malloc(sizeof(uspex_exec_struct));
	uspex_exec->job_id=g_slist_length(sysenv.uspex_calc_list);
	uspex_exec->have_result=FALSE;
	uspex_exec->job_uspex_exe=g_strdup_printf("%s",uspex_gui.calc.job_uspex_exe);
	uspex_exec->job_path=g_strdup_printf("%s",uspex_gui.calc.job_path);
	/*prepend to calc list*/
	sysenv.uspex_calc_list = g_slist_prepend (sysenv.uspex_calc_list,uspex_exec);
	/*launch uspex in a task*/
	GUI_LOCK(uspex_gui.button_save);
	GUI_LOCK(uspex_gui.button_exec);
	task_new("USPEX", &run_uspex_exec,uspex_exec,&cleanup_uspex_exec,uspex_exec,sysenv.active_model);
	/*when task is launched, close dialog*/
	uspex_cleanup();
	GUI_CLOSE(uspex_gui.window);
}
/********************************/
/* Quit uspex calculation dialog */
/********************************/
void quit_uspex_gui(GUI_OBJ *w, gpointer data){
	struct model_pak *model=data;
	uspex_cleanup();
	dialog_destroy(w,model);
}
/********************************/
/* USPEX calculation main dialog */
/********************************/
void gui_uspex_dialog(void){
/*launch the main interface*/
/*TODO: use distant connections*/
	gchar *title;
	gpointer dialog;
	GUI_OBJ *frame, *vbox, *hbox, *table;
	GUI_OBJ *notebook, *page, *button, *label;
	GUI_OBJ *separator;
	/* special */
	struct model_pak *data;
	struct core_pak *core;
	gint idx;
	gchar *tmp;
/* checks */
	data = sysenv.active_model;
	if (!data) return;
/* do we have a uspexrun.xml model?  */
	if (data->id == USPEX) {
		uspex_output_struct *uspex_output=data->uspex;
		uspex_calc_struct *uspex_calc;
		/*take GUI default from the opened USPEX output*/
		free_uspex_parameters(&(uspex_gui.calc));
		uspex_calc = read_uspex_parameters(_UO.calc->filename);
		if(uspex_calc){
			uspex_gui.have_output=TRUE;
			copy_uspex_parameters(uspex_calc,&(uspex_gui.calc));
			gui_uspex_init();
		}else{
			uspex_gui.have_output=FALSE;
			init_uspex_parameters(&(uspex_gui.calc));
			gui_uspex_init();
		}
		g_free(uspex_calc);
	} else {
		uspex_gui.have_output=FALSE;
		init_uspex_parameters(&(uspex_gui.calc));
		gui_uspex_init();
	}
/* dialog setup */
	title = g_strdup_printf("USPEX: %s", data->basename);
	dialog = dialog_request(CUSPEX, title, NULL, uspex_cleanup,data);
	g_free(title);
	uspex_gui.window = dialog_window(dialog);
/* --- Outside of notebook */
	GUI_FRAME_WINDOW(uspex_gui.window,frame);
	GUI_VBOX_FRAME(frame,vbox);
/* --- MODEL NAME */
	GUI_LINE_BOX(vbox,hbox);
	GUI_LABEL_BOX(hbox,label,"MODEL NAME");
	GUI_TEXT_ENTRY(hbox,uspex_gui.name,data->basename,TRUE);
GUI_TOOLTIP(uspex_gui.name,"The model name will be used in USPEX files\nas well as GDIS display.");
/* --- CONNECTED OUTPUT */
	GUI_LINE_BOX(vbox,hbox);
	GUI_LABEL_BOX(hbox,label,"CONNECTED OUTPUT");
if(uspex_gui.have_output){
	uspex_output_struct *uspex_output=data->uspex;
	tmp=g_strdup_printf("%s",_UO.calc->filename);
	GUI_TEXT_ENTRY(hbox,uspex_gui.file_entry,tmp,TRUE);
	g_free(tmp);
}else{
	GUI_TEXT_ENTRY(hbox,uspex_gui.file_entry,"",TRUE);

}
GUI_TOOLTIP(uspex_gui.file_entry,"Use the previous result of a USPEX calculation\nto fill the parameters of the current one.\nNOTE: all settings should be checked again. Carefully.");
	GUI_OPEN_BUTTON_BOX(hbox,button,load_parameters_dialog);
/* create frame in box */
	GUI_FRAME_WINDOW(uspex_gui.window,frame);
/* create notebook in frame */
	GUI_NOTE_FRAME(frame,notebook);

/*-----------------*/
/* page 1 -> SET-I */
/*-----------------*/
	GUI_PAGE_NOTE(notebook,page,"SET-I");
/* --- Type & System */
	GUI_FRAME_NOTE(page,frame,"Type & System");
/* create a table in the frame*/
	GUI_TABLE_FRAME(frame,table,5,3);
/* 1st line */
	GUI_COMBOBOX_TABLE(table,uspex_gui.calculationMethod,"calculationMethod: ",0,1,0,1);
	GUI_COMBOBOX_ADD(uspex_gui.calculationMethod,"USPEX");
	GUI_COMBOBOX_ADD(uspex_gui.calculationMethod,"META");
	GUI_COMBOBOX_ADD(uspex_gui.calculationMethod,"VCNEB");
	GUI_COMBOBOX_ADD(uspex_gui.calculationMethod,"PSO");
GUI_TOOLTIP(uspex_gui.calculationMethod,"calculationMethod: Ch. 4.1 DEFAULT: USPEX\nSet the method of calculation.");
	GUI_ENTRY_TABLE(table,uspex_gui.calculationType,uspex_gui.calc.calculationType,"%3i","calculationType:",1,2,0,1);
GUI_TOOLTIP(uspex_gui.calculationMethod,"calculationMethod: Ch. 4.1 DEFAULT: 300\nSets the dimensionality, molecularity, and variability\nof the calculation, can also be set individually.");
	GUI_SPIN_TABLE(table,uspex_gui._calctype_dim,uspex_gui._dim,spin_update_dim,"DIM",2,3,0,1);
GUI_TOOLTIP(uspex_gui._calctype_dim,"Set the dimension of the system.\n3, 2, 1, and 0 are equilvalent to 3D, 2D, 1D, and 0D.\n-2 correspond to 2D crystals.");
	GUI_CHECK_TABLE(table,button,uspex_gui.calc._calctype_mol,mol_toggle,"MOL",3,4,0,1);

	GUI_CHECK_TABLE(table,button,uspex_gui.calc._calctype_var,var_toggle,"VAR",4,5,0,1);
/* 2nd line */
	GUI_COMBOBOX_TABLE(table,uspex_gui.optType,"optType:",0,1,1,2);
	GUI_COMBOBOX_ADD(uspex_gui.optType,"1: MIN Enthalpy");
	GUI_COMBOBOX_ADD(uspex_gui.optType,"2: MIN Volume");
	GUI_COMBOBOX_ADD(uspex_gui.optType,"3: MAX Hardness");
	GUI_COMBOBOX_ADD(uspex_gui.optType,"4: MAX Degree of order");
	GUI_COMBOBOX_ADD(uspex_gui.optType,"5: MAX Structure difference");
	GUI_COMBOBOX_ADD(uspex_gui.optType,"6: MAX Dielect. suscept.");
	GUI_COMBOBOX_ADD(uspex_gui.optType,"7: MAX Band gap");
	GUI_COMBOBOX_ADD(uspex_gui.optType,"8: MAX Dielectric gap");
	GUI_COMBOBOX_ADD(uspex_gui.optType,"9: MAX Magnetization");
	GUI_COMBOBOX_ADD(uspex_gui.optType,"10: MAX Structure quasientropy");
	GUI_COMBOBOX_ADD(uspex_gui.optType,"1101: MAX Bulk modulus");
	GUI_COMBOBOX_ADD(uspex_gui.optType,"1102: MAX Shear modulus");
	GUI_COMBOBOX_ADD(uspex_gui.optType,"1103: MAX Young modulus");
	GUI_COMBOBOX_ADD(uspex_gui.optType,"1104: MAX Poisson Modulus");
	GUI_COMBOBOX_ADD(uspex_gui.optType,"1105: MAX Pugh ratio");
	GUI_COMBOBOX_ADD(uspex_gui.optType,"1106: MAX Vickers hardness");
	GUI_COMBOBOX_ADD(uspex_gui.optType,"1107: MAX Fracture toughness");
	GUI_COMBOBOX_ADD(uspex_gui.optType,"1108: MAX Debye temperature");
	GUI_COMBOBOX_ADD(uspex_gui.optType,"1109: MAX Sound velocity");
	GUI_COMBOBOX_ADD(uspex_gui.optType,"1110: MAX S-wave velocity");
	GUI_COMBOBOX_ADD(uspex_gui.optType,"1111: MAX P-wave velocity");
GUI_TOOLTIP(uspex_gui.optType,"optType: Ch. 4.1 DEFAULT: 1(Enthalpy)\nSelect the properties to optimize.");

/* initialize */
	GUI_COMBOBOX_SETUP(uspex_gui.calculationMethod,0,uspex_method_selected);
	GUI_COMBOBOX_SETUP(uspex_gui.optType,0,uspex_optimization_selected);
	mol_toggle();
	var_toggle();
/* --- end frame */

/*------------------*/
/* page 2 -> SET-II */
/*------------------*/
	GUI_PAGE_NOTE(notebook,page,"SET-II");
/* --- end frame */

/*---------------------*/
/* page 3 -> SPECIFICS */
/*---------------------*/
	GUI_PAGE_NOTE(notebook,page,"SPECIFICS");
/* --- end frame */

/*--------------------*/
/* page 4 -> ABINITIO */
/*--------------------*/
	GUI_PAGE_NOTE(notebook,page,"AB INITIO");
/* --- end frame */

/*----------------*/
/* page 5 -> EXEC */
/*----------------*/
	GUI_PAGE_NOTE(notebook,page,"EXEC");
/* --- end frame */

/* --- Outside of notebook */
	GUI_FRAME_WINDOW(uspex_gui.window,frame);
	GUI_VBOX_FRAME(frame,vbox);
/* Action buttons */
	GUI_SAVE_ACTION(uspex_gui.window,uspex_gui.button_save,save_uspex_calc,NULL);
	GUI_EXEC_ACTION(uspex_gui.window,uspex_gui.button_exec,uspex_exec_calc,NULL);
	GUI_CLOSE_ACTION(uspex_gui.window,button,quit_uspex_gui,dialog);
/* connect to signals */
	GUI_PAGE_CHANGE(notebook,uspex_gui_page_switch,NULL);
/* all done */
	uspex_gui_refresh();/*refresh once more*/
	GUI_SHOW(uspex_gui.window);/*display*/
	sysenv.refresh_dialog=TRUE;
}

