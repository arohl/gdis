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
extern struct elem_pak elements[];

/*globals variables*/
struct uspex_calc_gui uspex_gui;

/*************************/
/* initialize gui values */
/*************************/
void gui_uspex_init(struct model_pak *model){
	gint idx;
	gchar *line;
	uspex_gui.cur_page=USPEX_PAGE_SET_I;
	/*prepare local values*/
	uspex_gui._tmp_atom_typ=0;uspex_gui._tmp_atom_num=0;uspex_gui._tmp_atom_val=0;
	strcpy(uspex_gui._tmp_atom_sym,elements[uspex_gui._tmp_atom_typ].symbol);

	uspex_gui._tmp_specificTrans=NULL;
	if(uspex_gui.calc._nspetrans>0) {
		line=g_strdup("");
		for(idx=0;idx<uspex_gui.calc._nspetrans;idx++) line=g_strdup_printf("%s %i",line,uspex_gui.calc.specificTrans[idx]);
		uspex_gui._tmp_specificTrans=line;
		line=NULL;
	}
	/*get atom information if not available*/
	if(uspex_gui.calc.atomType==NULL){
		GSList *list;
		struct core_pak *core;
		/*first get the _nspecies*/
		list=find_unique(ELEMENT,model);
		uspex_gui.calc._nspecies=(gint)g_slist_length(list);
		uspex_gui.calc.atomType = g_malloc(uspex_gui.calc._nspecies*sizeof(gint));
		for (idx=0 ; list ; list=g_slist_next(list)){
			uspex_gui.calc.atomType[idx]=GPOINTER_TO_INT(list->data);
			idx++;
		}
		g_slist_free(list);
		/*fill the numType with information from current model*/
		if(uspex_gui.calc.numSpecies!=NULL) g_free(uspex_gui.calc.numSpecies);
		uspex_gui.calc._var_nspecies=1;/*we don't automatically vary formula at that point*/
		uspex_gui.calc.numSpecies = g_malloc(uspex_gui.calc._nspecies*sizeof(gint));
		for(idx=0;idx<uspex_gui.calc._nspecies;idx++) uspex_gui.calc.numSpecies[idx]=0;
		for (list=model->cores ; list ; list=g_slist_next(list)){
			idx=0;
			core=list->data;
			/*find corresponding species*/
			while(idx<uspex_gui.calc._nspecies) {
				if(core->atom_code==uspex_gui.calc.atomType[idx]) uspex_gui.calc.numSpecies[idx]++;
				idx++;
			}
		}
		if(uspex_gui.calc.valences!=NULL) g_free(uspex_gui.calc.valences);
		uspex_gui.calc.valences = g_malloc(uspex_gui.calc._nspecies*sizeof(gint));
		for(idx=0;idx<uspex_gui.calc._nspecies;idx++) uspex_gui.calc.valences[idx]=0;
	}
	/*we NEED to check optional valences*/
	if(uspex_gui.calc.valences==NULL){
		uspex_gui.calc.valences = g_malloc(uspex_gui.calc._nspecies*sizeof(gint));
		for(idx=0;idx<uspex_gui.calc._nspecies;idx++) uspex_gui.calc.valences[idx]=0;
	}
	uspex_gui.auto_bonds=(uspex_gui.calc.goodBonds==NULL);
	if(uspex_gui.calc.populationSize==0){
		/*take a default */
		uspex_gui.calc.populationSize=(2*model->num_atoms%10)*10;
		uspex_gui.calc.initialPopSize=uspex_gui.calc.populationSize;
	}
	if(uspex_gui.calc.stopCrit==0) {
		if(uspex_gui.calc._calctype_var) uspex_gui.calc.stopCrit=uspex_gui.calc.maxAt;
		else uspex_gui.calc.stopCrit=model->num_atoms;
	}
	if(uspex_gui.calc.keepBestHM==0) uspex_gui.calc.keepBestHM=(gint)(0.15*(uspex_gui.calc.populationSize));
	if(uspex_gui.calc.symmetries==NULL){
		switch(uspex_gui.calc._calctype_dim){
		case 0:
			uspex_gui.calc.symmetries=g_strdup("E C2 D2 C4 C3 C6 T S2 Ch1 Cv2 S4 S6 Ch3 Th Ch2 Ch4 D3 Ch6 O D4 Cv3 D6 Td Cv4 Dd3 Cv6 Oh C5 S5 S10 Cv5 Ch5 D5 Dd5 Dh5 I Ih");
			break;
		case 1:
			uspex_gui.calc.symmetries=g_strdup("");
			break;
		case 2:
			uspex_gui.calc.symmetries=g_strdup("2-17");
			break;
		case 3:
		default:
			uspex_gui.calc.symmetries=g_strdup("2-230");
		}
	}
	if((uspex_gui.calc.fracPerm==0)&&(uspex_gui.calc._nspecies>1)) uspex_gui.calc.fracPerm=0.1;
	if((uspex_gui.calc.fracRotMut==0)&&(uspex_gui.calc._calctype_mol)) uspex_gui.calc.fracRotMut=0.1;
	if((uspex_gui.calc.fracLatMut==0)&&(uspex_gui.calc.optRelaxType!=1)) uspex_gui.calc.fracLatMut=0.1;
	if((uspex_gui.calc.howManySwaps==0)&&(uspex_gui.calc._nspecies>1)&&(uspex_gui.calc.numSpecies!=NULL)){
		if(uspex_gui.calc.specificSwaps==NULL){
			gint i,j;
			for(i=0;i<(uspex_gui.calc._nspecies-1);i++)
				for(j=i+1;j<uspex_gui.calc._nspecies;j++)
					uspex_gui.calc.howManySwaps+=MIN(uspex_gui.calc.numSpecies[i],uspex_gui.calc.numSpecies[j]);
			uspex_gui.calc.howManySwaps*=0.5;
		}else{
			/*we need to count over only possible swaps*/
			gint i,j;
			gint n_swaps;
			gint *specificSwaps;
			gchar *ptr,*ptr2;
			/*populate specificSwaps*/
			/*1-count*/
			ptr=uspex_gui.calc.specificSwaps;n_swaps=0;
			while(*ptr==' ') ptr++;/*skip initial space (if any)*/
			ptr2=ptr;
			do{
				j=g_ascii_strtod(ptr,&ptr2);
				if(ptr2==ptr) break;
				n_swaps++;
				ptr=ptr2;
			}while(1);
			/*2-populate specificSwaps*/
			specificSwaps = g_malloc(n_swaps*sizeof(gint));
			ptr=uspex_gui.calc.specificSwaps;i=0;
			while(*ptr==' ') ptr++;/*skip initial space (if any)*/
			ptr2=ptr;
			do{
				j=g_ascii_strtod(ptr,&ptr2);
				if(ptr2==ptr) break;
				/*reg j*/
				specificSwaps[i] = j;
				ptr=ptr2;
				i++;
			}while(1);
			for(i=0;i<n_swaps-1;i++)
				for(j=i+1;j<n_swaps;j++)
					uspex_gui.calc.howManySwaps+=MIN(uspex_gui.calc.numSpecies[i],uspex_gui.calc.numSpecies[j]);
			uspex_gui.calc.howManySwaps*=0.5;
			g_free(specificSwaps);
		}
	}
	if(uspex_gui.calc.IonDistances==NULL){
		uspex_gui.calc.IonDistances = g_malloc(uspex_gui.calc._nspecies*uspex_gui.calc._nspecies*sizeof(gdouble));
		for(idx=0;idx<uspex_gui.calc._nspecies*uspex_gui.calc._nspecies;idx++) uspex_gui.calc.IonDistances[idx]=0.;/*no default*/
	}
	uspex_gui.is_dirty=TRUE;
}
/*****************************************************/
/*(re) populate atomType combo with atom information */
/*****************************************************/
void populate_atomType(void){
	gint idx,jdx;
	gchar *line;
	/*this is call on new atomType information, so wipe previous information*/
	GUI_COMBOBOX_WIPE(uspex_gui.atomType);
	GUI_COMBOBOX_WIPE(uspex_gui.goodBonds);
	for(idx=0;idx<uspex_gui.calc._nspecies;idx++){
		line=g_strdup_printf("%s(%i) (V=%i)",
			elements[uspex_gui.calc.atomType[idx]].symbol,uspex_gui.calc.numSpecies[idx],uspex_gui.calc.valences[idx]);
		GUI_COMBOBOX_ADD(uspex_gui.atomType,line);
		g_free(line);
	}
/*only update goodBonds if needed*/
if((uspex_gui.calc.goodBonds!=NULL)&&(uspex_gui.auto_bonds)){
	gchar *tmp;
	for(idx=0;idx<uspex_gui.calc._nspecies;idx++){
		tmp=g_strdup("");
		line=tmp;
		for(jdx=0;jdx<uspex_gui.calc._nspecies;jdx++){
			line=g_strdup_printf("%s %5f",tmp,uspex_gui.calc.goodBonds[jdx+idx*uspex_gui.calc._nspecies]);
			g_free(tmp);tmp=line;
		}
//		line=g_strdup_printf("%s\n",line);/* <- \n not needed*/
		GUI_COMBOBOX_ADD(uspex_gui.goodBonds,line);
	}
}
	/*add ending statement*/
	GUI_COMBOBOX_ADD(uspex_gui.atomType,"ADD ATOMTYPE");
	GUI_COMBOBOX_ADD(uspex_gui.goodBonds,"ADD GOODBOND");

}
/*************************/
/* Refresh SV parameters */
/*************************/
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
			test_calc = read_uspex_parameters(filename,uspex_gui.calc._nspecies);
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
/******************************************/
/* point to the path to run the USPEX job */
/******************************************/
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
	case 4://TPS
		uspex_gui.calc.calculationMethod = US_CM_TPS;break;
	case 5://MINHOP
		uspex_gui.calc.calculationMethod = US_CM_MINHOP;break;
	default://UNKNOWN
		uspex_gui.calc.calculationMethod = US_CM_UNKNOWN;
	}
}
/*****************************/
/* selecting calculationType */
/*****************************/
void uspex_type_selected(GUI_OBJ *w){
	gint index;
	GUI_COMBOBOX_GET(w,index);
	switch (index){
	case 1://s300
		uspex_gui.calc.calculationType=US_CT_s300;
		uspex_gui.calc._calctype_dim=3;
		uspex_gui.calc._calctype_mol=FALSE;
		uspex_gui.calc._calctype_var=FALSE;
		uspex_gui.calc._calctype_mag=TRUE;
		break;
	case 2://301
		uspex_gui.calc.calculationType=US_CT_301;
		uspex_gui.calc._calctype_dim=3;
		uspex_gui.calc._calctype_mol=FALSE;
		uspex_gui.calc._calctype_var=TRUE;
		uspex_gui.calc._calctype_mag=FALSE;
		break;
	case 3://s301
		uspex_gui.calc.calculationType=US_CT_s301;
		uspex_gui.calc._calctype_dim=3;
		uspex_gui.calc._calctype_mol=FALSE;
		uspex_gui.calc._calctype_var=TRUE;
		uspex_gui.calc._calctype_mag=TRUE;
		break;
	case 4://310
		uspex_gui.calc.calculationType=US_CT_310;
		uspex_gui.calc._calctype_dim=3;
		uspex_gui.calc._calctype_mol=TRUE;
		uspex_gui.calc._calctype_var=FALSE;
		uspex_gui.calc._calctype_mag=FALSE;
		break;
	case 5://311
		uspex_gui.calc.calculationType=US_CT_311;
		uspex_gui.calc._calctype_dim=3;
		uspex_gui.calc._calctype_mol=TRUE;
		uspex_gui.calc._calctype_var=TRUE;
		uspex_gui.calc._calctype_mag=FALSE;
		break;
	case 6://000
		uspex_gui.calc.calculationType=US_CT_000;
		uspex_gui.calc._calctype_dim=0;
		uspex_gui.calc._calctype_mol=FALSE;
		uspex_gui.calc._calctype_var=FALSE;
		uspex_gui.calc._calctype_mag=FALSE;
		break;
	case 7://s000
		uspex_gui.calc.calculationType=US_CT_s000;
		uspex_gui.calc._calctype_dim=0;
		uspex_gui.calc._calctype_mol=FALSE;
		uspex_gui.calc._calctype_var=FALSE;
		uspex_gui.calc._calctype_mag=TRUE;
		break;
	case 8://110
		uspex_gui.calc.calculationType=US_CT_110;
		uspex_gui.calc._calctype_dim=1;
		uspex_gui.calc._calctype_mol=TRUE;
		uspex_gui.calc._calctype_var=FALSE;
		uspex_gui.calc._calctype_mag=FALSE;
		break;
	case 9://200
		uspex_gui.calc.calculationType=US_CT_200;
		uspex_gui.calc._calctype_dim=2;
		uspex_gui.calc._calctype_mol=FALSE;
		uspex_gui.calc._calctype_var=FALSE;
		uspex_gui.calc._calctype_mag=FALSE;
		break;
	case 10://s200
		uspex_gui.calc.calculationType=US_CT_s200;
		uspex_gui.calc._calctype_dim=2;
		uspex_gui.calc._calctype_mol=FALSE;
		uspex_gui.calc._calctype_var=FALSE;
		uspex_gui.calc._calctype_mag=TRUE;
		break;
	case 11://201
		uspex_gui.calc.calculationType=US_CT_201;
		uspex_gui.calc._calctype_dim=2;
		uspex_gui.calc._calctype_mol=FALSE;
		uspex_gui.calc._calctype_var=TRUE;
		uspex_gui.calc._calctype_mag=FALSE;
		break;
	case 12://s201
		uspex_gui.calc.calculationType=US_CT_s201;
		uspex_gui.calc._calctype_dim=2;
		uspex_gui.calc._calctype_mol=FALSE;
		uspex_gui.calc._calctype_var=TRUE;
		uspex_gui.calc._calctype_mag=TRUE;
		break;
	case 13://-200
		uspex_gui.calc.calculationType=US_CT_m200;
		uspex_gui.calc._calctype_dim=-2;
		uspex_gui.calc._calctype_mol=FALSE;
		uspex_gui.calc._calctype_var=FALSE;
		uspex_gui.calc._calctype_mag=FALSE;
		break;
	case 14://-s200
		uspex_gui.calc.calculationType=US_CT_sm200;
		uspex_gui.calc._calctype_dim=-2;
		uspex_gui.calc._calctype_mol=FALSE;
		uspex_gui.calc._calctype_var=FALSE;
		uspex_gui.calc._calctype_mag=TRUE;
		break;
	case 0://300
	default://bad value -> default 300
		uspex_gui.calc.calculationType=US_CT_300;
		uspex_gui.calc._calctype_dim=3;
		uspex_gui.calc._calctype_mol=FALSE;
		uspex_gui.calc._calctype_var=FALSE;
		uspex_gui.calc._calctype_mag=FALSE;
	}
	/*now update uspex_gui._calctype_dim, uspex_gui._calctype_mol, and uspex_gui._calctype_var*/
	GUI_SPIN_SET(uspex_gui._calctype_dim,(gdouble)uspex_gui.calc._calctype_dim);
	if(uspex_gui.calc._calctype_mol) GUI_TOGGLE_ON(uspex_gui._calctype_mol);
	else GUI_TOGGLE_OFF(uspex_gui._calctype_mol);
	if(uspex_gui.calc._calctype_var) GUI_TOGGLE_ON(uspex_gui._calctype_var);
	else GUI_TOGGLE_OFF(uspex_gui._calctype_var);
	if(uspex_gui.calc._calctype_mag) GUI_TOGGLE_ON(uspex_gui._calctype_mag);
	else GUI_TOGGLE_OFF(uspex_gui._calctype_mag);
}
/**************************/
/* update calculationType */
/**************************/
void _update_calculationType(){
	gint i=0;
	switch(uspex_gui.calc._calctype_dim){
	case 3://300,301,310,311 - VER 10.1 s300 s301
		i=300;
		break;
	case 2://200,201 - VER 10.1 s200 s201
		i=200;
		uspex_gui.calc._calctype_mol=FALSE;
		break;
	case 1://110
		i=100;
		uspex_gui.calc._calctype_mag=FALSE;
		uspex_gui.calc._calctype_mol=TRUE;
		uspex_gui.calc._calctype_var=FALSE;
		break;
	case 0://000 - VER 10.1 s000
		i=0;
		uspex_gui.calc._calctype_mol=FALSE;
		uspex_gui.calc._calctype_var=FALSE;
		break;
	case -2://-200 - VER 10.1 -s200
		i=200;
		uspex_gui.calc._calctype_mol=FALSE;
		uspex_gui.calc._calctype_var=FALSE;
		break;
	default://should not happen
		i=300;/*reset to default value*/
		uspex_gui.calc._calctype_mag=FALSE;
		uspex_gui.calc._calctype_mol=FALSE;
		uspex_gui.calc._calctype_var=FALSE;
	}
	if(uspex_gui.calc._calctype_mag) i+=1000;
	if(uspex_gui.calc._calctype_mol) i+=10;
	if(uspex_gui.calc._calctype_var) i+=1;
	if(uspex_gui.calc._calctype_dim==-2) i*=-1;
	uspex_gui.calc.calculationType=i;
	/*reset toogles*/
	if(uspex_gui.calc._calctype_mol) GUI_TOGGLE_ON(uspex_gui._calctype_mol);
	else GUI_TOGGLE_OFF(uspex_gui._calctype_mol);
	if(uspex_gui.calc._calctype_var) GUI_TOGGLE_ON(uspex_gui._calctype_var);
	else GUI_TOGGLE_OFF(uspex_gui._calctype_var);
	if(uspex_gui.calc._calctype_mag) GUI_TOGGLE_ON(uspex_gui._calctype_mag);
	else GUI_TOGGLE_OFF(uspex_gui._calctype_mag);
	/*reset select*/
	switch(uspex_gui.calc.calculationType){
	case US_CT_s300:
		GUI_COMBOBOX_SET(uspex_gui.calculationType,1);break;
	case US_CT_301:
		GUI_COMBOBOX_SET(uspex_gui.calculationType,2);break;
	case US_CT_s301:
		GUI_COMBOBOX_SET(uspex_gui.calculationType,3);break;
	case US_CT_310:
		GUI_COMBOBOX_SET(uspex_gui.calculationType,4);break;
	case US_CT_311:
		GUI_COMBOBOX_SET(uspex_gui.calculationType,5);break;
	case US_CT_000:
		GUI_COMBOBOX_SET(uspex_gui.calculationType,6);break;
	case US_CT_s000:
		GUI_COMBOBOX_SET(uspex_gui.calculationType,7);break;
	case US_CT_110:
		GUI_COMBOBOX_SET(uspex_gui.calculationType,8);break;
	case US_CT_200:
		GUI_COMBOBOX_SET(uspex_gui.calculationType,9);break;
	case US_CT_s200:
		GUI_COMBOBOX_SET(uspex_gui.calculationType,10);break;
	case US_CT_201:
		GUI_COMBOBOX_SET(uspex_gui.calculationType,11);break;
	case US_CT_s201:
		GUI_COMBOBOX_SET(uspex_gui.calculationType,12);break;
	case US_CT_m200:
		GUI_COMBOBOX_SET(uspex_gui.calculationType,13);break;
	case US_CT_sm200:
		GUI_COMBOBOX_SET(uspex_gui.calculationType,14);break;
	case US_CT_300:
	default:
		GUI_COMBOBOX_SET(uspex_gui.calculationType,0);
	}
}
/************************/
/* update _calctype_dim */
/************************/
void spin_update_dim(){
	gint i=(gint)uspex_gui._dim;
	/*ban the -1 value*/
	if(i==-1){
		if(uspex_gui.calc._calctype_dim==0) i=-2;
		else i=0;
	}
	/*update calc value*/
	uspex_gui.calc._calctype_dim=i;
	/*update spin value*/
	uspex_gui._dim=(gdouble)i;
	GUI_SPIN_SET(uspex_gui._calctype_dim,uspex_gui._dim);
	/*done, now update calculationType*/
	_update_calculationType();
}
/*****************************************/
/* toggle mol -> update calculation type */
/*****************************************/
void mol_toggle(void){
	gint dim=uspex_gui.calc._calctype_dim*100;
	gint sgn;
	if(dim<0) {
		dim*=-1;
		sgn=-1;
	}else sgn=1;
	if(uspex_gui.calc._calctype_mol){
		if(uspex_gui.calc._calctype_var) uspex_gui.calc.calculationType = dim+10+1;
		else uspex_gui.calc.calculationType = dim+10;
	}else{
		if(uspex_gui.calc._calctype_mag){
			if(uspex_gui.calc._calctype_var) uspex_gui.calc.calculationType = 1000+dim+1;
			else uspex_gui.calc.calculationType = 1000+dim;
		}else{
			if(uspex_gui.calc._calctype_var) uspex_gui.calc.calculationType = dim+1;
			else uspex_gui.calc.calculationType = dim;
		}
	}
	uspex_gui.calc.calculationType*=sgn;
	/*done, now update calculationType*/
	_update_calculationType();
}
/*****************************************/
/* toggle var -> update calculation type */
/*****************************************/
void var_toggle(void){
	gint dim=uspex_gui.calc._calctype_dim*100;
	gint sgn; 
	if(dim<0) {
		dim*=-1;
		sgn=-1;
	}else sgn=1;
	if(uspex_gui.calc._calctype_var){
		if(uspex_gui.calc._calctype_mol) uspex_gui.calc.calculationType = dim+10+1;
		else 	if(uspex_gui.calc._calctype_mag) uspex_gui.calc.calculationType = 1000+dim+1;
			else uspex_gui.calc.calculationType = dim+1;
	}else{
		if(uspex_gui.calc._calctype_mol) uspex_gui.calc.calculationType = dim+10;
		else 	if(uspex_gui.calc._calctype_mag) uspex_gui.calc.calculationType = 1000+dim;
			else uspex_gui.calc.calculationType = dim;
	}
	uspex_gui.calc.calculationType*=sgn;
	/*done, now update calculationType*/
	_update_calculationType();
}
/*****************************************/
/* mag_toggle -> update calculation type */
/*****************************************/
void mag_toggle(void){
	gint dim=uspex_gui.calc._calctype_dim*100;
	gint sgn; 
	if(dim<0) {
		dim*=-1;
		sgn=-1;
	}else sgn=1;
	if(uspex_gui.calc._calctype_mag){
		if(uspex_gui.calc._calctype_var) uspex_gui.calc.calculationType = 1000+dim+1;
		else uspex_gui.calc.calculationType = 1000+dim;
	}else{
		if(uspex_gui.calc._calctype_mol){
			if(uspex_gui.calc._calctype_var) uspex_gui.calc.calculationType = dim+10+1;
			else uspex_gui.calc.calculationType = dim+10;
		}else{
			if(uspex_gui.calc._calctype_var) uspex_gui.calc.calculationType = dim+1;
			else uspex_gui.calc.calculationType = dim;
		}
	}
	uspex_gui.calc.calculationType*=sgn;
	/*done, now update calculationType*/
	_update_calculationType();
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
	case 10://11
		uspex_gui.calc.optType=US_OT_2R;break;
	case 11://14
		uspex_gui.calc.optType=US_OT_ZT;break;
	case 12://17
		uspex_gui.calc.optType=US_OT_Fphon;break;
	case 13://1101
		uspex_gui.calc.optType=US_OT_BULK_M;break;
	case 14://1102
		uspex_gui.calc.optType=US_OT_SHEAR_M;break;
	case 15://1103
		uspex_gui.calc.optType=US_OT_YOUNG_M;break;
	case 16://1104
		uspex_gui.calc.optType=US_OT_POISSON;break;
	case 17://1105
		uspex_gui.calc.optType=US_OT_PUGH_R;break;
	case 18://1106
		uspex_gui.calc.optType=US_OT_VICKERS_H;break;
	case 19://1107
		uspex_gui.calc.optType=US_OT_FRACTURE;break;
	case 20://1108
		uspex_gui.calc.optType=US_OT_DEBYE_T;break;
	case 21://1109
		uspex_gui.calc.optType=US_OT_SOUND_V;break;
	case 22://1110
		uspex_gui.calc.optType=US_OT_SWAVE_V;break;
	case 23://1111
		uspex_gui.calc.optType=US_OT_PWAVE_V;break;
	default:
		uspex_gui.calc.optType=US_OT_UNKNOWN;
	}
}
/**********************/
/* Selecting atomType */
/**********************/
void atomType_selected(GUI_OBJ *w){
        gchar *text;
        gchar sym[3];
        GUI_COMBOBOX_GET_TEXT(w,text);
	if (g_ascii_strcasecmp(text,"ADD ATOMTYPE") == 0) return;/*last selection*/
	sscanf(text,"%[^(](%i) (V=%i)",sym,&(uspex_gui._tmp_atom_num),&(uspex_gui._tmp_atom_val));
	sym[2]='\0';uspex_gui._tmp_atom_typ=elem_symbol_test(sym);
	strcpy(uspex_gui._tmp_atom_sym,sym);
	g_free(text);
	/*update entries*/
	GUI_ENTRY_TEXT(uspex_gui._atom_sym,uspex_gui._tmp_atom_sym);
	text=g_strdup_printf("%i",uspex_gui._tmp_atom_typ);
	GUI_ENTRY_TEXT(uspex_gui._atom_typ,text);
	g_free(text);text=g_strdup_printf("%i",uspex_gui._tmp_atom_num);
	GUI_ENTRY_TEXT(uspex_gui._atom_num,text);
	g_free(text);text=g_strdup_printf("%i",uspex_gui._tmp_atom_val);
	GUI_ENTRY_TEXT(uspex_gui._atom_val,text);
        g_free(text);
}
/*****************************************************************/
/* Change/ADD atomType, numSpecies, and valence of selected atom */
/*****************************************************************/
void apply_atom(){
	gint index;
	gchar *text;
	gboolean exists=FALSE;
	/*information*/
	gchar sym[3];
	gchar tmp[3];
	gint typ;
	gint num;
	gint val;
	/**/
	GUI_COMBOBOX_GET(uspex_gui.atomType,index);
	if(index==-1) {/*nothing selected -> should never happen*/
		GUI_COMBOBOX_SET(uspex_gui.atomType,0);
		return;
	}
	GUI_COMBOBOX_GET_TEXT(uspex_gui.atomType,text);
	/*mini-sync*/
	GUI_REG_VAL(uspex_gui._atom_sym,sym[0],"%s");
	GUI_REG_VAL(uspex_gui._atom_num,num,"%i");
	GUI_REG_VAL(uspex_gui._atom_val,val,"%i");
	typ = elem_symbol_test(sym);
	/*trick to have a valid sym*/
	strcpy(sym,elements[typ].symbol);
	if (g_ascii_strcasecmp(text,"ADD ATOMTYPE") == 0) {
		g_free(text);
		/*add a new atomType*/
		/*check if atomType does not exists*/
		index=0;
		GUI_COMBOBOX_SET(uspex_gui.atomType,0);
		GUI_COMBOBOX_GET_TEXT(uspex_gui.atomType,text);
		while(g_ascii_strcasecmp(text,"ADD ATOMTYPE") != 0){
			sscanf(text,"%[^(](%*i) (V=%*i)",tmp);
			uspex_gui._tmp_atom_typ=elem_symbol_test(tmp);
			if(uspex_gui._tmp_atom_typ==typ) exists=TRUE;
			g_free(text);
			index++;
			GUI_COMBOBOX_SET(uspex_gui.atomType,index);
			GUI_COMBOBOX_GET_TEXT(uspex_gui.atomType,text);
		}
		if(exists){
			text = g_strdup_printf("Atom type %s already exists! Can't add to atomType!\n",sym);
			gui_text_show(ERROR,text);
			g_free(text);
			return;
		}else{
			text = g_strdup_printf("%s(%i) (V=%i)",sym,num,val);
			GUI_COMBOBOX_ADD_TEXT(uspex_gui.atomType,index,text);
			GUI_COMBOBOX_SET(uspex_gui.atomType,index+1);/*select last*/
			g_free(text);
		}
	}else{
		/*get current typ*/
		sscanf(text,"%[^(](%*i) (V=%*i)",tmp);
		uspex_gui._tmp_atom_typ = elem_symbol_test(tmp);
		g_free(text);
		/*modify an atomType*/
		if(typ!=uspex_gui._tmp_atom_typ){
			gint index_bkup;
			/*type have changed, check if new type exists*/
			index_bkup=index;/*SAVE INDEX*/
			index=0;
			GUI_COMBOBOX_SET(uspex_gui.atomType,0);
			GUI_COMBOBOX_GET_TEXT(uspex_gui.atomType,text);
			while(g_ascii_strcasecmp(text,"ADD ATOMTYPE") != 0){
				sscanf(text,"%[^(](%*i) (V=%*i)",tmp);
				uspex_gui._tmp_atom_typ=elem_symbol_test(tmp);
				if(uspex_gui._tmp_atom_typ==typ) exists=TRUE;
				g_free(text);
				index++;
				GUI_COMBOBOX_SET(uspex_gui.atomType,index);
				GUI_COMBOBOX_GET_TEXT(uspex_gui.atomType,text);
			}
			if(exists){
				text = g_strdup_printf("Atom type %s already exists! Can't modify atomType!\n",sym);
				gui_text_show(ERROR,text);
				g_free(text);
				return;
			}
			index=index_bkup;/*RESTORE INDEX*/
		}
		/*just update values*/
		text = g_strdup_printf("%s(%i) (V=%i)",sym,num,val);
		GUI_COMBOBOX_ADD_TEXT(uspex_gui.atomType,index,text);
		GUI_COMBOBOX_DEL(uspex_gui.atomType,index+1);
		g_free(text);
	}
	/*select "ADD ATOM" (for convenience)*/
	GUI_COMBOBOX_SET(uspex_gui.atomType,index);
	GUI_COMBOBOX_GET_TEXT(uspex_gui.atomType,text);
	while(g_ascii_strcasecmp(text,"ADD ATOMTYPE") != 0){
		g_free(text);
		index++;
		GUI_COMBOBOX_SET(uspex_gui.atomType,index);
		GUI_COMBOBOX_GET_TEXT(uspex_gui.atomType,text);
	}
	g_free(text);
}
/****************************************/
/* Remove selected atomType information */
/****************************************/
void remove_atom(){
	gint index;
	gchar *text;
	GUI_COMBOBOX_GET(uspex_gui.atomType,index);
	if(index==-1) {/*nothing selected -> should never happen*/
		GUI_COMBOBOX_SET(uspex_gui.atomType,0);
		return;
	}
	GUI_COMBOBOX_GET_TEXT(uspex_gui.atomType,text);
	if (g_ascii_strcasecmp(text,"ADD ATOMTYPE") == 0) {
		/*can't delete this one*/
		g_free(text);
		return;
	}
	GUI_COMBOBOX_DEL(uspex_gui.atomType,index);
	if(index-1<0) index=1;
	GUI_COMBOBOX_SET(uspex_gui.atomType,index-1);
	g_free(text);
}
/*********************/
/* toggle auto_bonds */
/*********************/
void auto_bond_toggle(void){
	if(uspex_gui.auto_bonds){
		GUI_LOCK(uspex_gui.goodBonds);
		GUI_LOCK(uspex_gui._bond_d);
	}else{
		GUI_UNLOCK(uspex_gui.goodBonds);
		GUI_UNLOCK(uspex_gui._bond_d);
	}
}
/**********************/
/* Selecting goodBond */
/**********************/
void goodBonds_selected(GUI_OBJ *w){
        gchar *text;
	/**/
        GUI_COMBOBOX_GET_TEXT(w,text);
        if (g_ascii_strcasecmp(text,"ADD GOODBOND") == 0) return;/*last selection*/
	/*the line is going directly into _bond_d*/
	if(uspex_gui._tmp_bond_d!=NULL) g_free(uspex_gui._tmp_bond_d);
	uspex_gui._tmp_bond_d=g_strdup(text);
	GUI_ENTRY_TEXT(uspex_gui._bond_d,uspex_gui._tmp_bond_d);
	g_free(text);
}
/************************************/
/* Change/ADD goodBonds information */
/************************************/
void apply_bonds(){
	gint index;
	gchar *text;
	/**/
	GUI_COMBOBOX_GET(uspex_gui.goodBonds,index);
	GUI_COMBOBOX_GET_TEXT(uspex_gui.goodBonds,text);
	if(uspex_gui._tmp_bond_d!=NULL) g_free(uspex_gui._tmp_bond_d);
	GUI_ENTRY_GET_TEXT(uspex_gui._bond_d,uspex_gui._tmp_bond_d);
	GUI_COMBOBOX_ADD_TEXT(uspex_gui.goodBonds,index,uspex_gui._tmp_bond_d);
	if (g_ascii_strcasecmp(text,"ADD GOODBOND") == 0){
		/*new line*/
		GUI_COMBOBOX_SET(uspex_gui.goodBonds,index+1);
	}else{
		/*modify*/
		GUI_COMBOBOX_DEL(uspex_gui.goodBonds,index+1);
		g_free(text);
		/*select "ADD GOODBOND" (for convenience)*/
		GUI_COMBOBOX_SET(uspex_gui.goodBonds,index);
		GUI_COMBOBOX_GET_TEXT(uspex_gui.goodBonds,text);
		while(g_ascii_strcasecmp(text,"ADD GOODBOND") != 0){
			g_free(text);
			index++;
			GUI_COMBOBOX_SET(uspex_gui.goodBonds,index);
			GUI_COMBOBOX_GET_TEXT(uspex_gui.goodBonds,text);
		}
	}
	g_free(text);
}
/********************************/
/* Remove goodBonds information */
/********************************/
void remove_bonds(){
	gint index;
	gchar *text;
	GUI_COMBOBOX_GET(uspex_gui.goodBonds,index);
	if(index==-1) {/*nothing selected -> should never happen*/
		GUI_COMBOBOX_SET(uspex_gui.goodBonds,0);
		return;
	}
	GUI_COMBOBOX_GET_TEXT(uspex_gui.goodBonds,text);
	if (g_ascii_strcasecmp(text,"ADD GOODBOND") == 0){
		/*can't delete this one*/
		g_free(text);
		return;
	}
	GUI_COMBOBOX_DEL(uspex_gui.goodBonds,index);
	if(index-1<0) index=1;
	GUI_COMBOBOX_SET(uspex_gui.goodBonds,index-1);
	g_free(text);
}
/**********************/
/* toggle new_optType */
/**********************/
void opt_toggle(){
	if(uspex_gui.have_new_opt){
		GUI_UNLOCK(uspex_gui.new_optType);
	}else{
		GUI_LOCK(uspex_gui.new_optType);
	}
}
/**************************************/
/* Apply changes to IonDistances line */
/**************************************/
void apply_distances(void){
	gint index;
	gchar *ptr,*ptr2;
	gint i;
	gint ndistval;
	gdouble *distval=NULL;
	gdouble *tmp=NULL;
	gdouble dist;
	/**/
	GUI_COMBOBOX_GET(uspex_gui.IonDistances,index);
	GUI_ENTRY_GET_TEXT(uspex_gui._distances,uspex_gui._tmp_distances);
	/*get the number of distance values*/
	ptr=&(uspex_gui._tmp_distances[0]);
	while(*ptr==' ') ptr++;/*skip leading white space, if any*/
	ptr2=ptr;ndistval=0;
	do{
		dist=g_ascii_strtod(ptr,&ptr2);
		if(ptr2==ptr) break;
		ndistval++;
		/*update array*/
		tmp = g_malloc(ndistval*sizeof(gdouble));
		if(distval){
			for(i=0;i<ndistval-1;i++) tmp[i]=distval[i];
			g_free(distval);
		}
		tmp[ndistval-1]=dist;
		distval=tmp;
		ptr=ptr2;
	}while(1);
	/*now prepare a *valid* update for IonDistances line*/
	if(ndistval>=uspex_gui.calc._nspecies){
		/*use only nspecies values*/
		ptr2=g_strdup_printf("%4f",distval[0]);
		ptr=ptr2;
		for(i=1;i<uspex_gui.calc._nspecies;i++){
			ptr=g_strdup_printf("%s %4f",ptr2,distval[i]);
			g_free(ptr2);ptr2=ptr;
		}
	}else{
		/*use first ndistval values, pad with zero*/
		ptr2=g_strdup_printf("%4f",distval[0]);
		ptr=ptr2;
		for(i=1;i<ndistval;i++){
			ptr=g_strdup_printf("%s %4f",ptr2,distval[i]);
			g_free(ptr2);ptr2=ptr;
		}
		for(i=ndistval;i<uspex_gui.calc._nspecies;i++){
			ptr=g_strdup_printf("%s %4f",ptr2,0.);
			g_free(ptr2);ptr2=ptr;
		}
	}
	g_free(distval);
	GUI_COMBOBOX_ADD_TEXT(uspex_gui.IonDistances,index,ptr);
	GUI_COMBOBOX_DEL(uspex_gui.IonDistances,index+1);
	GUI_COMBOBOX_SET(uspex_gui.IonDistances,index);
	g_free(ptr);
}
/************************************/
/* Apply changes to MolCenters line */
/************************************/
void apply_centers(void){
	gint index;
	gchar *ptr,*ptr2;
	gint i;
	gint nmolval;
	gdouble *molval=NULL;
	gdouble *tmp=NULL;
	gdouble val;
	/**/
	if(uspex_gui.calc._nmolecules==0) return;
	GUI_COMBOBOX_GET(uspex_gui.MolCenters,index);
	GUI_ENTRY_GET_TEXT(uspex_gui._centers,uspex_gui._tmp_centers);
	/*get the number of molecular distance values*/
	ptr=&(uspex_gui._tmp_centers[0]);
	while(*ptr==' ') ptr++;/*skip leading white space, if any*/
	ptr2=ptr;nmolval=0;
	do{
		val=g_ascii_strtod(ptr,&ptr2);
		if(ptr2==ptr) break;
		nmolval++;
		/*update array*/
		tmp = g_malloc(nmolval*sizeof(gdouble));
		if(molval){
			for(i=0;i<nmolval-1;i++) tmp[i]=molval[i];
			g_free(molval);
		}
		tmp[nmolval-1]=val;
		molval=tmp;
		ptr=ptr2;
	}while(1);
	/*now prepare a *valid* update for MolCenters line*/
	if(nmolval>=uspex_gui.calc._nmolecules){
		/*use only nmolecules values*/
		ptr2=g_strdup_printf("%4f",molval[0]);
		ptr=ptr2;
		for(i=1;i<uspex_gui.calc._nmolecules;i++){
			ptr=g_strdup_printf("%s %4f",ptr2,molval[i]);
			g_free(ptr2);ptr2=ptr;
		}
	}else{
		/*use first nmolval values, pad with zero*/
		ptr2=g_strdup_printf("%4f",molval[0]);
		ptr=ptr2;
		for(i=1;i<nmolval;i++){
			ptr=g_strdup_printf("%s %4f",ptr2,molval[i]);
			g_free(ptr2);ptr2=ptr;
		}
		for(i=nmolval;i<uspex_gui.calc._nmolecules;i++){
			ptr=g_strdup_printf("%s %4f",ptr2,0.);
			g_free(ptr2);ptr2=ptr;
		}
	}
	g_free(molval);
	GUI_COMBOBOX_ADD_TEXT(uspex_gui.MolCenters,index,ptr);
	GUI_COMBOBOX_DEL(uspex_gui.MolCenters,index+1);
	GUI_COMBOBOX_SET(uspex_gui.MolCenters,index);
	g_free(ptr);
}
/*********************/
/* toggle auto_C_ion */
/*********************/
void toggle_auto_C_ion(void){
	gint i,j;
	gchar *text;
	gchar *tamp;
	if(uspex_gui.auto_C_ion){
		/*automatic*/
		GUI_LOCK(uspex_gui.IonDistances);
		GUI_LOCK(uspex_gui._distances);
	}else{
		/*manual*/
		GUI_UNLOCK(uspex_gui.IonDistances);
		GUI_UNLOCK(uspex_gui._distances);
		/*wipe IonDistances, rewrite*/
		GUI_COMBOBOX_WIPE(uspex_gui.IonDistances);
		for(i=0;i<uspex_gui.calc._nspecies;i++){
			tamp=g_strdup_printf("%4f",uspex_gui.calc.IonDistances[i*uspex_gui.calc._nspecies]);
			text=tamp;
			for(j=1;j<uspex_gui.calc._nspecies;j++) {
				text=g_strdup_printf("%s %4f",tamp,uspex_gui.calc.IonDistances[j+i*uspex_gui.calc._nspecies]);
				g_free(tamp);
				tamp=text;
			}
			GUI_COMBOBOX_ADD(uspex_gui.IonDistances,text);
			g_free(text);
		}
		GUI_COMBOBOX_SET(uspex_gui.IonDistances,0);
	}
}
/*********************************/
/* Refresh the Constraints frame */
/*********************************/
void refresh_constraints(void){
	gint i,j;
	gchar *text;
	gchar *tamp;
	if(uspex_gui.calc._nmolecules>0){
		/*refresh MolCenters*/
		GUI_UNLOCK(uspex_gui.MolCenters);
		GUI_UNLOCK(uspex_gui._centers);
		/*wipe MolCenters, rewrite*/
		GUI_COMBOBOX_WIPE(uspex_gui.MolCenters);
		for(i=0;i<uspex_gui.calc._nmolecules;i++){
			tamp=g_strdup_printf("%4f",uspex_gui.calc.MolCenters[i*uspex_gui.calc._nmolecules]);
			text=tamp;
			for(j=1;j<uspex_gui.calc._nmolecules;j++) {
				text=g_strdup_printf("%s %4f",tamp,uspex_gui.calc.MolCenters[j+i*uspex_gui.calc._nmolecules]);
				g_free(tamp);
				tamp=text;
			}
			GUI_COMBOBOX_ADD(uspex_gui.MolCenters,text);
			g_free(text);
		}
		GUI_COMBOBOX_SET(uspex_gui.MolCenters,0);
	}else{
		GUI_LOCK(uspex_gui.MolCenters);
		GUI_LOCK(uspex_gui._centers);
	}
	
}
/**************************/
/* Selecting IonDistances */
/**************************/
void uspex_IonDistances_selected(GUI_OBJ *w){
	gchar *text;
	GUI_COMBOBOX_GET_TEXT(w,text);
	GUI_ENTRY_TEXT(uspex_gui._distances,text);
	g_free(text);
}
/************************/
/* Selecting MolCenters */
/************************/
void uspex_MolCenters_selected(GUI_OBJ *w){
	gchar *text;
	GUI_COMBOBOX_GET_TEXT(w,text);
	GUI_ENTRY_TEXT(uspex_gui._centers,text);
	g_free(text);
}
/***************************************/
/* Apply changes to Lattivevalues line */
/***************************************/
void apply_latticevalue(void){
	gint index;
	gchar *ptr,*ptr2;
	gint i,format;
	gint nlval=0;
	gdouble *lval=NULL;
	gdouble *tmp=NULL;
	gdouble val;
	GUI_COMBOBOX_GET(uspex_gui.Latticevalues,index);
	GUI_ENTRY_GET_TEXT(uspex_gui._latticevalue,uspex_gui._tmp_latticevalue);
	/*get the number of lattice values*/
	ptr=&(uspex_gui._tmp_latticevalue[0]);
	while(*ptr==' ') ptr++;/*skip leading white space, if any*/
	ptr2=ptr;nlval=0;
	do{
		val=g_ascii_strtod(ptr,&ptr2);
		if(ptr2==ptr) break;
		nlval++;
		/*update array*/
		tmp = g_malloc(nlval*sizeof(gdouble));
		if(lval){
			for(i=0;i<nlval-1;i++) tmp[i]=lval[i];
			g_free(lval);
		}
		tmp[nlval-1]=val;
		lval=tmp;
		ptr=ptr2;
	}while(1);
	/*now prepare a *valid* update for Latticevalues line*/
	if(nlval<1) return;/*invalid line*/
	GUI_COMBOBOX_GET(uspex_gui._latticeformat,format);
	switch(format){
		case 1://lattice parameters
			switch(nlval){
			case 1:
				ptr=g_strdup_printf("%4f %4f %4f",lval[0],0.,0.);
				break;
			case 2:
				ptr=g_strdup_printf("%4f %4f %4f",lval[0],lval[1],0.);
				break;
			case 3:
			default://might skip bad/extra values
				ptr=g_strdup_printf("%4f %4f %4f",lval[0],lval[1],lval[2]);
				break;
			}
			break;
		case 2://crystal definition
			if(nlval>=6){/*drop every values over the 6th*/
				ptr=g_strdup_printf("%4f %4f %4f %4f %4f %4f",lval[0],lval[1],lval[2],lval[3],lval[4],lval[5]);
			}else if(nlval==4){/*special definition*/
				 ptr=g_strdup_printf("%4f %4f %4f %4f %4f %4f",lval[0],lval[1],0.,lval[2],lval[3],0.);
			}else{/*zero padded*/
				ptr2=g_strdup_printf("%4f",lval[0]);
				for(i=1;i<nlval;i++) {
					ptr=g_strdup_printf("%s %4f",ptr2,lval[i]);
					g_free(ptr2);ptr2=ptr;
				}
				for(i=nlval;i<6;i++) {
					ptr=g_strdup_printf("%s %4f",ptr2,0.);
					g_free(ptr2);ptr2=ptr;
				}
			}
			break;
		case 0:
		default:
			/*any nlval>1 is valid*/
			ptr2=g_strdup_printf("%4f",lval[0]);
			ptr=ptr2;
			for(i=1;i<nlval;i++) {
				ptr=g_strdup_printf("%s %4f",ptr2,lval[i]);
				g_free(ptr2);ptr2=ptr;
			}
	}
	g_free(lval);
	GUI_COMBOBOX_ADD_TEXT(uspex_gui.Latticevalues,index,ptr);
	GUI_COMBOBOX_DEL(uspex_gui.Latticevalues,index+1);
	GUI_COMBOBOX_SET(uspex_gui.Latticevalues,index);
	g_free(ptr);
}
/********************/
/* toggle auto_lval */
/********************/
void toggle_auto_lval(void){
	if(uspex_gui.auto_lval){
		/*automatic*/
		GUI_LOCK(uspex_gui.Latticevalues);
		GUI_LOCK(uspex_gui._latticevalue);
		GUI_LOCK(uspex_gui._latticeformat);
	}else{
		/*manual*/
		GUI_UNLOCK(uspex_gui.Latticevalues);
		GUI_UNLOCK(uspex_gui._latticevalue);
		GUI_UNLOCK(uspex_gui._latticeformat);
	}
}
/***************************/
/* Selecting latticeformat */
/***************************/
void uspex_latticeformat_selected(GUI_OBJ *w){
	gint index;
	gchar *text;
	GUI_COMBOBOX_GET(w,index);
	GUI_COMBOBOX_WIPE(uspex_gui.Latticevalues);/*wipe in any case*/
	switch (index){
	case 1://lattice parameters
		switch (uspex_gui.calc._nlatticevalues){
		case 1:
			text=g_strdup_printf("%4f %4f %4f",uspex_gui.calc.Latticevalues[0],0.,0.);
			GUI_COMBOBOX_ADD(uspex_gui.Latticevalues,text);
			g_free(text);text=g_strdup_printf("%4f %4f %4f",0.,0.,0.);
			GUI_COMBOBOX_ADD(uspex_gui.Latticevalues,text);
			GUI_COMBOBOX_ADD(uspex_gui.Latticevalues,text);
			g_free(text);
			break;
		case 2:
			text=g_strdup_printf("%4f %4f %4f",uspex_gui.calc.Latticevalues[0],uspex_gui.calc.Latticevalues[1],0.);
			GUI_COMBOBOX_ADD(uspex_gui.Latticevalues,text);
			g_free(text);text=g_strdup_printf("%4f %4f %4f",uspex_gui.calc.Latticevalues[2],uspex_gui.calc.Latticevalues[3],0.);
			GUI_COMBOBOX_ADD(uspex_gui.Latticevalues,text);
			g_free(text);text=g_strdup_printf("%4f %4f %4f",0.,0.,0.);
			GUI_COMBOBOX_ADD(uspex_gui.Latticevalues,text);
			g_free(text);
			break;
		case 3:
			text=g_strdup_printf("%4f %4f %4f",
				uspex_gui.calc.Latticevalues[0],uspex_gui.calc.Latticevalues[1],uspex_gui.calc.Latticevalues[2]);
			GUI_COMBOBOX_ADD(uspex_gui.Latticevalues,text);
			g_free(text);text=g_strdup_printf("%4f %4f %4f",
				uspex_gui.calc.Latticevalues[3],uspex_gui.calc.Latticevalues[4],uspex_gui.calc.Latticevalues[5]);
			GUI_COMBOBOX_ADD(uspex_gui.Latticevalues,text);
			g_free(text);text=g_strdup_printf("%4f %4f %4f",
				uspex_gui.calc.Latticevalues[6],uspex_gui.calc.Latticevalues[7],uspex_gui.calc.Latticevalues[8]);
			GUI_COMBOBOX_ADD(uspex_gui.Latticevalues,text);
			g_free(text);
			break;
		default:
			/*create an null lattice*/
			text=g_strdup_printf("%4f %4f %4f",0.,0.,0.);
			GUI_COMBOBOX_ADD(uspex_gui.Latticevalues,text);
			GUI_COMBOBOX_ADD(uspex_gui.Latticevalues,text);
			GUI_COMBOBOX_ADD(uspex_gui.Latticevalues,text);
			g_free(text);
			break;
		}
		break;
	case 2://crystal definition
		if(uspex_gui.calc._nlatticevalues==6){
			text=g_strdup_printf("%4f %4f %4f %4f %4f %4f",
				uspex_gui.calc.Latticevalues[0],uspex_gui.calc.Latticevalues[1],uspex_gui.calc.Latticevalues[2],
				uspex_gui.calc.Latticevalues[3],uspex_gui.calc.Latticevalues[4],uspex_gui.calc.Latticevalues[5]);
			GUI_COMBOBOX_ADD(uspex_gui.Latticevalues,text);
			g_free(text);
		}else if(uspex_gui.calc._nlatticevalues==4){
			text=g_strdup_printf("%4f %4f %4f %4f %4f %4f",
				uspex_gui.calc.Latticevalues[0],uspex_gui.calc.Latticevalues[1],0.,
				uspex_gui.calc.Latticevalues[3],uspex_gui.calc.Latticevalues[4],0.);
			GUI_COMBOBOX_ADD(uspex_gui.Latticevalues,text);
			g_free(text);
		}else{
			text=g_strdup_printf("%4f %4f %4f %4f %4f %4f",0.,0.,0.,0.,0.,0.);
			GUI_COMBOBOX_ADD(uspex_gui.Latticevalues,text);
			g_free(text);
		}
		break;
	case 0:
	default:
		/*ie. Volume(s)*/
		if(uspex_gui.calc._nlatticevalues==0){
			/*there was nothing but wants to add volume*/
			text=g_strdup_printf("%4f",0.);
			GUI_COMBOBOX_ADD(uspex_gui.Latticevalues,text);
			g_free(text);
		}else{
			gint i;
			gchar *tamp;
			tamp=g_strdup_printf("%4f",uspex_gui.calc.Latticevalues[0]);
			text=tamp;
			for (i=1;i<uspex_gui.calc._nlatticevalues;i++) {
				text=g_strdup_printf("%s %4f",tamp,uspex_gui.calc.Latticevalues[i]);
				g_free(tamp);tamp=text;
			}
			GUI_COMBOBOX_ADD(uspex_gui.Latticevalues,text);
			g_free(text);
		}
	}
	GUI_COMBOBOX_SET(uspex_gui.Latticevalues,0);
}
/***********************************/
/* Toggle spacegroup determination */
/***********************************/
void SG_toggle(void){
	if(uspex_gui.calc.doSpaceGroup){
		GUI_UNLOCK(uspex_gui.SymTolerance);
	}else{
		GUI_LOCK(uspex_gui.SymTolerance);
	}
}

/************************/
/* Switch notebook page */
/************************/
void uspex_gui_page_switch(GUI_NOTE *notebook,GUI_OBJ *page,guint page_num){
	gint from_page;
//	gchar *text;
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
//	FILE *save_fp;
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
//	GUI_OBJ *separator;
	/* special */
	struct model_pak *data;
//	struct core_pak *core;
//	gint idx;
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
		/*be careful to read Parameters.txt, specifically*/
		title = g_strdup_printf("%s%s",_UO.calc->path,"Parameters.txt");
		uspex_calc = read_uspex_parameters(title,data->num_species);
		g_free(title);
		if(uspex_calc){
			uspex_gui.have_output=TRUE;
			copy_uspex_parameters(uspex_calc,&(uspex_gui.calc));
			gui_uspex_init(data);
		}else{
			uspex_gui.have_output=FALSE;
			init_uspex_parameters(&(uspex_gui.calc));
			gui_uspex_init(data);
		}
		g_free(uspex_calc);
	} else {
		uspex_gui.have_output=FALSE;
		init_uspex_parameters(&(uspex_gui.calc));
		gui_uspex_init(data);
	}
/*bailout process*/

if(uspex_gui.calc.numSpecies==NULL){
	title=g_strdup("USPEX model is missing numSpecies information!\n");
	gui_text_show(ERROR,title);
	g_free(title);
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
	tmp=g_strdup_printf("%s%s",_UO.calc->path,"Parameters.txt");
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

/*------------------*/
/* page 1 -> SYSTEM */
/*------------------*/
	GUI_PAGE_NOTE(notebook,page,"SYSTEM");
/* create a table in the page*/
	GUI_TABLE_NOTE(page,table,17,4);
/* --- Type & System */
	GUI_LABEL_TABLE(table,"Type & System",0,4,0,1);
/* line 1 */
	GUI_COMBOBOX_TABLE(table,uspex_gui.calculationMethod,"Method: ",0,1,1,2);/*multiline*/
	GUI_COMBOBOX_ADD(uspex_gui.calculationMethod,"USPEX");
	GUI_COMBOBOX_ADD(uspex_gui.calculationMethod,"META");
	GUI_COMBOBOX_ADD(uspex_gui.calculationMethod,"VCNEB");
	GUI_COMBOBOX_ADD(uspex_gui.calculationMethod,"PSO");
	GUI_COMBOBOX_ADD(uspex_gui.calculationMethod,"TPS");
	GUI_COMBOBOX_ADD(uspex_gui.calculationMethod,"MINHOP");
GUI_TOOLTIP(uspex_gui.calculationMethod,"calculationMethod: Ch. 4.1 DEFAULT: USPEX\nSet the method of calculation.");
	GUI_COMBOBOX_TABLE(table,uspex_gui.calculationType,"Type: ",1,3,1,2);
	GUI_COMBOBOX_ADD(uspex_gui.calculationType,"300");
	GUI_COMBOBOX_ADD(uspex_gui.calculationType,"s300");
	GUI_COMBOBOX_ADD(uspex_gui.calculationType,"301");
	GUI_COMBOBOX_ADD(uspex_gui.calculationType,"s301");
	GUI_COMBOBOX_ADD(uspex_gui.calculationType,"310");
	GUI_COMBOBOX_ADD(uspex_gui.calculationType,"311");
	GUI_COMBOBOX_ADD(uspex_gui.calculationType,"000");
	GUI_COMBOBOX_ADD(uspex_gui.calculationType,"s000");
	GUI_COMBOBOX_ADD(uspex_gui.calculationType,"110");
	GUI_COMBOBOX_ADD(uspex_gui.calculationType,"200");
	GUI_COMBOBOX_ADD(uspex_gui.calculationType,"s200");
	GUI_COMBOBOX_ADD(uspex_gui.calculationType,"201");
	GUI_COMBOBOX_ADD(uspex_gui.calculationType,"s201");
	GUI_COMBOBOX_ADD(uspex_gui.calculationType,"-200");
	GUI_COMBOBOX_ADD(uspex_gui.calculationType,"-s200");
GUI_TOOLTIP(uspex_gui.calculationType,"calculationType: Ch. 4.1 DEFAULT: 300\nSets the dimensionality, molecularity, and variability\nof the calculation, can also be set individually.\n(311) and (110) calculations may not be supported.");
	GUI_CHECK_TABLE(table,uspex_gui._calctype_mag,uspex_gui.calc._calctype_mag,mag_toggle,"MAG",3,4,1,2);
GUI_TOOLTIP(uspex_gui._calctype_mag,"Set a magnetic calculation.");
/* line 2 */
	/*col 1: empty*/
	GUI_SPIN_TABLE(table,uspex_gui._calctype_dim,uspex_gui._dim,spin_update_dim,"DIM",1,2,2,3);
	GUI_SPIN_RANGE(uspex_gui._calctype_dim,-2.,3.);
GUI_TOOLTIP(uspex_gui._calctype_dim,"Set the dimension of the system.\n3, 2, 1, and 0 are equilvalent to 3D, 2D, 1D, and 0D.\n-2 correspond to 2D crystals.");
	GUI_CHECK_TABLE(table,uspex_gui._calctype_mol,uspex_gui.calc._calctype_mol,mol_toggle,"MOL",2,3,2,3);
GUI_TOOLTIP(uspex_gui._calctype_mol,"Set the molecularity of the system.");
	GUI_CHECK_TABLE(table,uspex_gui._calctype_var,uspex_gui.calc._calctype_var,var_toggle,"VAR",3,4,2,3);
GUI_TOOLTIP(uspex_gui._calctype_var,"Set the variability of chemical composition.");
/* line 3 */
	GUI_COMBOBOX_TABLE(table,uspex_gui.atomType,"atomType:",0,1,3,4);/*multiline*/
	GUI_COMBOBOX_ADD(uspex_gui.atomType,"ADD ATOMTYPE");
GUI_TOOLTIP(uspex_gui.atomType,"atomType: Ch. 4.1 DEFAULT: none\nThis list regroups several tags:\natomType, numSpecies, and valences.");
	GUI_ENTRY_TABLE(table,uspex_gui._atom_sym,uspex_gui._tmp_atom_sym,"%s","@Sym:",1,2,3,4);
GUI_TOOLTIP(uspex_gui._atom_sym,"atomSym: - DEFAULT: none\nAtomic symbol of current species.");
	GUI_ENTRY_TABLE(table,uspex_gui._atom_typ,uspex_gui._tmp_atom_typ,"%3i","@Z:",2,3,3,4);
GUI_TOOLTIP(uspex_gui._atom_typ,"atomTyp: - DEFAULT: none\nAtomic number of current species.");
	GUI_APPLY_BUTTON_TABLE(table,button,apply_atom,3,4,3,4);
/* line 4 */
	/*col 1: empty*/
	GUI_ENTRY_TABLE(table,uspex_gui._atom_num,uspex_gui._tmp_atom_num,"%3i","@Num:",1,2,4,5);
GUI_TOOLTIP(uspex_gui._atom_num,"atomNum: - DEFAULT: none\nNumber of atoms in current species.");
	GUI_ENTRY_TABLE(table,uspex_gui._atom_val,uspex_gui._tmp_atom_val,"%3i","@Val:",2,3,4,5);
GUI_TOOLTIP(uspex_gui._atom_val,"atomVal: - DEFAULT: auto\nValence of the current species.\nAutomatically determined if zero.");
	GUI_DELETE_BUTTON_TABLE(table,button,remove_atom,3,4,4,5);
/* line 5 */
	GUI_COMBOBOX_TABLE(table,uspex_gui.goodBonds,"goodBonds:",0,1,5,6);/*multiline*/
	GUI_COMBOBOX_ADD(uspex_gui.goodBonds,"ADD GOODBOND");
GUI_TOOLTIP(uspex_gui.goodBonds,"goodBonds: Ch. 4.1 DEFAULT: auto\nSet the minimum distance at which a bond is considered.");
	GUI_TEXT_TABLE(table,uspex_gui._bond_d,uspex_gui._tmp_bond_d,"Bonds: ",1,3,5,6);
GUI_TOOLTIP(uspex_gui._bond_d,"Minimum bond distance between selected and others species.");
	GUI_APPLY_BUTTON_TABLE(table,button,apply_bonds,3,4,5,6);
/* line 6 */
	/*col 1: empty*/
	/*col 2: empty*/
	GUI_CHECK_TABLE(table,button,uspex_gui.auto_bonds,auto_bond_toggle,"AUTO_BONDS",2,3,6,7);
GUI_TOOLTIP(button,"Automatically determine bonds (recommended).");
	GUI_DELETE_BUTTON_TABLE(table,button,remove_bonds,3,4,6,7);
/* line 7 */
        GUI_COMBOBOX_TABLE(table,uspex_gui.optType,"optType:",0,1,7,8);/*multiline*/
        GUI_COMBOBOX_ADD(uspex_gui.optType,"1: MIN Enthalpy (stable phases)");
        GUI_COMBOBOX_ADD(uspex_gui.optType,"2: MIN Volume (densest structure)");
        GUI_COMBOBOX_ADD(uspex_gui.optType,"3: MAX Hardness (hardest phase)");
        GUI_COMBOBOX_ADD(uspex_gui.optType,"4: MAX order (most order structure)");
        GUI_COMBOBOX_ADD(uspex_gui.optType,"5: MAX Structure average difference");
        GUI_COMBOBOX_ADD(uspex_gui.optType,"6: MAX Dielectric susceptibility");
        GUI_COMBOBOX_ADD(uspex_gui.optType,"7: MAX Band gap");
        GUI_COMBOBOX_ADD(uspex_gui.optType,"8: MAX electric energy storage capacity");
        GUI_COMBOBOX_ADD(uspex_gui.optType,"9: MAX Magnetization");
        GUI_COMBOBOX_ADD(uspex_gui.optType,"10: MAX Structure quasientropy");
	GUI_COMBOBOX_ADD(uspex_gui.optType,"11: MAX? L/H eigenvalue difference of refractive index");
	GUI_COMBOBOX_ADD(uspex_gui.optType,"14: MAX? ZT thermoelectric figure of merit");
	GUI_COMBOBOX_ADD(uspex_gui.optType,"17: MAX? free energy at finite temperature");
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
	GUI_TEXT_TABLE(table,uspex_gui.new_optType,uspex_gui._tmp_new_optType,"NEW optType: ",1,3,7,8);
GUI_TOOLTIP(uspex_gui.optType,"optType: Ch. 4.1 DEFAULT: MIN_enthalpy\nSelect the properties to optimize in the new\nVER 10.1 USPEX format.");
	GUI_CHECK_TABLE(table,button,uspex_gui.have_new_opt,opt_toggle,"NEW",3,4,7,8);
GUI_TOOLTIP(button,"Use the VER 10.1 optType format.\nMandatory for multiobjective optimization.");
/* line 8 */
	/*col 1: empty*/
	GUI_CHECK_TABLE(table,button,uspex_gui.calc.anti_opt,NULL,"ANTI-OPT",1,2,8,9);/*not calling anything*/
GUI_TOOLTIP(button,"anti-opt: - DEFAULT: FALSE\nIf set REVERSE the direction of optimization.\ni.e. MIN -> MAX & MAX -> MIN");
	GUI_CHECK_TABLE(table,button,uspex_gui.calc.checkMolecules,NULL,"ckMol",2,3,8,9);/*not calling anything*/
GUI_TOOLTIP(button,"checkMolecules: Ch. 4.1 DEFAULT: TRUE\nCheck and discard broken/merged molecules.");
	GUI_CHECK_TABLE(table,button,uspex_gui.calc.checkConnectivity,NULL,"ckCon",3,4,8,9);/*not calling anything*/
GUI_TOOLTIP(button,"checkConnectivity: Ch. 4.1 DEFAULT: FALSE\nCalculate hardness and add connectivity in softmutation.");
/* <- Ext Pressure, LDA+U, and magRation moved to p. III, II, and III, respectively.*/
/* initialize */
	GUI_COMBOBOX_SETUP(uspex_gui.calculationMethod,0,uspex_method_selected);
	GUI_COMBOBOX_SETUP(uspex_gui.calculationType,0,uspex_type_selected);
	GUI_COMBOBOX_SETUP(uspex_gui.optType,0,uspex_optimization_selected);
	GUI_LOCK(uspex_gui._atom_typ);
	populate_atomType();
	GUI_COMBOBOX_SETUP(uspex_gui.atomType,uspex_gui.calc._nspecies,atomType_selected);
	GUI_SPIN_SET(uspex_gui._calctype_dim,(gdouble) uspex_gui.calc._calctype_dim);
	mol_toggle();
	var_toggle();
	mag_toggle();
	GUI_COMBOBOX_SETUP(uspex_gui.goodBonds,0,goodBonds_selected);
	auto_bond_toggle();
	opt_toggle();
/* --- cell */
	GUI_LABEL_TABLE(table,"Cell",0,4,9,10);
/* line 10 */
        GUI_COMBOBOX_TABLE(table,uspex_gui.Latticevalues,"Lattice:",0,1,10,11);/*multiline*/
GUI_TOOLTIP(uspex_gui.Latticevalues,"Latticevalues: Ch. 4.6 DEFAULT: auto\nInitial volume of the unit cell, or known lattice parameters.");
        GUI_TEXT_TABLE(table,uspex_gui._latticevalue,uspex_gui._tmp_latticevalue,"VALUES:",1,3,10,11);
GUI_TOOLTIP(uspex_gui._latticevalue,"VALUES: values of the current line of Latticevalues.");
        GUI_APPLY_BUTTON_TABLE(table,button,apply_latticevalue,3,4,10,11);
/* line 11 */
        /*col 1: empty*/
        GUI_COMBOBOX_TABLE(table,uspex_gui._latticeformat,"FORMAT:",1,2,11,12);
        GUI_COMBOBOX_ADD(uspex_gui._latticeformat,"Volumes");
        GUI_COMBOBOX_ADD(uspex_gui._latticeformat,"Lattice");
        GUI_COMBOBOX_ADD(uspex_gui._latticeformat,"Crystal");
GUI_TOOLTIP(uspex_gui._latticeformat,"FORMAT: whether Lattice values correspond to a series of\nVolumes, lattive vectors, or crystallographic definition.");
        GUI_TEXT_TABLE(table,uspex_gui.splitInto,uspex_gui._tmp_splitInto,"split:",2,3,11,12);
GUI_TOOLTIP(uspex_gui.splitInto,"splitInto: Ch. 4.5 DEFAULT: 1\nNumber of identical subcells or pseudosubcells in the unitcell.");
        GUI_CHECK_TABLE(table,button,uspex_gui.auto_lval,toggle_auto_lval,"AUTO_LAT",3,4,11,12);
GUI_TOOLTIP(button,"AUTO_LAT: use automatic value for Latticevalues.");
/* Constraints */
	GUI_LABEL_TABLE(table,"Constraints",0,4,12,13);
/* line 13 */
	GUI_COMBOBOX_TABLE(table,uspex_gui.IonDistances,"Ion:",0,1,13,14);/*multiline*/
GUI_TOOLTIP(uspex_gui.IonDistances,"IonDistance: Ch. 4.5 DEFAULT: auto\nTriangular matrix of the minimum allowed\ninteratomic distances.");
	GUI_TEXT_TABLE(table,uspex_gui._distances,uspex_gui._tmp_distances,"DIST: ",1,3,13,14);
GUI_TOOLTIP(uspex_gui._distances,"distances: minimum allowed distance from the current atom to others.");
	GUI_APPLY_BUTTON_TABLE(table,button,apply_distances,3,4,13,14);
/* line 14 */
	/*col 1: empty*/
	GUI_ENTRY_TABLE(table,uspex_gui.minVectorLength,uspex_gui.calc.minVectorLength,"%.4f","MV:",1,2,14,15);
GUI_TOOLTIP(uspex_gui.minVectorLength,"minVectorLength: Ch. 4.5 DEFAULT: auto\nMinimum length of a new lattice parameter.");
	GUI_ENTRY_TABLE(table,uspex_gui.constraint_enhancement,uspex_gui.calc.constraint_enhancement,"%i","CE:",2,3,14,15);
GUI_TOOLTIP(uspex_gui.constraint_enhancement,"constraint_enhancement: Ch. 4.5 DEFAULT: 1\nApply CE times the IonDistance constraints.");
	GUI_CHECK_TABLE(table,button,uspex_gui.auto_C_ion,toggle_auto_C_ion,"AUTO_ION",3,4,14,15);
GUI_TOOLTIP(button,"AUTO_ION: use automatic values for ion constraints.");
/* line 15 */
	GUI_COMBOBOX_TABLE(table,uspex_gui.MolCenters,"Mol:",0,1,15,16);/*multiline*/
GUI_TOOLTIP(uspex_gui.MolCenters,"MolCenters: Ch. 4.5 DEFAULT: none\nTriangular matrix of the minimum allowed\ndistances between molecules centers.");
	GUI_TEXT_TABLE(table,uspex_gui._centers,uspex_gui._tmp_centers,"CENTER:",1,3,15,16);
GUI_TOOLTIP(uspex_gui._centers,"centers: minimum allowed distance from the current molecule center to others.");
	GUI_APPLY_BUTTON_TABLE(table,button,apply_centers,3,4,15,16);
/* line 16 */
	/*col 1: empty*/
	GUI_COMBOBOX_TABLE(table,uspex_gui._molModels,"Model:",1,3,16,17);/*NEW*/
	GUI_COMBOBOX_ADD(uspex_gui._molModels,"UNDER CONSTRUCTION");
GUI_TOOLTIP(uspex_gui._centers,"Associate each molecule with a GDIS model.");
	/*col 4: empty*/
/* initialize */
	GUI_COMBOBOX_SETUP(uspex_gui.IonDistances,0,uspex_IonDistances_selected);
	GUI_COMBOBOX_SETUP(uspex_gui.MolCenters,0,uspex_MolCenters_selected);
	refresh_constraints();
	toggle_auto_C_ion();
	GUI_COMBOBOX_SETUP(uspex_gui._latticeformat,1,uspex_latticeformat_selected);
	uspex_latticeformat_selected(uspex_gui._latticeformat);
/* --- end frame */

/*---------------------*/
/* page 2 -> STRUCTURE */
/*---------------------*/
	GUI_PAGE_NOTE(notebook,page,"STRUCTURE");
/* create a table in the page*/
        GUI_TABLE_NOTE(page,table,18,4);
/* --- Population & selection */
        GUI_LABEL_TABLE(table,"Population & selection",0,4,0,1);
/* line 1 */
	GUI_ENTRY_TABLE(table,uspex_gui.populationSize,uspex_gui.calc.populationSize,"%4i","SIZE:",0,1,1,2);
GUI_TOOLTIP(uspex_gui.populationSize,"populationSize: Ch. 4.2 DEFAULT: auto\nNumber of structures in each generation.");
	GUI_ENTRY_TABLE(table,uspex_gui.initialPopSize,uspex_gui.calc.initialPopSize,"%4i","INIT:",1,2,1,2);
GUI_TOOLTIP(uspex_gui.initialPopSize,"initialPopSize: Ch. 4.2 DEFAULT: populationSize\nNumber of structures in initial generation.");
	GUI_ENTRY_TABLE(table,uspex_gui.numGenerations,uspex_gui.calc.numGenerations,"%4i","NGEN:",2,3,1,2);
GUI_TOOLTIP(uspex_gui.numGenerations,"numGenerations: Ch. 4.2 DEFAULT: 100\nMaximum number of generations.");
	GUI_ENTRY_TABLE(table,uspex_gui.stopCrit,uspex_gui.calc.stopCrit,"%4i","STOP:",3,4,1,2);
GUI_TOOLTIP(uspex_gui.stopCrit,"stopCrit: Ch. 4.2 DEFAULT: auto\nMaximum number of generations.");
/* line 2 */
	GUI_ENTRY_TABLE(table,uspex_gui.mag_nm,uspex_gui.calc.magRatio[0],"%.4f","N.M.:",0,1,2,3);
GUI_TOOLTIP(uspex_gui.mag_nm,"magRatio: Ch. 4.1 DEFAULT: 0.1\n(magnetic calculation) Initial ratio of structures\nwith a non-magnetic order.");
	GUI_ENTRY_TABLE(table,uspex_gui.mag_fmls,uspex_gui.calc.magRatio[1],"%.4f","FM-LS:",1,2,2,3);
GUI_TOOLTIP(uspex_gui.mag_fmls,"magRatio: Ch. 4.1 DEFAULT: 0.1\n(magnetic calculation) Initial ratio of structures\nwith a low-spin ferromagnetic order.");
	GUI_ENTRY_TABLE(table,uspex_gui.mag_afml,uspex_gui.calc.magRatio[3],"%.4f","AFM-L:",2,3,2,3);
GUI_TOOLTIP(uspex_gui.mag_afml,"magRatio: Ch. 4.1 DEFAULT: 0.1\n(magnetic calculation) Initial ratio of structures\nwith a low-spin antiferromagnetic order.");
	GUI_ENTRY_TABLE(table,uspex_gui.mag_fmlh,uspex_gui.calc.magRatio[5],"%.4f","FM-LH:",3,4,2,3);
GUI_TOOLTIP(uspex_gui.mag_fmlh,"magRatio: Ch. 4.1 DEFAULT: 0.1\n(magnetic calculation) Initial ratio of structures\nwith a low/high spin mixed ferromagnetic order.");
/* line 3 */
	GUI_CHECK_TABLE(table,button,uspex_gui.calc._calctype_mag,mag_toggle,"MAG",0,1,3,4);
GUI_TOOLTIP(button,"Set a magnetic calculation.");
	GUI_ENTRY_TABLE(table,uspex_gui.mag_fmhs,uspex_gui.calc.magRatio[2],"%.4f","FM-HS:",1,2,3,4);
GUI_TOOLTIP(uspex_gui.mag_fmhs,"magRatio: Ch. 4.1 DEFAULT: 0.1\n(magnetic calculation) Initial ratio of structures\nwith a high-spin ferromagnetic.");
	GUI_ENTRY_TABLE(table,uspex_gui.mag_afmh,uspex_gui.calc.magRatio[4],"%.4f","AFM-H:",2,3,3,4);
GUI_TOOLTIP(uspex_gui.mag_afmh,"magRatio: Ch. 4.1 DEFAULT: 0.1\n(magnetic calculation) Initial ratio of structures\nwith a high-spin antiferromagnetic.");
	GUI_ENTRY_TABLE(table,uspex_gui.mag_aflh,uspex_gui.calc.magRatio[6],"%.4f","AF-LH:",3,4,3,4);
GUI_TOOLTIP(uspex_gui.mag_aflh,"magRatio: Ch. 4.1 DEFAULT: 0.1\n(magnetic calculation) Initial ratio of structures\nwith a low/high spin mixed antiferromagnetic.");
/* line 4 */
	GUI_ENTRY_TABLE(table,uspex_gui.bestFrac,uspex_gui.calc.bestFrac,"%.4f","Best:",0,1,4,5);
GUI_TOOLTIP(uspex_gui.bestFrac,"bestFrac: Ch. 4.3 DEFAULT: 0.7\nFraction of current generation used to generate the next.");
	GUI_ENTRY_TABLE(table,uspex_gui.keepBestHM,uspex_gui.calc.keepBestHM,"%3i","BestHM:",1,2,4,5);
GUI_TOOLTIP(uspex_gui.keepBestHM,"keepBestHM: Ch. 4.3 DEFAULT: auto\nNumber of best structures that will survive in next generation.");
	GUI_CHECK_TABLE(table,button,uspex_gui.calc.reoptOld,NULL,"reopt",2,3,4,5);/*not calling anything*/
GUI_TOOLTIP(button,"reoptOld: Ch. 4.3 DEFAULT: FALSE\nIf set surviving structure will be re-optimized.");
	GUI_ENTRY_TABLE(table,uspex_gui.fitLimit,uspex_gui.calc.fitLimit,"%.4f","BestFIT:",3,4,4,5);
GUI_TOOLTIP(uspex_gui.fitLimit,"fitLimit: Ch. 4.3 DEFAULT: none\nStop calculation when fitLimit fitness is reached.");
/* --- Structure & Variation */
	GUI_LABEL_TABLE(table,"Structure & Variation",0,4,5,6);
/* line 6 */
	GUI_TEXT_TABLE(table,uspex_gui.symmetries,uspex_gui.calc.symmetries,"symmetries: ",0,4,6,7);
GUI_TOOLTIP(uspex_gui.symmetries,"symmetries: Ch. 4.4 DEFAULT: auto\nPossible space group for crystals,\nplane group for 2D crystals/surface\nor point group for clusters.");
/* line 7 */
	GUI_ENTRY_TABLE(table,uspex_gui.fracGene,uspex_gui.calc.fracGene,"%.4f","Heredity:",0,1,7,8);
GUI_TOOLTIP(uspex_gui.fracGene,"fracGene: Ch. 4.4 DEFAULT: 0.5\nRatio of structures obtained by heredity.");
	GUI_ENTRY_TABLE(table,uspex_gui.fracRand,uspex_gui.calc.fracRand,"%.4f","Random:",1,2,7,8);
GUI_TOOLTIP(uspex_gui.fracRand,"fracRand: Ch. 4.4 DEFAULT: 0.2\nRatio of structures obtained randomly.");
	GUI_ENTRY_TABLE(table,uspex_gui.fracTopRand,uspex_gui.calc.fracTopRand,"%.4f","TOPRand:",2,3,7,8);
GUI_TOOLTIP(uspex_gui.fracTopRand,"fracTopRand: Ch. 4.4 DEFAULT: 0.2\nRatio of structures obtained by topological random generator.");
	GUI_ENTRY_TABLE(table,uspex_gui.fracPerm,uspex_gui.calc.fracPerm,"%.4f","Permutation:",3,4,7,8);
GUI_TOOLTIP(uspex_gui.fracPerm,"fracPerm: Ch. 4.4 DEFAULT: auto\nRatio of structures obtained by permutation.");
/* line 8 */
	GUI_ENTRY_TABLE(table,uspex_gui.fracAtomsMut,uspex_gui.calc.fracAtomsMut,"%.4f","AtmMut:",0,1,8,9);
GUI_TOOLTIP(uspex_gui.fracAtomsMut,"fracAtomsMut: Ch. 4.4 DEFAULT: 0.1\nRatio of structures obtained by softmutation.");
	GUI_ENTRY_TABLE(table,uspex_gui.fracRotMut,uspex_gui.calc.fracRotMut,"%.4f","RotMut:",1,2,8,9);
GUI_TOOLTIP(uspex_gui.fracRotMut,"fracRotMut: Ch. 4.4 DEFAULT: auto\nRatio of structures obtained by mutation of molecular orientation.");
	GUI_ENTRY_TABLE(table,uspex_gui.fracLatMut,uspex_gui.calc.fracLatMut,"%.4f","LatMut:",2,3,8,9);
GUI_TOOLTIP(uspex_gui.fracLatMut,"fracLatMut: Ch. 4.4 DEFAULT: auto\nRatio of structures obtained by lattice mutation.");
	GUI_ENTRY_TABLE(table,uspex_gui.fracSpinMut,uspex_gui.calc.fracSpinMut,"%.4f","SpinMut:",3,4,8,9);
GUI_TOOLTIP(uspex_gui.fracSpinMut,"fracSpinMut: Ch. 4.4 DEFAULT: 0.1\nRatio of structures obtained by spin mutation.");
/* line 9 */
	GUI_ENTRY_TABLE(table,uspex_gui.howManySwaps,uspex_gui.calc.howManySwaps,"%3i","NSwaps:",0,1,9,10);
GUI_TOOLTIP(uspex_gui.howManySwaps,"howManySwaps: Ch. 4.4 DEFAULT: auto\nNumber of pairwise swaps for permutation\ndistributed uniformly between [1,howManySwaps].");
	GUI_TEXT_TABLE(table,uspex_gui.specificSwaps,uspex_gui.calc.specificSwaps,"Swaps: ",1,4,9,10);
GUI_TOOLTIP(uspex_gui.specificSwaps,"specificSwaps: Ch. 4.4 DEFAULT: blank\nWich atoms are allow to swap during permutation.");
/* line 10 */
	GUI_ENTRY_TABLE(table,uspex_gui.mutationDegree,uspex_gui.calc.mutationDegree,"%.4f","mutationDegree:",0,1,10,11);
GUI_TOOLTIP(uspex_gui.mutationDegree,"mutationDegree: Ch. 4.4 DEFAULT: auto\nMaximum displacement in softmutation (Ang).");
	GUI_ENTRY_TABLE(table,uspex_gui.mutationRate,uspex_gui.calc.mutationRate,"%.4f","mutationRate:",1,2,10,11);
GUI_TOOLTIP(uspex_gui.mutationRate,"mutationRate: Ch. 4.4 DEFAULT: 0.5\nStd. dev. of epsilon in strain matrix for lattice mutation.");
	GUI_ENTRY_TABLE(table,uspex_gui.DisplaceInLatmutation,uspex_gui.calc.DisplaceInLatmutation,"%.4f","D_LATMUT:",2,3,10,11);
GUI_TOOLTIP(uspex_gui.DisplaceInLatmutation,"DisplaceInLatmutation: Ch. 4.4 DEFAULT: 1.0\nSets softmutation as part of lattice mutation\nand gives maximum displacement (Ang).");
	GUI_CHECK_TABLE(table,uspex_gui.AutoFrac,uspex_gui.calc.AutoFrac,NULL,"AutoFrac",3,4,10,11);/*not calling anything*/
GUI_TOOLTIP(uspex_gui.AutoFrac,"AutoFrac: Ch. 4.4 DEFAULT: FALSE\nIf set variation parameters will be optimized during run.");
/* --- Fingerprint & Antiseeds */
	GUI_LABEL_TABLE(table,"Fingerprint, antiseed, & spacegroup",0,4,11,12);
/* line 12 */
	GUI_LABEL_TABLE(table,"Fingerprint",0,1,12,13);
	GUI_ENTRY_TABLE(table,uspex_gui.RmaxFing,uspex_gui.calc.RmaxFing,"%.4f","RMax:",1,2,12,13);
GUI_TOOLTIP(uspex_gui.RmaxFing,"RmaxFing: Ch. 4.9 DEFAULT: 10.0\ncutoff distance (Ang.).");
	GUI_ENTRY_TABLE(table,uspex_gui.deltaFing,uspex_gui.calc.deltaFing,"%.4f","delta:",2,3,12,13);
GUI_TOOLTIP(uspex_gui.deltaFing,"deltaFing: Ch. 4.9 DEFAULT: 0.08\nDiscretization of the fingerprint function.");
	GUI_ENTRY_TABLE(table,uspex_gui.sigmaFing,uspex_gui.calc.sigmaFing,"%.4f","sigma:",3,4,12,13);
GUI_TOOLTIP(uspex_gui.sigmaFing,"sigmaFing: Ch. 4.9 DEFAULT: 0.03\nGaussian broadening of interatomic distance.");
/* line 13 */
	GUI_LABEL_TABLE(table,"Antiseed",0,1,13,14);
	GUI_ENTRY_TABLE(table,uspex_gui.antiSeedsActivation,uspex_gui.calc.antiSeedsActivation,"%4i","Activation:",1,2,13,14);
GUI_TOOLTIP(uspex_gui.antiSeedsActivation,"antiSeedsActivation: Ch. 4.10 DEFAULT: 5000\nGeneration at which antiseed is active.");
	GUI_ENTRY_TABLE(table,uspex_gui.antiSeedsMax,uspex_gui.calc.antiSeedsMax,"%.4f","Max:",2,3,13,14);
GUI_TOOLTIP(uspex_gui.antiSeedsMax,"antiSeedsMax: Ch. 4.10 DEFAULT: 0.000\nGaussian height in mean square deviation of the generation enthalpy\nAmong bestFrac structures, recommended value is 0.01.");
	GUI_ENTRY_TABLE(table,uspex_gui.antiSeedsSigma,uspex_gui.calc.antiSeedsSigma,"%.4f","sigma:",3,4,13,14);
GUI_TOOLTIP(uspex_gui.antiSeedsSigma,"antiSeedsSigma: Ch. 4.10 DEFAULT: 0.001\nGaussian width in average distances between generated structure\nAmong bestFrac structures, recommended calue is 0.005.");
/* line 14 */
	GUI_LABEL_TABLE(table,"Space group",0,1,14,15);
	/*col 2: empty*/
	GUI_CHECK_TABLE(table,uspex_gui.doSpaceGroup,uspex_gui.calc.doSpaceGroup,SG_toggle,"Active",2,3,14,15);
GUI_TOOLTIP(uspex_gui.doSpaceGroup,"doSpaceGroup: Ch. 4.11 DEFAULT: TRUE\nActivate space group determination.");
	GUI_ENTRY_TABLE(table,uspex_gui.SymTolerance,uspex_gui.calc.SymTolerance,"%.4f","TOL:",3,4,14,15);
GUI_TOOLTIP(uspex_gui.SymTolerance,"SymTolerance: Ch. 4.11 DEFAULT: 0.10\nPrecision for symmetry determination.");
/* initialize */
	SG_toggle();
/* --- Variable-composition */
	GUI_LABEL_TABLE(table,"Variable-composition",0,4,15,16);
/* line 16 */
	GUI_ENTRY_TABLE(table,uspex_gui.firstGeneMax,uspex_gui.calc.firstGeneMax,"%4i","1st_Gen:",0,1,16,17);
GUI_TOOLTIP(uspex_gui.firstGeneMax,"firstGeneMax: Ch. 5.5 DEFAULT: 11\nNumber of composition sampled in 1st generation.");
	GUI_ENTRY_TABLE(table,uspex_gui.minAt,uspex_gui.calc.minAt,"%4i","min@:",1,2,16,17);
GUI_TOOLTIP(uspex_gui.minAt,"minAt: Ch. 5.5 DEFAULT: none\nMinimum number of atoms (molecules) in unitcell\nfor the first generation.");
	GUI_ENTRY_TABLE(table,uspex_gui.maxAt,uspex_gui.calc.maxAt,"%4i","max@:",2,3,16,17);
GUI_TOOLTIP(uspex_gui.maxAt,"maxAt: Ch. 5.5 DEFAULT: none\nMaximum number of atoms (molecules) in unitcell\nfor the first generation.");
/*col 4: empty*/
/* line 17 */
	GUI_ENTRY_TABLE(table,uspex_gui.fracTrans,uspex_gui.calc.fracTrans,"%.4f","Trans:",0,1,17,18);
GUI_TOOLTIP(uspex_gui.fracTrans,"fracTrans: Ch. 5.5 DEFAULT: 0.1\nratio of structures obtained by transmutation.");
	GUI_ENTRY_TABLE(table,uspex_gui.howManyTrans,uspex_gui.calc.howManyTrans,"%.4f","NTrans:",1,2,17,18);
GUI_TOOLTIP(uspex_gui.howManyTrans,"howManyTrans: Ch. 5.5 DEFAULT: 0.2\nMaximum ratio of transmutated atoms in a structure.");
	GUI_ENTRY_TABLE(table,uspex_gui.specificTrans,uspex_gui._tmp_specificTrans,"%s","Specific:",2,4,17,18);
GUI_TOOLTIP(uspex_gui.specificTrans,"specificTrans: Ch. 5.5 DEFAULT: blank\nList of allowed transmutation.");

/*-----------------------*/
/* page 3 -> CALCULATION */
/*-----------------------*/
	GUI_PAGE_NOTE(notebook,page,"CALCULATION");
/* --- end frame */

/*--------------------*/
/* page 4 -> ADVANCED */
/*--------------------*/
	GUI_PAGE_NOTE(notebook,page,"ADVANCED");
/* --- end frame */

/*--------------------*/
/* page 5 -> SPECIFIC */
/*--------------------*/
	GUI_PAGE_NOTE(notebook,page,"SPECIFIC");
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

