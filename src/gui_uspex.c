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
#define GLIB_DISABLE_DEPRECATION_WARNINGS
#include <glib/gstdio.h>
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
#include "track.h"

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
	gchar *ptr;
	uspex_gui.cur_page=USPEX_PAGE_SYSTEM;
	/*prepare local values*/
	uspex_gui._tmp_atom_typ=0;uspex_gui._tmp_atom_num=0;uspex_gui._tmp_atom_val=0;
	strcpy(uspex_gui._tmp_atom_sym,elements[uspex_gui._tmp_atom_typ].symbol);
	if(uspex_gui.calc._nspecies<1) uspex_gui.calc._nspecies=1;
	if(uspex_gui.calc._num_opt_steps<1) uspex_gui.calc._num_opt_steps=1;
	uspex_gui._tmp_num_opt_steps=(gdouble)uspex_gui.calc._num_opt_steps;
	uspex_gui._tmp_curr_step=1.;
	if(uspex_gui.calc._isfixed==NULL) {
		uspex_gui.calc._isfixed=g_malloc0(uspex_gui.calc._num_opt_steps*sizeof(gboolean));/*malloc0 FIX*/
		for(idx=0;idx<uspex_gui.calc._num_opt_steps;idx++) uspex_gui.calc._isfixed[idx]=FALSE;
	}
	uspex_gui._tmp_isfixed=uspex_gui.calc._isfixed[0];
	if(uspex_gui.calc.abinitioCode==NULL) {
		uspex_gui.calc.abinitioCode=g_malloc(1*sizeof(gint));
		uspex_gui.calc.abinitioCode[0]=1;
	}
	/*get atom information if not available*/
	if(uspex_gui.calc.atomType==NULL){
		GSList *list;
		struct core_pak *core;
		/*first get the _nspecies*/
		list=find_unique(ELEMENT,model);
		uspex_gui.calc._nspecies=(gint)g_slist_length(list);
if(uspex_gui.calc._nspecies<1) {
	/*_nspecies is 0 which can happen in case of new_model*/
	uspex_gui.calc._nspecies=1;
	uspex_gui.calc.atomType = g_malloc0(uspex_gui.calc._nspecies*sizeof(gint));
}else{
		uspex_gui.calc.atomType = g_malloc(uspex_gui.calc._nspecies*sizeof(gint));
}
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
/*EX14 has valences, but no atomType \/(^^)\/ */
//		if(uspex_gui.calc.valences!=NULL) g_free(uspex_gui.calc.valences);
		if(uspex_gui.calc.valences==NULL){
			uspex_gui.calc.valences = g_malloc(uspex_gui.calc._nspecies*sizeof(gint));
			for(idx=0;idx<uspex_gui.calc._nspecies;idx++) uspex_gui.calc.valences[idx]=0;
		}
	}
	if(uspex_gui.calc.valences==NULL) uspex_gui.calc.valences = g_malloc0(uspex_gui.calc._nspecies*sizeof(gint));
	uspex_gui.have_ZT=FALSE;
	if(uspex_gui.calc.optType==US_OT_ZT){
		uspex_gui.have_ZT=TRUE;
	}else{
		if(uspex_gui.calc.new_optType!=NULL){
			if(find_in_string("ZT",uspex_gui.calc.new_optType) != NULL) uspex_gui.have_ZT=TRUE;
			if(find_in_string("14",uspex_gui.calc.new_optType) != NULL) uspex_gui.have_ZT=TRUE;
		}
	}
	/*we NEED to check optional valences*/
	if(uspex_gui.calc.valences==NULL){
		uspex_gui.calc.valences = g_malloc(uspex_gui.calc._nspecies*sizeof(gint));
		for(idx=0;idx<uspex_gui.calc._nspecies;idx++) uspex_gui.calc.valences[idx]=0;
	}
	uspex_gui.auto_bonds=(uspex_gui.calc.goodBonds==NULL);
	if(uspex_gui.calc._nmolecules<1) uspex_gui.calc._nmolecules=1;
	uspex_gui._tmp_num_mol=(gdouble)uspex_gui.calc._nmolecules;
	uspex_gui._tmp_curr_mol=1.;
	uspex_gui._tmp_mols_gdis=g_malloc(uspex_gui.calc._nmolecules*sizeof(gint));
	uspex_gui._tmp_mols_gulp=g_malloc(uspex_gui.calc._nmolecules*sizeof(gboolean));
	for(idx=0;idx<uspex_gui.calc._nmolecules;idx++){
		uspex_gui._tmp_mols_gdis[idx]=0;
		uspex_gui._tmp_mols_gulp[idx]=FALSE;
	}
	if(uspex_gui.calc._nsplits==0) {
		uspex_gui.calc._nsplits=1;
		uspex_gui.calc.splitInto=g_malloc(uspex_gui.calc._nsplits*sizeof(gint));
		uspex_gui.calc.splitInto[0]=1;
	}
	if(uspex_gui.calc._nspetrans<1) uspex_gui.calc._nspetrans=0;
	if(uspex_gui.calc._npickimg<1) uspex_gui.calc._npickimg=0;
	if(uspex_gui.calc.ldaU!=NULL){
		/*we have 1<n<_nspecies U values*/
		line=g_strdup("");
		for(idx=0;idx<uspex_gui.calc._nspecies;idx++) line=g_strdup_printf("%s %f",line,uspex_gui.calc.ldaU[idx]);
		uspex_gui._tmp_ldaU=line;
		line=NULL;
	}else{
		uspex_gui._tmp_ldaU=g_strdup("");
	}
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
			gint i,j;gdouble tmp_d=0.;/*FIX: b865cc*/
			gint n_swaps;
			gint *specificSwaps;
			gchar *ptr,*ptr2;
			/*populate specificSwaps*/
			/*1-count*/
			ptr=uspex_gui.calc.specificSwaps;n_swaps=0;
			while(*ptr==' ') ptr++;/*skip initial space (if any)*/
			ptr2=ptr;
			do{
				tmp_d+=g_ascii_strtod(ptr,&ptr2);/*FIX: b865cc*/
				if(ptr2==ptr) break;
				n_swaps++;
				ptr=ptr2;
			}while(1);
			tmp_d-=(gint)tmp_d;
			if((tmp_d>0.)||(n_swaps==0)) {
				/*erroneous data, give up specificSwaps information*/
				g_free(uspex_gui.calc.specificSwaps);
				uspex_gui.calc.specificSwaps=NULL;
			}else{
				/*2-populate specificSwaps*/
				specificSwaps = g_malloc(n_swaps*sizeof(gint));
				ptr=uspex_gui.calc.specificSwaps;
				while(*ptr==' ') ptr++;/*skip initial space (if any)*/
				ptr2=ptr;i=0;/*FIX: 736236*/
				do{
					tmp_d=g_ascii_strtod(ptr,&ptr2);
					if(ptr2==ptr) break;
					/*reg tmp_d*/
					specificSwaps[i] = (gint)tmp_d;
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
	}
	if(uspex_gui.calc.IonDistances==NULL) uspex_gui.auto_C_ion=TRUE;
	else uspex_gui.auto_C_ion=FALSE;
	if(uspex_gui.calc.Latticevalues==NULL){
		uspex_gui.calc._nlattice_line=0;
		uspex_gui.calc._nlattice_vals=0;
		uspex_gui.auto_C_lat=TRUE;
	}else{
		uspex_gui.auto_C_lat=FALSE;
	}
	if(uspex_gui.calc._calctype_dim==3) uspex_gui.calc.doSpaceGroup=TRUE;
	else uspex_gui.calc.doSpaceGroup=FALSE;
	/*calculation steps specific*/
	if(uspex_gui.have_specific){
		/*use job_path / Specific as a default.*/
		if(uspex_gui._tmp_spe_folder==NULL){
			uspex_gui._tmp_spe_folder=g_build_filename(uspex_gui.calc.path,"../Specific",NULL);
		}
	}
	if(uspex_gui.calc.KresolStart==NULL) {
		uspex_gui.calc.KresolStart=g_malloc(uspex_gui.calc._num_opt_steps*sizeof(gdouble));
		uspex_gui.calc.KresolStart[0]=0.2;
		for(idx=1;idx<uspex_gui.calc._num_opt_steps;idx++) 
			uspex_gui.calc.KresolStart[idx]=0.2-(gdouble)(idx)*(0.2-0.08)/(uspex_gui.calc._num_opt_steps-1);
	}
	if(uspex_gui.calc.vacuumSize==NULL) {
		uspex_gui.calc.vacuumSize=g_malloc(uspex_gui.calc._num_opt_steps*sizeof(gdouble));
		for(idx=0;idx<uspex_gui.calc._num_opt_steps;idx++) uspex_gui.calc.vacuumSize[idx]=10.;
	}
	uspex_gui._tmp_commandExecutable=g_malloc0(uspex_gui.calc._num_opt_steps*sizeof(gchar *));
//	for(idx=0;idx<uspex_gui.calc._num_opt_steps;idx++) uspex_gui._tmp_commandExecutable[idx]=NULL;
	if(uspex_gui.calc.commandExecutable!=NULL) {
		/*cut*/
		gint size;
		idx=0;/*for safety*/
		line=&(uspex_gui.calc.commandExecutable[0]);
		ptr=line;
		while(*line!='\0'){
			while(!g_ascii_isprint(*ptr)) ptr++;
			line=ptr;
			size=0;
			while(g_ascii_isprint(*ptr)) {
				size++;
				ptr++;
			}
			uspex_gui._tmp_commandExecutable[idx]=g_strndup(line,size);
			line=ptr;
			idx++;
		}
	}
	uspex_gui._tmp_ai_input=g_malloc(uspex_gui.calc._num_opt_steps*sizeof(gchar *));
	for(idx=0;idx<uspex_gui.calc._num_opt_steps;idx++) uspex_gui._tmp_ai_input[idx]=NULL;/*no default*/
	uspex_gui._tmp_ai_opt=g_malloc(uspex_gui.calc._num_opt_steps*sizeof(gchar *));
	for(idx=0;idx<uspex_gui.calc._num_opt_steps;idx++) uspex_gui._tmp_ai_opt[idx]=NULL;/*no default*/
	uspex_gui.auto_step=TRUE;/*auto step is ON by default*/
	uspex_gui._tmp_ai_lib_folder=NULL;
	uspex_gui._tmp_ai_lib_sel=NULL;
	uspex_gui._tmp_ai_lib_folder=g_malloc(uspex_gui.calc._num_opt_steps*sizeof(gchar *));
	uspex_gui._tmp_ai_lib_sel=g_malloc(uspex_gui.calc._num_opt_steps*sizeof(gchar *));
	for(idx=0;idx<uspex_gui.calc._num_opt_steps;idx++){
		uspex_gui._tmp_ai_lib_folder[idx]=NULL;
		uspex_gui._tmp_ai_lib_sel[idx]=NULL;
	}
	/**/
	uspex_gui.have_v1030=TRUE;
	uspex_gui.have_octave=TRUE;
	if(sysenv.uspex_path!=NULL) uspex_gui.calc.job_uspex_exe=g_strdup(sysenv.uspex_path);
	/*set everything for an update*/
	uspex_gui._tmp_nspecies=uspex_gui.calc._nspecies;/*BACKUP*/
	uspex_gui.is_dirty=TRUE;
}
/*****************************************************/
/*(re) populate atomType combo with atom information */
/*****************************************************/
void populate_atomType(void){
	gint idx,jdx;
	gint num;
	gchar *line;
	gint nspe,nval;
	gchar sym[3];
	/*this is call on new atomType information, so wipe previous information*/
	GUI_COMBOBOX_WIPE(uspex_gui.atomType);
	GUI_COMBOBOX_WIPE(uspex_gui.goodBonds);
/*there is a _BUG_ where numSpecies refers to number of molecules instead of number of elements!*/
if(uspex_gui.calc._calctype_mol) num=uspex_gui.calc._var_nspecies;
else num=uspex_gui.calc._nspecies;
	for(idx=0;idx<num;idx++){
		if(uspex_gui.calc.valences==NULL) nval=0;
		else nval=uspex_gui.calc.valences[idx];
		if(uspex_gui.calc.numSpecies==NULL) nspe=0;
		else nspe=uspex_gui.calc.numSpecies[idx];
		if(uspex_gui.calc.atomType==NULL) {
			sym[0]='X';
			sym[1]='X';
			sym[2]='\0';
		}else{
			sym[0]=elements[uspex_gui.calc.atomType[idx]].symbol[0];
			sym[1]=elements[uspex_gui.calc.atomType[idx]].symbol[1];
			sym[2]='\0';
		}
		line=g_strdup_printf("%s(%i) (V=%i)",sym,nspe,nval);
		GUI_COMBOBOX_ADD(uspex_gui.atomType,line);
		g_free(line);
	}
/*only update goodBonds if needed*/
if((uspex_gui.calc.goodBonds!=NULL)&&(uspex_gui.auto_bonds)){
	gchar *tmp;
	for(idx=0;idx<num;idx++){
		tmp=g_strdup("");
		line=tmp;
		for(jdx=0;jdx<uspex_gui.calc._nspecies;jdx++){
			line=g_strdup_printf("%s %5f",tmp,uspex_gui.calc.goodBonds[jdx+idx*uspex_gui.calc._nspecies]);
			g_free(tmp);tmp=line;
		}
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
			uspex_gui._tmp_nspecies=uspex_gui.calc._nspecies;/*BACKUP*/
			test_calc = read_uspex_parameters(filename,uspex_gui.calc._nspecies);
			if(test_calc == NULL) {
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
void load_uspex_exe(void){
	GUI_OBJ *file_chooser;
	gint have_answer;
	gchar *filename;
	gchar *text;
	/**/
	GUI_PREPARE_OPEN_DIALOG(uspex_gui.window,file_chooser,"Select USPEX Executable","*","uspex_exec");
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
/***********************************************/
/* check the possible format for latticevalues */
/***********************************************/
void check_latticevalues_format(void){
	gint index;
        if(uspex_gui.auto_C_lat) return;/*no need to care*/
	GUI_COMBOBOX_GET(uspex_gui._latticeformat,index);
        GUI_COMBOBOX_WIPE(uspex_gui._latticeformat);/*wipe in any case*/
        GUI_COMBOBOX_ADD(uspex_gui._latticeformat,"Volumes");/*this is always ok*/
        if((!uspex_gui.calc._calctype_var)&&((uspex_gui.calc._calctype_dim==3)||(uspex_gui.calc._calctype_dim==-2))){
                GUI_COMBOBOX_ADD(uspex_gui._latticeformat,"Lattice");
                GUI_COMBOBOX_ADD(uspex_gui._latticeformat,"Crystal");
        }else{
		if(index!=0) GUI_COMBOBOX_SET(uspex_gui._latticeformat,0);
	}
}
/***************************************************/
/* show/hide specific page + lock/unlock mechanism */
/***************************************************/
void update_specific(void){
	gint index;
	gboolean hide_specific;
	GUI_COMBOBOX_GET(uspex_gui.calculationMethod,index);
	/*lock by default*/
/*POPULATION*/
	GUI_LOCK(uspex_gui.populationSize);
	GUI_LOCK(uspex_gui.initialPopSize);
	GUI_LOCK(uspex_gui.numGenerations);
	GUI_LOCK(uspex_gui.stopCrit);
	GUI_LOCK(uspex_gui.mag_nm);
	GUI_LOCK(uspex_gui.mag_fmls);
	GUI_LOCK(uspex_gui.mag_fmhs);
	GUI_LOCK(uspex_gui.mag_afml);
	GUI_LOCK(uspex_gui.mag_afmh);
	GUI_LOCK(uspex_gui.mag_fmlh);
	GUI_LOCK(uspex_gui.mag_aflh);
/*SURVIVAL & SELECTION*/
	GUI_LOCK(uspex_gui.bestFrac);
	GUI_LOCK(uspex_gui.keepBestHM);
	GUI_LOCK(uspex_gui.reoptOld);
/*VARIATION OPERATORS*/
	GUI_LOCK(uspex_gui.symmetries);
	GUI_LOCK(uspex_gui.fracGene);
	GUI_LOCK(uspex_gui.fracRand);
	GUI_LOCK(uspex_gui.fracTopRand);
	GUI_LOCK(uspex_gui.fracPerm);
	GUI_LOCK(uspex_gui.fracAtomsMut);
	GUI_LOCK(uspex_gui.fracLatMut);
	GUI_LOCK(uspex_gui.fracSpinMut);
	GUI_LOCK(uspex_gui.howManySwaps);
	GUI_LOCK(uspex_gui.specificSwaps);
	GUI_LOCK(uspex_gui.mutationDegree);
	GUI_LOCK(uspex_gui.mutationRate);
	GUI_LOCK(uspex_gui.DisplaceInLatmutation);
	GUI_LOCK(uspex_gui.AutoFrac);
/*CELL*/
	GUI_LOCK(uspex_gui.Latticevalues);
	GUI_LOCK(uspex_gui._latticeformat);
	GUI_LOCK(uspex_gui._latticevalue);
	GUI_LOCK(uspex_gui.splitInto);
/*MOLECULAR*/
	GUI_LOCK(uspex_gui.checkConnectivity);
	GUI_LOCK(uspex_gui.checkMolecules);
	GUI_LOCK(uspex_gui.fracRotMut);
	GUI_LOCK(uspex_gui.MolCenters);
	GUI_LOCK(uspex_gui._centers);
	GUI_LOCK(uspex_gui._centers_button);
	GUI_LOCK(uspex_gui.mol_model);
	GUI_LOCK(uspex_gui.num_mol);
	GUI_LOCK(uspex_gui.mol_model_button);
	GUI_LOCK(uspex_gui.mol_gdis);
	GUI_LOCK(uspex_gui.curr_mol);
	GUI_LOCK(uspex_gui.mol_gulp);
	GUI_LOCK(uspex_gui.mol_apply_button);
/*BoltzTraP*/
	GUI_LOCK(uspex_gui.BoltzTraP_T_max);
	GUI_LOCK(uspex_gui.BoltzTraP_T_delta);
	GUI_LOCK(uspex_gui.BoltzTraP_T_efcut);
	GUI_LOCK(uspex_gui.TE_T_interest);
	GUI_LOCK(uspex_gui.TE_threshold);
	GUI_LOCK(uspex_gui.TE_goal);
	GUI_LOCK(uspex_gui.cmd_BoltzTraP);
	GUI_LOCK(uspex_gui.cmd_BoltzTraP_button);
/*CRYSTAL*/
	GUI_LOCK(uspex_gui.PhaseDiagram);
/*SURFACES*/
	GUI_LOCK(uspex_gui.thicknessS);
	GUI_LOCK(uspex_gui.thicknessB);
	GUI_LOCK(uspex_gui.reconstruct);
	GUI_LOCK(uspex_gui.StoichiometryStart);
	GUI_LOCK(uspex_gui.E_AB);
	GUI_LOCK(uspex_gui.Mu_A);
	GUI_LOCK(uspex_gui.Mu_B);
	GUI_LOCK(uspex_gui.substrate_model);
	GUI_LOCK(uspex_gui.substrate_model_button);
/*CLUSTERS*/
	GUI_LOCK(uspex_gui.numberparents);
/*VARCOMP*/
	GUI_LOCK(uspex_gui.numSpecies);
	GUI_LOCK(uspex_gui.blockSpecies);
	GUI_LOCK(uspex_gui.Species_apply_button);
	GUI_LOCK(uspex_gui.Species_delete_button);
	GUI_LOCK(uspex_gui.firstGeneMax);
	GUI_LOCK(uspex_gui.minAt);
	GUI_LOCK(uspex_gui.maxAt);
	GUI_LOCK(uspex_gui.fracTrans);
	GUI_LOCK(uspex_gui.howManyTrans);
	GUI_LOCK(uspex_gui.specificTrans);
/*META*/
	GUI_LOCK(uspex_gui.GaussianWidth);
	GUI_LOCK(uspex_gui.GaussianHeight);
	GUI_LOCK(uspex_gui.FullRelax);
	GUI_LOCK(uspex_gui.maxVectorLength);
	GUI_LOCK(uspex_gui.meta_model);
	GUI_LOCK(uspex_gui.meta_model_button);
/*VCNEB*/
	GUI_LOCK(uspex_gui.vcnebType);/* is *always* locked */
	GUI_LOCK(uspex_gui._vcnebtype_method);
	GUI_LOCK(uspex_gui._vcnebtype_img_num);
	GUI_LOCK(uspex_gui._vcnebtype_spring);
	GUI_LOCK(uspex_gui.numImages);
	GUI_LOCK(uspex_gui.numSteps);
	GUI_LOCK(uspex_gui.optReadImages);
	GUI_LOCK(uspex_gui.optimizerType);
	GUI_LOCK(uspex_gui.optRelaxType);
	GUI_LOCK(uspex_gui.dt);
	GUI_LOCK(uspex_gui.ConvThreshold);
	GUI_LOCK(uspex_gui.VarPathLength);
	GUI_LOCK(uspex_gui.K_min);
	GUI_LOCK(uspex_gui.K_max);
	GUI_LOCK(uspex_gui.Kconstant);
	GUI_LOCK(uspex_gui.optFreezing);
	GUI_LOCK(uspex_gui.optMethodCIDI);
	GUI_LOCK(uspex_gui.startCIDIStep);
	GUI_LOCK(uspex_gui.pickupImages);
	GUI_LOCK(uspex_gui.FormatType);
	GUI_LOCK(uspex_gui.PrintStep);
	GUI_LOCK(uspex_gui.img_model);
	GUI_LOCK(uspex_gui.img_model_button);
/*PSO*/
	GUI_LOCK(uspex_gui.PSO_softMut);
	GUI_LOCK(uspex_gui.PSO_BestStruc);
	GUI_LOCK(uspex_gui.PSO_BestEver);
/*TPS*/
	GUI_LOCK(uspex_gui.numIterations);
	GUI_LOCK(uspex_gui.speciesSymbol);
	GUI_LOCK(uspex_gui.mass);
	GUI_LOCK(uspex_gui.amplitudeShoot_AB);
	GUI_LOCK(uspex_gui.amplitudeShoot_BA);
	GUI_LOCK(uspex_gui.magnitudeShoot_success);
	GUI_LOCK(uspex_gui.magnitudeShoot_failure);
	GUI_LOCK(uspex_gui.shiftRatio);
	GUI_LOCK(uspex_gui.orderParaType);
	GUI_LOCK(uspex_gui.opCriteria_start);
	GUI_LOCK(uspex_gui.opCriteria_end);
	GUI_LOCK(uspex_gui.cmdOrderParameter);
	GUI_LOCK(uspex_gui.cmdOrderParameter_button);
	GUI_LOCK(uspex_gui.cmdEnthalpyTemperature);
	GUI_LOCK(uspex_gui.cmdEnthalpyTemperature_button);
	GUI_LOCK(uspex_gui.orderParameterFile);
	GUI_LOCK(uspex_gui.orderParameterFile_button);
	GUI_LOCK(uspex_gui.enthalpyTemperatureFile);
	GUI_LOCK(uspex_gui.enthalpyTemperatureFile_button);
	GUI_LOCK(uspex_gui.trajectoryFile);
	GUI_LOCK(uspex_gui.trajectoryFile_button);
	GUI_LOCK(uspex_gui.MDrestartFile);
	GUI_LOCK(uspex_gui.MDrestartFile_button);
	/*check method*/
	switch (index){
	case 0://USPEX
/*POPULATION*/
		GUI_UNLOCK(uspex_gui.populationSize);
		GUI_UNLOCK(uspex_gui.initialPopSize);
		GUI_UNLOCK(uspex_gui.numGenerations);
		GUI_UNLOCK(uspex_gui.stopCrit);
if(uspex_gui.calc._calctype_mag){
		GUI_UNLOCK(uspex_gui.mag_nm);
		GUI_UNLOCK(uspex_gui.mag_fmls);
		GUI_UNLOCK(uspex_gui.mag_fmhs);
		GUI_UNLOCK(uspex_gui.mag_afml);
		GUI_UNLOCK(uspex_gui.mag_afmh);
		GUI_UNLOCK(uspex_gui.mag_fmlh);
		GUI_UNLOCK(uspex_gui.mag_aflh);
}
/*SURVIVAL & SELECTION*/
		GUI_UNLOCK(uspex_gui.bestFrac);
		GUI_UNLOCK(uspex_gui.keepBestHM);
		GUI_UNLOCK(uspex_gui.reoptOld);
/*VARIATION OPERATORS*/
		GUI_UNLOCK(uspex_gui.symmetries);
		GUI_UNLOCK(uspex_gui.fracGene);
		GUI_UNLOCK(uspex_gui.fracRand);
		GUI_UNLOCK(uspex_gui.fracTopRand);
		GUI_UNLOCK(uspex_gui.fracPerm);
		GUI_UNLOCK(uspex_gui.fracAtomsMut);
		GUI_UNLOCK(uspex_gui.fracLatMut);
if(uspex_gui.calc._calctype_mag){
		GUI_UNLOCK(uspex_gui.fracSpinMut);
}
		GUI_UNLOCK(uspex_gui.howManySwaps);
		GUI_UNLOCK(uspex_gui.specificSwaps);
		GUI_UNLOCK(uspex_gui.mutationDegree);
		GUI_UNLOCK(uspex_gui.mutationRate);
		GUI_UNLOCK(uspex_gui.DisplaceInLatmutation);
		GUI_UNLOCK(uspex_gui.AutoFrac);
/*CELL*/
if(!uspex_gui.auto_C_lat){
		GUI_UNLOCK(uspex_gui.Latticevalues);
		GUI_UNLOCK(uspex_gui._latticeformat);
		GUI_UNLOCK(uspex_gui._latticevalue);
}
		GUI_UNLOCK(uspex_gui.splitInto);
/*MOLECULAR*/
if(uspex_gui.calc._calctype_mol){
		GUI_UNLOCK(uspex_gui.checkConnectivity);
		GUI_UNLOCK(uspex_gui.checkMolecules);
		GUI_UNLOCK(uspex_gui.fracRotMut);
		GUI_UNLOCK(uspex_gui.MolCenters);
		GUI_UNLOCK(uspex_gui._centers);
		GUI_UNLOCK(uspex_gui._centers_button);
		GUI_UNLOCK(uspex_gui.mol_model);
		GUI_LOCK(uspex_gui.num_mol);
		GUI_LOCK(uspex_gui.mol_model_button);
		GUI_LOCK(uspex_gui.mol_gdis);
		GUI_LOCK(uspex_gui.curr_mol);
		GUI_LOCK(uspex_gui.mol_gulp);
		GUI_LOCK(uspex_gui.mol_apply_button);
		GUI_COMBOBOX_SET(uspex_gui.mol_model,0);
}
/*BoltzTraP*/
if(uspex_gui.have_ZT){
		GUI_UNLOCK(uspex_gui.BoltzTraP_T_max);
		GUI_UNLOCK(uspex_gui.BoltzTraP_T_delta);
		GUI_UNLOCK(uspex_gui.BoltzTraP_T_efcut);
		GUI_UNLOCK(uspex_gui.TE_T_interest);
		GUI_UNLOCK(uspex_gui.TE_threshold);
		GUI_UNLOCK(uspex_gui.TE_goal);
		GUI_UNLOCK(uspex_gui.cmd_BoltzTraP);
		GUI_UNLOCK(uspex_gui.cmd_BoltzTraP_button);
}
/*CRYSTAL*/
if(uspex_gui.calc._calctype_dim==3){
		GUI_UNLOCK(uspex_gui.PhaseDiagram);
}
/*SURFACES*/
if(uspex_gui.calc._calctype_dim==2){
		GUI_UNLOCK(uspex_gui.thicknessS);
		GUI_UNLOCK(uspex_gui.thicknessB);
		GUI_UNLOCK(uspex_gui.reconstruct);
		GUI_UNLOCK(uspex_gui.StoichiometryStart);
		GUI_UNLOCK(uspex_gui.E_AB);
		GUI_UNLOCK(uspex_gui.Mu_A);
		GUI_UNLOCK(uspex_gui.Mu_B);
		GUI_UNLOCK(uspex_gui.substrate_model);
		GUI_UNLOCK(uspex_gui.substrate_model_button);
}
/*CLUSTERS*/
if(uspex_gui.calc._calctype_dim==0){
		GUI_UNLOCK(uspex_gui.numberparents);
}
/*VARCOMP*/
if(uspex_gui.calc._calctype_var){
		GUI_UNLOCK(uspex_gui.numSpecies);
		GUI_UNLOCK(uspex_gui.blockSpecies);
		GUI_UNLOCK(uspex_gui.Species_apply_button);
		GUI_UNLOCK(uspex_gui.Species_delete_button);
		GUI_UNLOCK(uspex_gui.firstGeneMax);
		GUI_UNLOCK(uspex_gui.minAt);
		GUI_UNLOCK(uspex_gui.maxAt);
		GUI_UNLOCK(uspex_gui.fracTrans);
		GUI_UNLOCK(uspex_gui.howManyTrans);
		GUI_UNLOCK(uspex_gui.specificTrans);
}
		break;
	case 1://META
/*POPULATION*/
		GUI_UNLOCK(uspex_gui.populationSize);
		GUI_UNLOCK(uspex_gui.initialPopSize);
		GUI_UNLOCK(uspex_gui.numGenerations);
		GUI_UNLOCK(uspex_gui.stopCrit);
if(uspex_gui.calc._calctype_mag){
		GUI_UNLOCK(uspex_gui.mag_nm);
		GUI_UNLOCK(uspex_gui.mag_fmls);
		GUI_UNLOCK(uspex_gui.mag_fmhs);
		GUI_UNLOCK(uspex_gui.mag_afml);
		GUI_UNLOCK(uspex_gui.mag_afmh);
		GUI_UNLOCK(uspex_gui.mag_fmlh);
		GUI_UNLOCK(uspex_gui.mag_aflh);
}
/*SURVIVAL & SELECTION*/
		GUI_UNLOCK(uspex_gui.bestFrac);
		GUI_UNLOCK(uspex_gui.keepBestHM);
		GUI_UNLOCK(uspex_gui.reoptOld);
/*VARIATION OPERATORS*/
		GUI_UNLOCK(uspex_gui.symmetries);
		GUI_UNLOCK(uspex_gui.fracGene);
		GUI_UNLOCK(uspex_gui.fracRand);
		GUI_UNLOCK(uspex_gui.fracTopRand);
		GUI_UNLOCK(uspex_gui.fracPerm);
		GUI_UNLOCK(uspex_gui.fracAtomsMut);
		GUI_UNLOCK(uspex_gui.fracLatMut);
if(uspex_gui.calc._calctype_mag){
		GUI_UNLOCK(uspex_gui.fracSpinMut);
}
		GUI_UNLOCK(uspex_gui.howManySwaps);
		GUI_UNLOCK(uspex_gui.specificSwaps);
		GUI_UNLOCK(uspex_gui.mutationDegree);
		GUI_UNLOCK(uspex_gui.mutationRate);
		GUI_UNLOCK(uspex_gui.DisplaceInLatmutation);
		GUI_UNLOCK(uspex_gui.AutoFrac);
/*CELL*/
if(!uspex_gui.auto_C_lat){
		GUI_UNLOCK(uspex_gui.Latticevalues);
		GUI_UNLOCK(uspex_gui._latticeformat);
		GUI_UNLOCK(uspex_gui._latticevalue);
}
		GUI_UNLOCK(uspex_gui.splitInto);
/*MOLECULAR*/
if(uspex_gui.calc._calctype_mol){
		GUI_UNLOCK(uspex_gui.checkConnectivity);
		GUI_UNLOCK(uspex_gui.checkMolecules);
		GUI_UNLOCK(uspex_gui.fracRotMut);
		GUI_UNLOCK(uspex_gui.MolCenters);
		GUI_UNLOCK(uspex_gui._centers);
		GUI_UNLOCK(uspex_gui._centers_button);
		GUI_UNLOCK(uspex_gui.mol_model);
		GUI_LOCK(uspex_gui.num_mol);
		GUI_LOCK(uspex_gui.mol_model_button);
		GUI_LOCK(uspex_gui.mol_gdis);
		GUI_LOCK(uspex_gui.curr_mol);
		GUI_LOCK(uspex_gui.mol_gulp);
		GUI_LOCK(uspex_gui.mol_apply_button);
		GUI_COMBOBOX_SET(uspex_gui.mol_model,0);
}
/*BoltzTraP*/
if(uspex_gui.have_ZT){
		GUI_UNLOCK(uspex_gui.BoltzTraP_T_max);
		GUI_UNLOCK(uspex_gui.BoltzTraP_T_delta);
		GUI_UNLOCK(uspex_gui.BoltzTraP_T_efcut);
		GUI_UNLOCK(uspex_gui.TE_T_interest);
		GUI_UNLOCK(uspex_gui.TE_threshold);
		GUI_UNLOCK(uspex_gui.TE_goal);
		GUI_UNLOCK(uspex_gui.cmd_BoltzTraP);
		GUI_UNLOCK(uspex_gui.cmd_BoltzTraP_button);
}
/*CRYSTAL*/
if(uspex_gui.calc._calctype_dim==3){
                GUI_UNLOCK(uspex_gui.PhaseDiagram);
}
/*SURFACES*/
if(uspex_gui.calc._calctype_dim==2){
		GUI_UNLOCK(uspex_gui.thicknessS);
		GUI_UNLOCK(uspex_gui.thicknessB);
		GUI_UNLOCK(uspex_gui.reconstruct);
		GUI_UNLOCK(uspex_gui.StoichiometryStart);
		GUI_UNLOCK(uspex_gui.E_AB);
		GUI_UNLOCK(uspex_gui.Mu_A);
		GUI_UNLOCK(uspex_gui.Mu_B);
		GUI_UNLOCK(uspex_gui.substrate_model);
		GUI_UNLOCK(uspex_gui.substrate_model_button);
}
/*CLUSTERS*/
if(uspex_gui.calc._calctype_dim==0){
                GUI_UNLOCK(uspex_gui.numberparents);
}
/*VARCOMP*/
if(uspex_gui.calc._calctype_var){
		GUI_UNLOCK(uspex_gui.numSpecies);
		GUI_UNLOCK(uspex_gui.blockSpecies);
		GUI_UNLOCK(uspex_gui.Species_apply_button);
		GUI_UNLOCK(uspex_gui.Species_delete_button);
		GUI_UNLOCK(uspex_gui.firstGeneMax);
		GUI_UNLOCK(uspex_gui.minAt);
		GUI_UNLOCK(uspex_gui.maxAt);
		GUI_UNLOCK(uspex_gui.fracTrans);
		GUI_UNLOCK(uspex_gui.howManyTrans);
		GUI_UNLOCK(uspex_gui.specificTrans);
}
/*META*/
		GUI_UNLOCK(uspex_gui.GaussianWidth);
		GUI_UNLOCK(uspex_gui.GaussianHeight);
		GUI_UNLOCK(uspex_gui.FullRelax);
		GUI_UNLOCK(uspex_gui.maxVectorLength);
		GUI_UNLOCK(uspex_gui.meta_model);
		GUI_UNLOCK(uspex_gui.meta_model_button);
		break;
	case 2://VCNEB
		GUI_UNLOCK(uspex_gui._vcnebtype_method);
		GUI_UNLOCK(uspex_gui._vcnebtype_img_num);
		GUI_UNLOCK(uspex_gui._vcnebtype_spring);
		GUI_UNLOCK(uspex_gui.numImages);
		GUI_UNLOCK(uspex_gui.numSteps);
		GUI_UNLOCK(uspex_gui.optReadImages);
		GUI_UNLOCK(uspex_gui.optimizerType);
		GUI_UNLOCK(uspex_gui.optRelaxType);
		GUI_UNLOCK(uspex_gui.dt);
		GUI_UNLOCK(uspex_gui.ConvThreshold);
		GUI_UNLOCK(uspex_gui.VarPathLength);
		GUI_UNLOCK(uspex_gui.K_min);
		GUI_UNLOCK(uspex_gui.K_max);
		GUI_UNLOCK(uspex_gui.Kconstant);
		GUI_UNLOCK(uspex_gui.optFreezing);
		GUI_UNLOCK(uspex_gui.optMethodCIDI);
		GUI_UNLOCK(uspex_gui.startCIDIStep);
		GUI_UNLOCK(uspex_gui.pickupImages);
		GUI_UNLOCK(uspex_gui.FormatType);
		GUI_UNLOCK(uspex_gui.PrintStep);
		GUI_UNLOCK(uspex_gui.img_model);
		GUI_UNLOCK(uspex_gui.img_model_button);
		break;
	case 3://PSO
		GUI_UNLOCK(uspex_gui.PSO_softMut);
		GUI_UNLOCK(uspex_gui.PSO_BestStruc);
		GUI_UNLOCK(uspex_gui.PSO_BestEver);
		break;
	case 4://TPS
		GUI_UNLOCK(uspex_gui.numIterations);
		GUI_UNLOCK(uspex_gui.speciesSymbol);
		GUI_UNLOCK(uspex_gui.mass);
		GUI_UNLOCK(uspex_gui.amplitudeShoot_AB);
		GUI_UNLOCK(uspex_gui.amplitudeShoot_BA);
		GUI_UNLOCK(uspex_gui.magnitudeShoot_success);
		GUI_UNLOCK(uspex_gui.magnitudeShoot_failure);
		GUI_UNLOCK(uspex_gui.shiftRatio);
		GUI_UNLOCK(uspex_gui.orderParaType);
		GUI_UNLOCK(uspex_gui.opCriteria_start);
		GUI_UNLOCK(uspex_gui.opCriteria_end);
		GUI_UNLOCK(uspex_gui.cmdOrderParameter);
		GUI_UNLOCK(uspex_gui.cmdOrderParameter_button);
		GUI_UNLOCK(uspex_gui.cmdEnthalpyTemperature);
		GUI_UNLOCK(uspex_gui.cmdEnthalpyTemperature_button);
		GUI_UNLOCK(uspex_gui.orderParameterFile);
		GUI_UNLOCK(uspex_gui.orderParameterFile_button);
		GUI_UNLOCK(uspex_gui.enthalpyTemperatureFile);
		GUI_UNLOCK(uspex_gui.enthalpyTemperatureFile_button);
		GUI_UNLOCK(uspex_gui.trajectoryFile);
		GUI_UNLOCK(uspex_gui.trajectoryFile_button);
		GUI_UNLOCK(uspex_gui.MDrestartFile);
		GUI_UNLOCK(uspex_gui.MDrestartFile_button);
		break;
	case 5:////MINHOP -- no specific interface
		break;
	case 6:///COPEX -- no specific interface
		break;
	default:
		break;
	}
	check_latticevalues_format();
	/*last, hide specifc page if needed*/
	hide_specific=TRUE;
	if((uspex_gui.calc.calculationMethod==US_CM_META)
	  ||(uspex_gui.calc.calculationMethod==US_CM_VCNEB)
	  ||(uspex_gui.calc.calculationMethod==US_CM_PSO))
		hide_specific=FALSE;
        hide_specific&=(uspex_gui.calc._calctype_mol==FALSE);
        hide_specific&=(uspex_gui.calc._calctype_dim!=2);
        if(!hide_specific) GUI_NOTE_PAGE_SHOW(uspex_gui.specific_page);
        else GUI_NOTE_PAGE_HIDE(uspex_gui.specific_page);
}
/*******************************/
/* selecting calculationMethod */
/*******************************/
void uspex_method_selected(GUI_OBJ *w){
	gint index;
	GUI_COMBOBOX_GET(w,index);
	/*consequences*/
	/*check method*/
	switch (index){
	case 0://USPEX
		uspex_gui.calc.calculationMethod = US_CM_USPEX;
		break;
	case 1://META
		uspex_gui.calc.calculationMethod = US_CM_META;
		break;
	case 2://VCNEB
		uspex_gui.calc.calculationMethod = US_CM_VCNEB;
		break;
	case 3://PSO
		uspex_gui.calc.calculationMethod = US_CM_PSO;
		break;
	case 4://TPS
		uspex_gui.calc.calculationMethod = US_CM_TPS;
		break;
	case 5://MINHOP -- no specific interface
		uspex_gui.calc.calculationMethod = US_CM_MINHOP;
		break;
	case 6://COPEX -- no specific interface
		uspex_gui.calc.calculationMethod = US_CM_COPEX;
		break;
	default://UNKNOWN
		uspex_gui.calc.calculationMethod = US_CM_UNKNOWN;
	}
	update_specific();
}
/*****************************/
/* selecting calculationType */
/*****************************/
void uspex_type_selected(GUI_OBJ *w){
	gint index;
	GUI_COMBOBOX_GET(w,index);
	/*check type*/
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
	case 8://001 - only in v10.1?
		uspex_gui.calc.calculationType=US_CT_001;
		uspex_gui.calc._calctype_dim=0;
		uspex_gui.calc._calctype_mol=FALSE;
		uspex_gui.calc._calctype_var=TRUE;
		uspex_gui.calc._calctype_mag=FALSE;
		break;
	case 9://110
		uspex_gui.calc.calculationType=US_CT_110;
		uspex_gui.calc._calctype_dim=1;
		uspex_gui.calc._calctype_mol=TRUE;
		uspex_gui.calc._calctype_var=FALSE;
		uspex_gui.calc._calctype_mag=FALSE;
		break;
	case 10://200
		uspex_gui.calc.calculationType=US_CT_200;
		uspex_gui.calc._calctype_dim=2;
		uspex_gui.calc._calctype_mol=FALSE;
		uspex_gui.calc._calctype_var=FALSE;
		uspex_gui.calc._calctype_mag=FALSE;
		break;
	case 11://s200
		uspex_gui.calc.calculationType=US_CT_s200;
		uspex_gui.calc._calctype_dim=2;
		uspex_gui.calc._calctype_mol=FALSE;
		uspex_gui.calc._calctype_var=FALSE;
		uspex_gui.calc._calctype_mag=TRUE;
		break;
	case 12://201
		uspex_gui.calc.calculationType=US_CT_201;
		uspex_gui.calc._calctype_dim=2;
		uspex_gui.calc._calctype_mol=FALSE;
		uspex_gui.calc._calctype_var=TRUE;
		uspex_gui.calc._calctype_mag=FALSE;
		break;
	case 13://s201
		uspex_gui.calc.calculationType=US_CT_s201;
		uspex_gui.calc._calctype_dim=2;
		uspex_gui.calc._calctype_mol=FALSE;
		uspex_gui.calc._calctype_var=TRUE;
		uspex_gui.calc._calctype_mag=TRUE;
		break;
	case 14://-200
		uspex_gui.calc.calculationType=US_CT_m200;
		uspex_gui.calc._calctype_dim=-2;
		uspex_gui.calc._calctype_mol=FALSE;
		uspex_gui.calc._calctype_var=FALSE;
		uspex_gui.calc._calctype_mag=FALSE;
		break;
	case 15://-s200
		uspex_gui.calc.calculationType=US_CT_sm200;
		uspex_gui.calc._calctype_dim=-2;
		uspex_gui.calc._calctype_mol=FALSE;
		uspex_gui.calc._calctype_var=FALSE;
		uspex_gui.calc._calctype_mag=TRUE;
		break;
	case 16://-201 v10.4
		uspex_gui.calc.calculationType=US_CT_m201;
		uspex_gui.calc._calctype_dim=-2;
		uspex_gui.calc._calctype_mol=FALSE;
		uspex_gui.calc._calctype_var=TRUE;
		uspex_gui.calc._calctype_mag=FALSE;
		break;
	case 17://-s201 v10.4
		uspex_gui.calc.calculationType=US_CT_sm201;
		uspex_gui.calc._calctype_dim=-2;
		uspex_gui.calc._calctype_mol=FALSE;
		uspex_gui.calc._calctype_var=TRUE;
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
	update_specific();
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
	case 0://000 - VER 10.1 s000 001
		i=0;
		uspex_gui.calc._calctype_mol=FALSE;
		break;
	case -2://-200 - VER 10.1 -s200
		i=200;
		uspex_gui.calc._calctype_mol=FALSE;
//		uspex_gui.calc._calctype_var=FALSE; v10.4 allows it
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
	case US_CT_001:
		GUI_COMBOBOX_SET(uspex_gui.calculationType,8);break;
	case US_CT_110:
		GUI_COMBOBOX_SET(uspex_gui.calculationType,9);break;
	case US_CT_200:
		GUI_COMBOBOX_SET(uspex_gui.calculationType,10);break;
	case US_CT_s200:
		GUI_COMBOBOX_SET(uspex_gui.calculationType,11);break;
	case US_CT_201:
		GUI_COMBOBOX_SET(uspex_gui.calculationType,12);break;
	case US_CT_s201:
		GUI_COMBOBOX_SET(uspex_gui.calculationType,13);break;
	case US_CT_m200:
		GUI_COMBOBOX_SET(uspex_gui.calculationType,14);break;
	case US_CT_sm200:
		GUI_COMBOBOX_SET(uspex_gui.calculationType,15);break;
	case US_CT_m201:
		GUI_COMBOBOX_SET(uspex_gui.calculationType,16);break;
	case US_CT_sm201:
		GUI_COMBOBOX_SET(uspex_gui.calculationType,17);break;
	case US_CT_300:
	default:
		GUI_COMBOBOX_SET(uspex_gui.calculationType,0);
	}
	update_specific();
}
/************************/
/* update _calctype_dim */
/************************/
void spin_update_dim(void){
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
	/*now enable/disable mag-related items*/
	if(uspex_gui.calc._calctype_mag){
		GUI_UNLOCK(uspex_gui.mag_nm);
		GUI_UNLOCK(uspex_gui.mag_fmls);
		GUI_UNLOCK(uspex_gui.mag_fmhs);
		GUI_UNLOCK(uspex_gui.mag_afml);
		GUI_UNLOCK(uspex_gui.mag_afmh);
		GUI_UNLOCK(uspex_gui.mag_fmlh);
		GUI_UNLOCK(uspex_gui.mag_aflh);
		GUI_UNLOCK(uspex_gui.fracSpinMut);
		
	}else{
		GUI_LOCK(uspex_gui.mag_nm);
		GUI_LOCK(uspex_gui.mag_fmls);
		GUI_LOCK(uspex_gui.mag_fmhs);
		GUI_LOCK(uspex_gui.mag_afml);
		GUI_LOCK(uspex_gui.mag_afmh);
		GUI_LOCK(uspex_gui.mag_fmlh);
		GUI_LOCK(uspex_gui.mag_aflh);
		GUI_LOCK(uspex_gui.fracSpinMut);
	}
}
/***********************************/
/* selecting optimization property */
/***********************************/
void uspex_optimization_selected(GUI_OBJ *w){
	gint index;
	GUI_COMBOBOX_GET(w,index);
	uspex_gui.have_ZT=FALSE;
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
	case 11://12
		uspex_gui.calc.optType=US_OT_HM;break;
	case 12://14
		uspex_gui.calc.optType=US_OT_ZT;
		uspex_gui.have_ZT=TRUE;
		break;
	case 13://17
		uspex_gui.calc.optType=US_OT_Fphon;break;
	case 14://1101
		uspex_gui.calc.optType=US_OT_BULK_M;break;
	case 15://1102
		uspex_gui.calc.optType=US_OT_SHEAR_M;break;
	case 16://1103
		uspex_gui.calc.optType=US_OT_YOUNG_M;break;
	case 17://1104
		uspex_gui.calc.optType=US_OT_POISSON;break;
	case 18://1105
		uspex_gui.calc.optType=US_OT_PUGH_R;break;
	case 19://1106
		uspex_gui.calc.optType=US_OT_VICKERS_H;break;
	case 20://1107
		uspex_gui.calc.optType=US_OT_FRACTURE;break;
	case 21://1108
		uspex_gui.calc.optType=US_OT_DEBYE_T;break;
	case 22://1109
		uspex_gui.calc.optType=US_OT_SOUND_V;break;
	case 23://1110
		uspex_gui.calc.optType=US_OT_SWAVE_V;break;
	case 24://1111
		uspex_gui.calc.optType=US_OT_PWAVE_V;break;
	default:
		uspex_gui.calc.optType=US_OT_UNKNOWN;
	}
	update_specific();
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
/************************************************************/
/* update numSpecies in case it is changed through atomType */
/************************************************************/
void update_numSpecies(void){
	gchar *text;
	gchar *tmp;
	gint idx;
	gint num;
	/**/
	if(uspex_gui.calc.numSpecies==NULL) return;/*FAIL*/
	if(uspex_gui.calc._calctype_var) return;/*we only update 1-line numSpecies here*/
GUI_UNLOCK(uspex_gui.numSpecies);
	GUI_COMBOBOX_WIPE(uspex_gui.numSpecies);
	text=g_strdup_printf("%i",uspex_gui.calc.numSpecies[0]);
/*there is a _BUG_ where numSpecies refers to number of molecules instead of number of elements!*/
if(uspex_gui.calc._calctype_mol) num=uspex_gui.calc._var_nspecies;
else num=uspex_gui.calc._nspecies;
	for(idx=1;idx<num;idx++) {
		tmp=g_strdup_printf("%s %i",text,uspex_gui.calc.numSpecies[idx]);
		g_free(text);
		text=tmp;
	}
	GUI_COMBOBOX_ADD(uspex_gui.numSpecies,text);
	GUI_COMBOBOX_ADD(uspex_gui.numSpecies,"ADD SPECIES BLOCK");
	GUI_COMBOBOX_SET(uspex_gui.numSpecies,0);
GUI_LOCK(uspex_gui.numSpecies);
}
/*****************************************************************/
/* Change/ADD atomType, numSpecies, and valence of selected atom */
/*****************************************************************/
void apply_atom(){
	gint index;
	gint idx;
	gchar *text;
	gboolean exists=FALSE;
	/*information*/
	gchar sym[3];
	gchar tmp[3];
	gint typ;
	gint num;
	gint val;
	gint *atomType;
	gint *numSpecies;
	/**/
	GUI_COMBOBOX_GET(uspex_gui.atomType,index);/*PUSH index*/
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
		idx=0;
		GUI_COMBOBOX_SET(uspex_gui.atomType,0);
		GUI_COMBOBOX_GET_TEXT(uspex_gui.atomType,text);
		while(g_ascii_strcasecmp(text,"ADD ATOMTYPE") != 0){
			sscanf(text,"%[^(](%*i) (V=%*i)",tmp);
			uspex_gui._tmp_atom_typ=elem_symbol_test(tmp);
			if(uspex_gui._tmp_atom_typ==typ) exists=TRUE;
			g_free(text);
			idx++;
			GUI_COMBOBOX_SET(uspex_gui.atomType,idx);
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
			g_free(text);
			/*synchronize: calc.atomType, numSpecies*/
			atomType=g_malloc((uspex_gui.calc._nspecies+1)*sizeof(gint));
			numSpecies=g_malloc((uspex_gui.calc._nspecies+1)*sizeof(gint));
			for(idx=0;idx<uspex_gui.calc._nspecies+1;idx++) {
				GUI_COMBOBOX_SET(uspex_gui.atomType,idx);
				GUI_COMBOBOX_GET_TEXT(uspex_gui.atomType,text);
				sscanf(text,"%[^(](%i) (V=%*i)",tmp,&num);
				atomType[idx]=elem_symbol_test(tmp);
				numSpecies[idx]=num;
				g_free(text);
			}
			if(uspex_gui.calc.atomType!=NULL) g_free(uspex_gui.calc.atomType);
			uspex_gui.calc.atomType=atomType;
			uspex_gui.calc._nspecies++;
			if(!uspex_gui.calc._calctype_var) {
				if(uspex_gui.calc.numSpecies!=NULL) g_free(uspex_gui.calc.numSpecies);
				uspex_gui.calc.numSpecies=numSpecies;
				update_numSpecies();
			} else g_free(numSpecies);/*FIX: 0259a1*/
			/*select last*/
			GUI_COMBOBOX_SET(uspex_gui.atomType,uspex_gui.calc._nspecies);
		}
		uspex_gui._tmp_nspecies=uspex_gui.calc._nspecies;/*BACKUP*/
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
		/*synchronize: calc.atomType, numSpecies*/
		uspex_gui.calc.atomType[index]=elem_symbol_test(sym);
		if(!uspex_gui.calc._calctype_var) {
			uspex_gui.calc.numSpecies[index]=num;
			update_numSpecies();
		}
		/*select modified atom*/
		GUI_COMBOBOX_SET(uspex_gui.atomType,index);/*PULL index*/
	}
	/*select "ADD ATOM" (for convenience)*/
	GUI_COMBOBOX_SET(uspex_gui.atomType,index);
}
/****************************************/
/* Remove selected atomType information */
/****************************************/
void remove_atom(){
	gint index;
	gchar *text;
	gint idx;
	gchar tmp[3];
	gint num;
	gint *atomType;
	gint *numSpecies;
	/**/
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
	g_free(text);
	/*synchronize: calc.atomType, numSpecies*/
	if(uspex_gui.calc._nspecies==1) return;/*this should never happen*/
	atomType=g_malloc((uspex_gui.calc._nspecies-1)*sizeof(gint));
	numSpecies=g_malloc((uspex_gui.calc._nspecies-1)*sizeof(gint));
	for(idx=0;idx<uspex_gui.calc._nspecies-1;idx++) {
		GUI_COMBOBOX_SET(uspex_gui.atomType,idx);
		GUI_COMBOBOX_GET_TEXT(uspex_gui.atomType,text);
		sscanf(text,"%[^(](%i) (V=%*i)",tmp,&num);
		atomType[idx]=elem_symbol_test(tmp);
		numSpecies[idx]=num;
		g_free(text);
	}
	if(uspex_gui.calc.atomType!=NULL) g_free(uspex_gui.calc.atomType);
	uspex_gui.calc.atomType=atomType;
	uspex_gui.calc._nspecies--;
	if(!uspex_gui.calc._calctype_var) {
		if(uspex_gui.calc.numSpecies!=NULL) g_free(uspex_gui.calc.numSpecies);
		uspex_gui.calc.numSpecies=numSpecies;
		update_numSpecies();
	} else g_free(numSpecies);/*FIX: 628b78*/
	/*select previous*/
	if(index-1<0) index=1;
	GUI_COMBOBOX_SET(uspex_gui.atomType,index-1);
	uspex_gui._tmp_nspecies=uspex_gui.calc._nspecies;/*BACKUP*/
}
/***************************/
/* select numSpecies block */
/***************************/
void uspex_numSpecies_selected(GUI_OBJ *w){
	gchar *text;
	GUI_COMBOBOX_GET_TEXT(w,text);
	if(g_ascii_strcasecmp(text,"ADD SPECIES BLOCK") == 0) return;/*no need to care about that one*/
	if(uspex_gui._tmp_blockSpecies!=NULL) g_free(uspex_gui._tmp_blockSpecies);
	uspex_gui._tmp_blockSpecies=g_strdup(text);
	GUI_ENTRY_TEXT(uspex_gui.blockSpecies,uspex_gui._tmp_blockSpecies);
}
/**************************/
/* SET initial numSpecies */
/**************************/
void set_numSpecies(void){
	gchar *text;
	gchar *tmp;
	gint idx;
	gint jdx;
	if(uspex_gui.calc._nspecies<1) return;
	if(uspex_gui.calc._var_nspecies>1){
		GUI_COMBOBOX_WIPE(uspex_gui.numSpecies);
		for(jdx=0;jdx<uspex_gui.calc._var_nspecies;jdx++){
			text=g_strdup_printf("%i",uspex_gui.calc.numSpecies[0+jdx*uspex_gui.calc._nspecies]);
			for(idx=1;idx<uspex_gui.calc._nspecies;idx++) {
				tmp=g_strdup_printf("%s %i",text,uspex_gui.calc.numSpecies[idx+jdx*uspex_gui.calc._nspecies]);
				g_free(text);
				text=tmp;
			}
			GUI_COMBOBOX_ADD_TEXT(uspex_gui.numSpecies,jdx,text);
		}
		GUI_COMBOBOX_ADD(uspex_gui.numSpecies,"ADD SPECIES BLOCK");
	}else update_numSpecies();
	/*always select the first one*/
	GUI_COMBOBOX_SET(uspex_gui.numSpecies,0);
}
/**********************************/
/* change/ADD block of numSpecies */
/**********************************/
void apply_block_species(void){
	gint index;
	gchar *text;
	if(!uspex_gui.calc._calctype_var) return;/*only varcomp allowed*/
	GUI_COMBOBOX_GET(uspex_gui.numSpecies,index);
	GUI_COMBOBOX_GET_TEXT(uspex_gui.numSpecies,text);
	if(uspex_gui._tmp_blockSpecies!=NULL) g_free(uspex_gui._tmp_blockSpecies);
	GUI_ENTRY_GET_TEXT(uspex_gui.blockSpecies,uspex_gui._tmp_blockSpecies);
	GUI_COMBOBOX_ADD_TEXT(uspex_gui.numSpecies,index,uspex_gui._tmp_blockSpecies);
	if(g_ascii_strcasecmp(text,"ADD SPECIES BLOCK") == 0){
		/*new block*/
		GUI_COMBOBOX_SET(uspex_gui.numSpecies,index+1);
	}else{
		/*modify block*/
		GUI_COMBOBOX_DEL(uspex_gui.numSpecies,index+1);
		GUI_COMBOBOX_SET(uspex_gui.numSpecies,index);
	}
	g_free(text);

}
/******************************/
/* delete block of numSpecies */
/******************************/
void delete_block_species(void){
	gint index;
	gchar *text;
	GUI_COMBOBOX_GET(uspex_gui.numSpecies,index);
	if(index==-1) {/*nothing selected -> should never happen*/
		GUI_COMBOBOX_SET(uspex_gui.numSpecies,0);
		return;
	}
	if(!uspex_gui.calc._calctype_var) return;/*only varcomp allowed*/
	GUI_COMBOBOX_GET_TEXT(uspex_gui.numSpecies,text);
	if (g_ascii_strcasecmp(text,"ADD SPECIES BLOCK") == 0){
		/*can't delete this one*/
		g_free(text);
		return;
	}
	GUI_COMBOBOX_DEL(uspex_gui.numSpecies,index);
	if(index-1<0) index=1;
	GUI_COMBOBOX_SET(uspex_gui.numSpecies,index-1);
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
		GUI_COMBOBOX_SET(uspex_gui.goodBonds,index);
	}
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
	uspex_gui.calc.have_new_optType=uspex_gui.have_new_optType;/*sync*/
	if(uspex_gui.have_new_optType){
		GUI_UNLOCK(uspex_gui.new_optType);
		GUI_LOCK(uspex_gui.optType);
	}else{
		GUI_LOCK(uspex_gui.new_optType);
		GUI_UNLOCK(uspex_gui.optType);
	}
}
/**************************************/
/* Apply changes to IonDistances line */
/**************************************/
void apply_distances(void){
	gint index;
	gchar *text;
	gchar *ptr;
	gchar *ptr2;
	gint idx;
	/**/
	GUI_COMBOBOX_GET(uspex_gui.IonDistances,index);
	GUI_ENTRY_GET_TEXT(uspex_gui._distances,text);
	/*there is _nspecies values per line (and _nspecies lines)*/
	ptr=&(text[0]);
	while((*ptr!='\0')&&(!g_ascii_isgraph(*ptr))) ptr++;
	ptr2=ptr;idx=0;
	while(*ptr!='\0'){
		uspex_gui.calc.IonDistances[index*(uspex_gui.calc._nspecies)+idx]=g_ascii_strtod(ptr,&ptr2);
		if(ptr2==ptr) break;
		ptr=ptr2+1;
		idx++;
		if(idx>=uspex_gui.calc._nspecies) break;
	}
	/*now refresch IonDistances line*/
	text=g_strdup_printf("%.4lf",uspex_gui.calc.IonDistances[index*(uspex_gui.calc._nspecies)]);
	for(idx=1;idx<uspex_gui.calc._nspecies;idx++){
		ptr=g_strdup_printf("%s %.4lf",text,uspex_gui.calc.IonDistances[index*(uspex_gui.calc._nspecies)+idx]);
		g_free(text);
		text=ptr;
	}
	GUI_COMBOBOX_ADD_TEXT(uspex_gui.IonDistances,index,text);
	GUI_COMBOBOX_DEL(uspex_gui.IonDistances,index+1);
	GUI_COMBOBOX_SET(uspex_gui.IonDistances,index);
	g_free(text);
}
/************************************/
/* Apply changes to MolCenters line */
/************************************/
void apply_centers(void){
	gint index;
	gchar *text;
	gchar *ptr;
	gchar *ptr2;
	gint idx;
	/**/
	if(uspex_gui.calc._nmolecules==0) return;
	GUI_COMBOBOX_GET(uspex_gui.MolCenters,index);
	GUI_ENTRY_GET_TEXT(uspex_gui._centers,text);
	/*there is _nmolecules values per line (and _nmolecules lines)*/
	ptr=&(text[0]);
	while((*ptr!='\0')&&(!g_ascii_isgraph(*ptr))) ptr++;
	ptr2=ptr;idx=0;
	while(*ptr!='\0'){
		uspex_gui.calc.MolCenters[index*(uspex_gui.calc._nmolecules)+idx]=g_ascii_strtod(ptr,&ptr2);
		if(ptr2==ptr) break;
		ptr=ptr2+1;
		idx++;
		if(idx>=uspex_gui.calc._nmolecules) break;
	}
	/*now refresh MolCenters line*/
	text=g_strdup_printf("%.4lf",uspex_gui.calc.MolCenters[index*(uspex_gui.calc._nmolecules)]);
	for(idx=1;idx<uspex_gui.calc._nmolecules;idx++){
		ptr=g_strdup_printf("%s %.4lf",text,uspex_gui.calc.MolCenters[index*(uspex_gui.calc._nmolecules)+idx]);
		g_free(text);
		text=ptr;
	}
	GUI_COMBOBOX_ADD_TEXT(uspex_gui.MolCenters,index,text);
	GUI_COMBOBOX_DEL(uspex_gui.MolCenters,index+1);
	GUI_COMBOBOX_SET(uspex_gui.MolCenters,index);
	g_free(text);
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
		/*NEW: create calc.IonDistances if not here...*/
		if(uspex_gui.calc.IonDistances==NULL){
			/*malloc0 to avoid a possible use of uninitialized value*/
			uspex_gui.calc.IonDistances = g_malloc0(uspex_gui.calc._nspecies*uspex_gui.calc._nspecies*sizeof(gdouble));
			for(i=0;i<uspex_gui.calc._nspecies*uspex_gui.calc._nspecies;i++) 
				uspex_gui.calc.IonDistances[i]=0.;/*no default*/
		}
		for(i=0;i<uspex_gui.calc._nspecies;i++){
			tamp=g_strdup_printf("%.4lf",uspex_gui.calc.IonDistances[i*uspex_gui.calc._nspecies]);
			text=tamp;
			for(j=1;j<uspex_gui.calc._nspecies;j++) {
				text=g_strdup_printf("%s %.4lf",tamp,uspex_gui.calc.IonDistances[j+i*uspex_gui.calc._nspecies]);
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
	if((uspex_gui.calc._nmolecules>0)&&(uspex_gui.calc.MolCenters!=NULL)){
		/*refresh MolCenters*/
		GUI_UNLOCK(uspex_gui.MolCenters);
		GUI_UNLOCK(uspex_gui._centers);
		/*wipe MolCenters, rewrite*/
		GUI_COMBOBOX_WIPE(uspex_gui.MolCenters);
		for(i=0;i<uspex_gui.calc._nmolecules;i++){
			tamp=g_strdup_printf("%.4lf",uspex_gui.calc.MolCenters[i*uspex_gui.calc._nmolecules]);
			text=tamp;
			for(j=1;j<uspex_gui.calc._nmolecules;j++) {
				text=g_strdup_printf("%s %.4lf",tamp,uspex_gui.calc.MolCenters[j+i*uspex_gui.calc._nmolecules]);
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
	gint nvals;
	gchar *ptr;
	gchar *ptr2;
	GUI_COMBOBOX_GET(uspex_gui.Latticevalues,index);
	if(uspex_gui._tmp_latticevalue!=NULL) g_free(uspex_gui._tmp_latticevalue);
	GUI_ENTRY_GET_TEXT(uspex_gui._latticevalue,uspex_gui._tmp_latticevalue);
	GUI_COMBOBOX_ADD_TEXT(uspex_gui.Latticevalues,index,uspex_gui._tmp_latticevalue);
	GUI_COMBOBOX_DEL(uspex_gui.Latticevalues,index+1);
	GUI_COMBOBOX_SET(uspex_gui.Latticevalues,index);
/*keep track of vals*/
	if(uspex_gui.calc._nlattice_line>1){
		/*matrix has to be rectangular*/
		ptr=uspex_gui._tmp_latticevalue;
		nvals=0;
		while((*ptr!='\0')&&(!g_ascii_isgraph(*ptr))) ptr++;
		while((*ptr!='\0')&&(nvals<uspex_gui.calc._nlattice_vals)){
			/**/
			uspex_gui.calc.Latticevalues[nvals+index*uspex_gui.calc._nlattice_line]=g_ascii_strtod(ptr,&ptr2);
			ptr=ptr2+1;
			while((*ptr!='\0')&&(!g_ascii_isgraph(*ptr))) ptr++;
			nvals++;
		}
	}else{
		/*we need to count the number of elements*/
		ptr=uspex_gui._tmp_latticevalue;
		nvals=0;
		while((*ptr!='\0')&&(!g_ascii_isgraph(*ptr))) ptr++;
		while((*ptr!='\0')&&(*ptr!='\n')){
			nvals++;
			while(g_ascii_isgraph(*ptr)) ptr++;/*to next space EOL*/
			while((*ptr!='\0')&&(!g_ascii_isgraph(*ptr))) ptr++;/*to next item*/
		}
		/*wipe and replace*/
		uspex_gui.calc._nlattice_vals=nvals;
		if(uspex_gui.calc.Latticevalues!=NULL) g_free(uspex_gui.calc.Latticevalues);/*FIX: 23188d*/
		uspex_gui.calc.Latticevalues=g_malloc0(nvals*sizeof(gdouble));
		ptr=uspex_gui._tmp_latticevalue;
		nvals=0;
		while((*ptr!='\0')&&(!g_ascii_isgraph(*ptr))) ptr++;
		while((*ptr!='\0')&&(nvals<uspex_gui.calc._nlattice_vals)){
			uspex_gui.calc.Latticevalues[nvals]=g_ascii_strtod(ptr,&ptr2);
			ptr=ptr2+1;
			while((*ptr!='\0')&&(!g_ascii_isgraph(*ptr))) ptr++;
			nvals++;
		}
	}
}
/***************************/
/* Selecting latticeformat */
/***************************/
void uspex_latticeformat_selected(GUI_OBJ *w){
	gint index;
	gint idx;
	gint jdx;
	gchar *text;
	gchar *tmp;
	GUI_COMBOBOX_GET(w,index);
	if(index==-1) return;/*nothing selected*/
	GUI_COMBOBOX_WIPE(uspex_gui.Latticevalues);/*wipe in any case*/
	/*new: we select based on current information*/
	if(uspex_gui.calc._calctype_var){
		/*only Volume*/
		if(uspex_gui.calc._nlattice_line!=1) uspex_gui.calc._nlattice_line=1;
		if(uspex_gui.calc._nlattice_vals>0){
			text=g_strdup_printf("%.4lf",uspex_gui.calc.Latticevalues[0]);
			for(idx=0;idx<uspex_gui.calc._nlattice_vals;idx++) {
				tmp=g_strdup_printf("%s %.4lf",text,uspex_gui.calc.Latticevalues[idx]);
				g_free(text);
				text=tmp;
			}
		}else{
			uspex_gui.calc._nlattice_vals=1;
			if(uspex_gui.calc.Latticevalues!=NULL) g_free(uspex_gui.calc.Latticevalues);
			uspex_gui.calc.Latticevalues=g_malloc(1*sizeof(gdouble));
			uspex_gui.calc.Latticevalues[0]=0.;
			text=g_strdup_printf("%.4lf",uspex_gui.calc.Latticevalues[0]);
		}
		GUI_COMBOBOX_ADD(uspex_gui.Latticevalues,text);
		g_free(text);
	}else{
		/*for now, switching from one type to another wipe everything TODO crystal <-> lattice*/
		switch (index){
		case 1://lattice parameters
			if((uspex_gui.calc._calctype_dim==3)&&(uspex_gui.calc._nlattice_line==3)&&(uspex_gui.calc._nlattice_vals==3)){
				/*do not wipe*/
			}else if((uspex_gui.calc._calctype_dim==2)&&(uspex_gui.calc._nlattice_line==2)&&(uspex_gui.calc._nlattice_vals==2)){
				/*do not wipe*/
			}else{
				/*wipe*/
				if(uspex_gui.calc._calctype_dim==3){
					uspex_gui.calc._nlattice_line=3;
					uspex_gui.calc._nlattice_vals=3;
				}else{
					uspex_gui.calc._nlattice_line=2;
					uspex_gui.calc._nlattice_vals=2;
				}
				if(uspex_gui.calc.Latticevalues!=NULL) g_free(uspex_gui.calc.Latticevalues);
uspex_gui.calc.Latticevalues=g_malloc0((uspex_gui.calc._nlattice_line*uspex_gui.calc._nlattice_vals)*sizeof(gdouble));/*FIX: 400a23*/
				for(idx=0;idx<uspex_gui.calc._nlattice_line*uspex_gui.calc._nlattice_vals;idx++)
					uspex_gui.calc.Latticevalues[idx]=0.;
			}
			/*fill up lines*/
			for(jdx=0;jdx<uspex_gui.calc._nlattice_line;jdx++){
				text=g_strdup_printf("%.4lf",uspex_gui.calc.Latticevalues[jdx*uspex_gui.calc._nlattice_vals]);
				for(idx=1;idx<uspex_gui.calc._nlattice_vals;idx++){
					tmp=g_strdup_printf("%s %.4lf",text,uspex_gui.calc.Latticevalues[idx+jdx*uspex_gui.calc._nlattice_vals]);
					g_free(text);
					text=tmp;
				}
				GUI_COMBOBOX_ADD(uspex_gui.Latticevalues,text);
				g_free(text);
			}
			break;
		case 2://crystal definition
			if((uspex_gui.calc._calctype_dim==3)&&(uspex_gui.calc._nlattice_line==1)&&(uspex_gui.calc._nlattice_vals==6)){
				/*do not wipe*/
			}else if((uspex_gui.calc._calctype_dim==2)&&(uspex_gui.calc._nlattice_line==1)&&(uspex_gui.calc._nlattice_vals==3)){
				/*do not wipe*/
			}else{
				/*wipe, also enforce "standard" crystal definition
				 *ie. a, b, c, alpha, beta, gamma	-- OVHPA*/
				uspex_gui.calc._nlattice_line=1;
				uspex_gui.calc._nlattice_vals=6;
				if(uspex_gui.calc.Latticevalues!=NULL) g_free(uspex_gui.calc.Latticevalues);
uspex_gui.calc.Latticevalues=g_malloc0((uspex_gui.calc._nlattice_line*uspex_gui.calc._nlattice_vals)*sizeof(gdouble));/*FIX: 018107*/
				for(idx=0;idx<uspex_gui.calc._nlattice_line*uspex_gui.calc._nlattice_vals;idx++)
					uspex_gui.calc.Latticevalues[idx]=0.;
			}
			/*fill up the line*/
			text=g_strdup_printf("%.4lf",uspex_gui.calc.Latticevalues[0]);
			for(idx=0;idx<uspex_gui.calc._nlattice_vals;idx++){
				tmp=g_strdup_printf("%s %.4lf",text,uspex_gui.calc.Latticevalues[idx]);
				g_free(text);
				text=tmp;
			}
			GUI_COMBOBOX_ADD(uspex_gui.Latticevalues,text);
			g_free(text);
			break;
		case 0://volume definition
			/* pass through */
		default:
			if(uspex_gui.calc._nlattice_line==1){
				/*do not wipe*/
			}else{
				/*wipe*/
				uspex_gui.calc._nlattice_line=1;
				uspex_gui.calc._nlattice_vals=1;
				if(uspex_gui.calc.Latticevalues!=NULL) g_free(uspex_gui.calc.Latticevalues);
uspex_gui.calc.Latticevalues=g_malloc((uspex_gui.calc._nlattice_line*uspex_gui.calc._nlattice_vals)*sizeof(gdouble));
				uspex_gui.calc.Latticevalues[0]=0.;
			}
			text=g_strdup_printf("%.4lf",uspex_gui.calc.Latticevalues[0]);
			for(idx=1;idx<uspex_gui.calc._nlattice_vals;idx++){
				tmp=g_strdup_printf("%s %.4lf",text,uspex_gui.calc.Latticevalues[idx]);
				g_free(text);
				text=tmp;
			}
			GUI_COMBOBOX_ADD(uspex_gui.Latticevalues,text);
			g_free(text);
		}
	}
	GUI_COMBOBOX_SET(uspex_gui.Latticevalues,0);/*always go to the first one*/
}
/**********************************************/
/* Guess which is the format of LatticeValues */
/**********************************************/
void set_lattice_format(){
	if(uspex_gui.auto_C_lat) return;/*no need*/
	if(uspex_gui.calc._calctype_var){
		GUI_COMBOBOX_SET(uspex_gui._latticeformat,0);
		return;
	}
	if(uspex_gui.calc._nlattice_line>1) {
		GUI_COMBOBOX_SET(uspex_gui._latticeformat,1);/*lattice*/
		return;
	}
	switch(uspex_gui.calc._nlattice_vals){
	case 3:
		if(uspex_gui.calc._calctype_dim==-2) 
			GUI_COMBOBOX_SET(uspex_gui._latticeformat,2);/*crystal*/
		else GUI_COMBOBOX_SET(uspex_gui._latticeformat,0);/*volume*/
		return;
	case 6:
		if((uspex_gui.calc._calctype_dim==-2)||(uspex_gui.calc._calctype_dim==3)) 
			GUI_COMBOBOX_SET(uspex_gui._latticeformat,2);/*crystal*/
		else GUI_COMBOBOX_SET(uspex_gui._latticeformat,0);/*volume*/
		return;
	default:
		GUI_COMBOBOX_SET(uspex_gui._latticeformat,0);/*Volume by default*/
	}
}
/********************/
/* toggle auto_C_lat */
/********************/
void toggle_auto_C_lat(void){
	if(uspex_gui.auto_C_lat){
		/*automatic*/
		GUI_LOCK(uspex_gui.Latticevalues);
		GUI_LOCK(uspex_gui._latticevalue);
		GUI_LOCK(uspex_gui._latticeformat);
	}else{
		/*manual*/
		/*NEW: create a base Latticevalues if none present*/
		if(uspex_gui.calc.Latticevalues==NULL) {
			if(uspex_gui.calc._nlattice_line==0) {
				uspex_gui.calc._nlattice_line=1;
				uspex_gui.calc._nlattice_vals=1;
			}
			uspex_gui.calc.Latticevalues=
				g_malloc((uspex_gui.calc._nlattice_line*uspex_gui.calc._nlattice_vals)*sizeof(gdouble));
			GUI_COMBOBOX_SET(uspex_gui._latticeformat,0);
			set_lattice_format();
			uspex_latticeformat_selected(uspex_gui._latticeformat);
		}
		GUI_UNLOCK(uspex_gui.Latticevalues);
		GUI_UNLOCK(uspex_gui._latticevalue);
		GUI_UNLOCK(uspex_gui._latticeformat);
	}
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
/*****************************************************************/
/* convert from _tmp_commandExecutable to calc.commandExecutable */
/*****************************************************************/
void register_commandExecutable(void){
	gchar *text;
	gint idx;
	gint jdx;
	gboolean repopulate=FALSE;
	/*1- wipe calc.commandExecutable*/
	if(uspex_gui.calc.commandExecutable!=NULL) g_free(uspex_gui.calc.commandExecutable);
	uspex_gui.calc.commandExecutable=NULL;
	/*2- register mandatory element*/
	if(uspex_gui._tmp_commandExecutable[0]==NULL) return;/*nothing to update at all*/
	uspex_gui.calc.commandExecutable=g_strdup(uspex_gui._tmp_commandExecutable[0]);
	for(idx=1;idx<uspex_gui.calc._num_opt_steps;idx++){
		if(uspex_gui._tmp_commandExecutable[idx]==NULL) continue;/*skip it*/
		text=g_strdup_printf("%s\n%s",uspex_gui.calc.commandExecutable,uspex_gui._tmp_commandExecutable[idx]);
		g_free(uspex_gui.calc.commandExecutable);
		uspex_gui.calc.commandExecutable=text;
		if(!uspex_gui.calc._isCmdList) repopulate=TRUE;
		uspex_gui.calc._isCmdList=TRUE;
	}
	if(repopulate){
		/*set each step explicitely*/
		//1bis- re-wipe
		g_free(uspex_gui.calc.commandExecutable);
		uspex_gui.calc.commandExecutable=NULL;
		uspex_gui.calc.commandExecutable=g_strdup(uspex_gui._tmp_commandExecutable[0]);
		for(idx=1;idx<uspex_gui.calc._num_opt_steps;idx++){
			if(uspex_gui._tmp_commandExecutable[idx]==NULL){
				/*not acceptable, have to be set*/
				if(uspex_gui.calc.abinitioCode[idx]==uspex_gui.calc.abinitioCode[idx-1]){
					/*same as previous one*/
					uspex_gui._tmp_commandExecutable[idx]=g_strdup(uspex_gui._tmp_commandExecutable[idx-1]);
				}else{
					/*find a similar one*/
					for(jdx=1;jdx<uspex_gui.calc._num_opt_steps;jdx++){
						if(jdx==idx) continue;
						if(uspex_gui.calc.abinitioCode[idx]==uspex_gui.calc.abinitioCode[jdx]){
							uspex_gui._tmp_commandExecutable[idx]=g_strdup(uspex_gui._tmp_commandExecutable[jdx]);
						}
					}
					/*if all fail*/
					if(uspex_gui._tmp_commandExecutable[idx]==NULL) uspex_gui._tmp_commandExecutable[idx]=g_strdup("N/A");
				}
			}
			text=g_strdup_printf("%s\n%s",uspex_gui.calc.commandExecutable,uspex_gui._tmp_commandExecutable[idx]);
			g_free(uspex_gui.calc.commandExecutable);
			uspex_gui.calc.commandExecutable=text;
		}
	}
}
/*****************************************/
/* sync the flavor of (pseudo-)potential */
/*****************************************/
void sync_ai_lib_flavor(void){
	gchar *path, *text;
	GSList *folders,*folder;
	gchar  *item;
	gchar sym[3];
	gint idx,sel;
	gint step=(gint)uspex_gui._tmp_curr_step-1;
	if(uspex_gui._tmp_ai_lib_folder[step]==NULL) {
		GUI_COMBOBOX_WIPE(uspex_gui.ai_lib_flavor);
		/*library folder is blank but there could be a flavor entered manually*/
		if(uspex_gui._tmp_ai_lib_sel[step]!=NULL){
			GUI_COMBOBOX_ADD(uspex_gui.ai_lib_flavor,uspex_gui._tmp_ai_lib_sel[step]);
		}
		GUI_COMBOBOX_ADD(uspex_gui.ai_lib_flavor,"N/A");/*terminator*/
		GUI_COMBOBOX_SET(uspex_gui.ai_lib_flavor,0);
		return;/*no folder selected*/
	}
	path=g_strdup(uspex_gui._tmp_ai_lib_folder[step]);
	folders=file_dir_list(path,FALSE);
	g_free(path);
	GUI_COMBOBOX_WIPE(uspex_gui.ai_lib_flavor);
	sym[2]='\0';
	for (folder=folders ; folder ; folder=g_slist_next(folder)){
		item=g_strdup_printf("%s",(char *)folder->data);
		switch (uspex_gui.calc.abinitioCode[step]){
		case 1://VASP
		case 2://SIESTA
		case 8://Q.E.
		case 10://ATK
		case 11://CASTEP
		case 18://Abinit
			/*per XX species: 
			 * 1- XX{*}/ folder for VASP PP,
			 * 2- XX.* files for the others. */
			for(idx=0;idx<uspex_gui.calc._nspecies;idx++){
#if DEBUG_POTCAR_FOLDER
	fprintf(stdout,"#DBG: looking for %s\n",sym);
#endif  
				sym[0]=elements[uspex_gui.calc.atomType[idx]].symbol[0];
				sym[1]=elements[uspex_gui.calc.atomType[idx]].symbol[1];
				if(sym[0]==item[0]){
					if(sym[1]==item[1]) GUI_COMBOBOX_ADD(uspex_gui.ai_lib_flavor,item);
					else{
						/*special case of single atoms*/
						if((!g_ascii_islower (sym[1]))
						 &&(!g_ascii_islower (item[1]))){
							GUI_COMBOBOX_ADD(uspex_gui.ai_lib_flavor,item);
						}
					}
				}
			}
			break;
		case 3://GULP
			/* check *.lib */
			if (find_in_string(".lib",item) != NULL) 
				GUI_COMBOBOX_ADD(uspex_gui.ai_lib_flavor,item);
			break;
		case 12://Tinker
			/* check *.prm */
			if (find_in_string(".prm",item) != NULL)
				GUI_COMBOBOX_ADD(uspex_gui.ai_lib_flavor,item);
			break;
		case 7://CP2K
			/* check * */
			/*add everything that is a file*/
			if(item[0] != '.')
				GUI_COMBOBOX_ADD(uspex_gui.ai_lib_flavor,item);
			break;
		case 15://DFTB
			/* Special: all A, B species skf have:
			 * A-A.skf, A-B.skf, B-A.skf, B-B.skf.
			 * Files will be copied automatically.*/
			if(folder==folders) GUI_COMBOBOX_ADD(uspex_gui.ai_lib_flavor,"AUTO");
			break;
		case 4://LAMMPS
		case 13://MOPAC
		case 14://BoltzTraP
		case 16://Gaussian
		case 17://N/A
		case 19://CRYSTAL
			/* No library is needed... (I think)
			 * So this setting imports all files
			 * from library directory.*/
		case 5://ORCA
		case 6://DMACRYS
		case 9://FHI-aims
			/* Some library could be needed, but
			 * I can't figure out format/case...
			 * so we also imports all files from
			 * library directory. (just in case)*/
			if(folder==folders) GUI_COMBOBOX_ADD(uspex_gui.ai_lib_flavor,"ALL");
			break;
		default:
			/*unsupported code?*/
			break;
		}
		g_free(item);
	}
	GUI_COMBOBOX_ADD(uspex_gui.ai_lib_flavor,"N/A");/*terminator*/
	g_slist_free(folders);
	/*determine if there was an element selected*/
	GUI_COMBOBOX_SET(uspex_gui.ai_lib_flavor,0);
	if(uspex_gui._tmp_ai_lib_sel[step]!=NULL){
		/*there is a selection, is it in the list?*/
		idx=0;sel=0;text=NULL;
		do{
			GUI_COMBOBOX_SET(uspex_gui.ai_lib_flavor,idx);
			GUI_COMBOBOX_GET_TEXT(uspex_gui.ai_lib_flavor,text);
			if(text==NULL) break;/*if nothing?*/
			if(find_in_string(text,"N/A")!=NULL) break;/*end here*/
			if(find_in_string(uspex_gui._tmp_ai_lib_sel[step],text) != NULL) sel=idx;
			g_free(text);
			idx++;
		}while(TRUE);/*endless loop*/
		GUI_COMBOBOX_SET(uspex_gui.ai_lib_flavor,sel);
	}

}
/*************************************/
/* sync the library flavor selection */
/*************************************/
void sync_ai_lib_sel(void){
	gchar *text;
	gint step=(gint)uspex_gui._tmp_curr_step-1;
	/**/
	if(uspex_gui._tmp_ai_lib_sel[step]==NULL)
		text=g_strdup_printf("N/A");
	else
		text=g_strdup_printf("%s",uspex_gui._tmp_ai_lib_sel[step]);
	GUI_UNLOCK(uspex_gui.ai_lib_sel);
	GUI_ENTRY_TEXT(uspex_gui.ai_lib_sel,text);
	GUI_LOCK(uspex_gui.ai_lib_sel);
	g_free(text);
}
/********************************/
/* set ai_library for last step */
/********************************/
void register_ai_library(void){
	gint idx;
	gint code;
	gint clone;
	gint tot=(gint)uspex_gui._tmp_num_opt_steps;
	gint step=(gint)uspex_gui._tmp_curr_step-1;
	/*determine if _tmp_ai_lib_folder needs an update*/
	if(uspex_gui._tmp_ai_lib_folder[step]!=NULL) return;
	/*search if we already have a setup for this ai_code*/
	code=uspex_gui.calc.abinitioCode[step];
	clone=step;
	for(idx=0;idx<tot;idx++)
		if((code==uspex_gui.calc.abinitioCode[idx])&&(idx!=step)) clone=idx;
	if(clone!=step){
		/*just clone previous similar step*/
		if(uspex_gui._tmp_ai_lib_folder[clone]!=NULL){
			uspex_gui._tmp_ai_lib_folder[step]=g_strdup(uspex_gui._tmp_ai_lib_folder[clone]);
			/*by default, selected flavor also propagate*/
			if(uspex_gui._tmp_ai_lib_sel[clone]!=NULL)
				uspex_gui._tmp_ai_lib_sel[step]=g_strdup(uspex_gui._tmp_ai_lib_sel[clone]);
		} else uspex_gui._tmp_ai_lib_folder[step]=NULL;/*already set*/
	}else{
		/*ensure NULL?*/
	}
	/*wipe & update uspex_gui.ai_lib_flavor*/
	sync_ai_lib_flavor();
	sync_ai_lib_sel();
}
/******************************/
/* toggle USE Specific folder */
/******************************/
void toggle_use_specific(void){
uspex_gui.have_specific=TRUE;
GUI_TOGGLE_ON(uspex_gui.use_specific);
GUI_TOGGLE_OFF(uspex_gui.set_specific);
GUI_UNLOCK(uspex_gui.spe_folder);
GUI_UNLOCK(uspex_gui.spe_folder_button);
GUI_LOCK(uspex_gui.ai_input);
GUI_LOCK(uspex_gui.ai_input_button);
GUI_LOCK(uspex_gui.ai_opt);
GUI_LOCK(uspex_gui.ai_opt_button);
GUI_LOCK(uspex_gui.ai_lib);
GUI_LOCK(uspex_gui.ai_lib_button);
GUI_LOCK(uspex_gui.ai_lib_flavor);
GUI_LOCK(uspex_gui.ai_lib_sel);
}
/**********************************/
/* toggle SET specific parameters */
/**********************************/
void toggle_set_specific(void){
uspex_gui.have_specific=FALSE;
GUI_TOGGLE_OFF(uspex_gui.use_specific);
GUI_TOGGLE_ON(uspex_gui.set_specific);
GUI_LOCK(uspex_gui.spe_folder);
GUI_LOCK(uspex_gui.spe_folder_button);
GUI_UNLOCK(uspex_gui.ai_input);
GUI_UNLOCK(uspex_gui.ai_input_button);
GUI_UNLOCK(uspex_gui.ai_opt);
GUI_UNLOCK(uspex_gui.ai_opt_button);
GUI_UNLOCK(uspex_gui.ai_lib);
GUI_UNLOCK(uspex_gui.ai_lib_button);
GUI_UNLOCK(uspex_gui.ai_lib_flavor);
GUI_UNLOCK(uspex_gui.ai_lib_sel);
}
/************************/
/* load Specific folder */
/************************/
void load_spe_folder_dialog(void){
        GUI_OBJ *file_chooser;
        gint have_answer;
        gchar *filename;
        gchar *text;
        /**/
        GUI_PREPARE_OPEN_FOLDER(uspex_gui.window,file_chooser,"Select the Specific folder");
        GUI_OPEN_DIALOG_RUN(file_chooser,have_answer,filename);
        if(have_answer){
                /*there should be no case where have_answer==TRUE AND filename==NULL, but just in case.*/
                if(filename) {
                        text=g_strdup_printf("%s",filename);
                        GUI_ENTRY_TEXT(uspex_gui.spe_folder,text);
                        g_free(text);
                        if(uspex_gui._tmp_spe_folder!=NULL) g_free(uspex_gui._tmp_spe_folder);
                        uspex_gui._tmp_spe_folder=g_strdup_printf("%s",filename);
                        g_free (filename);
                }
        }     
        GUI_KILL_OPEN_DIALOG(file_chooser);
}
/*************************************************/
/* change the total number of optimisation steps */
/*************************************************/
void spin_update_num_opt_steps(void){
	gint     *tmp_i;
	gdouble *tmp_d0;
	gdouble *tmp_d1;
	gboolean *tmp_b;
	gchar   **tmp_c;
	gchar **tmp_inp;
	gchar **tmp_opt;
	gchar **tmp_lib;
	gchar **tmp_sel;
	gint i=(gint)uspex_gui._tmp_num_opt_steps;
	gint idx;
	/**/
	tmp_i=g_malloc(i*sizeof(gint));
	tmp_d0=g_malloc(i*sizeof(gdouble));
	tmp_d1=g_malloc(i*sizeof(gdouble));
	tmp_b=g_malloc(i*sizeof(gboolean));
	tmp_c=g_malloc(i*sizeof(gchar *));
	tmp_inp=g_malloc(i*sizeof(gchar *));
	tmp_opt=g_malloc(i*sizeof(gchar *));
	tmp_lib=g_malloc(i*sizeof(gchar *));
	tmp_sel=g_malloc(i*sizeof(gchar *));
	if(i>uspex_gui.calc._num_opt_steps){
		/*increase*/
		for(idx=0;idx<uspex_gui.calc._num_opt_steps;idx++){
			tmp_i[idx]=uspex_gui.calc.abinitioCode[idx];
			tmp_b[idx]=uspex_gui.calc._isfixed[idx];
			tmp_d0[idx]=uspex_gui.calc.KresolStart[idx];
			tmp_d1[idx]=uspex_gui.calc.vacuumSize[idx];
			if(uspex_gui._tmp_commandExecutable[idx]!=NULL){
				tmp_c[idx]=g_strdup(uspex_gui._tmp_commandExecutable[idx]);
				g_free(uspex_gui._tmp_commandExecutable[idx]);
				uspex_gui._tmp_commandExecutable[idx]=NULL;
			}else tmp_c[idx]=NULL;
			if(uspex_gui._tmp_ai_input[idx]!=NULL){
				tmp_inp[idx]=g_strdup(uspex_gui._tmp_ai_input[idx]);
				g_free(uspex_gui._tmp_ai_input[idx]);
				uspex_gui._tmp_ai_input[idx]=NULL;
			}else tmp_inp[idx]=NULL;
			if(uspex_gui._tmp_ai_opt[idx]!=NULL){
				tmp_opt[idx]=g_strdup(uspex_gui._tmp_ai_opt[idx]);
				g_free(uspex_gui._tmp_ai_opt[idx]);
				uspex_gui._tmp_ai_opt[idx]=NULL;
			}else tmp_opt[idx]=NULL;
			if(uspex_gui._tmp_ai_lib_folder!=NULL){
				tmp_lib[idx]=g_strdup(uspex_gui._tmp_ai_lib_folder[idx]);
				g_free(uspex_gui._tmp_ai_lib_folder[idx]);
				uspex_gui._tmp_ai_lib_folder[idx]=NULL;
			}else tmp_lib[idx]=NULL;
			if(uspex_gui._tmp_ai_lib_sel!=NULL){
				tmp_sel[idx]=g_strdup(uspex_gui._tmp_ai_lib_sel[idx]);
				g_free(uspex_gui._tmp_ai_lib_sel[idx]);
				uspex_gui._tmp_ai_lib_sel[idx]=NULL;
			}else tmp_sel[idx]=NULL;
		}
		/*last step*/
		tmp_i[i-1]=1;
		tmp_b[i-1]=FALSE;
		tmp_d0[i-1]=0.2-(gdouble)(i-1)*(0.2-0.08)/(uspex_gui.calc._num_opt_steps);
		tmp_d1[i-1]=10.;
		tmp_c[i-1]=NULL;/*find a better default?*/
		tmp_inp[i-1]=g_strdup("N/A");
		tmp_opt[i-1]=g_strdup("N/A");
		tmp_lib[i-1]=NULL;/*updated later*/
		tmp_sel[i-1]=NULL;/*updated later*/
	}else{
		/*decrease*/
		for(idx=0;idx<i;idx++){
			tmp_i[idx]=uspex_gui.calc.abinitioCode[idx];
			tmp_b[idx]=uspex_gui.calc._isfixed[idx];
			tmp_d0[idx]=uspex_gui.calc.KresolStart[idx];
			tmp_d1[idx]=uspex_gui.calc.vacuumSize[idx];
			if(uspex_gui._tmp_commandExecutable[idx]!=NULL){
				tmp_c[idx]=g_strdup(uspex_gui._tmp_commandExecutable[idx]);
				g_free(uspex_gui._tmp_commandExecutable[idx]);
				uspex_gui._tmp_commandExecutable[idx]=NULL;
			}else tmp_c[idx]=NULL;
			if(uspex_gui._tmp_ai_input[idx]!=NULL){
				tmp_inp[idx]=g_strdup(uspex_gui._tmp_ai_input[idx]);
				g_free(uspex_gui._tmp_ai_input[idx]);
				uspex_gui._tmp_ai_input[idx]=NULL;
			}else tmp_inp[idx]=NULL;
			if(uspex_gui._tmp_ai_opt[idx]!=NULL){
				tmp_opt[idx]=g_strdup(uspex_gui._tmp_ai_opt[idx]);
				g_free(uspex_gui._tmp_ai_opt[idx]);
				uspex_gui._tmp_ai_opt[idx]=NULL;
			}else tmp_opt[idx]=NULL;
			if(uspex_gui._tmp_ai_lib_folder!=NULL){
				tmp_lib[idx]=g_strdup(uspex_gui._tmp_ai_lib_folder[idx]);
				g_free(uspex_gui._tmp_ai_lib_folder[idx]);
				uspex_gui._tmp_ai_lib_folder[idx]=NULL;
			}else tmp_lib[idx]=NULL;
			if(uspex_gui._tmp_ai_lib_sel!=NULL){
				tmp_sel[idx]=g_strdup(uspex_gui._tmp_ai_lib_sel[idx]);
				g_free(uspex_gui._tmp_ai_lib_sel[idx]);
				uspex_gui._tmp_ai_lib_sel[idx]=NULL;
			}else tmp_sel[idx]=NULL;
		}
	}
	g_free(uspex_gui.calc.abinitioCode);uspex_gui.calc.abinitioCode=tmp_i;
	g_free(uspex_gui.calc._isfixed);uspex_gui.calc._isfixed=tmp_b;
	g_free(uspex_gui.calc.KresolStart);uspex_gui.calc.KresolStart=tmp_d0;
	g_free(uspex_gui.calc.vacuumSize);uspex_gui.calc.vacuumSize=tmp_d1;
	g_free(uspex_gui._tmp_commandExecutable);uspex_gui._tmp_commandExecutable=tmp_c;
	g_free(uspex_gui._tmp_ai_input);uspex_gui._tmp_ai_input=tmp_inp;/*FIX d09393*/
	g_free(uspex_gui._tmp_ai_opt);uspex_gui._tmp_ai_opt=tmp_opt;/*FIX 4f5cad*/
	g_free(uspex_gui._tmp_ai_lib_folder);uspex_gui._tmp_ai_lib_folder=tmp_lib;
	g_free(uspex_gui._tmp_ai_lib_sel);uspex_gui._tmp_ai_lib_sel=tmp_sel;
	/*update _num_opt_steps*/
	uspex_gui.calc._num_opt_steps=i;
	register_commandExecutable();
	register_ai_library();
	GUI_SPIN_RANGE(uspex_gui._curr_step,1.,uspex_gui._tmp_num_opt_steps);
	if(uspex_gui._tmp_curr_step>i) uspex_gui._tmp_curr_step=(gdouble)i;
	GUI_SPIN_SET(uspex_gui._curr_step,uspex_gui._tmp_curr_step);
}
/***************************/
/* change the current step */
/***************************/
void spin_update_curr_step(void){
	gint i=(gint)uspex_gui._tmp_curr_step;
	gchar *text;
	/*update step information*/
	uspex_gui._tmp_isfixed=uspex_gui.calc._isfixed[i-1];
	if(uspex_gui._tmp_isfixed) GUI_TOGGLE_ON(uspex_gui._isfixed);
	else GUI_TOGGLE_OFF(uspex_gui._isfixed);
	GUI_COMBOBOX_SET(uspex_gui.abinitioCode,uspex_gui.calc.abinitioCode[i-1]);
	text=g_strdup_printf("%.4lf",uspex_gui.calc.KresolStart[i-1]);
	GUI_ENTRY_TEXT(uspex_gui.KresolStart,text);
	g_free(text);
	text=g_strdup_printf("%.4lf",uspex_gui.calc.vacuumSize[i-1]);
	GUI_ENTRY_TEXT(uspex_gui.vacuumSize,text);
	g_free(text);
	/*commandExecutable*/
	if(!uspex_gui.calc._isCmdList){
		/*only first command matter*/
		if(uspex_gui._tmp_commandExecutable[0]==NULL) text=g_strdup("N/A");
		else text=g_strdup(uspex_gui._tmp_commandExecutable[0]);
	}else{
		/*display ith command*/
		if(uspex_gui._tmp_commandExecutable[i-1]==NULL) text=g_strdup("N/A");
		else text=g_strdup(uspex_gui._tmp_commandExecutable[i-1]);
	}
	GUI_ENTRY_TEXT(uspex_gui.commandExecutable,text);
	g_free(text);
	/*uspex_gui._tmp_ai_input*/
	if(uspex_gui._tmp_ai_input[i-1]==NULL) text=g_strdup("N/A");
	else text=g_strdup(uspex_gui._tmp_ai_input[i-1]);
	GUI_ENTRY_TEXT(uspex_gui.ai_input,text);
	g_free(text);
	/*uspex_gui._tmp_ai_opt*/
	if(uspex_gui._tmp_ai_opt[i-1]==NULL) text=g_strdup("N/A");
	else text=g_strdup(uspex_gui._tmp_ai_opt[i-1]);
	GUI_ENTRY_TEXT(uspex_gui.ai_opt,text);
	g_free(text);
	/*uspex_gui._tmp_ai_lib_folder*/
	if(uspex_gui._tmp_ai_lib_folder[i-1]==NULL) text=g_strdup("N/A");
	else text=g_strdup(uspex_gui._tmp_ai_lib_folder[i-1]);
	GUI_ENTRY_TEXT(uspex_gui.ai_lib,text);
	g_free(text);
	/*uspex_gui.ai_lib_flavor*/
	sync_ai_lib_flavor();
	/*uspex_gui._tmp_ai_lib_sel*/
	sync_ai_lib_sel();
}
/********************/
/* toggle auto_step */
/********************/
void toggle_auto_step(void){
	if(uspex_gui.auto_step){
		GUI_LOCK(uspex_gui.ai_input);
		GUI_LOCK(uspex_gui.ai_input_button);
		GUI_LOCK(uspex_gui.ai_opt);
		GUI_LOCK(uspex_gui.ai_opt_button);
	}else{
		GUI_UNLOCK(uspex_gui.ai_input);
		GUI_UNLOCK(uspex_gui.ai_input_button);
		GUI_UNLOCK(uspex_gui.ai_opt);
		GUI_UNLOCK(uspex_gui.ai_opt_button);
	}
}
/********************************/
/* select AI code for this step */
/********************************/
void uspex_ai_selected(GUI_OBJ *w){
	gint idx;
	gint step=(gint)uspex_gui._tmp_curr_step;
	GUI_COMBOBOX_GET(w,uspex_gui.calc.abinitioCode[step-1]);
	/*propose some defaults command is missing*/
	if((uspex_gui._tmp_commandExecutable[step-1]==NULL)&&(uspex_gui.calc._isCmdList)){
		/*we don't have commandExecutable*/
		for(idx=0;idx<uspex_gui.calc._num_opt_steps;idx++){
			if(idx==(step-1)) continue;
			if(uspex_gui.calc.abinitioCode[idx]==uspex_gui.calc.abinitioCode[step-1])
				uspex_gui._tmp_commandExecutable[step-1]=g_strdup(uspex_gui._tmp_commandExecutable[idx]);
		}
		/*if all fail*/
		if(uspex_gui._tmp_commandExecutable[step-1]==NULL) uspex_gui._tmp_commandExecutable[step-1]=g_strdup("N/A");
	}
	uspex_gui.have_ZT=FALSE;
	for(idx=0;idx<uspex_gui.calc._num_opt_steps;idx++) uspex_gui.have_ZT|=(uspex_gui.calc.abinitioCode[idx]==14);
}
/****************************/
/* load input for this step */
/****************************/
void load_ai_input_dialog(void){
        GUI_OBJ *file_chooser;
        gint have_answer;
        gchar *filename;
        gchar *text;
	gint step=(gint)uspex_gui._tmp_curr_step;
        /**/
        GUI_PREPARE_OPEN_DIALOG(uspex_gui.window,file_chooser,"Select Specific INPUT_X file","*","INPUT_X");
        GUI_OPEN_DIALOG_RUN(file_chooser,have_answer,filename);
        if(have_answer){
                /*there should be no case where have_answer==TRUE AND filename==NULL, but just in case.*/
                if(filename) {
                        text=g_strdup_printf("%s",filename);
                        GUI_ENTRY_TEXT(uspex_gui.ai_input,text);
			if(uspex_gui._tmp_ai_input[step-1]!=NULL) g_free(uspex_gui._tmp_ai_input[step-1]);
			uspex_gui._tmp_ai_input[step-1]=NULL;
			uspex_gui._tmp_ai_input[step-1]=g_strdup(text);
                        g_free (filename);
			g_free(text);
                }
        }
        GUI_KILL_OPEN_DIALOG(file_chooser);
}
/*********************************************/
/* load option file (optional) for this step */
/*********************************************/
void load_ai_opt_dialog(void){
	GUI_OBJ *file_chooser;
	gint have_answer;
	gchar *filename;
	gchar *text;
	gint step=(gint)uspex_gui._tmp_curr_step;
	/**/
	GUI_PREPARE_OPEN_DIALOG(uspex_gui.window,file_chooser,"Select Specific OPTION_X file","*","OPTION_X");
	GUI_OPEN_DIALOG_RUN(file_chooser,have_answer,filename);
	if(have_answer){
		/*there should be no case where have_answer==TRUE AND filename==NULL, but just in case.*/
		if(filename) {
			text=g_strdup_printf("%s",filename);
			GUI_ENTRY_TEXT(uspex_gui.ai_opt,text);
			if(uspex_gui._tmp_ai_opt[step-1]!=NULL) g_free(uspex_gui._tmp_ai_opt[step-1]);
			uspex_gui._tmp_ai_opt[step-1]=NULL;
			uspex_gui._tmp_ai_opt[step-1]=g_strdup(text);
			g_free (filename);
			g_free(text);
		}
	}
	GUI_KILL_OPEN_DIALOG(file_chooser);
}
/*****************************************************/
/* define an ab initio code executable for this step */
/*****************************************************/
void load_abinitio_exe_dialog(void){
	GUI_OBJ *file_chooser;
	gint have_answer;
	gchar *filename;
	gchar *text;
	gint step=(gint)uspex_gui._tmp_curr_step;
	/**/
	GUI_PREPARE_OPEN_DIALOG(uspex_gui.window,file_chooser,"Select ab initio command","*","CommandExecutable");
	GUI_OPEN_DIALOG_RUN(file_chooser,have_answer,filename);
	if(have_answer){
		/*there should be no case where have_answer==TRUE AND filename==NULL, but just in case.*/
		if(filename) {
			text=g_strdup_printf("%s",filename);
			GUI_ENTRY_TEXT(uspex_gui.commandExecutable,text);
			if(uspex_gui._tmp_commandExecutable[step-1]!=NULL) g_free(uspex_gui._tmp_commandExecutable[step-1]);
			uspex_gui._tmp_commandExecutable[step-1]=NULL;
			uspex_gui._tmp_commandExecutable[step-1]=g_strdup(text);
			g_free (filename);
			g_free(text);
		}
	}
	GUI_KILL_OPEN_DIALOG(file_chooser);
}
/**********************************************************/
/* Open a (pseudo-)potential directory for this step/code */
/**********************************************************/
void load_ai_lib_dialog(void){
        GUI_OBJ *file_chooser;
        gint have_answer;
        gchar *filename;
        gchar *text;
	gint index=(gint)uspex_gui._tmp_curr_step-1;
	/* open folder */
        GUI_PREPARE_OPEN_FOLDER(uspex_gui.window,file_chooser,"Select the POTCAR folder");
        GUI_OPEN_DIALOG_RUN(file_chooser,have_answer,filename);
        if(have_answer){
                /*there should be no case where have_answer==TRUE AND filename==NULL, but just in case.*/
                if(filename) {
                        text=g_strdup_printf("%s",filename);
                        GUI_ENTRY_TEXT(uspex_gui.ai_lib,text);
                        g_free(text);
                        if(uspex_gui._tmp_ai_lib_folder[index]!=NULL) 
				g_free(uspex_gui._tmp_ai_lib_folder[index]);
                        uspex_gui._tmp_ai_lib_folder[index]=g_strdup_printf("%s",filename);
                        g_free (filename);
                        /*sync library flavor & sel*/
			sync_ai_lib_flavor();
			sync_ai_lib_sel();
                }
        }     
        GUI_KILL_OPEN_DIALOG(file_chooser);
}
/*************************************/
/* Apply (register) a library flavor */
/*************************************/
void apply_ai_lib_flavor(void){
	gint step=(gint)uspex_gui._tmp_curr_step-1;
	gint code,idx;
	gboolean is_ok;
	gchar *text;
	gchar *sel,*tmp,*tmp2;
	/**/
	GUI_COMBOBOX_GET_TEXT(uspex_gui.ai_lib_flavor,text);
	if(text==NULL) return;/*no entry*/
	code=uspex_gui.calc.abinitioCode[step];
	if(uspex_gui._tmp_ai_lib_sel[step]==NULL){
		/*depending on case, prepare a dummy selection of the flavors*/
		switch(code){
		case 1://VASP
		case 2://SIESTA
		case 8://Q.E.
		case 10://ATK
		case 11://CASTEP
		case 18://Abinit
			/*prepare a default species string*/
			sel=g_strdup(elements[uspex_gui.calc.atomType[0]].symbol);
			for(idx=1;idx<uspex_gui.calc._nspecies;idx++){
				tmp=g_strdup_printf("%s %s",sel,elements[uspex_gui.calc.atomType[idx]].symbol);
				g_free(sel);sel=tmp;
			}
			uspex_gui._tmp_ai_lib_sel[step]=sel;
			break;
		case 3://GULP
		case 12://Tinker
		case 7://CP2K
		case 15://DFTB
		case 4://LAMMPS
		case 13://MOPAC
		case 14://BoltzTraP
		case 16://Gaussian
		case 17://N/A
		case 19://CRYSTAL
		case 5://ORCA
		case 6://DMACRYS
		case 9://FHI-aims
		default:
			/*all other don't need such preparation*/
			break;
		}
	}
	/*now update uspex_gui._tmp_ai_lib_sel[step]*/
	switch(code){
	case 1://VASP
	case 2://SIESTA
	case 8://Q.E.
	case 10://ATK
	case 11://CASTEP
	case 18://Abinit
		/*change symbol chain with flavor*/
		sel=uspex_gui._tmp_ai_lib_sel[step];
		idx=0;
		is_ok=TRUE;
		do{
			while((sel[idx]!=text[0])&&(sel[idx]!='\0')) idx++;
			if(sel[idx]=='\0'){
				fprintf(stderr,"USPEX: species flavor not found!\n");
				return;/*this *is* bad, species was never found!*/
			}
			if(sel[idx+1]==text[1]) is_ok=FALSE;
			else{
				if((!g_ascii_islower (sel[idx+1]))
				 &&(!g_ascii_islower (text[1])))
					is_ok=FALSE;
				else idx++;
			}
		}while(is_ok);
		/*idx points to the first letter of species match*/
		if(idx==0){
			/*species chain is at start*/
			while((!g_ascii_isspace(sel[idx]))&&(sel[idx]!='\0')) idx++;
			if(sel[idx]=='\0'){
				/*only one species?*/
				g_free(uspex_gui._tmp_ai_lib_sel[step]);
				uspex_gui._tmp_ai_lib_sel[step]=g_strdup(text);
			}else{
				tmp=g_strdup_printf("%s %s",text,&(sel[idx+1]));
				g_free(uspex_gui._tmp_ai_lib_sel[step]);
				uspex_gui._tmp_ai_lib_sel[step]=tmp;
			}
		}else{
			/*species chain in the middle (or end)*/
			tmp=g_strndup(sel,idx-1);
			tmp2=g_strdup_printf("%s %s",tmp,text);
			g_free(tmp);tmp=tmp2;
			while((!g_ascii_isspace(sel[idx]))&&(sel[idx]!='\0')) idx++;
			if(sel[idx]!='\0'){
				/*add the remaining chain*/
				tmp2=g_strdup_printf("%s %s",tmp,&(sel[idx+1]));
				g_free(tmp);tmp=tmp2;
			}
			g_free(uspex_gui._tmp_ai_lib_sel[step]);
			uspex_gui._tmp_ai_lib_sel[step]=tmp;
		}
		break;
	case 3://GULP
	case 12://Tinker
	case 7://CP2K
	case 15://DFTB
	case 4://LAMMPS
	case 13://MOPAC
	case 14://BoltzTraP
	case 16://Gaussian
	case 17://N/A
	case 19://CRYSTAL
	case 5://ORCA
	case 6://DMACRYS
	case 9://FHI-aims
	default:
		/*something was selected, so put it here*/
		g_free(uspex_gui._tmp_ai_lib_sel[step]);
		uspex_gui._tmp_ai_lib_sel[step]=g_strdup(text);
		break;
	}
	g_free(text);
	sync_ai_lib_sel();
}
/**************************/
/* Apply step information */
/**************************/
void apply_step(void){
	gint i=(gint)uspex_gui._tmp_curr_step;
	/*REG step values*/
	uspex_gui.calc._isfixed[i-1]=uspex_gui._tmp_isfixed;
	GUI_COMBOBOX_GET(uspex_gui.abinitioCode,uspex_gui.calc.abinitioCode[i-1]);
	GUI_REG_VAL(uspex_gui.KresolStart,uspex_gui.calc.KresolStart[i-1],"%lf");
	GUI_REG_VAL(uspex_gui.vacuumSize,uspex_gui.calc.vacuumSize[i-1],"%lf");
	GUI_ENTRY_GET_TEXT(uspex_gui.commandExecutable,uspex_gui._tmp_commandExecutable[i-1]);
	register_commandExecutable();/*REG command*/
	apply_ai_lib_flavor();
	register_ai_library();
}
/******************************************************************/
/* sets a remote folder (in case it can be set by a local folder) */
/******************************************************************/
void load_remote_folder(void){
	GUI_OBJ *file_chooser;
	gint have_answer;
	gchar *filename;
	gchar *text;
	/**/
	GUI_PREPARE_OPEN_FOLDER(uspex_gui.window,file_chooser,"Select identical remote folder");
	GUI_OPEN_DIALOG_RUN(file_chooser,have_answer,filename);
	if(have_answer){
		/*there should be no case where have_answer==TRUE AND filename==NULL, but just in case.*/
		if(filename) {
			text=g_strdup_printf("%s",filename);
			GUI_ENTRY_TEXT(uspex_gui.remoteFolder,text);
			g_free(text);
			USPEX_REG_TEXT(remoteFolder);
			g_free (filename);
		}
	}
	GUI_KILL_OPEN_DIALOG(file_chooser);
}
/************************/
/* toggle USPEX v. 10.3 */
/************************/
void toggle_v1030(void){
	if(uspex_gui.have_v1030) {
		GUI_LOCK(uspex_gui.sel_octave);
		GUI_TOGGLE_ON(uspex_gui.sel_v1030_1);
		GUI_TOGGLE_ON(uspex_gui.sel_v1030_3);
		uspex_gui.have_new_optType=TRUE;
		/*old opt style is lost anyway*/
		GUI_TOGGLE_ON(uspex_gui.sel_new_opt);
		GUI_LOCK(uspex_gui.sel_new_opt);
	}else {
		GUI_UNLOCK(uspex_gui.sel_octave);
		GUI_TOGGLE_OFF(uspex_gui.sel_v1030_1);
		GUI_TOGGLE_OFF(uspex_gui.sel_v1030_3);
		uspex_gui.have_new_optType=FALSE;
		/*some previous version can have both*/
		GUI_UNLOCK(uspex_gui.sel_new_opt);
		/*this is the default anyway*/
		GUI_TOGGLE_OFF(uspex_gui.sel_new_opt);
	}
	opt_toggle();
}
/**************************/
/* select dynamicalBestHM */
/**************************/
void uspex_dyn_HM_selected(GUI_OBJ *w){
	gint index;
	GUI_COMBOBOX_GET(w,index);
	switch (index){
	case 0:
		uspex_gui.calc.dynamicalBestHM=0;
		break;
	case 1:
		uspex_gui.calc.dynamicalBestHM=1;
		break;
	case 2:
		/* fall through */
	default:
		uspex_gui.calc.dynamicalBestHM=2;
	}
}
/**********************/
/* select manyParents */
/**********************/
void uspex_manyParents_selected(GUI_OBJ *w){
	gint index;
	GUI_COMBOBOX_GET(w,index);
	switch (index){
	case 1:
		uspex_gui.calc.manyParents=1;
		break;
	case 2:
		uspex_gui.calc.manyParents=2;
		break;
	case 3:
		uspex_gui.calc.manyParents=3;
		break;
	case 0:
		/* fall through */
	default:
		uspex_gui.calc.manyParents=0;
	}
}
/******************/
/* select TE_goal */
/******************/
void uspex_TE_goal_selected(GUI_OBJ *w){
	gint index;
	GUI_COMBOBOX_GET(w,index);
	switch (index){
	case 1:
		uspex_gui.calc.TE_goal=US_BT_ZTxx;
		break;
	case 2:
		uspex_gui.calc.TE_goal=US_BT_ZTyy;
		break;
	case 3:
		uspex_gui.calc.TE_goal=US_BT_ZTzz;
		break;
	case 0:
		/* fall through */
	default:
		uspex_gui.calc.TE_goal=US_BT_ZT;
	}
}
/*************************************/
/* load the BoltzTraP script/command */
/*************************************/
void load_cmd_BoltzTraP(void){
	GUI_OBJ *file_chooser;
	gint have_answer;
	gchar *filename;
	gchar *text;
	/**/
	GUI_PREPARE_OPEN_DIALOG(uspex_gui.window,file_chooser,"Select BoltzTrap command","*","BoltzTraP");
	GUI_OPEN_DIALOG_RUN(file_chooser,have_answer,filename);
	if(have_answer){
		/*there should be no case where have_answer==TRUE AND filename==NULL, but just in case.*/
		if(filename) {
			text=g_strdup_printf("%s",filename);
			GUI_ENTRY_TEXT(uspex_gui.cmd_BoltzTraP,text);
			if(uspex_gui._tmp_cmd_BoltzTraP!=NULL) g_free(uspex_gui._tmp_cmd_BoltzTraP);
			uspex_gui._tmp_cmd_BoltzTraP=NULL;
			uspex_gui._tmp_cmd_BoltzTraP=g_strdup(text);
			g_free (filename);
			g_free(text);
		}
	}
	GUI_KILL_OPEN_DIALOG(file_chooser);
}
/*******************************/
/* load orderParameter command */
/*******************************/
void load_cmd_OP(void){
	GUI_OBJ *file_chooser;
	gint have_answer;
	gchar *filename;
	gchar *text;
	/**/
	GUI_PREPARE_OPEN_DIALOG(uspex_gui.window,file_chooser,"Select OrderParameter command","*","script");
	GUI_OPEN_DIALOG_RUN(file_chooser,have_answer,filename);
	if(have_answer){
		/*there should be no case where have_answer==TRUE AND filename==NULL, but just in case.*/
		if(filename) {
			text=g_strdup_printf("%s",filename);
			GUI_ENTRY_TEXT(uspex_gui.cmdOrderParameter,text);
			g_free(text);
			USPEX_REG_TEXT(cmdOrderParameter);
			g_free (filename);
		}
	}
	GUI_KILL_OPEN_DIALOG(file_chooser);
}
/************************************/
/* load enthalpyTemperature command */
/************************************/
void load_cmd_ET(void){
	GUI_OBJ *file_chooser;
	gint have_answer;
	gchar *filename;
	gchar *text;
	/**/
	GUI_PREPARE_OPEN_DIALOG(uspex_gui.window,file_chooser,"Select EnthalpyTemperature command","*","script");
	GUI_OPEN_DIALOG_RUN(file_chooser,have_answer,filename);
	if(have_answer){
		/*there should be no case where have_answer==TRUE AND filename==NULL, but just in case.*/
		if(filename) {
			text=g_strdup_printf("%s",filename);
			GUI_ENTRY_TEXT(uspex_gui.cmdEnthalpyTemperature,text);
			g_free(text);
			USPEX_REG_TEXT(cmdEnthalpyTemperature);
			g_free (filename);
		}
	}
	GUI_KILL_OPEN_DIALOG(file_chooser);
}
/***************************/
/*load orderParameter file */
/***************************/
void load_OP_file(void){
	GUI_OBJ *file_chooser;
	gint have_answer;
	gchar *filename;
	gchar *text;
	/**/
	GUI_PREPARE_OPEN_DIALOG(uspex_gui.window,file_chooser,"Select OP file","*.dat","data file");
	GUI_OPEN_DIALOG_RUN(file_chooser,have_answer,filename);
	if(have_answer){
		/*there should be no case where have_answer==TRUE AND filename==NULL, but just in case.*/
		if(filename) {
			text=g_strdup_printf("%s",filename);
			GUI_ENTRY_TEXT(uspex_gui.orderParameterFile,text);
			g_free(text);
			USPEX_REG_TEXT(orderParameterFile);
			g_free (filename);
		}
	}
	GUI_KILL_OPEN_DIALOG(file_chooser);
}
/********************************/
/*load enthalpyTemperature file */
/********************************/
void load_ET_file(){
	GUI_OBJ *file_chooser;
	gint have_answer;
	gchar *filename;
	gchar *text;
	/**/
	GUI_PREPARE_OPEN_DIALOG(uspex_gui.window,file_chooser,"Select ET file","*.dat","data file");
	GUI_OPEN_DIALOG_RUN(file_chooser,have_answer,filename);
	if(have_answer){
		/*there should be no case where have_answer==TRUE AND filename==NULL, but just in case.*/
		if(filename) {
			text=g_strdup_printf("%s",filename);
			GUI_ENTRY_TEXT(uspex_gui.enthalpyTemperatureFile,text);
			g_free(text);
			USPEX_REG_TEXT(enthalpyTemperatureFile);
			g_free (filename);
		}
	}
	GUI_KILL_OPEN_DIALOG(file_chooser);
}
/************************/
/* load trajectory file */
/************************/
void load_traj_file(){
	GUI_OBJ *file_chooser;
	gint have_answer;
	gchar *filename;
	gchar *text;
	/**/
	GUI_PREPARE_OPEN_DIALOG(uspex_gui.window,file_chooser,"Select MD trajectory file","*","traj file");
	GUI_OPEN_DIALOG_RUN(file_chooser,have_answer,filename);
	if(have_answer){
		/*there should be no case where have_answer==TRUE AND filename==NULL, but just in case.*/
		if(filename) {
			text=g_strdup_printf("%s",filename);
			GUI_ENTRY_TEXT(uspex_gui.trajectoryFile,text);
			g_free(text);
			USPEX_REG_TEXT(trajectoryFile);
			g_free (filename);
		}
	}
	GUI_KILL_OPEN_DIALOG(file_chooser);
}
/************************/
/* load MD restart file */
/************************/
void load_MDrestart_file(){
	GUI_OBJ *file_chooser;
	gint have_answer;
	gchar *filename;
	gchar *text;
	/**/
	GUI_PREPARE_OPEN_DIALOG(uspex_gui.window,file_chooser,"Select MD restart file","*","restart file");
	GUI_OPEN_DIALOG_RUN(file_chooser,have_answer,filename);
	if(have_answer){
		/*there should be no case where have_answer==TRUE AND filename==NULL, but just in case.*/
		if(filename) {
			text=g_strdup_printf("%s",filename);
			GUI_ENTRY_TEXT(uspex_gui.MDrestartFile,text);
			g_free(text);
			USPEX_REG_TEXT(MDrestartFile);
			g_free (filename);
		}
	}
	GUI_KILL_OPEN_DIALOG(file_chooser);
}
/********************/
/* select FullRelax */
/********************/
void uspex_relax_selected(GUI_OBJ *w){
	gint index;
	GUI_COMBOBOX_GET(w,index);
	switch (index){
	case 0:
		uspex_gui.calc.FullRelax=0;
		break;
	case 1:
		uspex_gui.calc.FullRelax=1;
		break;
	case 2:
		/* fall through */
	default:
		uspex_gui.calc.FullRelax=2;
	}
}
/********************************************/
/* init metadynamics model from GDIS models */
/********************************************/
void init_metadynamics_models(){
	gint index;
	gchar *text;
	GSList *list;
	struct model_pak *data;
	/*WIPE MODEL*/
	GUI_COMBOBOX_WIPE(uspex_gui.meta_model);
	GUI_COMBOBOX_ADD(uspex_gui.meta_model,"From VASP5 POSCAR file");
	index=1;
	for(list=sysenv.mal;list;list=g_slist_next(list)){
		data = list->data;
		/*to be registered a model must have at least 1 atom*/
		if(g_slist_length(data->cores)>0) {
			text=g_strdup_printf("%i: %s",index,data->basename);
			GUI_COMBOBOX_ADD(uspex_gui.meta_model,text);
			g_free(text);
		}
		index++;
	}
	GUI_COMBOBOX_SET(uspex_gui.meta_model,0);
}
/****************************************/
/* load metadynamics starting structure */
/****************************************/
void load_meta_start_file(void){
	GUI_OBJ *file_chooser;
	gint have_answer;
	gchar *filename;
	gchar *text;
	/**/
	GUI_PREPARE_OPEN_DIALOG(uspex_gui.window,file_chooser,"Select META starting point","POSCAR*","POSCAR_1");
	GUI_OPEN_DIALOG_RUN(file_chooser,have_answer,filename);
	if(have_answer){
	/*there should be no case where have_answer==TRUE AND filename==NULL, but just in case.*/
		if(filename) {
			text=g_strdup_printf("%s",filename);
			GUI_COMBOBOX_ADD_TEXT(uspex_gui.meta_model,0,text);
			GUI_COMBOBOX_DEL(uspex_gui.meta_model,1);
			GUI_COMBOBOX_SET(uspex_gui.meta_model,0);
			g_free(text);
			g_free (filename);
		}
	}
	GUI_KILL_OPEN_DIALOG(file_chooser);
}
/*****************************/
/* select metadynamics model */
/*****************************/
void uspex_meta_model_selected(GUI_OBJ *w){
	gint index;
	GUI_COMBOBOX_GET(w,index);
	switch (index){
	case 0:
		GUI_UNLOCK(uspex_gui.meta_model_button);
		break;
	default:
		GUI_LOCK(uspex_gui.meta_model_button);
	}
}
/*********************/
/* update VCNEB Type */
/*********************/
void update_vcnebType(void){
	gchar *text;
	gint type=uspex_gui.calc._vcnebtype_method*100;
	if(uspex_gui.calc.calculationMethod!=US_CM_VCNEB) return;
	if(uspex_gui.calc._vcnebtype_img_num) type+=10;
	if(uspex_gui.calc._vcnebtype_spring) type+=1;
	uspex_gui.calc.vcnebType=type;
	/*update uspex_gui.vcnebType */
	text=g_strdup_printf("%3i",type);
	GUI_ENTRY_TEXT(uspex_gui.vcnebType,text);
	g_free(text);
	/*consequences*/
	if(uspex_gui.calc._vcnebtype_img_num) GUI_UNLOCK(uspex_gui.VarPathLength);
	else GUI_LOCK(uspex_gui.VarPathLength);
	if(uspex_gui.calc._vcnebtype_spring){
		GUI_UNLOCK(uspex_gui.K_min);
		GUI_UNLOCK(uspex_gui.K_max);
		GUI_LOCK(uspex_gui.Kconstant);
	}else{
		GUI_LOCK(uspex_gui.K_min);
		GUI_LOCK(uspex_gui.K_max);
		GUI_UNLOCK(uspex_gui.Kconstant);
	}
}
/***********************/
/* select VCNEB method */
/***********************/
void uspex_vcneb_method_selected(GUI_OBJ *w){
	gint index;
	GUI_COMBOBOX_GET(w,index);
	switch (index){
	case 1:
		uspex_gui.calc._vcnebtype_method=2;
		break;
	case 0:
		/* fall through */
	default:
		uspex_gui.calc._vcnebtype_method=1;
	}
	update_vcnebType();
}
/******************************/
/* select VCNEB optReadImages */
/******************************/
void uspex_ReadImg_selected(GUI_OBJ *w){
	gint index;
	GUI_COMBOBOX_GET(w,index);
	switch (index){
	case 0:
		uspex_gui.calc.optReadImages=0;
		break;
	case 1:
		uspex_gui.calc.optReadImages=1;
		break;
	case 2:
		/* fall through */
	default:
		uspex_gui.calc.optReadImages=2;
	}
}
/*******************************/
/* select VCNEB optimizer type */
/*******************************/
void uspex_optimizerType_selected(GUI_OBJ *w){
	gint index;
	GUI_COMBOBOX_GET(w,index);
	switch (index){
	case 1:
		uspex_gui.calc.optimizerType=2;
	case 0:
		/* fall through */
	default:
		uspex_gui.calc.optimizerType=1;
	}
}
/********************************/
/* select VCNEB relaxation type */
/********************************/
void uspex_RelaxType_selected(GUI_OBJ *w){
	gint index;
	GUI_COMBOBOX_GET(w,index);
	switch (index){
	case 0:
		uspex_gui.calc.optRelaxType=1;
		break;
	case 1:
		uspex_gui.calc.optRelaxType=2;
		break;
	case 2:
		/* fall through */
	default:
		uspex_gui.calc.optRelaxType=3;
	}
}
/*****************************/
/* select VCNEB CI/DI method */
/*****************************/
void uspex_CIDI_selected(GUI_OBJ *w){
	gint index;
	GUI_COMBOBOX_GET(w,index);
	switch (index){
	case 1:
		uspex_gui.calc.optMethodCIDI=1;
		break;
	case 2:
		uspex_gui.calc.optMethodCIDI=-1;
		break;
	case 3:
		uspex_gui.calc.optMethodCIDI=2;
		break;
	case 0:
		/* fall through */
	default:
		uspex_gui.calc.optMethodCIDI=0;
	}
	/*consequences*/
	if(uspex_gui.calc.optMethodCIDI==0){
		GUI_LOCK(uspex_gui.startCIDIStep);
		GUI_LOCK(uspex_gui.pickupImages);
	}else{
		GUI_UNLOCK(uspex_gui.startCIDIStep);
		GUI_UNLOCK(uspex_gui.pickupImages);
	}
}
/***********************************/
/* select VCNEB PATH output format */
/***********************************/
void uspex_FormatType_selected(GUI_OBJ *w){
	gint index;
	GUI_COMBOBOX_GET(w,index);
	switch (index){
	case 0:
		uspex_gui.calc.FormatType=1;
		break;
	case 2:
		uspex_gui.calc.FormatType=3;
		break;
	case 1:
		/* fall through */
	default:
		uspex_gui.calc.FormatType=2;
	}
}
/*************************************/
/* initialize image models from GDIS */
/*************************************/
void init_img_gdis_model(void){
	gint index;
	gchar *text;
	GSList *list;
	struct model_pak *data;
	/*WIPE MODEL*/
	GUI_COMBOBOX_WIPE(uspex_gui.img_model);
	GUI_COMBOBOX_ADD(uspex_gui.img_model,"From VASP5 POSCAR file");
	index=1;
	for(list=sysenv.mal;list;list=g_slist_next(list)){
		data = list->data;
		/*to be registered, a model must 3D periodic... I think*/
		if(data->periodic==3){
			text=g_strdup_printf("%i: %s",index,data->basename);
			GUI_COMBOBOX_ADD(uspex_gui.img_model,text);
			g_free(text);
		}
		index++;
	}
	GUI_COMBOBOX_SET(uspex_gui.img_model,0);
}
/*****************************/
/* select VCNEB Images model */
/*****************************/
void uspex_img_model_selected(GUI_OBJ *w){
	gint index;
	GUI_COMBOBOX_GET(w,index);
	switch (index){
	case 0:
		GUI_UNLOCK(uspex_gui.img_model_button);
		break;
	default:
		GUI_LOCK(uspex_gui.img_model_button);
	}
}
/*************************/
/* load VCNEB Image file */
/*************************/
void load_img_model_file(void){
	GUI_OBJ *file_chooser;
	gint have_answer;
	gchar *filename;
	gchar *text;
	/**/
	GUI_PREPARE_OPEN_DIALOG(uspex_gui.window,file_chooser,"Select VCNEB image model","*","VASP5 format");
	GUI_OPEN_DIALOG_RUN(file_chooser,have_answer,filename);
	if(have_answer){
		/*there should be no case where have_answer==TRUE AND filename==NULL, but just in case.*/
		if(filename) {
			text=g_strdup_printf("%s",filename);
			GUI_COMBOBOX_ADD_TEXT(uspex_gui.img_model,0,text);
			GUI_COMBOBOX_DEL(uspex_gui.img_model,1);
			GUI_COMBOBOX_SET(uspex_gui.img_model,0);
			g_free(text);
			g_free (filename);
		}
	}
	GUI_KILL_OPEN_DIALOG(file_chooser);
}
/*****************************************/
/* initialize molecular models from GDIS */
/*****************************************/
void init_uspex_gdis_mol(void){
	gint index;
	gchar *text;
	GSList *list;
	struct model_pak *data;
	/*WIPE MODEL*/
	GUI_COMBOBOX_WIPE(uspex_gui.mol_gdis);
	index=0;
	for(list=sysenv.mal;list;list=g_slist_next(list)){
		data = list->data;
		/*to be registered, a model must have a molecule*/
		if(g_slist_length(data->moles)>0) {
			text=g_strdup_printf("%i: %s",index,data->basename);
			GUI_COMBOBOX_ADD(uspex_gui.mol_gdis,text);
			g_free(text);
		}
		index++;
	}
	GUI_COMBOBOX_SET(uspex_gui.mol_gdis,uspex_gui._tmp_mols_gdis[0]);
}
/**************************/
/* select Molecular model */
/**************************/
void uspex_mol_model_selected(GUI_OBJ *w){
	gint index;
	GUI_COMBOBOX_GET(w,index);
	switch (index){
	case 1:/*from folder*/
		GUI_UNLOCK(uspex_gui.num_mol);
		GUI_UNLOCK(uspex_gui.mol_model_button);
		GUI_LOCK(uspex_gui.mol_gdis);
		GUI_LOCK(uspex_gui.curr_mol);
		GUI_LOCK(uspex_gui.mol_gulp);
		GUI_LOCK(uspex_gui.mol_apply_button);
		break;
	case 2:/*from gdis*/
		GUI_UNLOCK(uspex_gui.num_mol);
		GUI_LOCK(uspex_gui.mol_model_button);
		GUI_UNLOCK(uspex_gui.mol_gdis);
		GUI_UNLOCK(uspex_gui.curr_mol);
		GUI_UNLOCK(uspex_gui.mol_gulp);
		GUI_UNLOCK(uspex_gui.mol_apply_button);
		break;
	case 0:/*already provided*/
		/* fall through */
	default:
		GUI_LOCK(uspex_gui.num_mol);
		GUI_LOCK(uspex_gui.mol_model_button);
		GUI_LOCK(uspex_gui.mol_gdis);
		GUI_LOCK(uspex_gui.curr_mol);
		GUI_LOCK(uspex_gui.mol_gulp);
		GUI_LOCK(uspex_gui.mol_apply_button);
	}
}
/*******************************/
/* load molecular model/folder */
/*******************************/
void load_mol_model_folder(void){
	GUI_OBJ *file_chooser;
	gint have_answer;
	gchar *filename;
	gchar *text;
	/**/
	GUI_PREPARE_OPEN_FOLDER(uspex_gui.window,file_chooser,"Select MOL_x folder");
	GUI_OPEN_DIALOG_RUN(file_chooser,have_answer,filename);
	if(have_answer){
		/*there should be no case where have_answer==TRUE AND filename==NULL, but just in case.*/
		if(filename) {
			text=g_strdup_printf("%s",filename);
			/*save that path as mol_model #1*/
			GUI_COMBOBOX_ADD_TEXT(uspex_gui.mol_model,1,text);
			GUI_COMBOBOX_DEL(uspex_gui.mol_model,2);
			GUI_COMBOBOX_SET(uspex_gui.mol_model,1);
			g_free(text);
			g_free(filename);
		}
	}
	GUI_KILL_OPEN_DIALOG(file_chooser);
}
/******************************/
/* update number of molecules */
/******************************/
void spin_update_num_mol(void){
	gint idx;
	gint i=uspex_gui._tmp_num_mol;
	gint *tmp_i;
	gboolean *tmp_b;
	/**/
	tmp_i=g_malloc(i*sizeof(gint));
	tmp_b=g_malloc(i*sizeof(gboolean));
	if(i<uspex_gui.calc._nmolecules){
		/*reduction*/
		for(idx=0;idx<i;idx++) {
			tmp_i[idx]=uspex_gui._tmp_mols_gdis[idx];
			tmp_b[idx]=uspex_gui._tmp_mols_gulp[idx];
		}
	}else{
		/*augmentation*/
		for(idx=0;idx<uspex_gui.calc._nmolecules;idx++) {
			tmp_i[idx]=uspex_gui._tmp_mols_gdis[idx];
			tmp_b[idx]=uspex_gui._tmp_mols_gulp[idx];
		}
		/*in case we increase more than one time <- should not happen*/
		for(idx=uspex_gui.calc._nmolecules;idx<i;idx++){
			tmp_i[idx]=0;
			tmp_b[idx]=FALSE;
		}
	}
	g_free(uspex_gui._tmp_mols_gdis);
	uspex_gui._tmp_mols_gdis=tmp_i;
	g_free(uspex_gui._tmp_mols_gulp);
	uspex_gui._tmp_mols_gulp=tmp_b;
	uspex_gui.calc._nmolecules=i;
	GUI_SPIN_RANGE(uspex_gui.curr_mol,1.,(gdouble)i);
}
/***********************/
/* Select a GDIS model */
/***********************/
void uspex_gdis_mol_selected(GUI_OBJ *w){
	/*nothing here for now*/
}
/***************************/
/* update current molecule */
/***************************/
void spin_update_curr_mol(void){
	gint i=uspex_gui._tmp_curr_mol;
	GUI_COMBOBOX_SET(uspex_gui.mol_gdis,uspex_gui._tmp_mols_gdis[i]);/*_BUG_ memory*/
	if(uspex_gui._tmp_mols_gulp[i]){
		GUI_TOGGLE_ON(uspex_gui.mol_gulp);
	}else{
		GUI_TOGGLE_OFF(uspex_gui.mol_gulp);
	}
}
/*****************/
/* apply gdis mol*/
/*****************/
void apply_gdis_mol(void){
	gint i=uspex_gui._tmp_curr_mol;
	gint index;
	GUI_COMBOBOX_GET(uspex_gui.mol_gdis,index);
	uspex_gui._tmp_mols_gdis[i]=index;
	uspex_gui._tmp_mols_gulp[i]=uspex_gui.mol_as_gulp;
}
/*****************************************/
/* init substrate model from GDIS models */
/*****************************************/
void init_substrate_models(){
	gint index;
	gchar *text;
	GSList *list;
	struct model_pak *data;
	/*WIPE MODEL*/
	GUI_COMBOBOX_WIPE(uspex_gui.substrate_model);
	GUI_COMBOBOX_ADD(uspex_gui.substrate_model,"From VASP5 POSCAR file");
	index=1;
	for(list=sysenv.mal;list;list=g_slist_next(list)){
		data = list->data;
		/*to be registered, a model must be 3D periodic*/
		if(data->periodic==3){
			text=g_strdup_printf("%i: %s",index,data->basename);
			GUI_COMBOBOX_ADD(uspex_gui.substrate_model,text);
			g_free(text);
		}
		index++;
	}
	GUI_COMBOBOX_SET(uspex_gui.substrate_model,0);
}
/**********************************/
/* select surface substrate model */
/**********************************/
void uspex_substrate_model_selected(GUI_OBJ *w){
	gint index;
	GUI_COMBOBOX_GET(w,index);
	switch (index){
	case 0:
		GUI_UNLOCK(uspex_gui.substrate_model_button);
		break;
	default:
		GUI_LOCK(uspex_gui.substrate_model_button);
	}
}
/********************************/
/* load surface substrate model */
/********************************/
void load_substrate_model_file(void){
	GUI_OBJ *file_chooser;
	gint have_answer;
	gchar *filename;
	gchar *text;
	/**/
	GUI_PREPARE_OPEN_DIALOG(uspex_gui.window,file_chooser,"Select Substrate structure","POSCAR_SUBSTRATE","VASP5 format");
	GUI_OPEN_DIALOG_RUN(file_chooser,have_answer,filename);
	if(have_answer){
		/*there should be no case where have_answer==TRUE AND filename==NULL, but just in case.*/
		if(filename) {
			text=g_strdup_printf("%s",filename);
			GUI_COMBOBOX_ADD_TEXT(uspex_gui.substrate_model,0,text);
			GUI_COMBOBOX_DEL(uspex_gui.substrate_model,1);
			GUI_COMBOBOX_SET(uspex_gui.substrate_model,0);
			g_free(text);
			g_free (filename);
		}
	}
	GUI_KILL_OPEN_DIALOG(file_chooser);
}
/************************/
/* Switch notebook page */
/************************/
void uspex_gui_page_switch(GUI_NOTE *notebook,GUI_OBJ *page,guint page_num){
	gint from_page;
	gint idx;
	/* some (few) values need to be updated from a page to another*/
	GUI_NOTE_PAGE_NUMBER(notebook,from_page);
	if(from_page==(gint)page_num) return;/*not moved*/
	uspex_gui.cur_page=page_num;
	GUI_LOCK(page);
	if(page_num==USPEX_PAGE_SYSTEM){
		/*sync uspex_gui._calctype_mag*/
		if(uspex_gui.calc._calctype_mag) {
			GUI_TOGGLE_ON(uspex_gui._calctype_mag);
		}else{
			GUI_TOGGLE_OFF(uspex_gui._calctype_mag);
		}
		auto_bond_toggle();
		toggle_auto_C_lat();
		toggle_auto_C_ion();
	} else if(page_num==USPEX_PAGE_STRUCTURES){
		/*sync uspex_gui._calctype_mag_2*/
		if(uspex_gui.calc._calctype_mag) {
			GUI_TOGGLE_ON(uspex_gui._calctype_mag_2);
		}else{
			GUI_TOGGLE_OFF(uspex_gui._calctype_mag_2);
		}
	} else if (page_num==USPEX_PAGE_CALCULATION){
		/**/
		toggle_auto_step();
		toggle_v1030();
	} else if (page_num==USPEX_PAGE_ADVANCED){
		/**/
	} else if (page_num==USPEX_PAGE_SPECIFIC){
		/*update meta_model_button*/
		GUI_COMBOBOX_GET(uspex_gui.meta_model,idx);
		if((uspex_gui.calc.calculationMethod==US_CM_META)&&(idx==0)) GUI_UNLOCK(uspex_gui.meta_model_button);
		else GUI_LOCK(uspex_gui.meta_model_button);
		/*update vcneb*/
		update_vcnebType();
	}
	GUI_UNLOCK(page);
}
/****************************************************/
/* convert current uspex_gui into uspex_calc_struct */
/****************************************************/
void uspex_gui_sync(){
	gint index;
	gint idx;
	gint jdx;
	gchar *text;
	gchar *ptr;
	gchar *ptr2;
	gint num;
	gint val;
	/* sync start here */
	USPEX_REG_TEXT(name);
	//calculationMethod is already sync
	//calculationType is already sync
	uspex_gui.calc._calctype_dim=(gint)uspex_gui._dim;
	//_calctype_mol is already sync
	//_calctype_var is already sync
	//_calctype_mag is already sync
	//optType is already sync
	//have_new_optType is already sync
	if(uspex_gui.have_new_optType) USPEX_REG_TEXT(new_optType);
	//anti_opt is already sync
	//n_species is already sync
	//atomType is already sync
	if(uspex_gui.calc._calctype_var){/*VARCOMP -> take numSpecies from numSpecies lines*/
		GUI_COMBOBOX_GET(uspex_gui.numSpecies,index);/*PUSH index*/
		/*count _var_nspecies*/
		idx=0;
		GUI_COMBOBOX_SET(uspex_gui.numSpecies,0);
		GUI_COMBOBOX_GET_TEXT(uspex_gui.numSpecies,text);
		while(g_ascii_strcasecmp(text,"ADD SPECIES BLOCK") != 0){
			g_free(text);
			idx++;
			GUI_COMBOBOX_SET(uspex_gui.numSpecies,idx);
			GUI_COMBOBOX_GET_TEXT(uspex_gui.numSpecies,text);
		}
		uspex_gui.calc._var_nspecies=idx;
		if(uspex_gui.calc.numSpecies!=NULL) g_free(uspex_gui.calc.numSpecies);
		uspex_gui.calc.numSpecies=g_malloc((uspex_gui.calc._nspecies*uspex_gui.calc._var_nspecies)*sizeof(gint));
		jdx=0;
		GUI_COMBOBOX_SET(uspex_gui.numSpecies,0);
		GUI_COMBOBOX_GET_TEXT(uspex_gui.numSpecies,text);
		while(g_ascii_strcasecmp(text,"ADD SPECIES BLOCK") != 0){
			ptr=text;
			while((*ptr!='\0')&&(!g_ascii_isgraph(*ptr))) ptr++;
			idx=0;
			while((idx<uspex_gui.calc._nspecies)&&(*ptr!='\0')){
				uspex_gui.calc.numSpecies[idx+jdx*uspex_gui.calc._nspecies]=(gint)g_ascii_strtoull(ptr,&ptr2,10);
				if(ptr2==NULL) break;/*_BUG_ fix #1*/
				if(*ptr2=='\0') break;/*_BUG_ fix #2*/
				ptr=ptr2+1;
				while((*ptr!='\0')&&(!g_ascii_isgraph(*ptr))) ptr++;/*_BUG_*/
				idx++;
			}
			jdx++;
			GUI_COMBOBOX_SET(uspex_gui.numSpecies,jdx);
			g_free(text);
			GUI_COMBOBOX_GET_TEXT(uspex_gui.numSpecies,text);
		}
		g_free(text);
		GUI_COMBOBOX_SET(uspex_gui.numSpecies,index);/*PULL index*/
		/*now get valences from atomType*/
		GUI_COMBOBOX_GET(uspex_gui.atomType,index);/*PUSH index*/
		GUI_LOCK(uspex_gui.atomType);
		/*re-dim*/
		uspex_gui.calc.valences=g_malloc((uspex_gui.calc._nspecies)*sizeof(gint));
		for(idx=0;idx<uspex_gui.calc._nspecies;idx++) {
			GUI_COMBOBOX_SET(uspex_gui.atomType,idx);
			GUI_COMBOBOX_GET_TEXT(uspex_gui.atomType,text);
			sscanf(text,"%*[^(](%*i) (V=%i)",&(val));
			uspex_gui.calc.valences[idx]=val;
			g_free(text);
		}
		GUI_UNLOCK(uspex_gui.atomType);
		GUI_COMBOBOX_SET(uspex_gui.atomType,index);/*PULL index*/
	}else {/*not VARCOMP -> numSpecies is a 1-line block, take it directly from atomType*/
		uspex_gui.calc._var_nspecies=1;
		GUI_COMBOBOX_GET(uspex_gui.atomType,index);/*PUSH index*/
		GUI_LOCK(uspex_gui.atomType);
		/*re-dim*/
		if(uspex_gui.calc.numSpecies!=NULL) g_free(uspex_gui.calc.numSpecies);
		if(uspex_gui.calc.valences!=NULL) g_free(uspex_gui.calc.valences);
		uspex_gui.calc.numSpecies=g_malloc((uspex_gui.calc._nspecies)*sizeof(gint));
		uspex_gui.calc.valences=g_malloc((uspex_gui.calc._nspecies)*sizeof(gint));
		for(idx=0;idx<uspex_gui.calc._nspecies;idx++) {
			GUI_COMBOBOX_SET(uspex_gui.atomType,idx);
			GUI_COMBOBOX_GET_TEXT(uspex_gui.atomType,text);
			sscanf(text,"%*[^(](%i) (V=%i)",&(num),&(val));
			uspex_gui.calc.numSpecies[idx]=num;
			uspex_gui.calc.valences[idx]=val;
			g_free(text);
		}
		GUI_UNLOCK(uspex_gui.atomType);
		GUI_COMBOBOX_SET(uspex_gui.atomType,index);/*PULL index*/
	}
	if(uspex_gui.calc._calctype_mag){
		GUI_REG_VAL(uspex_gui.mag_nm,uspex_gui.calc.magRatio[0],"%lf");
		GUI_REG_VAL(uspex_gui.mag_fmls,uspex_gui.calc.magRatio[1],"%lf");
		GUI_REG_VAL(uspex_gui.mag_fmhs,uspex_gui.calc.magRatio[2],"%lf");
		GUI_REG_VAL(uspex_gui.mag_afml,uspex_gui.calc.magRatio[3],"%lf");
		GUI_REG_VAL(uspex_gui.mag_afmh,uspex_gui.calc.magRatio[4],"%lf");
		GUI_REG_VAL(uspex_gui.mag_fmlh,uspex_gui.calc.magRatio[5],"%lf");
		GUI_REG_VAL(uspex_gui.mag_aflh,uspex_gui.calc.magRatio[6],"%lf");
	}
	GUI_ENTRY_GET_TEXT(uspex_gui.ldaU,text);
	num=0;ptr=text;
	while(*ptr!='\0') {
		if(g_ascii_isdigit(*ptr)) num++;
		ptr++;
	}
	if(num>0){/*we have some value to sync*/
		if(uspex_gui.calc.ldaU!=NULL) g_free(uspex_gui.calc.ldaU);
		uspex_gui.calc.ldaU=g_malloc((uspex_gui.calc._nspecies)*sizeof(gdouble));
		for(idx=0;idx<uspex_gui.calc._nspecies;idx++) uspex_gui.calc.ldaU[idx]=0.;
		ptr=text;
		while((*ptr!='\0')&&(!g_ascii_isgraph(*ptr))) ptr++;
		idx=0;
		while((idx<uspex_gui.calc._nspecies)&&(*ptr!='\0')){
			uspex_gui.calc.ldaU[idx]=g_ascii_strtod(ptr,&ptr2);
			ptr=ptr2+1;
			while((*ptr!='\0')&&(!g_ascii_isgraph(*ptr))) ptr++;
			idx++;
		}
	}
	g_free(text);
	USPEX_REG_VAL(ExternalPressure,"%lf");
if(!uspex_gui.auto_bonds){
	GUI_COMBOBOX_GET(uspex_gui.goodBonds,index);/*PUSH index*/
	GUI_COMBOBOX_SET(uspex_gui.goodBonds,0);
	GUI_COMBOBOX_GET_TEXT(uspex_gui.goodBonds,text);
	if((g_ascii_strcasecmp(text,"ADD GOODBOND") != 0)){
		/*there is something to update*/
		if(uspex_gui.calc.goodBonds!=NULL) g_free(uspex_gui.calc.goodBonds);
		uspex_gui.calc.goodBonds=g_malloc(uspex_gui.calc._nspecies*uspex_gui.calc._nspecies*sizeof(gdouble));
		jdx=0;/*FIX: 5c0648*/
		while((jdx<uspex_gui.calc._nspecies)&&(g_ascii_strcasecmp(text,"ADD GOODBOND") != 0)){
			ptr=text;
			while((*ptr!='\0')&&(!g_ascii_isgraph(*ptr))) ptr++;
			idx=0;
			while((idx<uspex_gui.calc._nspecies)&&(*ptr!='\0')){
				uspex_gui.calc.goodBonds[idx+jdx*uspex_gui.calc._nspecies]=g_ascii_strtod(ptr,&ptr2);
				ptr=ptr2+1;
				while((*ptr!='\0')&&(!g_ascii_isgraph(*ptr))) ptr++;
				idx++;
			}
			jdx++;
			GUI_COMBOBOX_SET(uspex_gui.goodBonds,jdx);
			g_free(text);
			GUI_COMBOBOX_GET_TEXT(uspex_gui.goodBonds,text);
		}
		g_free(text);
	}
	GUI_COMBOBOX_SET(uspex_gui.goodBonds,index);/*PULL index*/
}else{
	if(uspex_gui.calc.goodBonds!=NULL) g_free(uspex_gui.calc.goodBonds);
	uspex_gui.calc.goodBonds=NULL;
}
	//checkMolecules is already sync
	//checkConnectivity is already sync
	USPEX_REG_VAL(fitLimit,"%lf");
	USPEX_REG_VAL(populationSize,"%i");
	USPEX_REG_VAL(initialPopSize,"%i");
	USPEX_REG_VAL(numGenerations,"%i");
	USPEX_REG_VAL(stopCrit,"%i");
	USPEX_REG_VAL(bestFrac,"%lf");
	USPEX_REG_VAL(keepBestHM,"%i");
	//reoptOld is already sync
	USPEX_REG_TEXT(symmetries);
	USPEX_REG_VAL(fracGene,"%lf");
	USPEX_REG_VAL(fracRand,"%lf");
	USPEX_REG_VAL(fracTopRand,"%lf");
	USPEX_REG_VAL(fracPerm,"%lf");
	USPEX_REG_VAL(fracAtomsMut,"%lf");
	USPEX_REG_VAL(fracRotMut,"%lf");
	USPEX_REG_VAL(fracLatMut,"%lf");
	USPEX_REG_VAL(fracSpinMut,"%lf");
	USPEX_REG_VAL(howManySwaps,"%i");
	USPEX_REG_TEXT(specificSwaps);
	USPEX_REG_VAL(mutationDegree,"%lf");
	USPEX_REG_VAL(mutationRate,"%lf");
	USPEX_REG_VAL(DisplaceInLatmutation,"%lf");
	//AutoFrac is already sync
if(uspex_gui.auto_C_ion){
	if(uspex_gui.calc.IonDistances!=NULL) g_free(uspex_gui.calc.IonDistances);
	uspex_gui.calc.IonDistances=NULL;
}
	USPEX_REG_VAL(minVectorLength,"%lf");
	USPEX_REG_VAL(constraint_enhancement,"%i");
if(uspex_gui.calc._nmolecules==0){
	if(uspex_gui.calc.MolCenters!=NULL) g_free(uspex_gui.calc.MolCenters);
	uspex_gui.calc.MolCenters=NULL;
}
if(uspex_gui.auto_C_lat){
	if(uspex_gui.calc.Latticevalues!=NULL) g_free(uspex_gui.calc.Latticevalues);
	uspex_gui.calc.Latticevalues=NULL;
	uspex_gui.calc._nlattice_line=0;
	uspex_gui.calc._nlattice_vals=0;
}
	GUI_ENTRY_GET_TEXT(uspex_gui.splitInto,text);
	ptr=text;
	while((*ptr!='\0')&&(!g_ascii_isgraph(*ptr))) ptr++;
	if(uspex_gui.calc.splitInto!=NULL) g_free(uspex_gui.calc.splitInto);
	uspex_gui.calc.splitInto=NULL;
	uspex_gui.calc._nsplits=0;
	if(*ptr!='\0'){
		/*count number of splits*/
		while(*ptr!='\0'){
			if(g_ascii_strtoull(ptr,&ptr2,10)==G_MAXUINT64) break;/*FIX: 38832a*/
			if(ptr2==ptr) break;/*not a valid integer*/
			uspex_gui.calc._nsplits++;
			if(ptr2==NULL) break;/*_BUG_ fix #1*/
			if(*ptr2=='\0') break;/*_BUG_ fix #2*/
			ptr=ptr2+1;
			while((*ptr!='\0')&&(!g_ascii_isgraph(*ptr))) ptr++;/*_BUG_*/
		}
		/*get them*/
		uspex_gui.calc.splitInto=g_malloc(uspex_gui.calc._nsplits*sizeof(gint));
		for(idx=0;idx<uspex_gui.calc._nsplits;idx++) uspex_gui.calc.splitInto[idx]=0;
		ptr=text;idx=0;
		while((idx<uspex_gui.calc._nsplits)&&(*ptr!='\0')){
			uspex_gui.calc.splitInto[idx]=(gint)g_ascii_strtoull(ptr,&ptr2,10);
			if(ptr2==NULL) break;/*_BUG_ fix #1*/
			if(*ptr2=='\0') break;/*_BUG_ fix #2*/
			ptr=ptr2+1;
			while((*ptr!='\0')&&(!g_ascii_isgraph(*ptr))) ptr++;/*_BUG_*/
			idx++;
		}
	}
	//pickUpYN is already sync
	USPEX_REG_VAL(pickUpGen,"%i");
	USPEX_REG_VAL(pickUpFolder,"%i");
	uspex_gui.calc._num_opt_steps=(gint)uspex_gui._tmp_num_opt_steps;
	//abinitioCode is already sync
	//KresolStart is already sync
	//vacuumSize is already sync
	//_isfixed is already sync
        GUI_ENTRY_GET_TEXT(uspex_gui.numProcessors,text);
        ptr=text;
        while((*ptr!='\0')&&(!g_ascii_isgraph(*ptr))) ptr++;
        if(uspex_gui.calc.numProcessors!=NULL) g_free(uspex_gui.calc.numProcessors);
        uspex_gui.calc.numProcessors=NULL;
        if(*ptr!='\0') {
                /*There is _num_opt_steps numProcessors*/
                uspex_gui.calc.numProcessors=g_malloc(uspex_gui.calc._num_opt_steps*sizeof(gint));
                for(idx=0;idx<uspex_gui.calc._num_opt_steps;idx++) uspex_gui.calc.numProcessors[idx]=0;
                idx=0;
                while((idx<uspex_gui.calc._num_opt_steps)&&(*ptr!='\0')){
                        uspex_gui.calc.numProcessors[idx]=(gint)g_ascii_strtoull(ptr,&ptr2,10);
                        ptr=ptr2+1;
                        while((*ptr!='\0')&&(!g_ascii_isgraph(*ptr))) ptr++;
                        idx++;
                }
        }
	USPEX_REG_VAL(numParallelCalcs,"%i");
	register_commandExecutable();
	USPEX_REG_VAL(whichCluster,"%i");
	USPEX_REG_TEXT(remoteFolder);
	//PhaseDiagram is already sync
	USPEX_REG_VAL(RmaxFing,"%lf");
	USPEX_REG_VAL(deltaFing,"%lf");
	USPEX_REG_VAL(sigmaFing,"%lf");
	USPEX_REG_VAL(toleranceFing,"%lf");
	USPEX_REG_VAL(antiSeedsActivation,"%i");
	USPEX_REG_VAL(antiSeedsMax,"%lf");
	USPEX_REG_VAL(antiSeedsSigma,"%lf");
	//doSpaceGroup is already sync
	USPEX_REG_VAL(SymTolerance,"%lf");
	USPEX_REG_VAL(repeatForStatistics,"%i");
	USPEX_REG_VAL(stopFitness,"%lf");
	USPEX_REG_VAL(fixRndSeed,"%i");
	//collectForces is already sync
	//ordering_active is already sync
	//symmetrize is already sync
	GUI_ENTRY_GET_TEXT(uspex_gui.valenceElectr,text);
	ptr=text;
	while((*ptr!='\0')&&(!g_ascii_isgraph(*ptr))) ptr++;
	if(uspex_gui.calc.valenceElectr!=NULL) g_free(uspex_gui.calc.valenceElectr);
	uspex_gui.calc.valenceElectr=NULL;
	if(*ptr!='\0') {
		/*There is _nspecies valenceElectr*/
		uspex_gui.calc.valenceElectr=g_malloc(uspex_gui.calc._nspecies*sizeof(gint));
		for(idx=0;idx<uspex_gui.calc._nspecies;idx++) uspex_gui.calc.valenceElectr[idx]=0;
		idx=0;
		while((idx<uspex_gui.calc._nspecies)&&(*ptr!='\0')){
			uspex_gui.calc.valenceElectr[idx]=(gint)g_ascii_strtoull(ptr,&ptr2,10);
			ptr=ptr2+1;
			while((*ptr!='\0')&&(!g_ascii_isgraph(*ptr))) ptr++;
			idx++;
		}
	}
	USPEX_REG_VAL(percSliceShift,"%lf");
	GUI_COMBOBOX_GET(uspex_gui.dynamicalBestHM,index);
	uspex_gui.calc.dynamicalBestHM=index;
	USPEX_REG_TEXT(softMutOnly);
	USPEX_REG_VAL(maxDistHeredity,"%lf");
	GUI_COMBOBOX_GET(uspex_gui.manyParents,index);
	uspex_gui.calc.manyParents=index;
	USPEX_REG_VAL(minSlice,"%lf");
	USPEX_REG_VAL(maxSlice,"%lf");
	USPEX_REG_VAL(numberparents,"%i");
	//_nmolecules is already sync
	USPEX_REG_VAL(BoltzTraP_T_max,"%lf");
	USPEX_REG_VAL(BoltzTraP_T_delta,"%lf");
	USPEX_REG_VAL(BoltzTraP_T_efcut,"%lf");
	USPEX_REG_VAL(TE_T_interest,"%lf");
	USPEX_REG_VAL(TE_threshold,"%lf");
	//TE_goal is already sync
	USPEX_REG_VAL(thicknessS,"%lf");
	USPEX_REG_VAL(thicknessB,"%lf");
	USPEX_REG_VAL(reconstruct,"%i");
	GUI_ENTRY_GET_TEXT(uspex_gui.StoichiometryStart,text);
	ptr=text;
	while((*ptr!='\0')&&(!g_ascii_isgraph(*ptr))) ptr++;
	if(uspex_gui.calc.StoichiometryStart!=NULL) g_free(uspex_gui.calc.StoichiometryStart);
	uspex_gui.calc.StoichiometryStart=NULL;
	if(*ptr!='\0') {
		/*There is _nspecies StoichiometryStart*/
		uspex_gui.calc.StoichiometryStart=g_malloc(uspex_gui.calc._nspecies*sizeof(gint));
		for(idx=0;idx<uspex_gui.calc._nspecies;idx++) uspex_gui.calc.StoichiometryStart[idx]=0;
		idx=0;
		while((idx<uspex_gui.calc._nspecies)&&(*ptr!='\0')){
			uspex_gui.calc.StoichiometryStart[idx]=(gint)g_ascii_strtoull(ptr,&ptr2,10);
			ptr=ptr2+1;
			while((*ptr!='\0')&&(!g_ascii_isgraph(*ptr))) ptr++;
			idx++;
		}
	}
	USPEX_REG_VAL(E_AB,"%lf");
	USPEX_REG_VAL(Mu_A,"%lf");
	USPEX_REG_VAL(Mu_B,"%lf");
	USPEX_REG_VAL(firstGeneMax,"%i");
	USPEX_REG_VAL(minAt,"%i");
	USPEX_REG_VAL(maxAt,"%i");
	USPEX_REG_VAL(fracTrans,"%lf");
	USPEX_REG_VAL(howManyTrans,"%lf");
	GUI_ENTRY_GET_TEXT(uspex_gui.specificTrans,text);
	ptr=text;
	while((*ptr!='\0')&&(!g_ascii_isgraph(*ptr))) ptr++;
	if(uspex_gui.calc.specificTrans!=NULL) g_free(uspex_gui.calc.specificTrans);
	uspex_gui.calc.specificTrans=NULL;
	uspex_gui.calc._nspetrans=0;
	if(*ptr!='\0'){
		/*count number of specificTrans*/
		while(*ptr!='\0'){
			if(g_ascii_strtoull(ptr,&ptr2,10)==G_MAXUINT64) break;/*FIX: 73cf4c*/
			if(ptr2==ptr) break;/*not a valid integer*/
			uspex_gui.calc._nspetrans++;
			ptr=ptr2+1;
			while((*ptr!='\0')&&(!g_ascii_isgraph(*ptr))) ptr++;
		}
		/*get them*/
		uspex_gui.calc.specificTrans=g_malloc(uspex_gui.calc._nspetrans*sizeof(gint));
		for(idx=0;idx<uspex_gui.calc._nspetrans;idx++) uspex_gui.calc.specificTrans[idx]=0;
		ptr=text;idx=0;
		while((idx<uspex_gui.calc._nspetrans)&&(*ptr!='\0')){
			uspex_gui.calc.specificTrans[idx]=(gint)g_ascii_strtoull(ptr,&ptr2,10);
			ptr=ptr2+1;
			while((*ptr!='\0')&&(!g_ascii_isgraph(*ptr))) ptr++;
			idx++;
		}
	}
	//ExternalPressure is already sync
	USPEX_REG_VAL(GaussianWidth,"%lf");
	USPEX_REG_VAL(GaussianHeight,"%lf");
	//FullRelax is already sync
	USPEX_REG_VAL(maxVectorLength,"%lf");
	USPEX_REG_VAL(PSO_softMut,"%lf");
	USPEX_REG_VAL(PSO_BestStruc,"%lf");
	USPEX_REG_VAL(PSO_BestEver,"%lf");
	USPEX_REG_VAL(vcnebType,"%i");
	//_vcnebtype_method is already sync
	//_vcnebtype_img_num is already sync
	//_vcnebtype_spring is already sync
	USPEX_REG_VAL(numImages,"%i");
	USPEX_REG_VAL(numSteps,"%i");
	//optReadImages is already sync
	//optimizerType is already sync
	//optRelaxType is already sync
	USPEX_REG_VAL(dt,"%lf");
	USPEX_REG_VAL(ConvThreshold,"%lf");
	USPEX_REG_VAL(VarPathLength,"%lf");
	USPEX_REG_VAL(K_min,"%lf");
	USPEX_REG_VAL(K_max,"%lf");
	USPEX_REG_VAL(Kconstant,"%lf");
	//optFreezing is already sync
	//optMethodCIDI is already sync
	USPEX_REG_VAL(startCIDIStep,"%i");
	GUI_ENTRY_GET_TEXT(uspex_gui.pickupImages,text);
	ptr=text;
	while((*ptr!='\0')&&(!g_ascii_isgraph(*ptr))) ptr++;
	if(uspex_gui.calc.pickupImages!=NULL) g_free(uspex_gui.calc.pickupImages);
	uspex_gui.calc.pickupImages=NULL;
	uspex_gui.calc._npickimg=0;
	if(*ptr!='\0'){
		/*count number of pickupImages*/
		while(*ptr!='\0'){
			if(g_ascii_strtoull(ptr,&ptr2,10)==G_MAXUINT64) break;/*FIX: 61c045*/
			if(ptr2==ptr) break;/*not a valid integer*/
			uspex_gui.calc._npickimg++;
			ptr=ptr2+1;
			while((*ptr!='\0')&&(!g_ascii_isgraph(*ptr))) ptr++;
		}
		/*get them*/
		uspex_gui.calc.pickupImages=g_malloc(uspex_gui.calc._npickimg*sizeof(gint));
		for(idx=0;idx<uspex_gui.calc._npickimg;idx++) uspex_gui.calc.pickupImages[idx]=0;
		ptr=text;idx=0;
		while((idx<uspex_gui.calc._npickimg)&&(*ptr!='\0')){
			uspex_gui.calc.pickupImages[idx]=(gint)g_ascii_strtoull(ptr,&ptr2,10);
			ptr=ptr2+1;
			while((*ptr!='\0')&&(!g_ascii_isgraph(*ptr))) ptr++;
			idx++;
		}
	}
	//FormatType is already sync
	USPEX_REG_VAL(PrintStep,"%i");
	USPEX_REG_VAL(numIterations,"%i");
	USPEX_REG_TEXT(speciesSymbol);
	GUI_ENTRY_GET_TEXT(uspex_gui.mass,text);
	ptr=text;
	while((*ptr!='\0')&&(!g_ascii_isgraph(*ptr))) ptr++;
	if(uspex_gui.calc.mass!=NULL) g_free(uspex_gui.calc.mass);
	uspex_gui.calc.mass=NULL;
	if(*ptr!='\0') {
		/*there is _nspecies mass*/
		uspex_gui.calc.mass=g_malloc(uspex_gui.calc._nspecies*sizeof(gint));
		for(idx=0;idx<uspex_gui.calc._nspecies;idx++) uspex_gui.calc.mass[idx]=0.;
		idx=0;
		while((idx<uspex_gui.calc._nspecies)&&(*ptr!='\0')){
			uspex_gui.calc.mass[idx]=g_ascii_strtod(ptr,&ptr2);
			ptr=ptr2+1;
			while((*ptr!='\0')&&(!g_ascii_isgraph(*ptr))) ptr++;
			idx++;
		}
	}
	GUI_REG_VAL(uspex_gui.amplitudeShoot_AB,uspex_gui.calc.amplitudeShoot[0],"%lf");
	GUI_REG_VAL(uspex_gui.amplitudeShoot_BA,uspex_gui.calc.amplitudeShoot[1],"%lf");
	GUI_REG_VAL(uspex_gui.magnitudeShoot_success,uspex_gui.calc.magnitudeShoot[0],"%lf");
	GUI_REG_VAL(uspex_gui.magnitudeShoot_failure,uspex_gui.calc.magnitudeShoot[1],"%lf");
	USPEX_REG_VAL(shiftRatio,"%lf");
	//orderParaType is already sync
	GUI_REG_VAL(uspex_gui.opCriteria_start,uspex_gui.calc.opCriteria[0],"%lf");
	GUI_REG_VAL(uspex_gui.opCriteria_end,uspex_gui.calc.opCriteria[1],"%lf");
	USPEX_REG_TEXT(cmdOrderParameter);
	USPEX_REG_TEXT(cmdEnthalpyTemperature);
	USPEX_REG_TEXT(orderParameterFile);
	USPEX_REG_TEXT(enthalpyTemperatureFile);
	USPEX_REG_TEXT(trajectoryFile);
	USPEX_REG_TEXT(MDrestartFile);
	USPEX_REG_TEXT(job_path);
	USPEX_REG_TEXT(job_uspex_exe);
}
/***************************/
/* save current parameters */
/***************************/
gint save_uspex_calc(){
#define DEBUG_USPEX_SAVE 0
	gchar *filename;
	/*1-synchronize*/
	uspex_gui_sync();
	/*2-output*/
	//SAVE calculation parameters to INPUT.txt
	filename=g_strdup_printf("%s/INPUT.txt",uspex_gui.calc.job_path);
	if(dump_uspex_parameters(filename,&(uspex_gui.calc))<0){
		fprintf(stderr,"#ERR: while saving INPUT.txt!\n");
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
/***************************************/
/* copy Specific directory (if exists) */
/***************************************/
void uspex_copy_specific(uspex_exec_struct *uspex_exec){
GDir *src;
gchar *source;
gchar *target;
/**/
if(uspex_exec==NULL) return;
if(uspex_gui.calc.path==NULL) return;
source=g_build_filename(uspex_gui.calc.path,"../Specific/", NULL);
src=g_dir_open(source,0,NULL);
if(src==NULL) {
	fprintf(stderr,"Can't open Specific directory.\n");
	return;
}
target = g_build_filename((*uspex_exec).job_path,"Specific", NULL);
if(g_mkdir(target,0775)){
	fprintf(stderr,"Can't create %s directory.\n",target);
	g_dir_close(src);
	return;
}
g_dir_close(src);
dumb_dir_copy(source,target);
g_free(source);
g_free(target);
}
/******************************************/
/* Specific auto-settings for VASP: INCAR */
/******************************************/
void vasp_specific_step_incar(FILE *there,gint step_n){
/*PSO and USPEX method are supposed to take identical steps, COPEX method is undefined*/
gint tot=uspex_gui.calc._num_opt_steps;
gdouble sigma;
gboolean is_meta;
gchar *title;
is_meta=((uspex_gui.calc.calculationMethod==US_CM_META)||(uspex_gui.calc.calculationMethod==US_CM_MINHOP));
if((uspex_gui.calc.calculationMethod == US_CM_VCNEB)||(uspex_gui.calc.calculationMethod == US_CM_TPS)){
/*in that case, there should be only one step*/
	if(step_n!=1) {
		title = g_strdup_printf("USPEX: VC-NEB & TPS calculation only require 1 calculation step!\nUSPEX: STEP %i IGNORED!\n",step_n);
		gui_text_show(ERROR,title);
		g_free(title);
		return;
	}
	fprintf(there,"SYSTEM = GDIS_GEN_%i\n",step_n);
	fprintf(there,"! USPEX: VC-NEB / TPS STEP, AUTO-GENERATED!\n");
	fprintf(there,"ISTART = 0\n");
	fprintf(there,"ISMEAR = 1\n");
	fprintf(there,"SIGMA = 0.05\n");
	fprintf(there,"PREC = HIGH\n");/*see note [0]*/
	fprintf(there,"EDIFF = 1E-5\n");
	fprintf(there,"EDIFFG = 1E-4\n");/*useless*/
	fprintf(there,"NSW = 0\n");
	/*to get forces details*/
	fprintf(there,"ISIF = 2\n");
	fprintf(there,"IBRION = 2\n");
	fprintf(there,"POTIM = 0.02\n");
	return;
}
if(tot>3){
	/*always*/
	fprintf(there,"SYSTEM = GDIS_GEN_%i\n",step_n);
	if(is_meta) fprintf(there,"! USPEX: META / MINHOP STEP, AUTO-GENERATED!\n");
	else fprintf(there,"! USPEX: USPEX / PSO STEP, AUTO-GENERATED!\n");
	fprintf(there,"ISTART = 0\n");
	/*This is the 'recommended' smearing by USPEX for an insulator*/
	/*however, it should perform as well (maybe better) for metals*/
	/*(Maybe a gaussian smearing should be more universal) --OVHPA*/
	fprintf(there,"ISMEAR = 1\n");
	sigma = 0.1 - ((double) step_n - 1.0) * 0.05 / ((double) tot - 2.0);
	if(step_n==1){
		fprintf(there,"SIGMA = 0.1\n");
		fprintf(there,"PREC = LOW\n");
		fprintf(there,"EDIFF = 1E-2\n");
		fprintf(there,"EDIFFG = 1E-2\n");
		fprintf(there,"NSW = 40\n");/*see note [1]*/
		if(is_meta) {
			if(uspex_gui.calc._isfixed[step_n-1]) fprintf(there,"ISIF = 2\n");/*FIXED CELL*/
			else fprintf(there,"ISIF = 4\n");/*cell shape only*/
		}else fprintf(there,"ISIF = 4\n");/*cell shape only*/
		fprintf(there,"IBRION = 2\n");
		fprintf(there,"POTIM = 0.02\n");
		return;
	}
	if(step_n==tot){
		fprintf(there,"SIGMA = 0.05\n");
		fprintf(there,"PREC = HIGH\n");
		fprintf(there,"EDIFF = 1E-5\n");
		fprintf(there,"NSW = 0\n");
		return;
	}
	fprintf(there,"SIGMA = %f\n",sigma);
	switch(tot){
	case 6:/*6 steps calculation*/
		switch (step_n){
		case 2:
			fprintf(there,"PREC = NORMAL\n");
			fprintf(there,"EDIFF = 1E-3\n");
			fprintf(there,"EDIFFG = 1E-2\n");
			fprintf(there,"NSW = 80\n");/*see note [1]*/
			if(is_meta){
				if(uspex_gui.calc._isfixed[step_n-1]) fprintf(there,"ISIF = 2\n");/*FIXED CELL*/
				else fprintf(there,"ISIF = 4\n");/*second cell shape only*/
			} else fprintf(there,"ISIF = 4\n");/*second cell shape only*/
			fprintf(there,"IBRION = 1\n");
			fprintf(there,"POTIM = 0.3\n");
			return;
		case 3:
			fprintf(there,"PREC = NORMAL\n");
			fprintf(there,"EDIFF = 1E-3\n");
			fprintf(there,"EDIFFG = 1E-2\n");
			fprintf(there,"NSW = 50\n");/*see note [1]*/
			if(is_meta){
				if(uspex_gui.calc._isfixed[step_n-1]) fprintf(there,"ISIF = 2\n");/*FIXED CELL*/
				else fprintf(there,"ISIF = 3\n");/*optimize everything*/
			}else fprintf(there,"ISIF = 3\n");/*optimize everything*/
			fprintf(there,"IBRION = 2\n");
			fprintf(there,"POTIM = 0.02\n");
			return;
		case 4:
			fprintf(there,"PREC = HIGH\n");/*see note [2]*/
			fprintf(there,"EDIFF = 1E-4\n");
			fprintf(there,"EDIFFG = 1E-2\n");
			fprintf(there,"NSW = 40\n");/*see note [1]*/
			if(is_meta){
				if(uspex_gui.calc._isfixed[step_n-1]) fprintf(there,"ISIF = 2\n");/*FIXED CELL*/
				else fprintf(there,"ISIF = 3\n");/*optimize everything*/
			}else fprintf(there,"ISIF = 3\n");/*optimize everything*/
			fprintf(there,"IBRION = 2\n");
			fprintf(there,"POTIM = 0.02\n");
			return;
		case 5:
			fprintf(there,"PREC = HIGH\n");/*see note [2]*/
			fprintf(there,"EDIFF = 1E-5\n");
			fprintf(there,"EDIFFG = 1E-3\n");
			fprintf(there,"NSW = 120\n");/*see note [1]*/
			if(is_meta){
				/*metadynamics is slightly different here*/
				if(uspex_gui.calc._isfixed[step_n-1]) fprintf(there,"ISIF = 2\n");/*FIXED CELL*/
				else fprintf(there,"ISIF = 3\n");/*optimize everything*/
			}else fprintf(there,"ISIF = 2\n");/*optimize ion position*/
			fprintf(there,"IBRION = 2\n");
			fprintf(there,"POTIM = 0.02\n");
			return;
		case 1:/*already dealt with*/
		case 6:/*idem*/
		default:
			/*nope*/
			return;
		}
		break;
	case 5:/*5 steps calculation*/
		switch (step_n){
		case 2:
			fprintf(there,"PREC = NORMAL\n");
			fprintf(there,"EDIFF = 1E-3\n");
			fprintf(there,"EDIFFG = 1E-2\n");
			fprintf(there,"NSW = 85\n");/*see note [1]*/
			if(is_meta){
				if(uspex_gui.calc._isfixed[step_n-1]) fprintf(there,"ISIF = 2\n");/*FIXED CELL*/
				else fprintf(there,"ISIF = 3\n");/*optimize everything*/
			}else fprintf(there,"ISIF = 3\n");/*optimize everything*/
			fprintf(there,"IBRION = 2\n");
			fprintf(there,"POTIM = 0.02\n");
			return;
		case 3:
			fprintf(there,"PREC = HIGH\n");/*see note [2]*/
			fprintf(there,"EDIFF = 1E-4\n");
			fprintf(there,"EDIFFG = 1E-2\n");
			fprintf(there,"NSW = 55\n");/*see note [1]*/
			if(is_meta){
				if(uspex_gui.calc._isfixed[step_n-1]) fprintf(there,"ISIF = 2\n");/*FIXED CELL*/
				else fprintf(there,"ISIF = 3\n");/*optimize everything*/
			}else fprintf(there,"ISIF = 3\n");/*optimize everything*/
			fprintf(there,"IBRION = 2\n");
			fprintf(there,"POTIM = 0.02\n");
			return;
		case 4:
			fprintf(there,"PREC = HIGH\n");/*see note [2]*/
			fprintf(there,"EDIFF = 1E-5\n");
			fprintf(there,"EDIFFG = 1E-3\n");
			fprintf(there,"NSW = 120\n");/*see note [1]*/
			if(is_meta){
				/*metadynamics is slightly different here*/
				if(uspex_gui.calc._isfixed[step_n-1]) fprintf(there,"ISIF = 2\n");/*FIXED CELL*/
				else fprintf(there,"ISIF = 3\n");/*optimize everything*/
			}else fprintf(there,"ISIF = 2\n");/*optimize ion position*/
			fprintf(there,"IBRION = 2\n");
			fprintf(there,"POTIM = 0.02\n");
			return;
		case 1:/*already dealt with*/
		case 5:/*idem*/
		default:
			/*nope*/
			return;
		}
		break;
	case 4:
		switch (step_n){
		case 2:
			fprintf(there,"PREC = NORMAL\n");
			fprintf(there,"EDIFF = 1E-3\n");
			fprintf(there,"EDIFFG = 1E-2\n");
			fprintf(there,"NSW = 90\n");/*see note [1]*/
			if(is_meta){
				if(uspex_gui.calc._isfixed[step_n-1]) fprintf(there,"ISIF = 2\n");/*FIXED CELL*/
				else fprintf(there,"ISIF = 3\n");/*optimize everything*/
			}else fprintf(there,"ISIF = 3\n");/*optimize everything*/
			fprintf(there,"IBRION = 2\n");
			fprintf(there,"POTIM = 0.02\n");
			return;
		case 3:
			fprintf(there,"PREC = HIGH\n");/*see note [2]*/
			fprintf(there,"EDIFF = 1E-4\n");
			fprintf(there,"EDIFFG = 1E-3\n");
			fprintf(there,"NSW = 140\n");/*see note [1]*/
			if(is_meta){
				if(uspex_gui.calc._isfixed[step_n-1]) fprintf(there,"ISIF = 2\n");/*FIXED CELL*/
				else fprintf(there,"ISIF = 3\n");/*optimize everything*/
			}else fprintf(there,"ISIF = 3\n");/*optimize everything*/
			fprintf(there,"IBRION = 2\n");
			fprintf(there,"POTIM = 0.02\n");
			return;
		case 1:/*already dealt with*/
		case 4:/*idem*/
		default:
			/*nope*/
			return;
		}
		break;
	default:
		/*in case tot > 6*/
		title=g_strdup_printf("USPEX: Recommended number of steps for a USPEX JOB is 4~6\n");
		gui_text_show(ERROR,title);
		g_free(title);		
		switch (step_n){
		case 2:
			fprintf(there,"PREC = NORMAL\n");
			fprintf(there,"EDIFF = 1E-3\n");
			fprintf(there,"EDIFFG = 1E-2\n");
			fprintf(there,"NSW = 85\n");/*see note [1]*/
			if(is_meta){
				if(uspex_gui.calc._isfixed[step_n-1]) fprintf(there,"ISIF = 2\n");/*FIXED CELL*/
				else fprintf(there,"ISIF = 3\n");/*optimize everything*/
			}else fprintf(there,"ISIF = 3\n");/*optimize everything*/
			fprintf(there,"IBRION = 2\n");
			fprintf(there,"POTIM = 0.02\n");
			return;
		case 3:
			fprintf(there,"PREC = HIGH\n");/*see note [2]*/
			fprintf(there,"EDIFF = 1E-4\n");
			fprintf(there,"EDIFFG = 1E-2\n");
			fprintf(there,"NSW = 55\n");/*see note [1]*/
			if(is_meta){
				if(uspex_gui.calc._isfixed[step_n-1]) fprintf(there,"ISIF = 2\n");/*FIXED CELL*/
				else fprintf(there,"ISIF = 3\n");/*optimize everything*/
			}else fprintf(there,"ISIF = 3\n");/*optimize everything*/
			fprintf(there,"IBRION = 2\n");
			fprintf(there,"POTIM = 0.02\n");
			return;
		case 4:
			fprintf(there,"PREC = HIGH\n");/*see note [2]*/
			fprintf(there,"EDIFF = 1E-5\n");
			fprintf(there,"EDIFFG = 1E-3\n");
			fprintf(there,"NSW = 120\n");/*see note [1]*/
			if(is_meta){
				/*metadynamics is slightly different here*/
				if(uspex_gui.calc._isfixed[step_n-1]) fprintf(there,"ISIF = 2\n");/*FIXED CELL*/
				else fprintf(there,"ISIF = 3\n");/*optimize everything*/
			}else fprintf(there,"ISIF = 2\n");/*optimize ion position*/
			fprintf(there,"IBRION = 2\n");
			fprintf(there,"POTIM = 0.02\n");
			return;
		case 1:/*already dealt with*/
		default:
			/*nope*/
			break;
		}
		if(step_n<tot){
			/*deal with the extra steps*/
			fprintf(there,"! USPEX: EXTRA STEPS MIGHT NOT BE USEFUL!!\n");
			fprintf(there,"! USPEX: 4~6 VASP STEPS ARE RECOMMENDED...\n");
			fprintf(there,"PREC = HIGH\n");/*see note [2]*/
			fprintf(there,"EDIFF = 1E-5\n");
			if(step_n>6) fprintf(there,"EDIFFG = 1E-3\n");
			else fprintf(there,"EDIFFG = 1E-4\n");
			fprintf(there,"NSW = 50\n");/*JUST CYCLE 50 steps cycles of ISIF 2/3*/
			if(is_meta){
				/*metadynamics is slightly different here*/
				if(uspex_gui.calc._isfixed[step_n-1]) fprintf(there,"ISIF = 2\n");/*FIXED CELL*/
				else fprintf(there,"ISIF = 3\n");/*optimize everything*/
			}else {
				if(step_n%2) fprintf(there,"ISIF = 2\n");
				else fprintf(there,"ISIF = 3\n");
			}
			fprintf(there,"IBRION = 2\n");
			fprintf(there,"POTIM = 0.02\n");
			return;
		}
		break;
	}
}else{
	title=g_strdup_printf("USPEX: Recommended number of steps for a USPEX JOB is 4~6\n");
	gui_text_show(ERROR,title);
	g_free(title);
	switch (tot){
	case 3:
		fprintf(there,"SYSTEM = GDIS_GEN_%i\n",step_n);
		if(is_meta) fprintf(there,"! USPEX: META / MINHOP STEP, AUTO-GENERATED!\n");
		else fprintf(there,"! USPEX: USPEX / PSO STEP, AUTO-GENERATED!\n");
		/*This is an unsafe setting, make it clear*/
		fprintf(there,"! UNSAFE SETTING! PLEASE USE AT LEAST 4~6 VASP STEPS!\n");
		fprintf(there,"ISTART = 0\n");
		fprintf(there,"ISMEAR = 1\n");
		switch (step_n){
		case 1:
			/*see note [3]*/
			fprintf(there,"SIGMA = 0.1\n");
			fprintf(there,"PREC = NORMAL\n");
			fprintf(there,"NELMIN = 5\n");
			fprintf(there,"EDIFF = 1E-2\n");
			fprintf(there,"EDIFFG = -0.01\n");
			fprintf(there,"NSW = 40\n");
			if(is_meta){
				if(uspex_gui.calc._isfixed[step_n-1]) fprintf(there,"ISIF = 2\n");/*FIXED CELL*/
				else fprintf(there,"ISIF = 3\n");/*optimize everything*/
			}else fprintf(there,"ISIF = 3\n");
			fprintf(there,"IBRION = 3\n");
			fprintf(there,"SMASS = 0.3\n");
			fprintf(there,"POTIM = 0.03\n");
			return;
		case 2:
			fprintf(there,"SIGMA = 0.075\n");
			fprintf(there,"PREC = HIGH\n");/*see note [2]*/
			fprintf(there,"EDIFF = 1E-4\n");
			fprintf(there,"EDIFFG = 1E-3\n");
			fprintf(there,"NSW = 160\n");/*see note [1]*/
			if(is_meta){
				if(uspex_gui.calc._isfixed[step_n-1]) fprintf(there,"ISIF = 2\n");/*FIXED CELL*/
				else fprintf(there,"ISIF = 3\n");/*optimize everything*/
			}else fprintf(there,"ISIF = 3\n");/*optimize everything*/
			fprintf(there,"IBRION = 2\n");
			fprintf(there,"POTIM = 0.02\n");			
			return;
		case 3:
		default:
			fprintf(there,"SIGMA = 0.05\n");
			fprintf(there,"PREC = HIGH\n");
			fprintf(there,"EDIFF = 1E-5\n");
			fprintf(there,"NSW = 0\n");
			return;
		}
	case 2:/*2 steps should be for test/benchmark only*/
		fprintf(there,"SYSTEM = GDIS_GEN_%i\n",step_n);
		if(is_meta) fprintf(there,"! USPEX: META / MINHOP STEP, AUTO-GENERATED!\n");
		else fprintf(there,"! USPEX: USPEX / PSO STEP, AUTO-GENERATED!\n");
		/*This is an unsafe setting, make it clear*/
		fprintf(there,"! UNSAFE SETTING! PLEASE USE AT LEAST 4~6 VASP STEPS!\n");
		fprintf(there,"ISTART = 0\n");
		fprintf(there,"ISMEAR = 1\n");
		switch (step_n){
		case 1:
			/*see note [3]*/
			fprintf(there,"SIGMA = 0.1\n");
			fprintf(there,"PREC = NORMAL\n");
			fprintf(there,"NELMIN = 5\n");
			fprintf(there,"EDIFF = 1E-2\n");
			fprintf(there,"EDIFFG = -0.01\n");/*should almost never be reached in 40 steps*/
			fprintf(there,"NSW = 40\n");
			if(is_meta){
				if(uspex_gui.calc._isfixed[step_n-1]) fprintf(there,"ISIF = 2\n");/*FIXED CELL*/
				else fprintf(there,"ISIF = 3\n");/*optimize everything*/
			}else fprintf(there,"ISIF = 3\n");
			fprintf(there,"IBRION = 3\n");
			fprintf(there,"SMASS = 0.3\n");
			fprintf(there,"POTIM = 0.03\n");
			return;
		case 2:
		default:
			fprintf(there,"SIGMA = 0.05\n");
			fprintf(there,"PREC = HIGH\n");/*see note [2]*/
			fprintf(there,"EDIFF = 1E-5\n");
			fprintf(there,"EDIFFG = -0.01\n");
			fprintf(there,"NSW = 260\n");/*see note [1]*/
			if(is_meta){
				if(uspex_gui.calc._isfixed[step_n-1]) fprintf(there,"ISIF = 2\n");/*FIXED CELL*/
				else fprintf(there,"ISIF = 3\n");/*optimize everything*/
			}else fprintf(there,"ISIF = 3\n");/*optimize everything*/
			fprintf(there,"IBRION = 2\n");
			fprintf(there,"POTIM = 0.02\n");			
			return;
		}
	case 1:
	default:
		/*only one VASP step makes no sense!*/
		fprintf(there,"! USPEX: ONE VASP STEP MAKES NO SENSE!!!!\n");
		fprintf(there,"! USPEX: the possibility is left for both\n");
		fprintf(there,"! USPEX: TESTING and BENCHMARK, all other\n");
		fprintf(there,"! USPEX: cases should be considered WRONG\n");
		fprintf(there,"! USPEX: and with NO significance.--OVHPA\n");
		fprintf(there,"SIGMA = 0.75\n");
		fprintf(there,"PREC = HIGH\n");/*see note [2]*/
		fprintf(there,"NELMIN = 6\n");
		fprintf(there,"EDIFF = 1E-05\n");
		fprintf(there,"EDIFFG = 1E-05\n");
		fprintf(there,"NSW = 200\n");
		fprintf(there,"ISIF = 2\n");
		fprintf(there,"IBRION = 2\n");/*we do not relax cell/cell shape!*/
		fprintf(there,"POTIM = 0.15\n");
		/*we need an extra warning*/
		title=g_strdup_printf("USPEX: DOING a ONE-VASP-STEP USPEX JOB IS A TERRIBLE IDEA!\nUSPEX: results will be inconsistent and with NO significance!\n");
		gui_text_show(ERROR,title);
		g_free(title);
		return;
	}
}
/* NOTES:
[0] PREC = ACCURATE  is necessary to obtain a good estimate of the forces
and stress, however PREC = HIGH is a setting that will increase the ENCUT
by 30% (25% ENMAX) and that prove to be a better choice for systems under
a non-negligible external pressure.
[1] The number of necessary ionic steps have been determined by extensive
testing on _simple_ cases. They should provide a "good-enough" quality vs
speed convergence. However, it doesn't mean that they will be good enough
for all systems and setting: one need to check that, for most steps, good
convergence is achieved before maximum steps is reached.
  This is especially necessary for the penultimate step!
[2] The switch to HIGH precisions grant a 25% (30% in VASP documentation)
increase of ENMAX (ENCUT), amoung other acurqacy optimizations. While the
official documentation recommend to set ENCUT manually at all time, using
PREC=HIGH allow accurate forces description for all cases, including ones
under external pressure. The last step with PREC=Accurate ensure that the
energy is calculated properly. It can be replace with PREC=HIGH for cases
where forces are also required.
[3] When there is obviously not enough optimization steps (recommended is
4~6 steps), we use a damped algorithm with a damping factor, SMASS, and a
time step, POTIM, of 0.3 and 0.03, respectively. Both are chosen so small
that it should allow the geometry to be good enough for the real geometry
optimization, in the next step. However, these parameters are *not* given
as a recommendation and should never be used in any serious work.
  PLEASE USE AT LEAST 4~6 VASP STEPS!
                                                                 -- OVHPA
*/
}
/*******************************************/
/* Specific auto-settings for GULP: ginput */
/*******************************************/
void gulp_specific_step_ginput(FILE *there,gint step_n){
gint max_cyc, idx;
gchar *library=NULL;
if(step_n<1) return;
if(step_n>6) max_cyc=500;
else max_cyc=1000-100*step_n;
fprintf(there,"space\n");
fprintf(there,"1\n");
fprintf(there,"maxcyc %i\n",max_cyc);
/*we either take a library if one was defined or set default reaxff.lib*/
if(uspex_gui._tmp_ai_lib_sel[step_n-1]!=NULL){
	/*take the selected or defined library*/
	if(uspex_gui._tmp_ai_lib_folder[step_n-1]!=NULL)
		library=g_build_filename(uspex_gui._tmp_ai_lib_folder[step_n-1],uspex_gui._tmp_ai_lib_sel[step_n-1], NULL);
	else
		library=g_strdup(uspex_gui._tmp_ai_lib_sel[step_n-1]);
}else{
	/*see in other steps if a library was defined (take the last defined one)*/
	for(idx=uspex_gui.calc._num_opt_steps-1;(idx>=0)&&(library==NULL);idx--){
		if(uspex_gui.calc.abinitioCode[idx]!=3) continue;/*not GULP*/
		if(uspex_gui._tmp_ai_lib_folder[idx]==NULL) continue;
		else {
			/*we found the match*/
			if(uspex_gui._tmp_ai_lib_folder[idx]!=NULL)
				library=g_build_filename(uspex_gui._tmp_ai_lib_folder[idx],uspex_gui._tmp_ai_lib_sel[idx], NULL);
			else
				library=g_strdup(uspex_gui._tmp_ai_lib_sel[idx]);
			continue;/*will exit for*/
		}
	}
}
if(library!=NULL){
	fprintf(there,"library %s\n",library);
	g_free(library);
} else fprintf(there,"library reaxff.lib\n");/*the default choice*/
fprintf(there,"switch rfo 0.010\n");
}
/*********************************************/
/* Specific auto-settings for GULP: goptions */
/*********************************************/
void gulp_specific_step_goptions(FILE *there,gint step_n){
/*PSO and USPEX method are supposed to take identical steps, COPEX method is unknowned*/
if((uspex_gui.calc.calculationMethod == US_CM_VCNEB)||(uspex_gui.calc.calculationMethod == US_CM_TPS)){
/*VC-NEB & TPS*/
	if(step_n!=1) {
		gchar *title;
		title=g_strdup_printf("USPEX: VC-NEB & TPS calculation only require 1 calculation step!\nUSPEX: STEP %i IGNORED!\n",step_n);
		gui_text_show(ERROR,title);
		g_free(title);
		return;
	}
	fprintf(there,"opti prop conp stress\n");
	return;
}
if((uspex_gui.calc.calculationMethod==US_CM_META)||(uspex_gui.calc.calculationMethod==US_CM_MINHOP)){
/*META & MINHOP*/
	if(uspex_gui.calc._isfixed[step_n-1]){
		fprintf(there,"opti conj nosymmetry conv\n");
	}else{
		fprintf(there,"opti conj nosymmetry conp\n");
	}
	return;
}


if(step_n<1) return;
if(step_n==1) fprintf(there,"opti conj nosymmetry conv\n");
else fprintf(there,"opti conj nosymmetry conp\n");
}
/***************************************************/
/* create Specific directory or copy files into it */
/***************************************************/
void uspex_create_specific(uspex_exec_struct *uspex_exec){
GDir *src;
gchar *source;
gchar *target;
gchar *tmp;
gboolean hasend;
gchar *ptr,*sel;
gchar *filename;
FILE *tg;
gint step;
gint idx;
/* check if Specific directory exists*/
if(uspex_gui.calc.job_path==NULL) return;
source=g_build_filename(uspex_gui.calc.job_path,"/Specific/", NULL);
src=g_dir_open(source,0,NULL);
if(src==NULL) {
	/*if doesn't exists, create it*/
	if(g_mkdir(source,0775)){
		fprintf(stdout,"Can't create %s directory.\n",source);
		return;
	}
	src=g_dir_open(source,0,NULL);
	if(src==NULL) {
		fprintf(stdout,"Can't open newly created %s directory.\n",source);
		return;
	}
}
g_dir_close(src);
/*Specific is open, check what needs to be done*/
for(step=0;step<uspex_gui.calc._num_opt_steps;step++){
	if(uspex_gui.auto_step){
		/*use auto-generation (work in progress)*/
		switch (uspex_gui.calc.abinitioCode[step]){
		case 1:/*VASP*/
			/*1- create {Specific}/INCAR_X */
			tmp=g_strdup_printf("INCAR_%i",step+1);
			target=g_build_filename(source,tmp,NULL);
			tg=fopen(target,"wt");
			if(tg==NULL){
				fprintf(stderr,"USPEX: I/O error, can't create file %s\n",target);
				g_free(target);g_free(tmp);
				continue;/*and hope for the best*/
			}
			vasp_specific_step_incar(tg,step+1);
			fclose(tg);g_free(target);g_free(tmp);
			/*no option for VASP*/
			break;
		case 3:/*GULP*/
			/*1- create {Specific}/ginput_X */
			tmp=g_strdup_printf("ginput_%i",step+1);
			target=g_build_filename(source,tmp,NULL);
			tg=fopen(target,"wt");
			if(tg==NULL){
				fprintf(stderr,"USPEX: I/O error, can't create file %s\n",target);
				g_free(target);g_free(tmp);
				continue;/*and hope for the best*/
			}
			gulp_specific_step_ginput(tg,step+1);
			fclose(tg);g_free(target);g_free(tmp);
			/*2- create {Specific}/goptions_X */
			tmp=g_strdup_printf("goptions_%i",step+1);
			target=g_build_filename(source,tmp,NULL);
			tg=fopen(target,"wt");
			if(tg==NULL){
				fprintf(stderr,"USPEX: I/O error, can't create file %s\n",target);
				g_free(target);g_free(tmp);
				continue;/*and hope for the best*/
			}
			gulp_specific_step_goptions(tg,step+1);
			fclose(tg);g_free(target);g_free(tmp);
			break;
		case 2://SIESTA
		case 4://LAMMPS
		case 5://ORCA
		case 6://DMACRYS
		case 7://CP2K
		case 8://Q.E.
		case 9://FHI-aims
		case 10://ATK
		case 11://CASTEP
		case 12://Tinker
		case 13://MOPAC
		case 14://BoltzTraP
		case 15://DFTB
		case 16://Gaussian
		case 17://N/A
		case 18://Abinit
		case 19://CRYSTAL
		default:
			/*all of above are not ready (yet)*/
			break;
		}
	}else{
		/*we assume that everything has already been set by user*/
		/*1- copy input (if any)*/
		if(uspex_gui._tmp_ai_input[step]!=NULL){
			tg=fopen(uspex_gui._tmp_ai_input[step],"r");
			if(tg==NULL){
				fprintf(stderr,"USPEX: I/O error, can't create file %s\n",uspex_gui._tmp_ai_input[step]);
				continue;/*and hope for the best*/
			}
			fclose(tg);
			tmp=g_path_get_basename(uspex_gui._tmp_ai_input[step]);
			target=g_build_filename(source,tmp,NULL);
			if(!dumb_file_copy(uspex_gui._tmp_ai_input[step],target)){
				fprintf(stderr,"USPEX: I/O error, can't copy file %s to %s\n",uspex_gui._tmp_ai_input[step],target);
				g_free(tmp);g_free(target);
				continue;/*and hope for the best*/
			}
			g_free(tmp);g_free(target);
		}
		/*2- copy option (is any)*/
		if(uspex_gui._tmp_ai_opt[step]!=NULL){
			tg=fopen(uspex_gui._tmp_ai_opt[step],"r");
			if(tg==NULL){
				fprintf(stderr,"USPEX: I/O error, can't create file %s\n",uspex_gui._tmp_ai_opt[step]);
				continue;/*and hope for the best*/
			}
			fclose(tg);
			tmp=g_path_get_basename(uspex_gui._tmp_ai_opt[step]);
			target=g_build_filename(source,tmp,NULL);
			if(!dumb_file_copy(uspex_gui._tmp_ai_opt[step],target)){
				fprintf(stderr,"USPEX: I/O error, can't copy file %s to %s\n",uspex_gui._tmp_ai_opt[step],target);
				g_free(tmp);g_free(target);
				continue;/*and hope for the best*/
			}
			g_free(tmp);g_free(target);
		}
	}
	/*regardless of auto-generation setting, copy the libraries (when needed)*/
	if(uspex_gui._tmp_ai_lib_sel[step]!=NULL){
		switch (uspex_gui.calc.abinitioCode[step]){
		case 1:/*VASP*/
			/*special per-species*/
			if(uspex_gui._tmp_ai_lib_folder==NULL) continue;/*required*/
			sel=g_strdup(uspex_gui._tmp_ai_lib_sel[step]);
			hasend=FALSE;
			idx=0;
			while((sel[idx]==' ')||(sel[idx]=='\t')) {/*skip trailing space*/
				if(sel[idx]=='\0') {
					fprintf(stderr,"USPEX: ERROR while reading POTCAR flavors!\n");
					g_free(sel);
					sel=NULL;
					break;
				}
				idx++;
			}
			if(sel==NULL) continue;/*and hope for the best*/
			ptr=&(sel[idx]);
			while(!hasend){
				while((sel[idx]!=' ')){
					if(sel[idx]=='\0') {
						hasend=TRUE;
						break;  
					}
					idx++;
				}
				sel[idx]='\0';
				idx++;
				filename=g_strdup_printf("%s/%s/POTCAR",uspex_gui._tmp_ai_lib_folder[step],ptr);
				if(g_ascii_islower(ptr[1])) tmp=g_strdup_printf("POTCAR_%c%c",ptr[0],ptr[1]);
				else tmp=g_strdup_printf("POTCAR_%c",ptr[0]);
				target=g_build_filename(source,tmp,NULL);
				if(!dumb_file_copy(filename,target)){
					fprintf(stderr,"USPEX: I/O error, can't copy file %s to %s\n",filename,target);
					g_free(target);g_free(tmp);g_free(filename);
					continue;/*and hope for the best*/
				}
				g_free(target);g_free(tmp);g_free(filename);
				ptr=&(sel[idx]);
			}
			g_free(sel);
			break;
		case 3:/*GULP*/
			/*Since gulp can use the absolute path,*/
			/*there is no need to copy the library.*/
			break;
		case 7://CP2K
		case 12://Tinker
			//copy one library
		case 4://LAMMPS
		case 5://ORCA
		case 6://DMACRYS
		case 9://FHI-aims
		case 13://MOPAC
		case 14://BoltzTraP
		case 16://Gaussian
		case 17://N/A
		case 19://CRYSTAL
			//above LAMMPS~Gaussian either don't need
			//any library or use is unknown...
			target=g_build_filename(source,uspex_gui._tmp_ai_lib_sel[step],NULL);
			if(uspex_gui._tmp_ai_lib_folder[step]!=NULL)
				tmp=g_build_filename(uspex_gui._tmp_ai_lib_folder[step],uspex_gui._tmp_ai_lib_sel[step],NULL);
			else
				tmp=g_strdup(uspex_gui._tmp_ai_lib_sel[step]);
			if(!dumb_file_copy(tmp,target)){
				fprintf(stderr,"USPEX: I/O error, can't copy file %s to %s\n",tmp,target);
				g_free(target);g_free(tmp);
				continue;/*and hope for the best*/
			}
			g_free(target);g_free(tmp);
			break;
		case 2://SIESTA
		case 8://Q.E.
		case 10://ATK
		case 11://CASTEP
		case 18://Abinit
			//per-species potential
			break;
		case 15://DFTB
			//another special per-species
			break;
		default:
			//do nothing
			break;
		}
	}
}
g_free(source);
}
/*****************************************/
/* Execute or enqueue a uspex calculation */
/*****************************************/
void run_uspex_exec(uspex_exec_struct *uspex_exec,struct task_pak *task){
	/* Execute a uspex task TODO distant job */
	gchar *cmd;
	gchar *cwd;/*for push/pull*/

#if __WIN32
/*	At present running uspex in __WIN32 environment is impossible.
 *	However,  it should be possible to launch a (remote) job on a 
 *	distant server. See TODO */
	fprintf(stderr,"USPEX calculation can't be done in this environment.\n");return;
#else 
	/*direct launch: NEW - only use the USPEX own interface*/
	cmd = g_strdup_printf("%s -r > uspex.log",(*uspex_exec).job_uspex_exe);
#endif
	cwd=sysenv.cwd;/*push*/
	sysenv.cwd=g_strdup_printf("%s",(*uspex_exec).job_path);
	task->is_async = TRUE;
	task->status_file = g_build_filename((*uspex_exec).job_path,"uspex.log", NULL);
	task_async(cmd,&(task->pid));
	g_free(sysenv.cwd);
	sysenv.cwd=cwd;/*pull*/
	g_free(cmd);
//	(*uspex_exec).have_result=TRUE;
}
/*****************************/
/* cleanup task: load result */
/*****************************/
void cleanup_uspex_exec(uspex_exec_struct *uspex_exec){
	gchar *line;
	/*USPEX process ENDED*/
//	while(!(*uspex_exec).have_result) usleep(500*1000);/*wait for end of process*/
	line = g_strdup_printf("Completed USPEX job in directory %s\n",(*uspex_exec).job_path);
	gui_text_show(ITALIC,line);
	g_free(line);
}
/******************************/
/* Enqueue a uspex calculation */
/******************************/
void uspex_exec_calc(){
	FILE *fp;
	gint idx;
	gchar *filename;
	struct model_pak *result_model;
	uspex_exec_struct *uspex_exec;
	/*this will sync then enqueue a USPEX calculation*/
	if(save_uspex_calc()) return;/*sync and save all file*/
/*copy structure to the list*/
	uspex_exec=g_malloc(sizeof(uspex_exec_struct));
	uspex_exec->job_uspex_exe=g_strdup_printf("%s",uspex_gui.calc.job_uspex_exe);
	uspex_exec->job_path=g_strdup_printf("%s",uspex_gui.calc.job_path);
/*NEW: copy the Specific directory or create it*/
	if(uspex_gui.have_specific) uspex_copy_specific(uspex_exec);
	else uspex_create_specific(uspex_exec);
/*Get the future valid index*/
	idx=1;
	do{
		filename=g_strdup_printf("%s/results%i/OUTPUT.txt",(*uspex_exec).job_path,idx);
		fp=fopen(filename,"rt");
		if(fp==NULL) break;
		fclose(fp);
		idx++;
	}while(TRUE);
	(*uspex_exec).index=idx;
	filename=g_strdup_printf("%s/results%i/OUTPUT.txt",(*uspex_exec).job_path,(*uspex_exec).index);
	/*launch uspex in a task*/
	GUI_LOCK(uspex_gui.button_save);
	GUI_LOCK(uspex_gui.button_exec);
	task_new("USPEX", &run_uspex_exec,uspex_exec,&cleanup_uspex_exec,uspex_exec,sysenv.active_model);
	/*let's try something else.. we set a minimum model, and send it to tracking*/
	result_model=model_new();
	model_init(result_model);
	strcpy(result_model->filename,filename);
	result_model->uspex=NULL;
	result_model->track_me = TRUE;
	g_timeout_add_full(G_PRIORITY_DEFAULT,TRACKING_TIMEOUT,track_uspex,result_model,track_uspex_cleanup);
/*draw a dummy model*/
	g_free(result_model->basename);
if(GUI_ENTRY_LENGTH(uspex_gui.name)>0) GUI_ENTRY_GET_TEXT(uspex_gui.name,result_model->basename);
else result_model->basename = g_strdup("new_USPEX");
	sysenv.active_model = result_model;
	model_prep(result_model);
	result_model->mode = FREE;
	result_model->rmax = 5.0*RMAX_FUDGE;
	tree_model_add(result_model);
	tree_select_model(result_model);
	redraw_canvas(SINGLE);
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
/**************************************************/
/* refresh the GUI with value from uspex_gui.calc */
/**************************************************/
void uspex_gui_refresh(){
	gint i,j,num;
	gchar *line;
	gchar *ptr;
        /*refresh the whole GUI based on uspex_gui.calc information*/
	switch(uspex_gui.calc.calculationMethod){
		case US_CM_USPEX:
		GUI_COMBOBOX_SET(uspex_gui.calculationMethod,0);
		break;
		case US_CM_META:
		GUI_COMBOBOX_SET(uspex_gui.calculationMethod,1);
		break;
		case US_CM_VCNEB:
		GUI_COMBOBOX_SET(uspex_gui.calculationMethod,2);
		break;
		case US_CM_PSO:
		GUI_COMBOBOX_SET(uspex_gui.calculationMethod,3);
		break;
		case US_CM_TPS:
		GUI_COMBOBOX_SET(uspex_gui.calculationMethod,4);
		break;
		case US_CM_MINHOP:
		GUI_COMBOBOX_SET(uspex_gui.calculationMethod,5);
		break;
		case US_CM_COPEX:
		GUI_COMBOBOX_SET(uspex_gui.calculationMethod,6);
		break;
		case US_CM_UNKNOWN:
		default:
		GUI_COMBOBOX_SET(uspex_gui.calculationMethod,0);
	}
	switch(uspex_gui.calc.calculationType){
		case US_CT_300:
		GUI_COMBOBOX_SET(uspex_gui.calculationType,0);
		break;
		case US_CT_s300:
		GUI_COMBOBOX_SET(uspex_gui.calculationType,1);
		break;
		case US_CT_301:
		GUI_COMBOBOX_SET(uspex_gui.calculationType,2);
		break;
		case US_CT_s301:
		GUI_COMBOBOX_SET(uspex_gui.calculationType,3);
		break;
		case US_CT_310:
		GUI_COMBOBOX_SET(uspex_gui.calculationType,4);
		break;
		case US_CT_311:
		GUI_COMBOBOX_SET(uspex_gui.calculationType,5);
		break;
		case US_CT_000:
		GUI_COMBOBOX_SET(uspex_gui.calculationType,6);
		break;
		case US_CT_s000:
		GUI_COMBOBOX_SET(uspex_gui.calculationType,7);
		break;
		case US_CT_001:
		GUI_COMBOBOX_SET(uspex_gui.calculationType,8);
		break;
		case US_CT_110:
		GUI_COMBOBOX_SET(uspex_gui.calculationType,9);
		break;
		case US_CT_200:
		GUI_COMBOBOX_SET(uspex_gui.calculationType,10);
		break;
		case US_CT_s200:
		GUI_COMBOBOX_SET(uspex_gui.calculationType,11);
		break;
		case US_CT_201:
		GUI_COMBOBOX_SET(uspex_gui.calculationType,12);
		break;
		case US_CT_s201:
		GUI_COMBOBOX_SET(uspex_gui.calculationType,13);
		break;
		case US_CT_m200:
		GUI_COMBOBOX_SET(uspex_gui.calculationType,14);
		break;
		case US_CT_sm200:
		GUI_COMBOBOX_SET(uspex_gui.calculationType,15);
		break;
		case US_CT_m201:
		GUI_COMBOBOX_SET(uspex_gui.calculationType,16);
		break;
		case US_CT_sm201:
		GUI_COMBOBOX_SET(uspex_gui.calculationType,17);
		break;
		default:
		GUI_COMBOBOX_SET(uspex_gui.calculationType,0);
	}
	GUI_SPIN_SET(uspex_gui._calctype_dim,(gdouble) uspex_gui.calc._calctype_dim);
	if(uspex_gui.calc._calctype_mol) GUI_TOGGLE_ON(uspex_gui._calctype_mol);
	else GUI_TOGGLE_OFF(uspex_gui._calctype_mol);
	if(uspex_gui.calc._calctype_var) GUI_TOGGLE_ON(uspex_gui._calctype_var);
	else GUI_TOGGLE_OFF(uspex_gui._calctype_var);
	if(uspex_gui.calc._calctype_mag) {
		GUI_TOGGLE_ON(uspex_gui._calctype_mag);
		GUI_TOGGLE_ON(uspex_gui._calctype_mag_2);
	}else{
		GUI_TOGGLE_OFF(uspex_gui._calctype_mag);
		GUI_TOGGLE_OFF(uspex_gui._calctype_mag_2);
	}
	if(uspex_gui.calc._calctype_mag){
		GUI_ENTRY_SET(uspex_gui.mag_nm,uspex_gui.calc.magRatio[0],"%lf");
		GUI_ENTRY_SET(uspex_gui.mag_fmls,uspex_gui.calc.magRatio[1],"%lf");
		GUI_ENTRY_SET(uspex_gui.mag_fmhs,uspex_gui.calc.magRatio[2],"%lf");
		GUI_ENTRY_SET(uspex_gui.mag_afml,uspex_gui.calc.magRatio[3],"%lf");
		GUI_ENTRY_SET(uspex_gui.mag_afmh,uspex_gui.calc.magRatio[4],"%lf");
		GUI_ENTRY_SET(uspex_gui.mag_fmlh,uspex_gui.calc.magRatio[5],"%lf");
		GUI_ENTRY_SET(uspex_gui.mag_aflh,uspex_gui.calc.magRatio[6],"%lf");
	}
	switch(uspex_gui.calc.optType){
		case US_OT_ENTHALPY:
		GUI_COMBOBOX_SET(uspex_gui.optType,0);
		break;
		case US_OT_VOLUME:
		GUI_COMBOBOX_SET(uspex_gui.optType,1);
		break;
		case US_OT_HARDNESS:
		GUI_COMBOBOX_SET(uspex_gui.optType,2);
		break;
		case US_OT_ORDER:
		GUI_COMBOBOX_SET(uspex_gui.optType,3);
		break;
		case US_OT_DISTANCE:
		GUI_COMBOBOX_SET(uspex_gui.optType,4);
		break;
		case US_OT_DIELEC_S:
		GUI_COMBOBOX_SET(uspex_gui.optType,5);
		break;
		case US_OT_GAP:
		GUI_COMBOBOX_SET(uspex_gui.optType,6);
		break;
		case US_OT_DIELEC_GAP:
		GUI_COMBOBOX_SET(uspex_gui.optType,7);
		break;
		case US_OT_MAG:
		GUI_COMBOBOX_SET(uspex_gui.optType,8);
		break;
		case US_OT_QE:
		GUI_COMBOBOX_SET(uspex_gui.optType,9);
		break;
		case US_OT_2R:
		GUI_COMBOBOX_SET(uspex_gui.optType,10);
		break;
		case US_OT_HM:
		GUI_COMBOBOX_SET(uspex_gui.optType,11);
		break;
		case US_OT_ZT:
		GUI_COMBOBOX_SET(uspex_gui.optType,12);
		break;
		case US_OT_Fphon:
		GUI_COMBOBOX_SET(uspex_gui.optType,13);
		break;
		case US_OT_BULK_M:
		GUI_COMBOBOX_SET(uspex_gui.optType,14);
		break;
		case US_OT_SHEAR_M:
		GUI_COMBOBOX_SET(uspex_gui.optType,15);
		break;
		case US_OT_YOUNG_M:
		GUI_COMBOBOX_SET(uspex_gui.optType,16);
		break;
		case US_OT_POISSON:
		GUI_COMBOBOX_SET(uspex_gui.optType,17);
		break;
		case US_OT_PUGH_R:
		GUI_COMBOBOX_SET(uspex_gui.optType,18);
		break;
		case US_OT_VICKERS_H:
		GUI_COMBOBOX_SET(uspex_gui.optType,19);
		break;
		case US_OT_FRACTURE:
		GUI_COMBOBOX_SET(uspex_gui.optType,20);
		break;
		case US_OT_DEBYE_T:
		GUI_COMBOBOX_SET(uspex_gui.optType,21);
		break;
		case US_OT_SOUND_V:
		GUI_COMBOBOX_SET(uspex_gui.optType,22);
		break;
		case US_OT_SWAVE_V:
		GUI_COMBOBOX_SET(uspex_gui.optType,23);
		break;
		case US_OT_PWAVE_V:
		GUI_COMBOBOX_SET(uspex_gui.optType,24);
		break;
		case US_OT_UNKNOWN:
		default:
		GUI_COMBOBOX_SET(uspex_gui.optType,0);
	}
	if(uspex_gui.new_optType!=NULL){/*_BUG_ memory*/
		uspex_gui.have_new_optType=TRUE;
		if(uspex_gui._tmp_new_optType!=NULL) {
			g_free(uspex_gui._tmp_new_optType);
			uspex_gui._tmp_new_optType=NULL;
		}
		if(uspex_gui.calc.new_optType!=NULL){
			uspex_gui._tmp_new_optType=g_strdup(uspex_gui.calc.new_optType);
			GUI_ENTRY_TEXT(uspex_gui.new_optType,uspex_gui._tmp_new_optType);
		}
	}else{
		uspex_gui.have_new_optType=FALSE;
	}
	opt_toggle();
	populate_atomType();
	uspex_gui.auto_bonds=(uspex_gui.calc.goodBonds==NULL);
	if(!uspex_gui.calc._calctype_var) update_numSpecies();
	else{
		if(uspex_gui.calc.numSpecies!=NULL){
			/*update numSpecies BLOCKS*/
			GUI_UNLOCK(uspex_gui.numSpecies);
			GUI_COMBOBOX_WIPE(uspex_gui.numSpecies);
			if(uspex_gui.calc._calctype_mol) num=uspex_gui.calc._nmolecules;
			else num=uspex_gui.calc._nspecies;
			for(i=0;i<uspex_gui.calc._var_nspecies;i++){
				line=g_strdup(" ");
				for(j=0;j<num;j++){
					ptr=g_strdup_printf("%s %i",line,uspex_gui.calc.numSpecies[j+i*uspex_gui.calc._nmolecules]);
					g_free(line);
					line=ptr;
				}
				GUI_COMBOBOX_ADD(uspex_gui.numSpecies,line);
				g_free(line);
			}
			GUI_COMBOBOX_ADD(uspex_gui.numSpecies,"ADD SPECIES BLOCK");
			GUI_COMBOBOX_SET(uspex_gui.numSpecies,0);
			GUI_LOCK(uspex_gui.numSpecies);
		}
	}
set_numSpecies();/*fix _BUG_*/
	USPEX_SET_VAL(fitLimit,"%lf");
	if(uspex_gui.calc.ldaU!=NULL){
		line=g_strdup("");
		for(i=0;i<uspex_gui.calc._nspecies;i++) line=g_strdup_printf("%s %f",line,uspex_gui.calc.ldaU[i]);
		uspex_gui._tmp_ldaU=line;
		line=NULL;
	}else{
		uspex_gui._tmp_ldaU=g_strdup("");
	}
	GUI_ENTRY_TEXT(uspex_gui.ldaU,uspex_gui._tmp_ldaU);
	USPEX_SET_VAL(populationSize,"%i");
	USPEX_SET_VAL(initialPopSize,"%i");
	USPEX_SET_VAL(numGenerations,"%i");
	USPEX_SET_VAL(stopCrit,"%i");
	USPEX_SET_VAL(bestFrac,"%lf");
	USPEX_SET_VAL(keepBestHM,"%i");
	USPEX_SET_CHECK(reoptOld);
	if(uspex_gui.calc.symmetries!=NULL) USPEX_SET_TEXT(symmetries);
	USPEX_SET_VAL(fracGene,"%lf");
	USPEX_SET_VAL(fracRand,"%lf");
	USPEX_SET_VAL(fracTopRand,"%lf");
	USPEX_SET_VAL(fracPerm,"%lf");
	USPEX_SET_VAL(fracAtomsMut,"%lf");
	USPEX_SET_VAL(fracRotMut,"%lf");
	USPEX_SET_VAL(fracLatMut,"%lf");
	USPEX_SET_VAL(fracSpinMut,"%lf");
	USPEX_SET_VAL(howManySwaps,"%i");
	if(uspex_gui.calc.specificSwaps!=NULL) USPEX_SET_TEXT(specificSwaps);
	USPEX_SET_VAL(mutationDegree,"%lf");
	USPEX_SET_VAL(mutationRate,"%lf");
	USPEX_SET_VAL(DisplaceInLatmutation,"%lf");
	USPEX_SET_CHECK(AutoFrac);
/*TRY*/
	uspex_gui.auto_C_ion=(uspex_gui.calc.IonDistances==NULL);
	toggle_auto_C_ion();
	USPEX_SET_VAL(minVectorLength,"%lf");
	USPEX_SET_VAL(constraint_enhancement,"%i");
	refresh_constraints();
	uspex_gui.auto_C_lat=(uspex_gui.calc.Latticevalues==NULL);
	toggle_auto_C_lat();
	set_lattice_format();
	if((uspex_gui.calc._nsplits>0)&&(uspex_gui.calc.splitInto!=NULL)){
		line=g_strdup_printf("%i",uspex_gui.calc.splitInto[0]);
		for(i=1;i<uspex_gui.calc._nsplits;i++){
			ptr=g_strdup_printf("%s %i",line,uspex_gui.calc.splitInto[i]);
			g_free(line);
			line=ptr;
		}
		GUI_ENTRY_TEXT(uspex_gui.splitInto,line);
		g_free(line);
	}
	USPEX_SET_CHECK(pickUpYN);
	USPEX_SET_VAL(pickUpGen,"%i");
	USPEX_SET_VAL(pickUpFolder,"%i");
/*NEW: SET SPECIFIC*/
	if(uspex_gui.have_specific){
		toggle_use_specific();
		GUI_ENTRY_TEXT(uspex_gui.spe_folder,uspex_gui._tmp_spe_folder);
	}else{
		toggle_set_specific();
	}
num=(gint)uspex_gui._tmp_num_opt_steps;/*this is the PREVIOUS _num_opt_steps*/
	/*if we don't have abinitCode, set a default*/
	if(uspex_gui.calc.abinitioCode==NULL) {
		uspex_gui.calc.abinitioCode=g_malloc(uspex_gui.calc._num_opt_steps*sizeof(gint));
		for(i=1;i<uspex_gui.calc._num_opt_steps;i++) uspex_gui.calc.abinitioCode[i]=1;
	}
	/*same with KresolStart*/
	if(uspex_gui.calc.KresolStart==NULL) {
		uspex_gui.calc.KresolStart=g_malloc(uspex_gui.calc._num_opt_steps*sizeof(gdouble));
		uspex_gui.calc.KresolStart[0]=0.2;
		for(i=1;i<uspex_gui.calc._num_opt_steps;i++)
			uspex_gui.calc.KresolStart[i]=0.2-(gdouble)(i)*(0.2-0.08)/(uspex_gui.calc._num_opt_steps-1);
	}
	/*same with vacuumSize*/
	if(uspex_gui.calc.vacuumSize==NULL) {
		uspex_gui.calc.vacuumSize=g_malloc(uspex_gui.calc._num_opt_steps*sizeof(gdouble));
		for(i=0;i<uspex_gui.calc._num_opt_steps;i++) uspex_gui.calc.vacuumSize[i]=10.;
	}
	/*(re-)create _tmp_commandExecutable _tmp_ai_input _tmp_ai_opt*/
	if(uspex_gui._tmp_commandExecutable!=NULL) {
		for(i=0;i<num;i++){
			if(uspex_gui._tmp_commandExecutable[i]!=NULL) g_free(uspex_gui._tmp_commandExecutable[i]);
			uspex_gui._tmp_commandExecutable[i]=NULL;
		}
		g_free(uspex_gui._tmp_commandExecutable);
		uspex_gui._tmp_commandExecutable=NULL;
	}
	uspex_gui._tmp_commandExecutable=g_malloc0(uspex_gui.calc._num_opt_steps*sizeof(gchar *));
	if(uspex_gui._tmp_ai_input!=NULL) {
		for(i=0;i<num;i++) {
			if(uspex_gui._tmp_ai_input[i]!=NULL) g_free(uspex_gui._tmp_ai_input[i]);
			uspex_gui._tmp_ai_input[i]=NULL;
		}
		g_free(uspex_gui._tmp_ai_input);
		uspex_gui._tmp_ai_input=NULL;
	}
	uspex_gui._tmp_ai_input=g_malloc0(uspex_gui.calc._num_opt_steps*sizeof(gchar *));
	if(uspex_gui._tmp_ai_opt!=NULL){
		for(i=0;i<num;i++) {
			if(uspex_gui._tmp_ai_opt[i]!=NULL) g_free(uspex_gui._tmp_ai_opt[i]);
			uspex_gui._tmp_ai_opt[i]=NULL;
		}
		g_free(uspex_gui._tmp_ai_opt);
		uspex_gui._tmp_ai_opt=NULL;
	}
	uspex_gui._tmp_ai_opt=g_malloc0(uspex_gui.calc._num_opt_steps*sizeof(gchar *));
	/*redefine commandExecutable*/
	i=0;
	if(uspex_gui.calc.commandExecutable!=NULL) {
		gint size;
		line=&(uspex_gui.calc.commandExecutable[0]);
		ptr=line;
		while(*line!='\0'){
			while(!g_ascii_isprint(*ptr)) ptr++;
			line=ptr;
			size=0;
			while(g_ascii_isprint(*ptr)) {
				size++;
				ptr++;
			}
			uspex_gui._tmp_commandExecutable[i]=g_strndup(line,size);
			line=ptr;
			i++;
		}
	}
	uspex_gui.calc._isCmdList=(i>0);
	uspex_gui._tmp_num_opt_steps=(gdouble)uspex_gui.calc._num_opt_steps;
	GUI_SPIN_SET(uspex_gui._num_opt_steps,uspex_gui._tmp_num_opt_steps);
	uspex_gui._tmp_curr_step=1.;
	GUI_SPIN_SET(uspex_gui._curr_step,uspex_gui._tmp_curr_step);
	spin_update_curr_step();
	uspex_gui.auto_step=TRUE;
	toggle_auto_step();
        if(uspex_gui.calc.numProcessors==NULL){
                line=g_strdup("");
        }else{
                line=g_strdup_printf("%i",uspex_gui.calc.numProcessors[0]);
                for(i=1;i<uspex_gui.calc._num_opt_steps;i++){
                        ptr=g_strdup_printf("%s %i",line,uspex_gui.calc.numProcessors[i]);
                        g_free(line);
                        line=ptr;
                }
        }
        GUI_ENTRY_TEXT(uspex_gui.numProcessors,line);
        g_free(line);
	USPEX_SET_VAL(numParallelCalcs,"%i");
	USPEX_SET_VAL(whichCluster,"%i");
	if(uspex_gui.calc.remoteFolder!=NULL) USPEX_SET_TEXT(remoteFolder);
	USPEX_SET_CHECK(PhaseDiagram);
	USPEX_SET_VAL(RmaxFing,"%lf");
	USPEX_SET_VAL(deltaFing,"%lf");
	USPEX_SET_VAL(sigmaFing,"%lf");
	USPEX_SET_VAL(toleranceFing,"%lf");
	USPEX_SET_VAL(antiSeedsActivation,"%i");
	USPEX_SET_VAL(antiSeedsMax,"%lf");
	USPEX_SET_VAL(antiSeedsSigma,"%lf");
	USPEX_SET_CHECK(doSpaceGroup);
	SG_toggle();
	USPEX_SET_VAL(SymTolerance,"%lf");
	USPEX_SET_VAL(repeatForStatistics,"%i");
	USPEX_SET_VAL(stopFitness,"%lf");
	USPEX_SET_VAL(fixRndSeed,"%i");
	USPEX_SET_CHECK(collectForces);
	USPEX_SET_CHECK(ordering_active);
	USPEX_SET_CHECK(symmetrize);
	if(uspex_gui.calc.valenceElectr==NULL){
		line=g_strdup("");
	}else{
		line=g_strdup_printf("%i",uspex_gui.calc.valenceElectr[0]);
		for(i=1;i<uspex_gui.calc._nspecies;i++){
			ptr=g_strdup_printf("%s %i",line,uspex_gui.calc.valenceElectr[i]);
			g_free(line);
			line=ptr;
		}
	}
	GUI_ENTRY_TEXT(uspex_gui.valenceElectr,line);
	g_free(line);
	USPEX_SET_VAL(percSliceShift,"%lf");
	GUI_COMBOBOX_SET(uspex_gui.dynamicalBestHM,uspex_gui.calc.dynamicalBestHM);
	if(uspex_gui.calc.softMutOnly!=NULL) USPEX_SET_TEXT(softMutOnly);
	USPEX_SET_VAL(maxDistHeredity,"%lf");
	GUI_COMBOBOX_SET(uspex_gui.manyParents,uspex_gui.calc.manyParents);
	USPEX_SET_VAL(minSlice,"%lf");
	USPEX_SET_VAL(maxSlice,"%lf");
	USPEX_SET_VAL(numberparents,"%i");
	USPEX_SET_VAL(BoltzTraP_T_max,"%lf");
	USPEX_SET_VAL(BoltzTraP_T_delta,"%lf");
	USPEX_SET_VAL(BoltzTraP_T_efcut,"%lf");
	USPEX_SET_VAL(TE_T_interest,"%lf");
	USPEX_SET_VAL(TE_threshold,"%lf");
	uspex_gui.have_ZT=FALSE;
	if(uspex_gui.calc.optType==US_OT_ZT){
		uspex_gui.have_ZT=TRUE;
	}else{
		if(uspex_gui.calc.new_optType!=NULL){
			if(find_in_string("ZT",uspex_gui.calc.new_optType) != NULL) uspex_gui.have_ZT=TRUE;
			if(find_in_string("14",uspex_gui.calc.new_optType) != NULL) uspex_gui.have_ZT=TRUE;
		}
	}
	switch(uspex_gui.calc.TE_goal){
		case US_BT_ZTxx:
			GUI_COMBOBOX_SET(uspex_gui.TE_goal,1);
			break;
		case US_BT_ZTyy:
			GUI_COMBOBOX_SET(uspex_gui.TE_goal,2);
			break;
		case US_BT_ZTzz:
			GUI_COMBOBOX_SET(uspex_gui.TE_goal,3);
			break;
		case US_BT_ZT:
			/* fall through */
		default:
			GUI_COMBOBOX_SET(uspex_gui.TE_goal,0);
	}
	//cmd_BoltzTraP is left alone (for future use)
	USPEX_SET_VAL(thicknessS,"%lf");
	USPEX_SET_VAL(thicknessB,"%lf");
	USPEX_SET_VAL(reconstruct,"%i");
	if(uspex_gui.calc.StoichiometryStart==NULL){
		line=g_strdup("");
	}else{
		line=g_strdup_printf("%i",uspex_gui.calc.StoichiometryStart[0]);
		for(i=1;i<uspex_gui.calc._nspecies;i++){
			ptr=g_strdup_printf("%s %i",line,uspex_gui.calc.StoichiometryStart[i]);
			g_free(line);
			line=ptr;
		}
	}
	GUI_ENTRY_TEXT(uspex_gui.StoichiometryStart,line);
	g_free(line);
	USPEX_SET_VAL(E_AB,"%lf");
	USPEX_SET_VAL(Mu_A,"%lf");
	USPEX_SET_VAL(Mu_B,"%lf");
//	init_substrate_models();/*this is un-necessary*/
	USPEX_SET_VAL(firstGeneMax,"%i");
	USPEX_SET_VAL(minAt,"%i");
	USPEX_SET_VAL(maxAt,"%i");
	USPEX_SET_VAL(fracTrans,"%lf");
	USPEX_SET_VAL(howManyTrans,"%lf");
	if(uspex_gui.calc._nspetrans<1){
		line=g_strdup("");
	}else{
		line=g_strdup_printf("%i",uspex_gui.calc.specificTrans[0]);
		for(i=1;i<uspex_gui.calc._nspetrans;i++){
			ptr=g_strdup_printf("%s %i",line,uspex_gui.calc.specificTrans[i]);
			g_free(line);
			line=ptr;
		}
	}
	GUI_ENTRY_TEXT(uspex_gui.specificTrans,line);
	g_free(line);
	USPEX_SET_VAL(ExternalPressure,"%lf");
	USPEX_SET_VAL(GaussianWidth,"%lf");
	USPEX_SET_VAL(GaussianHeight,"%lf");
	GUI_COMBOBOX_SET(uspex_gui.FullRelax,uspex_gui.calc.FullRelax);
//	init_metadynamics_models();/*this is un-necessary*/
	USPEX_SET_VAL(maxVectorLength,"%lf");
	USPEX_SET_VAL(PSO_softMut,"%lf");
	USPEX_SET_VAL(PSO_BestStruc,"%lf");
	USPEX_SET_VAL(PSO_BestEver,"%lf");
	USPEX_SET_VAL(vcnebType,"%i");
	GUI_COMBOBOX_SET(uspex_gui._vcnebtype_method,uspex_gui.calc._vcnebtype_method-1);
	USPEX_SET_CHECK(_vcnebtype_img_num);
	USPEX_SET_CHECK(_vcnebtype_spring);
	USPEX_SET_VAL(numImages,"%i");
	USPEX_SET_VAL(numSteps,"%i");
	GUI_COMBOBOX_SET(uspex_gui.optReadImages,uspex_gui.calc.optReadImages);
	GUI_COMBOBOX_SET(uspex_gui.optimizerType,uspex_gui.calc.optimizerType-1);
	GUI_COMBOBOX_SET(uspex_gui.optRelaxType,uspex_gui.calc.optRelaxType-1);
	USPEX_SET_VAL(dt,"%lf");
	USPEX_SET_VAL(ConvThreshold,"%lf");
	USPEX_SET_VAL(VarPathLength,"%lf");
	USPEX_SET_VAL(K_min,"%lf");
	USPEX_SET_VAL(K_max,"%lf");
	USPEX_SET_VAL(Kconstant,"%lf");
	USPEX_SET_CHECK(optFreezing);
	GUI_COMBOBOX_SET(uspex_gui.optMethodCIDI,uspex_gui.calc.optMethodCIDI);
	USPEX_SET_VAL(startCIDIStep,"%i");
	if(uspex_gui.calc._npickimg<1){
		line=g_strdup("");
	}else{
		line=g_strdup_printf("%i",uspex_gui.calc.pickupImages[0]);
		for(i=1;i<uspex_gui.calc._npickimg;i++){
			ptr=g_strdup_printf("%s %i",line,uspex_gui.calc.pickupImages[i]);
			g_free(line);
			line=ptr;
		}
	}
	GUI_ENTRY_TEXT(uspex_gui.pickupImages,line);
	g_free(line);
	switch(uspex_gui.calc.FormatType){
		case 1:
		GUI_COMBOBOX_SET(uspex_gui.FormatType,0);
		break;
		case 3:
		GUI_COMBOBOX_SET(uspex_gui.FormatType,2);
		break;
		case 2:
			/* fall through */
		default:
		GUI_COMBOBOX_SET(uspex_gui.FormatType,1);
	}
//	init_img_gdis_model();/*this is un-necessary*/
	USPEX_SET_VAL(PrintStep,"%i");
	USPEX_SET_VAL(numIterations,"%i");
	if(uspex_gui.calc.speciesSymbol!=NULL) USPEX_SET_TEXT(speciesSymbol);
	if(uspex_gui.calc.mass==NULL){
		line=g_strdup("");
	}else{
		line=g_strdup_printf("%lf",uspex_gui.calc.mass[0]);
		for(i=1;i<uspex_gui.calc._nspecies;i++){
			ptr=g_strdup_printf("%s %lf",line,uspex_gui.calc.mass[i]);
			g_free(line);
			line=ptr;
		}
	}
	GUI_ENTRY_TEXT(uspex_gui.mass,line);
	g_free(line);
	GUI_ENTRY_SET(uspex_gui.amplitudeShoot_AB,uspex_gui.calc.amplitudeShoot[0],"%lf");
	GUI_ENTRY_SET(uspex_gui.amplitudeShoot_BA,uspex_gui.calc.amplitudeShoot[1],"%lf");
	GUI_ENTRY_SET(uspex_gui.magnitudeShoot_success,uspex_gui.calc.magnitudeShoot[0],"%lf");
	GUI_ENTRY_SET(uspex_gui.magnitudeShoot_failure,uspex_gui.calc.magnitudeShoot[1],"%lf");
	USPEX_SET_VAL(shiftRatio,"%lf");
	USPEX_SET_CHECK(orderParaType);
	GUI_ENTRY_SET(uspex_gui.opCriteria_end,uspex_gui.calc.opCriteria[1],"%lf");
	if(uspex_gui.calc.cmdOrderParameter!=NULL) USPEX_SET_TEXT(cmdOrderParameter);
	if(uspex_gui.calc.cmdEnthalpyTemperature!=NULL) USPEX_SET_TEXT(cmdEnthalpyTemperature);
	if(uspex_gui.calc.orderParameterFile!=NULL) USPEX_SET_TEXT(orderParameterFile);
	if(uspex_gui.calc.enthalpyTemperatureFile) USPEX_SET_TEXT(enthalpyTemperatureFile);
	if(uspex_gui.calc.trajectoryFile) USPEX_SET_TEXT(trajectoryFile);
	if(uspex_gui.calc.MDrestartFile) USPEX_SET_TEXT(MDrestartFile);
/*ADDON: molecular*/
	//mol_model <- nothing to do
	if(uspex_gui.calc._nmolecules < 1) uspex_gui.calc._nmolecules=1;
	uspex_gui._tmp_num_mol=(gdouble)uspex_gui.calc._nmolecules;
	uspex_gui._tmp_curr_mol=0.;
	GUI_SPIN_SET(uspex_gui.num_mol,uspex_gui._tmp_num_mol);
	GUI_SPIN_SET(uspex_gui.curr_mol,uspex_gui._tmp_curr_mol);
	if(uspex_gui._tmp_mols_gdis!=NULL) g_free(uspex_gui._tmp_mols_gdis);
	uspex_gui._tmp_mols_gdis=g_malloc0(uspex_gui.calc._nmolecules*sizeof(gint));
	if(uspex_gui._tmp_mols_gulp!=NULL) g_free(uspex_gui._tmp_mols_gulp);
	uspex_gui._tmp_mols_gulp=g_malloc0(uspex_gui.calc._nmolecules*sizeof(gboolean));
	uspex_gui.mol_as_gulp=FALSE;
	init_uspex_gdis_mol();/*this is un-necessary*/
	GUI_COMBOBOX_SET(uspex_gui.mol_model,0);
	spin_update_curr_mol();
/*SPECIFIC*/
	uspex_gui._tmp_nspecies=uspex_gui.calc._nspecies;/*BACKUP*/
	update_specific();
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
	/* special */
	uspex_output_struct *uspex_output;
	struct model_pak *data;
	gint idx;
	gchar *tmp;
/* checks */
	data = sysenv.active_model;
	if (!data) return;
/* do we have a Parameters.txt model?  */
	if(data->uspex!=NULL) uspex_output=(uspex_output_struct *)data->uspex;
	else uspex_output=NULL;
	if((uspex_output!=NULL)&&(data->id == USPEX)) {/*data->id check necessary?*/
		/*take GUI default from the opened USPEX output*/
		init_uspex_parameters(&(uspex_gui.calc));
		if(uspex_output->calc!=NULL){
			uspex_gui.have_output=TRUE;
			uspex_gui.have_specific=TRUE;
			copy_uspex_parameters(uspex_output->calc,&(uspex_gui.calc));
			gui_uspex_init(data);
		}else{
			uspex_gui.have_output=FALSE;
			uspex_gui.have_specific=FALSE;
			free_uspex_parameters(&(uspex_gui.calc));
			gui_uspex_init(data);
		}
	} else {
		uspex_gui.have_output=FALSE;
		uspex_gui.have_specific=FALSE;
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
/*update job_path <- default to current directory!*/
	if(uspex_gui.calc.job_path!=NULL) g_free(uspex_gui.calc.job_path);
	uspex_gui.calc.job_path=g_strdup_printf("%s",sysenv.cwd);
/* create frame in box */
	GUI_FRAME_WINDOW(uspex_gui.window,frame);
/* create notebook in frame */
	GUI_NOTE_FRAME(frame,notebook);

/*------------------*/
/* page 1 -> SYSTEM */
/*------------------*/
	GUI_PAGE_NOTE(notebook,page,"SYSTEM");
/* create a table in the page*/
	GUI_TABLE_NOTE(page,table,18,4);
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
	GUI_COMBOBOX_ADD(uspex_gui.calculationMethod,"COPEX");
	GUI_SPIN_TABLE(table,uspex_gui._calctype_dim,uspex_gui._dim,spin_update_dim,"DIM",1,2,1,2);
	GUI_SPIN_RANGE(uspex_gui._calctype_dim,-2.,3.);
GUI_TOOLTIP(uspex_gui._calctype_dim,"Set the dimension of the system.\n3, 2, 1, and 0 are equilvalent to 3D, 2D, 1D, and 0D.\n-2 correspond to 2D crystals.");
	GUI_ENTRY_TABLE(table,uspex_gui.ExternalPressure,uspex_gui.calc.ExternalPressure,"%.4f","ExtP:",2,3,1,2);
GUI_TOOLTIP(uspex_gui.ExternalPressure,"ExternalPressure: Ch. 4.1, 5.7 DEFAULT: none\nExternal pressure (GPa) for calculation.");
	GUI_CHECK_TABLE(table,uspex_gui.sel_v1030_1,uspex_gui.have_v1030,toggle_v1030,"Ver 10.4",3,4,1,2);
GUI_TOOLTIP(uspex_gui.sel_v1030_1,"Use USPEX v. 10.4 instead of 9.4.4.");
/* line 2 */
GUI_TOOLTIP(uspex_gui.calculationMethod,"calculationMethod: Ch. 4.1 DEFAULT: USPEX\nSet the method of calculation.");
	GUI_COMBOBOX_TABLE(table,uspex_gui.calculationType,"Type: ",0,1,2,3);
	GUI_COMBOBOX_ADD(uspex_gui.calculationType,"300");
	GUI_COMBOBOX_ADD(uspex_gui.calculationType,"s300");
	GUI_COMBOBOX_ADD(uspex_gui.calculationType,"301");
	GUI_COMBOBOX_ADD(uspex_gui.calculationType,"s301");
	GUI_COMBOBOX_ADD(uspex_gui.calculationType,"310");
	GUI_COMBOBOX_ADD(uspex_gui.calculationType,"311");
	GUI_COMBOBOX_ADD(uspex_gui.calculationType,"000");
	GUI_COMBOBOX_ADD(uspex_gui.calculationType,"s000");
	GUI_COMBOBOX_ADD(uspex_gui.calculationType,"001"); /*available only in v10.1?*/
	GUI_COMBOBOX_ADD(uspex_gui.calculationType,"110");
	GUI_COMBOBOX_ADD(uspex_gui.calculationType,"200");
	GUI_COMBOBOX_ADD(uspex_gui.calculationType,"s200");
	GUI_COMBOBOX_ADD(uspex_gui.calculationType,"201");
	GUI_COMBOBOX_ADD(uspex_gui.calculationType,"s201");
	GUI_COMBOBOX_ADD(uspex_gui.calculationType,"-200");
	GUI_COMBOBOX_ADD(uspex_gui.calculationType,"-s200");
	GUI_COMBOBOX_ADD(uspex_gui.calculationType,"-201");
	GUI_COMBOBOX_ADD(uspex_gui.calculationType,"-s201");
GUI_TOOLTIP(uspex_gui.calculationType,"calculationType: Ch. 4.1 DEFAULT: 300\nSets the dimensionality, molecularity, and variability\nof the calculation, can also be set individually.");
	GUI_CHECK_TABLE(table,uspex_gui._calctype_mag,uspex_gui.calc._calctype_mag,mag_toggle,"MAG",1,2,2,3);
GUI_TOOLTIP(uspex_gui._calctype_mag,"Set a magnetic calculation.");
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
	GUI_COMBOBOX_TABLE(table,uspex_gui.numSpecies,"numSpecies:",0,1,5,6);
	GUI_COMBOBOX_ADD(uspex_gui.numSpecies,"ADD SPECIES BLOCK");
GUI_TOOLTIP(uspex_gui.numSpecies,"numSpecies: Ch. 4.1 DEFAULT: none\nSpecifies the number of atoms of each types.\nCan be use to set blockSpecies for variable composition.");
	GUI_TEXT_TABLE(table,uspex_gui.blockSpecies,uspex_gui._tmp_blockSpecies,"Species: ",1,2,5,6);
GUI_TOOLTIP(uspex_gui.blockSpecies,"The number of atoms of each types for this block.");
//	GUI_2BUTTONS_TABLE(table,apply_block_species,delete_block_species,3,4,5,6);
	GUI_APPLY_BUTTON_TABLE(table,uspex_gui.Species_apply_button,apply_block_species,2,3,5,6);
	GUI_DELETE_BUTTON_TABLE(table,uspex_gui.Species_delete_button,delete_block_species,3,4,5,6);
/* line 6 */
	GUI_COMBOBOX_TABLE(table,uspex_gui.goodBonds,"goodBonds:",0,1,6,7);/*multiline*/
	GUI_COMBOBOX_ADD(uspex_gui.goodBonds,"ADD GOODBOND");
GUI_TOOLTIP(uspex_gui.goodBonds,"goodBonds: Ch. 4.1 DEFAULT: auto\nSet the minimum distance at which a bond is considered.");
	GUI_TEXT_TABLE(table,uspex_gui._bond_d,uspex_gui._tmp_bond_d,"Bonds: ",1,2,6,7);
GUI_TOOLTIP(uspex_gui._bond_d,"Minimum bond distance between selected and others species.");
	GUI_APPLY_BUTTON_TABLE(table,button,apply_bonds,2,3,6,7);
	GUI_DELETE_BUTTON_TABLE(table,button,remove_bonds,3,4,6,7);
/* line 7 */
	GUI_TEXT_TABLE(table,uspex_gui.ldaU,uspex_gui._tmp_ldaU,"lda+U:",0,3,7,8);
GUI_TOOLTIP(uspex_gui.ldaU,"ldaU: Ch. 4.1 DEFAULT: all 0\nHubbard U value (per atom) in L(S)DA+U method.");
	GUI_CHECK_TABLE(table,button,uspex_gui.auto_bonds,auto_bond_toggle,"AUTO_BONDS",3,4,7,8);
GUI_TOOLTIP(button,"Automatically determine bonds (recommended).");
/* line 8 */
        GUI_COMBOBOX_TABLE(table,uspex_gui.optType,"optType:",0,1,8,9);/*multiline*/
        GUI_COMBOBOX_ADD(uspex_gui.optType,"1: MIN Enthalpy (stable phases)");
        GUI_COMBOBOX_ADD(uspex_gui.optType,"2: MIN Volume (densest structure)");
        GUI_COMBOBOX_ADD(uspex_gui.optType,"3: MAX Hardness (hardest phase)");
        GUI_COMBOBOX_ADD(uspex_gui.optType,"4: MAX Order (most order structure)");
        GUI_COMBOBOX_ADD(uspex_gui.optType,"5: MAX Density");
        GUI_COMBOBOX_ADD(uspex_gui.optType,"6: MAX Dielectric susceptibility");
        GUI_COMBOBOX_ADD(uspex_gui.optType,"7: MAX Band gap");
        GUI_COMBOBOX_ADD(uspex_gui.optType,"8: MAX electric energy storage capacity");
        GUI_COMBOBOX_ADD(uspex_gui.optType,"9: MAX Magnetization");
        GUI_COMBOBOX_ADD(uspex_gui.optType,"10: MAX Structure quasientropy");
	GUI_COMBOBOX_ADD(uspex_gui.optType,"11: MAX L/H eigenvalue difference of refractive index");
	GUI_COMBOBOX_ADD(uspex_gui.optType,"12: MAX halfmetalicity parameter");
	GUI_COMBOBOX_ADD(uspex_gui.optType,"14: MAX ZT thermoelectric figure of merit");
	GUI_COMBOBOX_ADD(uspex_gui.optType,"17: MAX free energy at finite temperature");
        GUI_COMBOBOX_ADD(uspex_gui.optType,"1101: MAX Bulk modulus");
        GUI_COMBOBOX_ADD(uspex_gui.optType,"1102: MAX Shear modulus");
        GUI_COMBOBOX_ADD(uspex_gui.optType,"1103: MAX Young modulus");
        GUI_COMBOBOX_ADD(uspex_gui.optType,"1104: MAX Poisson ratio");
        GUI_COMBOBOX_ADD(uspex_gui.optType,"1105: MAX Pugh modulus ratio");
        GUI_COMBOBOX_ADD(uspex_gui.optType,"1106: MAX Vickers hardness");
        GUI_COMBOBOX_ADD(uspex_gui.optType,"1107: MAX Fracture toughness");
        GUI_COMBOBOX_ADD(uspex_gui.optType,"1108: MAX Debye temperature");
        GUI_COMBOBOX_ADD(uspex_gui.optType,"1109: MAX Sound velocity");
        GUI_COMBOBOX_ADD(uspex_gui.optType,"1110: MAX S-wave velocity");
        GUI_COMBOBOX_ADD(uspex_gui.optType,"1111: MAX P-wave velocity");
GUI_TOOLTIP(uspex_gui.optType,"optType: Ch. 4.1 DEFAULT: 1(Enthalpy)\nSelect the properties to optimize.");
	GUI_TEXT_TABLE(table,uspex_gui.new_optType,uspex_gui._tmp_new_optType,"NEW optType: ",1,3,8,9);
GUI_TOOLTIP(uspex_gui.new_optType,"optType: Ch. 4.1 DEFAULT: MIN_enthalpy\nSelect the properties to optimize in the new\nVER 10.1 USPEX format.");
	GUI_CHECK_TABLE(table,uspex_gui.sel_new_opt,uspex_gui.have_new_optType,opt_toggle,"NEW",3,4,8,9);
GUI_TOOLTIP(uspex_gui.sel_new_opt,"Use the NEW block optType format.\nMandatory for multiobjective optimization.");
	/*TODO: check on "change" to detect ZT*/
/* line 9 */
	/*col 1: empty*/
	GUI_CHECK_TABLE(table,button,uspex_gui.calc.anti_opt,NULL,"ANTI-OPT",1,2,9,10);/*not calling anything*/
GUI_TOOLTIP(button,"anti-opt: - DEFAULT: FALSE\nIf set REVERSE the direction of optimization.\ni.e. MIN -> MAX & MAX -> MIN");
	GUI_CHECK_TABLE(table,uspex_gui.checkMolecules,uspex_gui.calc.checkMolecules,NULL,"ckMol",2,3,9,10);/*not calling anything*/
GUI_TOOLTIP(uspex_gui.checkMolecules,"checkMolecules: Ch. 4.1 DEFAULT: TRUE\nCheck and discard broken/merged molecules.");
	GUI_CHECK_TABLE(table,uspex_gui.checkConnectivity,uspex_gui.calc.checkConnectivity,NULL,"ckCon",3,4,9,10);/*not calling anything*/
GUI_TOOLTIP(uspex_gui.checkConnectivity,"checkConnectivity: Ch. 4.1 DEFAULT: FALSE\nCalculate hardness and add connectivity in softmutation.");
/* <- Ext Pressure, LDA+U, and magRation moved to p. III, II, and III, respectively.*/
/* --- cell */
	GUI_LABEL_TABLE(table,"Cell",0,4,10,11);
/* line 11 */
        GUI_COMBOBOX_TABLE(table,uspex_gui.Latticevalues,"Lattice:",0,1,11,12);/*multiline*/
GUI_TOOLTIP(uspex_gui.Latticevalues,"Latticevalues: Ch. 4.6 DEFAULT: auto\nInitial volume of the unit cell, or known lattice parameters.");
        GUI_TEXT_TABLE(table,uspex_gui._latticevalue,uspex_gui._tmp_latticevalue,"VALUES:",1,3,11,12);
GUI_TOOLTIP(uspex_gui._latticevalue,"VALUES: values of the current line of Latticevalues.");
        GUI_APPLY_BUTTON_TABLE(table,button,apply_latticevalue,3,4,11,12);
/* line 12 */
        GUI_COMBOBOX_TABLE(table,uspex_gui._latticeformat,"FORMAT:",0,1,12,13);
	/*col 2: empty*/
        GUI_COMBOBOX_ADD(uspex_gui._latticeformat,"Volumes");
        GUI_COMBOBOX_ADD(uspex_gui._latticeformat,"Lattice");
        GUI_COMBOBOX_ADD(uspex_gui._latticeformat,"Crystal");
GUI_TOOLTIP(uspex_gui._latticeformat,"FORMAT: whether Lattice values correspond to a series of\nVolumes, lattive vectors, or crystallographic definition.");
	title=g_strdup_printf("%i",uspex_gui.calc.splitInto[0]);
	for(idx=1;idx<uspex_gui.calc._nsplits;idx++){
		tmp=g_strdup_printf("%s %i",title,uspex_gui.calc.splitInto[idx]);
		g_free(title);
		title=tmp;
	}
        GUI_TEXT_TABLE(table,uspex_gui.splitInto,title,"split:",1,3,12,13);
	g_free(title);
GUI_TOOLTIP(uspex_gui.splitInto,"splitInto: Ch. 4.4 DEFAULT: 1\nNumber of identical subcells or pseudosubcells in the unitcell.");
        GUI_CHECK_TABLE(table,button,uspex_gui.auto_C_lat,toggle_auto_C_lat,"AUTO_LAT",3,4,12,13);
GUI_TOOLTIP(button,"AUTO_LAT: use automatic value for Latticevalues.");
/* Constraints */
	GUI_LABEL_TABLE(table,"Constraints",0,4,13,14);
/* line 14 */
	GUI_COMBOBOX_TABLE(table,uspex_gui.IonDistances,"Ion:",0,1,14,15);/*multiline*/
GUI_TOOLTIP(uspex_gui.IonDistances,"IonDistance: Ch. 4.5 DEFAULT: auto\nTriangular matrix of the minimum allowed\ninteratomic distances.");
	GUI_TEXT_TABLE(table,uspex_gui._distances," ","DIST: ",1,3,14,15);
GUI_TOOLTIP(uspex_gui._distances,"distances: minimum allowed distance from the current atom to others.");
	GUI_APPLY_BUTTON_TABLE(table,button,apply_distances,3,4,14,15);
/* line 15 */
	/*col 1: empty*/
	GUI_ENTRY_TABLE(table,uspex_gui.minVectorLength,uspex_gui.calc.minVectorLength,"%.4f","MV:",1,2,15,16);
GUI_TOOLTIP(uspex_gui.minVectorLength,"minVectorLength: Ch. 4.5, 5.7 DEFAULT: auto\nMinimum length of a new lattice parameter.");
	GUI_ENTRY_TABLE(table,uspex_gui.constraint_enhancement,uspex_gui.calc.constraint_enhancement,"%i","CE:",2,3,15,16);
GUI_TOOLTIP(uspex_gui.constraint_enhancement,"constraint_enhancement: Ch. 4.5 DEFAULT: 1\nApply CE times the IonDistance constraints.");
	GUI_CHECK_TABLE(table,button,uspex_gui.auto_C_ion,toggle_auto_C_ion,"AUTO_ION",3,4,15,16);
GUI_TOOLTIP(button,"AUTO_ION: use automatic values for ion constraints.");
/* line 16 */
	GUI_COMBOBOX_TABLE(table,uspex_gui.MolCenters,"Mol:",0,1,16,17);/*multiline*/
GUI_TOOLTIP(uspex_gui.MolCenters,"MolCenters: Ch. 4.5 DEFAULT: none\nTriangular matrix of the minimum allowed\ndistances between molecules centers.");
	GUI_TEXT_TABLE(table,uspex_gui._centers," ","CENTER:",1,3,16,17);
GUI_TOOLTIP(uspex_gui._centers,"centers: minimum allowed distance from the current molecule center to others.");
	GUI_APPLY_BUTTON_TABLE(table,uspex_gui._centers_button,apply_centers,3,4,16,17);
/* line 17 */
/* --- end page */

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
GUI_TOOLTIP(uspex_gui.mag_fmls,"magRatio: Ch. 4.1 DEFAULT: 0.225\n(magnetic calculation) Initial ratio of structures\nwith a low-spin ferromagnetic order.");
	GUI_ENTRY_TABLE(table,uspex_gui.mag_afml,uspex_gui.calc.magRatio[3],"%.4f","AFM-L:",2,3,2,3);
GUI_TOOLTIP(uspex_gui.mag_afml,"magRatio: Ch. 4.1 DEFAULT: 0.225\n(magnetic calculation) Initial ratio of structures\nwith a low-spin antiferromagnetic order.");
	GUI_ENTRY_TABLE(table,uspex_gui.mag_fmlh,uspex_gui.calc.magRatio[5],"%.4f","FM-LH:",3,4,2,3);
GUI_TOOLTIP(uspex_gui.mag_fmlh,"magRatio: Ch. 4.1 DEFAULT: 0.0\n(magnetic calculation) Initial ratio of structures\nwith a low/high spin mixed ferromagnetic order.");
/* line 3 */
	GUI_CHECK_TABLE(table,uspex_gui._calctype_mag_2,uspex_gui.calc._calctype_mag,mag_toggle,"MAG",0,1,3,4);
GUI_TOOLTIP(uspex_gui._calctype_mag_2,"Set a magnetic calculation.");
	GUI_ENTRY_TABLE(table,uspex_gui.mag_fmhs,uspex_gui.calc.magRatio[2],"%.4f","FM-HS:",1,2,3,4);
GUI_TOOLTIP(uspex_gui.mag_fmhs,"magRatio: Ch. 4.1 DEFAULT: 0.225\n(magnetic calculation) Initial ratio of structures\nwith a high-spin ferromagnetic.");
	GUI_ENTRY_TABLE(table,uspex_gui.mag_afmh,uspex_gui.calc.magRatio[4],"%.4f","AFM-H:",2,3,3,4);
GUI_TOOLTIP(uspex_gui.mag_afmh,"magRatio: Ch. 4.1 DEFAULT: 0.225\n(magnetic calculation) Initial ratio of structures\nwith a high-spin antiferromagnetic.");
	GUI_ENTRY_TABLE(table,uspex_gui.mag_aflh,uspex_gui.calc.magRatio[6],"%.4f","AF-LH:",3,4,3,4);
GUI_TOOLTIP(uspex_gui.mag_aflh,"magRatio: Ch. 4.1 DEFAULT: 0.0\n(magnetic calculation) Initial ratio of structures\nwith a low/high spin mixed antiferromagnetic.");
/* line 4 */
	GUI_ENTRY_TABLE(table,uspex_gui.bestFrac,uspex_gui.calc.bestFrac,"%.4f","Best:",0,1,4,5);
GUI_TOOLTIP(uspex_gui.bestFrac,"bestFrac: Ch. 4.3 DEFAULT: 0.7\nFraction of current generation used to generate the next.");
	GUI_ENTRY_TABLE(table,uspex_gui.keepBestHM,uspex_gui.calc.keepBestHM,"%3i","BestHM:",1,2,4,5);
GUI_TOOLTIP(uspex_gui.keepBestHM,"keepBestHM: Ch. 4.3 DEFAULT: auto\nNumber of best structures that will survive in next generation.");
	GUI_CHECK_TABLE(table,uspex_gui.reoptOld,uspex_gui.calc.reoptOld,NULL,"reopt",2,3,4,5);/*not calling anything*/
GUI_TOOLTIP(uspex_gui.reoptOld,"reoptOld: Ch. 4.3 DEFAULT: FALSE\nIf set surviving structures will be re-optimized.");
	GUI_ENTRY_TABLE(table,uspex_gui.fitLimit,uspex_gui.calc.fitLimit,"%.4f","fitLimit:",3,4,4,5);
GUI_TOOLTIP(uspex_gui.fitLimit,"fitLimit: Ch. 4.3 DEFAULT: none\nStop calculation when fitLimit fitness is reached.");
/* --- Structure & Variation */
	GUI_LABEL_TABLE(table,"Structure & Variation",0,4,5,6);
/* line 6 */
	GUI_TEXT_TABLE(table,uspex_gui.symmetries,uspex_gui.calc.symmetries,"symmetries: ",0,4,6,7);
GUI_TOOLTIP(uspex_gui.symmetries,"symmetries: Ch. 4.4 DEFAULT: auto\nPossible space groups for crystals,\nplane groups for 2D crystals/surfaces\nor point groups for clusters.");
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
GUI_TOOLTIP(uspex_gui.mutationDegree,"mutationDegree: Ch. 4.13 DEFAULT: auto\nMaximum displacement in softmutation (Ang).");
	GUI_ENTRY_TABLE(table,uspex_gui.mutationRate,uspex_gui.calc.mutationRate,"%.4f","mutationRate:",1,2,10,11);
GUI_TOOLTIP(uspex_gui.mutationRate,"mutationRate: Ch. 4.13 DEFAULT: 0.5\nStd. dev. of epsilon in strain matrix for lattice mutation.");
	GUI_ENTRY_TABLE(table,uspex_gui.DisplaceInLatmutation,uspex_gui.calc.DisplaceInLatmutation,"%.4f","D_LATMUT:",2,3,10,11);
GUI_TOOLTIP(uspex_gui.DisplaceInLatmutation,"DisplaceInLatmutation: Ch. ?.? DEFAULT: 1.0\nSets softmutation as part of lattice mutation\nand gives maximum displacement (Ang).\nThis keyword has disapear in the 10.1 version of USPEX.");
	GUI_CHECK_TABLE(table,uspex_gui.AutoFrac,uspex_gui.calc.AutoFrac,NULL,"AutoFrac",3,4,10,11);/*not calling anything*/
GUI_TOOLTIP(uspex_gui.AutoFrac,"AutoFrac: Ch. 4.4 DEFAULT: FALSE\nIf set variation parameters will be optimized during run.");
/* --- Fingerprint & Antiseeds */
	GUI_LABEL_TABLE(table,"Fingerprint, antiseed, & spacegroup",0,4,11,12);
/* line 12 */
//	GUI_LABEL_TABLE(table,"Fingerprint",0,1,12,13);
	GUI_ENTRY_TABLE(table,uspex_gui.RmaxFing,uspex_gui.calc.RmaxFing,"%.4f","RMax:",0,1,12,13);
GUI_TOOLTIP(uspex_gui.RmaxFing,"RmaxFing: Ch. 4.9 DEFAULT: 10.0\nFingerprint cutoff distance (Ang.).");
	GUI_ENTRY_TABLE(table,uspex_gui.deltaFing,uspex_gui.calc.deltaFing,"%.4f","delta:",1,2,12,13);
GUI_TOOLTIP(uspex_gui.deltaFing,"deltaFing: Ch. 4.9 DEFAULT: 0.08\nDiscretization of the fingerprint function.");
	GUI_ENTRY_TABLE(table,uspex_gui.sigmaFing,uspex_gui.calc.sigmaFing,"%.4f","sigma:",2,3,12,13);
GUI_TOOLTIP(uspex_gui.sigmaFing,"sigmaFing: Ch. 4.9 DEFAULT: 0.03\nGaussian broadening of interatomic distance.");
	GUI_ENTRY_TABLE(table,uspex_gui.toleranceFing,uspex_gui.calc.toleranceFing,"%.4f","TOL:",3,4,12,13);
GUI_TOOLTIP(uspex_gui.toleranceFing,"toleranceFing: Ch. 4.9 DEFAULT: 0.008\nMinimal cosine distance for\n2 different structures.");
/* line 13 */
	GUI_LABEL_TABLE(table,"Antiseed",0,1,13,14);
	GUI_ENTRY_TABLE(table,uspex_gui.antiSeedsActivation,uspex_gui.calc.antiSeedsActivation,"%4i","Activation:",1,2,13,14);
GUI_TOOLTIP(uspex_gui.antiSeedsActivation,"antiSeedsActivation: Ch. 4.10 DEFAULT: 5000\nGeneration at which antiseed is active.");
	GUI_ENTRY_TABLE(table,uspex_gui.antiSeedsMax,uspex_gui.calc.antiSeedsMax,"%.4f","Max:",2,3,13,14);
GUI_TOOLTIP(uspex_gui.antiSeedsMax,"antiSeedsMax: Ch. 4.10 DEFAULT: 0.000\nGaussian height in mean square deviation of the generation enthalpy\namong bestFrac structures, recommended value is 0.01.");
	GUI_ENTRY_TABLE(table,uspex_gui.antiSeedsSigma,uspex_gui.calc.antiSeedsSigma,"%.4f","sigma:",3,4,13,14);
GUI_TOOLTIP(uspex_gui.antiSeedsSigma,"antiSeedsSigma: Ch. 4.10 DEFAULT: 0.001\nGaussian width in average distances between generated structure\namong bestFrac structures, recommended value is 0.005.");
/* line 14 */
	GUI_LABEL_TABLE(table,"Space group",0,1,14,15);
	/*col 2: empty*/
	GUI_CHECK_TABLE(table,uspex_gui.doSpaceGroup,uspex_gui.calc.doSpaceGroup,SG_toggle,"Active",2,3,14,15);
GUI_TOOLTIP(uspex_gui.doSpaceGroup,"doSpaceGroup: Ch. 4.11 DEFAULT: TRUE\nActivate space group determination.");
	GUI_ENTRY_TABLE(table,uspex_gui.SymTolerance,uspex_gui.calc.SymTolerance,"%.4f","TOL:",3,4,14,15);
GUI_TOOLTIP(uspex_gui.SymTolerance,"SymTolerance: Ch. 4.11 DEFAULT: 0.10\nPrecision for symmetry determination.");
/* --- Variable-composition */
	GUI_LABEL_TABLE(table,"Variable-composition",0,4,15,16);
/* line 16 */
	GUI_ENTRY_TABLE(table,uspex_gui.firstGeneMax,uspex_gui.calc.firstGeneMax,"%4i","1st_Gen:",0,1,16,17);
GUI_TOOLTIP(uspex_gui.firstGeneMax,"firstGeneMax: Ch. 5.1 DEFAULT: 11\nNumber of composition sampled in 1st generation.");
	GUI_ENTRY_TABLE(table,uspex_gui.minAt,uspex_gui.calc.minAt,"%4i","min@:",1,2,16,17);
GUI_TOOLTIP(uspex_gui.minAt,"minAt: Ch. 5.1 DEFAULT: none\nMinimum number of atoms (molecules) in unitcell\nfor the first generation.");
	GUI_ENTRY_TABLE(table,uspex_gui.maxAt,uspex_gui.calc.maxAt,"%4i","max@:",2,3,16,17);
GUI_TOOLTIP(uspex_gui.maxAt,"maxAt: Ch. 5.1 DEFAULT: none\nMaximum number of atoms (molecules) in unitcell\nfor the first generation.");
/*col 4: empty*/
/* line 17 */
	GUI_ENTRY_TABLE(table,uspex_gui.fracTrans,uspex_gui.calc.fracTrans,"%.4f","fTrans:",0,1,17,18);
GUI_TOOLTIP(uspex_gui.fracTrans,"fracTrans: Ch. 5.1 DEFAULT: 0.1\nFraction of structures obtained by transmutation.");
	GUI_ENTRY_TABLE(table,uspex_gui.howManyTrans,uspex_gui.calc.howManyTrans,"%.4f","rTrans:",1,2,17,18);
GUI_TOOLTIP(uspex_gui.howManyTrans,"howManyTrans: Ch. 5.1 DEFAULT: 0.2\nMaximum ratio of transmutated atoms in a structure.");
	if(uspex_gui.calc._nspetrans<1){
		title=g_strdup("");
	}else{
		title=g_strdup_printf("%i",uspex_gui.calc.specificTrans[0]);
		for(idx=1;idx<uspex_gui.calc._nspetrans;idx++){
			tmp=g_strdup_printf("%s %i",title,uspex_gui.calc.specificTrans[idx]);
			g_free(title);
			title=tmp;
		}
	}
	GUI_TEXT_TABLE(table,uspex_gui.specificTrans,title,"Trans: ",2,4,17,18);
	g_free(title);
GUI_TOOLTIP(uspex_gui.specificTrans,"specificTrans: Ch. 5.1 DEFAULT: blank\nList of allowed transmutation.");
/* --- end page */

/*-----------------------*/
/* page 3 -> CALCULATION */
/*-----------------------*/
	GUI_PAGE_NOTE(notebook,page,"CALCULATION");
/* create a table in the page*/
        GUI_TABLE_NOTE(page,table,18,4);
/* --- Ab initio */
        GUI_LABEL_TABLE(table,"Ab initio",0,4,0,1);
/* line 1 - NEW */
	GUI_RADIO2_TABLE(table,
		uspex_gui.use_specific,"USE Specific Folder",toggle_use_specific,
		uspex_gui.set_specific,"SET Specific Folder",toggle_set_specific,0,2,1,2);
GUI_TOOLTIP(uspex_gui.use_specific,"Use a previously defined Specific folder.");
GUI_TOOLTIP(uspex_gui.set_specific,"Set parameter to create a Specific folder.");
	GUI_TEXT_TABLE(table,uspex_gui.spe_folder,uspex_gui._tmp_spe_folder,"SPE:",2,3,1,2);
GUI_TOOLTIP(uspex_gui.spe_folder,"Specific folder.");
	GUI_OPEN_BUTTON_TABLE(table,uspex_gui.spe_folder_button,load_spe_folder_dialog,3,4,1,2);
/* line 2 */
	GUI_SPIN_TABLE(table,uspex_gui._num_opt_steps,uspex_gui._tmp_num_opt_steps,spin_update_num_opt_steps,"N_STEPS",0,1,2,3);
	GUI_SPIN_RANGE(uspex_gui._num_opt_steps,1.,(gdouble)USPEX_MAX_NUM_OPT_STEPS);
GUI_TOOLTIP(uspex_gui._num_opt_steps,"Set the number of optimisation steps.\nIncluding fixed-cells and full relaxation steps.");
	GUI_SPIN_TABLE(table,uspex_gui._curr_step,uspex_gui._tmp_curr_step,spin_update_curr_step,"STEPS #",1,2,2,3);
	GUI_SPIN_RANGE(uspex_gui._curr_step,1.,uspex_gui._tmp_num_opt_steps);
GUI_TOOLTIP(uspex_gui._curr_step,"Select current optimisation step number.");
	GUI_CHECK_TABLE(table,button,uspex_gui.auto_step,toggle_auto_step,"AUTO STEP",2,3,2,3);
GUI_TOOLTIP(button,"Let GDIS generate some \"Simplified\"\nINPUTS for this step, given\nthe calculation code.\nWARNING: GDIS use some general settings that may fail.");
	GUI_CHECK_TABLE(table,uspex_gui._isfixed,uspex_gui._tmp_isfixed,NULL,"Relax_fixed",3,4,2,3);/*not calling anything*/
GUI_TOOLTIP(uspex_gui._isfixed,"Set the current step as a fixed relaxation.\nUseful only with a mix between fixed and full relaxations.");
/* line 3 */
	GUI_COMBOBOX_TABLE(table,uspex_gui.abinitioCode,"CODE:",0,2,3,4);
	GUI_COMBOBOX_ADD(uspex_gui.abinitioCode,"0  - NONE");
	GUI_COMBOBOX_ADD(uspex_gui.abinitioCode,"1  - VASP");
	GUI_COMBOBOX_ADD(uspex_gui.abinitioCode,"2  - SIESTA");
	GUI_COMBOBOX_ADD(uspex_gui.abinitioCode,"3  - GULP");
	GUI_COMBOBOX_ADD(uspex_gui.abinitioCode,"4  - LAMMPS");
	GUI_COMBOBOX_ADD(uspex_gui.abinitioCode,"5  - ORCA");
	GUI_COMBOBOX_ADD(uspex_gui.abinitioCode,"6  - DMACRYS");
	GUI_COMBOBOX_ADD(uspex_gui.abinitioCode,"7  - CP2K");
	GUI_COMBOBOX_ADD(uspex_gui.abinitioCode,"8  - QE");
	GUI_COMBOBOX_ADD(uspex_gui.abinitioCode,"9  - FHI-aims");
	GUI_COMBOBOX_ADD(uspex_gui.abinitioCode,"10 - ATK");
	GUI_COMBOBOX_ADD(uspex_gui.abinitioCode,"11 - CASTEP");
	GUI_COMBOBOX_ADD(uspex_gui.abinitioCode,"12 - Tinker");
	GUI_COMBOBOX_ADD(uspex_gui.abinitioCode,"13 - MOPAC");
	GUI_COMBOBOX_ADD(uspex_gui.abinitioCode,"14 - BoltzTraP");
	GUI_COMBOBOX_ADD(uspex_gui.abinitioCode,"15 - DFTB");
	GUI_COMBOBOX_ADD(uspex_gui.abinitioCode,"16 - Gaussian");
	GUI_COMBOBOX_ADD(uspex_gui.abinitioCode,"17 - N/A");
	GUI_COMBOBOX_ADD(uspex_gui.abinitioCode,"18 - Abinit");
	GUI_COMBOBOX_ADD(uspex_gui.abinitioCode,"19 - CRYSTAL");
GUI_TOOLTIP(uspex_gui.abinitioCode,"abinitioCode: Ch. 4.8 DEFAULT: 1\nCode used for the calculations in this step.");
	GUI_ENTRY_TABLE(table,uspex_gui.KresolStart,uspex_gui.calc.KresolStart[0],"%.4f","Kresol:",2,3,3,4);
GUI_TOOLTIP(uspex_gui.KresolStart,"KresolStart: Ch. 4.8 DEFAULT: 0.2-0.08\nReciprocal-space resolution for this optimization step (2*pi/Ang).");
	GUI_ENTRY_TABLE(table,uspex_gui.vacuumSize,uspex_gui.calc.vacuumSize[0],"%.4f","vacuum:",3,4,3,4);
GUI_TOOLTIP(uspex_gui.vacuumSize,"vacuumSize: Ch. 4.8 DEFAULT: 10.0\nVacuum size (Ang.) between neighbor atoms from adjacent unit cells.\nFor cluster, 2D-crystal, and surfaces only.");
/* line 4 */
	GUI_TEXT_TABLE(table,uspex_gui.ai_input,uspex_gui._tmp_ai_input[0],"INP:",0,1,4,5);
GUI_TOOLTIP(uspex_gui.ai_input,"Input for the current calculation step.\nFile name as to end with *_N where N is the step number.");
	GUI_OPEN_BUTTON_TABLE(table,uspex_gui.ai_input_button,load_ai_input_dialog,1,2,4,5);
	GUI_TEXT_TABLE(table,uspex_gui.ai_opt,uspex_gui._tmp_ai_opt[0],"OPT:",2,3,4,5);
GUI_TOOLTIP(uspex_gui.ai_opt,"Complementary (optional) file for the current calculation step.\nWhen present, file name as to end with *_N where N is the step number.");
	GUI_OPEN_BUTTON_TABLE(table,uspex_gui.ai_opt_button,load_ai_opt_dialog,3,4,4,5);
/* line 5 */
	GUI_TEXT_TABLE(table,uspex_gui.ai_lib,uspex_gui._tmp_ai_lib_folder[0],"LIB:",0,1,5,6);
GUI_TOOLTIP(uspex_gui.ai_lib,"Set the library directory.\nCalculations often require pseudo-potentials for each species\nor a potential library file wich describe each species.\nEach different calculation step code require a library directory.");
	GUI_OPEN_BUTTON_TABLE(table,uspex_gui.ai_lib_button,load_ai_lib_dialog,1,2,5,6);
	GUI_COMBOBOX_TABLE(table,uspex_gui.ai_lib_flavor,"FLAVOR:",2,3,5,6);
GUI_TOOLTIP(uspex_gui.ai_lib_flavor,"Select flavor(s) that apply to current potential.\nIt is either a per-species (VASP)\nor per-potential (GULP) switch.");
	GUI_COMBOBOX_EDITABLE(uspex_gui.ai_lib_flavor);/*NEW: accept user-made library*/
	GUI_TEXT_TABLE(table,uspex_gui.ai_lib_sel,uspex_gui._tmp_ai_lib_sel[0],"lib:",3,4,5,6);
GUI_TOOLTIP(uspex_gui.ai_lib_sel,"The flavor of selected libraries.");
	GUI_LOCK(uspex_gui.ai_lib_sel);
/* line 6 */
	GUI_TEXT_TABLE(table,uspex_gui.commandExecutable,uspex_gui._tmp_commandExecutable[0],"EXE:",0,1,6,7);
GUI_TOOLTIP(uspex_gui.commandExecutable,"commandExecutable: Ch. 4.8 DEFAULT: none\nCurrent optimisation step executable of submission script.");
	GUI_OPEN_BUTTON_TABLE(table,button,load_abinitio_exe_dialog,1,2,6,7);
	/*2,3 blank*/
	GUI_APPLY_BUTTON_TABLE(table,button,apply_step,3,4,6,7);
GUI_TOOLTIP(button,"Apply changes to this step (only).");
/* --- USPEX launch */
	GUI_LABEL_TABLE(table,"USPEX launch",0,4,7,8);
/* line 8 */
	GUI_TEXT_TABLE(table,uspex_gui.job_uspex_exe,uspex_gui.calc.job_uspex_exe,"USPEX:",0,1,8,9);
GUI_TOOLTIP(uspex_gui.job_uspex_exe,"The USPEX script executable.");
	GUI_OPEN_BUTTON_TABLE(table,button,load_uspex_exe,1,2,8,9);
	GUI_ENTRY_TABLE(table,uspex_gui.whichCluster,uspex_gui.calc.whichCluster,"%i","Cluster:",2,3,8,9);
GUI_TOOLTIP(uspex_gui.whichCluster,"whichCluster: Ch. 4.8 DEFAULT: 0\nType of job submission, including:\n0 - no job-script;\n1 - local submission;\n2 - remote submission.\n>2 - specific supercomputer.");
	GUI_CHECK_TABLE(table,uspex_gui.PhaseDiagram,uspex_gui.calc.PhaseDiagram,NULL,"PhaseDiag",3,4,8,9);/*not calling anything*/
GUI_TOOLTIP(uspex_gui.PhaseDiagram,"PhaseDiagram: Ch. 4.8 DEFAULT: FALSE\nSet calculation of an estimated phase diagram\nfor calculationMethod 30X (300, 301, s300, and s301).");
/* line 9 */
	GUI_CHECK_TABLE(table,uspex_gui.sel_v1030_3,uspex_gui.have_v1030,toggle_v1030,"Ver 10.4",0,1,9,10);
GUI_TOOLTIP(uspex_gui.sel_v1030_3,"Use USPEX v. 10.4 instead of 9.4.4.");
        if(uspex_gui.calc.numProcessors==NULL){
                title=g_strdup("");
        }else{
                title=g_strdup_printf("%i",uspex_gui.calc.numProcessors[0]);
                for(idx=1;idx<uspex_gui.calc._num_opt_steps;idx++){
                        tmp=g_strdup_printf("%s %i",title,uspex_gui.calc.numProcessors[idx]);
                        g_free(title);
                        title=tmp;
                }
        }
	GUI_TEXT_TABLE(table,uspex_gui.numProcessors,title,"CPU:",1,2,9,10);
        g_free(title);
GUI_TOOLTIP(uspex_gui.numProcessors,"numProcessors: Ch. ? DEFAULT: 1\nNumber of processors used in each steps.");
	GUI_ENTRY_TABLE(table,uspex_gui.numParallelCalcs,uspex_gui.calc.numParallelCalcs,"%i","PAR:",2,3,9,10);
GUI_TOOLTIP(uspex_gui.numParallelCalcs,"numParallelCalcs: Ch. 4.8 DEFAULT: 1\nNumber of structure relaxations in parallel.");
	GUI_CHECK_TABLE(table,uspex_gui.sel_octave,uspex_gui.have_octave,NULL,"OCTAVE",3,4,9,10);/*not calling anything*/
GUI_TOOLTIP(uspex_gui.sel_octave,"Use octave instead of matlab.");
/* line 10 */
	GUI_TEXT_TABLE(table,uspex_gui.job_path,uspex_gui.calc.job_path,"Folder:",0,1,10,11);
GUI_TOOLTIP(uspex_gui.job_path,"Select USPEX calculation folder.\nThis folder is where run_uspex will be launched\nand result will be read by GDIS.");
	GUI_OPEN_BUTTON_TABLE(table,button,uspex_path_dialog,1,2,10,11);
	GUI_TEXT_TABLE(table,uspex_gui.remoteFolder,uspex_gui.calc.remoteFolder,"Remote:",2,3,10,11);
GUI_TOOLTIP(uspex_gui.remoteFolder,"remoteFolder: Ch. 4.8 DEFAULT: none\nWhen using remote submission, this hold the\nexecution folder on the distant computer.");
	GUI_OPEN_BUTTON_TABLE(table,button,load_remote_folder,3,4,10,11);
	/*TODO: disable remote when whichCluster<2*/
/* --- Restart */
	GUI_LABEL_TABLE(table,"Restart",0,4,11,12);
/* line 12 */
	GUI_CHECK_TABLE(table,uspex_gui.pickUpYN,uspex_gui.calc.pickUpYN,NULL,"RESTART",0,1,12,13);/*not calling anything*/
GUI_TOOLTIP(uspex_gui.pickUpYN,"pickUpYN: deprecated DEFAULT: 0\nSet to restart a calculation, deprecated on version <= 10.");
	GUI_ENTRY_TABLE(table,uspex_gui.pickUpGen,uspex_gui.calc.pickUpGen,"%i","GEN:",1,2,12,13);
GUI_TOOLTIP(uspex_gui.pickUpGen,"pickUpGen: Ch. 4.7 DEFAULT: 0\nSelect the generation at which USPEX should restart.\nA new calculation is started for 0 (default).");
	GUI_ENTRY_TABLE(table,uspex_gui.pickUpFolder,uspex_gui.calc.pickUpFolder,"%i","Folder:",2,3,12,13);
GUI_TOOLTIP(uspex_gui.pickUpFolder,"pickUpFolder: Ch. 4.7 DEFAULT: 0\nSelect the folder number from which to restart from.\nSetting N means restarting from a folder named \"resultN\".");
	GUI_CHECK_TABLE(table,button,uspex_gui.restart_cleanup,NULL,"CLEANUP",3,4,12,13);/*not calling anything*/
GUI_TOOLTIP(button,"Cleanup the calculation directory before restart,\nie. remove still_running, NOT_YET, etc.");
/* reserved for future use */
	for(idx=13;idx<18;idx++) GUI_LABEL_TABLE(table," ",0,4,idx,idx+1);
/* --- end page */

/*--------------------*/
/* page 4 -> ADVANCED */
/*--------------------*/
	GUI_PAGE_NOTE(notebook,page,"ADVANCED");
/* create a table in the page*/
	GUI_TABLE_NOTE(page,table,18,4);
/* --- Developers */
	GUI_LABEL_TABLE(table,"Developers",0,4,0,1);
/* line 1 */
	GUI_ENTRY_TABLE(table,uspex_gui.repeatForStatistics,uspex_gui.calc.repeatForStatistics,"%i","REPEAT:",0,1,1,2);
GUI_TOOLTIP(uspex_gui.repeatForStatistics,"repeatForStatistics: Ch. 4.12 DEFAULT: 1\nNumber of USPEX repeated job execution.\nFor default 1 value, no statistic is collected.");
	GUI_ENTRY_TABLE(table,uspex_gui.stopFitness,uspex_gui.calc.stopFitness,"%.4f","STOP_FIT:",1,2,1,2);
GUI_TOOLTIP(uspex_gui.stopFitness,"stopFitness: Ch. 4.12 DEFAULT: none\nThe fitness value at which USPEX calculation will stop.\nSee also fitLimit in STRUCTURE page.");
	GUI_ENTRY_TABLE(table,uspex_gui.fixRndSeed,uspex_gui.calc.fixRndSeed,"%i","RND_SEED:",2,3,1,2);
GUI_TOOLTIP(uspex_gui.fixRndSeed,"fixRndSeed: Ch. 4.12 DEFAULT: 0\nFor non-zero values, fix the random seed\nto check that for two USPEX calculations\nwith the same seed produce the same results.");
	GUI_CHECK_TABLE(table,uspex_gui.collectForces,uspex_gui.calc.collectForces,NULL,"CollectForces",3,4,1,2);/*not calling anything*/
GUI_TOOLTIP(uspex_gui.collectForces,"collectForces: Ch. 4.12 DEFAULT: FALSE\nCollect relaxation information from VASP calculation\n(energies, forces, positions, geometry, and stress).");
/* --- Seldom */
	GUI_LABEL_TABLE(table,"Seldom",0,4,2,3);
/* line 3 */
	GUI_CHECK_TABLE(table,uspex_gui.ordering_active,uspex_gui.calc.ordering_active,NULL,"ordering",0,1,3,4);/*not calling anything*/
GUI_TOOLTIP(uspex_gui.ordering_active,"ordering_active: Ch. 4.12 DEFAULT: TRUE\nSwitch the \"biasing of variation operators by local order parameters\".");
	GUI_CHECK_TABLE(table,uspex_gui.symmetrize,uspex_gui.calc.symmetrize,NULL,"symmetrize",1,2,3,4);/*not calling anything*/
GUI_TOOLTIP(uspex_gui.symmetrize,"symmetrize: Ch. 4.13 DEFAULT: FALSE\nTransform all structure to \"standard symmetry-adapted crystallografic setting\".");
        if(uspex_gui.calc.valenceElectr==NULL){
                title=g_strdup("");
        }else{
                title=g_strdup_printf("%i",uspex_gui.calc.valenceElectr[0]);
		for(idx=1;idx<uspex_gui.calc._nspecies;idx++){
			tmp=g_strdup_printf("%s %i",title,uspex_gui.calc.valenceElectr[idx]);
			g_free(title);
			title=tmp;
		}
        }
	GUI_TEXT_TABLE(table,uspex_gui.valenceElectr,title,"VALENCE:",2,4,3,4);
	g_free(title);
GUI_TOOLTIP(uspex_gui.valenceElectr,"valenceElectr: Ch. 4.13 DEFAULT: AUTO\nNumber of valence electrons.\nOverrides the tabulated values.");
/* line 4 */
	GUI_ENTRY_TABLE(table,uspex_gui.percSliceShift,uspex_gui.calc.percSliceShift,"%.4f","SliceShift:",0,1,4,5);
GUI_TOOLTIP(uspex_gui.percSliceShift,"percSliceShift: Ch. 4.13 DEFAULT: 1.0\nProbability of shifting slabs.");
	GUI_ENTRY_TABLE(table,uspex_gui.minSlice,uspex_gui.calc.minSlice,"%.4f","minSlice:",1,2,4,5);
GUI_TOOLTIP(uspex_gui.minSlice,"minSlice: Ch. 4.13 DEFAULT: N/A\nMinimal Slice thickness (Ang.).\nRecommended value is ~1 Ang.");
	GUI_COMBOBOX_TABLE(table,uspex_gui.dynamicalBestHM,"DYN_HM:",2,4,4,5);
	GUI_COMBOBOX_ADD(uspex_gui.dynamicalBestHM,"0 - NONE");
	GUI_COMBOBOX_ADD(uspex_gui.dynamicalBestHM,"1 - lowest Energy");
	GUI_COMBOBOX_ADD(uspex_gui.dynamicalBestHM,"2 - promote diversity");
GUI_TOOLTIP(uspex_gui.dynamicalBestHM,"dynamicalBestHM: Ch. ?.? DEFAULT: 2\nSpecify if and how the number of surviving structure will vary.\nThis keyword has disapear in the 10.1 version of USPEX.");
/* line 5 */
	/*col 1: empty*/
	GUI_ENTRY_TABLE(table,uspex_gui.maxSlice,uspex_gui.calc.maxSlice,"%.4f","maxSlice:",1,2,5,6);
GUI_TOOLTIP(uspex_gui.maxSlice,"maxSlice: Ch. 4.13 DEFAULT: N/A\nMaximum Slice thickness (Ang.).\nRecommended value is ~6 Ang.");
	GUI_TEXT_TABLE(table,uspex_gui.softMutOnly,uspex_gui.calc.softMutOnly,"SoftMut:",2,4,5,6);
GUI_TOOLTIP(uspex_gui.softMutOnly,"softMutOnly: Ch. ?.? DEFAULT: 0\nWhich/how many generations are produced from soft mutation only.\nThis keyword has disapear in the 10.1 version of USPEX.");
/* line 6 */
	GUI_ENTRY_TABLE(table,uspex_gui.maxDistHeredity,uspex_gui.calc.maxDistHeredity,"%.4f","DistHer:",0,1,6,7);
GUI_TOOLTIP(uspex_gui.maxDistHeredity,"maxDistHeredity: Ch. 4.13 DEFAULT: 0.5\nMax fingerprint distance between structures chosen for heredity.");
	GUI_ENTRY_TABLE(table,uspex_gui.numberparents,uspex_gui.calc.numberparents,"%i","NumP:",1,2,6,7);
GUI_TOOLTIP(uspex_gui.numberparents,"numberparents: Ch. 4.13 DEFAULT: 2\nNumber of parents chosen for heredity\nfor cluster calculation.");
	GUI_COMBOBOX_TABLE(table,uspex_gui.manyParents,"many_P:",2,4,6,7);
	GUI_COMBOBOX_ADD(uspex_gui.manyParents,"0 - 2 parents,    1 slice each");
	GUI_COMBOBOX_ADD(uspex_gui.manyParents,"1 - n structures, 1 slice each");
	GUI_COMBOBOX_ADD(uspex_gui.manyParents,"2 - 2 structures, n slices, independants");
	GUI_COMBOBOX_ADD(uspex_gui.manyParents,"3 - 2 structures, n slices, fixed offset");
GUI_TOOLTIP(uspex_gui.manyParents,"manyParents: Ch. 4.13 DEFAULT: 0\nDetermine if and how more slices\nand parents structures are chosen.\nFor large systems recommended\nsetting 3 may be beneficial.");
/* --- BoltzTraP */
	GUI_LABEL_TABLE(table,"BoltzTraP",0,4,7,8);
/* line 8 */
	GUI_COMBOBOX_TABLE(table,uspex_gui.TE_goal,"Goal:",0,1,8,9);
	GUI_COMBOBOX_ADD(uspex_gui.TE_goal,"ZT tensor trace");
	GUI_COMBOBOX_ADD(uspex_gui.TE_goal,"ZT_xx x component");
	GUI_COMBOBOX_ADD(uspex_gui.TE_goal,"ZT_yy y component");
	GUI_COMBOBOX_ADD(uspex_gui.TE_goal,"ZT_zz z component");
GUI_TOOLTIP(uspex_gui.TE_goal,"TE_goal: Ch. 5.2 DEFAULT: ZT\nSet the component of ZT to be optimized.");
	GUI_ENTRY_TABLE(table,uspex_gui.BoltzTraP_T_max,uspex_gui.calc.BoltzTraP_T_max,"%.4f","T_Max:",1,2,8,9);
GUI_TOOLTIP(uspex_gui.BoltzTraP_T_max,"BoltzTraP_T_max: Ch. 5.2 DEFAULT: 800.0\nMaximum BoltzTraP calculation temperature.");
	GUI_ENTRY_TABLE(table,uspex_gui.BoltzTraP_T_delta,uspex_gui.calc.BoltzTraP_T_delta,"%.4f","T_delta:",2,3,8,9);
GUI_TOOLTIP(uspex_gui.BoltzTraP_T_delta,"BoltzTraP_T_delta: Ch. 5.2 DEFAULT: 50.0\nBoltzTraP calculation temperature increment.");
	GUI_ENTRY_TABLE(table,uspex_gui.BoltzTraP_T_efcut,uspex_gui.calc.BoltzTraP_T_efcut,"%.4f","T_efcut:",3,4,8,9);
GUI_TOOLTIP(uspex_gui.BoltzTraP_T_efcut,"BoltzTraP_T_efcut: Ch. 5.2 DEFAULT: 0.15\nBoltzTraP calculation chemical potential (eV) interval.");
/* line 9 */
	GUI_TEXT_TABLE(table,uspex_gui.cmd_BoltzTraP,uspex_gui._tmp_cmd_BoltzTraP,"cmd:",0,1,9,10);
GUI_TOOLTIP(uspex_gui.cmd_BoltzTraP,"BoltzTraP software command/script (optional).");
GUI_OPEN_BUTTON_TABLE(table,uspex_gui.cmd_BoltzTraP_button,load_cmd_BoltzTraP,1,2,9,10);
	GUI_ENTRY_TABLE(table,uspex_gui.TE_T_interest,uspex_gui.calc.TE_T_interest,"%.4f","T_target:",2,3,9,10);
GUI_TOOLTIP(uspex_gui.TE_T_interest,"TE_T_interest: Ch. 5.2 DEFAULT: 300.0Y\nTarget temperature to optimize thermoelectric efficiency.");
	GUI_ENTRY_TABLE(table,uspex_gui.TE_threshold,uspex_gui.calc.TE_threshold,"%.4f","Threshold:",3,4,9,10);
GUI_TOOLTIP(uspex_gui.TE_threshold,"TE_threshold: Ch. 5.2 DEFAULT: 0.5\nStructures with a ZT figure of merit\nbelow this threshold will be discarded.");
/* --- Transition Path Sampling */
	GUI_LABEL_TABLE(table,"Transition Path Sampling",0,4,10,11);
/* line 11 */
	GUI_ENTRY_TABLE(table,uspex_gui.numIterations,uspex_gui.calc.numIterations,"%i","N_Iter:",0,1,11,12);
GUI_TOOLTIP(uspex_gui.numIterations,"numIterations: Ch. 6.3 DEFAULT: 1000\nMaximum number of TPS iterations.");
	/*col 2: empty*/
	GUI_ENTRY_TABLE(table,uspex_gui.shiftRatio,uspex_gui.calc.shiftRatio,"%.4f","rShift:",2,3,11,12);
GUI_TOOLTIP(uspex_gui.shiftRatio,"shiftRatio: Ch. 6.3 DEFAULT: 0.1\nFraction of shooter-after-shifter operations.");
	GUI_CHECK_TABLE(table,uspex_gui.orderParaType,uspex_gui.calc.orderParaType,NULL,"OP_TYPE",3,4,11,12);/*not calling anything*/
GUI_TOOLTIP(uspex_gui.orderParaType,"orderParaType: Ch. 6.3 DEFAULT: none\nSelect if the method of order parameter calculation is:\nthe fingerprint method (TRUE)\na user-defined method (FALSE).\nContrary to USPEX lack of default, TRUE is pre-selected.");
/* line 12 */
	GUI_TEXT_TABLE(table,uspex_gui.speciesSymbol,uspex_gui.calc.speciesSymbol,"SpeciesSymbols:",0,2,12,13);
GUI_TOOLTIP(uspex_gui.speciesSymbol,"speciesSymbol: Ch. 6.3 DEFAULT: none\nIdentity of all chemical species (atoms/molecules).");
	if(uspex_gui.calc.mass==NULL){
		title=g_strdup("");
	}else{
		title=g_strdup_printf("%lf",uspex_gui.calc.mass[0]);
		for(idx=1;idx<uspex_gui.calc._nspecies;idx++){
			tmp=g_strdup_printf("%s %lf",title,uspex_gui.calc.mass[idx]);
			g_free(title);
			title=tmp;
		}
	}
	GUI_TEXT_TABLE(table,uspex_gui.mass,title,"mass:",2,4,12,13);
	g_free(title);
GUI_TOOLTIP(uspex_gui.mass,"mass: Ch. 6.3 DEFAULT: Auto\nMass of each corresponding species.");
/* line 13 */
	GUI_ENTRY_TABLE(table,uspex_gui.amplitudeShoot_AB,uspex_gui.calc.amplitudeShoot[0],"%.4f","A(A->B):",0,1,13,14);
GUI_TOOLTIP(uspex_gui.amplitudeShoot_AB,"amplitudeShoot: Ch. 6.3 DEFAULT: 0.1\nMomentum amplitude for A->B shooting.");
	GUI_ENTRY_TABLE(table,uspex_gui.amplitudeShoot_BA,uspex_gui.calc.amplitudeShoot[1],"%.4f","A(B->A):",1,2,13,14);
GUI_TOOLTIP(uspex_gui.amplitudeShoot_BA,"amplitudeShoot: Ch. 6.3 DEFAULT: 0.1\nMomentum amplitude for B->A shooting.");
	GUI_ENTRY_TABLE(table,uspex_gui.magnitudeShoot_success,uspex_gui.calc.magnitudeShoot[0],"%.4f","M(success):",2,3,13,14);
GUI_TOOLTIP(uspex_gui.magnitudeShoot_success,"magnitudeShoot: Ch. 6.3 DEFAULT: 1.05\nAmplitude increasing factor on MD trajectory success.");
	GUI_ENTRY_TABLE(table,uspex_gui.magnitudeShoot_failure,uspex_gui.calc.magnitudeShoot[1],"%.4f","M(failure):",3,4,13,14);
GUI_TOOLTIP(uspex_gui.magnitudeShoot_failure,"magnitudeShoot: Ch. 6.3 DEFAULT: 1.05\nAmplitude decreasing factor on MD trajectory failure.");
/* line 14 */
	GUI_TEXT_TABLE(table,uspex_gui.cmdOrderParameter,uspex_gui.calc.cmdOrderParameter,"cmdOP:",0,2,14,15);
GUI_TOOLTIP(uspex_gui.cmdOrderParameter,"cmdOrderParameter: Ch. 6.3 DEFAULT: none\nUser-defined command for order parameter calculation.");
	GUI_OPEN_BUTTON_TABLE(table,uspex_gui.cmdOrderParameter_button,load_cmd_OP,2,3,14,15);
	GUI_ENTRY_TABLE(table,uspex_gui.opCriteria_start,uspex_gui.calc.opCriteria[0],"%.4f","SIM(start):",3,4,14,15);
GUI_TOOLTIP(uspex_gui.opCriteria_start,"opCriteria: Ch. 6.3 DEFAULT: none\nAllowable degree of similarity between starting states.");
/* line 15 */
	GUI_TEXT_TABLE(table,uspex_gui.cmdEnthalpyTemperature,uspex_gui.calc.cmdEnthalpyTemperature,"cmdET:",0,2,15,16);
GUI_TOOLTIP(uspex_gui.cmdEnthalpyTemperature,"cmdEnthalpyTemperature: Ch. 6.3 DEFAULT: none\nUser-defined command for enthalpy/temperature\nextraction from the MD results.");
	GUI_OPEN_BUTTON_TABLE(table,uspex_gui.cmdEnthalpyTemperature_button,load_cmd_ET,2,3,15,16);
	GUI_ENTRY_TABLE(table,uspex_gui.opCriteria_end,uspex_gui.calc.opCriteria[1],"%.4f","SIM(end):",3,4,15,16);
GUI_TOOLTIP(uspex_gui.opCriteria_end,"opCriteria: Ch. 6.3 DEFAULT: none\nAllowable degree of similarity between ending states.");
/* line 16 */
	GUI_TEXT_TABLE(table,uspex_gui.orderParameterFile,uspex_gui.calc.orderParameterFile,"OP_file:",0,1,16,17);
GUI_TOOLTIP(uspex_gui.orderParameterFile,"orderParameterFile: Ch. 6.3 DEFAULT: fp.dat\nOrder parameter history file.");
	GUI_OPEN_BUTTON_TABLE(table,uspex_gui.orderParameterFile_button,load_OP_file,1,2,16,17);
	GUI_TEXT_TABLE(table,uspex_gui.enthalpyTemperatureFile,uspex_gui.calc.enthalpyTemperatureFile,"ET_file:",2,3,16,17);
GUI_TOOLTIP(uspex_gui.enthalpyTemperatureFile,"enthalpyTemperatureFile: Ch. 6.3 DEFAULT: HT.dat\nEnthalpy and temperature history file.");
	GUI_OPEN_BUTTON_TABLE(table,uspex_gui.enthalpyTemperatureFile_button,load_ET_file,3,4,16,17);
/* line 17 */
	GUI_TEXT_TABLE(table,uspex_gui.trajectoryFile,uspex_gui.calc.trajectoryFile,"traj_file:",0,1,17,18);
GUI_TOOLTIP(uspex_gui.trajectoryFile,"trajectoryFile: Ch. 6.3 DEFAULT: traj.dat\nMD trajectory file.");
	GUI_OPEN_BUTTON_TABLE(table,uspex_gui.trajectoryFile_button,load_traj_file,1,2,17,18);
	GUI_TEXT_TABLE(table,uspex_gui.MDrestartFile,uspex_gui.calc.MDrestartFile,"MD_file:",2,3,17,18);
GUI_TOOLTIP(uspex_gui.MDrestartFile,"MDrestartFile: Ch. 6.3 DEFAULT: traj.restart\nMD restart file.");
	GUI_OPEN_BUTTON_TABLE(table,uspex_gui.MDrestartFile_button,load_MDrestart_file,3,4,17,18);
/* --- end page */

/*--------------------*/
/* page 5 -> SPECIFIC */
/*--------------------*/
	GUI_PAGE_NOTE(notebook,uspex_gui.specific_page,"SPECIFIC");
/* create a table in the page*/
	GUI_TABLE_NOTE(uspex_gui.specific_page,table,18,4);
/* --- Metadynamics */
	GUI_LABEL_TABLE(table,"Metadynamics",0,4,0,1);
/* line 1 */
	GUI_COMBOBOX_TABLE(table,uspex_gui.FullRelax,"Relax:",0,1,1,2);
	GUI_COMBOBOX_ADD(uspex_gui.FullRelax,"0 - No full relaxation (fix cells)");
	GUI_COMBOBOX_ADD(uspex_gui.FullRelax,"1 - Relax only the best structures");
	GUI_COMBOBOX_ADD(uspex_gui.FullRelax,"2 - Relax all different structures");
GUI_TOOLTIP(uspex_gui.FullRelax,"FullRelax: Ch. 5.7 DEFAULT: 2\nPerform full relaxation of which structures for analysis.\nRecommended value is 2.");
	GUI_ENTRY_TABLE(table,uspex_gui.maxVectorLength,uspex_gui.calc.maxVectorLength,"%.4f","MaxV:",1,2,1,2);
GUI_TOOLTIP(uspex_gui.maxVectorLength,"maxVectorLength: Ch. 5.7 DEFAULT: none\nAdd a correction force to keep cell length below this setting.");
	GUI_ENTRY_TABLE(table,uspex_gui.GaussianWidth,uspex_gui.calc.GaussianWidth,"%.4f","GaussW:",2,3,1,2);
GUI_TOOLTIP(uspex_gui.GaussianWidth,"GaussianWidth: Ch. 5.7 DEFAULT: AUTO\nWidth of Gaussian added to PES to accelerate phase transition.\nRecommended values is 0.10~0.15L (L = minimum cell length).");
	GUI_ENTRY_TABLE(table,uspex_gui.GaussianHeight,uspex_gui.calc.GaussianHeight,"%.4f","GaussH:",3,4,1,2);
GUI_TOOLTIP(uspex_gui.GaussianHeight,"GaussianHeight: Ch. 5.7 DEFAULT: AUTO\nHeight of Gaussian added to PES to accelerate phase transition.\nRecommended values is L.dh^2.G with L = average cell length,\ndh = GaussW, and G = shear modulus.");
/* line 2 */
	/*col 1: empty*/
//	GUI_ENTRY_TABLE(table,uspex_gui.ExternalPressure,uspex_gui.calc.ExternalPressure,"%.4f","ExtP:",0,1,2,3);
//GUI_TOOLTIP(uspex_gui.ExternalPressure,"ExternalPressure: Ch. 4.1, 5.7 DEFAULT: none\nExternal pressure (GPa) for calculation.");
	GUI_COMBOBOX_TABLE(table,uspex_gui.meta_model,"MODEL:",1,3,2,3);
	GUI_COMBOBOX_ADD(uspex_gui.meta_model,"From POSCAR FILE");
	GUI_COMBOBOX_ADD(uspex_gui.meta_model,"UNDER CONSTRUCTION");
GUI_TOOLTIP(uspex_gui.meta_model,"Select the model from which metadynamics is started.\nA good structure, relaxed at ExtP is necessary.");
	GUI_OPEN_BUTTON_TABLE(table,uspex_gui.meta_model_button,load_meta_start_file,3,4,2,3);
/* --- Particles Swarm optimization */
/* line 3 */
	/*col 1: empty*/
	GUI_LABEL_TABLE(table,"PSO:",0,1,3,4);
	GUI_ENTRY_TABLE(table,uspex_gui.PSO_softMut,uspex_gui.calc.PSO_softMut,"%.4f","SoftMut:",1,2,3,4);
GUI_TOOLTIP(uspex_gui.PSO_softMut,"PSO_softMut: Ch. 5.8 DEFAULT: 1\nSoft mutation weight.");
	GUI_ENTRY_TABLE(table,uspex_gui.PSO_BestStruc,uspex_gui.calc.PSO_BestStruc,"%.4f","BestStruct:",2,3,3,4);
GUI_TOOLTIP(uspex_gui.PSO_BestStruc,"PSO_BestStruc: Ch. 5.8 DEFAULT: 1\nWeight of heredity with best position of a given PSO particle.");
	GUI_ENTRY_TABLE(table,uspex_gui.PSO_BestEver,uspex_gui.calc.PSO_BestEver,"%.4f","BestEver:",3,4,3,4);
GUI_TOOLTIP(uspex_gui.PSO_BestEver,"PSO_BestEver: Ch. 5.8 DEFAULT: 1\nWeight of heredity with globally best PSO particle.");
/* --- Variable-cell nudged elastic band */
	GUI_LABEL_TABLE(table,"Variable-cell nudged elastic band",0,4,4,5);
/* line 5 */
	GUI_COMBOBOX_TABLE(table,uspex_gui._vcnebtype_method,"Method:",0,1,5,6);
	GUI_COMBOBOX_ADD(uspex_gui._vcnebtype_method,"1 - VC-NEB method");
	GUI_COMBOBOX_ADD(uspex_gui._vcnebtype_method,"2 - simple relaxation");
GUI_TOOLTIP(uspex_gui._vcnebtype_method,"Choose between VC-NEB method and simple structure relaxation.");
	GUI_ENTRY_TABLE(table,uspex_gui.vcnebType,uspex_gui.calc.vcnebType,"%3i","VC-NEB:",1,2,5,6);
GUI_TOOLTIP(uspex_gui.vcnebType,"vcnebType: Ch. 6.2 DEFAULT: 110\nType of VC-NEB calculation.");
GUI_LOCK(uspex_gui.vcnebType);/*not directly modifiable*/
	GUI_CHECK_TABLE(table,uspex_gui._vcnebtype_img_num,uspex_gui.calc._vcnebtype_img_num,update_vcnebType,"Var_Image",2,3,5,6);
GUI_TOOLTIP(uspex_gui._vcnebtype_img_num,"Set whether number of images should be kept fixed.");
	GUI_CHECK_TABLE(table,uspex_gui._vcnebtype_spring,uspex_gui.calc._vcnebtype_spring,update_vcnebType,"Var_Spring",3,4,5,6);
GUI_TOOLTIP(uspex_gui._vcnebtype_spring,"Set whether spring constants should be kept fixed.");
/* line 6 */
	GUI_COMBOBOX_TABLE(table,uspex_gui.optReadImages,"Img:",0,1,6,7);
	GUI_COMBOBOX_ADD(uspex_gui.optReadImages,"0 - All structures are needed");
	GUI_COMBOBOX_ADD(uspex_gui.optReadImages,"1 - Only initial and final");
	GUI_COMBOBOX_ADD(uspex_gui.optReadImages,"2 - Initial and final + intermediates");
GUI_TOOLTIP(uspex_gui.optReadImages,"optReadImages: Ch. 6.2 DEFAULT: 2\nSet the method for reading Images file.");
	GUI_ENTRY_TABLE(table,uspex_gui.numImages,uspex_gui.calc.numImages,"%i","N_Img:",1,2,6,7);
GUI_TOOLTIP(uspex_gui.numImages,"numImages: Ch. 6.2 DEFAULT: 9\nInitial number of images.");
	GUI_ENTRY_TABLE(table,uspex_gui.numSteps,uspex_gui.calc.numSteps,"%4i","N_Step:",2,3,6,7);
GUI_TOOLTIP(uspex_gui.numSteps,"numSteps: Ch. 6.2 DEFAULT: 600\nMaximum VC-NEB step iterations.\nA value of at least 500 is recommended.");
	GUI_CHECK_TABLE(table,uspex_gui.optFreezing,uspex_gui.calc.optFreezing,NULL,"Freeze_Img",3,4,6,7);/*not calling anything*/
GUI_TOOLTIP(uspex_gui.optFreezing,"optFreezing: Ch. 6.2 DEFAULT: FALSE\nActivate freezing of Image structure when ConvThreshold is reached.");
/* line 7 */
	GUI_COMBOBOX_TABLE(table,uspex_gui.optimizerType,"Opt:",0,1,7,8);
	GUI_COMBOBOX_ADD(uspex_gui.optimizerType,"1 - Steepest Descent");
	GUI_COMBOBOX_ADD(uspex_gui.optimizerType,"2 - Fast Inertial Relaxation Engine");
GUI_TOOLTIP(uspex_gui.optimizerType,"optimizerType: Ch. 6.2 DEFAULT: 1\nSelect the optimization algorithm (SD or FIRE).");
	GUI_ENTRY_TABLE(table,uspex_gui.dt,uspex_gui.calc.dt,"%.4f","dt:",1,2,7,8);
GUI_TOOLTIP(uspex_gui.dt,"dt: Ch. 6.2 DEFAULT: 0.05\nTime step for structure relaxation.");
	GUI_ENTRY_TABLE(table,uspex_gui.ConvThreshold,uspex_gui.calc.ConvThreshold,"%.4f","Conv:",2,3,7,8);
GUI_TOOLTIP(uspex_gui.ConvThreshold,"ConvThreshold: Ch. 6.2 DEFAULT: 0.003\nHalting condition (ev/Ang.) for RMS forces on images.");
	GUI_ENTRY_TABLE(table,uspex_gui.VarPathLength,uspex_gui.calc.VarPathLength,"%.4f","PathLength:",3,4,7,8);
GUI_TOOLTIP(uspex_gui.VarPathLength,"VarPathLength: Ch. 6.2 DEFAULT: AUTO\nCriterion to determine image creation/deletion for variable image method.\nif L>1.5*criterion image is added\nif L<0.5*criterion image is deleted\nL=length between two neighbor images.");
/* line 8 */
	GUI_COMBOBOX_TABLE(table,uspex_gui.optRelaxType,"Relax:",0,1,8,9);
	GUI_COMBOBOX_ADD(uspex_gui.optRelaxType,"1 - fixed cell, positions relaxed (=NEB)");
	GUI_COMBOBOX_ADD(uspex_gui.optRelaxType,"2 - cell lattice only (only for testing)");
	GUI_COMBOBOX_ADD(uspex_gui.optRelaxType,"3 - full, cell and positions relaxation.");
GUI_TOOLTIP(uspex_gui.optRelaxType,"optRelaxType: Ch. 6.2 DEFAULT: 3\nStructure relaxation mode.");
	GUI_ENTRY_TABLE(table,uspex_gui.K_min,uspex_gui.calc.K_min,"%.4f","Kmin:",1,2,8,9);
GUI_TOOLTIP(uspex_gui.K_min,"K_min: Ch. 6.2 DEFAULT: 5\nMinimum spring constant (eV/Ang.^2).");
	GUI_ENTRY_TABLE(table,uspex_gui.K_max,uspex_gui.calc.K_max,"%.4f","Kmax:",2,3,8,9);
GUI_TOOLTIP(uspex_gui.K_max,"K_max: Ch. 6.2 DEFAULT: 5\nMaximum spring constant (eV/Ang.^2).");
	GUI_ENTRY_TABLE(table,uspex_gui.Kconstant,uspex_gui.calc.Kconstant,"%.4f","Kcte:",3,4,8,9);
GUI_TOOLTIP(uspex_gui.Kconstant,"Kconstant: Ch. 6.2 DEFAULT: 5\nFixed spring constant (eV/Ang.^2).");
/* line 9 */
	GUI_COMBOBOX_TABLE(table,uspex_gui.optMethodCIDI,"CI/DI:",0,1,9,10);
	GUI_COMBOBOX_ADD(uspex_gui.optMethodCIDI," 0 - No CI/DI method will be used");
	GUI_COMBOBOX_ADD(uspex_gui.optMethodCIDI," 1 - single CI on highest energy TS");
	GUI_COMBOBOX_ADD(uspex_gui.optMethodCIDI,"-1 - Single DI on lowest  energy LM");
	GUI_COMBOBOX_ADD(uspex_gui.optMethodCIDI," 2 - Multi-CI/DI on provided  TS/LM");
GUI_TOOLTIP(uspex_gui.optMethodCIDI,"optMethodCIDI: Ch. 6.2 DEFAULT: 0\nOption for Climbing-Image (CI) and Descending-Image (DI).");
	GUI_ENTRY_TABLE(table,uspex_gui.startCIDIStep,uspex_gui.calc.startCIDIStep,"%i","Start CI/DI:",1,2,9,10);
GUI_TOOLTIP(uspex_gui.startCIDIStep,"startCIDIStep: Ch. 6.2 DEFAULT: 100\nStarting step for CI/DI method.");
	if(uspex_gui.calc._npickimg<1){
		title=g_strdup("");
	}else{
		title=g_strdup_printf("%i",uspex_gui.calc.pickupImages[0]);
		for(idx=1;idx<uspex_gui.calc._npickimg;idx++){
			tmp=g_strdup_printf("%s %i",title,uspex_gui.calc.pickupImages[idx]);
			g_free(title);
			title=tmp;
		}
	}
	GUI_TEXT_TABLE(table,uspex_gui.pickupImages,title,"Pickup:",2,3,9,10);
	g_free(title);
GUI_TOOLTIP(uspex_gui.pickupImages,"pickupImages: Ch. 6.2 DEFAULT: AUTO\nNumber/which images to be picked up for CI/DI method.");
	GUI_ENTRY_TABLE(table,uspex_gui.PrintStep,uspex_gui.calc.PrintStep,"%i","PrintStep:",3,4,9,10);
GUI_TOOLTIP(uspex_gui.PrintStep,"PrintStep: Ch. 6.2 DEFAULT: 1\nSave restart file every PrintStep times.");
/* line 10 */
	GUI_COMBOBOX_TABLE(table,uspex_gui.FormatType,"Format:",0,1,10,11);
	GUI_COMBOBOX_ADD(uspex_gui.FormatType,"1 - XCRYSTDENS format (.xsf)");
	GUI_COMBOBOX_ADD(uspex_gui.FormatType,"2 - VASP v5 (POSCAR) format.");
	GUI_COMBOBOX_ADD(uspex_gui.FormatType,"3 - XYZ format with lattice.");
GUI_TOOLTIP(uspex_gui.FormatType,"FormatType: Ch. 6.2 DEFAULT: 2\nFormat of structures in PATH output directory.");
	GUI_COMBOBOX_TABLE(table,uspex_gui.img_model,"Model:",1,3,10,11);
GUI_TOOLTIP(uspex_gui.img_model,"Select the original Images model.");
	GUI_OPEN_BUTTON_TABLE(table,uspex_gui.img_model_button,load_img_model_file,3,4,10,11);
/* --- Molecules */
	GUI_LABEL_TABLE(table,"Molecules",0,4,11,12);
/* line 12 */
	GUI_COMBOBOX_TABLE(table,uspex_gui.mol_model,"MODEL:",0,1,12,13);
	GUI_COMBOBOX_ADD(uspex_gui.mol_model,"Already provided");
	GUI_COMBOBOX_ADD(uspex_gui.mol_model,"From MOL_ folder");
	GUI_COMBOBOX_ADD(uspex_gui.mol_model,"From GDIS models");
GUI_TOOLTIP(uspex_gui.mol_model,"Select the model(s) for molecular calculation.");
	GUI_SPIN_TABLE(table,uspex_gui.num_mol,uspex_gui._tmp_num_mol,spin_update_num_mol,"N_MOLS",1,2,12,13);
GUI_TOOLTIP(uspex_gui.num_mol,"Select the number of molecules.");
	/*col 3: empty*/
	GUI_OPEN_BUTTON_TABLE(table,uspex_gui.mol_model_button,load_mol_model_folder,3,4,12,13);
/* line 13 */
	GUI_COMBOBOX_TABLE(table,uspex_gui.mol_gdis,"GDIS_MOL:",0,1,13,14);
GUI_TOOLTIP(uspex_gui.mol_gdis,"Select the molecule from GDIS model.");
	GUI_SPIN_TABLE(table,uspex_gui.curr_mol,uspex_gui._tmp_curr_mol,spin_update_curr_mol,"MOL_#",1,2,13,14);
	GUI_SPIN_RANGE(uspex_gui.curr_mol,1.,(gdouble)uspex_gui.calc._nmolecules);
GUI_TOOLTIP(uspex_gui.curr_mol,"Current molecule number.");
	GUI_CHECK_TABLE(table,uspex_gui.mol_gulp,uspex_gui.mol_as_gulp,NULL,"GULP_FORM",2,3,13,14);/*not calling anything*/
GUI_TOOLTIP(uspex_gui.mol_gulp,"Use GULP chemical labels and charge Zmatrix form.");
	GUI_APPLY_BUTTON_TABLE(table,uspex_gui.mol_apply_button,apply_gdis_mol,3,4,13,14);
/* --- Surfaces */
	GUI_LABEL_TABLE(table,"Surfaces",0,4,14,15);
/* line 15 */
	GUI_COMBOBOX_TABLE(table,uspex_gui.substrate_model,"MODEL:",1,3,15,16);
GUI_TOOLTIP(uspex_gui.substrate_model,"Select the model to use as a substrate\nie. without buffer, surface, and vacuum region.");
	GUI_OPEN_BUTTON_TABLE(table,uspex_gui.substrate_model_button,load_substrate_model_file,3,4,15,16);
/*line 16 */
	GUI_ENTRY_TABLE(table,uspex_gui.reconstruct,uspex_gui.calc.reconstruct,"%i","N_Surf:",1,2,16,17);
GUI_TOOLTIP(uspex_gui.reconstruct,"reconstruct: Ch. 5.4 DEFAULT: 1\nNumber of replication of the surface cell.");
	GUI_ENTRY_TABLE(table,uspex_gui.thicknessS,uspex_gui.calc.thicknessS,"%.4f","S_thick:",2,3,16,17);
GUI_TOOLTIP(uspex_gui.thicknessS,"thicknessS: Ch. 5.4 DEFAULT: 2.0\nThickness (Ang.) of the surface region.");
        GUI_ENTRY_TABLE(table,uspex_gui.thicknessB,uspex_gui.calc.thicknessB,"%.4f","B_thick:",3,4,16,17);
GUI_TOOLTIP(uspex_gui.thicknessB,"thicknessB: Ch. 5.4 DEFAULT: 3.0\nThickness (Ang.) of the buffer region.");
/* line 17 */
	if(uspex_gui.calc.StoichiometryStart==NULL){
		title=g_strdup("");
	}else{
		title=g_strdup_printf("%i",uspex_gui.calc.StoichiometryStart[0]);
		for(idx=1;idx<uspex_gui.calc._nspecies;idx++){
			tmp=g_strdup_printf("%s %i",title,uspex_gui.calc.StoichiometryStart[idx]);
			g_free(title);
			title=tmp;
		}
	}
	GUI_TEXT_TABLE(table,uspex_gui.StoichiometryStart,title,"Stoichio:",0,1,17,18);
	g_free(title);
GUI_TOOLTIP(uspex_gui.StoichiometryStart,"StoichiometryStart: Ch. 5.4 DEFAULT: ?\nDefine the initial stoichiometry of the BULK.");
	GUI_ENTRY_TABLE(table,uspex_gui.E_AB,uspex_gui.calc.E_AB,"%.4f","E_AB:",1,2,17,18);
GUI_TOOLTIP(uspex_gui.E_AB,"E_AB: Ch. 5.4 DEFAULT: ?\nDFT energy (eV/formula) of AmBn.");
	GUI_ENTRY_TABLE(table,uspex_gui.Mu_A,uspex_gui.calc.Mu_A,"%.4f","Mu_A:",2,3,17,18);
GUI_TOOLTIP(uspex_gui.Mu_A,"Mu_A: Ch. 5.4 DEFAULT: ?\nDFT energy (eV/atom) of elemental A.");
	GUI_ENTRY_TABLE(table,uspex_gui.Mu_B,uspex_gui.calc.Mu_B,"%.4f","Mu_B:",3,4,17,18);
GUI_TOOLTIP(uspex_gui.Mu_B,"Mu_B: Ch. 5.4 DEFAULT: ?\nDFT energy (eV/atom) of elemental B.");
/* --- end page */
/* initialize everything */
	GUI_COMBOBOX_SETUP(uspex_gui.calculationMethod,0,uspex_method_selected);
	GUI_COMBOBOX_SETUP(uspex_gui.calculationType,0,uspex_type_selected);
	GUI_COMBOBOX_SETUP(uspex_gui.optType,uspex_gui.calc.optType-1,uspex_optimization_selected);
	GUI_LOCK(uspex_gui._atom_typ);
	populate_atomType();
	GUI_COMBOBOX_SETUP(uspex_gui.atomType,uspex_gui.calc._nspecies,atomType_selected);
	GUI_SPIN_SET(uspex_gui._calctype_dim,(gdouble) uspex_gui.calc._calctype_dim);
	mol_toggle();
	var_toggle();
	GUI_COMBOBOX_SETUP(uspex_gui.goodBonds,0,goodBonds_selected);
	auto_bond_toggle();
	opt_toggle();
	uspex_method_selected(uspex_gui.calculationMethod);
	GUI_COMBOBOX_SETUP(uspex_gui.IonDistances,0,uspex_IonDistances_selected);
	GUI_COMBOBOX_SETUP(uspex_gui.MolCenters,0,uspex_MolCenters_selected);
	refresh_constraints();
	toggle_auto_C_ion();
	GUI_COMBOBOX_SETUP(uspex_gui._latticeformat,0,uspex_latticeformat_selected);
	set_lattice_format();
	toggle_auto_C_lat();
	SG_toggle();
	set_numSpecies();
	GUI_COMBOBOX_SETUP(uspex_gui.numSpecies,0,uspex_numSpecies_selected);
	uspex_numSpecies_selected(uspex_gui.numSpecies);
	toggle_v1030();
/*per optimization step*/
	GUI_COMBOBOX_SETUP(uspex_gui.abinitioCode,0,uspex_ai_selected);/*ai means ab initio, not...*/
	spin_update_curr_step();
	uspex_ai_selected(uspex_gui.abinitioCode);
	toggle_auto_step();
	mag_toggle();
	GUI_COMBOBOX_SETUP(uspex_gui.dynamicalBestHM,2,uspex_dyn_HM_selected);
	GUI_COMBOBOX_SETUP(uspex_gui.manyParents,0,uspex_manyParents_selected);
	GUI_COMBOBOX_SETUP(uspex_gui.TE_goal,0,uspex_TE_goal_selected);
	GUI_COMBOBOX_SETUP(uspex_gui.FullRelax,2,uspex_relax_selected);
	GUI_COMBOBOX_SETUP(uspex_gui.meta_model,0,uspex_meta_model_selected);
	init_metadynamics_models();
	GUI_COMBOBOX_SETUP(uspex_gui._vcnebtype_method,0,uspex_vcneb_method_selected);
	GUI_COMBOBOX_SETUP(uspex_gui.optReadImages,2,uspex_ReadImg_selected);
	GUI_COMBOBOX_SETUP(uspex_gui.optimizerType,0,uspex_optimizerType_selected);
	GUI_COMBOBOX_SETUP(uspex_gui.optRelaxType,2,uspex_RelaxType_selected);
	GUI_COMBOBOX_SETUP(uspex_gui.optMethodCIDI,0,uspex_CIDI_selected);
	GUI_COMBOBOX_SETUP(uspex_gui.FormatType,0,uspex_FormatType_selected);
	uspex_CIDI_selected(uspex_gui.optMethodCIDI);
	GUI_COMBOBOX_SETUP(uspex_gui.img_model,0,uspex_img_model_selected);
	init_img_gdis_model();
	update_vcnebType();
	GUI_COMBOBOX_SETUP(uspex_gui.mol_model,0,uspex_mol_model_selected);
	init_uspex_gdis_mol();
	GUI_COMBOBOX_SETUP(uspex_gui.mol_gdis,0,uspex_gdis_mol_selected);
	GUI_COMBOBOX_SETUP(uspex_gui.substrate_model,0,uspex_substrate_model_selected);
	init_substrate_models();
/* --- Outside of notebook */
	GUI_FRAME_WINDOW(uspex_gui.window,frame);
	GUI_VBOX_FRAME(frame,vbox);
/* Action buttons */
	GUI_SAVE_ACTION(uspex_gui.window,uspex_gui.button_save,save_uspex_calc,NULL);
	GUI_EXEC_ACTION(uspex_gui.window,uspex_gui.button_exec,uspex_exec_calc,NULL);
	GUI_CLOSE_ACTION(uspex_gui.window,button,quit_uspex_gui,dialog);
/* connect to signals */
	GUI_PAGE_CHANGE(notebook,uspex_gui_page_switch,NULL);
/* all done, we need to copy uspex parameters from model again */
if(uspex_gui.have_output) {
		copy_uspex_parameters(uspex_output->calc,&(uspex_gui.calc));
		toggle_use_specific();
}else{
		toggle_set_specific();
}
	uspex_gui_refresh();/*refresh once more*/
	GUI_SHOW(uspex_gui.window);/*display*/
	update_specific();
	sysenv.refresh_dialog=TRUE;
}

