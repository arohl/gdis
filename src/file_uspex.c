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

/* USPEX Parser:
 * The reading of USPEX result file is a little different:
 * 	1/ make a graph representing all structures energy levels vs.  generation numbers
 *		-> make each point on that graph clickable so that selecting a point will
 *		   open a new model containing the corresponding structure.
 *	2/ load the optimal structures as base model.
 *		-> each optimal structure, for each generation, is represented as a frame
 *		   in this model.
 * TODO:
 * 	1/ include other data plots:
 * 		- hardness / force vs. generation
 * 	2/ calculate some other properties...
 * */

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>

#include "gdis.h"
#include "coords.h"
#include "error.h"
#include "file.h"
#include "parse.h"
#include "matrix.h"
#include "zmatrix.h"
#include "zmatrix_pak.h"
#include "model.h"
#include "interface.h"
#include "graph.h"
#include "file_vasp.h"
#include "file_uspex.h"

/* main structures */
extern struct sysenv_pak sysenv;
extern struct elem_pak elements[];

#define DEBUG_USPEX_READ 0
#define _UO (*uspex_output)

gint read_parameter_uspex(gchar *filename, struct model_pak *model){
	FILE *vf;
        gchar *line=NULL;
	gchar *ptr;
	gint idx;
        uspex_output_struct *uspex_output=model->uspex;
	/* specific */
	gchar method;
	/*start*/
	_UO.nspecies=0;
        vf = fopen(filename, "rt");
        if (!vf) return 1;
        line = file_read_line(vf);
        while (line){
                if (find_in_string("calculationMethod",line) != NULL) {
                        sscanf(line,"%c%*s : calculationMethod %*s",&(method));
                        switch (method){
                                case 'U':
                                case 'u':
                                /*uspex*/
				_UO.method=US_CM_USPEX;
				property_add_ranked(2, "Calculation", "USPEX", model);
#if DEBUG_USPEX_READ
fprintf(stdout,"#DBG: USPEX calculationMethod=USPEX\n");
#endif
				break;
                                case 'M':
                                case 'm':
                                /*metadynamics*/
				_UO.method=US_CM_META;
#if DEBUG_USPEX_READ
fprintf(stdout,"#DBG: USPEX calculationMethod=META\n");
#endif
				break;
                                case 'V':
                                case 'v':
                                /*variable cell nudged elastic band*/
				_UO.method=US_CM_VCNEB;
#if DEBUG_USPEX_READ
fprintf(stdout,"#DBG: USPEX calculationMethod=VCNEB\n");
#endif
				break;
                                case 'P':
                                case 'p':
                                /*particle swarm optimization*/
				_UO.method=US_CM_PSO;
#if DEBUG_USPEX_READ
fprintf(stdout,"#DBG: USPEX calculationMethod=PSO\n");
#endif
				break;
                                default:
                                        g_free(line);
                                        line = g_strdup_printf("ERROR: unsupported USPEX method!\n");
                                        gui_text_show(ERROR, line);
                                        g_free(line);
					return 1;
                        }
                }
                if (find_in_string("calculationType",line) != NULL) {
                        sscanf(line,"%i : calculationType %*s",&(_UO.type));
#if DEBUG_USPEX_READ
fprintf(stdout,"#DBG: USPEX calculationType=%i\n",_UO.type);
#endif
                }
		if (find_in_string("optType",line) != NULL) {
			sscanf(line,"%i : optType %*s",&(_UO.opt_type));
#if DEBUG_USPEX_READ
fprintf(stdout,"#DBG: USPEX optType=%i\n",_UO.type);
#endif

		}
		if (find_in_string("atomType",line) !=NULL){
		  /*let's count number of species here*/
		  g_free(line);
		  line = file_read_line(vf);
		  ptr=&(line[0]);
		  while(*ptr==' ') ptr++;/*skip initial space (if any)*/
		  _UO.nspecies=1;/*there is at least one species*/
		  while(*ptr!='\0'){
		  	if(*ptr==' ') {
	 	 		_UO.nspecies++;ptr++;
				while(g_ascii_isgraph(*ptr)) ptr++;/*go to next space/end*/
 	 		}else ptr++;
  		  }
#if DEBUG_USPEX_READ
fprintf(stdout,"#DBG: USPEX number_of_species=%i\n",_UO.nspecies);
#endif
		  /*now let's find each species symbol*/
		  _UO.spe_Z=g_malloc(_UO.nspecies*sizeof(gint));
		  ptr=&(line[0]);idx=0;
		  while(*ptr==' ') ptr++;/*skip initial space (if any)*/
		  while((*ptr!='\n')&&(*ptr!='\0')){
		  	if(*ptr==' ') while(*ptr==' ') ptr++;/*skip space(s)*/
        		if(g_ascii_isdigit(*ptr)) _UO.spe_Z[idx]=g_ascii_strtod(ptr,NULL);/*user provided Z number*/
			else _UO.spe_Z[idx]=elem_symbol_test(ptr);/*user provided symbol*/
			if((_UO.spe_Z[idx]<=0)||(_UO.spe_Z[idx]>MAX_ELEMENTS-1)){/*X not allowed*/
				g_free(line);
				line = g_strdup_printf("ERROR: USPEX Parameters.txt file contain invalid atom definition!\n");
				gui_text_show(ERROR, line);
				g_free(line);
				return -1;
			}
#if DEBUG_USPEX_READ
fprintf(stdout,"#DBG: USPEX species[%i] Z=%i\n",idx,_UO.spe_Z[idx]);
#endif
			ptr++;idx++;
			while(g_ascii_isgraph(*ptr)) ptr++;/*go to next space/end*/
		    }
		}
	g_free(line);
	line = file_read_line(vf);
        }
/*other parameters will be filled later for gui_uspex */



/*end of Parameters.txt read*/
	return 0;
}
/************************************************************************/
/* Read Parameters.txt and fill keywords of uspex_calc_struct structure */
/*For reading Parameters.txt files, *not* preparing a model for display.*/
/************************************************************************/
gint read_uspex_parameters(gchar *filename,uspex_calc_struct *calc){
#define CALC (*calc)
	FILE *vf;
	long int vfpos;
	gchar *line=NULL;
	gchar *ptr;
//	gchar tmp[64];
	gchar c;
	gint i,j,k;
//	gdouble d;
	/*tests*/
	if(filename==NULL) return -1;
	if(calc==NULL) return -1;
	vf = fopen(filename, "rt");
	if (!vf) return -2;
	/**/
	CALC._nspecies=0;/*set hidden to 0*/
/*some lazy defines*/
#define __NSPECIES(line) do{\
	if(CALC._nspecies==0){\
		ptr=&(line[0]);\
		while(*ptr==' ') ptr++;\
		CALC._nspecies=1;\
		while(*ptr!='\0'){\
			if(*ptr==' ') {\
				CALC._nspecies++;ptr++;\
				while(g_ascii_isgraph(*ptr)) ptr++;\
			}else ptr++;\
		}\
	}\
}while(0)
#define __Q(a) #a

#define __GET_BOOL(value) if (find_in_string(__Q(value),line)!=NULL){\
	k=0;sscanf(line,"%i%*s",&(k));\
	CALC.value=(k==1);\
	g_free(line);line = file_read_line(vf);\
	continue;\
}
#define __GET_INT(value) if (find_in_string(__Q(value),line) != NULL) {\
	sscanf(line,"%i%*s",&(CALC.value));\
	g_free(line);line = file_read_line(vf);\
	continue;\
}
#define __GET_DOUBLE(value) if (find_in_string(__Q(value),line) != NULL) {\
	sscanf(line,"%lf%*s",&(CALC.value));\
	g_free(line);line = file_read_line(vf);\
	continue;\
}
#define __GET_STRING(value) if (find_in_string(__Q(value),line) != NULL) {\
	g_free(line);line = file_read_line(vf);\
	CALC.value=g_strdup(line);\
	g_free(line);line = file_read_line(vf);\
	g_free(line);line = file_read_line(vf);\
	continue;\
}
#define __GET_CHARS(value) if (find_in_string(__Q(value),line) != NULL) {\
	ptr = g_strdup_printf("%%s : %s",__Q(value));\
	sscanf(line,ptr,&(CALC.value));\
	g_free(ptr);g_free(line);line = file_read_line(vf);\
	continue;\
}

	line = file_read_line(vf);
	while (line){
		if (find_in_string("calculationMethod",line) != NULL) {
			c='\0';sscanf(line,"%c%*s",&(c));
			switch(c){
			case 'u':
			case 'U':/*USPEX method*/
				CALC.calculationMethod=US_CM_USPEX;
				break;
			case 'm':
			case 'M':
				CALC.calculationMethod=US_CM_META;
				break;
			case 'v':
			case 'V':
				CALC.calculationMethod=US_CM_VCNEB;
				break;
			case 'p':
			case 'P':
				/*PSO method is not supported yet*/
				CALC.calculationMethod=US_CM_PSO;
				break;
			default:
				CALC.calculationMethod=US_CM_UNKNOWN;
			}
			g_free(line);
			line = file_read_line(vf);
			continue;
		}
		if (find_in_string("calculationType",line) != NULL) {
			k=0;sscanf(line,"%i%*s",&(k));
			j=k;
			i=((int)(j/100))%10;
			if((i<0)||(i>3)) {
				g_free(line);
				line = file_read_line(vf);
				continue;
			}
			CALC._calctype_dim=i;
			j-=i*100;
			i=((int)(j/10))%10;
			if((i<0)||(i>1)){
				g_free(line);
				line = file_read_line(vf);
				continue;
			}
			CALC._calctype_mol=(i==1);
			j-=i*10;
			i=j%10;
			if((i<0)||(i>1)){
				g_free(line);
				line = file_read_line(vf);
				continue;
			}
			CALC._calctype_var=(i==1);
			CALC.calculationType=k;
			g_free(line);
			line = file_read_line(vf);
			continue;
		}
		if (find_in_string("optType",line) != NULL) {
			k=0;sscanf(line,"%i%*s",&(k));
			switch (k){
			case 1:
				CALC.optType=US_OT_ENTHALPY;
				break;
			case 2:
				CALC.optType=US_OT_VOLUME;
				break;
			case 3:
				CALC.optType=US_OT_HARDNESS;
				break;
			case 4:
				CALC.optType=US_OT_ORDER;
				break;
			case 5:
				CALC.optType=US_OT_DISTANCE;
				break;
			case 6:
				CALC.optType=US_OT_DIELEC_S;
				break;
			case 7:
				CALC.optType=US_OT_GAP;
				break;
			case 8:
				CALC.optType=US_OT_DIELEC_GAP;
				break;
			case 9:
				CALC.optType=US_OT_MAG;
				break;
			case 10:
				CALC.optType=US_OT_QE;
				break;
			case 1101:
				CALC.optType=US_OT_BULK_M;
				break;
			case 1102:
				CALC.optType=US_OT_SHEAR_M;
				break;
			case 1103:
				CALC.optType=US_OT_YOUNG_M;
				break;
			case 1104:
				CALC.optType=US_OT_POISSON;
				break;
			case 1105:
				CALC.optType=US_OT_PUGH_R;
				break;
			case 1106:
				CALC.optType=US_OT_VICKERS_H;
				break;
			case 1107:
				CALC.optType=US_OT_FRACTURE;
				break;
			case 1108:
				CALC.optType=US_OT_DEBYE_T;
				break;
			case 1109:
				CALC.optType=US_OT_SOUND_V;
				break;
			case 1110:
				CALC.optType=US_OT_SWAVE_V;
				break;
			case 1111:
				CALC.optType=US_OT_PWAVE_V;
				break;
			default:
				CALC.optType=US_OT_UNKNOWN;
			}
			g_free(line);
			line = file_read_line(vf);
			continue;
		}
		if (find_in_string("atomType",line) !=NULL){
			g_free(line);line = file_read_line(vf);/*go next line*/
			/*if no _nspecies, get it now*/
			__NSPECIES(line);
			/*look for atomType*/
			CALC.atomType = g_malloc(CALC._nspecies*sizeof(gint));
			for(i=0;i<CALC._nspecies;i++) CALC.atomType[i]=0;
			ptr=&(line[0]);i=0;
//			while(*ptr==' ') ptr++;/*skip initial space (if any)*/
			while((*ptr!='\n')&&(*ptr!='\0')){
				if(*ptr==' ') while(*ptr==' ') ptr++;/*skip space(s)*/
				if(g_ascii_isdigit(*ptr)) CALC.atomType[i]=g_ascii_strtod(ptr,NULL);/*user provided Z number*/
				else CALC.atomType[i]=elem_symbol_test(ptr);/*try user provided symbol*/
				if((CALC.atomType[i]<=0)||(CALC.atomType[i]>MAX_ELEMENTS-1)){/*invalid Z*/
					CALC.atomType[i]=0;/*allow it for now*/
				}
				ptr++;i++;
				while(g_ascii_isgraph(*ptr)) ptr++;/*go to next space/end*/
			}
			g_free(line);
			line = file_read_line(vf);/*this is the EndAtomType line*/
			g_free(line);
			line = file_read_line(vf);
			continue;
		}
		if (find_in_string("numSpecies",line) != NULL) {
			vfpos=ftell(vf);/* flag */
			g_free(line);line = file_read_line(vf);/*go next line*/
			/*if no _nspecies, get it now*/
			__NSPECIES(line);
			/*look for numSpecies*/
			/*1st get the total number of lines*/
			CALC._var_nspecies=0;
			do{
				g_free(line);line = file_read_line(vf);/*go next line*/
				CALC._var_nspecies++;
			}while(find_in_string("umSpecies",line) == NULL);
			fseek(vf,vfpos,SEEK_SET);/* rewind to flag */
			line = file_read_line(vf);/*first line of numSpecies*/
			CALC.numSpecies = g_malloc(CALC._nspecies*CALC._var_nspecies*sizeof(gint));
			for(i=0;i<CALC._nspecies*CALC._var_nspecies;i++) CALC.numSpecies[i]=0;
			for(j=0;j<CALC._var_nspecies;j++){
				ptr=&(line[0]);i=0;
				while((*ptr!='\n')&&(*ptr!='\0')){
					if(*ptr==' ') while(*ptr==' ') ptr++;/*skip space(s)*/
					sscanf(ptr,"%i%*s",&(CALC.numSpecies[i+j*CALC._nspecies]));
					ptr++;i++;
					while(g_ascii_isgraph(*ptr)) ptr++;/*go to next space/end*/
				}
				g_free(line);
				line = file_read_line(vf);
			}
			if (find_in_string("umSpecies",line) != NULL) {
				g_free(line);
				line = file_read_line(vf);
			}
			g_free(line);
			line = file_read_line(vf);
			continue;
		}
		__GET_DOUBLE(ExternalPressure)
		if (find_in_string("valences",line) != NULL) {
			g_free(line);line = file_read_line(vf);/*go next line*/
			/*if no _nspecies, get it now*/
			__NSPECIES(line);
			/*look for valences*/
			CALC.valences = g_malloc(CALC._nspecies*sizeof(gint));
			for(i=0;i<CALC._nspecies;i++) CALC.valences[i]=0;
			ptr=&(line[0]);i=0;
			while((*ptr!='\n')&&(*ptr!='\0')){
				if(*ptr==' ') while(*ptr==' ') ptr++;/*skip space(s)*/
				sscanf(ptr,"%i%*s",&(CALC.valences[i]));
				ptr++;i++;
				while(g_ascii_isgraph(*ptr)) ptr++;/*go to next space/end*/
			}
			g_free(line);
			line = file_read_line(vf);/*this is the EndValences line*/
			g_free(line);
			line = file_read_line(vf);
			continue;
		}
		if (find_in_string("goodBonds",line) != NULL) {
			g_free(line);line = file_read_line(vf);/*go next line*/
			/*if no _nspecies, get it now*/
			__NSPECIES(line);
			/*look for goodBonds*/
			CALC.goodBonds = g_malloc(CALC._nspecies*CALC._nspecies*sizeof(gdouble));
			for(i=0;i<CALC._nspecies*CALC._nspecies;i++) CALC.goodBonds[i]=0.;
			j=0;
			while((j<CALC._nspecies)&&(find_in_string("oodBonds",line)==NULL)){
				ptr=&(line[0]);i=0;
				while((*ptr!='\n')&&(*ptr!='\0')){
					if(*ptr==' ') while(*ptr==' ') ptr++;/*skip space(s)*/
					sscanf(ptr,"%lf%*s",&(CALC.goodBonds[i+j*CALC._nspecies]));
					ptr++;i++;
					while(g_ascii_isgraph(*ptr)) ptr++;/*go to next space/end*/
				}
				j++;
				g_free(line);
				line = file_read_line(vf);
			}
			if(find_in_string("oodBonds",line) != NULL){
				g_free(line);
				line = file_read_line(vf);
			}
			g_free(line);
			line = file_read_line(vf);
			continue;
		}
		__GET_BOOL(checkMolecules);
		__GET_BOOL(checkConnectivity);
		__GET_INT(populationSize);
		__GET_INT(initialPopSize);
		__GET_INT(numGenerations);
		__GET_INT(stopCrit);
		__GET_DOUBLE(bestFrac);
		__GET_INT(keepBestHM);
		__GET_BOOL(reoptOld);
		__GET_STRING(symmetries);
		__GET_DOUBLE(fracGene);
		__GET_DOUBLE(fracRand);
		__GET_DOUBLE(fracPerm);
		__GET_DOUBLE(fracAtomsMut);
		__GET_DOUBLE(fracRotMut);
		__GET_DOUBLE(fracLatMut);
		__GET_INT(howManySwaps);
		__GET_STRING(specificSwaps);
		__GET_DOUBLE(mutationDegree);
		__GET_DOUBLE(mutationRate);
		__GET_DOUBLE(DisplaceInLatmutation);
		__GET_BOOL(AutoFrac);
		__GET_BOOL(minVectorLength);
		if (find_in_string("IonDistances",line) != NULL) {
			g_free(line);line = file_read_line(vf);/*go next line*/
			/*if no _nspecies, get it now*/
			__NSPECIES(line);
			/*look for IonDistances*/
			CALC.IonDistances = g_malloc(CALC._nspecies*CALC._nspecies*sizeof(gdouble));
			for(i=0;i<CALC._nspecies*CALC._nspecies;i++) CALC.IonDistances[i]=0.;
			j=0;
			while((j<CALC._nspecies)&&(find_in_string("onDistances",line)==NULL)){
				ptr=&(line[0]);i=0;
				while((*ptr!='\n')&&(*ptr!='\0')){
					if(*ptr==' ') while(*ptr==' ') ptr++;/*skip space(s)*/
					sscanf(ptr,"%lf%*s",&(CALC.IonDistances[i+j*CALC._nspecies]));
					ptr++;i++;
					while(g_ascii_isgraph(*ptr)) ptr++;/*go to next space/end*/
				}
				j++;
				g_free(line);
				line = file_read_line(vf);
			}
			if(find_in_string("onDistances",line) != NULL){
				g_free(line);
				line = file_read_line(vf);
			}
			g_free(line);
			line = file_read_line(vf);
			continue;
		}
		__GET_INT(constraint_enhancement);
		if (find_in_string("MolCenters",line) != NULL) {
			g_free(line);line = file_read_line(vf);/*go next line*/
			/*count number of molecules*/
			ptr=&(line[0]);
			while(*ptr==' ') ptr++;
			CALC._nmolecules=1;
			while(*ptr!='\0'){
				if(*ptr==' ') {
					CALC._nmolecules++;
					while(g_ascii_isgraph(*ptr)) ptr++;
				}else ptr++;
			}
			CALC.MolCenters = g_malloc(CALC._nmolecules*CALC._nmolecules*sizeof(gdouble));
			for(i=0;i<CALC._nmolecules;i++) CALC.MolCenters[i]=0.;
			j=0;
			while((j<CALC._nmolecules)&&(find_in_string("EndMol",line)==NULL)){
				ptr=&(line[0]);i=0;
				while((*ptr!='\n')&&(*ptr!='\0')){
					if(*ptr==' ') while(*ptr==' ') ptr++;/*skip space(s)*/
					sscanf(ptr,"%lf%*s",&(CALC.MolCenters[i+j*CALC._nmolecules]));
					ptr++;i++;
					while(g_ascii_isgraph(*ptr)) ptr++;/*go to next space/end*/
				}
				j++;
				g_free(line);
				line = file_read_line(vf);
			}
			if(find_in_string("EndMol",line) != NULL){
				g_free(line);
				line = file_read_line(vf);
			}
			g_free(line);
			line = file_read_line(vf);
			continue;
		}
		if (find_in_string("Latticevalues",line) != NULL) {
			/*Latticevalues can be one or several volumes, one line crystal definition, or 3 lines lattice cell vectors...*/
			vfpos=ftell(vf);/* flag */
			g_free(line);line = file_read_line(vf);/*go next line*/
			/*count first line number of items*/
			ptr=&(line[0]);
			while(*ptr==' ') ptr++;
			CALC._nlatticevalues=1;
			while(*ptr!='\0'){
				if(*ptr==' ') {
					CALC._nlatticevalues++;
					while(g_ascii_isgraph(*ptr)) ptr++;
				}else ptr++;
			}
			/*depending on varcomp*/
			if(CALC._calctype_var){
				/*varcomp calculation: should be only one line*/
				g_free(line);line = file_read_line(vf);/*go next line*/
				if(find_in_string("Endvalues",line) != NULL){/*only volumes values*/
					CALC.Latticevalues = g_malloc(CALC._nlatticevalues*sizeof(gdouble));
					for(i=0;i<CALC._nlatticevalues;i++) CALC.Latticevalues[i]=0.;
					fseek(vf,vfpos,SEEK_SET);/* rewind to flag */
					g_free(line);line = file_read_line(vf);/*go next line*/
					ptr=&(line[0]);i=0;
					while((*ptr!='\n')&&(*ptr!='\0')){
						if(*ptr==' ') while(*ptr==' ') ptr++;/*skip space(s)*/
						sscanf(ptr,"%lf%*s",&(CALC.Latticevalues[i]));
						ptr++;i++;
						while(g_ascii_isgraph(*ptr)) ptr++;/*go to next space/end*/
					}
				}else{/*varcomp but more than one line -> _nlatticevalues <= 3 or bad setting*/
					if((CALC._nlatticevalues<=3)&&(CALC._nlatticevalues>1)){
						/*this may be a definition of 1, 2, or 3D cell, but why?*/
						CALC.Latticevalues = g_malloc(CALC._nlatticevalues*CALC._nlatticevalues*sizeof(gdouble));
						for(i=0;i<CALC._nlatticevalues*CALC._nlatticevalues;i++) CALC.Latticevalues[i]=0.;
						fseek(vf,vfpos,SEEK_SET);/* rewind to flag */
						g_free(line);line = file_read_line(vf);/*go next line*/
						ptr=&(line[0]);
						j=0;
						while((j<CALC._nlatticevalues)&&(find_in_string("Endvalues",line)==NULL)){
							ptr=&(line[0]);i=0;
							while((*ptr!='\n')&&(*ptr!='\0')){
								if(*ptr==' ') while(*ptr==' ') ptr++;/*skip space(s)*/
								sscanf(ptr,"%lf%*s",&(CALC.Latticevalues[i+j*CALC._nlatticevalues]));
								ptr++;i++;
								while(g_ascii_isgraph(*ptr)) ptr++;/*go to next space/end*/
							}
							j++;
							g_free(line);
							line = file_read_line(vf);
						}
					}else{/*this is a bad setting*/
						g_free(line);
						line = g_strdup_printf("ERROR: USPEX bad Latticevalues/varcomp setting in Parameters.txt!\n");
						gui_text_show(ERROR, line);
						g_free(line);
						g_free(line);line = file_read_line(vf);/*go next line*/
					}
				}
			}else{
				/*NOT varcomp calculation: should be a crystal definition*/
				if(CALC._nlatticevalues==6){/*crystal definition*/
					CALC.Latticevalues = g_malloc(CALC._nlatticevalues*sizeof(gdouble));
					for(i=0;i<CALC._nlatticevalues;i++) CALC.Latticevalues[i]=0.;
					fseek(vf,vfpos,SEEK_SET);/* rewind to flag */
					g_free(line);line = file_read_line(vf);/*go next line*/
					ptr=&(line[0]);i=0;
					while((*ptr!='\n')&&(*ptr!='\0')){
						if(*ptr==' ') while(*ptr==' ') ptr++;/*skip space(s)*/
						sscanf(ptr,"%lf%*s",&(CALC.Latticevalues[i]));
						ptr++;i++;
						 while(g_ascii_isgraph(*ptr)) ptr++;/*go to next space/end*/
					}
				}else if((CALC._nlatticevalues<=3)&&(CALC._nlatticevalues>1)){/*lattice vectors*/
					CALC.Latticevalues = g_malloc(CALC._nlatticevalues*CALC._nlatticevalues*sizeof(gdouble));
					for(i=0;i<CALC._nlatticevalues*CALC._nlatticevalues;i++) CALC.Latticevalues[i]=0.;
					fseek(vf,vfpos,SEEK_SET);/* rewind to flag */
					g_free(line);line = file_read_line(vf);/*go next line*/
					ptr=&(line[0]);
					j=0;
					while((j<CALC._nlatticevalues)&&(find_in_string("Endvalues",line)==NULL)){
						ptr=&(line[0]);i=0;
						while((*ptr!='\n')&&(*ptr!='\0')){
							if(*ptr==' ') while(*ptr==' ') ptr++;/*skip space(s)*/
							sscanf(ptr,"%lf%*s",&(CALC.Latticevalues[i+j*CALC._nlatticevalues]));
							ptr++;i++;
							while(g_ascii_isgraph(*ptr)) ptr++;/*go to next space/end*/
						}
						j++;
						g_free(line);
						line = file_read_line(vf);
					}
				}else if(CALC._nlatticevalues==1){/*only the volume value*/
					CALC.Latticevalues = g_malloc(sizeof(gdouble));
					CALC.Latticevalues[0]=0.;
					sscanf(ptr,"%lf%*s",&(CALC.Latticevalues[0]));
					g_free(line);
					line = file_read_line(vf);
				}else{/*bad setting*/
					g_free(line);
					line = g_strdup_printf("ERROR: USPEX bad Latticevalues setting in Parameters.txt!\n");
					gui_text_show(ERROR, line);
					g_free(line);
					g_free(line);line = file_read_line(vf);/*go next line*/
				}
			}
			if(find_in_string("Endvalues",line) != NULL){
				g_free(line);
				line = file_read_line(vf);
			}
			g_free(line);
			line = file_read_line(vf);
			continue;
		}
		if (find_in_string("splitInto",line) != NULL) {
			g_free(line);line = file_read_line(vf);/*go next line*/
			/*count number of splits*/
			ptr=&(line[0]);
			while(*ptr==' ') ptr++;
			CALC._nsplits=1;
			while(*ptr!='\0'){
				if(*ptr==' ') {
					CALC._nsplits++;
					while(g_ascii_isgraph(*ptr)) ptr++;
				}else ptr++;
			}
			/*look for splits*/
			CALC.splitInto = g_malloc(CALC._nsplits*sizeof(gint));
			for(i=0;i<CALC._nsplits;i++) CALC.splitInto[i]=0;
			ptr=&(line[0]);i=0;
			while((*ptr!='\n')&&(*ptr!='\0')){
				if(*ptr==' ') while(*ptr==' ') ptr++;/*skip space(s)*/
				sscanf(ptr,"%i%*s",&(CALC.splitInto[i]));
				ptr++;i++;
				while(g_ascii_isgraph(*ptr)) ptr++;/*go to next space/end*/
			}
			g_free(line);
			line = file_read_line(vf);/*this is the EndSplitInto line*/
			g_free(line);
			line = file_read_line(vf);
			continue;
		}
		if (find_in_string("abinitioCode",line) != NULL) {
			g_free(line);line = file_read_line(vf);/*go next line*/
			/*count number of optimization steps*/
			ptr=&(line[0]);
			while(*ptr==' ') ptr++;
			CALC._num_opt_steps=1;
			while(*ptr!='\0'){
				if(*ptr==' ') {
					CALC._num_opt_steps++;
					while(g_ascii_isgraph(*ptr)) ptr++;
				}else ptr++;
			}/*note that the METADYNAMICS parenthesis are ignored as well*/
			/*look for abinitioCode*/
			CALC.abinitioCode = g_malloc(CALC._num_opt_steps*sizeof(gint));
			for(i=0;i<CALC._nsplits;i++) CALC.abinitioCode[i]=0;
			ptr=&(line[0]);i=0;
			while((*ptr!='\n')&&(*ptr!='\0')){
				if((*ptr==' ')||(*ptr=='(')) while((*ptr==' ')||(*ptr=='(')) ptr++;/*skip space(s) & parenthesis*/
				sscanf(ptr,"%i%*s",&(CALC.abinitioCode[i]));
				ptr++;i++;
				while(g_ascii_isgraph(*ptr)) ptr++;/*go to next space/end (skip ')' if any)*/
			}
			g_free(line);
			line = file_read_line(vf);/*this is the ENDabinit line*/
			g_free(line);
			line = file_read_line(vf);
			continue;
		}
		if (find_in_string("KresolStart",line) != NULL) {
			g_free(line);line = file_read_line(vf);/*go next line*/
			/*in my understanding, there is either 1 or _num_opt_steps values of KresolStart*/
			CALC.KresolStart = g_malloc(CALC._num_opt_steps*sizeof(gdouble));
			for(i=0;i<CALC._nsplits;i++) CALC.KresolStart[i]=0.;
			ptr=&(line[0]);i=0;
			while((*ptr!='\n')&&(*ptr!='\0')){
				if(*ptr==' ') while(*ptr==' ') ptr++;/*skip space(s)*/
				sscanf(ptr,"%lf%*s",&(CALC.KresolStart[i]));
				ptr++;i++;
				while(g_ascii_isgraph(*ptr)) ptr++;/*go to next space/end*/
			}
			g_free(line);
			line = file_read_line(vf);/*this is the Kresolend line*/
			g_free(line);
			line = file_read_line(vf);
			continue;
		}
		if (find_in_string("vacuumSize",line) != NULL) {
			g_free(line);line = file_read_line(vf);/*go next line*/
			/*in my understanding, there is also either 1 or _num_opt_steps values of vacuumSize*/
			CALC.vacuumSize = g_malloc(CALC._num_opt_steps*sizeof(gdouble));
			for(i=0;i<CALC._nsplits;i++) CALC.vacuumSize[i]=0.;
			ptr=&(line[0]);i=0;
			while((*ptr!='\n')&&(*ptr!='\0')){
				if(*ptr==' ') while(*ptr==' ') ptr++;/*skip space(s)*/
				sscanf(ptr,"%lf%*s",&(CALC.vacuumSize[i]));
				ptr++;i++;
				while(g_ascii_isgraph(*ptr)) ptr++;/*go to next space/end*/
			}
			g_free(line);
			line = file_read_line(vf);/*this is the endVacuumSize line*/
			g_free(line);
			line = file_read_line(vf);
			continue;
		}
		__GET_INT(numParallelCalcs);
		if (find_in_string("commandExecutable",line) != NULL) {
			g_free(line);line = file_read_line(vf);/*go next line*/
			/*there is 1<x<_num_opt_steps lines of commandExecutable*/
			j=0;
			ptr = g_strdup("");
			while((j<CALC._num_opt_steps)&&(find_in_string("EndExecutable",ptr)==NULL)){
				CALC.commandExecutable = g_strdup_printf("%s\n%s",ptr,line);
				g_free(ptr);
				g_free(line);
				line = file_read_line(vf);
				ptr = CALC.commandExecutable;
				j++;
			}
			if(find_in_string("EndExecutable",line) != NULL){
				g_free(line);
				line = file_read_line(vf);
			}
			g_free(line);
			line = file_read_line(vf);
			continue;
		}
		__GET_INT(whichCluster);
		__GET_CHARS(remoteFolder);
		__GET_BOOL(PhaseDiagram);
		__GET_DOUBLE(RmaxFing);
		__GET_DOUBLE(deltaFing);
		__GET_DOUBLE(sigmaFing);
		__GET_INT(antiSeedsActivation);
		__GET_DOUBLE(antiSeedsMax);
		__GET_DOUBLE(antiSeedsSigma);
		__GET_BOOL(doSpaceGroup);
		if (find_in_string("SymTolerance",line) != NULL) {
			/*can be either a number or quantifier*/
			if(g_ascii_isalpha(line[0])){
				/*Either high, medium, or low*/
				switch(line[0]){
				case 'H':
				case 'h':
					CALC.SymTolerance=0.05;
					break;
				case 'M':
				case 'm':
					CALC.SymTolerance=0.10;
					break;
				case 'L':
				case 'l':
					CALC.SymTolerance=0.20;
					break;
				default:/*everything else is a bad setting*/
					g_free(line);
					line = g_strdup_printf("ERROR: USPEX bad SymTolerance setting in Parameters.txt!\n");
					gui_text_show(ERROR, line);
				}
			}else{
				/*get the number directly*/
				sscanf(line,"%lf%*s",&(CALC.SymTolerance));
			}
			g_free(line);
			line = file_read_line(vf);
			continue;
		}
		__GET_INT(repeatForStatistics);
		__GET_DOUBLE(stopFitness);
		__GET_BOOL(collectForces);
		__GET_BOOL(ordering_active);
		__GET_BOOL(symmetrize);
		if (find_in_string("valenceElectr",line) != NULL) {
			g_free(line);line = file_read_line(vf);/*go next line*/
			/*if no _nspecies, get it now*/
			__NSPECIES(line);
			/*look for valenceElectr*/
			CALC.valenceElectr = g_malloc(CALC._nspecies*sizeof(gint));
			for(i=0;i<CALC._nspecies;i++) CALC.valenceElectr[i]=0;
			ptr=&(line[0]);i=0;
			while((*ptr!='\n')&&(*ptr!='\0')){
				if(*ptr==' ') while(*ptr==' ') ptr++;/*skip space(s)*/
				sscanf(ptr,"%i%*s",&(CALC.valenceElectr[i]));
				ptr++;i++;
				while(g_ascii_isgraph(*ptr)) ptr++;/*go to next space/end*/
			}
			g_free(line);
			line = file_read_line(vf);/*this is the endValenceElectr line*/
			g_free(line);
			line = file_read_line(vf);
			continue;
		}
		__GET_DOUBLE(percSliceShift);
		__GET_INT(dynamicalBestHM);
		__GET_STRING(softMutOnly);/*TODO: rewrite as "smart" int array*/
		__GET_DOUBLE(maxDistHeredity);
		__GET_INT(manyParents);
		__GET_DOUBLE(minSlice);
		__GET_DOUBLE(maxSlice);
		__GET_INT(numberparents);
		__GET_DOUBLE(thicknessS);
		__GET_DOUBLE(thicknessB);
		__GET_INT(reconstruct);
		__GET_INT(firstGeneMax);
		__GET_INT(minAt);
		__GET_INT(maxAt);
		__GET_DOUBLE(fracTrans);
		__GET_DOUBLE(howManyTrans);
		if (find_in_string("specificTrans",line) != NULL) {
			g_free(line);line = file_read_line(vf);/*go next line*/
			/*count number of specificTrans*/
			ptr=&(line[0]);
			CALC._nspetrans=1;
			while(*ptr!='\0'){
				if(*ptr==' ') {
					CALC._nspetrans++;
					while(g_ascii_isgraph(*ptr)) ptr++;
				}else ptr++;
			}
			/*look for specificTrans*/
			CALC.specificTrans = g_malloc(CALC._nspetrans*sizeof(gint));
			for(i=0;i<CALC._nspetrans;i++) CALC.specificTrans[i]=0;
			ptr=&(line[0]);i=0;
			while((*ptr!='\n')&&(*ptr!='\0')){
				if(*ptr==' ') while(*ptr==' ') ptr++;/*skip space(s)*/
				sscanf(ptr,"%i%*s",&(CALC.specificTrans[i]));
				ptr++;i++;
				while(g_ascii_isgraph(*ptr)) ptr++;/*go to next space/end*/
			}
			g_free(line);
			line = file_read_line(vf);/*this is the EndTransSpecific line*/
			g_free(line);
			line = file_read_line(vf);
			continue;
		}
		__GET_DOUBLE(GaussianWidth);
		__GET_DOUBLE(GaussianHeight);
		__GET_INT(FullRelax);
		__GET_DOUBLE(maxVectorLength);
		__GET_DOUBLE(PSO_softMut);
		__GET_DOUBLE(PSO_BestStruc);
		__GET_DOUBLE(PSO_BestEver);
		if (find_in_string("vcnebType",line) != NULL) {
			k=0;sscanf(line,"%i%*s",&(k));
			j=k;
			i=((int)(j/100))%10;
			if((i<1)||(i>2)) {
				g_free(line);
				line = file_read_line(vf);
				continue;
			}
			CALC._vcnebtype_method=i;
			j-=i*100;
			i=((int)(j/10))%10;
			if((i<0)||(i>1)){
				g_free(line);
				line = file_read_line(vf);
				continue;
			}
			CALC._vcnebtype_img_num=(i==1);
			j-=i*10;
			i=j%10;
			if((i<0)||(i>1)){
				g_free(line);
				line = file_read_line(vf);
				continue;
			}
			CALC._vcnebtype_spring=(i==1);
			CALC.vcnebType=k;
			g_free(line);
			line = file_read_line(vf);
			continue;
		}
		__GET_INT(numImages);
		__GET_INT(numSteps);
		__GET_INT(optReadImages);
		__GET_BOOL(optimizerType);
		__GET_INT(optRelaxType);
		__GET_DOUBLE(dt);
		__GET_DOUBLE(ConvThreshold);
		__GET_DOUBLE(VarPathLength);
		__GET_DOUBLE(K_min);
		__GET_DOUBLE(K_max);
		__GET_DOUBLE(Kconstant);
		__GET_BOOL(optFreezing);
		__GET_INT(optMethodCIDI);
		__GET_INT(startCIDIStep);
		if (find_in_string("pickupImages",line) != NULL) {
			g_free(line);line = file_read_line(vf);/*go next line*/
			/*count number of picked-up images*/
			ptr=&(line[0]);
			CALC._npickimg=1;
			while(*ptr!='\0'){
				if(*ptr==' ') {
					CALC._npickimg++;
					while(g_ascii_isgraph(*ptr)) ptr++;
				}else ptr++;
			}
			/*look for pickupImages*/
			CALC.pickupImages = g_malloc(CALC._npickimg*sizeof(gint));
			for(i=0;i<CALC._npickimg;i++) CALC.pickupImages[i]=0;
			ptr=&(line[0]);i=0;
			while((*ptr!='\n')&&(*ptr!='\0')){
				if(*ptr==' ') while(*ptr==' ') ptr++;/*skip space(s)*/
				sscanf(ptr,"%i%*s",&(CALC.pickupImages[i]));
				ptr++;i++;
				while(g_ascii_isgraph(*ptr)) ptr++;/*go to next space/end*/
			}
			g_free(line);
			line = file_read_line(vf);/*this is the EndPickupImages line*/
			g_free(line);
			line = file_read_line(vf);
			continue;
		}
		__GET_INT(FormatType);
		__GET_INT(PrintStep);
	/*in case no match was done:*/
	g_free(line);
	line = file_read_line(vf);
	}

	return 0;
#undef __NSPECIES
#undef CALC
}

gint read_individuals_uspex(gchar *filename, struct model_pak *model){
        FILE *vf;
	long int vfpos;
        gchar *line=NULL;
	gchar tmp[64];
        uspex_output_struct *uspex_output=model->uspex;
        /* specific */
	gint idx=0;
	gchar *ptr;
	gint max_struct=0;
	gint red_index;
	/*since energy is sometimes eV, and sometimes eV/atom ?*/
	gboolean e_red=FALSE;
	/*start*/
        vf = fopen(filename, "rt");
        if (!vf) return 1;
	/*Have fitness? Better watch for each case...*/
	/*^^^^ I thought that calculationMethod = USPEX => fitness, but got defeated by EX19*/
	line = file_read_line(vf);
	if (find_in_string("Fitness",line) != NULL) _UO.have_fitness=TRUE;
	else _UO.have_fitness=FALSE;
	rewind(vf);
	/* deal with the problematic units */
	if(fetch_in_file(vf,"eV/atom")!=0) e_red=TRUE;
	rewind(vf);
        /* Skip header: FIX a _BUG_ when Individuals file sometimes has 2 line and sometimes only 1 */
        vfpos=ftell(vf);/* flag */
        line = file_read_line(vf);
        while(line){
                ptr=&(line[0]);
                while(*ptr==' ') ptr++;
                if(g_ascii_isdigit(*ptr)) {
                        fseek(vf,vfpos,SEEK_SET);/* rewind to flag */
                        break;
                }
                vfpos=ftell(vf);/* flag */
                g_free(line);
                line = file_read_line(vf);
        }
        if(line==NULL) return -1;
	/* META calculation have a gen=0 Individual
	  unfortunately, it is *not* in the gather-
	 -POSCARS file, so we ignore gen=0 for now. */
	if(_UO.method==US_CM_META){
		g_free(line);
		line = file_read_line(vf);
		vfpos=ftell(vf);/* flag */
	}
	/*do a 1st pass to count stuctures*/
	_UO.num_struct=0;
        while (line){
		_UO.num_struct++;
		sscanf(line," %*i %i %*s",&idx);
		if(idx>max_struct) max_struct=idx;
		g_free(line);
		line = file_read_line(vf);
	}
	_UO.num_struct--;
	fseek(vf,vfpos,SEEK_SET);/* rewind to flag */
	line = file_read_line(vf);
	/* now prepare arrays <- JOB=USPEX/300 example */
	if(_UO.num_struct==0) return 1;
	_UO.red_index=g_malloc(max_struct*sizeof(gint));
	for(idx=0;idx<max_struct;idx++) _UO.red_index[idx]=0;
	sscanf(line," %*[^[][%[^]]] %lf %*f %*s",
		&(tmp[0]),&(_UO.min_E));
	if(!e_red) _UO.min_E/=(gdouble)_UO.ind[0].natoms;
	_UO.max_E=_UO.min_E;
	idx=0;_UO.num_gen=1;
	while (line){
if(_UO.have_fitness){
	sscanf(line," %i %i%*[^[][%[^]]] %lf %lf %*f %lf %*s",
		&(_UO.ind[idx].gen),&(red_index),&(tmp[0]),&(_UO.ind[idx].energy),&(_UO.ind[idx].volume),&(_UO.ind[idx].fitness));
	if(isnan(_UO.ind[idx].fitness)) _UO.ind[idx].fitness=10000;
}else
	sscanf(line," %i %i%*[^[][%[^]]] %lf %lf %*s",
		&(_UO.ind[idx].gen),&(red_index),&(tmp[0]),&(_UO.ind[idx].energy),&(_UO.ind[idx].volume));
		if(_UO.ind[idx].gen<_UO.num_gen) _UO.ind[idx].gen++;
		if(_UO.ind[idx].gen>_UO.num_gen) _UO.num_gen=_UO.ind[idx].gen;
		_UO.red_index[red_index]=idx;
		if(!e_red) _UO.ind[idx].E=_UO.ind[idx].energy/_UO.ind[idx].natoms;
		else _UO.ind[idx].E=_UO.ind[idx].energy;
#if DEBUG_USPEX_READ
fprintf(stdout,"#DBG: USPEX[Individual]: GEN=%i STRUCT=%i e=%lf n_atom=%i E=%lf\n",_UO.ind[idx].gen,idx+1,_UO.ind[idx].energy,_UO.ind[idx].natoms,_UO.ind[idx].E);
#endif
		if(_UO.ind[idx].E<_UO.min_E) _UO.min_E=_UO.ind[idx].E;
		if(_UO.ind[idx].E>_UO.max_E) _UO.max_E=_UO.ind[idx].E;
		g_free(line);
		line = file_read_line(vf);
		idx++;
	}
	/*all done*/
	return 0;
}

gint read_output_uspex(gchar *filename, struct model_pak *model){
	gchar *line;
	FILE *vf;
	long int vfpos;
	gchar *ptr;
	gchar *ptr2;
        gchar *res_folder;
	gchar *aux_file;
	gint ix;
	gint idx;
	gint jdx;
	gint min,max,med;
	gint job;
	/* results */
	FILE *f_src;
	FILE *f_dest;
	gdouble *e;
	gdouble *c;
	gdouble *tag;
	gdouble *tmp;
	gdouble compo;
	gint n_compo;
	gint skip=0;
	gint gen=1;
	gint num;
	gint natoms;
	gint species_index;
	gdouble min_E;
	gdouble min_F=0.;
	gdouble max_E;
	gchar *atoms;
	uspex_output_struct *uspex_output;
	/* checks */
	g_return_val_if_fail(model != NULL, 1);
	g_return_val_if_fail(filename != NULL, 2);
	if (find_in_string("OUTPUT.txt",filename) == NULL){
		line = g_strdup_printf("ERROR: Please load USPEX OUTPUT.txt file!\n");
		gui_text_show(ERROR, line);
		g_free(line);
		return -1;
	}
	/*allocs*/
	if(model->uspex!=NULL) g_free(model->uspex);
	model->uspex=g_malloc(sizeof(uspex_output_struct));
	uspex_output=model->uspex;
	/* TODO: check that selected file is OUTPUT.txt */
	vf = fopen(filename, "rt");
	if (!vf) return 1;
	error_table_clear();
/* --- get result folder*/
	res_folder=g_strdup (filename);
	ptr=g_strrstr (res_folder,"OUTPUT.txt");
	*ptr='\0';/*cut file at folder*/
#if DEBUG_USPEX_READ
	fprintf(stdout,"#DBG: USPEX_FOLDER=%s\n",res_folder);
#endif
/* --- setup environment */
	sysenv.render.show_energy = TRUE;
	
/* --- read the system OUTPUT.txt (will be confirmed by Parameters.txt)*/
	if(fetch_in_file(vf,"Version")==0){
		/*can't get version number (10.0.0 example will fail here)*/
		/* -> just produce a warning... file might actually be processed correctly*/
		line = g_strdup_printf("WARNING: USPEX version unsupported!\n");
		gui_text_show(WARNING, line);
		g_free(line);
	}else{
		rewind(vf);
		line = file_read_line(vf);
		while(line) {
			if (find_in_string("Version",line) != NULL){
				sscanf(line,"| Version %i.%i.%i %*s",&(max),&(med),&(min));
				_UO.version=min+10*med+100*max;
				if(_UO.version!=944){
					g_free(line);
					line = g_strdup_printf("WARNING: USPEX version %i unsupported!\n",_UO.version);
					gui_text_show(WARNING, line);

				}
				g_free(line);
				break;
			}
			g_free(line);
			line = file_read_line(vf);
		}
	}
	rewind(vf);
	job=0;
	if(fetch_in_file(vf,"Block for system description")==0) {
/*since VCNEB will fail here due to a non-unified OUTPUT.txt format
 *we need to be a little more permissive... */ 
		rewind(vf);
		if(fetch_in_file(vf,"VCNEB")==0) goto uspex_fail;/*definitely wrong*/
		job=-1;
	}
if(job>-1){
	/*next line contains Dimensionality*/
	ix=0;
	line = file_read_line(vf);
	/*unless from a previous version and it contains a lot of '-' character*/
	if (find_in_string("Dimensionality",line) == NULL) {
		g_free(line);
		line = file_read_line(vf);
	}
	sscanf(line," Dimensionality : %i ",&(ix));
	switch (ix){/*to be use later*/
		case -2:
		/*we do not deal with 2D crystal yet*/
		case 0:
		case 1:
		case 2:
		case 3:
		_UO.dim=ix;
		break;
		default:
		line = g_strdup_printf("ERROR: reading USPEX OUTPUT.txt file!\n");
		gui_text_show(ERROR, line);
		g_free(line);
		goto uspex_fail;
	}
	job=ix*100;
	g_free(line);
	/*next line contains Molecular*/
	line = file_read_line(vf);
	sscanf(line," Molecular : %i %*s",&(ix));
	job+=(ix*10);
	if(ix==0) _UO.mol=FALSE;
	else if(ix==1) _UO.mol=TRUE;
	else{
		line = g_strdup_printf("ERROR: reading USPEX OUTPUT.txt file!\n");
		gui_text_show(ERROR, line);
		g_free(line);
		goto uspex_fail;
	}
	g_free(line);
	/*next line contains Variable Composition*/
	line = file_read_line(vf);
	sscanf(line," Variable Composition : %i %*s",&(ix));
	job+=ix;
	if(ix==0) _UO.var=FALSE;
	else if(ix==1) _UO.var=TRUE;
	else{
		line = g_strdup_printf("ERROR: reading USPEX OUTPUT.txt file!\n");
		gui_text_show(ERROR, line);
		g_free(line);
		goto uspex_fail;
	}
#if DEBUG_USPEX_READ
	fprintf(stdout,"#DBG: USPEX OUTPUT calculationType=%i\n",job);
#endif
	g_free(line);
}/*in case of VCNEB, we don't have a correct job number*/
	line = file_read_line(vf);
	num=0;
	while(line){
		if(find_in_string("types of atoms in the system",line) != NULL) {
			/**/
			sscanf(line," There are %i %*s",&(num));
		}
		g_free(line);
		line = file_read_line(vf);
	}
	g_free(line);
/* --- and close */
	fclose(vf);
/* --- READ Parameters.txt */
	aux_file = g_strdup_printf("%s%s",res_folder,"Parameters.txt");
	if(read_parameter_uspex(aux_file,model)!=0) {
		line = g_strdup_printf("ERROR: reading USPEX Parameter.txt file!\n");
		gui_text_show(ERROR, line);
		g_free(line);
		goto uspex_fail;
	}
/* --- CHECK type from Parameters.txt matches that of OUTPUT.TXT */
	if((job!=-1)&&(job!=_UO.type)){
/*since VCNEB will fail here due to a non-unified OUTPUT.txt format
 *we need to be a little more permissive... */
		line = g_strdup_printf("ERROR: inconsistent USPEX files (%i!=%i)!\n",job,_UO.type);
		gui_text_show(ERROR, line);
		g_free(line);
		goto uspex_fail;
	}
	g_free(aux_file);
	/*special case: no nspecies information*/
	if(_UO.nspecies==0) _UO.nspecies=num;
if((_UO.method==US_CM_USPEX)||(_UO.method==US_CM_META)){
/* --- if calculationMethod={USPEX,META} */
	model->basename=g_strdup_printf("uspex");
/* we have to open structure file first because of inconsistency in atom reporting in Individuals file*/
/* --- READ gatheredPOSCARS <- this is going to be the main model file */
/* ^^^ META: If gatheredPOSCARS_relaxed exists, read it instead */
        if(_UO.method==US_CM_META) {
                aux_file = g_strdup_printf("%s%s",res_folder,"gatheredPOSCARS_relaxed");
                vf = fopen(aux_file, "rt");
                if (!vf) aux_file = g_strdup_printf("%s%s",res_folder,"gatheredPOSCARS");
                else fclose(vf);
        }else{  
                aux_file = g_strdup_printf("%s%s",res_folder,"gatheredPOSCARS");
        }
        vf = fopen(aux_file, "rt");
        if (!vf) {
                if(_UO.ind!=NULL) g_free(_UO.ind);
                line = g_strdup_printf("ERROR: can't open USPEX gatheredPOSCARS file!\n");
                gui_text_show(ERROR, line);
                g_free(line);
                goto uspex_fail;
        }
	/*count number of structures*/
	_UO.num_struct=0;
	line = file_read_line(vf);
	while(!feof(vf)){
		if((line[0]=='E')&&(line[1]=='A')) _UO.num_struct++;
		g_free(line);
		line = file_read_line(vf);
	}
	rewind(vf);
        model->num_frames=0;
        vfpos=ftell(vf);/* flag */
        line = file_read_line(vf);
	_UO.ind=g_malloc((_UO.num_struct)*sizeof(uspex_individual));
	idx=0;ix=-1;
        while(!feof(vf)){
                if((line[0]=='E')&&(line[1]=='A')) {
			idx=0;ix++;
                        fseek(vf,vfpos,SEEK_SET);/* rewind to flag */
                        add_frame_offset(vf, model);
                        model->num_frames++;
                        g_free(line);
                        line = file_read_line(vf);
                }
		if(idx==5){/*calculate number of species and atoms*/
			_UO.ind[ix].atoms=g_malloc(_UO.nspecies*sizeof(gint));
			for(jdx=0;jdx<_UO.nspecies;jdx++) _UO.ind[ix].atoms[jdx]=0;/*init atoms*/
			_UO.ind[ix].natoms=0;
			/*get this structure nspecies*/
			ptr=&(line[0]);
			while(*ptr==' ') ptr++;/*skip initial space (if any)*/
			jdx=1;/*there is at least one species*/
			while(*ptr!='\0'){
				if(*ptr==' ') {
					jdx++;ptr++;
					while(*ptr==' ') ptr++;/*skip space*/
				}else ptr++;
			}
			if(jdx<_UO.nspecies){
				gchar *line2;/*double buffer*/
				/*we need to know which species is here*/
				line2 = file_read_line(vf);
				ptr=&(line[0]);
				ptr2=&(line2[0]);
				idx++;
				jdx=0;
				while(*ptr==' ') ptr++;/*skip initial space (if any)*/
				while(*ptr2==' ') ptr2++;/*skip initial space (if any)*/
				while((*ptr!='\0')&&(*ptr2!='\0')){
					/*get the correct jdx number*/
					jdx=0;/*we reset jdx to cope with 'OUT OF ORDER' poscar species, if any*/
					while((jdx<_UO.nspecies)&&(_UO.spe_Z[jdx]!=elem_symbol_test(ptr))) jdx++;
					_UO.ind[ix].atoms[jdx]=g_ascii_strtod(ptr2,NULL);
					_UO.ind[ix].natoms+=_UO.ind[ix].atoms[jdx];
					ptr++;ptr2++;
					while(g_ascii_isgraph(*ptr)) ptr++;/*go to next space/end*/
					while(g_ascii_isgraph(*ptr2)) ptr2++;/*go to next space/end*/
					ptr++;ptr2++;
					while(*ptr==' ') ptr++;/*skip space*/
					while(*ptr2==' ') ptr2++;/*skip space*/
				}
				g_free(line);
				line=line2;
			}else{
				/*get a new line (the number of each atoms)*/
				g_free(line);
				line = file_read_line(vf);
				idx++;
				ptr=&(line[0]);
				ptr2=ptr;jdx=0;
				do{/*get the number of each species (in the Parameters.txt order)*/
					_UO.ind[ix].atoms[jdx]=g_ascii_strtod(ptr,&ptr2);
					_UO.ind[ix].natoms+=_UO.ind[ix].atoms[jdx];
					if(ptr2==ptr) break;
					jdx++;
					ptr=ptr2;
				}while(1);
			}
		}
                vfpos=ftell(vf);/* flag */
		idx++;
                g_free(line);
                line = file_read_line(vf);
        }
        strcpy(model->filename,aux_file);// which means that we "forget" about OUTPUT.txt
        g_free(aux_file);

/* --- READ Individuals <- information about all structures */
/* ^^^ META: If Individuals_relaxed exists, read it instead */
	if(_UO.method==US_CM_META) {
		aux_file = g_strdup_printf("%s%s",res_folder,"Individuals_relaxed");
		vf = fopen(aux_file, "rt");
		if (!vf) aux_file = g_strdup_printf("%s%s",res_folder,"Individuals");
		else fclose(vf);
	}else{
		aux_file = g_strdup_printf("%s%s",res_folder,"Individuals");
	}
	if(read_individuals_uspex(aux_file,model)!=0) {
		line = g_strdup_printf("ERROR: reading USPEX Individuals file!\n");
		gui_text_show(ERROR, line);
		g_free(line);
		goto uspex_fail;
	}
	g_free(aux_file);
/* --- open the last frame */
	if(_UO.method==US_CM_META) {
		aux_file = g_strdup_printf("%s%s",res_folder,"gatheredPOSCARS_relaxed");
		vf = fopen(aux_file, "rt");
		if (!vf) aux_file = g_strdup_printf("%s%s",res_folder,"gatheredPOSCARS");
		else fclose(vf);
	}else{
		aux_file = g_strdup_printf("%s%s",res_folder,"gatheredPOSCARS");
	}
        vf = fopen(aux_file, "rt");
        if (!vf) {/*very unlikely, we just opened it*/
                if(_UO.ind!=NULL) g_free(_UO.ind);
                line = g_strdup_printf("ERROR: can't open USPEX gatheredPOSCARS file!\n");
                gui_text_show(ERROR, line);
                g_free(line);
                goto uspex_fail;
        }
	read_frame_uspex(vf,model);/*open first frame*/
        fclose(vf);

/* --- Prepare the summary graph */
	/*no need to load BESTIndividuals since we already know everything from Individuals*/
	_UO.graph=graph_new("ALL", model);
	_UO.graph_best=graph_new("BEST", model);
	if(_UO.num_struct>_UO.num_gen){
		skip=0;
		gen=1;
		_UO.num_best=_UO.num_gen;
		_UO.best_ind=g_malloc((1+_UO.num_best)*sizeof(gint));
		_UO.best_ind[0]=0;
		ix=0;
		while(gen<=_UO.num_gen){
		  while(_UO.ind[ix].gen<gen) ix++;/*reach generation gen*/
		  num=0;
		  skip=ix;
		  min_E=_UO.ind[ix].E;
		  if(_UO.have_fitness) min_F=_UO.ind[ix].fitness;
		  _UO.best_ind[gen]=ix;
		  while(_UO.ind[ix].gen==gen){
			if(_UO.have_fitness){
				if(_UO.ind[ix].fitness<min_F){
					min_F=_UO.ind[ix].fitness;
					_UO.best_ind[gen]=ix;
				}
			}else{
				if(_UO.ind[ix].E<min_E){
					min_E=_UO.ind[ix].E;
					_UO.best_ind[gen]=ix;
				}
			}
			num++;
			ix++;
		  }
		  e=g_malloc((1+num)*sizeof(gdouble));
		  ix=skip;
		  e[0]=(gdouble)num;
		  while(_UO.ind[ix].gen==gen){
			e[ix-skip+1]=_UO.ind[ix].E;
			ix++;
		  }
		  graph_add_borned_data(
			_UO.num_gen,e,0,_UO.num_gen,
			_UO.min_E-(_UO.max_E-_UO.min_E)*0.05,_UO.max_E+(_UO.max_E-_UO.min_E)*0.05,GRAPH_USPEX,_UO.graph);
		  gen++;
		  g_free(e);
		}
		/*prepare array*/
		e=g_malloc((1+_UO.num_gen)*sizeof(gdouble));
		e[0]=(gdouble)_UO.num_gen;
		max_E=_UO.ind[_UO.best_ind[1]].E;
		min_E=max_E;
		for(gen=0;gen<_UO.num_gen;gen++) {
			if(_UO.ind[_UO.best_ind[gen+1]].E>max_E) max_E=_UO.ind[_UO.best_ind[gen+1]].E;
			if(_UO.ind[_UO.best_ind[gen+1]].E<min_E) min_E=_UO.ind[_UO.best_ind[gen+1]].E;
			e[gen+1]=_UO.ind[_UO.best_ind[gen+1]].E;
#if DEBUG_USPEX_READ
fprintf(stdout,"#DBG: USPEX_BEST: GEN=%i STRUCT=%i E=%lf\n",gen+1,1+_UO.best_ind[gen+1],e[gen+1]);
#endif
		}
		if(_UO.method==US_CM_META) {/*shift best_ind by 1*/
			for(gen=1;gen<=_UO.num_gen;gen++) _UO.best_ind[gen]++;
		}
		graph_add_borned_data(
			_UO.num_gen,e,0,_UO.num_gen,
			min_E-(max_E-min_E)*0.05,max_E+(max_E-min_E)*0.05,GRAPH_USPEX_BEST,_UO.graph_best);
		g_free(e);
	}
	/*set ticks*/
	if(_UO.num_gen>15) ix=5;
	else ix=_UO.num_gen+1;
	graph_set_xticks(TRUE,ix,_UO.graph);
	graph_set_yticks(TRUE,5,_UO.graph);
	graph_set_xticks(TRUE,ix,_UO.graph_best);
	graph_set_yticks(TRUE,5,_UO.graph_best);
/* --- variable composition */
if(_UO.var){
/* --- NEW: add a simple composition % graph for each species */
if(_UO.nspecies>1){
n_compo=2;
ix=0;
for(species_index=0;species_index<_UO.nspecies;species_index++){
	c=g_malloc(2*sizeof(gdouble));/*suppose at least 2 different compositions*/
	c[0]=(gdouble)_UO.ind[0].atoms[species_index] / (gdouble)_UO.ind[0].natoms;
	c[1]=c[0];
	compo=c[0];
	ix=0;
	while((compo==c[0])&&(ix<_UO.num_struct)){
		compo=(gdouble)_UO.ind[ix].atoms[species_index] / (gdouble)_UO.ind[ix].natoms;
		if(compo!=c[0]){
			if(compo>c[0]) c[1]=compo;
			else {
				c[1]=c[0];c[0]=compo;
			}
		}
		ix++;
	}
	if(c[1]==c[0]) continue;/*in VARCOMP calculation, one species can be kept fixed*/
	for(ix=0;ix<_UO.num_struct;ix++){
		compo=(gdouble)_UO.ind[ix].atoms[species_index] / (gdouble)_UO.ind[ix].natoms;
		if(compo<c[0]) {/*we have a new first element*/
			tmp=g_malloc((n_compo+1)*sizeof(gdouble));
				tmp=g_malloc((n_compo+1)*sizeof(gdouble));
				tmp[0]=compo;
				memcpy(&(tmp[1]),c,n_compo*sizeof(gdouble));
				g_free(c);
				c=tmp;
				n_compo++;
				continue;
			}
			if(compo>c[n_compo-1]){/*we have a new last element*/
				tmp=g_malloc((n_compo+1)*sizeof(gdouble));
				tmp[n_compo]=compo;
				memcpy(&(tmp[0]),c,n_compo*sizeof(gdouble));
				g_free(c);
				c=tmp;
				n_compo++;
				continue;
			}
			for(jdx=1;jdx<n_compo;jdx++){
				if((compo>c[jdx-1])&&(compo<c[jdx])){
				/*insert a new element at position jdx*/
					tmp=g_malloc((n_compo+1)*sizeof(gdouble));
					tmp[jdx]=compo;
					memcpy(&(tmp[0]),c,(jdx)*sizeof(gdouble));
					memcpy(&(tmp[jdx+1]),&(c[jdx]),(n_compo-jdx)*sizeof(gdouble));
					g_free(c);
					c=tmp;
					n_compo++;
					jdx=n_compo;/*exit loop*/
					continue;
				}
				
			}
		}
		/*build the energy array*/
		e=g_malloc(n_compo*sizeof(gdouble));
		tag=g_malloc(n_compo*sizeof(gdouble));
		max_E=_UO.min_E;
		for(ix=0;ix<n_compo;ix++) e[ix]=_UO.max_E;
		for(ix=0;ix<_UO.num_struct;ix++){
			for(jdx=0;jdx<n_compo;jdx++){
				compo=(gdouble)_UO.ind[ix].atoms[species_index] / (gdouble)_UO.ind[ix].natoms;
				if(c[jdx]==compo){
					if(_UO.ind[ix].E < e[jdx]) {
						e[jdx] = _UO.ind[ix].E;
						tag[jdx] = (gdouble) ix;
						if(e[jdx]>max_E) max_E=e[jdx];
					}
				}
			}
		}
		/* prepare composition diagram*/
		line=g_strdup_printf("COMP_%s",elements[_UO.spe_Z[species_index]].symbol);
		_UO.graph_comp=graph_new(line, model);
                graph_add_borned_data(
			n_compo,tag,0,1,_UO.min_E-(max_E-_UO.min_E)*0.05,max_E+(max_E-_UO.min_E)*0.05,GRAPH_USPEX_2D,_UO.graph_comp);
		graph_add_borned_data(
			n_compo,c,0,1,_UO.min_E-(max_E-_UO.min_E)*0.05,max_E+(max_E-_UO.min_E)*0.05,GRAPH_USPEX_2D,_UO.graph_comp);
		graph_add_borned_data(
			n_compo,e,0,1,_UO.min_E-(max_E-_UO.min_E)*0.05,max_E+(max_E-_UO.min_E)*0.05,GRAPH_USPEX_2D,_UO.graph_comp);
                g_free(line);
		g_free(tag);
		g_free(c);
		g_free(e);
        	/*set ticks*/
        	graph_set_xticks(TRUE,3,_UO.graph_comp);
        	graph_set_yticks(TRUE,5,_UO.graph_comp);

	}

}/*variable composition with only 1 species...*/
}/*not a VARCOMP calculation*/

	/*refresh model*/
	tree_model_refresh(model);
        model->redraw = TRUE;
        /* always show this information */
        gui_text_show(ITALIC,g_strdup_printf("USPEX: %i structures detected.\n",model->num_frames));
        model_prep(model);


}else if(_UO.method==US_CM_VCNEB){
/* --- tentative at adding VCNEB, which format is different**/
/* ^^^ *BUT unfortunately, we NEED a vasp5 output (VCNEB stayed at VASP4) */
/*      thus we create a new file called transitionPath_POSCARs5 which is a
 *      translation of transitionPath_POSCARs for VASP5...*/
	atoms = g_strdup_printf("%s",elements[_UO.spe_Z[0]].symbol);
	for(idx=1;idx<_UO.nspecies;idx++){
		line = g_strdup_printf("%s %s",atoms,elements[_UO.spe_Z[idx]].symbol);
		g_free(atoms);
		atoms=line;
	}
/* open transitionPath_POSCARs and start to convert to transitionPath_POSCARs5 */
	aux_file = g_strdup_printf("%s%s",res_folder,"transitionPath_POSCARs");
	f_src = fopen(aux_file, "rt");
	if(!f_src) {
		line = g_strdup_printf("ERROR: can't open USPEX transitionPath_POSCARs file!\n");
		gui_text_show(ERROR, line);
		g_free(line);
		goto uspex_fail;
	}
	g_free(aux_file);
	aux_file = g_strdup_printf("%s%s",res_folder,"transitionPath_POSCARs5");
	f_dest = fopen(aux_file, "w");
	if(!f_dest) {
		vf=f_src;
		line = g_strdup_printf("ERROR: can't WRITE into USPEX transitionPath_POSCARs5 file!\n");
		gui_text_show(ERROR, line);
		g_free(line);
		goto uspex_fail;
	}
/*get the total number of atoms*/
	for(idx=0;idx<5;idx++){
		line = file_read_line(f_src);
		g_free(line);
	}
	natoms=1;/*there is at least one atom*/
	ptr=&(line[0]);
	do{
		natoms+=g_ascii_strtod(ptr,&ptr2);
		if(ptr2==ptr) break;
		ptr=ptr2;
	}while(1);
#if DEBUG_USPEX_READ
fprintf(stdout,"#DBG: USPEX[VCNEB]: NATOMS=%i ATOMS=%s\n",natoms,atoms);
#endif
	/*translate input VASP4 into proper VASP5*/
	rewind(f_src);
	idx=0;
	_UO.num_struct=0;
	line = file_read_line(f_src);
	while(line){
		fprintf(f_dest,"%s",line);
		if (find_in_string("Image",line) != NULL) {
			idx=0;
			_UO.num_struct++;
		}
		if(idx==4) fprintf(f_dest,"%s\n",atoms);
		idx++;
		g_free(line);
		line = file_read_line(f_src);
	}
	fclose(f_src);
	fclose(f_dest);
	fclose(vf);
	g_free(atoms);
	_UO.ind=g_malloc((_UO.num_struct)*sizeof(uspex_individual));
/* --- read energies from the Energy file */
	g_free(aux_file);
	aux_file = g_strdup_printf("%s%s",res_folder,"Energy");
	vf=fopen(aux_file,"rt");
	if (!vf) {
		g_free(_UO.ind);
		line = g_strdup_printf("ERROR: can't open USPEX Energy file!\n");
		gui_text_show(ERROR, line);
		g_free(line);
		goto uspex_fail;
	}
	/*first get to the first data line*/
	line = file_read_line(vf);
	ptr=&(line[0]);
	while((ptr)&&(*ptr==' ')) ptr++;/*skip blanks*/
	while(!g_ascii_isdigit(*ptr)){
		g_free(line);
		line = file_read_line(vf);
		ptr=&(line[0]);
		while((ptr)&&(*ptr==' ')) ptr++;/*skip blanks*/
	}/*line should now point to the first data line*/
	/*get first energy*/
	_UO.best_ind=g_malloc((1+_UO.num_struct)*sizeof(gint));
	e = g_malloc((_UO.num_struct+1)*sizeof(gdouble));
	e[0]=(gdouble) _UO.num_struct;/*first element is the size*/
	sscanf(line," %*i %lf %*s",&(e[1]));
	e[1] /= (gdouble)natoms;
	_UO.min_E=e[1];
	_UO.max_E=e[1];
	for(idx=0;idx<_UO.num_struct;idx++){
		if(!line) break;/*<- useful? */
		sscanf(line," %*i %lf %*s",&(_UO.ind[idx].energy));
		_UO.ind[idx].atoms=g_malloc(_UO.nspecies*sizeof(gint));
		_UO.ind[idx].natoms=natoms;/*<-constant*/
		_UO.ind[idx].E = _UO.ind[idx].energy / (gdouble)_UO.ind[idx].natoms;
		_UO.ind[idx].gen=idx;/*<- this is actually *not* true*/
		_UO.best_ind[idx+1]=idx;
		e[idx+1]=_UO.ind[idx].E;
		if(e[idx+1]<_UO.min_E) _UO.min_E=e[idx+1];
		if(e[idx+1]>_UO.max_E) _UO.max_E=e[idx+1];
#if DEBUG_USPEX_READ
fprintf(stdout,"#DBG: USPEX[VCNEB] Image %i: natoms=%i energy=%lf e=%lf\n",_UO.ind[idx].gen,_UO.ind[idx].natoms,_UO.ind[idx].energy,e[idx]);
#endif
		g_free(line);
		line = file_read_line(vf);
	}
	fclose(vf);
	g_free(aux_file);
	/*create the reduced index table*/
	_UO.red_index=g_malloc(_UO.num_struct*sizeof(gint));
	for(idx=0;idx<_UO.num_struct;idx++) _UO.red_index[idx]=idx;
	aux_file = g_strdup_printf("%s%s",res_folder,"transitionPath_POSCARs5");
/* --- reopen the newly created transitionPath_POSCARs5 file for reading <- this will be our new model file*/
	vf=fopen(aux_file,"rt");
        if (!vf) {/*very unlikely: we just create it*/
		g_free(_UO.ind);
		line = g_strdup_printf("ERROR: can't open USPEX transitionPath_POSCARs5 file!\n");
                gui_text_show(ERROR, line);
                g_free(line);
                goto uspex_fail;
        }
        model->num_frames=0;
        vfpos=ftell(vf);/* flag */
        line = file_read_line(vf);
	idx=0;
        while(!feof(vf)){
		if (find_in_string("Image",line) != NULL) {
			idx=0;
                        fseek(vf,vfpos,SEEK_SET);/* rewind to flag */
                        add_frame_offset(vf, model);
                        model->num_frames++;
                        g_free(line);
                        line = file_read_line(vf);
		}
                vfpos=ftell(vf);/* flag */
                g_free(line);
                line = file_read_line(vf);
		idx++;
        }
        strcpy(model->filename,aux_file);// which means that we "forget" about OUTPUT.txt
        g_free(aux_file);
        rewind(vf);
        read_frame_uspex(vf,model);/*open first frame*/
        fclose(vf);

	/*do a "best" graph with each image*/
        _UO.graph_best=graph_new("PATH", model);
	_UO.num_gen=_UO.num_struct;
        graph_add_borned_data(
                _UO.num_gen,e,0,_UO.num_gen,
                _UO.min_E-(_UO.max_E-_UO.min_E)*0.05,_UO.max_E+(_UO.max_E-_UO.min_E)*0.05,GRAPH_USPEX_BEST,_UO.graph_best);
        g_free(e);
        /*set ticks*/
        if(_UO.num_gen>15) ix=5;
        else ix=_UO.num_gen+1;
        graph_set_xticks(TRUE,ix,_UO.graph_best);
        graph_set_yticks(TRUE,5,_UO.graph_best);


        /*refresh model*/
        tree_model_refresh(model);
        model->redraw = TRUE;
        /* always show this information */
        gui_text_show(ITALIC,g_strdup_printf("USPEX: %i structures detected.\n",model->num_frames));
        model_prep(model);



}else{/*calculation method is not supported*/
	line = g_strdup_printf("ERROR: USPEX support is limited to calculationMethod={USPEX,META,VCNEB} for now!\n");
	gui_text_show(ERROR, line);
	g_free(line);
	goto uspex_fail;
}


	error_table_print_all();
	return 0;
uspex_fail:
	line = g_strdup_printf("ERROR: loading USPEX failed!\n");
	gui_text_show(ERROR, line);
	g_free(line);
	fclose(vf);
	return 3;
}
gint read_frame_uspex(FILE *vf, struct model_pak *model){
	gchar *line;
	uspex_output_struct *uspex_output=model->uspex;
	long int vfpos;/*to counter _BUG_ where loading a raw frame does not update model->cur_frame*/
	gchar *ptr;
	gint idx=0;
	/*read frame*/
        g_assert(vf != NULL);

vfpos=ftell(vf);/* flag */
/*_BUG_ model->cur_frame not updated, use frame title instead*/
line = file_read_line(vf);
if(line==NULL) return -1;
ptr=&(line[0]);
while(g_ascii_isalpha(*ptr)) ptr++;
sscanf(ptr,"%i%*s",&idx);
fseek(vf,vfpos,SEEK_SET);/* rewind to flag */
        if(vasp_load_poscar5(vf,model)<0) return 3;
/* in case of meta idx is *not* the real number of the structure */
/* now the index should always be _UO.red_index[idx-1] <- I THINK*/
	line=g_strdup_printf("%lf eV",_UO.ind[_UO.red_index[idx]].energy);
	property_add_ranked(3, "Energy", line, model);
	model->volume=_UO.ind[_UO.red_index[idx]].volume;
	g_free(line);
	line=g_strdup_printf("%i",idx);
	property_add_ranked(8, "Structure", line, model);
	g_free(line);
	return 0;
}

#undef _UO

