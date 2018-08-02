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
#define _UC (*uspex_calc)

gint read_parameter_uspex(gchar *filename, struct model_pak *model){
	FILE *vf;
        gchar *line=NULL;
	gchar *ptr;
	gint idx;
        uspex_output_struct *uspex_calc=model->uspex;
	/* specific */
	gchar method;
	/*start*/
	_UC.nspecies=0;
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
				_UC.method=US_CM_USPEX;
				property_add_ranked(2, "Calculation", "USPEX", model);
#if DEBUG_USPEX_READ
fprintf(stdout,"#DBG: USPEX calculationMethod=USPEX\n");
#endif
				break;
                                case 'M':
                                case 'm':
                                /*metadynamics*/
				_UC.method=US_CM_META;
#if DEBUG_USPEX_READ
fprintf(stdout,"#DBG: USPEX calculationMethod=META\n");
#endif
				break;
                                case 'V':
                                case 'v':
                                /*variable cell nudged elastic band*/
				_UC.method=US_CM_VCNEB;
#if DEBUG_USPEX_READ
fprintf(stdout,"#DBG: USPEX calculationMethod=VCNEB\n");
#endif
				break;
                                case 'P':
                                case 'p':
                                /*particle swarm optimization*/
				_UC.method=US_CM_PSO;
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
                        sscanf(line,"%i : calculationType %*s",&(_UC.type));
#if DEBUG_USPEX_READ
fprintf(stdout,"#DBG: USPEX calculationType=%i\n",_UC.type);
#endif
                }
		if (find_in_string("optType",line) != NULL) {
			sscanf(line,"%i : optType %*s",&(_UC.opt_type));
#if DEBUG_USPEX_READ
fprintf(stdout,"#DBG: USPEX optType=%i\n",_UC.type);
#endif

		}
		if (find_in_string("atomType",line) !=NULL){
		  /*let's count number of species here*/
		  g_free(line);
		  line = file_read_line(vf);
		  ptr=&(line[0]);
		  while(*ptr==' ') ptr++;/*skip initial space (if any)*/
		  _UC.nspecies=1;/*there is at lease one species*/
		  while(*ptr!='\0'){
		  	if(*ptr==' ') {
	 	 		_UC.nspecies++;ptr++;
				while(g_ascii_isgraph(*ptr)) ptr++;/*go to next space/end*/
 	 		}else ptr++;
  		  }
#if DEBUG_USPEX_READ
fprintf(stdout,"#DBG: USPEX number_of_species=%i\n",_UC.nspecies);
#endif
		  /*now let's find each species symbol*/
		  _UC.spe_Z=g_malloc(_UC.nspecies*sizeof(gint));
		  ptr=&(line[0]);idx=0;
		  while(*ptr==' ') ptr++;/*skip initial space (if any)*/
		  while((*ptr!='\n')&&(*ptr!='\0')){
		  	if(*ptr==' ') while(*ptr==' ') ptr++;/*skip space(s)*/
        		if(g_ascii_isdigit(*ptr)) _UC.spe_Z[idx]=g_ascii_strtod(ptr,NULL);/*user provided Z number*/
			else _UC.spe_Z[idx]=elem_symbol_test(ptr);/*user provided symbol*/
			if((_UC.spe_Z[idx]<=0)||(_UC.spe_Z[idx]>MAX_ELEMENTS-1)){/*X not allowed*/
				g_free(line);
				line = g_strdup_printf("ERROR: USPEX Parameters.txt file contain invalid atom definition!\n");
				gui_text_show(ERROR, line);
				g_free(line);
				return -1;
			}
#if DEBUG_USPEX_READ
fprintf(stdout,"#DBG: USPEX species[%i] Z=%i\n",idx,_UC.spe_Z[idx]);
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

gint read_individuals_uspex(gchar *filename, struct model_pak *model){
        FILE *vf;
	long int vfpos;
        gchar *line=NULL;
	gchar tmp[64];
        uspex_output_struct *uspex_calc=model->uspex;
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
	if (find_in_string("Fitness",line) != NULL) _UC.have_fitness=TRUE;
	else _UC.have_fitness=FALSE;
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
	if(_UC.method==US_CM_META){
		g_free(line);
		line = file_read_line(vf);
		vfpos=ftell(vf);/* flag */
	}
	/*do a 1st pass to count stuctures*/
	_UC.num_struct=0;
        while (line){
		_UC.num_struct++;
		sscanf(line," %*i %i %*s",&idx);
		if(idx>max_struct) max_struct=idx;
		g_free(line);
		line = file_read_line(vf);
	}
	_UC.num_struct--;
	fseek(vf,vfpos,SEEK_SET);/* rewind to flag */
	line = file_read_line(vf);
	/* now prepare arrays <- JOB=USPEX/300 example */
	if(_UC.num_struct==0) return 1;
	_UC.red_index=g_malloc(max_struct*sizeof(gint));
	for(idx=0;idx<max_struct;idx++) _UC.red_index[idx]=0;
	sscanf(line," %*[^[][%[^]]] %lf %*f %*s",
		&(tmp[0]),&(_UC.min_E));
	if(!e_red) _UC.min_E/=(gdouble)_UC.ind[0].natoms;
	_UC.max_E=_UC.min_E;
	idx=0;_UC.num_gen=1;
	while (line){
if(_UC.have_fitness){
	sscanf(line," %i %i%*[^[][%[^]]] %lf %lf %*f %lf %*s",
		&(_UC.ind[idx].gen),&(red_index),&(tmp[0]),&(_UC.ind[idx].energy),&(_UC.ind[idx].volume),&(_UC.ind[idx].fitness));
	if(isnan(_UC.ind[idx].fitness)) _UC.ind[idx].fitness=10000;
}else
	sscanf(line," %i %i%*[^[][%[^]]] %lf %lf %*s",
		&(_UC.ind[idx].gen),&(red_index),&(tmp[0]),&(_UC.ind[idx].energy),&(_UC.ind[idx].volume));
		if(_UC.ind[idx].gen<_UC.num_gen) _UC.ind[idx].gen++;
		if(_UC.ind[idx].gen>_UC.num_gen) _UC.num_gen=_UC.ind[idx].gen;
		_UC.red_index[red_index]=idx;
		if(!e_red) _UC.ind[idx].E=_UC.ind[idx].energy/_UC.ind[idx].natoms;
		else _UC.ind[idx].E=_UC.ind[idx].energy;
#if DEBUG_USPEX_READ
fprintf(stdout,"#DBG: USPEX[Individual]: GEN=%i STRUCT=%i e=%lf n_atom=%i E=%lf\n",_UC.ind[idx].gen,idx+1,_UC.ind[idx].energy,_UC.ind[idx].natoms,_UC.ind[idx].E);
#endif
		if(_UC.ind[idx].E<_UC.min_E) _UC.min_E=_UC.ind[idx].E;
		if(_UC.ind[idx].E>_UC.max_E) _UC.max_E=_UC.ind[idx].E;
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
	uspex_output_struct *uspex_calc;
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
	uspex_calc=model->uspex;
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
				_UC.version=min+10*med+100*max;
				if(_UC.version!=944){
					g_free(line);
					line = g_strdup_printf("WARNING: USPEX version %i unsupported!\n",_UC.version);
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
		_UC.dim=ix;
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
	if(ix==0) _UC.mol=FALSE;
	else if(ix==1) _UC.mol=TRUE;
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
	if(ix==0) _UC.var=FALSE;
	else if(ix==1) _UC.var=TRUE;
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
	if((job!=-1)&&(job!=_UC.type)){
/*since VCNEB will fail here due to a non-unified OUTPUT.txt format
 *we need to be a little more permissive... */
		line = g_strdup_printf("ERROR: inconsistent USPEX files (%i!=%i)!\n",job,_UC.type);
		gui_text_show(ERROR, line);
		g_free(line);
		goto uspex_fail;
	}
	g_free(aux_file);
	/*special case: no nspecies information*/
	if(_UC.nspecies==0) _UC.nspecies=num;
if((_UC.method==US_CM_USPEX)||(_UC.method==US_CM_META)){
/* --- if calculationMethod={USPEX,META} */
	model->basename=g_strdup_printf("uspex");
/* we have to open structure file first because of inconsistency in atom reporting in Individuals file*/
/* --- READ gatheredPOSCARS <- this is going to be the main model file */
/* ^^^ META: If gatheredPOSCARS_relaxed exists, read it instead */
        if(_UC.method==US_CM_META) {
                aux_file = g_strdup_printf("%s%s",res_folder,"gatheredPOSCARS_relaxed");
                vf = fopen(aux_file, "rt");
                if (!vf) aux_file = g_strdup_printf("%s%s",res_folder,"gatheredPOSCARS");
                else fclose(vf);
        }else{  
                aux_file = g_strdup_printf("%s%s",res_folder,"gatheredPOSCARS");
        }
        vf = fopen(aux_file, "rt");
        if (!vf) {
                if(_UC.ind!=NULL) g_free(_UC.ind);
                line = g_strdup_printf("ERROR: can't open USPEX gatheredPOSCARS file!\n");
                gui_text_show(ERROR, line);
                g_free(line);
                goto uspex_fail;
        }
	/*count number of structures*/
	_UC.num_struct=0;
	line = file_read_line(vf);
	while(!feof(vf)){
		if((line[0]=='E')&&(line[1]=='A')) _UC.num_struct++;
		g_free(line);
		line = file_read_line(vf);
	}
	rewind(vf);
        model->num_frames=0;
        vfpos=ftell(vf);/* flag */
        line = file_read_line(vf);
	_UC.ind=g_malloc((_UC.num_struct)*sizeof(uspex_individual));
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
			_UC.ind[ix].atoms=g_malloc(_UC.nspecies*sizeof(gint));
			for(jdx=0;jdx<_UC.nspecies;jdx++) _UC.ind[ix].atoms[jdx]=0;/*init atoms*/
			_UC.ind[ix].natoms=0;
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
			if(jdx<_UC.nspecies){
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
					while((jdx<_UC.nspecies)&&(_UC.spe_Z[jdx]!=elem_symbol_test(ptr))) jdx++;
					_UC.ind[ix].atoms[jdx]=g_ascii_strtod(ptr2,NULL);
					_UC.ind[ix].natoms+=_UC.ind[ix].atoms[jdx];
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
					_UC.ind[ix].atoms[jdx]=g_ascii_strtod(ptr,&ptr2);
					_UC.ind[ix].natoms+=_UC.ind[ix].atoms[jdx];
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
	if(_UC.method==US_CM_META) {
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
	if(_UC.method==US_CM_META) {
		aux_file = g_strdup_printf("%s%s",res_folder,"gatheredPOSCARS_relaxed");
		vf = fopen(aux_file, "rt");
		if (!vf) aux_file = g_strdup_printf("%s%s",res_folder,"gatheredPOSCARS");
		else fclose(vf);
	}else{
		aux_file = g_strdup_printf("%s%s",res_folder,"gatheredPOSCARS");
	}
        vf = fopen(aux_file, "rt");
        if (!vf) {/*very unlikely, we just opened it*/
                if(_UC.ind!=NULL) g_free(_UC.ind);
                line = g_strdup_printf("ERROR: can't open USPEX gatheredPOSCARS file!\n");
                gui_text_show(ERROR, line);
                g_free(line);
                goto uspex_fail;
        }
	read_frame_uspex(vf,model);/*open first frame*/
        fclose(vf);

/* --- Prepare the summary graph */
	/*no need to load BESTIndividuals since we already know everything from Individuals*/
	_UC.graph=graph_new("ALL", model);
	_UC.graph_best=graph_new("BEST", model);
	if(_UC.num_struct>_UC.num_gen){
		skip=0;
		gen=1;
		_UC.num_best=_UC.num_gen;
		_UC.best_ind=g_malloc((1+_UC.num_best)*sizeof(gint));
		_UC.best_ind[0]=0;
		ix=0;
		while(gen<=_UC.num_gen){
		  while(_UC.ind[ix].gen<gen) ix++;/*reach generation gen*/
		  num=0;
		  skip=ix;
		  min_E=_UC.ind[ix].E;
		  if(_UC.have_fitness) min_F=_UC.ind[ix].fitness;
		  _UC.best_ind[gen]=ix;
		  while(_UC.ind[ix].gen==gen){
			if(_UC.have_fitness){
				if(_UC.ind[ix].fitness<min_F){
					min_F=_UC.ind[ix].fitness;
					_UC.best_ind[gen]=ix;
				}
			}else{
				if(_UC.ind[ix].E<min_E){
					min_E=_UC.ind[ix].E;
					_UC.best_ind[gen]=ix;
				}
			}
			num++;
			ix++;
		  }
		  e=g_malloc((1+num)*sizeof(gdouble));
		  ix=skip;
		  e[0]=(gdouble)num;
		  while(_UC.ind[ix].gen==gen){
			e[ix-skip+1]=_UC.ind[ix].E;
			ix++;
		  }
		  graph_add_borned_data(
			_UC.num_gen,e,0,_UC.num_gen,
			_UC.min_E-(_UC.max_E-_UC.min_E)*0.05,_UC.max_E+(_UC.max_E-_UC.min_E)*0.05,GRAPH_USPEX,_UC.graph);
		  gen++;
		  g_free(e);
		}
		/*prepare array*/
		e=g_malloc((1+_UC.num_gen)*sizeof(gdouble));
		e[0]=(gdouble)_UC.num_gen;
		max_E=_UC.ind[_UC.best_ind[1]].E;
		min_E=max_E;
		for(gen=0;gen<_UC.num_gen;gen++) {
			if(_UC.ind[_UC.best_ind[gen+1]].E>max_E) max_E=_UC.ind[_UC.best_ind[gen+1]].E;
			if(_UC.ind[_UC.best_ind[gen+1]].E<min_E) min_E=_UC.ind[_UC.best_ind[gen+1]].E;
			e[gen+1]=_UC.ind[_UC.best_ind[gen+1]].E;
#if DEBUG_USPEX_READ
fprintf(stdout,"#DBG: USPEX_BEST: GEN=%i STRUCT=%i E=%lf\n",gen+1,1+_UC.best_ind[gen+1],e[gen+1]);
#endif
		}
		if(_UC.method==US_CM_META) {/*shift best_ind by 1*/
			for(gen=1;gen<=_UC.num_gen;gen++) _UC.best_ind[gen]++;
		}
		graph_add_borned_data(
			_UC.num_gen,e,0,_UC.num_gen,
			min_E-(max_E-min_E)*0.05,max_E+(max_E-min_E)*0.05,GRAPH_USPEX_BEST,_UC.graph_best);
		g_free(e);
	}
	/*set ticks*/
	if(_UC.num_gen>15) ix=5;
	else ix=_UC.num_gen+1;
	graph_set_xticks(TRUE,ix,_UC.graph);
	graph_set_yticks(TRUE,5,_UC.graph);
	graph_set_xticks(TRUE,ix,_UC.graph_best);
	graph_set_yticks(TRUE,5,_UC.graph_best);
/* --- variable composition */
if(_UC.var){
/* --- NEW: add a simple composition % graph for each species */
if(_UC.nspecies>1){
n_compo=2;
ix=0;
for(species_index=0;species_index<_UC.nspecies;species_index++){
	c=g_malloc(2*sizeof(gdouble));/*suppose at least 2 different compositions*/
	c[0]=(gdouble)_UC.ind[0].atoms[species_index] / (gdouble)_UC.ind[0].natoms;
	c[1]=c[0];
	compo=c[0];
	ix=0;
	while((compo==c[0])&&(ix<_UC.num_struct)){
		compo=(gdouble)_UC.ind[ix].atoms[species_index] / (gdouble)_UC.ind[ix].natoms;
		if(compo!=c[0]){
			if(compo>c[0]) c[1]=compo;
			else {
				c[1]=c[0];c[0]=compo;
			}
		}
		ix++;
	}
	if(c[1]==c[0]) continue;/*in VARCOMP calculation, one species can be kept fixed*/
	for(ix=0;ix<_UC.num_struct;ix++){
		compo=(gdouble)_UC.ind[ix].atoms[species_index] / (gdouble)_UC.ind[ix].natoms;
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
		max_E=_UC.min_E;
		for(ix=0;ix<n_compo;ix++) e[ix]=_UC.max_E;
		for(ix=0;ix<_UC.num_struct;ix++){
			for(jdx=0;jdx<n_compo;jdx++){
				compo=(gdouble)_UC.ind[ix].atoms[species_index] / (gdouble)_UC.ind[ix].natoms;
				if(c[jdx]==compo){
					if(_UC.ind[ix].E < e[jdx]) {
						e[jdx] = _UC.ind[ix].E;
						tag[jdx] = (gdouble) ix;
						if(e[jdx]>max_E) max_E=e[jdx];
					}
				}
			}
		}
		/* prepare composition diagram*/
		line=g_strdup_printf("COMP_%s",elements[_UC.spe_Z[species_index]].symbol);
		_UC.graph_comp=graph_new(line, model);
                graph_add_borned_data(
			n_compo,tag,0,1,_UC.min_E-(max_E-_UC.min_E)*0.05,max_E+(max_E-_UC.min_E)*0.05,GRAPH_USPEX_2D,_UC.graph_comp);
		graph_add_borned_data(
			n_compo,c,0,1,_UC.min_E-(max_E-_UC.min_E)*0.05,max_E+(max_E-_UC.min_E)*0.05,GRAPH_USPEX_2D,_UC.graph_comp);
		graph_add_borned_data(
			n_compo,e,0,1,_UC.min_E-(max_E-_UC.min_E)*0.05,max_E+(max_E-_UC.min_E)*0.05,GRAPH_USPEX_2D,_UC.graph_comp);
                g_free(line);
		g_free(tag);
		g_free(c);
		g_free(e);
        	/*set ticks*/
        	graph_set_xticks(TRUE,3,_UC.graph_comp);
        	graph_set_yticks(TRUE,5,_UC.graph_comp);

	}

}/*variable composition with only 1 species...*/
}/*not a VARCOMP calculation*/

	/*refresh model*/
	tree_model_refresh(model);
        model->redraw = TRUE;
        /* always show this information */
        gui_text_show(ITALIC,g_strdup_printf("USPEX: %i structures detected.\n",model->num_frames));
        model_prep(model);


}else if(_UC.method==US_CM_VCNEB){
/* --- tentative at adding VCNEB, which format is different**/
/* ^^^ *BUT unfortunately, we NEED a vasp5 output (VCNEB stayed at VASP4) */
/*      thus we create a new file called transitionPath_POSCARs5 which is a
 *      translation of transitionPath_POSCARs for VASP5...*/
	atoms = g_strdup_printf("%s",elements[_UC.spe_Z[0]].symbol);
	for(idx=1;idx<_UC.nspecies;idx++){
		line = g_strdup_printf("%s %s",atoms,elements[_UC.spe_Z[idx]].symbol);
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
	_UC.num_struct=0;
	line = file_read_line(f_src);
	while(line){
		fprintf(f_dest,"%s",line);
		if (find_in_string("Image",line) != NULL) {
			idx=0;
			_UC.num_struct++;
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
	_UC.ind=g_malloc((_UC.num_struct)*sizeof(uspex_individual));
/* --- read energies from the Energy file */
	g_free(aux_file);
	aux_file = g_strdup_printf("%s%s",res_folder,"Energy");
	vf=fopen(aux_file,"rt");
	if (!vf) {
		g_free(_UC.ind);
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
	_UC.best_ind=g_malloc((1+_UC.num_struct)*sizeof(gint));
	e = g_malloc((_UC.num_struct+1)*sizeof(gdouble));
	e[0]=(gdouble) _UC.num_struct;/*first element is the size*/
	sscanf(line," %*i %lf %*s",&(e[1]));
	e[1] /= (gdouble)natoms;
	_UC.min_E=e[1];
	_UC.max_E=e[1];
	for(idx=0;idx<_UC.num_struct;idx++){
		if(!line) break;/*<- useful? */
		sscanf(line," %*i %lf %*s",&(_UC.ind[idx].energy));
		_UC.ind[idx].atoms=g_malloc(_UC.nspecies*sizeof(gint));
		_UC.ind[idx].natoms=natoms;/*<-constant*/
		_UC.ind[idx].E = _UC.ind[idx].energy / (gdouble)_UC.ind[idx].natoms;
		_UC.ind[idx].gen=idx;/*<- this is actually *not* true*/
		_UC.best_ind[idx+1]=idx;
		e[idx+1]=_UC.ind[idx].E;
		if(e[idx+1]<_UC.min_E) _UC.min_E=e[idx+1];
		if(e[idx+1]>_UC.max_E) _UC.max_E=e[idx+1];
#if DEBUG_USPEX_READ
fprintf(stdout,"#DBG: USPEX[VCNEB] Image %i: natoms=%i energy=%lf e=%lf\n",_UC.ind[idx].gen,_UC.ind[idx].natoms,_UC.ind[idx].energy,e[idx]);
#endif
		g_free(line);
		line = file_read_line(vf);
	}
	fclose(vf);
	g_free(aux_file);
	/*create the reduced index table*/
	_UC.red_index=g_malloc(_UC.num_struct*sizeof(gint));
	for(idx=0;idx<_UC.num_struct;idx++) _UC.red_index[idx]=idx;
	aux_file = g_strdup_printf("%s%s",res_folder,"transitionPath_POSCARs5");
/* --- reopen the newly created transitionPath_POSCARs5 file for reading <- this will be our new model file*/
	vf=fopen(aux_file,"rt");
        if (!vf) {/*very unlikely: we just create it*/
		g_free(_UC.ind);
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
        _UC.graph_best=graph_new("PATH", model);
	_UC.num_gen=_UC.num_struct;
        graph_add_borned_data(
                _UC.num_gen,e,0,_UC.num_gen,
                _UC.min_E-(_UC.max_E-_UC.min_E)*0.05,_UC.max_E+(_UC.max_E-_UC.min_E)*0.05,GRAPH_USPEX_BEST,_UC.graph_best);
        g_free(e);
        /*set ticks*/
        if(_UC.num_gen>15) ix=5;
        else ix=_UC.num_gen+1;
        graph_set_xticks(TRUE,ix,_UC.graph_best);
        graph_set_yticks(TRUE,5,_UC.graph_best);


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
	uspex_output_struct *uspex_calc=model->uspex;
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
/* now the index should always be _UC.red_index[idx-1] <- I THINK*/
	line=g_strdup_printf("%lf eV",_UC.ind[_UC.red_index[idx]].energy);
	property_add_ranked(3, "Energy", line, model);
	model->volume=_UC.ind[_UC.red_index[idx]].volume;
	g_free(line);
	line=g_strdup_printf("%i",idx);
	property_add_ranked(8, "Structure", line, model);
	g_free(line);
	return 0;
}

#undef _UC

