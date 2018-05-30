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
 *	1/ load the optimal structure as base model.
 *		-> each optimal structure, for each generation, is represented as a frame
 *		   in this model.
 *	2/ make a graph representing all structures energy levels vs.  generation numbers
 *		-> make each point on that graph clickable so that selecting a point will
 *		   open a new model containing the corresponding structure.
 * TODO for later:
 * 	1/ include other data plots:
 * 		- hardness vs. generation
 * 		- energy vs. stoechiometry
 * 		- energy vs. image number
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
        uspex_calc_struct *uspex_calc=model->uspex;
	/* specific */
	gchar method;
	/*start*/
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
        uspex_calc_struct *uspex_calc=model->uspex;
        /* specific */
	gint idx=0;
	gchar *ptr;
	gchar *ptr2;
	gint spe;
	/*start*/
        vf = fopen(filename, "rt");
        if (!vf) return 1;
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
	/**/
	/*do a 1st pass to count stuctures*/
	_UC.num_struct=0;
        while (line){
		_UC.num_struct++;
		g_free(line);
		line = file_read_line(vf);
	}
	_UC.num_struct--;
	fseek(vf,vfpos,SEEK_SET);/* rewind to flag */
	line = file_read_line(vf);
	/* now prepare arrays <- JOB=USPEX/300 example */
	if(_UC.num_struct==0) return 1;
	_UC.ind=g_malloc(_UC.num_struct*sizeof(uspex_individual));
	sscanf(line," %*i %*i %*[a-zA-Z] [%[^]]] %lf %*f %*s",
		&(tmp[0]),&(_UC.min_E));
	idx=0;ptr=&(tmp[0]);ptr2=ptr;
	_UC.nspecies=0;
	do{/*get number of atoms*/
		idx+=g_ascii_strtod(ptr,&ptr2);
		_UC.nspecies++;
		if(ptr2==ptr) break;
		ptr=ptr2;
	}while(1);
	_UC.nspecies--;
	_UC.min_E/=(gdouble)idx;
	_UC.max_E=_UC.min_E;
	idx=0;_UC.num_gen=1;
	while (line){
		sscanf(line," %i %*i %*[a-zA-Z] [%[^]]] %lf %lf %*s",
			&(_UC.ind[idx].gen),&(tmp[0]),&(_UC.ind[idx].energy),&(_UC.ind[idx].volume));
		if(_UC.ind[idx].gen>_UC.num_gen) _UC.num_gen=_UC.ind[idx].gen;
		/*calculate number of atoms*/
		_UC.ind[idx].natoms=0;ptr=&(tmp[0]);ptr2=ptr;
		spe=0;_UC.ind[idx].atoms=g_malloc(_UC.nspecies*sizeof(gint));
		do{
			_UC.ind[idx].atoms[spe]=g_ascii_strtod(ptr,&ptr2);
			_UC.ind[idx].natoms+=_UC.ind[idx].atoms[spe];
			if(ptr2==ptr) break;
			spe++;
			ptr=ptr2;
		}while(1);
		_UC.ind[idx].E=_UC.ind[idx].energy/_UC.ind[idx].natoms;
		if(_UC.ind[idx].E<_UC.min_E) _UC.min_E=_UC.ind[idx].E;
		if(_UC.ind[idx].E>_UC.max_E) _UC.max_E=_UC.ind[idx].E;
		g_free(line);
		line = file_read_line(vf);
		idx++;
	}

	/*all done*/
	return 0;
}

gint read_best_individuals_uspex(gchar *filename, struct model_pak *model){
        FILE *vf;
        long int vfpos;
        gchar *line=NULL;
	gchar *ptr;
	gint idx;
        uspex_calc_struct *uspex_calc=model->uspex;
        /*start*/
        vf = fopen(filename, "rt");
        if (!vf) return 1;
	/* Skip header */
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
	/*get the number of best structure (should be the same as generation number, but...)*/
	_UC.num_best=0;
	line = file_read_line(vf);
	while(line){
		_UC.num_best++;
		g_free(line);
		line = file_read_line(vf);
	}
	fseek(vf,vfpos,SEEK_SET);/* rewind to flag */
	/*now alloc*/
	_UC.best_ind=g_malloc(_UC.num_best*sizeof(gint));
	line = file_read_line(vf);
	idx=0;
	while(line){
		sscanf(line," %*i %i %*s",&(_UC.best_ind[idx]));
		g_free(line);
		line = file_read_line(vf);
		idx++;
	}
	/*that's all*/
	return 0;
}

gint read_output_uspex(gchar *filename, struct model_pak *model){
	gchar *line=NULL;
	FILE *vf;
	long int vfpos;
	gchar *ptr;
        gchar *res_folder;
	gchar *aux_file;
	/* of interest */
	gint ix=0;
	gint job;
	uspex_calc_struct *uspex_calc;
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
	model->uspex=g_malloc(sizeof(uspex_calc_struct));
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
	
/* --- read the system parameters (AS A CHECK)*/
	if(fetch_in_file(vf,"Block for system description")==0) goto uspex_fail;
	/*next line contains Dimensionality*/
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
	if(job!=_UC.type) {
		line = g_strdup_printf("ERROR: inconsistent USPEX files (%i!=%i)!\n",job,_UC.type);
		gui_text_show(ERROR, line);
		g_free(line);
		goto uspex_fail;
	}
/* --- if calculationMethod=USPEX */
	g_free(aux_file);
if(_UC.method==US_CM_USPEX){
	model->basename=g_strdup_printf("uspex");
/* --- READ Individuals <- information about all structures */
	aux_file = g_strdup_printf("%s%s",res_folder,"Individuals");
	if(read_individuals_uspex(aux_file,model)!=0) {
		line = g_strdup_printf("ERROR: reading USPEX Individuals file!\n");
		gui_text_show(ERROR, line);
		g_free(line);
		goto uspex_fail;
	}
/* --- READ gatheredPOSCARS <- this is going to be the main model file */
	aux_file = g_strdup_printf("%s%s",res_folder,"gatheredPOSCARS");
        vf = fopen(aux_file, "rt");
        if (!vf) {
                if(_UC.ind!=NULL) g_free(_UC.ind);
                line = g_strdup_printf("ERROR: can't open USPEX gatheredPOSCARS file!\n");
                gui_text_show(ERROR, line);
                g_free(line);
                goto uspex_fail;
        }
        model->num_frames=0;
	vfpos=ftell(vf);/* flag */
	line = file_read_line(vf);
        while(!feof(vf)){
		if((line[0]=='E')&&(line[1]=='A')) {
			fseek(vf,vfpos,SEEK_SET);/* rewind to flag */
			add_frame_offset(vf, model);
			model->num_frames++;
			g_free(line);
			line = file_read_line(vf);
		}
		vfpos=ftell(vf);/* flag */
		g_free(line);
		line = file_read_line(vf);
        }
        strcpy(model->filename,aux_file);// which means that we "forget" about OUTPUT.txt
	rewind(vf);
	read_frame_uspex(vf,model);/*open first frame*/
        fclose(vf);

/* --- Prepare the summary graph */
	/*no need to load BESTIndividuals since we already know everything from Individuals*/
	_UC.graph=graph_new("ALL", model);
	_UC.graph_best=graph_new("BEST", model);
	if(_UC.num_struct>_UC.num_gen){
		gdouble *e;
		gint skip=0;
		gint gen=1;
		gint num;
		gdouble min_E;
		gdouble max_E;
		_UC.num_best=_UC.num_gen;
		_UC.best_ind=g_malloc((1+_UC.num_best)*sizeof(gint));
		_UC.best_ind[0]=0;
		ix=0;
		while(gen<=_UC.num_gen){
		  while(_UC.ind[ix].gen<gen) ix++;/*reach generation gen*/
		  num=0;
		  skip=ix;
		  min_E=_UC.ind[ix].E;
		  _UC.best_ind[gen]=ix;
		  while(_UC.ind[ix].gen==gen){
			if(_UC.ind[ix].E<min_E){
				min_E=_UC.ind[ix].E;
				_UC.best_ind[gen]=ix;
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
		for(gen=0;gen<_UC.num_gen;gen++) {
			if(_UC.ind[_UC.best_ind[gen+1]].E>max_E) max_E=_UC.ind[_UC.best_ind[gen+1]].E;
			e[gen+1]=_UC.ind[_UC.best_ind[gen+1]].E;
#if DEBUG_USPEX_READ
fprintf(stdout,"#DBG: USPEX_BEST: GEN=%i STRUCT=%i E=%lf\n",gen+1,1+_UC.best_ind[gen+1],e[gen+1]);
#endif
		}
		graph_add_borned_data(
			_UC.num_gen,e,0,_UC.num_gen,
			_UC.min_E-(max_E-_UC.min_E)*0.05,max_E+(max_E-_UC.min_E)*0.05,GRAPH_USPEX_BEST,_UC.graph_best);
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
	/*add a graph with binary composition, if binary*/
	if(_UC.nspecies==2){
		gdouble *e;
		gdouble *c;
		gdouble *tag;
		gdouble *tmp;
		gdouble compo;
		gdouble max_E;
		gint n_compo;/*number of different compositions*/
		gint jdx;
		n_compo=2;
		ix=0;
		/*check for 2 different compositions <- useful?*/
		c=g_malloc(2*sizeof(gdouble));
		c[0]=(gdouble)_UC.ind[0].atoms[1] / (gdouble)_UC.ind[0].natoms;
		c[1]=c[0];
		do{
			compo=(gdouble)_UC.ind[ix].atoms[1] / (gdouble)_UC.ind[ix].natoms;
			if(compo!=c[0]){
				if(compo>c[0]) c[1]=compo;
				else {
					c[1]=c[0];c[0]=compo;
				}
			}
			ix++;
		}while((compo==c[0])&&(ix<_UC.num_struct));
if(c[1]!=c[0]){
		/*we have 2 different composition*/
		/*build the compo array*/
		for(ix=0;ix<_UC.num_struct;ix++){
			compo=(gdouble)_UC.ind[ix].atoms[1] / (gdouble)_UC.ind[ix].natoms;
			if(compo<c[0]) {/*we have a new first element*/
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
				compo=(gdouble)_UC.ind[ix].atoms[1] / (gdouble)_UC.ind[ix].natoms;
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
		_UC.graph_comp=graph_new("COMP", model);
                graph_add_borned_data(
			n_compo,tag,0,1,_UC.min_E-(max_E-_UC.min_E)*0.05,max_E+(max_E-_UC.min_E)*0.05,GRAPH_USPEX_2D,_UC.graph_comp);
		graph_add_borned_data(
			n_compo,c,0,1,_UC.min_E-(max_E-_UC.min_E)*0.05,max_E+(max_E-_UC.min_E)*0.05,GRAPH_USPEX_2D,_UC.graph_comp);
		graph_add_borned_data(
			n_compo,e,0,1,_UC.min_E-(max_E-_UC.min_E)*0.05,max_E+(max_E-_UC.min_E)*0.05,GRAPH_USPEX_2D,_UC.graph_comp);
                g_free(tag);
		g_free(c);
		g_free(e);
        	/*set ticks*/
        	graph_set_xticks(TRUE,3,_UC.graph_comp);
        	graph_set_yticks(TRUE,5,_UC.graph_comp);
}
/*2 species, but only one composition?*/
	}
}

	/*refresh model*/
	tree_model_refresh(model);
        model->redraw = TRUE;
        /* always show this information */
        gui_text_show(ITALIC,g_strdup_printf("USPEX: %i structures detected.\n",model->num_frames));
        model_prep(model);


}else{/*calculation is not USPEX*/
	line = g_strdup_printf("ERROR: USPEX support is limited to calculationMethod=USPEX for now!\n");
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
	uspex_calc_struct *uspex_calc=model->uspex;
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
	line=g_strdup_printf("%lf eV",_UC.ind[idx-1].energy);
	property_add_ranked(3, "Energy", line, model);
	model->volume=_UC.ind[idx-1].volume;
	g_free(line);
	line=g_strdup_printf("%i",idx);
	property_add_ranked(8, "Structure", line, model);
	g_free(line);
	return 0;
}

#undef _UC

