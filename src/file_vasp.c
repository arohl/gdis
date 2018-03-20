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

/* simple VASP xml Parser, _not_ using xml parser (libxml2 was buggy). */

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

enum {VASP_DEFAULT, VASP_XML, VASP_TRIKS, VASP_NLOOPS};

/* main structures */
extern struct sysenv_pak sysenv;
extern struct elem_pak elements[];
/* TODO: use proper read_token */
#define find_in_string(a,b) strstr(b,a)

long int fetch_in_file(FILE *vf,const gchar *target){
	gchar *line;
	line = file_read_line(vf);
	while (line){
		if (strstr(line,target) != NULL) break;
		g_free(line);
		line = file_read_line(vf);
	}
	if(feof(vf)) return 0;
	return ftell(vf);
}

int vasp_xml_read_header(FILE *vf){
/* read the vasp xml header */
	gchar *line;
	int isok=0;
	line = file_read_line(vf);
	while (line){
		if (find_in_string("/generator",line) != NULL) break;
		if (find_in_string("program",line) != NULL) isok+=10*(find_in_string("vasp", line)!=NULL);
		/* TODO: read version*/
		g_free(line);
		line = file_read_line(vf);
	}
	isok*=(int)(line!=0);
	g_free(line);
	return (isok);
/* return value:
 * 0 => wrong file / file type
 * 1 => not vasp ?
 * 10 => ok
*/
}

int vasp_xml_read_org_incar(FILE *vf){
/* TODO: include original INCAR */
	gchar *line;
	int isok;
	line = file_read_line(vf);
	while (line) {/*just skip for now*/
		if (find_in_string("/incar",line) != NULL) break;
		g_free(line);
		line = file_read_line(vf);
	}
	isok=(int)(line != 0);
	g_free(line);
	return (isok);
}

int vasp_xml_read_kpoints(FILE *vf,struct model_pak *model){
/* get number of kpoints and their distance to gamma */
	gchar *line;
	gint nkpts;
	gdouble xo,yo,zo;
	gdouble xi,yi,zi;
	gdouble dist=0.;
	long int vfpos=ftell(vf);/*flag begining of kpoints*/
	line = file_read_line(vf);
/*unfortunately, there is no NKPTS entry in vasprun.xml!*/
	if(fetch_in_file(vf,"kpointlist")==0) {
		fseek(vf,vfpos,SEEK_SET);/* rewind to flag */
		return -1;/*broken kpoints entry*/
	}
	vfpos=ftell(vf);/*flag begining of kpointlist*/
	nkpts=0;
	line = file_read_line(vf);
	while (line) {/*first get the number of kpoints*/
	        if (find_in_string("/varray",line) != NULL) break;
		nkpts++;
		g_free(line);
		line = file_read_line(vf);
	}
	if (line == NULL) {
		fseek(vf,vfpos,SEEK_SET);/* rewind to flag */
		return -1;/* Incomplete kpoint list? */
	}
	model->nkpoints=nkpts;
	if(model->kpts_d!=NULL) g_free(model->kpts_d);
	model->kpts_d=g_malloc((model->nkpoints)*sizeof(gdouble));
	fseek(vf,vfpos,SEEK_SET);/* rewind to flag */
	/*now get kpoints distances*/
	nkpts=0;
	xo=0.;
	yo=0.;
	zo=0.;
	line = file_read_line(vf);
	while (line) {/*now get kpoints distances*/
		if (find_in_string("/varray",line) != NULL) break;
		sscanf(line," <v> %lf %lf %lf </v> ",&xi,&yi,&zi);
		dist+=sqrt((xi-xo)*(xi-xo)+(yi-yo)*(yi-yo)+(zi-zo)*(zi-zo));
		model->kpts_d[nkpts]=dist;
//fprintf(stdout,"#DBG: kpt[%i]={%lf,%lf,%lf} d=%lf\n",nkpts,xi,yi,zi,model->kpts_d[nkpts]);
		xo=xi;yo=yi;zo=zi;
		nkpts++;
		g_free(line);
		line = file_read_line(vf);
	}
	vfpos=ftell(vf);/*flag end of kpointlist*/
	/*all done for now! TODO: get kpoints coordinates and weights*/
	if(fetch_in_file(vf,"</kpoints>")==0) {
		fseek(vf,vfpos,SEEK_SET);/* rewind to flag */
		return -1;/*broken kpoints entry*/
	}
	return 0;
}

int vasp_xml_read_incar(FILE *vf,struct model_pak *model){
/* This is the 'full length' INCAR (mostly ignored). */
	gchar *line;
	int isok=0;
	int ii=0;
	int ispin=0;
	line = file_read_line(vf);
	while (line) {
	        if (find_in_string("/parameters",line) != NULL) break;
		if (find_in_string("SYSTEM",line) != NULL) {
			g_free(model->basename);
			model->basename=g_malloc((strlen(line)-35)*sizeof(gchar));
			sscanf(line," <i type=\"string\" name=\"SYSTEM\"> %s",model->basename);
			/*get rid of the </i> part*/
			while((model->basename[ii] != '<')&&(model->basename[ii] != '\0')) ii++;
			model->basename[ii]='\0';
		}
		/*fill needed properties*/
		if (find_in_string("NBANDS",line) != NULL)
			sscanf(line," <i type=\"int\" name=\"NBANDS\"> %i</i> ",&(model->nbands));
		if (find_in_string("ISPIN",line) != NULL){
			sscanf(line," <i type=\"int\" name=\"ISPIN\"> %i</i> ",&ispin);
			if(ispin>1) model->spin_polarized=TRUE;
		}
		if (find_in_string("NEDOS",line) != NULL)
			sscanf(line," <i type=\"int\" name=\"NEDOS\"> %i</i> ",&(model->ndos));
		g_free(line);
		line = file_read_line(vf);
	}
	isok=(int)(line != 0);
	g_free(line);
	return (isok);
}

int vasp_xml_read_atominfo(FILE *vf, struct model_pak *model){
/* read current information about atoms */
	gchar *line;
	int idx=0;
	int natom=0;
	int ii;
	gchar *label=NULL;
	gint number;
	gdouble charge;
	struct core_pak *core;
	/* Goto end of first <set> array */
	if(fetch_in_file(vf,"</set>")==0) return -1;
	/* Goto the begining of 2nd <set> array */
	if(fetch_in_file(vf,"<set>")==0) return -1;
	/* Construct the species list */
	line = file_read_line(vf);
	label = g_malloc(3*sizeof(gchar));
	while (line) {
		if (find_in_string("</set>",line) != NULL) break;
		sscanf(line," <rc><c>  %i</c><c>%2c</c><c> %*f</c><c> %lf</c><c> %*s",&number,label,&charge);
		label[2]='\0';
		for(ii=0;ii<number;ii++){
			core=new_core(label,model);
			core->charge=charge;
			model->cores=g_slist_append(model->cores,core);
		}
		idx++;
		natom+=number;
		g_free(line);
		line = file_read_line(vf);
	}
	free(label);
	if (line == NULL) return -1;/* Incomplete atom definition? */
	g_free(line);
	/* set total number of atoms */
	model->vasp.num_atoms=natom;
        model->vasp.num_species=idx;
	model->expected_cores=natom;
        model->expected_shells=0;
	model->num_atoms=natom;
        /* exit at </atominfo> */
	if(fetch_in_file(vf,"</atominfo>")==0) return -1;
	return 0;
}

int vasp_xml_read_energy(FILE *vf, struct model_pak *model){
	gchar *line;
	long int vfpos=ftell(vf);
	if(fetch_in_file(vf,"<energy>")==0) {
		/*we didnt't get it until EOF, which is normal when using finalpos*/
		rewind(vf);
		/*find the last valid <structure>*/
		while(fetch_in_file(vf,"<structure>")!=0) vfpos=ftell(vf);/*flag*/
		fseek(vf,vfpos,SEEK_SET);/* rewind to flag */
		if(fetch_in_file(vf,"<energy>")==0) return -1;/*still no <energy> tag?*/
	}
	if(fetch_in_file(vf,"e_fr_energy")==0) return -1;/*no energy information in <energy>?*/
	line = file_read_line(vf);/*next line is e_wo_entrp*/
	if (find_in_string("e_wo_entrp",line) != NULL) {
		sscanf(line," <i name=\"e_wo_entrp\"> %lf </i>",&model->vasp.energy);
		sprintf(line,"%lf eV",model->vasp.energy);
		property_add_ranked(3, "Energy", line, model);
	} else {
		return -1;
	}
	g_free(line);
	return 0;
}

int vasp_xml_read_frequency(FILE *vf, struct model_pak *model){
/*factor is actually:
 * factor = sqrt((1.60217733e-19/1e-20)/1.6605402e-27)/(2.*PI*2.99792458e+10)
 * calculated with quadmath precision. It convert eigenvalue of hessian to
 * frequency, in cm^{-1} 
*/
	gdouble factor=521.4708336735473879;
	gchar *line;
	long int vfpos;
	gint idx,jdx;

	if(fetch_in_file(vf,"<dynmat>")==0) return -1;/*we don't have a dynmat array*/
	vfpos=ftell(vf);
	if(fetch_in_file(vf,"eigenvalues")==0) return -1;/* no eigenvalues -> no frequency */
	fseek(vf,vfpos,SEEK_SET);/* rewind to flag */
	fetch_in_file(vf,"</varray>");
	line = file_read_line(vf);/*next line is eigenvalues*/
	/* count eigenvalues*/
	idx=0;
	while(line[idx]!='>') idx++;
	/*from here the number of eigenvalues is equal to the number of space characters :)*/
	model->nfreq=0;
	while(line[idx]!='\0') {
		if(line[idx]==' ') model->nfreq++;
		idx++;
	}
	if(model->nfreq>0) model->have_frequency=TRUE;
	else return -1;/*an eigenvalue array but no data?*/
	model->freq=g_malloc(model->nfreq*(sizeof(gdouble)));
	/*now, scan for each value TODO: do both at once?*/
	idx=0;jdx=0;
	while(line[idx]!='>') idx++;
	idx++;
	while(idx++){
		sscanf(&line[idx],"%lE %*s",&(model->freq[jdx]));
		model->freq[jdx]=factor*sqrt(fabs(model->freq[jdx]));
		while((line[idx]!='\0')&&(line[idx]!=' ')) idx++;/*go to next number, if any*/
		if(line[idx]=='\0') break;
		jdx++;
	}
	g_free(line);
	/*note that there is no frequency intensity information in vaspxml 
 * 	  and AFAIK there is no easy way to obtain such information. OVHPA*/
	return 0;
}


int vasp_xml_read_pos(FILE *vf,struct model_pak *model){
/* read the atoms positions */
	gchar *line;
	int idx=0;
	struct core_pak *core;
	/* Goto crystal varray "basis" */
	if(fetch_in_file(vf,"basis")==0) return -1;
	/* Get 3x3 vector lattice */
	line = file_read_line(vf);
	while (line) {
		if (find_in_string("/varray",line) != NULL) break;
		sscanf(line," <v> %lf %lf %lf </v>",&model->latmat[idx],&model->latmat[idx+3],&model->latmat[idx+6]);
		idx+=1;
		g_free(line);
		line = file_read_line(vf);
	}
	if (line == NULL) return -1;/*incomplete basis definition*/
	g_free(line);
	/*Always true in VASP*/
	model->fractional=TRUE;
	model->coord_units=ANGSTROM;
	model->construct_pbc = TRUE;
	model->periodic = 3;
	/* next line should be volume */
	line = file_read_line(vf);
	if (find_in_string("volume",line) != NULL) {
		sscanf(line," <i name=\"volume\"> %lf </i>",&model->volume);
	}
	g_free(line);
	/* TODO: gdouble rlatmat[9] can be filled */
	/* TODO: find the difference between ilatmat and rlatmat */
	/* Goto to "positions" */
	if(fetch_in_file(vf,"positions")==0) return -1;
	/* fill every atom position */
	line = file_read_line(vf);
	idx=0;core=g_slist_nth_data(model->cores,0);
	while (line) {
		if (find_in_string("/varray",line) != NULL) break;
		sscanf(line," <v> %lf %lf %lf </v>",&core->x[0],&core->x[1],&core->x[2]);
		idx++;
		core=g_slist_nth_data(model->cores,idx);
		g_free(line);
		line = file_read_line(vf);
	}
	if (line == NULL) return -1;/*incomplete positions*/
	g_free(line);
	if (idx != model->vasp.num_atoms) /* This should not happen */
		fprintf(stderr,"WARNING: Expecting %i atoms but got %i!\n",model->vasp.num_atoms,idx);
	/* look for energies */
	vasp_xml_read_energy(vf,model);
	return 0;
}
int vasp_xml_read_bands(FILE *vf,struct model_pak *model){
/* get the (non-projected) bandstructure */
        gchar *line;
        gint idx;
	gint ik,ib;
        long int vfpos=ftell(vf);
        /*start*/
#define DEBUG_BAND 0
	if(model->nkpoints<2) return 1;/*there is only one k-point*/
	if(model->nbands<2) return 1;/*there is only one band (how is that even possible?)*/
	if(fetch_in_file(vf,"<eigenvalues>")==0) {
		/*when there is no band information after <dos>
 * 		it is still possible to find some before.. */
		rewind(vf);
		if(fetch_in_file(vf,"<eigenvalues>")==0) {
			fseek(vf,vfpos,SEEK_SET);/* rewind to flag */
			return -1;/* no band information */
		}
	}
	/*NOTE: the projected information should be an <array> *AFTER* </eigenvalues>*/
	/*there is nbands * nkpoints * nspin values of energy to get*/
/*spin up*/
	if(fetch_in_file(vf,"kpoint 1")==0) {
		fseek(vf,vfpos,SEEK_SET);/* rewind to flag */
		return -1;/*no values*/
	}
        if(model->band_up != NULL) g_free(model->band_up);
        /*allocation TODO: g_malloc_try?*/
        model->band_up=g_malloc((model->nkpoints*model->nbands)*sizeof(gdouble));
	ik=0;ib=0;
	idx=0;
	line = file_read_line(vf);
	while(line){
		if(ib>(model->nbands-1)) {
			ik++;
			ib=0;
		}
		if(ik>(model->nkpoints-1)) break;
		if(find_in_string("</set>",line) != NULL) {
			g_free(line);
			line = file_read_line(vf);
			g_free(line);
			line = file_read_line(vf);
		}
		if(find_in_string("set comment",line) != NULL) break;/*next component*/
		if(find_in_string("</array>",line) != NULL) break;/*end of array*/
		idx=ik*(model->nbands)+ib;/*current index*/
		sscanf(line," <r> %lf %*s </r> ",&(model->band_up[idx]));
#if DEBUG_BAND
fprintf(stdout,"#DBG: spin=  up kpt=%i band=%i eval=%lf\n",ik,ib,model->band_up[idx]);
#endif
		ib++;
		g_free(line);
		line = file_read_line(vf);
	}
	if(line == NULL) return -1;/*incomplete set*/
	if(!model->spin_polarized) return 0;/*only one component*/
/*spin down*/
        if(fetch_in_file(vf,"kpoint 1")==0) {
                fseek(vf,vfpos,SEEK_SET);/* rewind to flag */
                return -1;/*no values*/
        }
        if(model->band_down != NULL) g_free(model->band_down);
        model->band_down=g_malloc((model->nkpoints*model->nbands)*sizeof(gdouble));
        ik=0;ib=0;
        idx=0;
        line = file_read_line(vf);
        while(line){
                if(ib>(model->nbands-1)) {
                        ik++;
                        ib=0;
                }
                if(ik>(model->nkpoints-1)) break;
                if(find_in_string("</set>",line) != NULL) {
                        g_free(line);
                        line = file_read_line(vf);
			g_free(line);
			line = file_read_line(vf);
                }
                if(find_in_string("set comment",line) != NULL) break;/*next (very unexpected) component*/
                if(find_in_string("</array>",line) != NULL) break;/*end of array*/
                idx=ik*(model->nbands)+ib;/*current index*/
                sscanf(line," <r> %lf %*s </r> ",&(model->band_down[idx]));
#if DEBUG_BAND
fprintf(stdout,"#DBG: spin=down kpt=%i band=%i eval=%lf\n",ik,ib,model->band_up[idx]);
#endif
                ib++;
                g_free(line);
                line = file_read_line(vf);
        }
        if(line == NULL) return -1;/*incomplete set*/
	return 0;
}



int vasp_xml_read_dos(FILE *vf, struct model_pak *model){
	/* get the final density of states (total dos only) */
	gchar *line;
	gint idx;
	/*start*/
	if(fetch_in_file(vf,"<dos>")==0) return -1;/* no dos information */
	line = file_read_line(vf);
	sscanf(line," <i name=\"efermi\"> %lf </i> ",&(model->efermi));
#define DEBUG_DOS 0
#if DEBUG_DOS
fprintf(stdout,"#DBG: FERMI ENERGY: %lf\n",model->efermi);
#endif
	/* the next line should contain <total> */
	line = file_read_line(vf);
	if(find_in_string("</set>",line) == NULL) {
		/*there is a pb: we have a dos, but not a total one*/
		/*solution: rewind and look for <total> keyword*/
		rewind(vf);
		if(fetch_in_file(vf,"<total>")==0) return -1;/*no total dos info?*/
	}
	/*go to first spin component*/
	if(fetch_in_file(vf,"spin 1")==0) return -1;/*no component*/
	if(model->dos_eval!=NULL) g_free(model->dos_eval);
	if(model->dos_spin_up!=NULL) g_free(model->dos_spin_up);
	model->dos_eval=g_malloc(model->ndos*sizeof(gdouble));
	model->dos_spin_up=g_malloc(model->ndos*sizeof(gdouble));
	line = file_read_line(vf);
	idx=0;
	while(line){
		if(idx>(model->ndos-1)) break;
		if(find_in_string("</set>",line) != NULL) break;
		sscanf(line," <r> %lf %lf %*s",&(model->dos_eval[idx]),&(model->dos_spin_up[idx]));
#if DEBUG_DOS
fprintf(stdout,"#DBG: CATCH SPIN UP: %i %lf \t %lf\n",idx,model->dos_eval[idx],model->dos_spin_up[idx]);
#endif
		idx++;
		g_free(line);
		line = file_read_line(vf);
	}
	if(!model->spin_polarized) return 0;/*only one component*/
	if(model->dos_spin_down!=NULL) g_free(model->dos_spin_down);
	model->dos_spin_down=g_malloc(model->ndos*sizeof(gdouble));
	line = file_read_line(vf);
	idx=0;
	while(line){
		if(idx>(model->ndos-1)) break;
		if(find_in_string("</set>",line) != NULL) break;
		sscanf(line," <r> %*f %lf %*s",&(model->dos_spin_down[idx]));
#if DEBUG_DOS
fprintf(stdout,"#DBG: CATCH SPIN DOWN: %i %lf \t %lf\n",idx,model->dos_eval[idx],model->dos_spin_down[idx]);
#endif
		idx++;
		g_free(line);
		line = file_read_line(vf);
	}
	return 0;
}
/* here will be some more features */
/* general parser / writter */
gint read_xml_vasp(gchar *filename, struct model_pak *model){
/* READER init for VASP XML
 * check the current status:
 * vaspxml =>
 * 	singlepoint calculation -> display finalpos (no frame)
 * 	partial calculation w/ nframe=1 -> display the 1st frame (no frame)
 * 	partial calculation w/ nframe>1 -> display the last frame (w/ frame)
 * 	normal calculation -> display finalpos (w/ frame)
*/
	int isok=0;
	int num_frames=1;
	FILE *vf;
	long int vfpos;
	gchar *line;
	g_return_val_if_fail(model != NULL, 1);
	g_return_val_if_fail(filename != NULL, 2);
	vf = fopen(filename, "rt");
	if (!vf) return 1;
	error_table_clear();
	/* some defaults */
	sysenv.render.show_energy = TRUE;
	/* start reading */
	line = file_read_line(vf);
	/* the first xml tag ie <?xml version="x.x" encoding="ISO-xxxx-x"?> */
	if (find_in_string("xml",line) == NULL) return 3;/* not even an xml file */
	g_free(line);
	/* <generator> tag */
	if(fetch_in_file(vf,"<generator>")==0) return 3;
	isok=vasp_xml_read_header(vf);
	if (isok == 10) gui_text_show(STANDARD,g_strdup_printf("VASP detected\n"));
	else {
		if (isok == 1) gui_text_show(WARNING,g_strdup_printf("not generated by vasp, will try to read however.\n"));
		if (isok == 0) return 3;/* not a valid vaspxml */
	}
vfpos=ftell(vf);/*flag*/
	/* <incar> tag - if none, rewind and ignore */
	if(fetch_in_file(vf,"<incar>")==0) fseek(vf,vfpos,SEEK_SET);
	else vasp_xml_read_org_incar(vf);
/* Starting vasp 5.4.X with X>1 primitive_cell and primitive_index are provided here. */
/* For now this information is discarded. */
	/* <kpoints> tag - if none, rewind and ignore */	
	if(fetch_in_file(vf,"<kpoints>")==0) fseek(vf,vfpos,SEEK_SET);
	else vasp_xml_read_kpoints(vf,model);
	/* <parameters> tag - if none, rewind and ignore */
	if(fetch_in_file(vf,"<parameters>")==0) fseek(vf,vfpos,SEEK_SET);
	else vasp_xml_read_incar(vf,model);
	/* <atominfo> tag - mandatory */
	if(fetch_in_file(vf,"<atominfo>")==0) return 3;
	if(vasp_xml_read_atominfo(vf,model)<0) return 3;
	/* <structure name="initialpos" > tag */
vfpos=ftell(vf);/*flag*/
	/* Counting the # of frames */
	if(fetch_in_file(vf,"initialpos")==0) return 3;
	line = file_read_line(vf);
	while (line){
		if (find_in_string("/calculation",line) != NULL) {
			add_frame_offset(vf, model);
			num_frames++;
		}
		if (find_in_string("/modeling",line) != NULL) break;
		g_free(line);
		line = file_read_line(vf);
	}
	g_free(line);
	fseek(vf,vfpos,SEEK_SET);/* rewind to flag */
	vasp_xml_read_frequency(vf,model);/* Read Frequency */
	fseek(vf,vfpos,SEEK_SET);/* rewind to flag */
	vasp_xml_read_dos(vf,model);/* Read DOS */
	vasp_xml_read_bands(vf,model);/* Read bands */
	fseek(vf,vfpos,SEEK_SET);/* rewind to flag */
	if (num_frames == 1){
		model->animation=FALSE;
		if(vasp_xml_read_pos(vf,model)<0) return 3;
	} else { /* we have num_frames-1 frames */
		model->cur_frame=1;
		if (num_frames <= 2){/* special cases, only initialpos and finalpos or single point calculation */
			num_frames=1;
			model->animation=FALSE;/* get rid of frame display */
			if(fetch_in_file(vf,"finalpos")==0) {
				fseek(vf,vfpos,SEEK_SET);/* rewind to flag */
				/*display at least the initial position, with a warning*/
				if(fetch_in_file(vf,"initialpos")==0) return 3;
				line = g_strdup_printf("WARNING: incomplete calculation!\n");
				gui_text_show(ERROR, line);
				g_free(line);
			}
			if(vasp_xml_read_pos(vf,model)<0) return 3;
		} else {
			model->animation=TRUE;
			num_frames=num_frames-2;/* because the finalpos is already part of the ionic steps */
			if(fetch_in_file(vf,"finalpos")==0) {
				/*imcomplete calculation, go the the last valid <structure>*/
				fseek(vf,vfpos,SEEK_SET);/* rewind to flag */
				while(fetch_in_file(vf,"<structure>")!=0) vfpos=ftell(vf);/*flag*/
				fseek(vf,vfpos,SEEK_SET);/* rewind to flag */
				line = g_strdup_printf("WARNING: incomplete calculation!\n");
				gui_text_show(ERROR, line);
				g_free(line);
			}
			if(vasp_xml_read_pos(vf,model)<0) return 3;
			model->cur_frame = num_frames - 1;
		}
	}
	/* at the end of file, or </modeling> tag */
	model->num_frames = num_frames;
	model->redraw = TRUE;
	/* always show this information */
	gui_text_show(ITALIC,g_strdup_printf("-> %i frames detected.\n",num_frames));
	strcpy(model->filename, filename);/* strcpy is ok? */
	fflush(stdout);
	model_prep(model);
	error_table_print_all();
	fclose(vf);
	return(0);
}
/* simplified frame reading */
gint read_xml_vasp_frame(FILE *vf, struct model_pak *model){
	g_assert(vf != NULL);
	if(fetch_in_file(vf,"scstep")==0) return 3;
	if(vasp_xml_read_pos(vf,model)<0) return 3;
	return 0;
}

