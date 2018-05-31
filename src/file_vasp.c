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
#include "edit.h"
#include "error.h"
#include "file.h"
#include "parse.h"
#include "matrix.h"
#include "zmatrix.h"
#include "zmatrix_pak.h"
#include "model.h"
#include "interface.h"
#include "file_vasp.h"

enum {VASP_DEFAULT, VASP_XML, VASP_TRIKS, VASP_NLOOPS};

/* main structures */
extern struct sysenv_pak sysenv;
extern struct elem_pak elements[];

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
/**********************************************************/
/* populate vasp_calc_struct with values from vasprun.xml */
/**********************************************************/
int vasp_xml_load_calc(FILE *vf,vasp_calc_struct *calc){
#define CALC (*calc)
#define DEBUG_XML2CALC 0
#define XC_REG_DOUBLE(field,value) do{\
        if (find_in_string(field,line) != NULL) {\
                sscanf(line," <i name=\"%*[^\"]\"> %[^<]</i>",tamp);\
                value=g_ascii_strtod(tamp,NULL);\
        }\
}while(0)
#define XC_REG_INT(field,value) do{\
	if (find_in_string(field,line) != NULL) {\
		sscanf(line," <i type=\"int\" name=\"%*[^\"]\"> %[^<]</i>",tamp);\
		value=(gint)g_ascii_strtoll(tamp,NULL,10);\
	}\
}while(0)
#define XC_REG_BOOL(field,value) do{\
	if (find_in_string(field,line) != NULL) {\
		sscanf(line," <i type=\"logical\" name=\"%*[^\"]\"> %s",tamp);\
		value=(tamp[0]=='T');\
	}\
}while(0)
#define XC_REG_TEXT(field,value) do{\
	if (find_in_string(field,line) != NULL) {\
		if(value!=NULL) g_free(value);\
		sscanf(line," <v name=\"%*[^\"]\"> %[^<]</v>",tamp);\
		value=g_strdup_printf("%s",tamp);\
	}\
}while(0)

	gchar *line;
	gchar tamp[512];/*seems large but some double arrays can become quite huge (magmom)*/
	if(vf==NULL) return -1;
	if(calc==NULL) return -1;
	/* first rewind*/
	rewind(vf);
	if(fetch_in_file(vf,"parameters")==0) return -1;
	/* we are now in parameters: general*/
	line = file_read_line(vf);
	while(line){
		if (find_in_string("</separator>",line) != NULL) break;/*end of general*/
		if (find_in_string("SYSTEM",line) != NULL){
			/*register SYSTEM*/
			if(CALC.name!=NULL) g_free(CALC.name);
			sscanf(line," <i type=\"string\" name=\"SYSTEM\"> %[^<]</i>",tamp);
			CALC.name=g_strdup_printf("%s",tamp);
#if DEBUG_XML2CALC
	fprintf(stdout,"#DBG XML2CALC: SYSTEM=%s\n",CALC.name);
#endif
		}
		//Skipped: LCOMPAT
		g_free(line);
		line = file_read_line(vf);
	}
	if(line==NULL) return -1;
	line = file_read_line(vf);/*electronic*/
	while(line){
		if (find_in_string("</separator>",line) != NULL) break;/*end of electronic*/
		if (find_in_string("PREC",line) != NULL) {
			/*register PREC*/
			sscanf(line," <i type=\"string\" name=\"PREC\">%s",tamp);
			tamp[0]=g_ascii_tolower(tamp[0]);
			switch (tamp[0]){
			case 's':
				CALC.prec=VP_SINGLE;break;
			case 'a':
				CALC.prec=VP_ACCURATE;break;
			case 'h':
				CALC.prec=VP_HIGH;break;
			case 'm':
				CALC.prec=VP_MED;break;
			case 'l':
				CALC.prec=VP_LOW;break;
			case 'n':
			default:
				CALC.prec=VP_NORM;
			}
#if DEBUG_XML2CALC
        fprintf(stdout,"#DBG XML2CALC: PREC=%c\n",tamp[0]);
#endif
		}
		XC_REG_DOUBLE("ENMAX",CALC.encut);
		XC_REG_DOUBLE("ENAUG",CALC.enaug);
		XC_REG_DOUBLE("EDIFF",CALC.ediff);
		XC_REG_DOUBLE("NELECT",CALC.nelect);
		//Skipped: EREF (unknown)
		XC_REG_DOUBLE("SIGMA",CALC.sigma);
		XC_REG_DOUBLE("KSPACING",CALC.kspacing);
		XC_REG_INT("IALGO",CALC.ialgo);
		XC_REG_INT("IWAVPR",CALC.iwavpr);
		XC_REG_INT("NBANDS",CALC.nbands);
		//Skipped: TURBO (unknown)
		//Skipped: IRESTART (unknown)
		//Skipped: NREBOOT (unknown)
		//Skipped: NMIN (unknown)
		XC_REG_INT("ISMEAR",CALC.ismear);
		XC_REG_BOOL("KGAMMA",CALC.kgamma);
		/*read next line*/
		g_free(line);
		line = file_read_line(vf);
	}
#if DEBUG_XML2CALC
	fprintf(stdout,"#DBG XML2CALC: ENCUT=%lf\n",CALC.encut);
	fprintf(stdout,"#DBG XML2CALC: ENAUG=%lf\n",CALC.enaug);
	fprintf(stdout,"#DBG XML2CALC: EDIFF=%lE\n",CALC.ediff);
	fprintf(stdout,"#DBG XML2CALC: IALGO=%i\n",CALC.ialgo);
	fprintf(stdout,"#DBG XML2CALC: IWAVPR=%i\n",CALC.iwavpr);
	fprintf(stdout,"#DBG XML2CALC: NBANDS=%i\n",CALC.nbands);
	fprintf(stdout,"#DBG XML2CALC: NELECT=%lf\n",CALC.nelect);
	fprintf(stdout,"#DBG XML2CALC: ISMEAR=%i\n",CALC.ismear);
	fprintf(stdout,"#DBG XML2CALC: SIGMA=%lf\n",CALC.sigma);
	fprintf(stdout,"#DBG XML2CALC: KSPACING=%lf\n",CALC.kspacing);
if(CALC.kgamma) fprintf(stdout,"#DBG XML2CALC: KGAMMA=.TRUE.\n");
else 		fprintf(stdout,"#DBG XML2CALC: SIGMA=.FALSE.\n");
#endif
	line = file_read_line(vf);/*electronic projector*/
	while(line){
		if (find_in_string("</separator>",line) != NULL) break;/*end of electronic projector*/
		//LREAL: specia
		if (find_in_string("LREAL",line) != NULL) {/*LREAL is still a boolean in vasprun.xml*/
			sscanf(line," <i type=\"logical\" name=\"LREAL\"> %s",tamp);
			if(tamp[0]=='T') CALC.lreal=VLR_TRUE;
			else CALC.lreal=VLR_FALSE;
		}
		XC_REG_TEXT("ROPT",CALC.ropt);
		XC_REG_INT("LMAXPAW",CALC.lmaxpaw);
		XC_REG_INT("LMAXMIX",CALC.lmaxmix);
		//Skipped: NLSPLINE (unknown)
		/*read next line*/
		g_free(line);
		line = file_read_line(vf);
	}
#if DEBUG_XML2CALC
if(CALC.lreal==VLR_TRUE)	fprintf(stdout,"#DBG XML2CALC: LREAL=.TRUE.\n");
else				fprintf(stdout,"#DBG XML2CALC: LREAL=.FALSE.\n");
	fprintf(stdout,"#DBG XML2CALC: ROPT=%s\n",CALC.ropt);
	fprintf(stdout,"#DBG XML2CALC: LMAXPAW=%i\n",CALC.lmaxpaw);
	fprintf(stdout,"#DBG XML2CALC: LMAXMIX=%i\n",CALC.lmaxmix);
#endif
	line = file_read_line(vf);/*electronic startup*/
	while(line){
		if (find_in_string("</separator>",line) != NULL) break;/*end of electronic startup*/
		XC_REG_INT("ISTART",CALC.istart);
		XC_REG_INT("ICHARG",CALC.icharg);
		//INIWAV:special
		if (find_in_string("INIWAV",line) != NULL) {
			sscanf(line," <i type=\"int\" name=\"INIWAV\"> %[^<]</i>",tamp);
			CALC.iniwav=(g_ascii_digit_value(tamp[0])==1);
		}
		/*read next line*/
		g_free(line);
		line = file_read_line(vf);
	}
#if DEBUG_XML2CALC
	fprintf(stdout,"#DBG XML2CALC: ISTART=%i\n",CALC.istart);
	fprintf(stdout,"#DBG XML2CALC: ICHARG=%i\n",CALC.icharg);
if(CALC.iniwav) fprintf(stdout,"#DBG XML2CALC: INIWAV=1\n");
else 		fprintf(stdout,"#DBG XML2CALC: INIWAV=0\n");
#endif
	line = file_read_line(vf);/*electronic spin*/
	while(line){
		if (find_in_string("</separator>",line) != NULL) break;/*end of electronic spin*/
		//ISPIN: special 
		if (find_in_string("ISPIN",line) != NULL) {
			sscanf(line," <i type=\"int\" name=\"ISPIN\"> %[^<]</i>",tamp);
			CALC.ispin=(g_ascii_digit_value(tamp[0])==2);
		}
		XC_REG_BOOL("LNONCOLLINEAR",CALC.non_collinear);
		XC_REG_TEXT("MAGMOM",CALC.magmom);
		XC_REG_DOUBLE("NUPDOWN",CALC.nupdown);
		XC_REG_BOOL("LSORBIT",CALC.lsorbit);
		XC_REG_TEXT("SAXIS",CALC.saxis);
		//Skipped: LSPIRAL (unknown)
		//Skipped: QSPIRAL (unknown)
		//Skipped: LZEROZ (unknown)
		/*read next line*/
		g_free(line);
		line = file_read_line(vf);
	}
#if DEBUG_XML2CALC
if(CALC.ispin)		fprintf(stdout,"#DBG XML2CALC: ISPIN=2\n");
else 			fprintf(stdout,"#DBG XML2CALC: ISPIN=1\n");
if(CALC.non_collinear) 	fprintf(stdout,"#DBG XML2CALC: LNONCOLLINEAR=.TRUE.\n");
else 			fprintf(stdout,"#DBG XML2CALC: LNONCOLLINEAR=.FALSE.\n");
	fprintf(stdout,"#DBG XML2CALC: MAGMOM=%s\n",CALC.magmom);
	fprintf(stdout,"#DBG XML2CALC: NUPDOWN=%i\n",CALC.nupdown);
if(CALC.lsorbit)	fprintf(stdout,"#DBG XML2CALC: LSORBIT=.TRUE.\n");
else			fprintf(stdout,"#DBG XML2CALC: LSORBIT=.FALSE.\n");
	fprintf(stdout,"#DBG XML2CALC: SAXIS=%s\n",CALC.saxis);
#endif
	line = file_read_line(vf);/*electronic exchange-correlation*/
	while(line){
		if (find_in_string("</separator>",line) != NULL) break;/*end of electronic exchange-correlation*/
		XC_REG_BOOL("LASPH",CALC.lasph);
		XC_REG_BOOL("LMETAGGA",CALC.lmetagga);/*FIXME: when LMETAGGA=TRUE METAGGA and tags are give in INCAR ONLY*/
		/*read next line*/
		g_free(line);
		line = file_read_line(vf);
	}
#if DEBUG_XML2CALC
if(CALC.lasph)		fprintf(stdout,"#DBG XML2CALC: LASPH=.TRUE.\n");
else			fprintf(stdout,"#DBG XML2CALC: LASPH=.FALSE.\n");
if(CALC.lmetagga)	fprintf(stdout,"#DBG XML2CALC: LMETAGGA=.TRUE.\n");
else			fprintf(stdout,"#DBG XML2CALC: LMETAGGA=.FALSE.\n");
#endif
	line = file_read_line(vf);/*electronic convergence*/
	while(line){
		if (find_in_string("</separator>",line) != NULL) break;/*end of electronic convergence*/
		XC_REG_INT("NELM",CALC.nelm);
		XC_REG_INT("NELMDL",CALC.nelmdl);
		XC_REG_INT("NELMIN",CALC.nelmin);
		//Skipped: ENINI (unknown)
		XC_REG_BOOL("LDIAG",CALC.ldiag);
		//Skipped: SUBROT (adv.)
		//Skipped: WEIMIN (adv.)
		//Skipped: EBREAK (adv.)
		//Skipped: DEPER (adv.)
		//Skipped: NRMM (unknown)
		//TIME: special (because it is called vtime in vasp_calc_struct)
		if (find_in_string("TIME",line) != NULL) {
			sscanf(line," <i name=\"%*[^\"]\"> %[^<]</i>",tamp);
			CALC.vtime=g_ascii_strtod(tamp,NULL);
		}
		/*read next line*/
		g_free(line);
		line = file_read_line(vf);
	}
#if DEBUG_XML2CALC
	fprintf(stdout,"#DBG XML2CALC: NELM=%i\n",CALC.nelm);
	fprintf(stdout,"#DBG XML2CALC: NELMDL=%i\n",CALC.nelmdl);
	fprintf(stdout,"#DBG XML2CALC: NELMIN=%i\n",CALC.nelmin);
if(CALC.ldiag) 	fprintf(stdout,"#DBG XML2CALC: LDIAG=.TRUE.\n");
else		fprintf(stdout,"#DBG XML2CALC: LDIAG=.FALSE.\n");
	fprintf(stdout,"#DBG XML2CALC: TIME=%lf\n",CALC.vtime);
#endif
	line = file_read_line(vf);/*because there is 2 end </separator>*/
	g_free(line);
	line = file_read_line(vf);/*electronic mixer*/
	while(line){
		if (find_in_string("</separator>",line) != NULL) break;/*end of electronic mixer*/
		XC_REG_DOUBLE("AMIX",CALC.amix);
		XC_REG_DOUBLE("BMIX",CALC.bmix);
		XC_REG_DOUBLE("AMIN",CALC.amin);
		XC_REG_DOUBLE("AMIX_MAG",CALC.amix_mag);
		XC_REG_DOUBLE("BMIX_MAG",CALC.bmix_mag);
		XC_REG_INT("IMIX",CALC.imix);
		//Skipped: MIXFIRST (unknown)
		XC_REG_INT("MAXMIX",CALC.maxmix);
		XC_REG_DOUBLE("WC",CALC.wc);
		XC_REG_INT("INIMIX",CALC.inimix);
		XC_REG_INT("MIXPRE",CALC.mixpre);
		//Skipped: MREMOVE (unknown)
		/*read next line*/
		g_free(line);
		line = file_read_line(vf);
	}
#if DEBUG_XML2CALC
	fprintf(stdout,"#DBG XML2CALC: AMIX=%lf\n",CALC.amix);
	fprintf(stdout,"#DBG XML2CALC: BMIX=%lf\n",CALC.bmix);
	fprintf(stdout,"#DBG XML2CALC: AMIN=%lf\n",CALC.amin);
	fprintf(stdout,"#DBG XML2CALC: AMIX_MAG=%lf\n",CALC.amix_mag);
	fprintf(stdout,"#DBG XML2CALC: BMIX_MAG=%lf\n",CALC.bmix_mag);
	fprintf(stdout,"#DBG XML2CALC: IMIX=%i\n",CALC.imix);
	fprintf(stdout,"#DBG XML2CALC: MAXMIX=%i\n",CALC.maxmix);
	fprintf(stdout,"#DBG XML2CALC: WC=%lf\n",CALC.wc);
	fprintf(stdout,"#DBG XML2CALC: INIMIX=%i\n",CALC.inimix);
	fprintf(stdout,"#DBG XML2CALC: MIXPRE=%i\n",CALC.mixpre);
#endif
	line = file_read_line(vf);/*because there is 2 end </separator>*/
	g_free(line);
	line = file_read_line(vf);/*electronic dipolcorrection*/
	while(line){
		if (find_in_string("</separator>",line) != NULL) break;/*end of electronic dipolcorrection*/
		XC_REG_BOOL("LDIPOL",CALC.ldipol);
		XC_REG_BOOL("LMONO",CALC.lmono);
		XC_REG_INT("IDIPOL",CALC.idipol);
		XC_REG_DOUBLE("EPSILON",CALC.epsilon);
		XC_REG_TEXT("DIPOL",CALC.dipol);
		XC_REG_DOUBLE("EFIELD",CALC.efield);
		/*read next line*/
		g_free(line);
		line = file_read_line(vf);
	}
#if DEBUG_XML2CALC
if(CALC.ldipol) fprintf(stdout,"#DBG XML2CALC: LDIPOL=.TRUE.\n");
else 		fprintf(stdout,"#DBG XML2CALC: LDIPOL=.FALSE.\n");
if(CALC.lmono) 	fprintf(stdout,"#DBG XML2CALC: LMONO=.TRUE.\n");
else 		fprintf(stdout,"#DBG XML2CALC: LMONO=.FALSE.\n");
	fprintf(stdout,"#DBG XML2CALC: IDIPOL=%i\n",CALC.idipol);
	fprintf(stdout,"#DBG XML2CALC: EPSILON=%lf\n",CALC.epsilon);
	fprintf(stdout,"#DBG XML2CALC: DIPOL=%s\n",CALC.dipol);
	fprintf(stdout,"#DBG XML2CALC: EFIELD=%lf\n",CALC.efield);
#endif
	line = file_read_line(vf);/*because there is 2 end </separator>*/
	g_free(line);
	line = file_read_line(vf);/*grids*/
	while(line){
		if (find_in_string("</separator>",line) != NULL) break;/*end of grids*/
		XC_REG_INT("NGX\"",CALC.ngx);
		XC_REG_INT("NGY\"",CALC.ngy);
		XC_REG_INT("NGZ\"",CALC.ngz);
		XC_REG_INT("NGXF",CALC.ngxf);
		XC_REG_INT("NGYF",CALC.ngyf);
		XC_REG_INT("NGZF",CALC.ngzf);
		XC_REG_BOOL("ADDGRID",CALC.addgrid);
		/*read next line*/
		g_free(line);
		line = file_read_line(vf);
	}
#if DEBUG_XML2CALC
	fprintf(stdout,"#DBG XML2CALC: NGX=%i\n",CALC.ngx);
	fprintf(stdout,"#DBG XML2CALC: NGY=%i\n",CALC.ngy);
	fprintf(stdout,"#DBG XML2CALC: NGZ=%i\n",CALC.ngz);
	fprintf(stdout,"#DBG XML2CALC: NGXF=%i\n",CALC.ngxf);
	fprintf(stdout,"#DBG XML2CALC: NGYF=%i\n",CALC.ngyf);
	fprintf(stdout,"#DBG XML2CALC: NGZF=%i\n",CALC.ngzf);
if(CALC.addgrid) 	fprintf(stdout,"#DBG XML2CALC: ADDGRID=.TRUE.\n");
else 			fprintf(stdout,"#DBG XML2CALC: ADDGRID=.FALSE.\n");
#endif
	line = file_read_line(vf);/*ionic*/
	while(line){
		if (find_in_string("</separator>",line) != NULL) break;/*end of ionic*/
		XC_REG_INT("NSW",CALC.nsw);
		XC_REG_INT("IBRION\"",CALC.ibrion);
		XC_REG_INT("ISIF",CALC.isif);
		XC_REG_DOUBLE("PSTRESS",CALC.pstress);
		XC_REG_DOUBLE("EDIFFG",CALC.ediffg);
		XC_REG_INT("NFREE",CALC.nfree);
		XC_REG_DOUBLE("POTIM",CALC.potim);
		XC_REG_DOUBLE("SMASS",CALC.smass);
		//Skipped: SCALEE (unknown)
		/*read next line*/
		g_free(line);
		line = file_read_line(vf);
	}
#if DEBUG_XML2CALC
	fprintf(stdout,"#DBG XML2CALC: NSW=%i\n",CALC.nsw);
	fprintf(stdout,"#DBG XML2CALC: IBRION=%i\n",CALC.ibrion);
	fprintf(stdout,"#DBG XML2CALC: ISIF=%i\n",CALC.isif);
	fprintf(stdout,"#DBG XML2CALC: PSTRESS=%lf\n",CALC.pstress);
	fprintf(stdout,"#DBG XML2CALC: EDIFFG=%lE\n",CALC.ediffg);
	fprintf(stdout,"#DBG XML2CALC: NFREE=%i\n",CALC.nfree);
	fprintf(stdout,"#DBG XML2CALC: POTIM=%lf\n",CALC.potim);
	fprintf(stdout,"#DBG XML2CALC: SMASS=%lf\n",CALC.smass);
#endif
	line = file_read_line(vf);/*ionic md*/
	while(line){
		if (find_in_string("</separator>",line) != NULL) break;/*end of ionic md*/
		XC_REG_DOUBLE("TEBEG",CALC.tebeg);
		XC_REG_DOUBLE("TEEND",CALC.teend);
		XC_REG_INT("NBLOCK",CALC.nblock);
		XC_REG_INT("KBLOCK",CALC.kblock);
		XC_REG_INT("NPACO",CALC.npaco);
		XC_REG_DOUBLE("APACO",CALC.apaco);
		/*read next line*/
		g_free(line);
		line = file_read_line(vf);
	}
#if DEBUG_XML2CALC
	fprintf(stdout,"#DBG XML2CALC: TEBEG=%lf\n",CALC.tebeg);
	fprintf(stdout,"#DBG XML2CALC: TEEND=%lf\n",CALC.teend);
	fprintf(stdout,"#DBG XML2CALC: NBLOCK=%i\n",CALC.nblock);
	fprintf(stdout,"#DBG XML2CALC: KBLOCK=%i\n",CALC.kblock);
	fprintf(stdout,"#DBG XML2CALC: NPACO=%i\n",CALC.npaco);
	fprintf(stdout,"#DBG XML2CALC: APACO=%lf\n",CALC.apaco);
#endif
	line = file_read_line(vf);/*symmetry*/
	while(line){
		if (find_in_string("</separator>",line) != NULL) break;/*end of symmetry*/
		XC_REG_INT("ISYM",CALC.isym);
		XC_REG_DOUBLE("SYMPREC",CALC.sym_prec);
		/*read next line*/
		g_free(line);
		line = file_read_line(vf);
	}
#if DEBUG_XML2CALC
	fprintf(stdout,"#DBG XML2CALC: ISYM=%i\n",CALC.isym);
	fprintf(stdout,"#DBG XML2CALC: SYMPREC=%lE\n",CALC.sym_prec);
#endif
	line = file_read_line(vf);/*dos*/
	while(line){
		if (find_in_string("</separator>",line) != NULL) break;/*end of dos*/
		XC_REG_BOOL("LORBIT",CALC.lorbit);
		XC_REG_TEXT("RWIGS",CALC.rwigs);
		XC_REG_INT("NEDOS",CALC.nedos);
		XC_REG_DOUBLE("EMIN",CALC.emin);
		XC_REG_DOUBLE("EMAX",CALC.emax);
		XC_REG_DOUBLE("EFERMI",CALC.efermi);//actually (unknown)
		/*read next line*/
		g_free(line);
		line = file_read_line(vf);
	}
#if DEBUG_XML2CALC
if(CALC.lorbit) fprintf(stdout,"#DBG XML2CALC: LORBIT=.TRUE.\n");
else 		fprintf(stdout,"#DBG XML2CALC: LORBIT=.FALSE.\n");
	fprintf(stdout,"#DBG XML2CALC: RWIGS=%s\n",CALC.rwigs);
	fprintf(stdout,"#DBG XML2CALC: NEDOS=%lf\n",CALC.nedos);
	fprintf(stdout,"#DBG XML2CALC: EMIN=%lf\n",CALC.emin);
	fprintf(stdout,"#DBG XML2CALC: EMAX=%lf\n",CALC.emax);
	fprintf(stdout,"#DBG XML2CALC: EFERMI=%lf\n",CALC.efermi);
#endif
	line = file_read_line(vf);/*writing*/
	while(line){
		if (find_in_string("</separator>",line) != NULL) break;/*end of writing*/
		XC_REG_INT("NWRITE",CALC.nwrite);
		XC_REG_BOOL("LWAVE",CALC.lwave);
		XC_REG_BOOL("LCHARG",CALC.lcharg);
		//Skipped: LPARD (adv.)
		XC_REG_BOOL("LVTOT",CALC.lvtot);
		XC_REG_BOOL("LVHAR",CALC.lvhar);
		XC_REG_BOOL("LELF",CALC.lelf);
		XC_REG_BOOL("LOPTICS",CALC.loptics);/*why here?*/
		//Skipped: STM (adv.)
		/*read next line*/
		g_free(line);
		line = file_read_line(vf);
	}
#if DEBUG_XML2CALC
	fprintf(stdout,"#DBG XML2CALC: NWRITE=%i\n",CALC.nwrite);
if(CALC.lwave) 	fprintf(stdout,"#DBG XML2CALC: LWAVE=.TRUE.\n");
else 		fprintf(stdout,"#DBG XML2CALC: LWAVE=.FALSE.\n");
if(CALC.lcharg)	fprintf(stdout,"#DBG XML2CALC: LCHARG=.TRUE.\n");
else		fprintf(stdout,"#DBG XML2CALC: LCHARG=.FALSE.\n");
if(CALC.lvtot) 	fprintf(stdout,"#DBG XML2CALC: LVTOT=.TRUE.\n");
else		fprintf(stdout,"#DBG XML2CALC: LVTOT=.FALSE.\n");
if(CALC.lvhar) 	fprintf(stdout,"#DBG XML2CALC: LVHAR=.TRUE.\n");
else		fprintf(stdout,"#DBG XML2CALC: LVHAR=.FALSE.\n");
if(CALC.lelf) 	fprintf(stdout,"#DBG XML2CALC: LELF=.TRUE.\n");
else		fprintf(stdout,"#DBG XML2CALC: LELF=.FALSE.\n");
if(CALC.loptics)
		fprintf(stdout,"#DBG XML2CALC: LOPTICS=.TRUE.\n");
else		fprintf(stdout,"#DBG XML2CALC: LOPTICS=.FALSE.\n");
#endif
	line = file_read_line(vf);/*performance*/
	while(line){
		if (find_in_string("</separator>",line) != NULL) break;/*end of performance*/
		XC_REG_INT("NPAR",CALC.npar);
		XC_REG_INT("NSIM",CALC.nsim);
		//Skipped: NBLK (adv.)
		XC_REG_BOOL("LPLANE",CALC.lplane);
		XC_REG_BOOL("LSCALAPACK",CALC.lscalapack);
		//Skipped: LSCAAWARE (unknown)
		XC_REG_BOOL("LSCALU",CALC.lscalu);
		//Skipped: LASYNC (unknown)
		//Skipped: LORBITALREAL (unknown)
		/*read next line*/
		g_free(line);
		line = file_read_line(vf);
	}
	
#if DEBUG_XML2CALC
	fprintf(stdout,"#DBG XML2CALC: NPAR=%i\n",CALC.npar);
	fprintf(stdout,"#DBG XML2CALC: NSIM=%i\n",CALC.nsim);
if(CALC.lplane)	fprintf(stdout,"#DBG XML2CALC: LPLANE=.TRUE.\n");
else		fprintf(stdout,"#DBG XML2CALC: LPLANE=.FALSE.\n");
if(CALC.lscalapack)
		fprintf(stdout,"#DBG XML2CALC: LSCALAPACK=.TRUE.\n");
else		fprintf(stdout,"#DBG XML2CALC: LSCALAPACK=.FALSE.\n");
if(CALC.lscalu)	fprintf(stdout,"#DBG XML2CALC: LSCALU=.TRUE.\n");
else		fprintf(stdout,"#DBG XML2CALC: LSCALU=.FALSE.\n");
#endif
	line = file_read_line(vf);/*miscellaneous*/
	while(line){/*skip the section*/
		if (find_in_string("</separator>",line) != NULL) break;/*end of miscellaneous*/
		/*read next line*/
		g_free(line);
		line = file_read_line(vf);
	}
	line = file_read_line(vf);/*UNNAMED section (ie. we are reading GGA_COMPAT already)*/
	while(line){
		if (find_in_string("<separator",line) != NULL) break;
			/*end of UNNAMED section: the section does not have ending separator*/
		XC_REG_BOOL("GGA_COMPAT",CALC.gga_compat);
		//Skipped: LBERRY (adv.)
		//Skipped: ICORELEVEL (adv.)
		XC_REG_BOOL("LDAU",CALC.ldau);
		//Skipped: I_CONSTRAINED_M (unknown)
		/*read next line*/
		g_free(line);
		line = file_read_line(vf);
	}
#if DEBUG_XML2CALC
if(CALC.gga_compat)
		fprintf(stdout,"#DBG XML2CALC: GGA_COMPAT=.TRUE.\n");
else		fprintf(stdout,"#DBG XML2CALC: GGA_COMPAT=.FALSE.\n");
if(CALC.ldau)	fprintf(stdout,"#DBG XML2CALC: LDAU=.TRUE.\n");
else		fprintf(stdout,"#DBG XML2CALC: LDAU=.FALSE.\n");
#endif
	/*jump directly to next section: electronic exchange-correlation*/
	while(line){
		if (find_in_string("</separator>",line) != NULL) break;/*end of electronic exchange-correlation*/
		//GGA: special
                if (find_in_string("GGA\"",line) != NULL) {
                        /*register PREC*/
                        sscanf(line," <i type=\"string\" name=\"GGA\">%s",tamp);
                        tamp[0]=g_ascii_tolower(tamp[0]);
			tamp[1]=g_ascii_tolower(tamp[1]);
                        switch (tamp[0]){
                        case '9':
                                CALC.gga=VG_91;break;
                        case 'R':
                                CALC.gga=VG_RP;break;
                        case 'A':
                                CALC.gga=VG_AM;break;
                        case 'P':
                        default:
				CALC.gga=VG_PE;
				if(tamp[1]=='s') CALC.gga=VG_PS;
                        }
#if DEBUG_XML2CALC
        fprintf(stdout,"#DBG XML2CALC: GGA=%c%c\n",tamp[0],tamp[1]);
#endif
                }
		XC_REG_BOOL("VOSKOWN",CALC.voskown);
		//Everything else is skipped and reserved for (adv.)
		/*read next line*/
		g_free(line);
		line = file_read_line(vf);

	}	
#if DEBUG_XML2CALC
if(CALC.voskown)
		fprintf(stdout,"#DBG XML2CALC: VOSKOWN=.TRUE.\n");
else		fprintf(stdout,"#DBG XML2CALC: VOSKOWN=.FALSE.\n");
#endif
	line = file_read_line(vf);/*vdW DFT*/
	while(line){
		if (find_in_string("</separator>",line) != NULL) break;/*end of vdW DFT*/
		XC_REG_BOOL("LUSE_VDW",CALC.use_vdw);
		XC_REG_DOUBLE("Zab_VDW",CALC.zab_vdw);
		XC_REG_DOUBLE("PARAM1",CALC.param1_vdw);
		XC_REG_DOUBLE("PARAM2",CALC.param2_vdw);
		XC_REG_DOUBLE("PARAM3",CALC.param3_vdw);
		/*read next line*/
		g_free(line);
		line = file_read_line(vf);
	}
#if DEBUG_XML2CALC
if(CALC.use_vdw)
		fprintf(stdout,"#DBG XML2CALC: LUSE_VDW=.TRUE.\n");
else		fprintf(stdout,"#DBG XML2CALC: LUSE_VDW=.FALSE.\n");
	fprintf(stdout,"#DBG XML2CALC: Zab_VDW=%lf\n",CALC.zab_vdw);
	fprintf(stdout,"#DBG XML2CALC: PARAM1=%lf\n",CALC.param1_vdw);
	fprintf(stdout,"#DBG XML2CALC: PARAM2=%lf\n",CALC.param2_vdw);
	fprintf(stdout,"#DBG XML2CALC: PARAM3=%lf\n",CALC.param3_vdw);
#endif
	line = file_read_line(vf);/*model GW*/
	while(line){
		if (find_in_string("</separator>",line) != NULL) break;/*end of model GW*/
		//Everything is skipped and reserved for (adv.)
		/*read next line*/
		g_free(line);
		line = file_read_line(vf);
	}
	line = file_read_line(vf);/*linear response parameters*/
	while(line){
		if (find_in_string("</separator>",line) != NULL) break;/*end of linear response parameters*/
		XC_REG_BOOL("LEPSILON",CALC.lepsilon);
		XC_REG_BOOL("LRPA",CALC.lrpa);
		XC_REG_BOOL("LNABLA",CALC.lnabla);
		//Skipped: LVEL (unknown)
		//Skipped: KINTER (unknown)
		XC_REG_DOUBLE("CSHIFT",CALC.cshift);
		//Skipped: OMEGAMAX (unknown)
		//Skipped: DEG_THRESHOLD (unknown)
		/*read next line*/
		g_free(line);
		line = file_read_line(vf);
	}
#if DEBUG_XML2CALC
if(CALC.lepsilon)
		fprintf(stdout,"#DBG XML2CALC: LEPSILON=.TRUE.\n");
else		fprintf(stdout,"#DBG XML2CALC: LEPSILON=.FALSE.\n");
if(CALC.lrpa)	fprintf(stdout,"#DBG XML2CALC: LRPA=.TRUE.\n");
else 		fprintf(stdout,"#DBG XML2CALC: LRPA=.FALSE.\n");
if(CALC.lnabla)	fprintf(stdout,"#DBG XML2CALC: LNABLA=.TRUE.\n");
else		fprintf(stdout,"#DBG XML2CALC: LNABLA=.FALSE.\n");
	fprintf(stdout,"#DBG XML2CALC: CSHIFT=%lf\n",CALC.cshift);
#endif

	//Everything else is skipped and reserved for (adv.)

	/* TODO (adv.)*/

	return 0;
#undef CALC
}
/****************************************************/
/* setup vasp_calc_struct (using file if available) */
/****************************************************/
gint vasprun_update(gchar *filename,vasp_calc_struct *calc){
#define CALC (*calc)
	struct model_pak *data=sysenv.active_model;
        GSList *list2;
        struct core_pak *core;
        /*Common setup for xml or regular models*/
/*4-POSCAR*/
        /* unexposed */
        CALC.species_symbols=NULL;
        CALC.species_numbers=NULL;
        /* exposed */
        CALC.poscar_x=0.;
        CALC.poscar_y=0.;
        CALC.poscar_z=0.;
	/* from model */
/* This fix a BUG specific to QEOUT, where the last position can have fractional=FALSE, while coordinate is fractional*/
if(data->id==QE_OUT){
	FILE *fp=fopen(data->filename, "r");
        if (!fp){
                gui_text_show(ERROR,g_strdup_printf("I/O ERROR: can't open .qeout file!\n"));
                return -1;
        }
        read_raw_frame(fp,data->cur_frame,data);
        fclose(fp);
}
	CALC.poscar_a0=1.0;
	if(data->fractional) CALC.poscar_direct=TRUE;
	else CALC.poscar_direct=FALSE;
	if(data->periodic<1) {/*not periodic: create a big cubic box*/
		gdouble xmax=0,ymax=0,zmax=0;
		gdouble xmin=0,ymin=0,zmin=0;
		for (list2=data->cores ; list2 ; list2=g_slist_next(list2)){
			core=list2->data;
			if(xmax<(core->x[0])) xmax=core->x[0];
			if(ymax<(core->x[1])) ymax=core->x[1];
			if(zmax<(core->x[2])) zmax=core->x[2];
			if(xmin>(core->x[0])) xmin=core->x[0];
			if(ymin>(core->x[1])) ymin=core->x[1];
			if(zmin>(core->x[2])) zmin=core->x[2];
		}
		CALC.poscar_ux=10.+(xmax-xmin);
		CALC.poscar_uy=0.;
		CALC.poscar_uz=0.;
		CALC.poscar_vx=0.;
		CALC.poscar_vy=10.+(ymax-ymin);
		CALC.poscar_vz=0.;
		CALC.poscar_wx=0.;
		CALC.poscar_wy=0.;
		CALC.poscar_wz=10.+(zmax-zmin);
		CALC.poscar_direct=FALSE;/*always for 0D*/
	} else if(data->periodic<2){/* 1D-periodic: add 10. unit interspace */
		gdouble xmax=0,ymax=0;
		gdouble xmin=0,ymin=0;
		for (list2=data->cores ; list2 ; list2=g_slist_next(list2)){
			core=list2->data;
			if(xmax<(core->x[0])) xmax=core->x[0];
			if(ymax<(core->x[1])) ymax=core->x[1];
			if(xmin>(core->x[0])) xmin=core->x[0];
			if(ymin>(core->x[1])) ymin=core->x[1];

		}
		for (list2=data->cores ; list2 ; list2=g_slist_next(list2)){
			core=list2->data;/*TODO: confirm 1D behavior (find example)*/
			core->x[0]=core->x[0]/(10.+(xmax-xmin));
			core->x[1]=core->x[1]/(10.+(ymax-ymin));
		}
		CALC.poscar_ux=10.+(xmax-xmin);
		CALC.poscar_uy=0.;
		CALC.poscar_uz=0.;
		CALC.poscar_vx=0.;
		CALC.poscar_vy=10.+(ymax-ymin);
		CALC.poscar_vz=0.;
		CALC.poscar_wx=0.;
		CALC.poscar_wy=0.;
		CALC.poscar_wz=data->latmat[0];
	} else if(data->periodic<3){/* 2D-periodic: add a 10. unit vacuum layer */
		gdouble zmax=0.;
		gdouble zmin=0.;
		for (list2=data->cores ; list2 ; list2=g_slist_next(list2)){
			core=list2->data;
			if(zmax<(core->x[2])) zmax=core->x[2];
			if(zmin>(core->x[2])) zmin=core->x[2];
		}
		for (list2=data->cores ; list2 ; list2=g_slist_next(list2)){
			core=list2->data;
			core->x[2]=core->x[2]/(10.+(zmax-zmin));
		}
		CALC.poscar_ux=data->latmat[0];
		CALC.poscar_uy=data->latmat[3];
		CALC.poscar_uz=data->latmat[6];
		CALC.poscar_vx=data->latmat[1];
		CALC.poscar_vy=data->latmat[4];
		CALC.poscar_vz=data->latmat[7];
		CALC.poscar_wx=0.;
		CALC.poscar_wy=0.;
		CALC.poscar_wz=10.+(zmax-zmin);
	} else {/* 3D-periodic: Nothing much to do */
		CALC.poscar_ux=data->latmat[0];
		CALC.poscar_uy=data->latmat[3];
		CALC.poscar_uz=data->latmat[6];
		CALC.poscar_vx=data->latmat[1];
		CALC.poscar_vy=data->latmat[4];
		CALC.poscar_vz=data->latmat[7];
		CALC.poscar_wx=data->latmat[2];
		CALC.poscar_wy=data->latmat[5];
		CALC.poscar_wz=data->latmat[8];
	}
	CALC.poscar_sd=TRUE;
/*5-KPOINTS*/
        CALC.kpoints_nkpts=0;
/*...*/
        CALC.tetra_a=1;
        CALC.tetra_b=1;
        CALC.tetra_c=1;
        CALC.tetra_d=1;
/*6-POTCAR*/
        CALC.potcar_folder=NULL;
        CALC.potcar_file=NULL;
        CALC.potcar_species=NULL;
        CALC.potcar_species_flavor=NULL;
/*7-RUN*/
        CALC.job_vasp_exe=g_strdup_printf("%s",sysenv.vasp_path);
        CALC.job_mpirun=g_strdup_printf("%s",sysenv.mpirun_path);
        CALC.job_path=g_strdup_printf("%s",sysenv.cwd);
        CALC.job_nproc=1.;
        CALC.ncore=1.;
        CALC.kpar=1.;

/*specific part starts here*/

	if(filename){

		/*load most values from provided file*/
		/*TODO: add vaspxml reader*/
		FILE *fp=fopen(filename,"r");
		if(!fp) vasprun_update(NULL,calc);
		if(vasp_xml_load_calc(fp,calc)!=0) vasprun_update(NULL,calc);
		fclose(fp);
		CALC.use_prec=TRUE;
		CALC.auto_elec=TRUE;
		CALC.auto_mixer=TRUE;
		CALC.auto_grid=TRUE;

		/*POTCAR is unrealated to vasprun.xml*/
                CALC.potcar_folder=NULL;
                CALC.potcar_file=NULL;
                CALC.potcar_species=NULL;
                CALC.potcar_species_flavor=NULL;

		/*this might be caught in INCAR*/
		CALC.algo=VA_IALGO;
		CALC.ncore=1.;
		CALC.kpar=1.;
	}else{
		/*use some default values*/
/*name*/
		CALC.name=NULL;
/*1-general*/
		CALC.prec=VP_NORM;
		CALC.use_prec=TRUE;
		CALC.encut=0.0;/*ie not set*/
		CALC.enaug=0.0;/*ie not set*/
		CALC.ediff=1e-4;
		CALC.algo=VA_NORM;
		CALC.ialgo=VIA_KOSUGI;
		CALC.ldiag=TRUE;
		CALC.auto_elec=TRUE;
		CALC.nbands=0;/*ie not set*/
		CALC.nelect=0.0;/*ie not set*/
		CALC.iwavpr=10;
		CALC.nsim=4;
		CALC.vtime=0.4;
/*2-smearing*/
		CALC.ismear=1;
		CALC.sigma=0.2;
		CALC.fermwe=NULL;
		CALC.fermdo=NULL;
		CALC.kspacing=0.5;
		CALC.kgamma=TRUE;
/*3-projector*/
	        CALC.lreal=VLR_FALSE;
	        CALC.ropt=NULL;
	        CALC.addgrid=FALSE;
	        CALC.lmaxmix=2;
	        CALC.lmaxpaw=-100;/*ie don't show*/
/*4-startup*/
		CALC.istart=-1;/*ie not set*/
		CALC.icharg=-1;/*ie not set*/
		CALC.iniwav=TRUE;
/*5-spin*/
		CALC.ispin=FALSE;
		CALC.non_collinear=FALSE;
		CALC.magmom=NULL;
		CALC.nupdown=-1.0;/*ie not set*/
		CALC.lsorbit=FALSE;
		CALC.saxis=NULL;
		CALC.gga_compat=TRUE;
/*6-xc*/
		CALC.gga=VG_PE;
		CALC.voskown=FALSE;
		CALC.lasph=FALSE;
		CALC.lmaxtau=0;
		CALC.lmixtau=FALSE;
		CALC.lmetagga=FALSE;
		CALC.mgga=VMG_MBJ;/*default is none actually*/
		CALC.cmbj=NULL;
		CALC.cmbja=-0.012;
		CALC.cmbjb=1.023;
		/*vdw part*/
/*7-convergence*/
		CALC.nelm=60;
		CALC.nelmdl=-12;
		CALC.nelmin=2;
/*8-mixer*/
		CALC.auto_mixer=TRUE;
		CALC.imix=VM_4;
		CALC.amix=0.4;
		CALC.bmix=1.0;
		CALC.amin=0.1;
		CALC.amix_mag=1.6;
		CALC.bmix_mag=1.0;
		CALC.maxmix=-45;
		CALC.wc=1000.;
		CALC.inimix=VIM_1;
		CALC.mixpre=VMP_1;
/*9-dipole*/
		CALC.idipol=VID_0;
		CALC.ldipol=FALSE;
		CALC.lmono=FALSE;
		CALC.epsilon=1.0;/*assuming dielectric constant for vacuum*/
		CALC.dipol=NULL;
		CALC.efield=0.0;
/*10-dos*/
		CALC.lorbit=0;
		CALC.nedos=301;
		CALC.emin=0.0;/*ie not set*/
		CALC.emax=0.0;/*ie not set*/
		CALC.efermi=0.0;/*undocument*/
		CALC.rwigs=NULL;
/*11-linear response*/
		CALC.loptics=FALSE;
		CALC.lepsilon=FALSE;
		CALC.lrpa=FALSE;
		CALC.lnabla=FALSE;
		CALC.cshift=0.1;
		CALC.lcalceps=FALSE;/*use for hybride-DFT TODO (adv.)*/
/*ionic*/
/*0-grid*/
		CALC.auto_grid=TRUE;/*no default values for grids integers*/
/*1-general*/
		CALC.nsw=0;
		CALC.ibrion=-1;
		CALC.isif=2;
		CALC.pstress=0.0;
		CALC.ediffg=CALC.ediff*10.0;
		CALC.potim=0.5;
		CALC.nfree=2;/*actually there is no default*/
/*2-md*/
		CALC.smass=-3.0;
		CALC.npaco=256;
		CALC.apaco=16.;/*an integer in Ang?*/
		CALC.tebeg=0.0;
		CALC.teend=CALC.tebeg;
		CALC.nblock=1;
		CALC.kblock=CALC.nsw;
/*3-symmetry*/
		CALC.isym=2;/*let's take the PAW default*/
		CALC.sym_prec=1e-5;
/* PERFS */
		CALC.lplane=TRUE;
		CALC.lscalu=TRUE;
		CALC.lscalapack=TRUE;
		/*output*/
		CALC.nwrite=2;
		CALC.lwave=TRUE;
		CALC.lcharg=TRUE;
		CALC.lvtot=FALSE;
		CALC.lvhar=FALSE;
		CALC.lelf=FALSE;
	}
	return 0;
#undef CALC
}
/****************************/
/* free vasp_calc structure */
/****************************/
void vasprun_free(vasp_calc_struct *calc){
#define CALC (*calc)
	if(calc==NULL) return;
	if(CALC.name!=NULL) g_free(CALC.name);
	if(CALC.fermwe!=NULL) g_free(CALC.fermwe);
	if(CALC.fermdo!=NULL) g_free(CALC.fermdo);
	if(CALC.ropt!=NULL) g_free(CALC.ropt);
	if(CALC.magmom!=NULL) g_free(CALC.magmom);
	if(CALC.saxis!=NULL) g_free(CALC.saxis);
	if(CALC.cmbj!=NULL) g_free(CALC.cmbj);
	if(CALC.ldaul!=NULL) g_free(CALC.ldaul);
	if(CALC.ldauu!=NULL) g_free(CALC.ldauu);
	if(CALC.ldauj!=NULL) g_free(CALC.ldauj);
	if(CALC.dipol!=NULL) g_free(CALC.dipol);
	if(CALC.species_symbols!=NULL) g_free(CALC.species_symbols);
	if(CALC.species_numbers!=NULL) g_free(CALC.species_numbers);
	if(CALC.potcar_folder!=NULL) g_free(CALC.potcar_folder);
	if(CALC.potcar_file!=NULL) g_free(CALC.potcar_file);
	if(CALC.potcar_species!=NULL) g_free(CALC.potcar_species);
	if(CALC.potcar_species_flavor!=NULL) g_free(CALC.potcar_species_flavor);
	if(CALC.rwigs!=NULL) g_free(CALC.rwigs);
	if(CALC.job_vasp_exe!=NULL) g_free(CALC.job_vasp_exe);
	if(CALC.job_mpirun!=NULL) g_free(CALC.job_mpirun);
	if(CALC.job_path!=NULL) g_free(CALC.job_path);
//result_model should not be freed ;)
#undef CALC
}
/***********************************/
/* convert vasp structure to INCAR */
/***********************************/
void vasp_calc_to_incar(FILE *output,vasp_calc_struct calc){
/*skipping default parts*/
	fprintf(output,"!INCAR GENERATED BY GDIS %4.2f.%d (C) %d\n",VERSION,PATCH,YEAR);
	if(calc.name!=NULL) fprintf(output,"SYSTEM=%s\n",calc.name);
	fprintf(output,"PREC=");
	switch (calc.prec){
		case VP_SINGLE:
			fprintf(output,"SINGLE\n");break;
		case VP_ACCURATE:
			fprintf(output,"ACCURATE\n");break;
		case VP_HIGH:
			fprintf(output,"HIGH\n");break;
		case VP_MED:
			fprintf(output,"MEDIUM\n");break;
		case VP_LOW:
			fprintf(output,"LOW\n");break;
		case VP_NORM:
		default:
			fprintf(output,"NORMAL\n");
	}
	if(calc.use_prec) fprintf(output,"!USING PREC SETTINGS FOR ENCUT and ENAUG\n");
	else{
		if(calc.encut!=0.0) fprintf(output,"ENCUT=%lf\n",calc.encut);
		if(calc.enaug!=0.0) fprintf(output,"ENAUG=%lf\n",calc.enaug);
	}
	if(calc.ediff!=1E-4) fprintf(output,"EDIFF=%lE\n",calc.ediff);
	if(calc.algo!=VA_IALGO){
		fprintf(output,"ALGO=");
		switch (calc.algo){
		case VA_NORM:
			fprintf(output,"NORMAL\n");break;
		case VA_VERYFAST:
			fprintf(output,"VERYFAST\n");break;
		case VA_FAST:
			fprintf(output,"FAST\n");break;
		case VA_CONJ:
			fprintf(output,"CONJ\n");break;
		case VA_ALL:
			fprintf(output,"ALL\n");break;
		case VA_DAMPED:
			fprintf(output,"DAMPED\n");break;
		case VA_SUBROT:
			fprintf(output,"SUBROT\n");break;
		case VA_EIGEN:
			fprintf(output,"EIGENVAL\n");break;
		case VA_NONE:
			fprintf(output,"NONE\n");break;
		case VA_NOTHING:
			fprintf(output,"NOTHING\n");break;
		case VA_EXACT:
			fprintf(output,"EXACT\n");break;
		case VA_DIAG:
			fprintf(output,"DIAG\n");break;
		case VA_IALGO:
		default:
			/*should never happen*/
			fprintf(output,"NORMAL\n");
		}
	} else fprintf(output,"IALGO=%i\n",(gint)calc.ialgo);
	if(!(calc.ldiag)) fprintf(output,"LDIAG=.FALSE.\n");
	if(!(calc.auto_elec)){
		if(calc.nbands!=0) fprintf(output,"NBANDS=%i\n",calc.nbands);
		if(calc.nelect!=0.0) fprintf(output,"NELECT=%lf\n",calc.nelect);
		if(calc.iwavpr!=10) fprintf(output,"IWAVPR=%i\n",calc.iwavpr);
		if(calc.nsim!=4) fprintf(output,"NSIM=%i\n",calc.nsim);
		if(calc.vtime!=0.4) fprintf(output,"VTIME=%lf\n",calc.vtime);
	}
	if(calc.ismear!=1) fprintf(output,"ISMEAR=%i\n",calc.ismear);
	if(calc.sigma!=0.2) fprintf(output,"SIGMA=%lf\n",calc.sigma);
	if(calc.fermwe!=NULL) fprintf(output,"FERMWE=%s\n",calc.fermwe);
	if(calc.fermdo!=NULL) fprintf(output,"FERMDO=%s\n",calc.fermdo);
	if(calc.kspacing!=0.5) fprintf(output,"KSPACING=%lf\n",calc.kspacing);
	if(!calc.kgamma) fprintf(output,"KGAMMA=.FALSE.\n");
	/*always write LREAL tag*/
	fprintf(output,"LREAL=");
	switch (calc.lreal){
	case VLR_AUTO:
		fprintf(output,"AUTO\n");break;
	case VLR_ON:
		fprintf(output,"ON\n");break;
	case VLR_TRUE:
		fprintf(output,".TRUE.\n");break;
	case VLR_FALSE:
	default:
		fprintf(output,".FALSE.\n");
	}
	if(calc.ropt!=NULL) fprintf(output,"ROPT=%s\n",calc.ropt);
	if(calc.addgrid) fprintf(output,"ADDGRID=.TRUE.\n");
	if(calc.lmaxmix!=2) fprintf(output,"LMAXMIX=%i\n",calc.lmaxmix);
	if(calc.lmaxpaw!=-100) fprintf(output,"LMAXPAW=%i\n",calc.lmaxpaw);
        if(calc.istart!=-1) fprintf(output,"ISTART=%i\n",calc.istart);
        if(calc.icharg!=-1) fprintf(output,"ICHARG=%i\n",calc.icharg);
        if(!calc.iniwav) fprintf(output,"INIWAV=0\n");
	if(calc.ispin) fprintf(output,"ISPIN=2\n");
	if(calc.non_collinear) fprintf(output,"LNONCOLLINEAR=.TRUE.\n");
	if(calc.magmom!=NULL) fprintf(output,"MAGMOM=%s\n",calc.magmom);
	if(calc.nupdown!=-1) fprintf(output,"NUPDOWN=%lf\n",calc.nupdown);
	if(calc.lsorbit) fprintf(output,"LSORBIT=.TRUE.\n");
	if(calc.saxis!=NULL) fprintf(output,"MAGMOM=%s\n",calc.saxis);
	if(!calc.gga_compat) fprintf(output,"GGA_COMPAT=.FALSE.\n");
	/*always write GGA tag*/
	fprintf(output,"GGA=");
	switch (calc.gga){
	case VG_91:
		fprintf(output,"91\n");break;
	case VG_RP:
		fprintf(output,"RP\n");break;
	case VG_AM:
		fprintf(output,"AM\n");break;
	case VG_PS:
		fprintf(output,"PS\n");break;
	case VG_PE:
	default:
		fprintf(output,"PE\n");
	}
	if(calc.voskown) fprintf(output,"VOSKOWN=1\n");
	if(calc.lasph) fprintf(output,"LASPH=.TRUE.\n");
	if(calc.lmaxtau!=0) fprintf(output,"LMAXTAU=%i\n",calc.lmaxtau);
	if(calc.lmixtau) fprintf(output,"LMIXTAU=.TRUE.\n");
	if(calc.lmetagga){
		fprintf(output,"METAGGA=");
		switch (calc.mgga){
		case VMG_TPSS:
			fprintf(output,"TPSS\n");break;
		case VMG_RTPSS:
			fprintf(output,"RTPSS\n");break;
		case VMG_M06L:
			fprintf(output,"M06L\n");break;
		case VMG_MBJ:
		default:
			fprintf(output,"MBJ\n");
		}
		if(calc.mgga==VMG_MBJ){/*print parameters as well*/
			if(calc.cmbj!=NULL) fprintf(output,"CMBJ=%s\n",calc.cmbj);
			if(calc.cmbja!=-0.012) fprintf(output,"CMBJA=%lf\n",calc.cmbja);
			if(calc.cmbjb!=1.023) fprintf(output,"CMBJB=%lf\n",calc.cmbjb);
		}
	}
	if(calc.ldau){
		fprintf(output,"LDAU=.TRUE.\n");
		fprintf(output,"LDAUTYPE=%i\n",(gint)calc.ldau_type);/*always write*/
		fprintf(output,"LDAUPRINT=%i\n",(gint)calc.ldau_output);/*always write*/
		if(calc.ldaul!=NULL) fprintf(output,"LDAUL=%s\n",calc.ldaul);
		if(calc.ldauu!=NULL) fprintf(output,"LDAUU=%s\n",calc.ldauu);
		if(calc.ldauj!=NULL) fprintf(output,"LDAUJ=%s\n",calc.ldauj);
	}
	if(calc.nelm!=60) fprintf(output,"NELM=%i\n",calc.nelm);
	if(calc.nelmdl!=-12) fprintf(output,"NELMDL=%i\n",calc.nelmdl);
	if(calc.nelmin!=2) fprintf(output,"NELMIN=%i\n",calc.nelmin);
	if(!calc.auto_mixer){
		fprintf(output,"IMIX=%i\n",(gint)calc.imix);/*always write*/
		if(calc.amix!=0.4) fprintf(output,"AMIX=%lf\n",calc.amix);
		if(calc.bmix!=1.0) fprintf(output,"BMIX=%lf\n",calc.bmix);
		if(calc.amin!=0.1) fprintf(output,"AMIN=%lf\n",calc.amin);
		if(calc.amix_mag!=1.6) fprintf(output,"AMIX_MAG=%lf\n",calc.amix_mag);
		if(calc.bmix_mag!=1.0) fprintf(output,"BMIX_MAG=%lf\n",calc.bmix_mag);
		if(calc.maxmix!=-45) fprintf(output,"MAXMIX=%i\n",calc.maxmix);
		if(calc.wc!=1000.0) fprintf(output,"WC=%lf\n",calc.wc);
		if(calc.inimix!=VIM_1) fprintf(output,"INIMIX=%i\n",(gint)calc.inimix);
		if(calc.mixpre!=VMP_1) fprintf(output,"MIXPRE=%i\n",(gint)calc.mixpre);
	}
	if(calc.idipol!=VID_0){
		if(calc.ldipol) fprintf(output,"LDIPOL=.TRUE.\n");
		if(calc.lmono) fprintf(output,"LMONO=.TRUE.\n");
		if(calc.epsilon!=1.0) fprintf(output,"EPSILON=%lf\n",calc.epsilon);
		if(calc.dipol!=NULL) fprintf(output,"DIPOL=%s\n",calc.dipol);
		if((calc.efield!=0.0)&&(calc.idipol!=VID_4)) fprintf(output,"EFIELD=%lf\n",calc.efield);
	}
	if(calc.lorbit!=0) fprintf(output,"LORBIT=%i\n",calc.lorbit);
	if(calc.nedos!=301) fprintf(output,"NEDOS=%i\n",calc.nedos);
	if(calc.emin!=0.) fprintf(output,"EMIN=%lf\n",calc.emin);
	if(calc.emax!=0.) fprintf(output,"EMAX=%lf\n",calc.emax);
	if(calc.efermi!=0.) fprintf(output,"EFERMI=%lf\n",calc.efermi);
	if(calc.rwigs!=NULL) fprintf(output,"RWIGS=%s\n",calc.rwigs);
	if(calc.loptics) fprintf(output,"LOPTICS=.TRUE.\n");
	if(calc.lepsilon) fprintf(output,"LEPSILON=.TRUE.\n");
	if(calc.lrpa) fprintf(output,"LRPA=.TRUE.\n");
	if(calc.lnabla) fprintf(output,"LNABLA=.TRUE.\n");
	if(calc.lcalceps) fprintf(output,"LCALCEPS=.TRUE.\n");
	if(calc.cshift!=0.1) fprintf(output,"CSHIFT=%lf\n",calc.cshift);
	if(!calc.auto_grid){
		if(calc.ngx>0) fprintf(output,"NGX=%i\n",(gint)calc.ngx);
		if(calc.ngy>0) fprintf(output,"NGY=%i\n",(gint)calc.ngy);
		if(calc.ngz>0) fprintf(output,"NGZ=%i\n",(gint)calc.ngz);
		if(calc.ngxf>0) fprintf(output,"NGXF=%i\n",(gint)calc.ngxf);
		if(calc.ngyf>0) fprintf(output,"NGYF=%i\n",(gint)calc.ngyf);
		if(calc.ngzf>0) fprintf(output,"NGZF=%i\n",(gint)calc.ngzf);
	}
	if(calc.nsw!=0) fprintf(output,"NSW=%i\n",calc.nsw);
	if(calc.ibrion!=-1) fprintf(output,"IBRION=%i\n",calc.ibrion);
	if(calc.isif!=2) fprintf(output,"ISIF=%i\n",calc.isif);
	if(calc.pstress!=0.0) fprintf(output,"PSTRESS=%lf\n",calc.pstress);
	if(calc.ediffg!=(calc.ediff*10.0)) fprintf(output,"EDIFFG=%lE\n",calc.ediffg);
	if(calc.potim!=0.5) fprintf(output,"POTIM=%lf\n",calc.potim);
	if(calc.nfree!=2) fprintf(output,"NFREE=%i\n",calc.nfree);
	if(calc.ibrion==0){
		if(calc.smass!=-3.0) fprintf(output,"SMASS=%lf\n",calc.smass);
		if(calc.npaco!=256) fprintf(output,"NPACO=%i\n",calc.npaco);
		if(calc.apaco!=16.) fprintf(output,"APACO=%lf\n",calc.apaco);
		if(calc.tebeg!=0.0) fprintf(output,"TEBEG=%lf\n",calc.tebeg);
		if(calc.teend!=calc.tebeg) fprintf(output,"TEEND=%lf\n",calc.teend);
		if(calc.nblock!=1) fprintf(output,"NBLOCK=%i\n",calc.nblock);
		if(calc.kblock!=calc.nsw) fprintf(output,"KBLOCK=%i\n",calc.kblock);
	}
	if(calc.isym!=2) fprintf(output,"ISYM=%i\n",calc.isym);
	if((calc.isym>0)&&(calc.sym_prec!=1e-5)) fprintf(output,"SYM_PREC=%lE\n",calc.sym_prec);
	if(calc.ncore>1.) fprintf(output,"NCORE=%i\n",(gint)calc.ncore);
	if(calc.kpar>1.) fprintf(output,"KPAR=%i\n",(gint)calc.kpar);
	if(!calc.lplane) fprintf(output,"LPLANE=.FALSE.\n");
	if(!calc.lscalu) fprintf(output,"LSCALU=.FALSE.\n");
	if(!calc.lscalapack) fprintf(output,"LSCALAPACK=.FALSE.\n");
	if(!calc.lwave) fprintf(output,"LWAVE=.FALSE.\n");
	if(!calc.lcharg) fprintf(output,"LCHARG=.FALSE.\n");
	if(calc.lvtot) fprintf(output,"LVTOT=.TRUE.\n");
	if(calc.lvhar) fprintf(output,"LVHAR=.TRUE.\n");
	if(calc.lelf) fprintf(output,"LELF=.TRUE.\n");

	/*TODO (adv.)*/
	
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
        model->num_species=idx;
	model->expected_cores=natom;
        model->expected_shells=0;
	model->num_atoms=natom;
        /* exit at </atominfo> */
	if(fetch_in_file(vf,"</atominfo>")==0) return -1;
	return 0;
}

int vasp_xml_read_energy(FILE *vf, struct model_pak *model){
	gchar *line;
	gdouble energy=0.;
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
		sscanf(line," <i name=\"e_wo_entrp\"> %lf </i>",&energy);
		sprintf(line,"%lf eV",energy);
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
	if (idx != model->num_atoms) /* This should not happen */
		fprintf(stderr,"WARNING: Expecting %i atoms but got %i!\n",model->num_atoms,idx);
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
	gdouble init;
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
	/*The first line is special: when integration is not 0, it contain a huge value:
 * 	  while this value makes intergration of DOS (ie. nb of electron) correct, it is
 * 	  not a correct representation of DOS at that point... (FIXED by setting to 0)*/
	sscanf(line," <r> %lf %lf %lf %*s",&(model->dos_eval[idx]),&(model->dos_spin_up[idx]),&init);
	if(init > 0.0) model->dos_spin_up[0]=0.0;
	line = file_read_line(vf);
	idx++;
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
/*same as above*/
	sscanf(line," <r> %lf %lf %lf %*s",&(model->dos_eval[idx]),&(model->dos_spin_down[idx]),&init);
	if(init > 0.0) model->dos_spin_up[0]=0.0;
	line = file_read_line(vf);
	idx++;
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
	if(vasp_xml_read_dos(vf,model)){;/* Read DOS */
		model->ndos=0;/*didn't work*/
	}
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
/*********************/
/* helpers functions */
/*********************/
gint vasp_load_poscar5(FILE *vf,struct model_pak *model){
	gint idx,ix;
	gchar *line;
	gchar *label;
	gchar *spec;
	gchar *name;
	gdouble a0;
	/*atom determination*/
	gchar *ptr;
	gchar *ptr2;
	gchar *ptr3;
	gchar sym[3];
	struct core_pak *core;
	gboolean is_direct;
	/*read a vasp5 formated POSCAR*/
	if(vf==NULL) return 1;
	if(model==NULL) return 1;
	if(feof(vf)) return 1;
	/*we are good to go*/
        line = file_read_line(vf);/*title*/
	if(line==NULL) return 1;
/*add local structure title to properties?*/
	line = file_read_line(vf);/*lattice parameter*/
	if(line==NULL) return 1;
	sscanf(line,"%lf%*s",&(a0));
	idx=0;
	line = file_read_line(vf);/*basis*/
        while (idx<3) {
                sscanf(line," %lf %lf %lf%*s",&model->latmat[idx],&model->latmat[idx+3],&model->latmat[idx+6]);
		/*multiply everything by lattice parameter so a0 -> 1*/
		model->latmat[idx]*=a0;
		model->latmat[idx+3]*=a0;
		model->latmat[idx+6]*=a0;
                idx+=1;
                g_free(line);
                line = file_read_line(vf);
		if(line == NULL) return -1;/*incomplete basis*/
        }
	label=g_strdup(line);/*atomic symbols*/
	g_free(line);
	line = file_read_line(vf);/*number of each species*/
	idx=0;ptr2=&(line[0]);ptr3=label;
	sym[2]='\0';/*always*/
	model->num_atoms=0;
/*FIX _BUG_ core list grow!*/
core_delete_all(model);
	spec=NULL;
	name=g_strdup_printf("%c",'\0');
        do{
		ptr=ptr2;
		ix=g_ascii_strtod(ptr,&ptr2);
		if(ptr2==ptr) break;
		while(*ptr3==' ') ptr3++;
		sym[0]=*ptr3;
		ptr3++;
		if((*ptr3==' ')||(*ptr3=='\0')||(*ptr3=='\n')) sym[1]='\0';
		else sym[1]=*ptr3;
		ptr3++;
		for(idx=0;idx<ix;idx++){
                        core=new_core(sym,model);
                        core->charge=0.;/*no such information on POSCAR*/
                        model->cores=g_slist_append(model->cores,core);
		}
		g_free(spec);
		spec=g_strdup_printf("%s%s(%i)",name,sym,ix);
		g_free(name);
		name=g_strdup(spec);
                model->num_atoms+=ix;
        }while(1);
property_add_ranked(7, "Formula", name, model);
	g_free(line);
	/*Always true in VASP*/
        model->fractional=TRUE;
        model->coord_units=ANGSTROM;
        model->construct_pbc = TRUE;
        model->periodic = 3;
	line = file_read_line(vf);/*direct/cartesian switch*/
	if((line[0]=='d')||(line[0]=='D')) is_direct=TRUE;
	else if((line[0]=='c')||(line[0]=='C')||(line[0]=='k')||(line[0]=='K')) is_direct=FALSE;
	else {/*no info: bailout*/
		core_delete_all(model);
		return -2;
	}
        line = file_read_line(vf);
	if((line[0]=='s')||(line[0]=='S')) {
		/*Selective switch: skip*/
		g_free(line);
		line = file_read_line(vf);
	}
	/*now start registering atoms*/
	for(idx=0;idx<(model->num_atoms-1);idx++){
		core=g_slist_nth_data(model->cores,idx);
                sscanf(line," %lf %lf %lf%*s",&core->x[0],&core->x[1],&core->x[2]);
		if(!is_direct){
			core->x[0]/=(model->latmat[0]+model->latmat[1]+model->latmat[2]);
			core->x[1]/=(model->latmat[3]+model->latmat[4]+model->latmat[5]);
			core->x[2]/=(model->latmat[6]+model->latmat[7]+model->latmat[8]);
		}
                g_free(line);
                line = file_read_line(vf);
		if(line==NULL){
			/*problem reading atoms: bailout*/
			core_delete_all(model);
			return -2;
		}
        }
	/*do the last one outside of loop to avoid file_read_linegoing too far */
	core=g_slist_nth_data(model->cores,model->num_atoms-1);
	sscanf(line," %lf %lf %lf%*s",&core->x[0],&core->x[1],&core->x[2]);
	if(!is_direct){
		core->x[0]/=(model->latmat[0]+model->latmat[1]+model->latmat[2]);
		core->x[1]/=(model->latmat[3]+model->latmat[4]+model->latmat[5]);
		core->x[2]/=(model->latmat[6]+model->latmat[7]+model->latmat[8]);
	}
	g_free(line);
	/*should be all done*/
	return 0;
}


