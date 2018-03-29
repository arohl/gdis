/*
Copyright (C) 2004 by Andrew Lloyd Rohl

andrew@ivec.org

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

#include <math.h>
#include <stdio.h>
#include <string.h>

#include "gdis.h"
#include "coords.h"
#include "file.h"
#include "parse.h"
#include "model.h"
#include "interface.h"
#include "matrix.h"

/* main structures */
extern struct sysenv_pak sysenv;
extern struct elem_pak elements[];

/* VZ routines to read .cell files */
/*************************************************/
/* read a .cell block into the model structure     */
/*************************************************/
gint read_cell_block(FILE *fp, struct model_pak *data)
{
    gchar line[LINELEN], *atom_name, *element, **buf;
    struct core_pak *core;
    GSList *clist;
    gint num_tokens;
    
    clist = data->cores;
    /* loop while there's data */
    while (!fgetline(fp, line))
    {
        
        /* 3-D structure? */
        if (g_strrstr(line, "BLOCK LATTICE_ABC") != 0)
        {
            printf("Here\n");
            if (fgetline(fp, line)) {
                gui_text_show(ERROR, "unexpected end of file skipping to gradients\n");
                return(2);
            }
            data->periodic = 3;
            buf=tokenize(line,&num_tokens);
            printf("entries %d\n",num_tokens);
            data->pbc[0] = str_to_float(*(buf+0));
            data->pbc[1] = str_to_float(*(buf+1));
            data->pbc[2] = str_to_float(*(buf+2));
            if (fgetline(fp, line)) {
                gui_text_show(ERROR, "unexpected end of file skipping to gradients\n");
                return(2);
            }
            g_strfreev(buf);
            buf=tokenize(line,&num_tokens);
            printf("entries %d\n",num_tokens);
            data->pbc[3] = str_to_float(*(buf+0))* PI/180.0;
            data->pbc[4] = str_to_float(*(buf+1)) * PI/180.0;
            data->pbc[5] = str_to_float(*(buf+2)) * PI/180.0;
            printf("al be ga %f %f %f\n",data->pbc[3],data->pbc[4],data->pbc[5]);
            if (fgetline(fp, line)) {
                gui_text_show(ERROR, "unexpected end of file skipping to gradients\n");
                return(2);
            } /*omit next line*/
        }
        else if (g_strrstr(line, "BLOCK POSITIONS_") != 0) {
            
            if (g_strrstr(line, "BLOCK POSITIONS_ABS") != 0)
                data->fractional = FALSE;
            else
                data->fractional = TRUE;
            
            if (fgetline(fp, line)) {
                gui_text_show(ERROR, "unexpected end of file skipping to gradients\n");
                return(2);
            } /* go to next line */
            
            while (g_strrstr(line, "ENDBLOCK POSITIONS_") == 0) {
                buf = tokenize(line,&num_tokens);
                element = g_strdup(*(buf+0));
                atom_name = g_strdup(element);
                if (clist)
                {
                    core = clist->data;
                    clist = g_slist_next(clist);
                }
                else
                {
                    core = core_new(element, atom_name, data);
                    data->cores = g_slist_append(data->cores, core);
                }
                g_free(atom_name);
                
                core->x[0] = str_to_float(*(buf+1));
                core->x[1] = str_to_float(*(buf+2));
                core->x[2] = str_to_float(*(buf+3));
                g_strfreev(buf);
                if (fgetline(fp, line)) {
                    gui_text_show(ERROR, "unexpected end of file skipping to gradients\n");
                    return(2);
                }
            }
        }
    }
    return(-1);
}

/*********************/
/* PDB frame read    */
/*********************/
gint read_cell_frame(FILE *fp, struct model_pak *data)
{
    /* frame overwrite */
    read_cell_block(fp, data);
    return(0);
}

/******************/
/* PDB reading    */
/******************/
#define DEBUG_READ_PDB 0
gint read_cell(gchar *filename, struct model_pak *data)
{
    gint frame=0;
    FILE *fp;
    
    /* checks */
    g_return_val_if_fail(data != NULL, 1);
    g_return_val_if_fail(filename != NULL, 2);
    
    /* initialisers */
    data->periodic = 0;
    
    fp = fopen(filename, "rt");
    if (!fp)
        return(3);
    add_frame_offset(fp, data);
    
    while (read_cell_block(fp, data) == 0)
    {
        add_frame_offset(fp, data);
        /* increment counter */
        frame++;
    }
    
    /* model setup */
    strcpy(data->filename, filename);
    g_free(data->basename);
    data->basename = parse_strip(filename);
    
    data->num_frames = data->cur_frame = frame;
    data->cur_frame--;
    
    model_prep(data);
    return(0);
}

/*******************************************/
/* read single CASTEP output configuration */
/*******************************************/
#define DEBUG_READ_CASTEP_OUT 0
gint read_castep_out_block(FILE *fp, struct model_pak *model)
{
gint i;
gint num_tokens;
gint no_grad = 0;
gint last_frame = FALSE;
gdouble grad, max_grad = 0.0, rms_grad = 0.0;
gdouble pressure;
gchar **buff, line[LINELEN], *ext, *text;
GString *title, *grad_string, *pressure_string;
GSList *clist;
struct core_pak *core;

clist = model->cores;

/* read to end of iteration */
while (TRUE)
  {
  if (fgetline(fp, line))
    {
    gui_text_show(ERROR, "unexpected end of file reading to end of iteration\n");
    return(2);
    }
    

  if (g_ascii_strncasecmp(line, "==========", 10) == 0)
      break;
      
  if (g_ascii_strncasecmp(line, "Writing model to", 16) == 0)
      break;
  
/* read cell vectors? */
  if (g_strrstr(line, "Real Lattice") != NULL)
    {
    /* read in cell vectors */
    /* NB: gdis wants transposed matrix */
    for (i=0; i<3; i++)
      {
      if (fgetline(fp, line))
        {
        gui_text_show(ERROR, "unexpected end of file reading cell vectors\n");
        return(2);
        }
      buff = tokenize(line, &num_tokens);
      model->latmat[0+i] = str_to_float(*(buff+0));
      model->latmat[3+i] = str_to_float(*(buff+1));
      model->latmat[6+i] = str_to_float(*(buff+2));
      g_strfreev(buff);
      }
    }


gdouble* shifts[2000];
int shift_counter=0;
gdouble* aniso[2000];
gdouble* asym[2000];
gdouble* cq[2000];
gdouble* efgasym[2000];

/* read coordinates */
  if (g_strrstr(line, "x  Element    Atom        Fractional coordinates of atoms  x") != NULL)
    {

int skip=0;
  for (skip=0; skip<3; skip++)
    if (fgetline(fp, line))
      {
      gui_text_show(ERROR, "unexpected end of file skipping to coordinates\n");
      return(2);
      }


    buff = tokenize(line, &num_tokens);
    while (num_tokens == 7)
      {
      if (clist)
        {
        core = (struct core_pak *) clist->data;
        clist = g_slist_next(clist);
        }
      else
        {
	gchar* atom_name=g_strdup(*(buff+1));
/* Next line is corrected by Anne-Christine Uldry. Works if number of spicies >100 */
	g_strlcat(atom_name,*(buff+2),5);
/*	g_strlcat(atom_name,*(buff+2),4);*/
        core = core_new(*(buff+1),atom_name, model);
        model->cores = g_slist_append(model->cores, core);
        }
      core->x[0] = str_to_float(*(buff+3));
      core->x[1] = str_to_float(*(buff+4));
      core->x[2] = str_to_float(*(buff+5));

      shifts[shift_counter]=&core->atom_nmr_shift;
      aniso[shift_counter]=&core->atom_nmr_aniso;
      asym[shift_counter]=&core->atom_nmr_asym;
      cq[shift_counter]=&core->atom_nmr_cq;
      efgasym[shift_counter]=&core->atom_nmr_efgasym;
      shift_counter++;

      #if DEBUG_READ_CASTEP_OUT
      printf("new coords %f %f %f\n", core->x[0],  core->x[1], core->x[2]);
      #endif

    /* get next line */
      g_strfreev(buff);
      if (fgetline(fp, line))
        {
        gui_text_show(ERROR, "unexpected end of file reading coordinates\n");
        return(2);
        }
      buff = tokenize(line, &num_tokens);
      }
    g_strfreev(buff);
    clist = model->cores;
    }


/*Iso chemical shifts */
  if (g_strrstr(line, "|  Species   Ion    Iso(ppm)") != NULL)
  {
    if (fgetline(fp, line))
      {
      gui_text_show(ERROR, "unexpected end of file skipping to coordinates\n");
      return(2);
      }
    buff = tokenize(line, &num_tokens);
    int i=0;
    while ((num_tokens == 9) || (num_tokens == 7))
    {
/*	if (i<=shift_counter)*/
	*shifts[i]=str_to_float(*(buff+3));
	*aniso[i]=str_to_float(*(buff+4));
	*asym[i]=str_to_float(*(buff+5));
	if (num_tokens == 9) { //read EFG parameters if present
		*cq[i]=str_to_float(*(buff+6));
		*efgasym[i]=str_to_float(*(buff+7));
	}
	i++;
    /* get next line */
      g_strfreev(buff);
      if (fgetline(fp, line))
        {
        gui_text_show(ERROR, "unexpected end of file reading coordinates\n");
        return(2);
        }
      buff = tokenize(line, &num_tokens);
      }
    g_strfreev(buff);
    }
/* EFG only computation*/
  if (g_strrstr(line, " |  Species   Ion             Cq(MHz)") != NULL)
  {
    if (fgetline(fp, line))
      {
      gui_text_show(ERROR, "unexpected end of file skipping to coordinates\n");
      return(2);
      }
    buff = tokenize(line, &num_tokens);
    int i=0;
    while (num_tokens == 6)
    {
	*cq[i]=str_to_float(*(buff+3));
	*efgasym[i]=str_to_float(*(buff+4));
	i++;
    /* get next line */
      g_strfreev(buff);
      if (fgetline(fp, line))
        {
        gui_text_show(ERROR, "unexpected end of file reading coordinates\n");
        return(2);
        }
      buff = tokenize(line, &num_tokens);
      }
    g_strfreev(buff);
    }



/* SCF or Final Energy */
  if (g_strrstr(line, "Number of kpoints used") != NULL)
    {
    buff = g_strsplit(line, "=", 2);
    g_strchug(g_strchomp(*(buff+1)));
    property_add_ranked(10, "K-points", *(buff+1), model);
    g_strfreev(buff);
    }

  if (g_strrstr(line, "Files used for pseudopotentials:") != NULL)
    {
    if (fgetline(fp, line))
      {
      gui_text_show(ERROR, "unexpected end of file reading pseudotential file names\n");
      return(2);
      }
      buff = tokenize(line, &num_tokens);
      ext = find_char(*(buff+1), '.', LAST);
      ext++;
      if (g_strrstr(ext, "usp") != NULL)
        property_add_ranked(10, "Pseudopotentials", "Ultrasoft", model);
      else if (g_strrstr(ext, "recpot") != NULL)
        property_add_ranked(10, "Pseudopotentials", "Norm conserving", model);
      else
        property_add_ranked(10, "Pseudopotentials", "Unknown", model);
    }
  if (g_ascii_strncasecmp(line, "Final energy", 12) == 0)
    {
    buff = g_strsplit(line, "=", 2);
    text = format_value_and_units(*(buff+1), 5);
    property_add_ranked(3, "Energy", text, model);
    g_free(text);
    model->castep.energy = str_to_float(*(buff+1));
    model->castep.have_energy = TRUE;
    g_strfreev(buff);
    }
  if (g_strrstr(line, "Final Enthalpy") != NULL)
    {
    last_frame = TRUE;
    model->castep.min_ok = TRUE;
    buff = g_strsplit(line, "=", 2);
    text = format_value_and_units(*(buff+1), 5);
    property_add_ranked(3, "Energy", text, model);
    g_free(text);
    model->castep.energy = str_to_float(*(buff+1));
    model->castep.have_energy = TRUE;
    g_strfreev(buff);
    }

  /* gradients */
  if (g_strrstr(line, "* -") != NULL)
    {
    for (i=0; i<2; i++)
      {
      if (fgetline(fp, line))
        {
        gui_text_show(ERROR, "unexpected end of file skipping to gradients\n");
        return(2);
        }
      }
    /* read gradients */
    if (fgetline(fp, line))
      {
      gui_text_show(ERROR, "unexpected end of file reading gradients\n");
      return(2);
      }
    buff = tokenize(line, &num_tokens);
    while (num_tokens == 7)
      {
      for (i=0; i<3; i++)
        {
        grad = str_to_float(*(buff+3+i));
        no_grad++;
        rms_grad = grad*grad;
        max_grad = MAX(fabs(grad), max_grad);
      }
      /* get next line */
      g_strfreev(buff);
      if (fgetline(fp, line))
        {
        gui_text_show(ERROR, "unexpected end of file reading gradients\n");
        return(2);
        }
      buff = tokenize(line, &num_tokens);
      }
    g_strfreev(buff);
    model->castep.rms_grad = sqrt(rms_grad/no_grad);
    model->castep.have_rms_grad = TRUE;
    model->castep.max_grad = max_grad;
    model->castep.have_max_grad = TRUE;
    grad_string = g_string_new("");
    g_string_append_printf(grad_string, "%.5f %s", model->castep.rms_grad, property_lookup("Gradient Units", model));
    property_add_ranked(4, "RMS Gradient", grad_string->str, model);
    g_string_free(grad_string, TRUE); 

    #if DEBUG_READ_CASTEP_OUT
    printf("RMS gradient: %f Max gradient: %f\n", model->castep.rms_grad, model->castep.max_grad);
    #endif
    }
      
  if (g_strrstr(line, "Pressure") != NULL)
    {
    /* read pressure */
    buff = tokenize(line, &num_tokens);
    pressure = str_to_float(*(buff+2));
    pressure_string = g_string_new("");
    g_string_append_printf(pressure_string, "%.4f %s", pressure, property_lookup("Pressure Units", model));
    property_add_ranked(5, "Pressure", pressure_string->str, model);
    g_string_free(pressure_string, TRUE);
    g_strfreev(buff);
    }
  }
  
g_free(model->title);
title = g_string_new("");
g_string_append_printf(title, "E");
g_string_append_printf(title, "(DFT)");
/* TODO read energy units from CASTEP file */
g_string_append_printf(title, " = %.5f eV", model->castep.energy);
g_string_append_printf(title, ", grad = %.5f", model->castep.rms_grad);
model->title = g_strdup(title->str);
g_string_free(title, TRUE); 

return(0);
}

/*******************************/
/* CASTEP output frame reading */
/*******************************/
gint read_castep_out_frame(FILE *fp, struct model_pak *model)
{
    return(read_castep_out_block(fp, model));
}



/********************************/
/* Read in a CASTEP output file */
/********************************/
gint read_castep_out(gchar *filename, struct model_pak *model)
{
gint frame;
gchar **buff, line[LINELEN], *text;
gboolean band_structure;
FILE *fp;

if (g_strrstr(filename, ".cell") != NULL)
	return read_cell(filename, model);


fp = fopen(filename, "rt");
if (!fp)
  return(1);

model->construct_pbc = TRUE;
model->periodic = 3;
model->fractional = TRUE;

band_structure=FALSE;/*switch for bandstructure reading*/

frame=0;
while (!fgetline(fp, line))
  {
  /* find relevent units */
  if (g_strrstr(line, "force unit") != NULL)
    {
    buff = g_strsplit(line, ":", 2);
    g_strstrip(*(buff+1));
    property_add_ranked(0, "Gradient Units", *(buff+1), model);
    g_strfreev(buff);
    }
  if (g_strrstr(line, "pressure unit") != NULL)
    {
    buff = g_strsplit(line, ":", 2);
    g_strstrip(*(buff+1));
    property_add_ranked(0, "Pressure Units", *(buff+1), model);
    g_strfreev(buff);
    }
  /* store selected parameters in hash table */
  if (g_ascii_strncasecmp(line, " type of calculation", 20) == 0)
    {
    buff = g_strsplit(line, ":", 2);
    g_strstrip(*(buff+1));
    if (g_ascii_strncasecmp(*(buff+1), "geometry optimization", 21) == 0)
      property_add_ranked(2, "Calculation", "Optimisation", model);
    else if (g_ascii_strncasecmp(*(buff+1), "single point energy", 19) == 0)
      property_add_ranked(2, "Calculation", "Single Point", model);
    else if (g_ascii_strncasecmp(*(buff+1), "band structure", 14) == 0){
      property_add_ranked(2, "Calculation", "Band Structure", model);
      band_structure=TRUE;
    }else
      property_add_ranked(2, "Calculation", *(buff+1), model);
    g_strfreev(buff);
    }
  if (g_ascii_strncasecmp(line, " using functional", 17) == 0)
    {
    buff = g_strsplit(line, ":", 2);
    g_strstrip(*(buff+1));
    if (g_ascii_strncasecmp(*(buff+1), "Perdew Burke Ernzerhof", 22) == 0)
      property_add_ranked(7, "Functional", "PBE", model);
    else if (g_ascii_strncasecmp(*(buff+1), "Perdew Wang (1991)", 18) == 0)
      property_add_ranked(7, "Functional", "PW91", model);
    else
      property_add_ranked(7, "Functional", *(buff+1), model);
    g_strfreev(buff);
    }
  if (g_ascii_strncasecmp(line, " plane wave basis set cut-off", 29) == 0)
    {
    buff = g_strsplit(line, ":", 2);
    g_strstrip(*(buff+1));
    text = format_value_and_units(*(buff+1), 1);
    property_add_ranked(6, "Basis Cutoff", text, model);
    g_free(text);
    }  
  if (g_ascii_strncasecmp(line, " size of standard grid", 22) == 0)
    {
    buff = g_strsplit(line, ":", 2);
    g_strstrip(*(buff+1));
    text = format_value_and_units(*(buff+1), 2);
    property_add_ranked(10, "Grid", text, model);
    g_free(text);
    }  
 
  /* read coordinates */
	if ((g_strrstr(line, "Unit Cell") != NULL) || (g_strrstr(line, "Cell Contents") != NULL))
		{
		/* go through all frames to count them */
		add_frame_offset(fp, model);
		read_castep_out_block(fp, model);
		frame++;
    }
  }

fclose(fp);/*for some reason it seems that castep files are not closed after read?*/

if(band_structure) {/*try to read band structure information, if any.*/
	gint ispin=0;
	gint ik=0;
	gint ib=0;
	gdouble xi,yi,zi;
	gdouble *x,*y,*z;
	gdouble ev,evmin,evmax;
	gdouble dist=0.;
	gchar *ptr=NULL;
	gchar *band_filename=g_strdup (filename);
	ptr=g_strrstr (band_filename,".castep");
	if(ptr==NULL) goto castep_done;
	sprintf(ptr,".bands");
	fp = fopen(band_filename, "rt");
	if(!fp) goto castep_done;
	/*get number of kpoints*/
	if(fgetline(fp, line)) goto castep_done;
	sscanf(line,"Number of k-points %d",&(model->nkpoints));
	if(fgetline(fp, line)) goto castep_done;
	sscanf(line,"Number of spin components %d",&ispin);
	if(ispin>1) model->spin_polarized=TRUE;
	if(fgetline(fp, line)) goto castep_done;
	/*skip number of electrons*/
	if(fgetline(fp, line)) goto castep_done;
	sscanf(line,"Number of eigenvalues %d",&(model->nbands));
	if(fgetline(fp, line)) goto castep_done;
	sscanf(line,"Fermi energy (%*[^)]) %lf",&(model->efermi));
	while (!fgetline(fp, line)) if (g_strrstr(line, "K-point") != NULL) break;
	if(!line) goto castep_done;
	model->kpts_d=g_malloc((model->nkpoints)*sizeof(gdouble));
	x=g_malloc((model->nkpoints)*sizeof(gdouble));
	y=g_malloc((model->nkpoints)*sizeof(gdouble));
	z=g_malloc((model->nkpoints)*sizeof(gdouble));
	model->band_up=g_malloc((model->nkpoints*model->nbands)*sizeof(gdouble));
	if(model->spin_polarized) model->band_down=g_malloc((model->nkpoints*model->nbands)*sizeof(gdouble));
	evmin=0.;evmax=0.;
	while(line){
		sscanf(line,"K-point  %i %lf %lf %lf %*f",&ik,&xi,&yi,&zi);
		ik--;
		/*kpts are *not* in order*/
		x[ik]=xi;
		y[ik]=yi;
		z[ik]=zi;
		/*read component 1*/
		if(fgetline(fp, line)) break;
		/* line should contain "Spin component 1", skipping */
		ib=0;
		while(ib<model->nbands) {
			if(fgetline(fp, line)) break;
			sscanf(line," %lf",&(ev));
			if(ev<evmin) evmin=ev;
			if(ev>evmax) evmax=ev;
			model->band_up[ik*(model->nbands)+ib]=ev;
			ib++;
		}
		if(fgetline(fp, line)) break;
		if(model->spin_polarized){/*read component 2*/
			/* line should contain "Spin component 2", skipping */
			ib=0;
			while(ib<model->nbands) {
				if(fgetline(fp, line)) break;
				sscanf(line," %lf",&(ev));
				if(ev<evmin) evmin=ev;
				if(ev>evmax) evmax=ev;
				model->band_down[ik*(model->nbands)+ib]=ev;
				ib++;
			}
		}
	}
	xi=x[0];yi=y[0];zi=z[0];dist=0.;
	for(ik=0;ik<model->nkpoints;ik++){
		/*calculate distances*/
		dist+=sqrt((x[ik]-xi)*(x[ik]-xi)+(y[ik]-yi)*(y[ik]-yi)+(z[ik]-zi)*(z[ik]-zi));
		model->kpts_d[ik]=dist;
		xi=x[ik];yi=y[ik];zi=z[ik];
	}
	g_free(x);
	g_free(y);
	g_free(z);
	fclose(fp);
	g_free(band_filename);
/* BAND reading is done, now let's calculate DOS
 * because CASTEP does not provide it directly..*/
  if(model->nbands>2)
  {
/*we want DOS_NB DOS using a gaussian smearing of SIGMA (in au)*/
#define DOS_NB 1000
#define SIGMA 0.005
gdouble emin=((int)(evmin*10.))/10.;
gdouble emax=((int)(0.5+evmax*10.))/10.;
gdouble dos;
gdouble dos_dn;
gdouble de;
gdouble delta;
int i;
	model->ndos=DOS_NB;
	model->dos_eval=g_malloc(DOS_NB*sizeof(gdouble));
	model->dos_spin_up=g_malloc(DOS_NB*sizeof(gdouble));
	if(model->spin_polarized) model->dos_spin_down=g_malloc(DOS_NB*sizeof(gdouble));
	for(i=0;i<DOS_NB;++i){
		model->dos_eval[i]=emin+i*(emax-emin)/DOS_NB;
		dos=0.;dos_dn=0.;
		for(ik=0;ik<model->nkpoints;ik++)
			for(ib=0;ib<model->nbands;ib++){
				de=model->dos_eval[i]-model->band_up[ik*(model->nbands)+ib];
				de=de/SIGMA;
				delta=exp(-1.0*de*de)/(sqrt(PI));/*GAUSSIAN SMEARING*/
				dos+=delta/SIGMA;
				if(model->spin_polarized){/*process spin down*/
					de=model->dos_eval[i]-model->band_down[ik*(model->nbands)+ib];
					de=de/SIGMA;
					delta=exp(-1.0*de*de)/(sqrt(PI));/*GAUSSIAN SMEARING*/
					dos_dn+=delta/SIGMA;
				}
			}
		model->dos_spin_up[i]=dos;
		if(model->spin_polarized) model->dos_spin_down[i]=dos_dn;
	}
  }
  }
/* done */

castep_done:
strcpy(model->filename, filename);
g_free(model->basename);
model->basename = parse_strip(filename);
model->num_frames = model->cur_frame = frame;
model->cur_frame--;

model_prep(model);

return(0);
}


/* VEZ write in CASTEP cell file */

gint write_castep_cell(gchar *filename, struct model_pak *model)
{
gdouble c;
gdouble x[3];
GSList *list;
struct core_pak *core;
gint i = 0;
FILE *fp;

/* checks */
g_return_val_if_fail(model != NULL, 1);
g_return_val_if_fail(filename != NULL, 2);

/* open the file */
fp = fopen(filename, "wt");
if (!fp)
  return(3);

if (model->periodic)
  {
  if (model->periodic == 2)
    /* saving a surface as a 3D model - make depth large enough to fit everything */
    c =  2.0*model->rmax;
  else
    c = model->pbc[2];

  if ((180.0*model->pbc[3]/PI<0.0) || (180.0*model->pbc[3]/PI>180.0))
	gui_text_show(ERROR, "Cell angle alpha is not in 0-180 degree range, please correct manually\n");
  if ((180.0*model->pbc[4]/PI<0.0) || (180.0*model->pbc[4]/PI>180.0))
	gui_text_show(ERROR, "Cell angle beta is not in 0-180 degree range, please correct manually\n");
  if ((180.0*model->pbc[5]/PI<0.0) || (180.0*model->pbc[5]/PI>180.0))
	gui_text_show(ERROR, "Cell angle gamma is not in 0-180 degree range, please correct manually\n");


fprintf(fp,"%%BLOCK LATTICE_ABC\n");
fprintf(fp,"%f  %f   %f\n",model->pbc[0], model->pbc[1], c);
fprintf(fp,"%f  %f   %f\n",fabs(180.0*model->pbc[3]/PI),fabs(180.0*model->pbc[4]/PI),fabs(180.0*model->pbc[5]/PI));
fprintf(fp,"%%ENDBLOCK LATTICE_ABC\n\n");

}

/*
The ATOM records
*/

fprintf(fp,"%%BLOCK POSITIONS_ABS\n");

for (list=model->cores ; list ; list=g_slist_next(list))
  {
  core = list->data;

/* everything is cartesian after latmat mult */
  ARR3SET(x, core->x);
  vecmat(model->latmat, x);
  i++;
  fprintf(fp,"%s   %f   %f   %f\n",elements[core->atom_code].symbol,x[0], x[1], x[2]);

/*     fprintf(fp,"ATOM  %5d %-4.4s %-3.3s  %c %3d %8.3f%8.3f%8.3f%6.2f%6.2f      %-4s%2s%-2s\n",
              i, core->atom_label, core->res_name, core->chain, core->res_no, x[0], x[1], x[2], core->sof, 0.0, "", elements[core->atom_code].symbol, "");*/
}

fprintf(fp,"%%ENDBLOCK POSITIONS_ABS\n\n");


/* fixed atoms */
gint constrains_counter=1;
gint atomic_counter[200];
int qq;
for (qq=0; qq<200; qq++)
	atomic_counter[qq]=1;

fprintf(fp,"%%BLOCK IONIC_CONSTRAINTS\n");
for (list=model->cores ; list ; list=g_slist_next(list))
  {
	core = list->data;
	if (core->ghost) {
		fprintf(fp,"%d %s %d 1 0 0\n",constrains_counter, elements[core->atom_code].symbol, atomic_counter[core->atom_code]);
		constrains_counter++;
		fprintf(fp,"%d %s %d 0 1 0\n",constrains_counter, elements[core->atom_code].symbol, atomic_counter[core->atom_code]);
		constrains_counter++;
		fprintf(fp,"%d %s %d 0 0 1\n",constrains_counter, elements[core->atom_code].symbol, atomic_counter[core->atom_code]);
		constrains_counter++;
	}
	atomic_counter[core->atom_code]++;
  }
fprintf(fp,"%%ENDBLOCK IONIC_CONSTRAINTS\n\n");


fprintf(fp,"KPOINTS_MP_SPACING 0.05\n\nFIX_ALL_CELL : TRUE\n\nFIX_COM : FALSE\n\nSYMMETRY_GENERATE\n");

fclose(fp);
return(0);
}


