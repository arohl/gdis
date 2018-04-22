/*
Copyright (C) 2003 by Andrew Lloyd Rohl

a.rohl@curtin.edu.au

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
#include "model.h"
#include "file.h"
#include "parse.h"
#include "matrix.h"
#include "interface.h"

/* main structures */
extern struct sysenv_pak sysenv;
extern struct elem_pak elements[];


/*
The ATOM records present the atomic coordinates for standard residues. They
also present the occupancy and temperature factor for each atom. Heterogen
coordinates use the HETATM record type. The element symbol is always present
on each ATOM record; segment identifier and charge are optional.
COLUMNS        DATA TYPE       FIELD         DEFINITION
---------------------------------------------------------------------------------
 1 -  6        Record name     "ATOM  "
 7 - 11        Integer         serial        Atom serial number.
13 - 16        Atom            name          Atom name.
17             Character       altLoc        Alternate location indicator.
18 - 20        Residue name    resName       Residue name.
22             Character       chainID       Chain identifier.
23 - 26        Integer         resSeq        Residue sequence number.
27             AChar           iCode         Code for insertion of residues.
31 - 38        Real(8.3)       x             Orthogonal coordinates for X in Angstroms.
39 - 46        Real(8.3)       y             Orthogonal coordinates for Y in Angstroms.
47 - 54        Real(8.3)       z             Orthogonal coordinates for Z in Angstroms.
55 - 60        Real(6.2)       occupancy     Occupancy.
61 - 66        Real(6.2)       tempFactor    Temperature factor.
73 - 76        LString(4)      segID         Segment identifier, left-justified.
77 - 78        LString(2)      element       Element symbol, right-justified.
79 - 80        LString(2)      charge        Charge on the atom.
*/

/*
Note that in common with GPTA from Paolo Raiteri, we have used columns 81-100 for accurate charges.
*/

/* FIXME - deal with res_name , chain values of NULL */
void write_pdb_core(FILE *fp, gint i, struct core_pak *core, struct model_pak *model)
{
gdouble x[3];
gchar *res_name;

ARR3SET(x, core->x);
vecmat(model->latmat, x);

if (!core->res_name)
  res_name = "UNK";
else
  res_name = core->res_name;

/* NEW - enforce string length minimums */
fprintf(fp,"ATOM  %5d %-4.4s %3.3s %c %3d    %8.3f%8.3f%8.3f%6.2f%6.2f      %-4.4s%2.2s%-2.2s",
            i, core->atom_label, res_name, core->chain, core->res_no,
            x[0], x[1], x[2], core->sof, 0.0, "", elements[core->atom_code].symbol, "");
if (!core->lookup_charge)
  fprintf(fp,"%12.8f", core->charge);
fprintf(fp,"\n");
}


/******************/
/* PDB writing */
/******************/
gint write_pdb(gchar *filename, struct model_pak *model)
{
gdouble c;
GSList *list, *list2;
struct bond_pak *bond;
struct core_pak *core, *core2;
gint i;
FILE *fp;

/* checks */
g_return_val_if_fail(model != NULL, 1);
g_return_val_if_fail(filename != NULL, 2);

/* open the file */
fp = fopen(filename, "wt");
if (!fp)
  return(3);

/* TODO - write multiple frames if it has them */

/*
  The CRYST1 record presents the unit cell parameters, space group, and Z value.
  If the structure was not determined by crystallographic means, CRYST1 simply 
  defines a unit cube. 
  COLUMNS       DATA TYPE      FIELD         DEFINITION
  -------------------------------------------------------------
    1 -  6       Record name    "CRYST1"
    7 - 15       Real(9.3)      a             a (Angstroms).
   16 - 24       Real(9.3)      b             b (Angstroms).
   25 - 33       Real(9.3)      c             c (Angstroms).
   34 - 40       Real(7.2)      alpha         alpha (degrees).
   41 - 47       Real(7.2)      beta          beta (degrees).
   48 - 54       Real(7.2)      gamma         gamma (degrees).
   56 - 66       LString        sGroup        Space group.
   67 - 70       Integer        z             Z value.
*/

if (model->periodic)
  {
  if (model->periodic == 2)
    /* saving a surface as a 3D model - make depth large enough to fit everything */
    c =  2.0*model->rmax;
  else
    c = model->pbc[2];

if (model->sginfo.spacename)
  fprintf(fp,"CRYST1%9.3f%9.3f%9.3f%7.2f%7.2f%7.2f %-11s%4d\n",
              model->pbc[0], model->pbc[1], c,
              180.0*model->pbc[3]/PI,
              180.0*model->pbc[4]/PI,
              180.0*model->pbc[5]/PI,
              model->sginfo.spacename, 0);
else
  fprintf(fp,"CRYST1%9.3f%9.3f%9.3f%7.2f%7.2f%7.2f %-11s%4d\n",
              model->pbc[0], model->pbc[1], c,
              180.0*model->pbc[3]/PI,
              180.0*model->pbc[4]/PI,
              180.0*model->pbc[5]/PI,
              "P1", 0);

  }

i=1;
for (list=model->cores ; list ; list=g_slist_next(list))
  {
  write_pdb_core(fp, i++, list->data, model);
  }

/* CURRENT */
/* connectivity */
i=1;
for (list=model->cores ; list ; list=g_slist_next(list))
  {
  core = list->data;

  if (core->bonds)
    {
fprintf(fp, "CONECT");

fprintf(fp, " %4d", i);

    for (list2=core->bonds ; list2 ; list2=g_slist_next(list2))
      {
      bond = list2->data;

      if (bond->status == NORMAL)
        {
        if (connect_split(bond))
          continue;

        if (core == bond->atom1)
          core2 = bond->atom2;
        else
          core2 = bond->atom1;
        fprintf(fp, " %4d", 1+g_slist_index(model->cores, core2));
        }
      }

fprintf(fp, "\n");

    }
  i++;
  }

fprintf(fp, "END\n");

fclose(fp);
return(0);
}

/***********************************************/
/* read a line of FORTRAN output padding to    */
/* PDBLINELEN characters (with nulls) if       */ 
/* required                                    */
/***********************************************/
#define PDBLINELEN 100
gint fort_read_line(FILE *fp, gchar **line)
{
gint ret_val;
gchar the_line[PDBLINELEN+1];

ret_val = fgetline(fp, the_line);
if (!ret_val)
  {
  *line = g_strndup(the_line, PDBLINELEN);
  }
return(ret_val);
}

/**********************************************/
/* read a gint from a line of FORTRAN output  */
/**********************************************/
gint fort_read_gint(char *line, gint start_col, gint end_col)
{
gchar *field;
gdouble value;

g_return_val_if_fail(start_col > 0, 0.0);
g_return_val_if_fail(start_col <= end_col, 0.0);
field = g_strndup(line+start_col-1, end_col-start_col+1);
value = (gint) str_to_float(field);
g_free(field);
return(value);
}

/*************************************************/
/* read a gdouble from a line of FORTRAN output  */
/*************************************************/
gdouble fort_read_gdouble(char *line, gint start_col, gint end_col)
{
gchar *field;
gdouble value;

g_return_val_if_fail(start_col > 0, 0.0);
g_return_val_if_fail(start_col <= end_col, 0.0);
field = g_strndup(line+start_col-1, end_col-start_col+1);
value = str_to_float(field);
g_free(field);
return(value);
}

/*************************************************/
/* read a string from a line of FORTRAN output  */
/*************************************************/
gchar *fort_read_string(char *line, gint start_col, gint end_col)
{
  gchar *field;
  
  g_return_val_if_fail(start_col > 0, "");
  g_return_val_if_fail(start_col <= end_col, "");
  field = g_strndup(line+start_col-1, end_col-start_col+1);
  return(field);
}

/*************************************************/
/* read a character from a line of FORTRAN output  */
/*************************************************/
gchar fort_read_char(char *line, gint start_col)
{
  g_return_val_if_fail(start_col > 0, "");
  return(line[start_col-1]);
}

/*************************************************/
/* read a PDB block into the model structure     */
/*************************************************/
gint read_pdb_block(FILE *fp, struct model_pak *data)
{
gchar *line, *atom_name, *res_name, *element;
struct core_pak *core;
GSList *clist;

clist = data->cores;

/* loop while there's data */
while (!fort_read_line(fp, &line))
  {
  /* 3-D structure? */
  if (g_ascii_strncasecmp("CRYST1", line, 6) == 0)
      /*
        The CRYST1 record presents the unit cell parameters, space group, and Z value.
        If the structure was not determined by crystallographic means, CRYST1 simply 
        defines a unit cube. 
        COLUMNS       DATA TYPE      FIELD         DEFINITION
        -------------------------------------------------------------
          1 -  6       Record name    "CRYST1"
          7 - 15       Real(9.3)      a             a (Angstroms).
         16 - 24       Real(9.3)      b             b (Angstroms).
         25 - 33       Real(9.3)      c             c (Angstroms).
         34 - 40       Real(7.2)      alpha         alpha (degrees).
         41 - 47       Real(7.2)      beta          beta (degrees).
         48 - 54       Real(7.2)      gamma         gamma (degrees).
         56 - 66       LString        sGroup        Space group.
         67 - 70       Integer        z             Z value.
    */
    {
    data->periodic = 3;
    data->pbc[0] = fort_read_gdouble(line, 7, 15);
    data->pbc[1] = fort_read_gdouble(line, 16, 24);
    data->pbc[2] = fort_read_gdouble(line, 25, 33);
    data->pbc[3] = fort_read_gdouble(line, 34, 40) * PI/180.0;
    data->pbc[4] = fort_read_gdouble(line, 41, 47) * PI/180.0;
    data->pbc[5] = fort_read_gdouble(line, 48, 54) * PI/180.0;
    data->sginfo.spacename = fort_read_string(line, 56, 66);
    
    /* indicate that name should used in lookup (IF a spacegroup was found) */
    if (strlen(data->sginfo.spacename))
      data->sginfo.spacenum = -1;
    else
      {
      g_free(data->sginfo.spacename);
      data->sginfo.spacename = g_strdup("P 1");
      }
    }
  else if ((g_ascii_strncasecmp("ATOM", line, 4) == 0) || (g_ascii_strncasecmp("HETATM", line, 6) == 0)) 
    {
    /*
    The ATOM records present the atomic coordinates for standard residues. They
    also present the occupancy and temperature factor for each atom. Heterogen
    coordinates use the HETATM record type. The element symbol is always present
    on each ATOM record; segment identifier and charge are optional.
    COLUMNS        DATA TYPE       FIELD         DEFINITION
    ---------------------------------------------------------------------------------
     1 -  6        Record name     "ATOM  "
     7 - 11        Integer         serial        Atom serial number.
    13 - 16        Atom            name          Atom name.
    17             Character       altLoc        Alternate location indicator.
    18 - 20        Residue name    resName       Residue name.
    22             Character       chainID       Chain identifier.
    23 - 26        Integer         resSeq        Residue sequence number.
    27             AChar           iCode         Code for insertion of residues.
    31 - 38        Real(8.3)       x             Orthogonal coordinates for X in Angstroms.
    39 - 46        Real(8.3)       y             Orthogonal coordinates for Y in Angstroms.
    47 - 54        Real(8.3)       z             Orthogonal coordinates for Z in Angstroms.
    55 - 60        Real(6.2)       occupancy     Occupancy.
    61 - 66        Real(6.2)       tempFactor    Temperature factor.
    73 - 76        LString(4)      segID         Segment identifier, left-justified.
    77 - 78        LString(2)      element       Element symbol, right-justified.
    79 - 80        LString(2)      charge        Charge on the atom.

    Note that in common with GPTA from Paolo Raiteri, we have used columns 81-100 for accurate charges.
    */
    atom_name = fort_read_string(line, 13, 16);
    if (clist)
      {
      core = clist->data;
      clist = g_slist_next(clist);
      }
    else
      {
      element = fort_read_string(line, 77, 78);
      element = g_strchug(element);
/*    printf("\"%s\"\n", element);*/
      if (*element != '\0')
        core = core_new(element, atom_name, data);
      else
        core = new_core(atom_name, data);
      data->cores = g_slist_append(data->cores, core);
      }
    g_free(atom_name);
    
    g_free(core->res_name);
    res_name = fort_read_string(line, 18, 20);
    core->res_name = g_strdup(res_name);
    g_free(res_name);
  
    core->chain = fort_read_char(line, 22);
    if ((core->chain >= 'A') && (core->chain <= 'Z'))
      core->region = core->chain - 'A';

    core->res_no = fort_read_gint(line, 23, 26);
    
    core->x[0] = fort_read_gdouble(line, 31, 38);
    core->x[1] = fort_read_gdouble(line, 39, 46);
    core->x[2] = fort_read_gdouble(line, 47, 54);
    
    core->sof = fort_read_gdouble(line, 55, 60);
    if (line[81] != '\0')
      {
      core->lookup_charge = FALSE;
      core->charge = fort_read_gdouble(line, 81, 100);
      }
    }
  else if (g_ascii_strncasecmp("ENDMDL", line, 6) == 0)
    {
    /*
    The ENDMDL records are paired with MODEL records to group individual structures
    found in a coordinate entry. 
    COLUMNS         DATA TYPE        FIELD           DEFINITION
    ------------------------------------------------------------------
      1 -  6         Record name      "ENDMDL"
    */
    return(0);
    }
  else if (g_ascii_strncasecmp("END", line, 3) == 0)
    {
/* end of a frame? */
    return(0);
    }

  g_free(line);
  }
  return(-1);
}

/*********************/
/* PDB frame read    */
/*********************/
gint read_pdb_frame(FILE *fp, struct model_pak *data)
{
/* frame overwrite */
read_pdb_block(fp, data);
return(0);
}

/******************/
/* PDB reading    */
/******************/
#define DEBUG_READ_PDB 0
gint read_pdb(gchar *filename, struct model_pak *data)
{
gint frame=0;
FILE *fp;

/* checks */
g_return_val_if_fail(data != NULL, 1);
g_return_val_if_fail(filename != NULL, 2);

/* initialisers */
data->periodic = 0;
data->protein = TRUE;
data->fractional = FALSE; 

fp = fopen(filename, "rt");
if (!fp)
  return(3);
add_frame_offset(fp, data);

while (read_pdb_block(fp, data) == 0)
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

                            /*
                                The MODEL record specifies the model serial number when multiple structures
                                are presented in a single coordinate entry, as is often the case with structures
                                determined by NMR. 
                                COLUMNS       DATA TYPE      FIELD         DEFINITION
                                ----------------------------------------------------------------------
                                1 -  6       Record name    "MODEL "
                                11 - 14       Integer        serial        Model serial number.
                            */





//	 		if ((line80.startsWith("ATOM") || line80.startsWith("HETATM"))) {
	 		    /*
				The ATOM records present the atomic coordinates for standard residues. They
				also present the occupancy and temperature factor for each atom. Heterogen
				coordinates use the HETATM record type. The element symbol is always present
				on each ATOM record; segment identifier and charge are optional.
				COLUMNS        DATA TYPE       FIELD         DEFINITION
				---------------------------------------------------------------------------------
				 1 -  6        Record name     "ATOM  "
				 7 - 11        Integer         serial        Atom serial number.
				13 - 16        Atom            name          Atom name.
				17             Character       altLoc        Alternate location indicator.
				18 - 20        Residue name    resName       Residue name.
				22             Character       chainID       Chain identifier.
				23 - 26        Integer         resSeq        Residue sequence number.
				27             AChar           iCode         Code for insertion of residues.
				31 - 38        Real(8.3)       x             Orthogonal coordinates for X in Angstroms.
				39 - 46        Real(8.3)       y             Orthogonal coordinates for Y in Angstroms.
				47 - 54        Real(8.3)       z             Orthogonal coordinates for Z in Angstroms.
				55 - 60        Real(6.2)       occupancy     Occupancy.
				61 - 66        Real(6.2)       tempFactor    Temperature factor.
				73 - 76        LString(4)      segID         Segment identifier, left-justified.
				77 - 78        LString(2)      element       Element symbol, right-justified.
				79 - 80        LString(2)      charge        Charge on the atom.
			    */

/*	 			serial = Integer.valueOf(line80.substring(6, 11).trim()).intValue();
	 			if (countOnly) {
	 				if (serial > nAtoms)
	 					nAtoms = serial;
	 			}
	 			else {
	 				name = line80.substring(12, 16).trim();
	    			x = Double.valueOf(line80.substring(30, 38).trim()).doubleValue();
	    			y = Double.valueOf(line80.substring(38, 46).trim()).doubleValue();
	    			z = Double.valueOf(line80.substring(46, 54).trim()).doubleValue();
	 				atomlist[serial] = new Atom (name, x, y, z, blackBackground);
	 				sortindex[serial] = serial;
	 			}
	 		} */

                            /*
                                The ENDMDL records are paired with MODEL records to group individual structures
                                found in a coordinate entry. 
                                COLUMNS         DATA TYPE        FIELD           DEFINITION
                                ------------------------------------------------------------------
                                 1 -  6         Record name      "ENDMDL"
                            */

				/*
				The CONECT records specify connectivity between atoms for which coordinates
				are supplied. The connectivity is described using the atom serial number as
				found in the entry. CONECT records are mandatory for HET groups (excluding
				water) and for other bonds not specified in the standard residue
				connectivity table which involve atoms in standard residues (see Appendix 4
				for the list of standard residues).
				COLUMNS         DATA TYPE        FIELD           DEFINITION
				---------------------------------------------------------------------------------
				 1 -  6         Record name      "CONECT"
				 7 - 11         Integer          serial          Atom serial number
				12 - 16         Integer          serial          Serial number of bonded atom
				17 - 21         Integer          serial          Serial number of bonded atom
				22 - 26         Integer          serial          Serial number of bonded atom
				27 - 31         Integer          serial          Serial number of bonded atom
				32 - 36         Integer          serial          Serial number of hydrogen bonded atom
				37 - 41         Integer          serial          Serial number of hydrogen bonded atom
				42 - 46         Integer          serial          Serial number of salt bridged atom
				47 - 51         Integer          serial          Serial number of hydrogen bonded atom
				52 - 56         Integer          serial          Serial number of hydrogen bonded atom
				57 - 61         Integer          serial          Serial number of salt bridged atom
				*/
