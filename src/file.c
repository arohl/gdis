/*
Copyright (C) 2003 by Sean David Fleming

sean@ivec.org

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
#include <strings.h>
#include <ctype.h>
#include <stdlib.h>
#include <time.h>
#include <unistd.h>
#define GTK_DISABLE_DEPRECATED
#include <glib/gstdio.h>
#include <sys/stat.h>
#include <sys/types.h>

#include "gdis.h"
#include "coords.h"
#include "error.h"
#include "file.h"
#include "parse.h"
#include "matrix.h"
#include "model.h"
#include "space.h"
#include "zone.h"
#include "render.h"
#include "select.h"
#include "gui_shorts.h"
#include "interface.h"
#include "dialog.h"
#include "opengl.h"

#define DEBUG_MORE 0
#define MAX_KEYS 15

/* main structures */
extern struct sysenv_pak sysenv;
extern struct elem_pak elements[];

/***************************************/
/* setup the recognized file type list */
/***************************************/
#define DEBUG_FILE_INIT 0
void file_init(void)
{
GSList *rlist=NULL;
struct file_pak *file_data;

/* NEW - build a recognized image read format list */
/* also used in picture writing as it avoids overlap */
/* problems when getting the file_data structure by extension */ 
#define PICTURE_SUPPORT 0
#if PICTURE_SUPPORT
for (list=gdk_pixbuf_get_formats() ; list ; list=g_slist_next(list))
  {
  gint i;
  gchar **ext;

  ext = gdk_pixbuf_format_get_extensions(list->data);

  i = 0;
  while (*(ext+i))
    {
/*
    if (gdk_pixbuf_format_is_writable(list->data))
      wlist = g_slist_prepend(wlist, g_strdup(*(ext+i)));
*/
    rlist = g_slist_prepend(rlist, g_strdup(*(ext+i)));
    i++;
    }
  }
#endif

#if DEBUG_FILE_INIT
printf("read: ");
for (list=rlist ; list ; list=g_slist_next(list))
  {
  printf("[%s] ", (gchar *) list->data);
  }
printf("\n");
#endif

/* build the recognized file type list */
sysenv.file_list = NULL;

/* supported file type */
file_data = g_malloc(sizeof(struct file_pak));
file_data->id = DATA;                           /* unique identifier */ 
file_data->group = DATA;                        /* used to group inp/out types */
file_data->menu = TRUE;                         /* include in menu listing */
file_data->label = g_strdup("All known types"); /* text info for the user */
file_data->ext = NULL;                          /* extension matching */
file_data->write_file = NULL;                   /* file creation */
file_data->read_file = NULL;                    /* file reading */
file_data->read_frame = NULL;                   /* frame reading */
sysenv.file_list = g_slist_prepend(sysenv.file_list, file_data);

/* supported file type */
file_data = g_malloc(sizeof(struct file_pak));
file_data->id = ABINIT_OUT;
file_data->group = ABINIT;
file_data->menu = FALSE;
file_data->label = g_strdup("ABINIT output");
file_data->ext = NULL;
file_data->ext = g_slist_prepend(file_data->ext, "about");
file_data->ext = g_slist_prepend(file_data->ext, "abot");
file_data->write_file = NULL;
file_data->read_file = read_about;
file_data->read_frame = read_about_frame;
sysenv.file_list = g_slist_prepend(sysenv.file_list, file_data);

/* supported file type */
file_data = g_malloc(sizeof(struct file_pak));
file_data->id = MD_ANALYSIS;
file_data->group = MD_ANALYSIS;
file_data->menu = FALSE;
file_data->label = g_strdup("Analysis");
file_data->ext = NULL;
file_data->ext = NULL;
/* import/export handled in analysis as a special case */
file_data->write_file = NULL;
file_data->read_file = NULL;
file_data->read_frame = NULL;                   /* frame reading */
sysenv.file_list = g_slist_prepend(sysenv.file_list, file_data);

/* supported file type  */
file_data = g_malloc(sizeof(struct file_pak));
file_data->id = BIOGRAF;
file_data->group = BIOGRAF;
file_data->menu = TRUE;
file_data->label = g_strdup("Biograf");                 
file_data->ext = NULL;
file_data->ext = g_slist_prepend(file_data->ext, "bgf");
file_data->write_file = write_bgf;
file_data->read_file = read_bgf;
file_data->read_frame = NULL;                   /* frame reading */
sysenv.file_list = g_slist_prepend(sysenv.file_list, file_data);  

/* supported file type */
file_data = g_malloc(sizeof(struct file_pak));
file_data->id = BIOSYM;
file_data->group = BIOSYM;
file_data->menu = TRUE;
file_data->label = g_strdup("Biosym");
file_data->ext = NULL;
file_data->ext = g_slist_prepend(file_data->ext, "car");
file_data->ext = g_slist_prepend(file_data->ext, "cor");
file_data->ext = g_slist_prepend(file_data->ext, "arc");
file_data->write_file = write_arc;
file_data->read_file = read_arc;
file_data->read_frame = read_arc_frame;
sysenv.file_list = g_slist_prepend(sysenv.file_list, file_data);

/* supported file type */
file_data = g_malloc(sizeof(struct file_pak));
file_data->id = CASTEP;
file_data->group = CASTEP;
file_data->menu = TRUE;
file_data->label = g_strdup("CASTEP");
file_data->ext = NULL;
file_data->ext = g_slist_append(file_data->ext, "cell");
file_data->write_file = write_castep_cell;
file_data->read_file = NULL;
file_data->read_frame = NULL;
sysenv.file_list = g_slist_prepend(sysenv.file_list, file_data);

/* supported file type */
file_data = g_malloc(sizeof(struct file_pak));
file_data->id = CASTEP_OUT;
file_data->group = CASTEP;
file_data->menu = FALSE;
file_data->label = g_strdup("CASTEP output");
file_data->ext = NULL;
file_data->ext = g_slist_append(file_data->ext, "castep");
file_data->write_file = NULL;
file_data->read_file = read_castep_out;
file_data->read_frame = read_castep_out_frame;
sysenv.file_list = g_slist_prepend(sysenv.file_list, file_data);

/* supported file type */
file_data = g_malloc(sizeof(struct file_pak));
file_data->id = CIF;
file_data->group = CIF;
file_data->menu = TRUE;
file_data->label = g_strdup("CIF");
file_data->ext = NULL;
file_data->ext = g_slist_prepend(file_data->ext, "cif");
file_data->write_file = write_cif;
file_data->read_file = read_cif;
file_data->read_frame = NULL;                   /* frame reading */
sysenv.file_list = g_slist_prepend(sysenv.file_list, file_data);

/* NEW: monty crystal graph file format */
/* supported file type  */
file_data = g_malloc(sizeof(struct file_pak));
file_data->id = CRYSTAL_GRAPH;
file_data->group = CRYSTAL_GRAPH;
file_data->menu = TRUE;
file_data->label = g_strdup("Crystal Graph Files");                 
file_data->ext = NULL;
file_data->ext = g_slist_prepend(file_data->ext, "cgf");
file_data->write_file = write_cgf;
file_data->read_file = read_cgf;
file_data->read_frame = NULL;                   /* frame reading */
sysenv.file_list = g_slist_prepend(sysenv.file_list, file_data);

/* supported file type */
file_data = g_malloc(sizeof(struct file_pak));
file_data->id = CSSR;
file_data->group = CSSR;
file_data->menu = TRUE;
file_data->label = g_strdup("CSSR");
file_data->ext = NULL;
file_data->ext = g_slist_append(file_data->ext, "cssr");
file_data->write_file = write_cssr;
file_data->read_file = read_cssr;
sysenv.file_list = g_slist_prepend(sysenv.file_list, file_data);

/* supported file type */
file_data = g_malloc(sizeof(struct file_pak));
file_data->id = DIFFAX_INP;
file_data->group = DIFFAX_INP;
file_data->menu = FALSE;
file_data->label = g_strdup("DIFFaX");
file_data->ext = NULL;
file_data->ext = g_slist_prepend(file_data->ext, "dfx");
file_data->write_file = write_diffax;
file_data->read_file = read_diffax;
file_data->read_frame = NULL;                   /* frame reading */
sysenv.file_list = g_slist_prepend(sysenv.file_list, file_data);

/* supported file type */
file_data = g_malloc(sizeof(struct file_pak));
file_data->id = DLP;
file_data->group = DLP;
file_data->menu = TRUE;
file_data->label = g_strdup("DLP");
file_data->ext = NULL;
file_data->ext = g_slist_prepend(file_data->ext, "dlp");
file_data->write_file = write_dlp;
file_data->read_file = read_dlp;
file_data->read_frame = NULL;
sysenv.file_list = g_slist_prepend(sysenv.file_list, file_data);

/* supported file type */
file_data = g_malloc(sizeof(struct file_pak));
file_data->id = DMOL_INPUT;
file_data->group = DMOL_INPUT;
file_data->menu = TRUE;
file_data->label = g_strdup("DMOL");
file_data->ext = NULL;
file_data->ext = g_slist_prepend(file_data->ext, "dmol");
file_data->write_file = write_dmol;
file_data->read_file = read_dmol;
file_data->read_frame = read_dmol_frame;
sysenv.file_list = g_slist_prepend(sysenv.file_list, file_data);

/* supported file type */
file_data = g_malloc(sizeof(struct file_pak));
file_data->id = DLPOLY;
file_data->group = DLPOLY;
file_data->menu = TRUE;
file_data->label = g_strdup("DL_POLY");
file_data->ext = NULL;
file_data->ext = g_slist_append(file_data->ext, "");
file_data->ext = g_slist_append(file_data->ext, "dlpoly");
file_data->write_file = write_dlpoly;
file_data->read_file = read_dlpoly;
file_data->read_frame = read_dlpoly_frame;
sysenv.file_list = g_slist_prepend(sysenv.file_list, file_data);

/* supported file type */
file_data = g_malloc(sizeof(struct file_pak));
file_data->id = GAMESS;
file_data->group = GAMESS;
file_data->menu = TRUE;
file_data->label = g_strdup("GAMESS");
file_data->ext = NULL;
file_data->ext = g_slist_append(file_data->ext, "inp");
file_data->write_file = write_gms;
file_data->read_file = read_gms;
file_data->read_frame = NULL;                   /* frame reading */
sysenv.file_list = g_slist_prepend(sysenv.file_list, file_data);

/* supported file type */
file_data = g_malloc(sizeof(struct file_pak));
file_data->id = GAMESS_OUT;
file_data->group = GAMESS;
file_data->menu = FALSE;
file_data->label = g_strdup("GAMESS Output");
file_data->ext = NULL;
file_data->ext = g_slist_prepend(file_data->ext, "gmout");
file_data->ext = g_slist_prepend(file_data->ext, "gmot");
file_data->write_file = NULL;
file_data->read_file = read_gms_out;
file_data->read_frame = read_gms_out_frame;
sysenv.file_list = g_slist_prepend(sysenv.file_list, file_data);

/* supported file type */
file_data = g_malloc(sizeof(struct file_pak));
file_data->id = GAUSS;
file_data->group = GAUSS;
file_data->menu = TRUE;
file_data->label = g_strdup("GAUSSIAN");
file_data->ext = NULL;
file_data->ext = g_slist_append(file_data->ext, "com");
file_data->write_file = write_gauss;
file_data->read_file = NULL;
file_data->read_frame = NULL;
sysenv.file_list = g_slist_prepend(sysenv.file_list, file_data);

/* supported file type */
file_data = g_malloc(sizeof(struct file_pak));
file_data->id = GAUSS_OUT;
file_data->group = GAUSS;
file_data->menu = FALSE;
file_data->label = g_strdup("GAUSSIAN output");
file_data->ext = NULL;
file_data->ext = g_slist_append(file_data->ext, "log");
file_data->write_file = NULL;
file_data->read_file = read_gauss_out;
file_data->read_frame = read_gauss_out_frame;
sysenv.file_list = g_slist_prepend(sysenv.file_list, file_data);

/* supported file type */
file_data = g_malloc(sizeof(struct file_pak));
file_data->id = QE;
file_data->group = QE;
file_data->menu = TRUE;
file_data->label = g_strdup("QE");
file_data->ext = NULL;
file_data->ext = g_slist_append(file_data->ext, "qein");
file_data->write_file = write_qe;
file_data->read_file = read_qe;
file_data->read_frame = NULL;
sysenv.file_list = g_slist_prepend(sysenv.file_list, file_data);

/* supported file type */
file_data = g_malloc(sizeof(struct file_pak));
file_data->id = QE_OUT;
file_data->group = QE;
file_data->menu = FALSE;
file_data->label = g_strdup("QE output");
file_data->ext = NULL;
file_data->ext = g_slist_append(file_data->ext, "qeout");
file_data->write_file = NULL;
file_data->read_file = read_qe_out;
file_data->read_frame = read_qe_out_frame;
sysenv.file_list = g_slist_prepend(sysenv.file_list, file_data);
    

/* supported file type */
file_data = g_malloc(sizeof(struct file_pak));
file_data->id = GROMACS;
file_data->group = GROMACS;
file_data->menu = FALSE;
file_data->label = g_strdup("GROMACS input");
file_data->ext = NULL;
file_data->ext = g_slist_append(file_data->ext, "gro");
file_data->write_file = write_gromacs;
file_data->read_file = read_gromacs_gro;
file_data->read_frame = NULL;
sysenv.file_list = g_slist_prepend(sysenv.file_list, file_data);

/* supported file type */
file_data = g_malloc(sizeof(struct file_pak));
file_data->id = MORPH;
file_data->group = MORPH;
file_data->menu = TRUE;
file_data->label = g_strdup("GDIS Morphology");
file_data->ext = NULL;
file_data->ext = g_slist_append(file_data->ext, "gmf");
file_data->write_file = write_gmf;
file_data->read_file = read_gmf;
file_data->read_frame = NULL;                   /* frame reading */
sysenv.file_list = g_slist_prepend(sysenv.file_list, file_data);

/* supported file type */
file_data = g_malloc(sizeof(struct file_pak));
file_data->id = GEOMVIEW_OFF;
file_data->group = GEOMVIEW_OFF;
file_data->menu = TRUE;
file_data->label = g_strdup("Geomview OFF");
file_data->ext = NULL;
file_data->ext = g_slist_prepend(file_data->ext, "off");
file_data->write_file = NULL;
file_data->read_file = read_off;
file_data->read_frame = NULL;
sysenv.file_list = g_slist_prepend(sysenv.file_list, file_data);

/* supported file type */
file_data = g_malloc(sizeof(struct file_pak));
file_data->id = GULP;
file_data->group = GULP;
file_data->menu = TRUE;
file_data->label = g_strdup("GULP");
file_data->ext = NULL;
file_data->ext = g_slist_prepend(file_data->ext, "gin");
file_data->ext = g_slist_prepend(file_data->ext, "res");
file_data->write_file = write_gulp;
file_data->read_file = read_gulp;
file_data->read_frame = NULL;                   /* frame reading */
sysenv.file_list = g_slist_prepend(sysenv.file_list, file_data);

/* supported file type */
file_data = g_malloc(sizeof(struct file_pak));
file_data->id = GULPOUT;
file_data->group = GULP;
file_data->menu = FALSE;
file_data->label = g_strdup("GULP output");
file_data->ext = NULL;
file_data->ext = g_slist_prepend(file_data->ext, "got");
file_data->ext = g_slist_prepend(file_data->ext, "gout");
file_data->write_file = NULL;
file_data->read_file = read_gulp_output;
file_data->read_frame = NULL;                   /* frame reading */
sysenv.file_list = g_slist_prepend(sysenv.file_list, file_data);

/* supported file type */
file_data = g_malloc(sizeof(struct file_pak));
file_data->id = MARVIN;
file_data->group = MARVIN;
file_data->menu = TRUE;
file_data->label = g_strdup("MARVIN");
file_data->ext = NULL;
file_data->ext = g_slist_prepend(file_data->ext, "mar");
file_data->ext = g_slist_prepend(file_data->ext, "mvn");
file_data->ext = g_slist_prepend(file_data->ext, "mvn-r");
file_data->write_file = write_marvin;
file_data->read_file = read_marvin;
file_data->read_frame = NULL;                   /* frame reading */
sysenv.file_list = g_slist_prepend(sysenv.file_list, file_data);

/* supported file type */
file_data = g_malloc(sizeof(struct file_pak));
file_data->id = MVNOUT;
file_data->group = MARVIN;
file_data->menu = FALSE;
file_data->label = g_strdup("MARVIN output");
file_data->ext = NULL;
file_data->ext = g_slist_prepend(file_data->ext, "mot");
file_data->ext = g_slist_prepend(file_data->ext, "mvout");
file_data->ext = g_slist_prepend(file_data->ext, "mvnout");
file_data->write_file = NULL;
file_data->read_file = read_mvnout;
file_data->read_frame = NULL;                   /* frame reading */
sysenv.file_list = g_slist_prepend(sysenv.file_list, file_data);

/* supported file type */
file_data = g_malloc(sizeof(struct file_pak));
file_data->id = META_DATA;
file_data->group = META_DATA;
file_data->menu = FALSE;
file_data->label = g_strdup("Meta data");
file_data->ext = NULL;
file_data->ext = g_slist_prepend(file_data->ext, "meta");
file_data->write_file = write_meta;
file_data->read_file = NULL;
file_data->read_frame = NULL;
sysenv.file_list = g_slist_prepend(sysenv.file_list, file_data);
  
/* supported file type */
file_data = g_malloc(sizeof(struct file_pak));
file_data->id = MOL2;
file_data->group = MOL2;
file_data->menu = TRUE;
file_data->label = g_strdup("Tripos mol2");
file_data->ext = NULL;
file_data->ext = g_slist_append(file_data->ext, "mol2");
file_data->write_file = NULL;
file_data->read_file = read_mol2;
file_data->read_frame = NULL; 
sysenv.file_list = g_slist_prepend(sysenv.file_list, file_data);

/* supported file type */
file_data = g_malloc(sizeof(struct file_pak));
file_data->id = NWCHEM;
file_data->group = NWCHEM;
file_data->menu = TRUE;
file_data->label = g_strdup("NWChem");
file_data->ext = NULL;
file_data->ext = g_slist_prepend(file_data->ext, "nwin");
file_data->ext = g_slist_prepend(file_data->ext, "nw");
file_data->write_file = NULL;
file_data->read_file = read_nw;
file_data->read_frame = NULL;                   /* frame reading */
sysenv.file_list = g_slist_prepend(sysenv.file_list, file_data);

/* supported file type */
file_data = g_malloc(sizeof(struct file_pak));
file_data->id = NWCHEM_OUT;
file_data->group = NWCHEM;
file_data->menu = FALSE;
file_data->label = g_strdup("NWChem output");
file_data->ext = NULL;
file_data->ext = g_slist_prepend(file_data->ext, "nwout");
file_data->ext = g_slist_prepend(file_data->ext, "nwot");
file_data->ext = g_slist_prepend(file_data->ext, "nwo");
file_data->write_file = NULL;
file_data->read_file = read_nwout;
file_data->read_frame = read_nwout_frame;
sysenv.file_list = g_slist_prepend(sysenv.file_list, file_data);

/* supported file type */
file_data = g_malloc(sizeof(struct file_pak));
file_data->id = PDB;
file_data->group = PDB;
file_data->menu = TRUE;
file_data->label = g_strdup("PDB");
file_data->ext = NULL;
file_data->ext = g_slist_append(file_data->ext, "pdb");
file_data->write_file = write_pdb;
file_data->read_file = read_pdb;
file_data->read_frame = read_pdb_frame;
sysenv.file_list = g_slist_prepend(sysenv.file_list, file_data);

/* supported file type */
file_data = g_malloc(sizeof(struct file_pak));
file_data->id = PICTURE;
file_data->group = PICTURE;
file_data->menu = FALSE;
/* TODO - check supported image types with glib - add to string */
file_data->label = g_strdup("Picture (jpg)");
file_data->ext = rlist;
file_data->write_file = NULL;
file_data->read_file = NULL;
file_data->read_frame = NULL;                   /* frame reading */
sysenv.file_list = g_slist_prepend(sysenv.file_list, file_data);

/* supported file type */
file_data = g_malloc(sizeof(struct file_pak));
file_data->id = POVRAY;
file_data->group = POVRAY;
file_data->menu = FALSE;
file_data->label = g_strdup("POVRay");
file_data->ext = NULL;
file_data->ext = g_slist_append(file_data->ext, "pov");
file_data->write_file = write_povray;
file_data->read_file = NULL;
file_data->read_frame = NULL;                   /* frame reading */
sysenv.file_list = g_slist_prepend(sysenv.file_list, file_data);

/* supported file type */
file_data = g_malloc(sizeof(struct file_pak));
file_data->id = PROJECT;
file_data->group = PROJECT;
file_data->menu = FALSE;
file_data->label = g_strdup("Project");
file_data->ext = NULL;
file_data->ext = g_slist_append(file_data->ext, "pcf");
file_data->write_file = NULL;
file_data->read_file = project_read;
file_data->read_frame = NULL;
sysenv.file_list = g_slist_prepend(sysenv.file_list, file_data);

/* supported file type */
file_data = g_malloc(sizeof(struct file_pak));
file_data->id = CEL;
file_data->group = CEL;
file_data->menu = TRUE;
file_data->label = g_strdup("PowderCell");
file_data->ext = NULL;
file_data->ext = g_slist_append(file_data->ext, "cel");
file_data->write_file = NULL;
file_data->read_file = read_cel;
file_data->read_frame = NULL;
sysenv.file_list = g_slist_prepend(sysenv.file_list, file_data);

/* supported file type */
file_data = g_malloc(sizeof(struct file_pak));
file_data->id = RIETICA;
file_data->group = RIETICA;
file_data->menu = TRUE;
file_data->label = g_strdup("Rietica");
file_data->ext = NULL;
file_data->ext = g_slist_prepend(file_data->ext, "inp");
file_data->write_file = NULL;
file_data->read_file = read_rietica;
file_data->read_frame = NULL;                   /* frame reading */
sysenv.file_list = g_slist_prepend(sysenv.file_list, file_data);

/* supported file type */
file_data = g_malloc(sizeof(struct file_pak));
file_data->id = FDF;
file_data->group = FDF;
file_data->menu = TRUE;
file_data->label = g_strdup("SIESTA");
file_data->ext = NULL;
file_data->ext = g_slist_prepend(file_data->ext, "fdf");
file_data->write_file = write_fdf;
file_data->read_file = read_fdf;
file_data->read_frame = NULL;
sysenv.file_list = g_slist_prepend(sysenv.file_list, file_data);

/* supported file type */
file_data = g_malloc(sizeof(struct file_pak));
file_data->id = SIESTA_OUT;
file_data->group = FDF;
file_data->menu = FALSE;
file_data->label = g_strdup("SIESTA output");
file_data->ext = NULL;
file_data->ext = g_slist_prepend(file_data->ext, "sout");
file_data->ext = g_slist_prepend(file_data->ext, "sot");
file_data->write_file = NULL;
file_data->read_file = read_sout;
file_data->read_frame = read_sout_frame;
sysenv.file_list = g_slist_prepend(sysenv.file_list, file_data);

/* supported file type */
file_data = g_malloc(sizeof(struct file_pak));
file_data->id = VASP;
file_data->group = VASP;
file_data->menu = TRUE;
file_data->label = g_strdup("VASP xml");
file_data->ext = NULL;
file_data->ext = g_slist_prepend(file_data->ext, "xml");
file_data->write_file = NULL;
file_data->read_file = read_xml_vasp;
file_data->read_frame = read_xml_vasp_frame;
sysenv.file_list = g_slist_prepend(sysenv.file_list, file_data);

/* supported file type */
file_data = g_malloc(sizeof(struct file_pak));
file_data->id = USPEX;
file_data->group = USPEX;
file_data->menu = TRUE;
file_data->label = g_strdup("USPEX OUTPUT");
file_data->ext = NULL;
file_data->ext = g_slist_prepend(file_data->ext, "txt");
file_data->write_file = NULL;
file_data->read_file = read_output_uspex;
file_data->read_frame = read_frame_uspex;
sysenv.file_list = g_slist_prepend(sysenv.file_list, file_data);

/* supported file type */
file_data = g_malloc(sizeof(struct file_pak));
file_data->id = TEXT;
file_data->group = TEXT;
file_data->menu = FALSE;
file_data->label = g_strdup("Text file");
file_data->ext = NULL;
file_data->ext = g_slist_prepend(file_data->ext, "txt");
file_data->write_file = NULL;
file_data->read_file = NULL;
file_data->read_frame = NULL;
sysenv.file_list = g_slist_prepend(sysenv.file_list, file_data);

/* supported file type */
file_data = g_malloc(sizeof(struct file_pak));
file_data->id = XML;
file_data->group = XML;
file_data->menu = TRUE;
file_data->label = g_strdup("XML");
file_data->ext = NULL;
file_data->ext = g_slist_prepend(file_data->ext, "xml");
file_data->write_file = write_xml;
file_data->read_file = read_xml;
file_data->read_frame = NULL;                   /* frame reading */
sysenv.file_list = g_slist_prepend(sysenv.file_list, file_data);

/* supported file type */
file_data = g_malloc(sizeof(struct file_pak));
file_data->id = XTL;
file_data->group = XTL;
file_data->menu = TRUE;
file_data->label = g_strdup("XTL");
file_data->ext = NULL;
file_data->ext = g_slist_append(file_data->ext, "xtl");
file_data->write_file = write_xtl;
file_data->read_file = read_xtl;
file_data->read_frame = NULL;                   /* frame reading */
sysenv.file_list = g_slist_prepend(sysenv.file_list, file_data);

/* supported file type */
file_data = g_malloc(sizeof(struct file_pak));
file_data->id = XYZ;
file_data->group = XYZ;
file_data->menu = TRUE;
file_data->label = g_strdup("XYZ");
file_data->ext = NULL;
/*
file_data->ext = g_slist_append(file_data->ext, "ani");
file_data->ext = g_slist_append(file_data->ext, "xmol");
*/
file_data->ext = g_slist_append(file_data->ext, "xyz");
file_data->write_file = write_xyz;
file_data->read_file = read_xyz;
file_data->read_frame = read_xyz_frame;
sysenv.file_list = g_slist_prepend(sysenv.file_list, file_data);

sysenv.file_list = g_slist_reverse(sysenv.file_list);
}

/******************************/
/* free file recognition list */
/******************************/
void file_free(void)
{
GSList *list;

/* free data embedded in the file structure */
for (list=sysenv.file_list ; list ; list=g_slist_next(list))
  {
  struct file_pak *file_data = list->data;

  g_slist_free(file_data->ext);
  g_free(file_data->label);
  }

/* free the list and the file structure pointers */
free_slist(sysenv.file_list);
}

/************************/
/* current version info */
/************************/
void gdis_blurb(FILE *fp)
{
fprintf(fp, "Created by GDIS version %4.2f.%d\n", VERSION, PATCH);
}

/*********************************************/
/* alphabetically sort the directory listing */
/*********************************************/
gint alpha_slist_sort(gpointer ptr1, gpointer ptr2)
{
return(g_ascii_strcasecmp((gchar *) ptr1, (gchar *) ptr2));
}

/*********************************************/
/* an all-platform directory listing routine */
/*********************************************/
GSList *file_dir_list(const gchar *path, gint sort)
{
const gchar *name;
GDir *dir;
GSList *files=NULL;

/* ensure we can go up a directory */
files = g_slist_prepend(files, g_strdup(".."));

/* build the directory list */
dir = g_dir_open(path, 0, NULL);
name = g_dir_read_name(dir);
while (name)
  {
  if (g_ascii_strncasecmp(".", name, 1) != 0)
    files = g_slist_prepend(files, g_strdup(name));

  name = g_dir_read_name(dir);
  }
g_dir_close(dir);

if (sort)
  files = g_slist_sort(files, (gpointer) alpha_slist_sort);

return(files);
}

/*********************************************************************/
/* get a directory listing matched against a supplied (glob) pattern */
/*********************************************************************/
#define DEBUG_FILE_PATTERN 0
GSList *file_dir_list_pattern(const gchar *dir, const gchar *pattern)
{
gchar *name;
GPatternSpec *ps;
GSList *list, *files;

g_assert(dir != NULL);

files = file_dir_list(dir, FALSE);

if (pattern)
  ps = g_pattern_spec_new(pattern);
else
  return(files);

list = files;
while (list)
  {
  name = list->data;
  list = g_slist_next(list);

  if (!g_pattern_match_string(ps, name))
    {
    files = g_slist_remove(files, name);
    g_free(name);
    }
  }

g_pattern_spec_free(ps);

#if DEBUG_FILE_PATTERN
printf("Found %d matches\n", g_slist_length(files));
#endif

return(files);
}

/*******************************************************/
/* allocate a new string containing the file extension */
/*******************************************************/
gchar *file_extension_get(const gchar *name)
{
gint i,n;
gchar *ext;

if (!name)
  return(NULL);

/* search for '.' character in reverse order */
/*i =*/ n = strlen(name);/*FIX 4ff5be*/
/* minimum size, avoids troublesome cases (eg "..") */
if (n > 2)
  {
  for (i=n-1 ; i-- ; )
    {
    if (name[i] == '.')
      {
      ext = g_strndup(name+i+1, n-i-1);
/* ignore any .# extension */
      if (!str_is_float(ext))
        return(ext);
      g_free(ext);
/* mark the .# as the new end */
      n = i;
      }
    }
  }
return(NULL);
}

/************************************************/
/* routine to determine if a file is recognized */
/************************************************/
/* returns pointer to file info if found, NULL otherwise */
#define DEBUG_GET_FILE_INFO 0
struct file_pak *get_file_info(gpointer ptr, gint type)
{
gint code=-1;
gchar *text=NULL, *ext;
GSList *file, *ext_list;
struct file_pak *file_data;

/* checks */
g_return_val_if_fail(ptr != NULL, NULL);

/* init for search */
switch(type)
  {
  case BY_LABEL:
    text = g_strdup(ptr);

#if DEBUG_GET_FILE_INFO
printf("Searching for type [%s]\n", text);
#endif
    break;

  case BY_EXTENSION:
/* fix to get around the quotes needed to properly pass win32 filesnames with spaces */
#if __WIN32
    {
    gchar *tmp = g_shell_unquote(ptr, NULL);
    text = file_extension_get(tmp);
    g_free(tmp); 
    }
#else
    text = file_extension_get(ptr);
#endif

    if (!text)
      return(NULL);

#if DEBUG_GET_FILE_INFO
printf("Searching for extension [%s]\n", text);
#endif
    break;

  case BY_FILE_ID:
    code = GPOINTER_TO_INT(ptr);
#if DEBUG_GET_FILE_INFO
printf("Searching for code [%d]\n", code);
#endif
    break;
  }

/* search */
for (file=sysenv.file_list ; file ; file=g_slist_next(file))
  {
  file_data = file->data;

  switch (type)
    {
/* compare to all extensions in list */
    case BY_EXTENSION: 
/* go through all extensions listed under this file type */
      for (ext_list=file_data->ext ; ext_list ; ext_list=g_slist_next(ext_list))
        {
        ext = ext_list->data;
        if (strlen(text) == strlen(ext))
          {
          if (g_ascii_strcasecmp(text, ext) == 0)
            {
#if DEBUG_GET_FILE_INFO
printf("Matched: %s\n", file_data->label);
#endif
            g_free(text);
            return(file_data);
            }
          }
        }
      break;

/* compare with label */
    case BY_LABEL:
      if (strlen(text) != strlen(file_data->label))
        break;
      if (g_ascii_strcasecmp(text, file_data->label) == 0)
        {
#if DEBUG_GET_FILE_INFO
printf("Matched: %s\n", file_data->label);
#endif
        g_free(text);
        return(file_data);
        }
      break;

    case BY_FILE_ID:
      if (code == file_data->id)
        {
#if DEBUG_GET_FILE_INFO
printf("Matched: %s\n", file_data->label);
#endif
        return(file_data);
        }
      break;


    default:
      printf("get_file_info() error: bad search type.\n");
    }
  }

if (text)
  g_free(text);

return(NULL);
}

/**************************/
/* detect valid filetypes */
/**************************/
gint file_extension_valid(const gchar *name)
{
gchar *ext;
GSList *item, *list;
struct file_pak *file;

/* locate the extension */
 /*FIXME: can find_char be replaced by standard strchr()/strrchr() functions?
   - MW */
ext = find_char(name, '.', LAST);
if (ext)
  ext++;

if (sysenv.file_type == DATA) 
  /* compare against all extension types - any match is allowed */
  {
  if (!ext) /* files without extensions are matched only in strict case */
    return(FALSE);  
  for (list=sysenv.file_list ; list ; list=g_slist_next(list))
    {
    file = list->data; 
    for (item=file->ext ; item ; item=g_slist_next(item))
      if (g_ascii_strcasecmp(item->data, ext) == 0)
        return(TRUE);
    }
  }
else
  /* strict case - file must match the dialog filter */
  {
  for (list=sysenv.file_list ; list ; list=g_slist_next(list))
    {
    file = list->data; 
    if (sysenv.file_type == file->group)
        for (item=file->ext ; item ; item=g_slist_next(item))
            if ((ext && g_ascii_strcasecmp(item->data, ext) == 0)
                || (!ext && strlen(item->data) == 0))
              return(TRUE);
    }
  }

return(FALSE);
}

/*******************************************************/
/* get an unused BASENAME (of supplied extension type) */
/*******************************************************/
/* TODO - supply a basename+ext & this inserts _? until new? */
#define DEBUG_GUN 0
gchar *gun(const gchar *ext)
{
gint i;
gchar *name;
GString *filename;
FILE *fp;

/* seek a file name that doesn't exist (avoid background overwriting) */
filename = g_string_new(NULL);
i=0;
do
  {
  g_string_printf(filename,"dummy_%d.%s", i, ext);//deprecated g_string_sprintf
#if DEBUG_GUN
printf("testing: %s\n", filename->str);
#endif
  i++;
  }
while (g_file_test(filename->str, G_FILE_TEST_EXISTS));

/* create the file to prevent another process from taking it */
/* b4 the current caller gets around to writing anything to it */
fp = fopen(filename->str, "wt");
if (fp)
  {
  fprintf(fp, "locked.\n");
  fclose(fp);
  }
else
  {
  printf("Fatal error in gun()\n");
  return(NULL);
  }

name = g_string_free(filename, FALSE);

#if DEBUG_GUN
printf("using base: %s\n", name);
#endif

return(name);
}

/**************************************************************/
/* correct numbers in binary files with reverse byte ordering */
/**************************************************************/
void swap_bytes(void *ptr, const gint size)
{
gint i,j;
gchar tmp;

/*
printf("start: ");
for (i=0 ; i<size ; i++)
  printf("%d ", *((char *)(ptr+i)));
printf("\n");
*/
j=size-1;
for (i=0 ; i<size/2 ; i++)
  {
  tmp = *((gchar *)(ptr+j));
  *((gchar *)(ptr+j)) = *((gchar *)(ptr+i));
  *((gchar *)(ptr+i)) = tmp;
  j--;
  }
/*
printf(" stop: ");
for (i=0 ; i<size ; i++)
  printf("%d ", *((char *)(ptr+i)));
printf("\n");
*/
}

/**************************************/
/* get a non trivial line from a file */
/**************************************/
gint fgetline(FILE *fp, gchar *line)
{
gint i, linlen;

for(;;)
  {
/* get a line */
  if (fgets(line, LINELEN/2, fp) == NULL)
    return(1);
  linlen = strlen(line);
/* only treated as data if not a comment and non empty */
  if (line[0] != '#' && linlen)
    {
/* ampersand concatenation */
/* TODO - extra var in fgetline() (eg mode) that requests this */
/* TODO - handle multiple line concatenation */
    for (i=linlen ; i-- ; )
      {
      if (line[i] == '&')
        {
        if (fgets(&line[linlen], LINELEN/2, fp) == NULL)
          return(1);
        break;
        }
      }
    break;
    }
  }

/* all clear */
return(0);
}
/****************************/
/* reach a string in a file */
/****************************/
long int fetch_in_file(FILE *vf,const gchar *target){
	gchar *line;
	line = file_read_line(vf);
	while (line){
		if (find_in_string(target,line) != NULL) break;
		g_free(line);
		line = file_read_line(vf);
	}
	if(feof(vf)) return 0;
	return ftell(vf);
}
/*****************************************/
/* read in and return a line of any size */
/*****************************************/
gint skip_comment = TRUE;
void file_skip_comment(gint skip)
{
skip_comment = skip;
}
/* returns NULL on EOF */
gchar *file_read_line(FILE *fp)
{
gboolean skip;
gchar c;
GString *buff;

/* checks */
g_assert(fp != NULL);
c = fgetc(fp);
if (feof(fp))
  return(NULL);

/* read single chars into an expandable buffer */
buff = g_string_new(NULL);
while (!feof(fp))
  {
  skip = FALSE;
  if (c == '#' && skip_comment)
    skip = TRUE;

/* read a complete line */
  while (!feof(fp))
    {
    g_string_append_c(buff, c);
    if (c == '\n')
      break;
    if (c == '\r')
      {
      gchar d;

      d = fgetc(fp);
      if (d == '\n')
        {
/*
        printf("skipping MSDOS crap.\n");
*/
        }
      else
        ungetc(d, fp);

      break;
      }
    c = fgetc(fp);
    }
/* ignore? (keep reading) */
  if (skip)
    {
    g_string_assign(buff, "");
    c = fgetc(fp);
    }
  else
    break;
  }

/* create a properly terminated line of text */
g_string_append_c(buff, '\0');

/* free the GString, but not the text */
return(g_string_free(buff, FALSE));
}

/**********************************************************/
/* format value as a float printed to dp places plus unit */
/**********************************************************/
gchar *format_value_and_units(gchar *string, gint dp)
{
gint num_tokens;
#ifdef _BUG_
GString *text, *format_string;
#else
gchar *format_str;
#endif //_BUG_ (OVHPA_1)
gchar *text_str, *units, **buff;
gdouble float_val;

float_val = str_to_float(string);
buff = tokenize(string, &num_tokens);
if (num_tokens < 2)
  units = "";
else
  units = *(buff+1);
#ifdef _BUG_
format_string = g_string_new("");
g_string_append_printf(format_string, "%%.%df %%s", dp);
text = g_string_new("");
g_string_append_printf(text, format_string->str, float_val, units);
text_str = text->str;
g_string_free(format_string, TRUE); 
g_string_free(text, FALSE);
g_strfreev(buff);
#else
format_str=malloc(32*sizeof(char));
text_str=calloc(MAX_VALUE_SIZE,sizeof(char));
sprintf(format_str,"%%.%df %%s", dp);
sprintf(text_str,format_str,float_val,units);
free(format_str);
g_strfreev(buff);
#endif //_BUG_ (OVHPA_1)
return (text_str);
}

/************************************************/
/* read in a raw frame (animation and analysis) */
/************************************************/
#define DEBUG_READ_RAW 0
gint read_raw_frame(FILE *fp, gint n, struct model_pak *model)
{
gint flag=99;
gchar *filename;
fpos_t *offset;
FILE *fp2;
struct file_pak *file;

g_assert(fp != NULL);
g_assert(model != NULL);

/*in case model->id is invalid (-1) file is set to NULL*/
file = get_file_info(GINT_TO_POINTER(model->id), BY_FILE_ID);

#if DEBUG_READ_RAW
printf("reading frame: %d\n", n);
#endif

/* position the file pointer */
if (model->id != GULP)
  {
  offset = g_list_nth_data(model->frame_list, n);

#if DEBUG_READ_RAW
printf("offset = %d\n", offset);
#endif

  if (!offset)
    return(1);
  if (fsetpos(fp, offset))
    {
    printf("Error positioning file pointer.\n");
    return(2);
    }
/* use supplied routine (if available) */
  if(file==NULL) {
#if DEBUG_READ_RAW
fprintf(stdout,"WRONG model id -> can't process frame!\n");
#endif
	return 5;
}
  if (file->read_frame) 
    flag = file->read_frame(fp, model);
  else
    gui_text_show(ERROR, "No animation read routine.\n");
  }
else
  {
/* GULP exception */
  if (!n)
    {
/* NEW - 0th frame exception - load original GULP coords */
    if (model->cores)
      g_assert(model->gulp.orig_cores != NULL);
    if (model->shels)
      g_assert(model->gulp.orig_shells != NULL);

/* free connectivity lists */
    free_slist(model->bonds);
    free_slist(model->moles);
    model->bonds = NULL;
    model->moles = NULL;

/* restore coords */
    free_slist(model->cores);
    free_slist(model->shels);
    model->cores = dup_core_list(model->gulp.orig_cores);
    model->shels = dup_shell_list(model->gulp.orig_shells);

/* restore lattice */
/*
      model->fractional = model->gulp.orig_fractional;
*/
/* FIXME - I don't know why this works (and the above doesn't)... */
    model->fractional = TRUE;
    model->construct_pbc = model->gulp.orig_construct_pbc;
    memcpy(model->latmat, model->gulp.orig_latmat, 9*sizeof(gdouble));
    memcpy(model->pbc, model->gulp.orig_pbc, 6*sizeof(gdouble));

/* turn space group lookup on */
    model->sginfo.lookup = TRUE;
 
    flag = 0;
    }
  else
    {
/* turn off space group lookup */
/* trg files contain ALL atoms & symmetry may be broken */
    model->sginfo.lookup = FALSE;

    filename = g_strdup_printf("%s/%s", sysenv.cwd, model->gulp.trj_file);
/* NB: GULP binary trajectory files - MUST specify "rb" */
/* FIXME - some GULP traj files may be ascii I think? */
    fp2 = fopen(filename, "rb");
    if (!fp2)
      {
      printf("Failed to open: %s\n", filename);
      g_free(filename);
      return(3);
      }

/* FIXME - ugly hack to get the nth frame */
#define CRAP 0
#if CRAP
    read_trj_header(fp2, model);
    for (i=0 ; i<n-1 ; i++)
      read_trj_frame(fp2, model, FALSE);
    read_trj_frame(fp2, model, TRUE);
    fclose(fp2);
#else

/* CURRENT - implementing fgetpos/fsetpos method */
  offset = g_list_nth_data(model->frame_list, n);
  if (fsetpos(fp2, offset))
    {
    printf("Error positioning file pointer.\n");
    g_free(filename);
    return(2);
    }

  read_trj_frame(fp2, model, TRUE);

#endif



    g_free(filename);
    flag = 0;
    }
  }
return(flag);
}

/***********************************/
/* handler for normal file loading */
/***********************************/
#define DEBUG_FILE_LOAD 0
void file_load(gchar *filename, struct model_pak *mdata)
{
gint j;
gint flag, status=-1;
gchar *fullname;
GSList *list;
struct model_pak *data;
struct file_pak *file_data;

#if DEBUG_FILE_LOAD
printf("loading: [%s] into: %p\n", filename, mdata);
#endif

/* FIXME - new error reporting - doesnt work for gulp */
error_table_clear();
error_table_enable();

/* get the file structure required to parse the file */
switch (sysenv.file_type)
  {
/* cases for which duplicate extensions occur */
  case RIETICA:
  case GAMESS:
  case DLPOLY:
    file_data = get_file_info(GINT_TO_POINTER(sysenv.file_type), BY_FILE_ID);
    break;

/* cases for which extensions must be used to determine type (eg .gin/.got) */
  default:
    file_data = get_file_info((gpointer *) filename, BY_EXTENSION);
  }
if (!file_data)
  return;

#if DEBUG_FILE_LOAD
printf("Using file load routine for: [%s] files.\n", file_data->label);
#endif

if (mdata)
  data = mdata;
else
  {
/* get the new model number */
  data = model_new();
  if (!data)
    {
    gui_text_show(ERROR, "Model memory allocation failed.\n");
    return;
    }
  }

data->id = file_data->id;

/* read file if exists, else try prepending current working directory */
if (file_data->read_file)
  {
  if (g_path_is_absolute(filename))
    {
    status = file_data->read_file(filename, data);
    }
  else
    {
    fullname = g_build_filename(sysenv.cwd, filename, NULL);
    status = file_data->read_file(fullname, data);
    g_free(fullname);
    }
  }
else
  gui_text_show(ERROR, "No read routine for this type. ");

/* NEW - display error (if any) */
if ((data->error_file_read)->len)
  {
  gui_text_show(ERROR, (data->error_file_read)->str);
  }

/* check for successful file load */
if (status)
  {
  gui_text_show(ERROR, "Load failed.\n");

printf("Load failed, error code: %d\n", status);

  model_delete(data);
  }
else
  {
/* we don't know how many new models were loaded so */
/* scan through them all & check for initialization */
  for (list=sysenv.mal ; list ; list=g_slist_next(list))
    {
    data = list->data;

/* skip if already on the tree */
    if (data->grafted)
      continue;


/* NEW - big model display mode exceptions */
/* CURRENT - this doesn't work - something not set up? */
/* TODO - automaticaly choose an appropriate zone size */
/*
if (data->num_atoms > 100000)
  {
  gpointer za;

  core_render_mode_set(ZONE, data->cores);

  za = zone_make(10.0, data);
  zone_display_init(za, data);
  zone_free(za);
  }
else
*/
  {
  if (data->num_atoms > 1000)
    core_render_mode_set(STICK, data->cores);
  }

/* surfaces are always conv by default */
    if (data->periodic == 2)
      data->gulp.method = CONV;
/* not on tree - must have just been loaded */
    tree_model_add(data);

/* create gulp supercells */
    flag=0;
    for (j=0 ; j<3 ; j++)
      {
      if (data->gulp.super[j] > 1)
        {
        data->image_limit[2*j+1] = data->gulp.super[j];
        data->expected_cores *= data->gulp.super[j];
        data->expected_shells *= data->gulp.super[j];
        flag++;
        }
      }
    if (flag)
      {
      space_make_images(CREATE, data);
      space_make_supercell(data);
      model_prep(data);
      }

/* NEW - store initial frame, as trajectory animation will overwrite */
if (data->gulp.trj_file && !data->gulp.orig_cores)
  {
/* FIXME - saving fractional type here instead of at the end of read_gulp() stuffs things up  - why? */
/* NB: we need to save the coords here (and not earlier) in case a supercell was created */

  memcpy(data->gulp.orig_pbc, data->pbc, 6*sizeof(gdouble));
  memcpy(data->gulp.orig_latmat, data->latmat, 9*sizeof(gdouble));
  data->gulp.orig_cores = dup_core_list(data->cores);
  data->gulp.orig_shells = dup_shell_list(data->shels);
  }

    sysenv.active_model = data;
    }
  }

/* finished - only now we destroy the file dialog */
dialog_destroy_type(FILE_SELECT);
tree_select_active();

error_table_print_all();
}

/****************/
/* load dialog */
/***************/
#ifdef WITH_GUI
void file_load_dialog(void)
{
/* revert to data filetypes for special cases */
switch (sysenv.file_type)
  {
  case MD_ANALYSIS:
  case PICTURE:
  case PROJECT:
  case GEOMVIEW_OFF:
    sysenv.file_type = DATA;
  }

file_dialog("Load file", NULL, FILE_LOAD, (gpointer) file_load, sysenv.file_type);
}
#endif

/**********************************/
/* save file with name */
/**********************************/
gint file_save_as(gchar *filename, struct model_pak *model)
{
gint ret;
#ifdef UNUSED_BUT_SET
gint id;
#endif
gchar *name;
struct file_pak *file_data;

/* setup & checks */
if (!model)
  model = sysenv.active_model;

if (!filename || !model)
  return(1);

file_data = get_file_info((gpointer *) filename, BY_EXTENSION);

if (file_data)
  {
/* use extension */
#ifdef UNUSED_BUT_SET
  id = file_data->id;
#endif
  strcpy(model->filename, filename);
  g_free(model->basename);
  model->basename = parse_strip(filename);
  }
else
  {
  gchar *ext;

/* no extension - attempt to use filter type */
  file_data = get_file_info(GINT_TO_POINTER(sysenv.file_type), BY_FILE_ID);
  if (!file_data)
    return(2);
  ext = g_slist_nth_data(file_data->ext, 0);
  if (!ext)
    return(2);

#ifdef UNUSED_BUT_SET
  id = file_data->id;
#endif
  g_free(model->basename);
  model->basename = parse_strip(filename);

  name = g_strdup_printf("%s.%s", model->basename, ext);
  strcpy(model->filename, name);
  g_free(name);
  }

/* file info indicates routine to call */
if (file_data->write_file == NULL)
  {
  printf("No write routine for this type.\n");
  return(3);
  }
else
  {
/* build the full name */
  name = g_build_filename(sysenv.cwd, model->filename, NULL);
  ret = file_data->write_file(name, model);
  g_free(name);
  }

/* update */
if (ret)
  gui_text_show(ERROR, "Save failed.\n");
else
  {
  gui_text_show(STANDARD, "Saved model.\n");
  dialog_destroy_type(FILE_SELECT);
  tree_model_refresh(model);
  redraw_canvas(SINGLE);
  }

return(ret);
}

/*************************************/
/* save active model using same name */
/*************************************/
void file_save(void)
{
gchar *name;
struct model_pak *model;

model = sysenv.active_model;
if (!model)
  return;

/* strip the path, as file_save_as() prepends the cwd */
name = g_path_get_basename(model->filename);
file_save_as(name, model);
g_free(name);
}

/******************/
/* save as dialog */
/******************/
#ifdef WITH_GUI
void file_save_dialog(void)
{
gchar *text=NULL;
struct model_pak *model;

model = sysenv.active_model;
if (model)
  text = model->basename;

file_dialog("File save", text, FILE_SAVE, (gpointer) file_save_as,
                                                sysenv.file_type);
}
#endif

/*************************************/
/* get the size of a file (in bytes) */
/*************************************/
gint file_byte_size(const gchar *filename)
{
struct stat buff;

stat(filename, &buff);

return(buff.st_size);
}

/************************************/
/* search for an executable program */
/************************************/
/* TODO - put elsewhere ... look for a program ... ask for a path (can cancel) if not found */
gchar *file_find_program(const gchar *name)
{
gchar *path;

path = g_find_program_in_path(name);

return(path);
}

/***********************************/
/* change to a specified directory */
/***********************************/
#define DEBUG_SET_PATH 0
gint set_path(const gchar *txt)
{
gint i, n, status=0;
gchar **buff, *text;
GString *path;

/* split by directory separator */
text = g_strdup(txt);
g_strstrip(text);
buff = g_strsplit(text, DIR_SEP, 0);
path = g_string_new(NULL);

/* find the number of tokens */
n=0;
while(*(buff+n))
  {
#if DEBUG_SET_PATH
printf("[%s] ", *(buff+n));
#endif
  n++;
  }
#if DEBUG_SET_PATH
printf(": found %d tokens.\n", n);
#endif

/* truncate token list if parent directory selected */
if (n > 1)
  if (g_ascii_strncasecmp("..",*(buff+n-1),2) == 0)
    n -= 2;

/* build new path */
/* this is a bit fiddly, as we don't want a dir_sep on the end */
/* EXCEPT if it's the root directory (blame windows for this) */
g_string_printf(path, "%s", *buff);//deprecated g_string_sprintf
i=1;
while (i < n)
  {
  g_string_append_printf(path, "%s%s", DIR_SEP, *(buff+i));//deprecated g_string_sprintfa
  i++;
  }
if (n < 2)
  g_string_append_printf(path, "%s", DIR_SEP);//deprecated g_string_sprintfa

#if DEBUG_SET_PATH
printf("testing path [%s] ... \n", path->str); 
#endif

if (g_file_test(path->str, G_FILE_TEST_IS_DIR))
  {
  g_free(sysenv.cwd);
  sysenv.cwd = g_strdup(path->str);
  }
else
  {
#if DEBUG_SET_PATH
printf(" Not a directory.\n");
#endif
  status++;
  }

g_free(text);
g_strfreev(buff);
g_string_free(path, TRUE);

return(status);
}
/************************************/
/* copy a text file f_src to f_dest */
/************************************/
#define DEBUG_DUMB_COPY 0
gboolean dumb_file_copy(gchar *f_src,gchar *f_dest){
/*This is a very dumb version which reads each line of a file src, and copy in into destination */
/*A better and portable way to do that would be to use the g_file_copy from GIO library --OVHPA */
gchar *buffer;
//GError *error=NULL;
gsize length;
#if DEBUG_DUMB_COPY
fprintf(stdout,"copy %s into %s ... ",f_src,f_dest);
#endif
/* error can be set even though g_file_get_contents succeed!
if(g_file_get_contents (f_src,&buffer,&length,&error)){
if(error!=NULL) {//who would set error on a success?
#if DEBUG_DUMB_COPY
	fprintf(stdout,"error#1=%s\n",error->message);
#endif
	g_error_free (error);
	g_free(error);
}
*/
if(g_file_get_contents (f_src,&buffer,&length,NULL)){
//	if(g_file_set_contents (f_dest,buffer,(gssize) length,&error)){
	if(g_file_set_contents (f_dest,buffer,(gssize) length,NULL)){
#if DEBUG_DUMB_COPY
fprintf(stdout,"SUCCESS!\n");
#endif
		return TRUE;
	}
}
#if DEBUG_DUMB_COPY
fprintf(stdout,"FAILED!\n");
#endif
//fprintf(stderr,"Error during copy of %s into %s\nError code %i: %s",f_src,f_dest,error->code, error->message);
fprintf(stderr,"Error during copy of %s into %s\n",f_src,f_dest);
return FALSE;
}
/*************************************************************/
/* copy all files from a src directory to a target directory */
/* src and target must i) exists and ii) be readable, target */
/* should be iii) writable too. Uses dumb_file_copy. --OVHPA */
/*************************************************************/
gboolean dumb_dir_copy(gchar *src, gchar *dest){
gboolean is_ok=FALSE;
GDir *d_src;
const gchar *current=NULL;
gchar *f_src;
gchar *f_dest;
d_src=g_dir_open(src,0,NULL);
current=g_dir_read_name(d_src);
while(current){
	f_src=g_build_filename(src,current,NULL);
	f_dest=g_build_filename(dest,current,NULL);
	if(g_file_test(f_src,G_FILE_TEST_IS_DIR)){
		if(!g_mkdir(f_dest,0775)) {
			fprintf(stderr,"COPY aborted: Can't create directory %s\n",f_dest);
			g_free(f_src);
			g_free(f_dest);
			return FALSE;
		}
		is_ok=dumb_dir_copy(f_src,f_dest);/*recursive*/
	}else{
		is_ok=dumb_file_copy(f_src,f_dest);
	}
	if(is_ok==FALSE) {
		fprintf(stdout,"Interrupted copy!\n");
		g_free(f_src);
		g_free(f_dest);
		return FALSE;
	}
	g_free(f_src);
	g_free(f_dest);
	current=g_dir_read_name(d_src);
}
return TRUE;
}




