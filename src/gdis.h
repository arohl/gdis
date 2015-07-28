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
#define JJM_DEBUG
#define VERSION 0.91
#define PATCH 0
#define YEAR 2015

/********************************************/
/*** Debugging switches					 ****/
/********************************************/

/*** Molecular surface debugging ***/
#define DEBUG_MOLSURF 1

#if DEBUG_MOLSURF
#define DEBUG_CREATE_TRIANGLES 1
#define DEBUG_MS_CUBE 1
#endif

/***Hirshfeld surface debugging ***/
#define DEBUG_HIRSHFELD 1

#if DEBUG_HIRSHFELD
#define DEBUG_DE 1
#define DEBUG_SHAPE 1
#define DEBUG_HFS_NORMALS 1
#define DEBUG_HFS_CALC_ENV 1
#endif

/**** end debugging switches ***/

#define MAX_ELEMENTS 110
#define MAX_TOKENS 100
#define MAX_DISPLAYED 9
#define MAX_MODELS 14
#define ELEM_LABEL_SIZE 7
#define LABEL_SIZE 7
#define TAIL_SIZE 80
#define DEFAULT_NA 30
#define DEFAULT_NL 10
#define DEFAULT_SS 10
#define BOND_FUDGE 0.1
#define RMAX_FUDGE 1.1
/*#define ROTATE_SCALE 1.0*/
#define ROTATE_SCALE 0.05
#define SCALE_MAG 10
#define MESH_POINTS_MIN 2
#define MESH_POINTS_MAX 20
#define MESH_SPACE_MIN 0.1
#define MESH_SPACE_MAX 10.0
/* 1st is roughly (no accounting for latmat) in Angs */
#define POSITION_TOLERANCE 0.01
#define FRACTION_TOLERANCE 0.000001

/* parsing */
#define FILELEN 512
#define LINELEN  300

/* maths */
#define PI 3.1415926535897932384626433832795
#define D2R 0.01745329252
#define R2D 57.295779513
#define AU2ANG 0.529177249

/* fundamental constants, etc (C. Fisher 2004) */
#define ELCHARGE 1.60217733e-19
#define AMU 1.6605402e-27
#define AVOGADRO 6.0221367e23

/* location of data file - TODO - IMPROVE ie rc file/scan etc */
#ifdef __WIN32
#define INIT_FILE "gdis.rc"
#else
#define INIT_FILE ".gdisrc"
#endif

#define ELEM_FILE "gdis.elements"
#define LIBRARY_FILE "gdis.library"

/*
#ifdef __WIN32
#define GFONT "-adobe-helvetica-medium-r-normal--*-120-*-*-*-*-*-*"
#define DFONT "-adobe-helvetica-medium-r-normal--*-160-*-*-*-*-*-*"
#else
#endif
*/
#define GTK_FONT "Sans 10"
/*#define GL_FONT "Sans 14"*/
#define GL_FONT "Helvetica 14"

/* camera placement for POVray */
#define FRAMEZ 2003

/* non-contact molecular surface types */
#define REENTRANT -1
#define HOLE -2
#define VALLEY -3

/* convenience macros */
#define SPIN_FVAL gtk_spin_button_get_value_as_float
#define SPIN_IVAL gtk_spin_button_get_value_as_int

#define P3VEC(s,v) (printf("%s (%.8lf, %.8lf, %.8lf)\n",s,*v,*(v+1),*(v+2)))
#define P4VEC(s,v) (printf("%s (%.8lf, %.8lf, %.8lf, %.8lf)\n",s,*v,*(v+1),*(v+2),*(v+3)))

#define P3MAT(s,m) {printf("%s\n|%lf, %lf, %lf|\n|%lf, %lf, %lf|\n|%lf, %lf, %lf|\n",\
s,*m,*(m+1),*(m+2),*(m+3),*(m+4),*(m+5),*(m+6),*(m+7),*(m+8));}

#define SUBLINE(pix,gc,sub,x1,y1,x2,y2) (gdk_draw_line(pix, gc, \
x1+sysenv.subcenx[sub], y1+sysenv.subceny[sub], \
x2+sysenv.subcenx[sub], y2+sysenv.subceny[sub]))

#define SUBTEXT(pix,gc,sub,x1,y1,text) (gdk_draw_text(pix, sysenv.dfont, gc, \
x1+sysenv.subcenx[sub], y1+sysenv.subceny[sub], \
text, strlen(text)))

#define SUBARC(pix,gc,sub,x,y,w,h,a1,a2) (gdk_draw_arc(pix, gc, FALSE,\
x+sysenv.subcenx[sub], y+sysenv.subceny[sub], w, h, a1, a2))

#define FSUBARC(pix,gc,sub,x,y,w,h,a1,a2) (gdk_draw_arc(pix, gc, TRUE,\
x+sysenv.subcenx[sub], y+sysenv.subceny[sub], w, h, a1, a2))

/* linked list loops */
#define FOR_SLIST(item, list) for (item=list ; (item=g_slist_next(item)) ; )
#define FOR_LIST(item, list) for (item=list ; (item=g_list_next(item)) ; )


/* file parsing keyword/option stuff */
enum 
{
E_SINGLE, E_OPTIMIZE, DOCK, FREE_ENERGY, MD, CONP, CONV,
MOLE, MOLMEC, MOLQ, NOBUILD, NOFLAGS, 
NAME, CELL, CART, FRAC, SFRAC, PFRAC,
SURFACE_CELL, SURFACE_VECTORS, POLYMER_CELL, POLYMER_VECTOR, LATTICE_VECTORS,
SPACE, ORIGIN, SPECIES, MAXCYC, ENERGY, TOTAL_ENERGY, SBULK_ENERGY,
GULP_POTENTIAL, GULP_SINGLE_LINE_POTENTIAL, GULP_NOAUTOBOND, GULP_CONNECT,
SWITCH_OFF, SWITCH_ON, CYCLE, GNORM, BFGS_OPT, RFO_OPT, CONJ_OPT, UNIT_HESSIAN,
TEMPERATURE, PRESSURE, ENSEMBLE,
SHRINK, QEQ, ZSISA, COMPARE, NOSYM, FIX, PHONON, EIGEN,
SUPER_CELL, NVE, NVT, NPT, ENDFORCE,
TIMESTEP, EQUILIBRATION, PRODUCTION, SAMPLE, WRITE, PRINT, 

/* FIXME - separate out solvation to its own sub-structure */
GULP_SOLVATION_NONE, GULP_SOLVATION_COSMO, GULP_SOLVATION_COSMIC,
GULP_COSMO_SHAPE, GULP_COSMO_POINTS, GULP_COSMO_SEGMENTS, GULP_COSMO_SMOOTHING,
GULP_COSMO_SOLVENT_EPSILON, GULP_COSMO_SOLVENT_RADIUS, GULP_COSMO_SOLVENT_RMAX,

KPOINTS, GULP_LIBRARY, IGNORE, ERONGI,
TITLE, END, GULP_OBSERVABLES, TEMP_FILE, DUMP_FILE, LIB_FILE, OUTPUT, 
NO_ESURF, NO_EATT, D_HKL,
KEYWORD, 
ELEMENT, SYMBOL, NUMBER, LABEL, WEIGHT, COVALENT, IONIC, VDW, CHARGE, COLOUR,
GULP_VARIABLES,
LABEL_NORMAL, LABEL_GHOST, LABEL_FF,
COORD_X, COORD_Y, COORD_Z,
CORE_GROWTH_SLICE, CORE_REGION, CORE_TRANSLATE, CORE_FF,
CIF_AUDIT, CIF_CREATION_DATE, CIF_DATA_START, CIF_LOOP_START,
CIF_DATABASE_CODE, CIF_CHEMICAL_NAME, CIF_MINERAL_NAME,
CIF_CELL_A, CIF_CELL_B, CIF_CELL_C, CIF_CELL_ALPHA, CIF_CELL_BETA,
CIF_CELL_GAMMA, CIF_SPACE_NAME, CIF_SPACE_NUM, CIF_REFINE,
CIF_ATOM_SITE, CIF_LABEL, CIF_TYPE_SYMBOL,
CIF_SOF, CIF_EQUIV_SITE, CIF_EQUIV_POS,
CIF_CART_X, CIF_CART_Y, CIF_CART_Z, CIF_FRAC_X, CIF_FRAC_Y, CIF_FRAC_Z,
MARVIN_BASIS, MARVIN_COORD,
DIPOLE_TOLERANCE,
/* Moldy parameters */
CONTROL_FILE, SYSSPEC_FILE, SAVE_FILE, RESTART_FILE, PRINTINT, NDUMPS,
DELTAT, STEPS, CUTOFF, ALPHA, KCUTOFF,
SCALEINT, SCALEEND, DUMPSTART, DUMPINT, AVSTART, AVINT, ROLLINT,
RDFSTART, RDFINT, RDFOUT, RDFLIM, NBINS,
BACKINT, BACKFILE, TTMASS, RTMASS,
MASSPARM, SUBCELL, ACCURACY, NMOLS, DENSITY,
DREIDING, GASTEIGER,
/* gdis config file parsing */
GDIS_ELEM_START, GDIS_END
};

/* Moldy unit types */
enum { TIME_UNIT, LENGTH_UNIT, MASS_UNIT, CHARGE_UNIT };

/* misc specifiers */
enum {UP, DOWN};
enum {DEFAULT, MODIFIED};
enum {ALL, SINGLE};
enum {APPEND, INSERT, REPLACE};
enum {LOW, MEDIUM, HIGH};
enum {PLANE_XY, PLANE_YZ, PLANE_XZ};
enum {CHECK, CONDENSE};
enum {EXEC_PATH};
enum {SELECTED, ANY};
enum {BY_FILE_ID, BY_EXTENSION, BY_LABEL};
enum {STANDARD, WARNING, ERROR, BOLD, ITALIC};
enum {CARTESIAN, OTHER, ANGSTROM, BOHR, DEGREES, RADIANS};



/* surface types */
enum {SOLSURF, MOLSURF};

/* some colour names */
enum {BLACK, WHITE, RED, GREEN, BLUE, YELLOW, GOLD, ORANGE,
      SEA_GREEN, AQUAMARINE, NAVY_BLUE, LIGHT_TAN, LEMON, PURPLE};

/* data type specifiers */
enum {MODEL, CORE, BOND, UBOND, HBOND, MOL, SHELL, ELEM,
      ATOM_LABEL, ATOM_TYPE, ELEM_MOL, FRAGMENT, REGION, GROWTH_SLICE,
      TRANSLATE, OCCUPANCY, VELOCITY,
      GEOM, DIST, ANGLE, DIHEDRAL, SELECT,
      ZEOLITE, RIBBON};

/* coords_init() call modes */
enum 
{
INIT_COORDS, CENT_COORDS, REDO_COORDS,
MESH_ON, MESH_OFF, BOX_ON, BOX_OFF,
};

/* crystal lattice types */
enum {XS_UNKNOWN, XS_TRICLINIC, XS_MONOCLINIC, XS_ORTHOROMBIC,
      XS_TETRAGONAL, XS_TRIGONAL, XS_HEXAGONAL, XS_CUBIC};

/* unit cell facets */
enum {A_LO, A_HI, B_LO, B_HI, C_LO, C_HI};

/* general update call modes */
enum {INITIAL, REFRESH, RESTORE, CREATE};

/* atom status flags */
#define NORMAL 1
#define DELETED 2
#define HIDDEN 4
#define SELECT 8
#define PRUNED 16
#define ZEOL_HIDDEN 32
#define OFF_SCREEN 64

/* drawing mode */
enum {DRAW_OFF, DRAW_GDK, DRAW_GL};

/* creator grid */
enum {MESH_POINTS, MESH_SPACING};

/* marvin region types */
enum {REGION1A, REGION2A, REGION1B, REGION2B, NONE};

/* selection modes */
enum {CLEAN, START, UPDATE, STOP, ASSIGN, RECALL, RELEASE};

/*****************/
/* MAIN INCLUDES */
/*****************/
#include <glib.h>
#include <gtk/gtk.h>
#include "pak.h"

#if GTK_EXTRA
#include <gtkextra/gtksheet.h>
#endif

