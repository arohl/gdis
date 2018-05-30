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

/* header for USPEX file */

typedef enum {/*released USPEX calculation methods as of 2018 (5)*/
	US_CM_USPEX,	/*USPEX*/
	US_CM_META,	/*metadynamics*/
	US_CM_VCNEB,	/*VC-NEB*/
	US_CM_PSO,	/*PSO*/
} uspex_method;

typedef struct {
        gint gen;		/*generation of this individual*/
	gint natoms;		/*number of atoms*/
	gint *atoms;		/*number of atoms per species*/
        gdouble energy;		/*total energy*/
	gdouble E;		/*reduce energy / atoms*/
        gdouble volume;		/*total volume*/
} uspex_individual;

/* global calculation structure */
typedef struct {
/*name*/
	gchar *name;
/*system parameters*/
	uspex_method method;
	gint type;
	/* in detail */
        gint dim;
        gboolean mol;
        gboolean var;
	/* optimization */
	gint opt_type;

/*structures*/
	gint num_gen;
	gint num_struct;
	gint nspecies;          /*number of species*/
	gdouble min_E;
	gdouble max_E;
	uspex_individual *ind;
	gint num_best;
	gint *best_ind;
/*interpretation*/
	gpointer graph;
	gpointer graph_best;
	gpointer graph_comp;
} uspex_calc_struct;
/* execution structure */
typedef struct {
	int job_id;
	gboolean have_result;
	/*job related*/
	gchar *job_uspex_exe;
	gchar *job_path;
} uspex_exec_struct;

/*methods of interest*/


