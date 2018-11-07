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

/* simple USPEX calcul interface */
/* DEFINES: sync uspex_gui and uspex_gui.calc */
#define USPEX_REG_VAL(value,format) do{\
	GUI_REG_VAL(uspex_gui.value,uspex_gui.calc.value,format);\
}while(0)
#define USPEX_REG_TEXT(value) do{\
	if(uspex_gui.calc.value!=NULL) g_free(uspex_gui.calc.value);\
	if(GUI_ENTRY_LENGHT(uspex_gui.value)>0) GUI_ENTRY_GET_TEXT(uspex_gui.value,uspex_gui.calc.value);\
	else uspex_gui.calc.value=NULL;\
}while(0)
/*page numbers*/
#define USPEX_PAGE_SET_I 1
#define USPEX_PAGE_SET_II 2
#define USPEX_PAGE_SPECIFICS 3
#define USPEX_PAGE_ABINITIO 4
#define USPEX_PAGE_EXEC 5
/* gui structure */
struct uspex_calc_gui{
	/*window information*/
	GUI_OBJ *window;
	/*connection to calculation parameters*/
	uspex_calc_struct calc;
	/*actual GUI*/
	gint cur_page;
	gboolean is_dirty;
	GUI_OBJ *name;
	GUI_OBJ *file_entry;
	gboolean have_output;
/*4.1 Type of run & System*/
	GUI_OBJ *calculationMethod;
	GUI_OBJ *calculationType;
	GUI_OBJ *_calctype_dim;
	gdouble _dim;/*absurd spin on double*/
	GUI_OBJ *_calctype_mol;
	GUI_OBJ *_calctype_var;
	GUI_OBJ *_calctype_mag;			/*VER 10.1*/
	GUI_OBJ *optType;
	gboolean have_new_opt;			/*VER 10.1*/
	GUI_OBJ *new_optType;			/*VER 10.1*/
	gchar *_tmp_new_optType;		/*VER 10.1*/
	/*atoms definition: apply/remove*/
	GUI_OBJ *atomType;
	GUI_OBJ *_atom_sym;
	GUI_OBJ *_atom_typ;
	GUI_OBJ *_atom_num;
	GUI_OBJ *_atom_val;
	/*temporary values*/
	gchar _tmp_atom_sym[3];
	gint _tmp_atom_typ;
	gint _tmp_atom_num;
	gint _tmp_atom_val;
	/*bond definition: apply/auto*/
	GUI_OBJ *goodBonds;
	GUI_OBJ *_bond_d;
	gchar *_tmp_bond_d;
	gboolean auto_bonds;
	/**/
	//checkMolecules is auto-sync
	//checkConnectivity is auto-sync
	GUI_OBJ *fitLimit;			/*VER 10.1*/
/*4.2 Population*/
	GUI_OBJ *populationSize;
	GUI_OBJ *initialPopSize;
	GUI_OBJ *numGenerations;
	GUI_OBJ *stopCrit;
	GUI_OBJ *mag_nm;			/*VER 10.1*/
	GUI_OBJ *mag_fmls;			/*VER 10.1*/
	GUI_OBJ *mag_fmhs;			/*VER 10.1*/
	GUI_OBJ *mag_afml;			/*VER 10.1*/
	GUI_OBJ *mag_afmh;			/*VER 10.1*/
	GUI_OBJ *mag_fmlh;			/*VER 10.1*/
	GUI_OBJ *mag_aflh;			/*VER 10.1*/
/*4.3 Survival of the fittest & Selection*/
	GUI_OBJ *bestFrac;
	GUI_OBJ *keepBestHM;
	GUI_OBJ *reoptOld;
/*4.4 Structure generation and variation operators*/
	GUI_OBJ *symmetries;
	GUI_OBJ *fracGene;
	GUI_OBJ *fracRand;
	GUI_OBJ *fracTopRand;			/*VER 10.1*/
	GUI_OBJ *fracPerm;
	GUI_OBJ *fracAtomsMut;
	GUI_OBJ *fracRotMut;
	GUI_OBJ *fracLatMut;
	GUI_OBJ *fracSpinMut;			/*VER 10.1*/
	GUI_OBJ *howManySwaps;
	GUI_OBJ *specificSwaps;
	GUI_OBJ *mutationDegree;
	GUI_OBJ *mutationRate;
	GUI_OBJ *DisplaceInLatmutation;
	GUI_OBJ *AutoFrac;
/*4.5 Constrains*/
	GUI_OBJ *IonDistances;
	GUI_OBJ *_distances;
	gchar *_tmp_distances;
//	GUI_OBJ *curr_distance;
	gint _curr_distance;
	GUI_OBJ *minVectorLength;
	GUI_OBJ *MolCenters;
	GUI_OBJ *_centers;
	gchar *_tmp_centers;
//	GUI_OBJ *curr_center;
	gint _curr_center;
	GUI_OBJ *constraint_enhancement;
	gboolean auto_C_ion;
	gboolean auto_C_lat;
	GUI_OBJ *_molModels;
/*4.6 Cell*/
	GUI_OBJ *Latticevalues;
	GUI_OBJ *_latticeformat;
	GUI_OBJ *_latticevalue;
	gchar *_tmp_latticevalue;
	GUI_OBJ *splitInto;
	gchar *_tmp_splitInto;
	gboolean auto_lval;
/*4.9 fingerprint*/
	GUI_OBJ *RmaxFing;
	GUI_OBJ *deltaFing;
	GUI_OBJ *sigmaFing;
/*4.10 Antiseed*/
	GUI_OBJ *antiSeedsActivation;
	GUI_OBJ *antiSeedsMax;
	GUI_OBJ *antiSeedsSigma;
/*4.11 spacegroup*/
	GUI_OBJ *doSpaceGroup;
	GUI_OBJ *SymTolerance;


/*5.5 variable composition*/
	GUI_OBJ *firstGeneMax;
	GUI_OBJ *minAt;
	GUI_OBJ *maxAt;
	GUI_OBJ *fracTrans;
	GUI_OBJ *howManyTrans;
	gchar *_tmp_specificTrans;
	GUI_OBJ *specificTrans;

/* CALCUL */
	GUI_OBJ *job_uspex_exe;
	GUI_OBJ *job_path;
	gboolean have_result;
	gint index;
/*buttons*/
	GUI_OBJ *button_save;
	GUI_OBJ *button_exec;
};

/*methods of interest*/

void uspex_gui_refresh();
