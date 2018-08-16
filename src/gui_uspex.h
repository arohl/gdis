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
	GUI_OBJ *optType;
	/*atoms definition: apply/remove*/
	GUI_OBJ *atomType;
	GUI_OBJ *_atom_sym;
	GUI_OBJ *_atom_typ;
	GUI_OBJ *_atom_num;
	GUI_OBJ *_atom_val;
	/*temporary value*/
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
/*4.2 Population*/
	GUI_OBJ *populationSize;
	GUI_OBJ *initialPopSize;
	GUI_OBJ *numGenerations;
	GUI_OBJ *stopCrit;
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
