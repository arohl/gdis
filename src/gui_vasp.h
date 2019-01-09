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

/* simple VASP calcul interface */
/* DEFINES: sync vasp_gui and vasp_gui.calc */
#define VASP_REG_VAL(value,format) do{\
	GUI_REG_VAL(vasp_gui.value,vasp_gui.calc.value,format);\
}while(0)
#define VASP_REG_TEXT(value) do{\
	if(vasp_gui.calc.value!=NULL) g_free(vasp_gui.calc.value);\
	if(GUI_ENTRY_LENGTH(vasp_gui.value)>0) GUI_ENTRY_GET_TEXT(vasp_gui.value,vasp_gui.calc.value);\
	else vasp_gui.calc.value=NULL;\
}while(0)
/*page numbers*/
#define VASP_PAGE_SIMPLIFIED 0
#define VASP_PAGE_CONVERGENCE 1
#define VASP_PAGE_ELECTRONIC 2
#define VASP_PAGE_ELEC2 3
#define VASP_PAGE_IONIC 4
#define VASP_PAGE_KPOINTS 5
#define VASP_PAGE_POTCAR 6
#define VASP_PAGE_EXEC 7
/* gui structure */
struct vasp_calc_gui{
	/*just a recall*/
	GUI_OBJ *window;
	vasp_calc_struct calc;
	/*actual GUI*/
/*simple interface*/
	GUI_OBJ *simple_calcul;
	GUI_OBJ *simple_system;
	gboolean simple_rgeom;
	GUI_OBJ *simple_dim;
	gdouble dimension;
	GUI_OBJ *simple_kgrid;
	GUI_OBJ *simple_poscar;
	GUI_OBJ *simple_species;
	GUI_OBJ *simple_potcar;
	GUI_OBJ *simple_potcar_button;
	GUI_OBJ *simple_np;
	GUI_OBJ *simple_ncore;
	GUI_OBJ *simple_kpar;
	GUI_OBJ *simple_apply;
	GUI_OBJ *simple_message;
	GUI_TEXTVIEW_BUFFER *simple_message_buff;
/*full interface*/
	gint cur_page;
	GUI_OBJ *name;
	GUI_OBJ *file_entry;
	GUI_OBJ *prec;
	GUI_OBJ *encut;
	GUI_OBJ *enaug;
	GUI_OBJ *ediff;
	GUI_OBJ *algo;
	GUI_OBJ *ialgo;
	GUI_OBJ *nsim;
	GUI_OBJ *vtime;
	GUI_OBJ *nbands;
	GUI_OBJ *nelect;
	GUI_OBJ *iwavpr;
	GUI_OBJ *ismear;
	GUI_OBJ *ismear_3;/*intentional duplicate*/
	GUI_OBJ *sigma;
	GUI_OBJ *fermwe;
	GUI_OBJ *fermdo;
	GUI_OBJ *kgamma;
	GUI_OBJ *kgamma_3;/*intentional duplicates*/
	GUI_OBJ *kspacing;
	GUI_OBJ *kspacing_3;/*intentional duplicates*/
	GUI_OBJ *lreal;
	GUI_OBJ *ropt;
	GUI_OBJ *lmaxmix;
	GUI_OBJ *lmaxpaw;
	GUI_OBJ *istart;
	GUI_OBJ *icharg;
	GUI_OBJ *nupdown;
	GUI_OBJ *magmom;
	GUI_OBJ *lmaxtau;
	GUI_OBJ *lnoncoll;
	GUI_OBJ *saxis;
	GUI_OBJ *gga;
	GUI_OBJ *lsorbit;
	GUI_OBJ *metagga;
	GUI_OBJ *cmbj;
	GUI_OBJ *cmbja;
	GUI_OBJ *cmbjb;
	GUI_OBJ *ldau;
	GUI_OBJ *ldau_print;
	GUI_OBJ *ldaul;
	GUI_OBJ *ldauu;
	GUI_OBJ *ldauj;
/*mixer*/
	GUI_OBJ *nelm;
	GUI_OBJ *nelmdl;
	GUI_OBJ *nelmin;
	GUI_OBJ *mixer;
	GUI_OBJ *amix;
	GUI_OBJ *bmix;
	GUI_OBJ *amin;
	GUI_OBJ *amix_mag;
	GUI_OBJ *bmix_mag;
	GUI_OBJ *maxmix;
	GUI_OBJ *wc;
	GUI_OBJ *inimix;
	GUI_OBJ *mixpre;
/* dipol */
	GUI_OBJ *idipol;
	GUI_OBJ *ldipol;
	GUI_OBJ *lmono;
	GUI_OBJ *dipol;
	GUI_OBJ *epsilon;
	GUI_OBJ *efield;
/* dos */
	GUI_OBJ *lorbit;
	GUI_OBJ *have_paw;
	GUI_OBJ *nedos;
	GUI_OBJ *emin;
	GUI_OBJ *emax;
	GUI_OBJ *efermi;
	GUI_OBJ *rwigs;
/* linear response */
	GUI_OBJ *loptics;
	GUI_OBJ *lepsilon;
	GUI_OBJ *lrpa;
	GUI_OBJ *lnabla;
	GUI_OBJ *lcalceps;
	GUI_OBJ *cshift;
/* grid */
	GUI_OBJ *ngx;
	GUI_OBJ *ngy;
	GUI_OBJ *ngz;
	GUI_OBJ *ngxf;
	GUI_OBJ *ngyf;
	GUI_OBJ *ngzf;
/*ionic*/
	GUI_OBJ *nsw;
	GUI_OBJ *ibrion;
	GUI_OBJ *isif;
	GUI_OBJ *relax_ions;
	GUI_OBJ *relax_shape;
	GUI_OBJ *relax_volume;
	GUI_OBJ *ediffg;
	GUI_OBJ *pstress;
	GUI_OBJ *nfree;
	GUI_OBJ *potim;
/*md*/
	GUI_OBJ *tebeg;
	GUI_OBJ *teend;
	GUI_OBJ *smass;
	GUI_OBJ *nblock;
	GUI_OBJ *kblock;
	GUI_OBJ *npaco;
	GUI_OBJ *apaco;
/*symmetry*/
	GUI_OBJ *isym;
	GUI_OBJ *sym_prec;
/* POSCAR */
	GUI_OBJ *poscar_free;
	GUI_OBJ *poscar_direct;
	GUI_OBJ *poscar_a0;
	GUI_OBJ *poscar_ux;
        GUI_OBJ *poscar_uy;
        GUI_OBJ *poscar_uz;
        GUI_OBJ *poscar_vx;
        GUI_OBJ *poscar_vy;
        GUI_OBJ *poscar_vz;
        GUI_OBJ *poscar_wx;
        GUI_OBJ *poscar_wy;
        GUI_OBJ *poscar_wz;
	GUI_OBJ *poscar_atoms;
	GUI_OBJ *poscar_index;
	GUI_OBJ *poscar_symbol;
	GUI_OBJ *poscar_x;
	GUI_OBJ *poscar_y;
	GUI_OBJ *poscar_z;
	GUI_OBJ *poscar_tx;
	GUI_OBJ *poscar_ty;
	GUI_OBJ *poscar_tz;
/* KPOINTS */
	GUI_OBJ *kpoints_nkpts;
	GUI_OBJ *kpoints_mode;
	GUI_OBJ *kpoints_cart;
	GUI_OBJ *kpoints_kx;
        GUI_OBJ *kpoints_ky;
        GUI_OBJ *kpoints_kz;
	GUI_OBJ *kpoints_sx;
	GUI_OBJ *kpoints_sy;
	GUI_OBJ *kpoints_sz;
	GUI_OBJ *kpoints_coord;
	GUI_OBJ *have_tetra;
	GUI_OBJ *kpoints_kpts;
	GUI_OBJ *kpoints_i;
	GUI_OBJ *kpoints_x;
	GUI_OBJ *kpoints_y;
	GUI_OBJ *kpoints_z;
	GUI_OBJ *kpoints_w;
	GUI_OBJ *kpoints_tetra;
	GUI_OBJ *tetra;
	GUI_OBJ *tetra_total;
	GUI_OBJ *tetra_volume;
	GUI_OBJ *tetra_i;
	GUI_OBJ *tetra_w;
	GUI_OBJ *tetra_a;
	GUI_OBJ *tetra_b;
	GUI_OBJ *tetra_c;
	GUI_OBJ *tetra_d;
/* POTCAR */
	GUI_OBJ *poscar_species;
	GUI_OBJ *species_flavor;
	GUI_OBJ *species_button;
	gboolean have_potcar_folder;
	GUI_OBJ *potcar_select_file;
	GUI_OBJ *potcar_file;
	GUI_OBJ *potcar_file_button;
	GUI_OBJ *potcar_select_folder;
	GUI_OBJ *potcar_folder;
	GUI_OBJ *potcar_folder_button;
	GUI_OBJ *potcar_species;
	GUI_OBJ *potcar_species_flavor;
/* PERF */
	GUI_OBJ *ncore;
	GUI_OBJ *kpar;
/*switches*/
	gboolean have_xml;
	gboolean poscar_dirty;
	gboolean rions;
	gboolean rshape;
	gboolean rvolume;
	gboolean have_result;
/* CALCUL */
	GUI_OBJ *job_vasp_exe;/*will be taken directly from gdis in the future*/
	GUI_OBJ *job_mpirun;
	GUI_OBJ *job_path;
	GUI_OBJ *job_nproc;
/*buttons*/
	GUI_OBJ *button_save;
	GUI_OBJ *button_exec;
};

/*methods of interest*/

void vasp_gui_refresh();
void vasp_poscar_sync();
