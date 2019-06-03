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

/* header for VASP file and calcul */

typedef enum {
        VP_NORM=0,
        VP_SINGLE,
        VP_ACCURATE,
        VP_HIGH=-1,
        VP_MED=-2,
        VP_LOW=-3
} vasp_prec;

typedef enum {
	VA_IALGO=-1,
        VA_NORM=0,
        VA_VERYFAST,
        VA_FAST,
        VA_CONJ,
        VA_ALL,
        VA_DAMPED,
        VA_SUBROT,
        VA_EIGEN,
        VA_NONE,
        VA_NOTHING,
        VA_EXACT,
        VA_DIAG,
} vasp_algo;

typedef enum {
        VIA_OE_FIXED=2,
        VIA_O_FIXED=3,
        VIA_SUBROT=4,
        VIA_STEEP=5,
        VIA_CG=6,
        VIA_PSTEEP=7,
        VIA_PCG=8,
        VIA_KOSUGI=38,
        VIA_ESTEEP=44,
        VIA_RMMP=46,
        VIA_PRMM=48,
        VIA_VAR_DAMP=53,
        VIA_VAR_QUENCH=54,
        VIA_VAR_PCG=58,
        VIA_EXACT=90
} vasp_ialgo;

typedef enum {
	VLR_AUTO,
	VLR_ON,
	VLR_TRUE,
	VLR_FALSE,
} vasp_lreal;

typedef enum {
	VG_91,		/*PW91*/
	VG_PE,		/*PBE*/
	VG_RP,		/*rPBE*/
	VG_AM,		/*AM05*/
	VG_PS,		/*PBEsol*/
} vasp_gga;

typedef enum {
	VMG_TPSS,
	VMG_RTPSS,
	VMG_M06L,
	VMG_MBJ,
} vasp_mgga;

typedef enum {
	VU_1=1,		/*LSDA+U Liechtenstein*/
	VU_2=2,		/*LSDA+U Dudarev*/
	VU_4=4,		/*LDA+U Liechtenstein*/
} vasp_ldau;

typedef enum {
	VUO_0=0,	/*silent*/
	VUO_1=1,	/*occupency matrix*/
	VUO_2=2,	/*full*/
} vasp_ldau_output;

typedef enum {
	VM_0=0,		/*no mixer*/
	VM_1=1,		/*Kerker mixer*/
	VM_2=2,		/*Tchebychev mixer*/
	VM_4=4,		/*Broyden mixer*/
} vasp_mixer;

typedef enum {
	VIM_0=0,	/*linear*/
	VIM_1=1,	/*Kerker*/
} vasp_inimix;

typedef enum {
	VMP_0=0,	/*none*/
	VMP_1=1,	/*inverse Kerker factor 20*/
	VMP_2=2,	/*inverse Kerker factor 200*/
} vasp_mixpre;

typedef enum {
	VID_0=0,	/*no calculation*/
	VID_1=1,	/*dipole in the lattice u axis*/
	VID_2=2,	/*dipole in the lattice v axis*/
	VID_3=3,	/*dipole in the lattice w axis*/
	VID_4=4,	/*dipole in the all lattice axis*/
} vasp_idipol;

typedef enum {
	VPF_FIXED,	/*all atoms fixed ie. F   F   F*/
	VPF_FREE,	/*all atoms free  ie. T   T   T*/
	VPF_MAN,	/*manual selection*/
} vasp_poscar_free;

typedef enum {
	VKP_MAN,	/*Manual entry*/
	VKP_LINE,	/*Line mode (bandstructure calculation)*/
	VKP_AUTO,	/*Automatic M&P generation*/
	VKP_GAMMA,	/*Gamma centered M&P generation*/
	VKP_MP,		/*Monkhorst-Pack classic generation*/
	VKP_BASIS,	/*Advanced Basis Set definition*/
} vasp_kpoints_mode;

/* global calcualtion structure */
typedef struct {
/*name*/
	gchar *name;
/*I-electronic*/
/*I-1-general*/
	vasp_prec prec;
	gboolean use_prec;
	gdouble encut;
	gdouble enaug;
	gdouble ediff;
	vasp_algo algo;
	vasp_ialgo ialgo;
	gboolean ldiag;
	gboolean auto_elec;
	gint nbands;
	gdouble nelect;
	gint iwavpr;
	gint nsim;
	gdouble vtime;
/*I-2-smearing*/
	gint ismear;
	gdouble sigma;
	gchar *fermwe;
	gchar *fermdo;
	gdouble kspacing;
	gboolean kgamma;
/*I-3-projector*/
	vasp_lreal lreal;
	gchar *ropt;
	gboolean addgrid;
	gint lmaxmix;
	gint lmaxpaw;
/*I-4-startup*/
	gint istart;
	gint icharg;
	gboolean iniwav;
/*I-5-spin*/
	gboolean ispin;
	gboolean non_collinear;
	gchar *magmom;
	gdouble nupdown;
	/*advanced*/
	gboolean lsorbit;
	gchar *saxis;
	gboolean gga_compat;
/*I-6-xc*/
	vasp_gga gga;
	gboolean voskown;
	gboolean lasph;
	gint lmaxtau;
	gboolean lmixtau;
	gboolean lmetagga;
	vasp_mgga mgga;
	gchar *cmbj;
	gdouble cmbja;
	gdouble cmbjb;
	gboolean ldau;
	vasp_ldau ldau_type;
	vasp_ldau_output ldau_output;
	gchar *ldaul;
	gchar *ldauu;
	gchar *ldauj;
	/* vdw -> advanced */
/* TODO: NOT READY YET */
	gboolean use_vdw;
	gdouble zab_vdw;
	gdouble param1_vdw;
	gdouble param2_vdw;
	gdouble param3_vdw;
/*I-7-convergence*/
	gint nelm;
	gint nelmdl;
	gint nelmin;
/*I-8-mixer*/
	gboolean auto_mixer;
	gdouble amix;
	gdouble bmix;
	gdouble amin;
	gdouble amix_mag;
	gdouble bmix_mag;
	/*advanced*/
	vasp_mixer imix;
	gint maxmix;
	gdouble wc;
	vasp_inimix inimix;
	vasp_mixpre mixpre;
/*I-9-dipole*/
	gboolean ldipol;
	gboolean lmono;
	vasp_idipol idipol;
	gdouble epsilon;
	gchar *dipol;
	gdouble efield;
/*I-10-grid*/
	gboolean auto_grid;
	gint ngx;
	gint ngy;
	gint ngz;
	gint ngxf;
	gint ngyf;
	gint ngzf;
/*II-ionic*/
/*II-1-general*/
	gint nsw;
	gint ibrion;
	gint isif;
	gdouble pstress;
	gdouble ediffg;
	gint nfree;
	gdouble potim;
/*II-2-MD*/
	gdouble tebeg;
	gdouble teend;
	gdouble smass;
	gint nblock;
	gint kblock;
	gint npaco;
	gdouble apaco;
/*II-3-symmetry*/
	gint isym;
	gdouble sym_prec;
/*II-4-POSCAR*/
	gboolean poscar_sd;/*selective dynamics*/
	gboolean poscar_direct;/*direct coord.*/
	vasp_poscar_free poscar_free;
	gdouble poscar_a0;
	gdouble poscar_ux;
        gdouble poscar_uy;
        gdouble poscar_uz;
        gdouble poscar_vx;
        gdouble poscar_vy;
        gdouble poscar_vz;
        gdouble poscar_wx;
        gdouble poscar_wy;
        gdouble poscar_wz;
	gint poscar_index;
	gchar poscar_symbol[3];
/*unexposed species and species number*/
	gint atoms_total;
	gint electron_total;
	gint species_total;
	gchar *species_symbols;
	gchar *species_numbers;
/*back to regular values*/
	gdouble poscar_x;
	gdouble poscar_y;
	gdouble poscar_z;
	gboolean poscar_tx;
	gboolean poscar_ty;
	gboolean poscar_tz;
/*III-KPOINTS*/
	gint kpoints_nkpts;
	vasp_kpoints_mode kpoints_mode;
	gboolean kpoints_cart;
	gboolean kpoints_tetra;
	gdouble kpoints_kx;
	gdouble kpoints_ky;
	gdouble kpoints_kz;
	gdouble kpoints_sx;
	gdouble kpoints_sy;
	gdouble kpoints_sz;

        gint kpoints_i;
        gdouble kpoints_x;
        gdouble kpoints_y;
        gdouble kpoints_z;
        gdouble kpoints_w;
	
	gint tetra_total;
	gdouble tetra_volume;
	gint tetra_i;
	gdouble tetra_w;
	gint tetra_a;
	gint tetra_b;
	gint tetra_c;
	gint tetra_d;
/*IV-POTCAR*/
	gchar *potcar_folder;
	gchar *potcar_file;
	gchar *potcar_species;
	gchar *potcar_species_flavor;
/*IV-advanced*/
/*IV-1-DOS*/
	gint lorbit;
	gboolean have_paw;
	gchar *rwigs;
	gint nedos;
	gdouble emin;
	gdouble emax;
	gdouble efermi;
/*IV-2-linear response*/
	gboolean loptics;
	gboolean lcalceps;
	gboolean lepsilon;
	gboolean lrpa;
	gboolean lnabla;
	gdouble cshift;
/* (adv.) */
/*V-1-performance*/
	gint npar;
	gdouble ncore;
	gdouble kpar;
	gboolean lplane;
	gboolean lscalu;
	gboolean lscalapack;
/*V-2-writting*/
	gint nwrite;
	gboolean lwave;
	gboolean lcharg;
	gboolean lvtot;
	gboolean lvhar;
	gboolean lelf;

/*VI-run*/
	gchar *job_vasp_exe;
	gchar *job_mpirun;
	gchar *job_path;
	gdouble job_nproc;
/* RESULTS */
} vasp_calc_struct;
typedef enum {
	VASP_SINGLE=0,	/*single-point electronic calculation*/
	VASP_RUN=1,	/*currently running calculation*/
	VASP_OPT=2,	/*geometry optimization calculation*/
	VASP_ELEC=4,	/*electronic structure calculation*/
	VASP_FREQ=8,	/*vibrational modes calculation*/
} vasp_calc_type;
/* output structure */
typedef struct {
/*name*/
	gchar *name;
	gint version;
//	vasp_calc_struct *calc;/*calculation structure*/
	/*extra information*/
	vasp_calc_type calc_type;
	gint n_scf;
	gdouble *E;
	gdouble *V;
	gdouble *F;
	gdouble *P;
	long int last_pos;
/*NEW: attach frame-dependent graphs*/
	gpointer graph_energy;
	gpointer graph_volume;
	gpointer graph_forces;
	gpointer graph_stress;
} vasp_output_struct;
/*execution structure*/
typedef struct {
	gchar *job_vasp_exe;
	gchar *job_mpirun;
	gchar *job_path;
	gdouble job_nproc;
} vasp_exec_struct;

/*methods of interest*/

int vasp_xml_load_calc(FILE *vf,vasp_calc_struct *calc);
void vasp_calc_to_incar(FILE *output,vasp_calc_struct calc);
gint vasprun_update(gchar *filename,vasp_calc_struct *calc);
void vasprun_free(vasp_calc_struct *calc);
void vasp_out_free(vasp_output_struct * vo);

gint vasp_load_poscar5(FILE *vf,struct model_pak *model);

