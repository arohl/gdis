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
/* DEFINES: sync vasp_gui and vasp_calc */
#define VASP_REG_VAL(value,format) do{\
        tamp=g_strdup_printf("%s",gtk_entry_get_text(GTK_ENTRY(vasp_gui.value)));\
        sscanf(tamp,format,&(vasp_calc.value));\
        g_free(tamp);\
}while(0)
#define VASP_REG_TEXT(value) do{\
if(gtk_entry_get_text_length(GTK_ENTRY(vasp_gui.value))>0){\
        if(vasp_calc.value!=NULL) g_free(vasp_calc.value);\
        vasp_calc.value=g_strdup_printf("%s",gtk_entry_get_text(GTK_ENTRY(vasp_gui.value)));\
}else{\
        if(vasp_calc.value!=NULL) g_free(vasp_calc.value);\
        vasp_calc.value=NULL;\
}\
}while(0)
/* DEFINES: interface */
#define _SPACE 0
#define VASP_NEW_LINE() do{\
        hbox = gtk_hbox_new(FALSE, 0);\
        gtk_container_set_border_width(GTK_CONTAINER(hbox), _SPACE);\
        gtk_box_pack_start(GTK_BOX(vbox), hbox, FALSE, FALSE, 0);\
}while(0)
#define VASP_NEW_LABEL(text) do{\
        label = gtk_label_new(text);\
        gtk_box_pack_start(GTK_BOX(hbox), label, FALSE, FALSE, 0);\
}while(0)
#define VASP_TEXT_ENTRY(name,value,wilde) do{\
        name=gtk_entry_new();\
        if(value!=NULL) gtk_entry_set_text(GTK_ENTRY(name),g_strdup_printf("%s",value));\
        gtk_box_pack_start(GTK_BOX(hbox), name, wilde, wilde, 0);\
}while(0)
#define VASP_ENTRY(name,value,format,size) do{\
        name=gtk_entry_new();\
        gtk_entry_set_text(GTK_ENTRY(name),g_strdup_printf(format,value));\
        gtk_entry_set_width_chars (GTK_ENTRY(name),size);\
        gtk_box_pack_start(GTK_BOX(hbox), name, FALSE, FALSE, 0);\
}while(0)
#define VASP_NEW_SEPARATOR() do{\
        separator=gtk_hseparator_new();\
        gtk_box_pack_start(GTK_BOX(hbox), separator, TRUE, FALSE, 0);\
}while(0)
#define VASP_REG_COMBO(name,default_text,function) do{\
        name = gtk_combo_new();\
        gtk_entry_set_editable(GTK_ENTRY(GTK_COMBO(name)->entry), FALSE);\
        gtk_combo_set_popdown_strings(GTK_COMBO(name), list);\
        g_list_free(list);\
        gtk_box_pack_start(GTK_BOX(hbox), name, FALSE, FALSE, 0);\
        gtk_entry_set_text(GTK_ENTRY(GTK_COMBO(name)->entry),default_text);\
        g_signal_connect(GTK_OBJECT(GTK_COMBO(name)->entry), "changed",GTK_SIGNAL_FUNC(function),data);\
}while(0)
/*page numbers*/
#define VASP_PAGE_SIMPLIFIED 0
#define VASP_PAGE_CONVERGENCE 1
#define VASP_PAGE_ELECTRONIC 2
#define VASP_PAGE_IONIC 3
#define VASP_PAGE_KPOINTS 4
#define VASP_PAGE_POTCAR 5
#define VASP_PAGE_EXEC 6

/* gui structure */
struct vasp_calc_gui{
	/*just a recall*/
	GtkWidget *window;
	/*actual GUI*/
	GtkWidget *spinner;
	GtkWidget *name;
	GtkWidget *file_entry;
	GtkWidget *encut;
	GtkWidget *enaug;
	GtkWidget *ediff;
	GtkWidget *ialgo;
	GtkWidget *nsim;
	GtkWidget *vtime;
	GtkWidget *nbands;
	GtkWidget *nelect;
	GtkWidget *iwavpr;
	GtkWidget *ismear;
	GtkWidget *ismear_3;/*intentional duplicate*/
	GtkWidget *sigma;
	GtkWidget *fermwe;
	GtkWidget *fermdo;
	GtkWidget *kgamma;
	GtkWidget *kgamma_3;/*intentional duplicates*/
	GtkWidget *kspacing;
	GtkWidget *kspacing_3;/*intentional duplicates*/
	GtkWidget *ropt;
	GtkWidget *lmaxmix;
	GtkWidget *lmaxpaw;
	GtkWidget *istart;
	GtkWidget *icharg;
	GtkWidget *nupdown;
	GtkWidget *magmom;
	GtkWidget *lmaxtau;
	GtkWidget *lnoncoll;
	GtkWidget *saxis;
	GtkWidget *lsorbit;
	GtkWidget *metagga;
	GtkWidget *cmbj;
	GtkWidget *cmbja;
	GtkWidget *cmbjb;
	GtkWidget *ldau;
	GtkWidget *ldau_print;
	GtkWidget *ldaul;
	GtkWidget *ldauu;
	GtkWidget *ldauj;
/* ... */
/*mixer*/
	GtkWidget *nelm;
	GtkWidget *nelmdl;
	GtkWidget *nelmin;
	GtkWidget *mixer;
	GtkWidget *amix;
	GtkWidget *bmix;
	GtkWidget *amin;
	GtkWidget *amix_mag;
	GtkWidget *bmix_mag;
	GtkWidget *maxmix;
	GtkWidget *wc;
	GtkWidget *inimix;
	GtkWidget *mixpre;
/* ... */
	GtkWidget *idipol;
	GtkWidget *ldipol;
	GtkWidget *lmono;
	GtkWidget *dipol;
	GtkWidget *epsilon;
	GtkWidget *efield;

	GtkWidget *ngx;
	GtkWidget *ngy;
	GtkWidget *ngz;
	GtkWidget *ngxf;
	GtkWidget *ngyf;
	GtkWidget *ngzf;
/*ionic*/
	GtkWidget *nsw;
	GtkWidget *ibrion;
	GtkWidget *isif;
	GtkWidget *relax_ions;
	GtkWidget *relax_shape;
	GtkWidget *relax_volume;
	GtkWidget *ediffg;
	GtkWidget *pstress;
	GtkWidget *nfree;
	GtkWidget *potim;
/*md*/
	GtkWidget *tebeg;
	GtkWidget *teend;
	GtkWidget *smass;
	GtkWidget *nblock;
	GtkWidget *kblock;
	GtkWidget *npaco;
	GtkWidget *apaco;
/*symmetry*/
	GtkWidget *isym;
	GtkWidget *sym_prec;
/* POSCAR */
	GtkWidget *poscar_free;
	GtkWidget *poscar_direct;
	GtkWidget *poscar_a0;
	GtkWidget *poscar_ux;
        GtkWidget *poscar_uy;
        GtkWidget *poscar_uz;
        GtkWidget *poscar_vx;
        GtkWidget *poscar_vy;
        GtkWidget *poscar_vz;
        GtkWidget *poscar_wx;
        GtkWidget *poscar_wy;
        GtkWidget *poscar_wz;
	GtkWidget *poscar_atoms;
	GtkWidget *poscar_index;
	GtkWidget *poscar_symbol;
	GtkWidget *poscar_x;
	GtkWidget *poscar_y;
	GtkWidget *poscar_z;
	GtkWidget *poscar_tx;
	GtkWidget *poscar_ty;
	GtkWidget *poscar_tz;
/* KPOINTS */
	GtkWidget *kpoints_nkpts;
	GtkWidget *kpoints_mode;
	GtkWidget *kpoints_cart;
	GtkWidget *kpoints_kx;
        GtkWidget *kpoints_ky;
        GtkWidget *kpoints_kz;
	GtkWidget *kpoints_sx;
	GtkWidget *kpoints_sy;
	GtkWidget *kpoints_sz;
	GtkWidget *kpoints_coord;
	GtkWidget *have_tetra;
	GtkWidget *kpoints_kpts;
	GtkWidget *kpoints_i;
	GtkWidget *kpoints_x;
	GtkWidget *kpoints_y;
	GtkWidget *kpoints_z;
	GtkWidget *kpoints_w;
	GtkWidget *kpoints_tetra;
	GtkWidget *tetra;
	GtkWidget *tetra_total;
	GtkWidget *tetra_volume;
	GtkWidget *tetra_i;
	GtkWidget *tetra_w;
	GtkWidget *tetra_a;
	GtkWidget *tetra_b;
	GtkWidget *tetra_c;
	GtkWidget *tetra_d;
/* POTCAR */
	GtkWidget *species;
	GtkWidget *potcar_folder;
	GtkWidget *potcar_file;
	GtkWidget *potcar_species;
/* PERF */
	GtkWidget *ncore;
	GtkWidget *kpar;
/*switches*/
	gboolean have_xml;
	gboolean poscar_dirty;
	gboolean rions;
	gboolean rshape;
	gboolean rvolume;
	gboolean have_result;
/* CALCUL */
	GtkWidget *job_vasp_exe;/*will be taken directly from gdis in the future*/
	GtkWidget *job_mpirun;
	GtkWidget *job_path;
	GtkWidget *job_nproc;
/* RESULT */
	struct model_pak *result_model;
};
/* calcul structure*/
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
	gint nelect;
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
	gint nupdown;
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
/* TODO: NOT READY YET
	gboolean use_vdw;
	gdouble zab_vdw;
	gdouble param1_vdw;
	gdouble param2_vdw;
	gdouble param3_vdw;
*/
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
	gint smass;
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

/*IV-advanced*/
/*IV-1-DOS*/
	gboolean lorbit;
	gdouble *rwigs;
	gint nedos;
	gdouble emin;
	gdouble emax;
	gdouble efermi;
/*IV-2-linear response*/
	gboolean lepsilon;
	gboolean lrpa;
	gboolean lnabla;
	gboolean lvel;
	gint kinter;
	gdouble cshift;
	gdouble omegamax;
	gdouble deg_thershold;
/*IV-3-magnetization*/
	gboolean nucind;
	gdouble magpos[3];
	gboolean lnicsall;
	gboolean orbitalmag;
	gboolean lmagbloch;
	gboolean lchimag;
	gboolean lgauge;
	gint magatom;
	gdouble magdipol[3];
	gdouble avecconst[3];
/*V-performance*/
	gdouble ncore;
	gdouble kpar;
	gboolean lplane;
	gboolean lscalu;
	gboolean lscalapack;
/**/
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

} vasp_calc_struct;


/*methods of interest*/


void vasp_calc_to_incar(FILE *output,vasp_calc_struct vasp_calc);







