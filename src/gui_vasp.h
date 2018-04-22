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
        tamp=g_strdup_printf("%s",gtk_entry_get_text(GTK_ENTRY(vasp_gui.value)));\
        sscanf(tamp,format,&(vasp_gui.calc.value));\
        g_free(tamp);\
}while(0)
#define VASP_REG_TEXT(value) do{\
if(gtk_entry_get_text_length(GTK_ENTRY(vasp_gui.value))>0){\
        if(vasp_gui.calc.value!=NULL) g_free(vasp_gui.calc.value);\
        vasp_gui.calc.value=g_strdup_printf("%s",gtk_entry_get_text(GTK_ENTRY(vasp_gui.value)));\
}else{\
        if(vasp_gui.calc.value!=NULL) g_free(vasp_gui.calc.value);\
        vasp_gui.calc.value=NULL;\
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
/* DEFINES: table interface */
#define VASP_TEXT_TABLE(name,value,caption,l,r,t,b) do{\
        hbox = gtk_hbox_new(FALSE, 0);\
        VASP_NEW_LABEL(caption);\
	VASP_TEXT_ENTRY(name,value,TRUE);\
        gtk_table_attach_defaults(GTK_TABLE(table),hbox,l,r,t,b);\
}while(0)
#define VASP_ENTRY_TABLE(name,value,format,caption,l,r,t,b) do{\
	hbox = gtk_hbox_new(FALSE, 0);\
	VASP_NEW_LABEL(caption);\
	VASP_NEW_SEPARATOR();\
	VASP_ENTRY(name,value,format,8);\
	gtk_table_attach_defaults(GTK_TABLE(table),hbox,l,r,t,b);\
}while(0)
#define VASP_CHECK_TABLE(name,value,function,caption,l,r,t,b) do{\
	hbox = gtk_hbox_new(FALSE, 0);\
	name = gui_direct_check(caption,&(value),function,NULL,hbox);\
	gtk_table_attach_defaults(GTK_TABLE(table),hbox,l,r,t,b);\
}while(0)
#define VASP_COMBO_TABLE(name,default_text,function,caption,l,r,t,b) do{\
	hbox = gtk_hbox_new(FALSE, 0);\
	VASP_NEW_LABEL(caption);\
	VASP_NEW_SEPARATOR();\
	VASP_REG_COMBO(name,default_text,function);\
	gtk_table_attach_defaults(GTK_TABLE(table),hbox,l,r,t,b);\
}while(0)
#define VASP_SEPARATOR_TABLE(l,r,t,b) do{\
        separator=gtk_hseparator_new();\
	gtk_table_attach_defaults(GTK_TABLE(table),separator,l,r,t,b);\
}while(0)
#define VASP_COMBOBOX_TABLE(name,caption,l,r,t,b) do{\
        hbox = gtk_hbox_new(FALSE, 0);\
        label = gtk_label_new(caption);\
        gtk_box_pack_start(GTK_BOX(hbox),label,FALSE,FALSE,0);\
        name=gtk_combo_box_text_new_with_entry();\
        gtk_entry_set_editable(GTK_ENTRY(gtk_bin_get_child(GTK_BIN(name))),FALSE);\
        gtk_box_pack_start(GTK_BOX(hbox),name,TRUE,TRUE,0);\
        gtk_table_attach_defaults(GTK_TABLE(table),hbox,l,r,t,b);\
}while(0)
#define VASP_COMBOBOX_ADD(combobox,text) do{\
	gtk_combo_box_text_append_text(GTK_COMBO_BOX_TEXT(combobox),text);\
}while(0)
#define VASP_COMBOBOX_SETUP(combobox,default_value,function) do{\
	gtk_combo_box_set_active(GTK_COMBO_BOX(combobox),default_value);\
	g_signal_connect(GTK_COMBO_BOX_TEXT(combobox),"changed",GTK_SIGNAL_FUNC(function),data);\
}while(0)
/*default spin is always [1,100]*/
#define VASP_SPIN_TABLE(name,value,function,caption,l,r,t,b) do{\
	hbox = gtk_hbox_new(FALSE, 0);\
	name = gui_direct_spin(caption,&(value),1.0,100.0,1.0,function,NULL,hbox);\
	gtk_table_attach_defaults(GTK_TABLE(table),hbox,l,r,t,b);\
}while(0)
#define VASP_BUTTON_TABLE(name,type,function,l,r,t,b) do{\
        name=gtk_button_new_from_stock(type);\
	gtk_table_attach_defaults(GTK_TABLE(table),name,l,r,t,b);\
        g_signal_connect(GTK_OBJECT(name),"clicked",GTK_SIGNAL_FUNC(function),NULL);\
}while(0)
#define VASP_2BUTTONS_TABLE(name_1,name_2,function_1,function_2,l,r,t,b) do{\
	hbox = gtk_hbox_new(FALSE, 0);\
	name_1 = gui_stock_button(GTK_STOCK_APPLY,function_1,NULL,hbox);\
	name_2 = gui_stock_button(GTK_STOCK_DELETE,function_2,NULL,hbox);\
	gtk_table_attach_defaults(GTK_TABLE(table),hbox,l,r,t,b);\
}while(0)
/*left aligned labels*/
#define VASP_LABEL_TABLE(caption,l,r,t,b) do{\
	hbox = gtk_hbox_new(FALSE, 0);\
	label = gtk_label_new(caption);\
	gtk_box_pack_start(GTK_BOX(hbox),label,FALSE,FALSE,0);\
	VASP_NEW_SEPARATOR();\
	gtk_table_attach_defaults(GTK_TABLE(table),hbox,l,r,t,b);\
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
	GtkWidget *window;
	vasp_calc_struct calc;
	/*actual GUI*/
/*simple interface*/
	GtkWidget *simple_calcul;
	GtkWidget *simple_system;
	gboolean simple_rgeom;
	GtkWidget *simple_dim;
	gdouble dimension;
	GtkWidget *simple_kgrid;
	GtkWidget *simple_poscar;
	GtkWidget *simple_species;
	GtkWidget *simple_potcar;
	GtkWidget *simple_potcar_button;
	GtkWidget *simple_np;
	GtkWidget *simple_ncore;
	GtkWidget *simple_kpar;
	GtkWidget *simple_apply;
	GtkWidget *simple_message;
	GtkTextBuffer *simple_message_buff;
/*full interface*/
	gint cur_page;
	GtkWidget *spinner;
	GtkWidget *name;
	GtkWidget *file_entry;
	GtkWidget *prec;
	GtkWidget *encut;
	GtkWidget *enaug;
	GtkWidget *ediff;
	GtkWidget *algo;
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
	GtkWidget *lreal;
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
	GtkWidget *gga;
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
/* dipol */
	GtkWidget *idipol;
	GtkWidget *ldipol;
	GtkWidget *lmono;
	GtkWidget *dipol;
	GtkWidget *epsilon;
	GtkWidget *efield;
/* dos */
	GtkWidget *lorbit;
	GtkWidget *have_paw;
	GtkWidget *nedos;
	GtkWidget *emin;
	GtkWidget *emax;
	GtkWidget *efermi;
	GtkWidget *rwigs;
/* linear response */
	GtkWidget *loptics;
	GtkWidget *lepsilon;
	GtkWidget *lrpa;
	GtkWidget *lnabla;
	GtkWidget *lcalceps;
	GtkWidget *cshift;
/* grid */
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
	GtkWidget *poscar_species;
	GtkWidget *species_flavor;
	GtkWidget *species_button;
	gboolean have_potcar_folder;
	GtkWidget *potcar_select_file;
	GtkWidget *potcar_file;
	GtkWidget *potcar_file_button;
	GtkWidget *potcar_select_folder;
	GtkWidget *potcar_folder;
	GtkWidget *potcar_folder_button;
	GtkWidget *potcar_species;
	GtkWidget *potcar_species_flavor;
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
/*buttons*/
	GtkWidget *button_save;
	GtkWidget *button_exec;
/* RESULT */
	gint link_job_id;
};

/*methods of interest*/

void vasp_gui_refresh();
void vasp_poscar_sync();
