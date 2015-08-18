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

/**********************************/
/* global unified rendering setup */
/**********************************/
struct render_pak
{
/* main setup */
gdouble width;
gdouble height;
gint camera[3];
gdouble vp_dist;        /* vanishing point for perspective proj */
gchar filename[FILELEN];
/* flags & control */
gint mode;              /* default render mode for atoms */
gint perspective;
gint antialias;
gint shadowless;
gint fog;
gint axes;
gint atype;
gint show_energy;
gboolean wire_model;
gboolean wire_surface;
gboolean wire_show_hidden;
gint animate;
gint animate_type;       /* file output */
gint no_povray_exec;
gint no_keep_tempfiles;
gchar *animate_file;
gdouble delay;
gdouble mpeg_quality;

/* TODO - put all the stereo stuff here */
gboolean stereo_quadbuffer;
gboolean stereo_use_frustum;
gboolean stereo_left;
gboolean stereo_right;
gdouble stereo_eye_offset;
gdouble stereo_parallax;

/* composite object quality */
gint auto_quality;            /* switch auto quality on/off */
gdouble sphere_quality;
gdouble cylinder_quality;
gdouble ribbon_quality;
gint halos;
gint fast_rotation;

/* light list */
gint num_lights;
GSList *light_list;

/* old lighting */
gdouble ambience;
gdouble diffuse;
gdouble specular;

gdouble fog_density;
gdouble fog_start;

/* geometry */
gdouble ball_rad;
gdouble stick_rad;
gdouble stick_thickness;
gdouble line_thickness;
gdouble frame_thickness;
gdouble geom_line_width;
gdouble cpk_scale;
gdouble ribbon_thickness;
gdouble ribbon_curvature;
gdouble phonon_scaling;

/* exp grid size for zones */
gdouble zone_size;

/* molsurf */
gdouble ms_grid_size;

/* spatial highlighting */
gdouble shl_strength;
gdouble shl_size;
/* atom/bond highlighting */
gdouble ahl_strength;
gdouble ahl_size;

/* material properties*/
gdouble ref_index;
gdouble transmit;
gdouble ghost_opacity;
gdouble fg_colour[3];
gdouble bg_colour[3];
gdouble morph_colour[3];
gdouble rsurf_colour[3];         /* re-entrant surface */
gdouble label_colour[3];         /* geometry labels */
gdouble title_colour[3];         /* axes titles */
gdouble ribbon_colour[3];
gchar *morph_finish;
};

/***********************/
/* top level structure */
/***********************/
struct sysenv_pak
{
/* location info */
gchar *cwd;
gchar *init_file;
gchar *elem_file;

gchar *gdis_path;

gchar *babel_exe;
gchar *babel_path;

gchar *gulp_exe;
gchar *gulp_path;

gchar *povray_exe;
gchar *povray_path;

gchar *convert_exe;
gchar *convert_path;

gchar *viewer_exe;
gchar *viewer_path;

gchar *gamess_exe;
gchar *gamess_path;

gchar *monty_exe;
gchar *monty_path;

/* filetype filters */
gint file_type;
gint babel_type;
gint write_gdisrc;
int fortran_buffer;

/* centralized info on all supported file types */
GSList *file_list;
GSList *module_list;
GSList *projects;

/* database allocation size */
gint num_elements;
GSList *elements;
GHashTable *sfc_table;
GHashTable *library;
GHashTable *manual;
GHashTable *image_table;
GtkTreeStore *tree_store;
GtkWidget *tree;

/* approximate GTK font width */
gint gtk_fontsize;

/* OpenGL */
gchar gl_fontname[LINELEN];
GdkFont *gl_font;
gpointer glconfig;
gint snapshot;
gchar *snapshot_filename;

/* LEFT pane hideable boxes */
gint mtb_on;
gint mpb_on;
gint msb_on;
gint apb_on;
gint pib_on;
gint lmb_on;
gint mcb_on;

/* model related */
gpointer active_model;
GSList *mal;
gint num_trees;
gint num_displayed;
gint displayed[MAX_DISPLAYED];
gint select_mode;
/* NEW - perform the periodic bond calc, slower when it's done, but */
/* gives more info & speeds up subsequent stuff (eg unfrag) */
gint calc_pbonds;

/* drawing */
gboolean refresh_canvas;
gboolean refresh_dialog;
gboolean refresh_tree;
gboolean refresh_properties;
gboolean refresh_text;

gboolean stereo;                /* stereo on/off */
gboolean stereo_windowed;       /* is capable of windowed stereo */
gboolean stereo_fullscreen;     /* actual stereo mode (windowed/fullscreen) */
gint fps;
gint canvas;
gint moving;
gint roving;
gint canvas_rows;
gint canvas_cols;
GSList *canvas_list;
GdkVisual *visual;
GdkColormap *colourmap;
GtkWidget *display_box;
GtkWidget *glarea;
GtkWidget *surfaces;
GtkWidget *mpane;
GtkWidget *tpane;
/* screen dimensions */
gint x;
gint y;
gint width;
gint height;
/* model pane dimensions */
gint tree_width;
gint tray_height;
gint tree_divider;
/* for OpenGL - the smaller of width & height */
gint size;    /* pixel size */
gdouble rsize;  /* real world size */
gdouble aspect; /* aspect ratio */
/* size of colourmap */
gint depth;
/* number of subwindows */
gint subx;
gint suby;
/* pixel dimensions */
gint subwidth;
gint subheight;
/* smallest of the above two */
gint subscale;
/* fraction dimensions */
gdouble fsubwidth;
gdouble fsubheight;
/* list of subwindow centres */
gint subcenx[MAX_DISPLAYED];
gint subceny[MAX_DISPLAYED];
/* dialog list */
GSList *dialog_list;
/* copy active select */
struct model_pak *select_source;
/* rendering */
struct render_pak render;

/* task handling */
GSList *task_list;
GtkWidget *task_label;
gint max_threads;
GThreadPool *thread_pool;

/* CURRENT - task successor (remote job execution) */
GSList *host_list;

};

/********************************/
/* element properties structure */
/********************************/
struct elem_pak
{
gchar *symbol;
gchar *name;
gint number;         /* NB: keep! elem_pak is used for more than elem table */
gdouble weight;
gdouble cova;
gdouble vdw;
gdouble charge;
gdouble shell_charge;
gdouble colour[3];
};

/**********************/
/* spatial vector pak */
/**********************/
struct vec_pak
{
gdouble colour[3];
/* coords */
gdouble x[3];
gdouble rx[3];
/* normal */
gdouble n[3];
gdouble rn[3];
/* NEW - associated data/properties */
gpointer data;
};

/***********************/
/* special object pack */
/***********************/
struct object_pak
{
gint id;
gint type;          /* plane, vector, ribbon */
gpointer *data;     /* object definition data */
};

/**********************/
/* general point pack */
/**********************/
struct point_pak
{
gint type;
gint num_points;
gdouble *x[3];
};

/*****************************************************************/
/* structure for bulk computation of valid cuts & their energies */
/*****************************************************************/
struct shift_pak
{
gint locked;
gdouble shift;
gint region[2];
gdouble eatt[2];
gdouble esurf[2];
gdouble bbpa;
gboolean dipole_computed;
gdouble dipole;
gdouble gnorm;
gchar *procfile;
};

/********************************/
/* general plane data structure */
/********************************/
struct plane_pak
{
/* is this a symmetry generated plane or not */
gint primary;
gpointer parent;

/* XRD info */
gint multiplicity;
/* real, imaginary structure factor components */
gdouble f[2];

/* raw direction (should be miller, but use float just in case) */
gdouble m[3];
/* normalized/latmat'd -ve of the miller indices */
gdouble x[3];
/* hkl label */
gint index[3];
/* unrotated normal */
gdouble norm[3];
/* NEW - the surface and depth repeat vectors */
gdouble lattice[9];
/* face ranking data */
gdouble eatt_shift;
gdouble esurf_shift;
gdouble bbpa_shift;
gdouble dhkl;
gdouble eatt[2];
gdouble esurf[2];
gdouble bbpa;            /* broken bonds per unit area */
gdouble area;
/* attached shift list */
GSList *shifts;
/* drawing control flags */
gint present;
gint visible;
/* hull vertices */
GSList *vertices;
/* hkl label position */
gdouble rx[3];
};

/********************/
/* space group info */
/********************/
struct space_pak
{
gint lookup;
gint spacenum;
gint lattice;
gint pointgroup;
gint cellchoice;
gint inversion;
gint order;
gchar centering;
gchar *spacename;
gchar *latticename;
/* signs & ordering of x,y,z general positions */
gdouble **matrix;
/* const offset to general positions */
gdouble **offset;
};

/*****************/
/* symmetry info */
/*****************/
struct symop_pak
{
gint type;
gint order;
gdouble point[3];
gdouble norm[3];
gdouble matrix[9];
};

struct symmetry_pak
{
gint num_symops;
struct symop_pak *symops;
gint num_items;
gchar **items;
/* redundant? - pass around via attached data to callbacks? */
GtkWidget *summary;
GtkWidget *pg_entry;
gchar *pg_name;
};

/*************/
/* gulp data */
/*************/
struct gulp_pak
{
gint no_exec;
gint run;
gint method;
gint optimiser;
gint optimiser2;
gint switch_type;
gdouble switch_value;
gint unit_hessian;
gint qeq;
gint free;
gint zsisa;
gint prop;
gint phonon;
gint eigen;
gint compare;
gint nosym;
gint fix;
gint noautobond;
gint cycles;
gint ensemble;
gint coulomb;
gint rigid;
gchar *rigid_move;
gdouble maxcyc;
gint super[3];
gint no_esurf;     /* old gulp's can't handle sbulkenergy keyword */
gdouble energy;
gdouble esurf[2];   /* first - unre, second - re */
/*
gint num_kpoints;
gdouble kpoints[3];
*/
gchar *kpoints;
gchar *esurf_units;
gint no_eatt;      /* old gulp's can't handle dhkl keyword */
gdouble eatt[2];
gchar *eatt_units;
gdouble sbulkenergy;
gdouble sdipole;
/* electrostatics */
gdouble qsum;
gdouble gnorm;
/* dynamics info */
/*
gdouble frame_time;
gdouble frame_ke;
gdouble frame_pe;
gdouble frame_temp;
*/
/* HACK - store original structure so we can recall after reading trajectory */
gint orig_fractional;
gint orig_construct_pbc;
gdouble orig_pbc[6];
gdouble orig_latmat[9];
GSList *orig_cores;
GSList *orig_shells;
/* electrostatic potential calcs */
GSList *epot_vecs;
GSList *epot_vals;
gdouble epot_min;
gdouble epot_max;
/* widgets that need to be updated */
GtkWidget *energy_entry;
GtkWidget *sdipole_entry;
GtkWidget *esurf_entry;
GtkWidget *eatt_entry;
GtkWidget *sbe_entry;
/* HACK - store pointer to GULP dialog so can destroy after read in GULP output */
GtkWidget *dialog;

/* element/potentials data */
gchar *potentials;
gchar *elements;
gchar *species;

/* dynamics control */
gchar *timestep;
gchar *equilibration;
gchar *production;
gchar *sample;
gchar *write;
gchar *pressure;
gchar *temperature;

/* COSMO control */
gint solvation_model;
gint cosmo_sas;
gint cosmo_shape;
gint cosmo_shape_index[2];
gint cosmo_points;
gint cosmo_segments;

gdouble cosmo_smoothing;
gdouble cosmo_solvent_epsilon;
gdouble cosmo_solvent_radius;
gdouble cosmo_solvent_delta;
gdouble cosmo_solvent_rmax;

/* file output options */
gint print_charge;

/* filename control */
gchar *libfile;
gchar *temp_file;
gchar *out_file;
gchar *dump_file;
gchar *trj_file;
gchar *mov_file;

/* stuff gdis didn't understand */
gint output_extra_keywords;
gint output_extra;
gchar *extra_keywords;
gchar *extra;
};

/***************/
/* gamess data */
/***************/

typedef enum {
  GMS_ANGS, GMS_BOHR
} GMSUnits;

typedef enum {
  GMS_RUN, GMS_CHECK, GMS_DEBUG
} GMSExeType;

typedef enum {
  GMS_USER, GMS_MNDO, GMS_AM1, GMS_PM3, GMS_MINI, GMS_MIDI,
  GMS_STO, GMS_N21, GMS_N31, GMS_N311, GMS_P6311G,
  GMS_DZV, GMS_DH, GMS_TZV, GMS_MC
} GMSBasisType;

typedef enum {
  /* single coordinate options */
  GMS_ENERGY, GMS_GRADIENT, GMS_HESSIAN, GMS_PROP, GMS_MOROKUMA, GMS_TRANSITN, GMS_FFIELD, GMS_TDHF, GMS_MAKEFP,
  /*  multiple geometry options */
  GMS_OPTIMIZE, GMS_TRUDGE, GMS_SADPOINT, GMS_IRC, GMS_VSCF, GMS_DRC, GMS_GLOBOP, GMS_GRADEXTR, GMS_SURFACE
} GMSRunType;

typedef enum {
  GMS_RHF, GMS_UHF, GMS_ROHF
} GMSScfType;

typedef enum {
  GMS_QA, GMS_NR, GMS_RFO, GMS_SCHLEGEL
} GMSOptType;

struct GMS_keyword_pak
{
gchar *label;
gint id;
};

struct basis_pak
{
gchar *label;
GMSBasisType basis;
gint ngauss;
};

struct gamess_pak
{
GMSUnits units;
GMSExeType exe_type;
GMSRunType run_type;
GMSScfType scf_type;
GMSBasisType basis;
GMSOptType opt_type;
gboolean dft;
gint dft_functional;
gint ngauss;
gint converged;
gint have_CI;
gint have_CC;
gint wide_output;
gint MP_level;
gdouble total_charge;
gdouble multiplicity;
gdouble time_limit;
gdouble mwords;
gdouble num_p;
gdouble num_d;
gdouble num_f;
gint have_heavy_diffuse;
gint have_hydrogen_diffuse;
gdouble nstep;
gdouble maxit;
gdouble energy;
gdouble max_grad;
gdouble rms_grad;
gint have_energy;
gint have_max_grad;
gint have_rms_grad;
gchar *title;
gchar *srcfile;
gchar *temp_file;
gchar *out_file;

/* widgets that need to be updated */
GtkWidget *energy_entry;
};

/***************/
/* siesta data */
/***************/
struct siesta_pak
{
gint num_atoms;
gint have_energy;
gint have_max_grad;
gdouble energy;
gdouble max_grad;
gdouble epot_min;
gdouble epot_max;
gdouble eden_scale;
/* cartesian + spin grid */
gdouble cell[9];
gdouble icell[9];
gint grid[4];
gdouble *eden;
gdouble *edos;
gdouble *epot;

/*
 * Unit Details
 */

//Mass
enum siesta_mass_measurement_type { KG_MASS, G_MASS, AMU_MASS } mixing_weight_units;

//Length
enum siesta_length_measurement_type { M_LEN, CM_LEN, NM_LEN, ANG_LEN, BOHR_LEN } kgrid_cutoff_units, md_max_cg_displacement_units, md_precondition_variable_cell_units, md_fc_displ_units;

//Temp
enum siesta_temperature_measurement_type { K_TEMP, C_TEMP, F_TEMP }  md_inital_temperature_units, electronic_temperature_units, md_target_temperature_units;

//Time
enum siesta_time_measurement_type { S_TIME, FS_TIME, PS_TIME, NS_TIME } timestep_units, md_length_time_step_units, md_tau_relax_units;

//Energy
enum siesta_energy_measurement_type { J_ENRG, ERG_ENRG, EV_ENRG, MEV_ENRG, RY_ENRG, MRY_ENRG, HARTREE_ENRG, KO_ENRG, KCALMOLE_ENRG, MHARTREE_ENRG, KJMOL_ENRG, HZ_ENRG, THZ_ENRG, CM_ENRG, CMCM_ENRG } energy_shift_units, mesh_cutoff_units;

//Force
enum siesta_force_measurement_type { N_FORCE, EVANG_FORCE, RYBOHR_FORCE } md_max_force_tol_units;

//Pressure
enum siesta_pressure_measurement_type { P_PRES, MPA_PRES, GPA_PRES, ATM_PRES, BAR_PRES, MBAR_PRES, RYBOHR3_PRES, EVANG3_PRES } pressure_units, md_max_stress_tol_units, md_bulk_modulus_units, md_target_pressure_units;
//enum siesta_charge_measurement_type { C_CHRG, E_CHRG };
//enum siesta_dipole_measurement_type { CMM_DIPOLE, D_DIPOLE, DEBYE_DIPOLE, EBOHR_DIPOLE, EANG_DIPOLE };

//moment of inertia
enum siesta_mominert_measurement_type { KGM_MOMINERT, RYFS_MOMINERT } md_nose_mass_units, md_parrinello_rahman_mass_units;

//efield - what ever that is
//enum siesta_efield_measurement_type { VM_EFIELD, VMM_EFIELD, VANG_EFIELD, VBOHR_EFIELD, RYBOHRE_EFIELD, HARBOHRE_EFIELD };

enum siesta_xc_functional_type { GGA_XCFUNC, LDA_XCFUNC, LSD_XCFUNC } xc_functional_type;

enum siesta_xc_author_type { CA_XCAUTH, PZ_XCAUTH, PW_XCAUTH, PBE_XCAUTH, RPBE_XCAUTH, LYP_XCAUTH } xc_author_type;

enum siesta_solution_method_type { DIAGON_SOLMETH, ORDERN_SOLMETH } solution_method_type;

/*
 * Electronic Structure
 * Quality?
 */
gdouble split_zeta_norm;                        //PAO.SplitNorm
gdouble energy_shift;                           //PAO.EnergyShift
gdouble mesh_cutoff;
gdouble electronic_temperature;                  //ElectronicTemperature
gboolean spin_polarised;

//PAO.BasisType
enum siesta_pao_basis_type { SPLIT_PAOBT, SPLITGAUSS_PAOBT, NODES_PAOBT, NONODES_PAOBT } basis_type;

//PAO.BasisSize
enum { SZ_ZETA, DZ_ZETA, SZP_ZETA, DZP_ZETA, CUSTOM_ZETA } basis_set;

//CUSTOM_ZETA'S
gdouble custom_zeta;
gdouble custom_zeta_polarisation;
gboolean custom_zeta_checkbox;
GtkWidget * custom_zeta_frame;          //HACK HACK HACK

gboolean harris_functional;

// *If Periodic?
gboolean is_periodic;
gdouble kgrid_cutoff;

/*
 * SCF
 */
gboolean diag_divide_and_conquer;               //Diag.DivideAndConquer
gdouble eta_value;
gdouble no_of_itterations;
gdouble localisation_radius;

gdouble no_of_cycles;

gboolean pulay_mixing;                          //if pulay
gdouble mixing_weight;                          //DM.MixingWeight
gdouble no_of_pulay_matrices;                   //DM.NumberPulay

/*
 * Geometry
 */
enum { SINGLE_POINT, OPTIMISATION, MOLECULAR_DYNAMICS, PHONON_CALCULATION } run_type;

gdouble number_of_steps;
enum  {CONSTANT_PRESSURE, CONSTANT_VOLUME } constant_component;
gdouble therostat_parameter;
gdouble timestep;

GtkWidget * geom_notebook_hack;          //HACK HACK HACK

gdouble finite_diff_step_size;

gdouble pressure;

/*
 * MD run stuff
 */
enum siesta_md_run_type { CG_MDRUN, VERLET_MDRUN, NOSE_MDRUN, PARRINELLOPAHMAN_MDRUN, NOSEPARRINELLOPAHMAN_MDRUN, ANNEAL_MDRUN, FC_MDRUN, PHONON_MDRUN } md_type_of_run;

gboolean md_variable_cell;              //MD.VariableCell
gint md_num_cg_steps;                   //MD.NumCGsteps
gdouble md_max_cg_displacement;         //MD.MaxCGDispl
gdouble md_max_force_tol;               //MD.MaxForceTol
gdouble md_precondition_variable_cell;  //MD.PreconditionVariableCell
gdouble md_max_stress_tol;              //MD.MaxStressTol
gdouble md_inital_time_step;            //MD.InitialTimeStep
gdouble md_final_time_step;             //MD.FinalTimeStep
gdouble md_length_time_step;            //MD.LengthTimeStep
gboolean md_quench;                     //MD.Quench


gdouble md_inital_temperature;          //MD.InitialTemperature
gdouble md_target_temperature;          //MD.TargetTemperature
gdouble md_target_pressure;             //MD.TargetPressure
gdouble md_target_stress_xx;
gdouble md_target_stress_yy;
gdouble md_target_stress_zz;
gdouble md_target_stress_xy;
gdouble md_target_stress_xz;
gdouble md_target_stress_yz;


gdouble md_nose_mass;                   //MD.NoseMass
gdouble md_parrinello_rahman_mass;      //MD.ParrinelloRahmanMass

enum siesta_md_anneal_type { TEMPERATURE_ANNEAL, PRESSURE_ANNEAL, TEMP_AND_PRESS_ANNEAL } md_anneal_option;
//MD.AnnealOption

gdouble md_tau_relax;                   //MD.TauRelax
gdouble md_bulk_modulus;                //MD.BulkModulus

gdouble md_fc_displ;                    //MD.FCDispl
gint md_fc_first;                       //MD.FCfirst
gint md_fc_last;                        //MD.FClast
//TOO HARD
//MD.ATforPhonon


/*
 * FILE I/O
 */

gboolean use_saved_data;
gboolean long_output;
GtkWidget * long_output_widget;
gboolean write_density_matrix;

gboolean file_output_write_coor_init;       //WriteCoorInitial
gboolean file_output_write_coor_step;       //WriteCoorStep
gboolean file_output_write_forces;          //WriteForces
gboolean file_output_write_kpoints;         //WriteKpoints
gboolean file_output_write_eigenvalues;     //WriteEigenvalues
gboolean file_output_write_kbands;          //WriteKbands
gboolean file_output_write_bands;           //WriteBands
gboolean file_output_write_wavefunctions;   //WriteWaveFunctions
gdouble  file_output_write_mullikenpop;     //WriteMullikenPop

gboolean file_output_write_dm;              //WriteDM
gboolean file_output_write_coor_xmol;       //WriteCoorXmol
gboolean file_output_write_coor_cerius;     //WriteCoorCerius
gboolean file_output_write_md_xmol;         //WriteMDXmol
gboolean file_output_write_md_history;      //WriteMDhistory


//Animate frequencies.
gboolean show_vect;
GtkWidget * siesta_phonon_slider;
GtkWidget * freq_text_box;
gchar * freq_disp_str;
gchar * current_frequency_string;
gint current_animation;
gint num_animations;
gdouble current_frequency;
gboolean toggle_animate;
gboolean vibration_calc_complete;
gboolean animation_paused;
gdouble custom_scaler;

/*
MAT * eigen_xyz_atom_mat;
VEC * eigen_values;
*/
gpointer eigen_xyz_atom_mat;
gpointer eigen_values;


gint * sorted_eig_values;

gdouble atoms_per_job;


gdouble density_of_states;
gdouble density_on_mesh;
gdouble electrostatic_pot_on_mesh;

gchar * modelfilename;
};

/***************/
/* Monty2 data */
/***************/
struct monty_pak
{
  /* crystal graph images */
  gdouble image_x;
  gdouble image_y;
  gdouble image_z;
  
  /* &INPUT */
  gchar *hkls;
  gchar *supersaturations;
  gchar *input_surface;
  gchar *input_dir;
  gchar *input_cgf;
  gchar *energy_unit;
  gchar *esolv;
  gboolean run_nucleation;
  gboolean run_diffusion;
  gdouble nucleus_size;
  
  /* &OUTPUT */
  gchar *output_dirs;
  gchar *output_extension;
  gboolean write_surface;
  gboolean write_xyz;    /* -xmm */
  gboolean write_matlab; /* -xyz */
  gboolean write_msi;
  gboolean write_cube;
  
  /* &MODEL */
  gboolean spiral;
  gdouble xsteps;
  gdouble ysteps;
  gdouble kinetics;
  gdouble temperature;
  
  /* &MONITOR */
  gboolean monitor_height;
  gboolean monitor_energy;
  gboolean monitor_hhcorr;
  gboolean multi_frame_xyz;
  gboolean monitor_diffusion_profile;
  
  /* &RUN */
  gchar * random_seed;
  gdouble rows;
  gdouble cols;
  gdouble layers;
  gdouble increment;
  gdouble relax;
  gdouble cycles;
  gdouble moves;  
};

/***************/
/* abinit data */
/***************/
struct abinit_pak
{
gdouble energy;
gdouble max_grad;
gdouble rms_grad;
};

/***************/
/* NWChem data */
/***************/
struct nwchem_pak
{
gdouble energy;
gdouble max_grad;
gdouble rms_grad;
gint min_ok;
gint have_energy;
gint have_max_grad;
gint have_rms_grad;
};

/***************/
/* CASTEP data */
/***************/
struct castep_pak
{
  gdouble energy;
  gdouble max_grad;
  gdouble rms_grad;
  gint min_ok;
  gint have_energy;
  gint have_max_grad;
  gint have_rms_grad;
};

/*************************/
/* Quantum ESPRESSO data */
/*************************/
struct qe_species_pak
{
  gchar *label;
  gdouble mass;
  gchar *pseudo;
};

struct qe_pak
{
  gdouble energy;
  gdouble max_grad;
  gdouble rms_grad;
  gint min_ok;
  gint have_energy;
  gint have_max_grad;
  gint have_rms_grad;
  gdouble have_species;
  GSList *species;
  
  /* SYSTEM namelist */
  gchar *input_dft;
  gdouble ecutwfc;
  GSList *unparsed_system;
};

/**********************/
/* surface generation */
/**********************/
struct surface_pak
{
struct model_pak *model;
gint optimise;         /* relaxed/unrelaxed */
gint converge_eatt;    /* T - eatt, F - esurf */
gint converge_r1;
gint converge_r2;
gint include_polar;    /* calc shifts - keep polar surfaces */
gint ignore_bonding;   /* gensurf - ignore bonding */
gint create_surface;
gint true_cell;
gint keep_atom_order;
gint ignore_symmetry;
gint bonds_full;     /* NEW - total bonds in source unit cell */ 
gint bonds_cut;      /* NEW - cleaved bonds in surface */ 
gdouble miller[3];
gdouble shift;
gdouble dspacing;
gdouble depth;
gdouble depth_vec[3];
gdouble region[2];
gdouble dipole_tolerance;
};

/************/
/* XRD data */
/************/
struct diffract_pak
{
gint radiation;
gdouble wavelength;
gdouble theta[3];
gdouble dhkl_min;
gint broadening;
gdouble asym;
gdouble u, v, w;
};

/************************/
/* model data structure */
/************************/
struct model_pak
{
/* main info */
gint id;
gboolean locked;         /* prevent deletion (eg by dependent child process) */
gboolean protein;        /* cope with ugly PDB element labelling */
gint mode;
gint state;              /* the mode state (always reset when changing modes) */
gint grafted;            /* indicates if it's been put on the model tree yet */
gint redraw;             /* indicates a redraw is required */
gint redraw_count;       /* current number of redraws */
gint redraw_current;     /* current redraw time (micro-sec) */
gint redraw_cumulative;  /* cumulative redraw time (micro-sec) */
gint redraw_time;        /* average redraw time (micro-sec) */
gint fractional;
guint periodic;
guint num_frames;
guint cur_frame;
/* external animation frame info */
guint header_size;
guint frame_size;
guint file_size;
guint expected_cores;
guint expected_shells;
gboolean trj_swap;
gboolean build_molecules;
gboolean build_hydrogen;
gboolean build_polyhedra;
gboolean build_zeolite;
gint coord_units;      /* currently for converting bohr */
gint construct_pbc;
gint done_pbonds;   
gint cbu;              /* continuous update of bonds (# atoms related) */
gint colour_scheme;
gint has_sof;
gint sof_colourize;
gint atom_info;
gint shell_info;
gint default_render_mode;
gdouble rmax;
gdouble zoom;

/* display flags for the drawing area */
gint show_names;
gint show_title;
gint show_frame_number;
gint show_charge;
gint show_atom_charges;
gint show_atom_labels;
gint show_atom_types;
gint show_atom_index;
gint show_geom_labels;
gint show_cores;
gint show_bonds;
gint show_hbonds;
gint show_links;
gint show_shells;
gint show_axes;
gint show_cell;
gint show_cell_images;
gint show_cell_lengths;
gint show_waypoints;
gint show_selection_labels;

/*VZ*/
gint show_nmr_shifts;
gint show_nmr_csa;
gint show_nmr_efg;

/* ghosts */
gint num_atoms;
gint num_asym;
gint num_shells;
gint num_ghosts;
gint num_geom;

/* main lists */
GSList *cores;
GSList *shels;
GSList *ghosts;
GSList *bonds;
GSList *ubonds;
GSList *moles;
GSList *images;
GSList *selection;
GSList *unique_atom_list;
GSList *elements;
GSList *ff_list;
GSList *measure_list;
GSList *layer_list;
GSList *graph_list;
GSList *picture_list;
GSList *transform_list;
GSList *frame_data_list;
GSList *waypoint_list;
GSList *undo_list;
gpointer graph_active;
gpointer picture_active;
gpointer project;
gpointer zone_array;
gpointer zmatrix;

/* property hash table and ranked list */
GHashTable *property_table;
GSList *property_list;

/* special object lists */
GSList *ribbons;
GSList *spatial;
GSList *phonons;
GSList *ir_list;
GSList *raman_list;
gint num_phonons;
gint current_phonon;
gint pulse_count;
gint pulse_direction;
gdouble pulse_max;
GtkWidget *phonon_slider;
gboolean show_eigenvectors;
gboolean phonon_movie;
gchar *phonon_movie_name;
gchar *phonon_movie_type;

/* electrostatic scale */
gint epot_autoscale;
gdouble epot_min;
gdouble epot_max;
gint epot_div;

/* camera */
gpointer camera;
gpointer camera_default;

/* model geometric and display data */
gdouble offset[3];
gdouble pbc[6];
gdouble centroid[3];
gdouble cell_angles[3];
gdouble area;
gdouble volume;
/* translation matrix/vectors & its inverse */
gdouble latmat[9];
gdouble ilatmat[9];
/* reciprocal lattice vectors */
gdouble rlatmat[9];
/* rotation/translation transformation matrix */
gdouble display_lattice[16];
/* special object flags */
gint axes_type;
gint mesh_on, box_on, asym_on;
gboolean ms_colour_scale;
gint ms_colour_method;

/* special object data */
struct vec_pak xlat[6];
struct vec_pak axes[6];
struct vec_pak cell[8];
gint select_box[4];

/* morphology data */
gint morph_type;                /* BFDH, etc. */
gint morph_label;               /* flag - hkl labelling */
gint num_vertices;
gint num_planes;                /* total number of faces */
gint num_facets;                /* number of present facets in plane list */
gint sculpt_shift_use;          /* attempt to tailor sculpted morphology to shifts */
/* morphology */
GSList *planes;

/* space group */
struct space_pak sginfo;
/* symmetry info */
struct symmetry_pak symmetry;
/* periodic image numbers - start/stop multiples */
gint num_images;
gdouble image_limit[6];

/* region display (marvin only) */
gint region_max;
gint region[4];
gint region_empty[4];

/* text info */
gchar *title;
gchar filename[FILELEN];
gchar *basename;
GString *error_file_read;

/* energetics */
struct gulp_pak gulp;
struct gamess_pak gamess;
struct siesta_pak siesta;
struct abinit_pak abinit;
struct nwchem_pak nwchem;
struct castep_pak castep;
struct qe_pak qe;
struct diffract_pak diffract;
/* surface config */
struct surface_pak surface;
/* NEW: Monty settings */
struct monty_pak monty;

/* element data */
/*
gint num_elem;
struct elem_pak *elements;
*/
/* animation sequences */
gint animation;
gint animating;
gint anim_confine;
gboolean anim_fix;
gboolean anim_noscale;
gboolean anim_loop;
gdouble anim_speed;
gdouble anim_step;
gpointer analysis;
GList *frame_list;

/* current frame data */
gdouble frame_time;
gdouble frame_ke;
gdouble frame_pe;
gdouble frame_temp;

FILE *afp;
};

/**********************/
/* Supported file pak */
/**********************/
struct file_pak
{
gint id;                      /* unique identifier */
gint group;                   /* used to group inp/out types */
gint menu;                    /* (true/false) include in file menu listing? */
gchar *label;                 /* text info for the user */
GSList *ext;                  /* extension matching list */
gint (*read_file) (gchar *, struct model_pak *);     /* read routine */
gint (*write_file) (gchar *, struct model_pak *);    /* write routine */
gint (*read_frame) (FILE *, struct model_pak *);     /* animation read */
};

/*******************/
/* grid scheduling */
/*******************/
struct grid_config_pak
{
    gint id;                //unique ID
    gboolean active;        //just incase
    gchar * grid_name;      //a name for the grid
    gdouble max_cpus;       //MAX cpus
    gdouble max_memory;     //MAX memory
    gdouble max_walltime;   //walltime request
    gdouble max_jobfs;      //walltime request
    
    gdouble num_cpus;       //default number of CPUS
    gdouble memory;         //vmem request
    gdouble walltime;       //walltime request
    gdouble jobfs;          //file system request
    gchar * username;
    gchar * project_code;   //
    gchar * grid_address;   //hostname
    
    enum grid_access_method_type { SSH_ACCESS, MAGIC_ACCESS } access_method;
    
    gchar * queueing_system;
};

/*********************************/
/* widget convenience structures */
/*********************************/
struct callback_pak
{
gchar *label;
gint id;
};

/*******************/
/* main prototypes */
/*******************/
struct model_pak *model_ptr(gint, gint);

void sys_free(void);

void get_eledata(void);
void put_elem_data(struct elem_pak *, struct model_pak *);
gint get_elem_data(gint, struct elem_pak *, struct model_pak *);

void make_axes(struct model_pak *);
void make_cell(struct model_pak *);
void make_mesh(struct model_pak *);

void gui_render_dialog(void);
void gui_mdi_dialog(void);
void mdi_model_create(gint);

void write_dat(struct model_pak *);

gdouble ran2(void);
gint fill(void);

void update_box(gint, gint, struct model_pak *, gint);

gint elem_symbol_test(const gchar *);
gint elem_number_test(const gchar *);
gint elem_test(const gchar *);

gint pdb_elem_type(const gchar *);

void set_rgb_colour(GdkGC *, gint *);

gint write_gdisrc(void);
gint read_gdisrc(void);

gint event_render_modify(GtkWidget *, gpointer *); 

void delay(gint);

gint read_elem_data(FILE *, gint);
gint write_elem_data(FILE *);

void cmd_init(gint, gchar **);

GSList *slist_gchar_dup(GSList *);
void free_slist(GSList *);
void free_list(GList *);

