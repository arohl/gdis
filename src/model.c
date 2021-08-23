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
#include <stdlib.h>
#include <time.h>
#include <unistd.h>
#include <sys/stat.h>

#include "gdis.h"
#include "coords.h"
#include "matrix.h"
#include "edit.h"
#include "error.h"
#include "file.h"
#include "graph.h"
#include "morph.h"
#include "model.h"
#include "measure.h"
#include "project.h"
#include "analysis.h"
#include "render.h"
#include "opengl.h"
#include "select.h"
#include "space.h"
#include "surface.h"
#include "spatial.h"
#include "interface.h"
#include "dialog.h"
#include "zone.h"
#include "zmatrix.h"
#include "undo.h"

/* the main pak structures */
extern struct sysenv_pak sysenv;
extern struct elem_pak elements[];
/* main model database */
struct model_pak *models[MAX_MODELS];

/***********************************/
/* unified initialization/cleaning */
/***********************************/
void model_init(struct model_pak *data)
{
static gint n=0;

/* NB: all the important defines should go here */
/* ie those values that might cause a core dump due to */
/* code that requires they at least have a value assigned  */
data->id = -1;
data->locked = FALSE;
data->protein = FALSE;
data->num_atoms = 0;
data->num_asym = 0;
data->num_shells = 0;
data->num_species = 0;
data->num_geom = 0;
data->num_images = 0;
data->num_frames = 1;
data->cur_frame = 0;
data->expected_cores = 0;
data->expected_shells = 0;
data->header_size = 0;
data->frame_size = 0;
data->file_size = 0;
data->trj_swap = FALSE;
data->selection = NULL;
data->atom_info = -1;
data->shell_info = -1;
data->sof_colourize = FALSE;
data->construct_pbc = FALSE;
data->redraw = FALSE;
data->redraw_count = 0;
data->redraw_cumulative = 0;
data->redraw_time = 0;
data->colour_scheme = ELEM;
data->rmax = RMAX_FUDGE;
data->silent=FALSE;
data->track_me=FALSE;
data->t_next=NULL;
data->track_nb=0;
data->snapshot_eps=FALSE;
data->eps_file=NULL;

/* linked list init */
data->cbu = TRUE;
data->cores = NULL;
data->shels = NULL;
data->bonds = NULL;
data->ubonds = NULL;
data->moles = NULL;
data->planes = NULL;
data->ghosts = NULL;
data->images = NULL;
data->camera = NULL;
data->camera_default = NULL;
data->ff_list = NULL;
data->layer_list = NULL;
data->measure_list = NULL;
data->graph_list = NULL;
data->picture_list = NULL;
data->transform_list = NULL;
data->waypoint_list = NULL;
data->frame_data_list = NULL;
data->graph_active = NULL;
data->graph_ui_active = NULL;
data->picture_active = NULL;
data->project = NULL;
data->zone_array = NULL;
data->zmatrix = zmat_new();

/* display flags for the drawing area */
data->show_names = FALSE;
data->show_title = TRUE;
data->show_frame_number = TRUE;
data->show_atom_charges = FALSE;
data->show_atom_labels = FALSE;
data->show_atom_types = FALSE;
data->show_atom_index = FALSE;
data->show_geom_labels = TRUE;
data->show_cores = TRUE;
data->show_core_charges = FALSE;
data->show_shell_charges = FALSE;
data->show_shells = FALSE;
data->show_bonds = TRUE;
data->show_hbonds = FALSE;
data->show_links = FALSE;
data->show_axes = TRUE;
data->show_cell = TRUE;
data->show_cell_images = TRUE;
data->show_cell_lengths = FALSE;
data->show_waypoints = TRUE;
data->show_selection_labels = FALSE;

/*VZ*/
data->show_nmr_shifts = FALSE;
data->show_nmr_csa = FALSE;
data->show_nmr_efg = FALSE;

/* NEW */
data->property_table = g_hash_table_new_full(g_str_hash, g_str_equal, g_free, property_free);
data->property_list = NULL;

undo_init(data);

data->ribbons = NULL;
data->spatial = NULL;
data->phonons = NULL;
data->ir_list = NULL;
data->raman_list = NULL;
data->phonon_slider = NULL;
data->current_phonon = -1;
data->pulse_count = 0;
data->pulse_direction = 1;
data->pulse_max = 10.0;
data->num_phonons = 0;
data->show_eigenvectors = FALSE;
data->phonon_movie = FALSE;
data->phonon_movie_name = g_strdup("phonon_movie");
data->phonon_movie_type = g_strdup("gif");
data->gulp.phonon = FALSE;
data->gulp.kpoints = NULL;

/* VdW surface calc */
data->ms_colour_scale = FALSE;

data->epot_autoscale = TRUE;
data->epot_min = 999999999.9;
data->epot_max = -999999999.9;
data->epot_div = 11;

/* init flags */
data->periodic = 0;
data->grafted = FALSE;
data->fractional = FALSE;
data->done_pbonds = FALSE;
data->axes_type = CARTESIAN;
data->box_on = FALSE;
data->asym_on = FALSE;
data->coord_units = CARTESIAN;
data->default_render_mode = BALL_STICK;

/* bonding modes */
data->build_molecules = TRUE;
data->build_hydrogen = FALSE;
data->build_polyhedra = FALSE;
data->build_zeolite = FALSE;

/* symmetry */
space_init(&data->sginfo);
data->symmetry.num_symops = 0;
data->symmetry.symops = g_malloc(sizeof(struct symop_pak));
data->symmetry.num_items = 1;
data->symmetry.items = g_malloc(2*sizeof(gchar *));
*(data->symmetry.items) = g_strdup("none");
*((data->symmetry.items)+1) = NULL;
data->symmetry.pg_entry = NULL;
data->symmetry.summary = NULL;
data->symmetry.pg_name = g_strdup("unknown");
VEC3SET(&data->pbc[0], 1.0, 1.0, 1.0);
VEC3SET(&data->pbc[3], 0.5*G_PI, 0.5*G_PI, 0.5*G_PI);

/* periodic images */
data->image_limit[0] = 0.0;
data->image_limit[1] = 1.0;
data->image_limit[2] = 0.0;
data->image_limit[3] = 1.0;
data->image_limit[4] = 0.0;
data->image_limit[5] = 1.0;

data->region_max = 0;
/* existence */
data->region_empty[REGION1A] = FALSE;
data->region_empty[REGION2A] = TRUE;
data->region_empty[REGION1B] = TRUE;
data->region_empty[REGION2B] = TRUE;
/* lighting */
data->region[REGION1A] = TRUE;
data->region[REGION2A] = FALSE;
data->region[REGION1B] = FALSE;
data->region[REGION2B] = FALSE;
/* morphology */
data->morph_label = FALSE;/*FIX valgrind complaint --OVHPA*/
data->morph_type = DHKL;
data->num_vertices = 0;
data->num_planes = 0;
data->sculpt_shift_use = FALSE;

data->basename = g_strdup_printf("model_%d", g_slist_index(sysenv.mal, data));
strcpy(data->filename, data->basename);
n++;

/* gulp template init */
/* NEW - COSMO init ... TODO - include all GULP inits */
gulp_init(&data->gulp);

/* GAMESS defaults */
gamess_init(&data->gamess);

data->abinit.energy = data->abinit.max_grad = data->abinit.rms_grad = 0.0;

data->nwchem.energy = data->nwchem.max_grad = data->nwchem.rms_grad = 0.0;
data->nwchem.have_energy = data->nwchem.have_max_grad = data->nwchem.have_rms_grad = data->nwchem.min_ok = FALSE;

data->castep.energy = data->castep.max_grad = data->castep.rms_grad = 0.0;
data->castep.have_energy = data->castep.have_max_grad = data->castep.have_rms_grad = data->castep.min_ok = FALSE;
  
data->qe.energy = data->qe.max_grad = data->qe.rms_grad = 0.0;
data->qe.have_energy = data->qe.have_max_grad = data->qe.have_rms_grad = data->qe.min_ok = FALSE;
data->qe.have_species = FALSE;
data->qe.species = NULL;
  
data->qe.input_dft = "";
data->qe.ecutwfc = 75.0;
data->qe.unparsed_system = NULL;

/* SIESTA defaults */
data->siesta.num_atoms = 0;
data->siesta.energy = data->siesta.max_grad = 0.0;
data->siesta.have_energy = data->siesta.have_max_grad = 0.0;
data->siesta.epot_min = 0.0;
data->siesta.epot_max = 0.0;
data->siesta.eden_scale = 1.0;
data->siesta.eden = NULL;
data->siesta.edos = NULL;
data->siesta.epot = NULL;
data->siesta.freq_disp_str = NULL;
data->siesta.eigen_values = NULL;
data->siesta.eigen_xyz_atom_mat = NULL;


/* SIESTA - TerryRankine - terryrankine */
data->siesta.split_zeta_norm = 0.15;         //double
data->siesta.energy_shift = 0.02;            //double
data->siesta.energy_shift_units = RY_ENRG;
data->siesta.mesh_cutoff = 100.0;            //double
data->siesta.mesh_cutoff_units = RY_ENRG;

data->siesta.xc_author_type = PZ_XCAUTH;
data->siesta.xc_functional_type = LDA_XCFUNC;


data->siesta.electronic_temperature = 300.0;  //double
data->siesta.electronic_temperature_units = K_TEMP;

data->siesta.spin_polarised = FALSE;


data->siesta.basis_type = SPLIT_PAOBT;
data->siesta.basis_set = DZP_ZETA;

data->siesta.custom_zeta = 3.0;
data->siesta.custom_zeta_polarisation= 1.0;

data->siesta.harris_functional = FALSE;


data->siesta.is_periodic = TRUE;
data->siesta.kgrid_cutoff = 0.0;               //double
data->siesta.kgrid_cutoff_units = BOHR_LEN;

data->siesta.eta_value = 0.0;                  //double
data->siesta.no_of_itterations = 1000.0;
data->siesta.localisation_radius = 8.0;        //double

data->siesta.no_of_cycles = 50.0;              //MaxSCFIterations
data->siesta.mixing_weight = 0.25;             //double
data->siesta.pulay_mixing = FALSE;
data->siesta.no_of_pulay_matrices = 0.0;       //double
data->siesta.diag_divide_and_conquer = TRUE;   //Diag.DivideAndConquer



data->siesta.run_type = SINGLE_POINT;

data->siesta.number_of_steps =0.0;
data->siesta.constant_component = CONSTANT_PRESSURE;

data->siesta.md_inital_temperature = 0.0;
data->siesta.md_inital_temperature_units = K_TEMP;
data->siesta.md_target_temperature = 0.0;
data->siesta.md_target_temperature_units = K_TEMP;



data->siesta.therostat_parameter = 1.0;

data->siesta.finite_diff_step_size = 1.0;

//Termination conditions
data->siesta.md_max_cg_displacement = 0.2;   //MD.MaxCGDispl
data->siesta.md_max_cg_displacement_units = BOHR_LEN;
data->siesta.md_max_force_tol = 0.04;        //MD.MaxForceTol
data->siesta.md_max_force_tol_units = EVANG_FORCE;
data->siesta.md_max_stress_tol = 1.0;        //MD.MaxStressTol
data->siesta.md_max_stress_tol_units = GPA_PRES;

//siesta MD options
data->siesta.md_type_of_run = VERLET_MDRUN;

data->siesta.md_variable_cell = FALSE;
data->siesta.md_precondition_variable_cell = 5.0;
data->siesta.md_precondition_variable_cell_units = ANG_LEN;


data->siesta.md_inital_time_step = 1;
data->siesta.md_final_time_step = 1;
data->siesta.timestep = 1.0;
data->siesta.timestep_units = FS_TIME;


data->siesta.md_target_temperature = 0.0;
data->siesta.md_target_temperature_units = K_TEMP;

data->siesta.md_quench = FALSE;

data->siesta.md_target_pressure = 0.0;
data->siesta.md_target_pressure_units = GPA_PRES;
data->siesta.md_target_stress_xx = -1.0;
data->siesta.md_target_stress_yy = -1.0;
data->siesta.md_target_stress_zz = -1.0;
data->siesta.md_target_stress_xy = 0.0;
data->siesta.md_target_stress_xz = 0.0;
data->siesta.md_target_stress_yz = 0.0;

data->siesta.md_nose_mass = 100;
data->siesta.md_nose_mass_units = RYFS_MOMINERT;

data->siesta.md_parrinello_rahman_mass = 100;
data->siesta.md_parrinello_rahman_mass_units = RYFS_MOMINERT;

data->siesta.md_anneal_option = TEMP_AND_PRESS_ANNEAL;

data->siesta.md_tau_relax = 100;
data->siesta.md_tau_relax_units = FS_TIME;

data->siesta.md_bulk_modulus = 100.0;
data->siesta.md_bulk_modulus_units = RYBOHR3_PRES;

data->siesta.md_fc_displ = 0.04;
data->siesta.md_fc_displ_units = BOHR_LEN;

data->siesta.md_fc_first = 1.0;
data->siesta.md_fc_last = data->num_atoms;

data->siesta.pressure = 0.0;                   //double

data->siesta.use_saved_data = FALSE;
data->siesta.long_output = FALSE;
data->siesta.write_density_matrix = FALSE;

data->siesta.density_of_states = 1.0;          //ASK
data->siesta.density_on_mesh = 1.0;            //ASK
data->siesta.electrostatic_pot_on_mesh = 0.0;  //ASK

data->siesta.file_output_write_coor_init = TRUE;        //WriteCoorInital
data->siesta.file_output_write_coor_step = FALSE;       //WriteCoorStep
data->siesta.file_output_write_forces = FALSE;          //WriteForces
data->siesta.file_output_write_kpoints = FALSE;         //WriteKpoints
data->siesta.file_output_write_eigenvalues = FALSE;     //WriteEigenvalues
data->siesta.file_output_write_kbands = FALSE;          //WriteKbands
data->siesta.file_output_write_bands = FALSE;           //WriteBands
data->siesta.file_output_write_wavefunctions = FALSE;   //WriteWaveFunctions
data->siesta.file_output_write_mullikenpop = 0.0;       //WriteMullikenPop

data->siesta.file_output_write_dm = TRUE;               //WriteDM
data->siesta.file_output_write_coor_xmol = FALSE;       //WriteCoorXmol
data->siesta.file_output_write_coor_cerius = FALSE;     //WriteCoorCerius
data->siesta.file_output_write_md_xmol = FALSE;         //WriteMDXmol
data->siesta.file_output_write_md_history = TRUE;       //WriteMDhistory


data->siesta.vibration_calc_complete = FALSE;
data->siesta.custom_scaler = 0.3;

data->siesta.modelfilename = g_strdup(data->filename);

/* end TerryRankine */



/* Monty defaults */

/* crystal graph images */
data->monty.image_x = 1.0;
data->monty.image_y = 1.0;
data->monty.image_z = 1.0;

/* Monty &INPUT */
data->monty.hkls = g_strdup("[0 0 1]");
data->monty.supersaturations = g_strdup("10.0 0.0 1.0");
data->monty.input_surface = g_strdup("");
data->monty.input_dir = g_strdup_printf(".%s", DIR_SEP);
data->monty.input_cgf = g_strdup("");
data->monty.energy_unit = g_strdup("kcal/mol");
data->monty.esolv = g_strdup("0.0");
data->monty.run_diffusion = FALSE;
data->monty.run_nucleation = FALSE;
data->monty.nucleus_size = 50.0;


/* Monty &OUTPUT */
data->monty.output_dirs = g_strdup_printf(".%s001%s", DIR_SEP, DIR_SEP);
data->monty.output_extension = g_strdup(".monty2");
data->monty.write_surface = TRUE;
data->monty.write_xyz = FALSE;    /* -xmm */
data->monty.write_matlab = FALSE; /* -xyz */
data->monty.write_msi = FALSE;
data->monty.write_cube = FALSE;

/* Monty &MODEL */
data->monty.spiral = FALSE;
data->monty.xsteps = 0.0;
data->monty.ysteps = 0.0;
data->monty.kinetics = 0.0;
data->monty.temperature = 300.0;

/* Monty &MONITOR */
data->monty.monitor_height = TRUE;
data->monty.monitor_energy = FALSE;
data->monty.monitor_hhcorr = FALSE;
data->monty.multi_frame_xyz = FALSE;
data->monty.monitor_diffusion_profile = FALSE;

/* Monty &RUN */
data->monty.random_seed = g_strdup("0");
data->monty.rows = 40.0;
data->monty.cols = 40.0;
data->monty.layers = 20.0;
data->monty.increment = 5.0;
data->monty.relax = 100000.0;
data->monty.cycles = 100.0;
data->monty.moves = 10000.0;

/* END Monty setup */

/* diffraction defaults */
VEC3SET(data->diffract.theta, 0.0, 90.0, 0.1);
data->diffract.wavelength = 1.54180;
data->diffract.asym = 0.18;
data->diffract.u = 0.0;
data->diffract.v = 0.0;
data->diffract.w = 0.0;

/* surface creation init */
data->surface.optimise = FALSE;
data->surface.converge_eatt = FALSE;
data->surface.converge_r1 = TRUE;
data->surface.converge_r2 = TRUE;
data->surface.include_polar = FALSE;
data->surface.ignore_bonding = FALSE;
data->surface.create_surface = FALSE;
data->surface.true_cell = FALSE;
data->surface.keep_atom_order = FALSE;
data->surface.ignore_symmetry = FALSE;
data->surface.bonds_full = 0;
data->surface.bonds_cut = 0;
data->surface.miller[0] = 1;
data->surface.miller[1] = 0;
data->surface.miller[2] = 0;
data->surface.shift = 0.0;
data->surface.dspacing = 0.0;
data->surface.region[0] = 1;
data->surface.region[1] = 1;
data->surface.dipole_tolerance=0.005;
VEC3SET(data->surface.depth_vec, 0.0, 0.0, 0.0);

matrix_identity(data->latmat);
matrix_identity(data->ilatmat);
matrix_identity(data->rlatmat);

/* plots */
data->plots=FALSE;
data->plot=NULL;

/* density of states (DOS) */
data->spin_polarized=FALSE;
data->efermi=0.0;
data->ndos=0;
data->dos_eval=NULL;
data->dos_spin_up=NULL;
data->dos_spin_down=NULL;

/* bandstructure */
data->nbands=0;
data->nkpoints=0;
data->kpts_d=NULL;
data->band_up=NULL;
data->band_down=NULL;

/* frequency */
data->have_frequency = FALSE;
data->nfreq=0;
data->freq=NULL;
data->freq_intens=NULL;

/* uspex */
data->uspex=NULL;
/* vasp */
data->vasp=NULL;

/* no non-default element data for the model (yet) */
data->elements = NULL;
data->unique_atom_list = NULL;
/* file pointer */
data->animation = FALSE;
data->animating = FALSE;
data->anim_fix = FALSE;
data->anim_noscale = FALSE;
data->anim_loop = FALSE;
data->afp = NULL;
data->anim_confine = PBC_CONFINE_ATOMS;
data->anim_speed = 5.0;
data->anim_step = 1.0;
data->frame_list = NULL;
data->title = NULL;
data->analysis = analysis_new();

/* NEW */
data->error_file_read = g_string_new(NULL);
}

/*********************/
/* duplicate a model */
/*********************/
/* NB: this is primarily meant for threads to make a copy that can be */
/* manipulated without fear of deletion - consequently, much of the display */
/* related data does not need to be transferred (mainly want the coordinates) */
gpointer model_dup(struct model_pak *source)
{
struct model_pak *model;

/* checks */
g_assert(source != NULL);
if (source->num_frames > 1)
  {
printf("Not properly dealt with yet...\n");
/* TODO - if points to a file - ok, since the file (hopefully won't be deleted) */
  return(NULL);
  }

/* allocation */
model = g_malloc(sizeof(struct model_pak));
model_init(model);

/* duplication */
model->cores = dup_core_list(source->cores);
model->shels = dup_shell_list(source->shels);

/* TODO - have a lattice pointer that can be memcpy'd in one hit */
model->periodic = source->periodic;
model->fractional = source->fractional;

memcpy(model->latmat, source->latmat, 9*sizeof(gdouble));
memcpy(model->ilatmat, source->ilatmat, 9*sizeof(gdouble));
memcpy(model->rlatmat, source->rlatmat, 9*sizeof(gdouble));
memcpy(model->pbc, source->pbc, 6*sizeof(gdouble));

g_free(model->sginfo.spacename);
model->sginfo.spacename = source->sginfo.spacename;
space_lookup(model);

return(model);
}

/************************************/
/* standard preparation for display */
/************************************/
#define DEBUG_PREP_MODEL 0
gint model_prep(struct model_pak *data)
{
g_return_val_if_fail(data != NULL, 1);

#if DEBUG_PREP_MODEL
printf("start prep: %p\n", data);
#endif

error_table_clear();

/* init the original values */
/* FIXME - what if there are asymmetric duplicates? */
data->num_asym = g_slist_length(data->cores);

/* create the all important lattice matrix before anything else */
matrix_lattice_init(data);

/* convert to angstroms if required */
coords_init_units(data);

/* convert input cartesian coords to fractional */
if (!data->fractional)
  coords_make_fractional(data);

/* initialize the spatial partitioning */
zone_init(data);

/* very large model exception test */
if (data->num_asym < 100000)
  {
/* remove any repeats */
  delete_duplicate_cores(data);
  }
else
  {
  data->build_molecules = FALSE;
  }

/* NB: core-shell links are required before doing the space group stuff */
shell_make_links(data);

/* mark the largest region label */
data->region_max = region_max(data);

/* space group lookup/full cell generation */
if (data->periodic == 3)
  {
  if (!data->sginfo.spacenum)
    data->sginfo.spacenum=1;
  if (space_lookup(data))
    gui_text_show(ERROR, "Error in Space Group lookup.\n");

  data->axes_type = OTHER;
  }

/* multi-frame file is always an animation */
if (data->num_frames > 1)
  data->animation = TRUE;

/* calculate how many frames an associated gulp trg file has */
if (data->animation)
  {
  guint n, header_size;
  gdouble frame_size, num_frames;
  gchar *filename;
  struct stat buff;
  FILE *fp;
#define FRAME_EPS 0.001

/* calculate frame size for num_frame calculation */
  if (data->id == GULP && !data->file_size)
    {
    header_size = 6*sizeof(int) + sizeof(double);
/* sizeof data */
    data->expected_cores = g_slist_length(data->cores);
    data->expected_shells = g_slist_length(data->shels);
    n = data->expected_cores + data->expected_shells;
/* cope with the supercell cell multiplication */
    n *= data->gulp.super[0] * data->gulp.super[1] * data->gulp.super[2];

/* NB: a bit messy, but the supercell isn't created until later */
/* (see file_load() - after the entire file has been processed) */

/* sizeof RECORD */
    frame_size = 14 * sizeof(int);
    frame_size += (4 + 6*n) * sizeof(double);

if (data->gulp.ensemble == NPT)
  {
  frame_size += 4 * sizeof(int);

/* nstrains (cell velocities) */
  switch (data->periodic)
    {
    case 3:
      n=6;
      break;
    case 2:
      n=3;
      break;
    case 1:
      n=1;
      break;
    default:
      n=0;
    }

/* cell vectors */
  n += 9;

  frame_size += n * sizeof(double);
  }

/* stat the file (MUST have the full path) */
    filename = g_strdup_printf("%s/%s", sysenv.cwd, data->gulp.trj_file);
/* get num frames (as a float - for debugging) */
    if (stat(filename, &buff))
      {
      printf("Warning: bad or missing GULP trajectory file.\n");
      data->animation = FALSE;
      }
    else
      {
/* attempt to read the header */
fp = fopen(filename, "r");
if (fp)
  read_trj_header(fp, data);
else
  printf("Failed to open %s.\n", filename);

      data->file_size = buff.st_size;
      num_frames = data->file_size - header_size;
      num_frames /= frame_size;
      data->num_frames = (num_frames + FRAME_EPS);
      data->header_size = header_size;
      data->frame_size = frame_size;

#if DEBUG_PREP_MODEL
printf("  num cores: %d\n", data->expected_cores);
printf(" num shells: %d\n", data->expected_shells);
printf("header size: %d\n", data->header_size);
printf(" frame size: %d\n", data->frame_size);
printf("  file size: %d\n", data->file_size);
printf(" num frames: %d (%f)\n", data->num_frames, num_frames);
#endif

if ((num_frames - data->num_frames) > FRAME_EPS)
  {
  printf("Bad GULP trajectory file.\n");
  }
fclose(fp);

/* NEW - attempt to mark the start of all frames (for faster random access) */
      mark_trj_frames(data);

      }
    g_free(filename);
    }
  }

/* init */
data->num_atoms = g_slist_length(data->cores);
data->num_shells = g_slist_length(data->shels);
data->mode = FREE;

#if DEBUG_PREP_MODEL
printf(" num atoms: %d\n", data->num_atoms);
printf("num shells: %d\n", data->num_shells);
#endif

/* initial center/coord calc */
coords_init(INIT_COORDS, data);

/* connectivity */
connect_bonds(data);
connect_molecules(data);

/* assign colours */
model_colour_scheme(data->colour_scheme, data);

/* final periodicity setup */
if (data->periodic == 2)
  {
/* order by z, unless input core order should be kept */
  if (!data->surface.keep_atom_order)
    sort_coords(data);
  }

/* final recenter/coord calc (after molecule unfragmenting/surface sorting) */
coords_init(CENT_COORDS, data);

/* set up the camera - should be done after all coord centering etc. */
if (sysenv.canvas)
  camera_init(data);

#if DEBUG_PREP_MODEL
printf("end prep: %p\n", data);
#endif

error_table_print_all();

return(0);
}

/*****************************************/
/* calculate density of 3D cell in g/cm3 */
/*****************************************/
gdouble model_density(struct model_pak *model)
{
gdouble density, volume;
GSList *list;
struct core_pak *core;

density = 0.0;
volume = model->volume;

if ( volume != 0.0 )
  {
  /* loop over cores */
  for (list=model->cores ; list ; list=g_slist_next(list))
    {
    core = list->data;
    if( core->status & (DELETED | HIDDEN))
      continue;
    density += atom_mass(core);
    }
  density /= (AVOGADRO*volume*1e-24);
  }

return(density);
}

/*************************************************/
/* store the current position of the file stream */
/*************************************************/
#define DEBUG_ADD_FRAME_OFFSET 0
gint add_frame_offset(FILE *fp, struct model_pak *data)
{
fpos_t *offset;

g_assert(fp != NULL);
g_assert(data != NULL);

offset = g_malloc(sizeof(fpos_t));
/* prepend & reverse? */
if (!fgetpos(fp, offset))
  data->frame_list = g_list_append(data->frame_list, offset); 
else
  {
  g_free(offset);
  return(1);
  }

#if DEBUG_ADD_FRAME_OFFSET
printf("stored: %p\n", offset);
#endif

return(0);
}

/*********************************************************/
/* convenience routine for duplicating a list of strings */
/*********************************************************/
GSList *slist_gchar_dup(GSList *src)
{
GSList *list, *dest=NULL;

for (list=src ; list ; list=g_slist_next(list))
  dest = g_slist_prepend(dest, g_strdup(list->data));

if (dest)
  dest = g_slist_reverse(dest);

return(dest);
}

/*************************************************************/
/* generic call for an slist containing one freeable pointer */
/*************************************************************/
void free_slist(GSList *list)
{
GSList *item;

if (!list)
  return;

for (item=list ; item ; item=g_slist_next(item))
  g_free(item->data);

g_slist_free(list);
list=NULL;
}

/***********************************************************/
/* generic call for a list containing one freeable pointer */
/***********************************************************/
void free_list(GList *list)
{
GList *item;

for (item=list ; item ; item=g_list_next(item))
  g_free(item->data);

g_list_free(list);
list=NULL;
}

/*************************************************/
/* TODO - put stuff like this in a gamess.c file */
/*************************************************/
void gamess_data_free(struct model_pak *model)
{
g_free(model->gamess.title);
g_free(model->gamess.temp_file);
g_free(model->gamess.out_file);
}

/***********************************/
/* free the entire model structure */
/***********************************/
#define DEBUG_FREE_MODEL 0
void model_free(struct model_pak *data)
{
/* check */
g_return_if_fail(data != NULL);

#if DEBUG_FREE_MODEL
printf("freeing string data...\n");
#endif

g_free(data->basename);
g_free(data->symmetry.pg_name);
g_free(data->title);

space_free(&data->sginfo);

/*free vasp pointer*/
if(data->vasp!=NULL){
	free_vasp_out(data->vasp);
	g_free(data->vasp);
	data->vasp=NULL;
}
/*free uspex pointer*/
if(data->uspex!=NULL){
	free_uspex_out(data->uspex);
	g_free(data->uspex);
	data->uspex=NULL;
}

#if DEBUG_FREE_MODEL
printf("freeing zone data...\n");
#endif

zone_free(data->zone_array);

#if DEBUG_FREE_MODEL
printf("freeing coord data...\n");
#endif

/* free lists with associated data */
free_core_list(data);
free_slist(data->images);
free_list(data->frame_list);
free_slist(data->layer_list);
free_slist(data->picture_list);
graph_free_list(data);
free_slist(data->phonons);
free_slist(data->ir_list);
free_slist(data->raman_list);
free_slist(data->waypoint_list);
free_slist(data->transform_list);
free_slist(data->ff_list);

/* camera */
g_free(data->camera_default);

/* destroy the property table and list */
g_hash_table_destroy(data->property_table);
g_slist_free(data->property_list);

/* external package configurations */
gamess_data_free(data);
gulp_files_free(data);
gulp_data_free(data);

#if DEBUG_FREE_MODEL
printf("freeing symmetry data...\n");
#endif
g_free(data->symmetry.symops);
g_strfreev(data->symmetry.items);

#if DEBUG_FREE_MODEL
printf("freeing measurements...\n");
#endif
meas_prune_model(data);
measure_free_all(data);

#if DEBUG_FREE_MODEL
printf("freeing vertices... (%d)\n", data->num_vertices);
#endif

plane_data_free(data->planes);
g_slist_free(data->planes);
data->planes = NULL;
data->num_planes = 0;

/* spatial objects */
spatial_destroy_all(data);

#if DEBUG_FREE_MODEL
printf("freeing element data...\n");
#endif
free_slist(data->elements);
data->elements = NULL;

/*
#if DEBUG_FREE_MODEL
printf("closing animation stream...\n");
#endif
if (data->afp)
  fclose(data->afp);
*/

analysis_free(data->analysis);
data->analysis = NULL;

g_string_free(data->error_file_read, TRUE);
}

/************************/
/* delete a given model */
/************************/
void model_delete(struct model_pak *model)
{
GSList *list;
struct canvas_pak *canvas;

/* checks */
if (!model)
  return;
if (model->locked)
  {
  gui_text_show(ERROR, "Model is locked.\n");
  return;
  }

if (!g_slist_find(sysenv.mal, model))
  {
  printf("WARNING: unregistered model.\n");
  }
else
  {
/* remove reference from the canvas list */
  for (list=sysenv.canvas_list ; list ; list=g_slist_next(list))
    {
    canvas = list->data;
    if (canvas->model == model)
      canvas->model = NULL;
    }

  sysenv.mal = g_slist_remove(sysenv.mal, model);
  }

#ifdef WITH_GUI
/* destroy any associated dialogs */
dialog_destroy_model(model);
#endif

/* free the model's pointers */
model_free(model);
sysenv.mal = g_slist_remove(sysenv.mal, model);
g_free(model);

/* update */
sysenv.refresh_dialog=TRUE;
canvas_shuffle();
redraw_canvas(ALL);
}

/**********************************/
/* update the model/canvas layout */
/**********************************/
/* TODO - rename - eg model_active_shuffle() */
void canvas_shuffle(void)
{
gint c, m;
GSList *clist, *mlist;
struct canvas_pak *canvas;

/* find active model in canvas list */
c = 0;
mlist = NULL;
for (clist=sysenv.canvas_list ; clist ; clist=g_slist_next(clist))
  {
  canvas = clist->data;
  if (sysenv.active_model)
    {

    if (!canvas->model &&!mlist)
      {
      canvas->model = sysenv.active_model;
      }

    if (canvas->model == sysenv.active_model)
      {
      m = g_slist_index(sysenv.mal, canvas->model);
      if (m < c)
        {
/* not enough prior models to display active model at current canvas poisiton */
/* start at active model */
/* TODO - start at low enough model number so that active model is displayed */
        mlist = g_slist_find(sysenv.mal, sysenv.active_model);
        }
      else
        {
/* we have enough prior models to display active model at current canvas poisiton */
        mlist = g_slist_nth(sysenv.mal, m-c);
        }
      }

    }
  c++;
  }

/* active model not in the canvas list */
if (sysenv.active_model && !mlist)
  {
/* start at active model */
  mlist = g_slist_find(sysenv.mal, sysenv.active_model);
  }

/* fill the canvas list */
if (mlist)
  {
  for (clist=sysenv.canvas_list ; clist ; clist=g_slist_next(clist))
    {
    canvas = clist->data;
    if (mlist)
      {
      canvas->model = mlist->data;
      mlist = g_slist_next(mlist);
      }
    else
      canvas->model = NULL;
    }
  }
}

/*******************************/
/* replacement for ASSIGN call */
/*******************************/
struct model_pak *model_new(void)
{
struct model_pak *model;

model = g_malloc(sizeof(struct model_pak));
sysenv.mal = g_slist_append(sysenv.mal, model);

model_init(model);
return(model);
}

/*******************************************/
/* get pointer to a requested model's data */
/*******************************************/
#define DEBUG_PTR 0
struct model_pak *model_ptr(gint model, gint mode)
{
struct model_pak *ptr=NULL;

/* how were we called? */
switch (mode)
  {
  case RECALL:
/* FIXME - there are an awful lot of calls to this; check for routines */
/* that pass the model number, when they could pass the model ptr instead */
#if DEBUG_PTR
printf("Request for model %d [%d,%d], ptr %p\n", model, 0, sysenv.num_models, ptr);
#endif
    ptr = (struct model_pak *) g_slist_nth_data(sysenv.mal, model);
    break;
  }
return(ptr);
}

/***********************************************/
/* model property list modification primitives */
/***********************************************/

struct property_pak
{
guint rank;
gchar *label;
/* _BUG_OVHPA_1
gchar *value;
 * The gchar *value; part is updated many times
 * during the reading of a file or file frames.
 * Each rewrite of its value involves using the
 * g_malloc (via g_strdup)/g_free couple which,
 * in a threaded environment, will fail.
 * This is due to the g_malloc "smart" re-alloc
 * which is using memory unreachable by current
 * thread.
 * FIX: use static allocation for value.
 */
gchar value[MAX_VALUE_SIZE];//_BUG_OVHPA_1
};

/*************************/
/* destruction primitive */
/*************************/
void property_free(gpointer ptr_property)
{
struct property_pak *p = ptr_property;

/* NB: the label/key is free'd elsewhere */
//g_free(p->value);//_BUG_OVHPA_1
g_free(p);
}

/*********************/
/* ranking primitive */
/*********************/
gint property_ranking(gpointer ptr_p1, gpointer ptr_p2)
{
struct property_pak *p1 = ptr_p1, *p2 = ptr_p2;

if (p1->rank < p2->rank)
  return(-1);
if (p1->rank > p2->rank)
  return(1);
return(0);
}

/***************************************/
/* add a pre-ranked new model property */
/***************************************/
/* NB: property is not displayed if rank = 0 */
void property_add_ranked(guint rank,
                         const gchar *key,
                         const gchar *value,
                         struct model_pak *model)
{
struct property_pak *p;

g_assert(model != NULL);

/* check if label already exists */
p = g_hash_table_lookup(model->property_table, key);
if (p)
  { //_BUG_OVHPA_1
#define DEBUG_PROPERTY 0
	if(value==NULL) return; /* refuse to update with a null value */
	if(value[0]=='\0') return; /* refuse to update with no value */
	if(g_strlcpy(p->value,value,MAX_VALUE_SIZE)>=MAX_VALUE_SIZE)
		fprintf(stderr,"_BUG_ the MAX_VALUE_SIZE (=%i) must be increased and code recompiled!\n",MAX_VALUE_SIZE);
	p->rank = rank;
	model->property_list = g_slist_remove(model->property_list, p);
#if DEBUG_PROPERTY
fprintf(stdout,"#DBG: update old %i: %s %s -> ",p->rank,p->label,p->value);
fprintf(stdout,"new %i: %s %s\n",rank,key,value);
#endif //DEBUG_PROPERTY
  }
else
  {
/* create new property */
  p = g_malloc(sizeof(struct property_pak));
  if(g_strlcpy(p->value,value,MAX_VALUE_SIZE)>=MAX_VALUE_SIZE)
		fprintf(stderr,"_BUG_ the MAX_VALUE_SIZE (=%i) must be increase and code recompiled!\n",MAX_VALUE_SIZE);
  p->label = g_strdup(key);
  p->rank = rank;
  g_hash_table_replace(model->property_table, p->label, p);
#if DEBUG_PROPERTY
fprintf(stdout,"#DBG: create new %i: %s %s\n",rank,key,value);
#endif
  }

/* update sorted property list */
model->property_list = g_slist_insert_sorted(model->property_list, p, (gpointer) property_ranking);
}

/****************************************/
/* refresh the model content properties */
/****************************************/
void model_content_refresh(struct model_pak *model)
{
gint n;
gdouble density;
gchar *text;

g_assert(model != NULL);

calc_emp(model);
if( model->periodic == 3 )
  {
  density = model_density(model);
  text = g_strdup_printf("%6.3f g/cm3", density);
  property_add_ranked(1, "Density", text, model);
  g_free(text);
  }
else
  property_add_ranked(0, "Density", "dummy", model);

if( model->periodic == 2 )
  {
  text = g_strdup_printf("%6.3f D", fabs(model->gulp.sdipole)<1e-3?0.0:model->gulp.sdipole);
  property_add_ranked(1, "Dipole (z)", text, model);
  g_free(text);
  }
else
  property_add_ranked(0, "Dipole (z)", "dummy", model);

text = g_strdup_printf("%6.3f e", fabs(model->gulp.qsum)<1e-3?0.0:model->gulp.qsum);
property_add_ranked(1, "Total charge", text, model);
g_free(text);

text = g_strdup_printf("%d", g_slist_length(model->moles));
property_add_ranked(1, "Total molecules", text, model);
g_free(text);

/* can be confusing as it also counts hydrogen bonds */
/*
text = g_strdup_printf("%d", g_slist_length(model->bonds));
property_add_ranked(1, "Total bonds", text, model);
g_free(text);
*/

n = g_slist_length(model->shels);
if (n)
  {
  text = g_strdup_printf("%d", n);
  property_add_ranked(1, "Total shells", text, model);
  g_free(text);
  }
else
  property_add_ranked(0, "Total shells", "dummy", model);

text = g_strdup_printf("%d", g_slist_length(model->cores));
property_add_ranked(1, "Total atoms", text, model);
g_free(text);
gui_refresh(GUI_MODEL_PROPERTIES);
}

/****************************************/
/* property value extraction primitives */
/****************************************/
gchar *property_lookup(gchar *key, struct model_pak *model)
{
struct property_pak *p;

g_assert(model != NULL);

p = g_hash_table_lookup(model->property_table, key);
if (p)
  return(g_strdup (p->value));//_BUG_OVHPA_1
return(NULL);
}

guint property_rank(gpointer ptr_property)
{
struct property_pak *p = ptr_property;

return(p->rank);
}

gchar *property_label(gpointer ptr_property)
{
struct property_pak *p = ptr_property;

return(p->label);
}

gchar *property_value(gpointer ptr_property)
{
struct property_pak *p = ptr_property;

return(p->value);
}

