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

#include "gdis.h"
#include "coords.h"
#include "error.h"
#include "file.h"
#include "parse.h"
#include "matrix.h"
#include "zmatrix.h"
#include "zmatrix_pak.h"
#include "model.h"
#include "interface.h"

gint set_true_false();
gint set_energy_units();
gint set_temperature_units();
gint set_length_units();
gint set_force_units();
gint set_pressure_units();
gint set_time_units();
gint set_mominert_units();

gint print_file_energy_units();
gint print_file_mass_units();
gint print_file_length_units();
gint print_file_time_units();
gint print_file_temperature_units();
gint print_file_pressure_units();
gint print_file_force_units();
gint print_file_mominert_units();

enum {SIESTA_DEFAULT, SIESTA_SPECIES, SIESTA_LATTICE_PARAMETERS, SIESTA_LATTICE_VECTORS,
      SIESTA_COORDS, SIESTA_ZMATRIX};

/* main structures */
extern struct sysenv_pak sysenv;
extern struct elem_pak elements[];


/**********************************/
/* obtain a core's species number */
/**********************************/
gint fdf_species_index(gchar *atom_label, GSList *species_list)
{
gint i;
struct species_pak *species;
GSList *list;

/* find corresponding number of this element in the species block */
i=1;
for (list=species_list ; list ; list=g_slist_next(list))
  {
  species = list->data;

  if (g_ascii_strcasecmp(atom_label, species->label) == 0)
    return(i);

  i++;
  }
return(0);
}

/*****************************************************/
/* generate the list of unique species, SIESTA style */
/*****************************************************/
GSList *fdf_species_build(struct model_pak *model)
{
gint code;
struct species_pak *species_data;
GSList *list, *normal_list, *ghost_list, *species_list;

/* init the normal, ghost, and overall species lists */
normal_list = find_unique(LABEL_NORMAL, model);
ghost_list = find_unique(LABEL_GHOST, model);
species_list = NULL;

for (list=normal_list ; list ; list=g_slist_next(list))
  {
  code = elem_symbol_test(list->data);

/* update the species identifier list */
  species_data = g_malloc(sizeof(struct species_pak));
  species_data->label = list->data;
  species_data->number = code;
  species_list = g_slist_prepend(species_list, species_data);
  }

for (list=ghost_list ; list ; list=g_slist_next(list))
  {
  code = elem_symbol_test(list->data);

/* update the species identifier list */
  species_data = g_malloc(sizeof(struct species_pak));
  species_data->label = list->data;
  species_data->number = -code;
  species_list = g_slist_prepend(species_list, species_data);
  }
species_list = g_slist_reverse(species_list);

/* parse the zmatrix entries and make the type match the (possible new) species order */
if (model->zmatrix)
  {
  struct zmat_pak *zmatrix = model->zmatrix;

/* scan zval entries - make sure type value reflects element order */
  for (list=zmatrix->zlines ; list ; list=g_slist_next(list))
    {
    struct zval_pak *zval = list->data;
    zval->type = fdf_species_index(zval->elem, species_list);
    }
  }

return(species_list);
}

/****************/
/* file writing */
/****************/
gint write_fdf(gchar *filename, struct model_pak *model)
{
gint i;
gdouble x[3], depth;
GSList *list, *clist, *species_list;
struct core_pak *core;
struct species_pak *species_data;
FILE *fp;

/* checks */
g_return_val_if_fail(model != NULL, 1);
g_return_val_if_fail(filename != NULL, 2);

/* open the file */
fp = fopen(filename,"wt");
if (!fp)
  return(3);

/* print header */
fprintf(fp, "# ");
gdis_blurb(fp);
fprintf(fp, "#\n\n");
fprintf(fp, "SystemLabel      %s\n\n", model->basename);
fprintf(fp, "NumberOfAtoms    %d\n\n", g_slist_length(model->cores));

/* init the normal, ghost, and overall species lists */
species_list = fdf_species_build(model);

fprintf(fp, "NumberOfSpecies  %d\n", g_slist_length(species_list));

/* write the species (unique atom types) block */
fprintf(fp, "%%block ChemicalSpeciesLabel\n");
i=1;
for (list=species_list ; list ; list=g_slist_next(list))
  {
  species_data = list->data;

  fprintf(fp, "  %3d  %3d  %s\n", i, species_data->number, species_data->label);
  i++;
  }
fprintf(fp, "%%endblock ChemicalSpeciesLabel\n\n");


/* write the lattice data */
if (model->periodic)
  {
  fprintf(fp, "LatticeConstant 1.0 Ang\n");
  if (model->construct_pbc)
    {
/* NB: siesta matrices are transposed wrt gdis */
    fprintf(fp, "%%block LatticeVectors\n");
    fprintf(fp,"%15.10f %15.10f %15.10f\n",
                model->latmat[0], model->latmat[3], model->latmat[6]);
    fprintf(fp,"%15.10f %15.10f %15.10f\n",
                model->latmat[1], model->latmat[4], model->latmat[7]);
    fprintf(fp,"%15.10f %15.10f %15.10f\n",
                model->latmat[2], model->latmat[5], model->latmat[8]);
    fprintf(fp, "%%endblock LatticeVectors\n\n");
    }
  else
    {
    fprintf(fp, "%%block LatticeParameters\n");
/* if it's a surface - print Dhkl (mult by region sizes) */
    if (model->periodic == 2)
      {
/* get depth info */
      depth = (model->surface.region[0]+model->surface.region[1])
                                          *model->surface.depth;
/* no depth info - make it large enough to fit everything */
      if (depth < POSITION_TOLERANCE)
        depth = 2.0*model->rmax;

      fprintf(fp, "  %f  %f  %f", model->pbc[0], model->pbc[1], depth);
      }
    else
      fprintf(fp, "  %f  %f  %f", model->pbc[0], model->pbc[1], model->pbc[2]);

    fprintf(fp, "  %f  %f  %f\n", R2D*model->pbc[3], R2D*model->pbc[4], R2D*model->pbc[5]);
    fprintf(fp, "%%endblock LatticeParameters\n\n");
    }
  }


/* write the zmatrix data */
/*
if (model->zmatrix)
*/
if (zmat_entries_get(model->zmatrix))
  {
  fprintf(fp, "ZM.UnitsLength %s\n", zmat_distance_units_get(model->zmatrix));
  fprintf(fp, "ZM.UnitsAngle %s\n\n", zmat_angle_units_get(model->zmatrix));

  fprintf(fp, "%%block Zmatrix\n");
 
  zmat_coord_print(fp, species_list, model);
 
  zmat_mol_print(fp, model->zmatrix);
 
  zmat_var_print(fp, model->zmatrix);

  zmat_const_print(fp, model->zmatrix);
 
  fprintf(fp, "%%endblock Zmatrix\n");
  }
else
  {

/* write the atoms - for surfaces (and isolated molecules), coords must be cartesian */
if (model->fractional)
  fprintf(fp, "AtomicCoordinatesFormat Fractional\n");
else
  fprintf(fp, "AtomicCoordinatesFormat NotScaledCartesianAng\n");

fprintf(fp, "%%block AtomicCoordinatesAndAtomicSpecies\n");

for (clist=model->cores ; clist ; clist=g_slist_next(clist))
  {
  core = clist->data;
  if (core->status & DELETED)
    continue;

  ARR3SET(x, core->x);

/* NB: want fractional if 3D periodic, otherwise cartesian */
  if (!model->fractional)
    vecmat(model->latmat, x);

/* find corresponding number of this element in the species block */
  i=1;
  for (list=species_list ; list ; list=g_slist_next(list))
    {
    species_data = list->data;

    if (g_ascii_strcasecmp(core->atom_label, species_data->label) == 0)
      {
      if (core->ghost)
        {
        if (species_data->number == -core->atom_code)
          break;
        }
      else
        {
        if (species_data->number == core->atom_code)
          break;
        }
      }
    i++;
    }

/* found? */
  if (list)
    fprintf(fp,"  %14.9f  %14.9f  %14.9f    %d\n", x[0], x[1], x[2], i);
  else
    printf("write_fdf() error: bad species type.\n");
  }
fprintf(fp, "%%endblock AtomicCoordinatesAndAtomicSpecies\n\n");

  }

free_slist(species_list);

/*
 * Siesta option writer
 */

//  PAO.BasisType
switch (model->siesta.basis_type)
{
    case SPLIT_PAOBT:
    {
        fprintf(fp, "PAO.BasisType    split\n");
        break;
    }
    case SPLITGAUSS_PAOBT:
    {
        fprintf(fp, "PAO.BasisType    splitgauss\n");
        break;
    }
    case NODES_PAOBT:
    {
        fprintf(fp, "PAO.BasisType    nodes\n");
        break;
    }
    case NONODES_PAOBT:
    {
        fprintf(fp, "PAO.BasisType    nonodes\n");
        break;
    }
    default:
    {
            //split!
        fprintf(fp, "PAO.BasisType    split\n");
        break;
    }
}



switch (model->siesta.basis_set)
{
    case SZ_ZETA :
    {
        fprintf(fp,"PAO.BasisSize    SZ\n");
        break;
    }
    case SZP_ZETA :
    {
        fprintf(fp,"PAO.BasisSize    SZP\n");
        break;
    }
    case DZ_ZETA :
    {
        fprintf(fp,"PAO.BasisSize    DZ\n");
        break;
    }
    case DZP_ZETA :
    {
        fprintf(fp,"PAO.BasisSize    DZP\n");
        break;
    }
    case CUSTOM_ZETA :
    {
        //TODO - CUSTOM Zeta Levels? - should set a tag somewhere
        fprintf(fp,"PAO.BasisSize    DZP\n");
        break;
    }
}


//siesta_solution_method_type
switch (model->siesta.solution_method_type)
{
    case DIAGON_SOLMETH: fprintf(fp, "SolutionMethod     diagon\n"); break;
    case ORDERN_SOLMETH: fprintf(fp, "SolutionMethod     OrderN\n"); break;
    default:
    //should never get here....
        error_table_entry("Inconsistancy model->siesta.solution_method_type.\n");
        fprintf(fp, "SolutionMethod     diagon\n");
        break;
}


//     if (model->siesta.split_zeta_norm != 0.15)
fprintf(fp, "PAO.SplitNorm    %f\n", model->siesta.split_zeta_norm);


fprintf(fp, "PAO.EnergyShift    %f ", model->siesta.energy_shift);
if (model->siesta.energy_shift_units)
  print_file_energy_units(fp, model->siesta.energy_shift_units);
else
  fprintf(fp, "\n");


//Default is FALSE.
if (model->siesta.harris_functional == TRUE)
    fprintf(fp, "Harris_functional    true\n");
else
    fprintf(fp, "Harris_functional    false\n");

//Default == LDA
switch (model->siesta.xc_functional_type)
{
    case GGA_XCFUNC: fprintf(fp, "XC.functional    GGA\n"); break;
    case LDA_XCFUNC: fprintf(fp, "XC.functional    LDA\n"); break;
    case LSD_XCFUNC: fprintf(fp, "XC.functional    LSD\n"); break;
    default:
        //should never get here
        error_table_entry("Inconsistancy model->siesta.xc_functional.\n");
        break;
}

switch (model->siesta.xc_author_type)
{
    case CA_XCAUTH:    fprintf(fp, "XC.Authors    CA\n"); break;
    case PZ_XCAUTH:    fprintf(fp, "XC.Authors    PZ\n"); break;
    case PW_XCAUTH:    fprintf(fp, "XC.Authors    PW92\n"); break;
    case PBE_XCAUTH:    fprintf(fp, "XC.Authors    PBE\n"); break;
    case RPBE_XCAUTH:    fprintf(fp, "XC.Authors    RPBE\n"); break;
    case LYP_XCAUTH:    fprintf(fp, "XC.Authors    LYP\n"); break;
    default:
        //should never get here
        error_table_entry("Inconsistancy model->siesta.xc_authors.\n");
        break;
}

if (model->siesta.spin_polarised == TRUE)
    fprintf (fp, "SpinPolarized    true\n");
else
    fprintf (fp, "SpinPolarized    false\n");



//MeshCutoff
fprintf (fp, "MeshCutoff    %f", model->siesta.mesh_cutoff);
if (print_file_energy_units(fp, model->siesta.mesh_cutoff_units) != 0)
{
    error_table_entry("Warning - Units missing MeshCutoff\n");
}

//kgrid_cutoff

fprintf(fp, "kgrid_cutoff    %f", model->siesta.kgrid_cutoff);
if (print_file_length_units(fp, model->siesta.kgrid_cutoff_units) !=0)
{
    error_table_entry("Warning - Units missing kgrid_cutoff\n");
}

fprintf(fp, "ElectronicTemperature    %f", model->siesta.electronic_temperature);
if (model->siesta.electronic_temperature_units != -1)
{
    if (print_file_temperature_units(fp, model->siesta.electronic_temperature_units) !=0)
    {
        error_table_entry("Warning - Units missing electronic_temperature_units\n");
    }
}


fprintf(fp, "MaxSCFIterations    %d\n", (int) model->siesta.no_of_cycles);
fprintf(fp, "DM.NumberPulay    %d\n", (int) model->siesta.no_of_pulay_matrices);
fprintf(fp, "DM.MixingWeight    %f\n", model->siesta.mixing_weight);

if (model->siesta.pulay_mixing)
{
    //what does this map to?
}


switch (model->siesta.run_type)
{
    case SINGLE_POINT:
        fprintf(fp, "MD.TypeOfRun    CG\n");
        break;
    case OPTIMISATION: fprintf(fp, "MD.TypeOfRun    Verlet\n"); break;
    case MOLECULAR_DYNAMICS: fprintf(fp, "MD.TypeOfRun    Verlet\n"); break;
    case PHONON_CALCULATION: 
        switch (model->siesta.md_type_of_run)
        {
            case FC_MDRUN:      fprintf(fp, "MD.TypeOfRun    FC\n"); break;
            case PHONON_MDRUN:  fprintf(fp, "MD.TypeOfRun    Phonon\n"); break;
            default:
            break;
        }
        break;        
                

    //MORE OPTIONS HERE > see manual
    default :
        fprintf(fp," \n");
        error_table_entry("Warning - MD.TypeOfRun *NOT* written out\n");
        break;
}


if (model->siesta.md_variable_cell)
    fprintf(fp, "MD.VariableCell    true\n");
else
    fprintf(fp, "MD.VariableCell    false\n");

if (model->siesta.md_num_cg_steps != -1.0)
{
    fprintf(fp, "MD.NumCGsteps    %d\n", (gint) model->siesta.md_num_cg_steps);
}

fprintf(fp, "MD.MaxCGDispl    %f", model->siesta.md_max_cg_displacement);
if (print_file_length_units(fp, model->siesta.md_max_cg_displacement_units) != 0)
{
    error_table_entry("Warning - Units missing MD.MaxCGDispl\n");
}


fprintf(fp, "MD.PreconditionVariableCell    %f", model->siesta.md_precondition_variable_cell);
if (print_file_length_units(fp, model->siesta.md_precondition_variable_cell_units ) != 0)
{
    error_table_entry("Warning - Units missing MD.PreconditionVariableCell\n");
}


fprintf(fp, "MD.MaxStressTol    %f", model->siesta.md_max_stress_tol);
if (print_file_pressure_units(fp, model->siesta.md_max_stress_tol_units) != 0)
{
    error_table_entry("Warning - Units missing MD.MaxStressTol\n");
}

fprintf(fp, "\n%%block MD.TargetStress\n");
fprintf(fp, "%f %f %f %f %f %f\n", model->siesta.md_target_stress_xx,
                                   model->siesta.md_target_stress_yy,
                                   model->siesta.md_target_stress_zz,
                                   model->siesta.md_target_stress_xy,
                                   model->siesta.md_target_stress_xz,
                                   model->siesta.md_target_stress_yz);
fprintf(fp, "%%endblock MD.TargetStress\n\n");



fprintf(fp, "MD.MaxForceTol    %f", model->siesta.md_max_force_tol);
if (print_file_force_units(fp, model->siesta.md_max_force_tol_units) != 0)
{
    error_table_entry("Warning - Units missing MD.MaxForceTol\n");
}

//run_type is not CG_MDRUN
fprintf(fp, "MD.InitialTimeStep    %d\n", (gint) model->siesta.md_inital_time_step);
fprintf(fp, "MD.FinalTimeStep    %d\n", (gint) model->siesta.md_final_time_step);
fprintf(fp, "MD.LengthTimeStep    %d", (gint) model->siesta.timestep);
if (print_file_time_units(fp, model->siesta.timestep_units) != 0)
{
    error_table_entry("Warning - Units missing MD.LengthTimeStep\n");
}

fprintf(fp, "MD.InitialTemperature    %f", model->siesta.md_inital_temperature);
if (print_file_temperature_units(fp, model->siesta.md_inital_temperature_units) != 0)
{
    error_table_entry("Warning - Units missing MD.InitialTemperature\n");
}

if (model->siesta.md_quench)
    fprintf(fp, "MD.Quench    true\n");
else
    fprintf(fp, "MD.Quench    false\n");

fprintf(fp, "MD.TargetTemperature    %f", model->siesta.md_inital_temperature);
if (print_file_temperature_units(fp, model->siesta.md_inital_temperature_units) !=0)
{
    error_table_entry("Warning - Units missing MD.TargetTemperature\n");
}


fprintf(fp, "MD.NoseMass    %f", model->siesta.md_nose_mass);
if (print_file_mominert_units(fp, model->siesta.md_nose_mass_units) != 0)
{
    error_table_entry("Warning - Units missing MD.NoseMass\n");
}

fprintf(fp, "MD.ParrinelloRahmanMass    %f", model->siesta.md_parrinello_rahman_mass);
if (print_file_mominert_units(fp, model->siesta.md_parrinello_rahman_mass_units) != 0)
{
    error_table_entry("Warning - Units missing MD.ParrinelloRahmanMass\n");
}

switch (model->siesta.md_anneal_option)
{
    case TEMPERATURE_ANNEAL:
        fprintf(fp, "MD.AnnealOption    Temperature\n");
        break;
    case PRESSURE_ANNEAL:
        fprintf(fp, "MD.AnnealOption    Pressure\n");
        break;
    case TEMP_AND_PRESS_ANNEAL:
        fprintf(fp, "MD.AnnealOption    TemperatureandPressure\n");
        break;
    default:
        fprintf(fp, "MD.AnnealOption    TemperatureandPressure\n");
        break;
}

fprintf(fp, "MD.TauRelax    %f", model->siesta.md_tau_relax);
if (print_file_time_units(fp, model->siesta.md_tau_relax_units) !=0)
{
    error_table_entry("Warning - Units missing MD.TauRelax\n");
}

fprintf(fp, "MD.BulkModulus    %f", model->siesta.md_bulk_modulus);
if (print_file_pressure_units(fp, model->siesta.md_bulk_modulus_units) !=0)
{
    error_table_entry("Warning - Units missing MD.BulkModulus\n");
}


fprintf(fp, "MD.TargetPressure    %f", model->siesta.md_target_pressure);
if (print_file_pressure_units(fp, model->siesta.md_target_pressure_units) !=0)
{
    error_table_entry("Warning - Units missing MD.BulkModulus\n");
}


fprintf(fp, "MD.FCDispl    %f", model->siesta.md_fc_displ);
if (print_file_length_units(fp, model->siesta.md_fc_displ_units) !=0)
{
    error_table_entry("Warning - Units missing MD.FCDispl\n");
}
fprintf(fp, "MD.FCfirst    %d\n", model->siesta.md_fc_first );
fprintf(fp, "MD.FClast    %d\n\n", model->siesta.md_fc_last);


if (model->siesta.use_saved_data)
    fprintf(fp, "UseSaveData\t\ttrue\n");
else
    fprintf(fp, "UseSaveData\t\tfalse\n");


if (model->siesta.file_output_write_coor_init)          //WriteCoorInital
    fprintf(fp, "WriteCoorInital    true\n");
else
    fprintf(fp, "WriteCoorInital    false\n");

if (model->siesta.file_output_write_coor_step)          //WriteCoorStep
    fprintf(fp, "WriteCoorStep    true\n");
else
    fprintf(fp, "WriteCoorStep    false\n");

if (model->siesta.file_output_write_forces)             //WriteForces
    fprintf(fp, "WriteForces\t\ttrue\n");
else
    fprintf(fp, "WriteForces\t\tfalse\n");

if (model->siesta.file_output_write_kpoints)            //WriteKpoints
    fprintf(fp, "WriteKpoints\t\ttrue\n");
else
    fprintf(fp, "WriteKpoints\t\tfalse\n");

if (model->siesta.file_output_write_eigenvalues)        //WriteEigenvalues
    fprintf(fp, "WriteEigenvalues\t\ttrue\n");
else
    fprintf(fp, "WriteEigenvalues\t\tfalse\n");

if (model->siesta.file_output_write_kbands)             //WriteKbands
    fprintf(fp, "WriteKbands\t\ttrue\n");
else
    fprintf(fp, "WriteKbands\t\tfalse\n");

if (model->siesta.file_output_write_bands)              //WriteBands
    fprintf(fp, "WriteBands\t\ttrue\n");
else
    fprintf(fp, "WriteBands\t\tfalse\n");

if (model->siesta.file_output_write_wavefunctions)      //WriteWaveFunctions
    fprintf(fp, "WriteWaveFunctions\t\ttrue\n");
else
    fprintf(fp, "WriteWaveFunctions\t\tfalse\n");

if (model->siesta.file_output_write_mullikenpop != -1)         //WriteMullikenPop
    fprintf(fp, "WriteMullikenPop\t\t%d\n", (int) model->siesta.file_output_write_mullikenpop);

if (model->siesta.file_output_write_dm)                 //WriteDM
    fprintf(fp, "WriteDM\t\ttrue\n");
else
    fprintf(fp, "WriteDM\t\tfalse\n");

if (model->siesta.file_output_write_coor_xmol)          //WriteCoorXmol
    fprintf(fp, "WriteCoorXmol\t\ttrue\n");
else
    fprintf(fp, "WriteCoorXmol\t\tfalse\n");

if (model->siesta.file_output_write_coor_cerius)        //WriteCoorCerius
    fprintf(fp, "WriteCoorCerius\t\ttrue\n");
else
    fprintf(fp, "WriteCoorCerius\t\tfalse\n");

if (model->siesta.file_output_write_md_xmol)            //WriteMDXmol
    fprintf(fp, "WriteMDXmol\t\ttrue\n");
else
    fprintf(fp, "WriteMDXmol\t\tfalse\n");

if (model->siesta.file_output_write_md_history)         //WriteMDhistory
    fprintf(fp, "WriteMDhistory\t\ttrue\n");
else
    fprintf(fp, "WriteMDhistory\t\tfalse\n");

/*
g_slist_free(normal_list);
g_slist_free(ghost_list);
*/

fclose(fp);
return(0);
}

/****************************************************/
/* attempt to read header at a particular precision */
/****************************************************/
gint siesta_header_read(FILE *fp, gint precision, gdouble *cell, gint *grid)
{
gint i, swap=FALSE;
int grid_int[4];
/*long grid_long[4];*/
/*float cell_float[9];*/
double cell_double[9];

/* standard return values */
for (i=9 ; i-- ; )
  cell[i] = 0.0;
for (i=4 ; i-- ; )
  grid[i] = 0;

/* precision dependent read */
switch (precision)
  {
  case 1:
    READ_RECORD;
    fread(cell_double, sizeof(double), 9, fp);
    READ_RECORD;

    READ_RECORD;
    fread(&grid_int, sizeof(int), 4, fp);
    READ_RECORD;

/* byte ordering check */
    for (i=9 ; i-- ; )
      {
      if (fabs(cell_double[i]) > 1000000.0)
        {
        swap = TRUE;
        break;
        }
      }
    if (swap)
      {
      for (i=9 ; i-- ; )
        swap_bytes(&cell_double[i], sizeof(double));
      for (i=4 ; i-- ; )
        swap_bytes(&grid_int[i], sizeof(int));
      }

/* pass back */
    for (i=9 ; i-- ; )
      cell[i] = cell_double[i];
    for (i=4 ; i-- ; )
      grid[i] = grid_int[i];
    break;

  default:
    printf("FIXME - unsuported precision.\n");
  }

return(swap);
}

/******************************************/
/* check for additional siesta data files */
/******************************************/
#define DEBUG_SIESTA_DENSITY 0
void siesta_density_files(const gchar *source, struct model_pak *model)
{
gint grid[4];
gint i, j, k, n, size, file_size, type, swap;
gchar *text, *name;
float *f;
gdouble cell[9], icell[9], *ptr;
FILE *fp;

g_assert(source != NULL);
g_assert(model != NULL);

#if DEBUG_SIESTA_DENSITY
printf("source: %s\n", source);
#endif

/* search for siesta output types */
/* in same directory & with same basename - different extension */
for (type=0 ; type<2 ; type++)
  {
/* remove extension */
  name = text = NULL;
  n = strlen(source);
  for (i=n ; i-- ; )
    {
    if (*(source+i) == '.')
      {
      text = g_strndup(source, i);
      break;
      }
    }
/* should never get a name with no extension */
  if (!text)
    return;

/* build current type name */
  switch (type)
    {
    case 0:
      name = g_strdup_printf("%s.RHO", text);
      break;

    case 1:
      name = g_strdup_printf("%s.VH", text);
      break;
    }
  g_free(text);
  if (!name)
    continue;

  fp = fopen(name, "rt");
  if (!fp)
    {
#if DEBUG_SIESTA_DENSITY
printf("Could not find: %s\n", name);
#endif
    continue;
    }

  file_size = file_byte_size(name);
  g_free(name);

  swap = siesta_header_read(fp, 1, cell, grid);

  memcpy(icell, cell, 9*sizeof(gdouble));
  matrix_invert(icell);

#if DEBUG_SIESTA_DENSITY
printf("swap bytes: %d\n", swap);
P3MAT(" cell: ", cell);
P3MAT("icell: ", icell);
printf("grid: %d x %d x %d : %d\n", grid[0], grid[1], grid[2], grid[3]);
#endif

/* NB: FORTRAN has a 4 byte separator either side of each record (single write() call) */
/*
real*8 cell(3,3)
real*4 rho(nmesh(1)*nmesh(2)*nmesh(3),nspin)
integer mesh(3)
integer nspin

    write(iu) cell
    write(iu) mesh, nspin
    ind = 0
    do is  = 1,nspin
       do iz = 1,mesh(3)
         do iy = 1,mesh(2)
            write(iu) rho(ind+ix,is),ix = 1,mesh(1)
            ind = ind + mesh(1)
         enddo
       enddo
    enddo
*/

/* data array size */
  n = grid[1] * grid[2] * grid[3];
  size = 8 + 4*grid[0];
  size *= n;
/* cell + mesh size */
  size += 104;

#if DEBUG_SIESTA_DENSITY
printf("Data size: %d bytes\n", size);
printf("File size: %d bytes\n", file_size);
#endif

  if (size != file_size)
    {
/* TODO - retry at different precision */
    gui_text_show(ERROR, "Failed to parse binary file.\n");
    continue;
    }

  memcpy(model->siesta.cell, cell, 9*sizeof(gdouble));
  memcpy(model->siesta.icell, icell, 9*sizeof(gdouble));
  memcpy(model->siesta.grid, grid, 4*sizeof(gint));

/* CURRENT - test for a scaling factor (ie au output instead of angstrom) */
matmat(model->latmat, icell);

if (icell[0] > 0.2)
  {
#if DEBUG_SIESTA_DENSITY
printf("density grid in au, scale = %f\n", icell[0]);
#endif
  model->siesta.eden_scale = icell[0];
  }

/* CURRENT - read function data */
  switch (type)
    {
    case 0:
      model->siesta.eden = g_malloc(n*grid[0]*sizeof(gdouble));
      ptr = model->siesta.eden;
      break;

    case 1:
      model->siesta.epot = g_malloc(n*grid[0]*sizeof(gdouble));
      ptr = model->siesta.epot;
      break;

    default:
      ptr = NULL;
    }

  f = g_malloc(grid[0]*sizeof(float));
  k=0;
  for (i=0 ; i<n ; i++)
    {
    READ_RECORD;
    fread(f, sizeof(float), grid[0], fp);
    READ_RECORD;

    if (ptr)
      {
      if (swap)
        {
        for (j=0 ; j<grid[0] ; j++)
          {
          swap_bytes(f+j, sizeof(float));
          *(ptr+k++) = *(f+j);
          }
        }
      else
        {
        for (j=0 ; j<grid[0] ; j++)
          *(ptr+k++) = *(f+j);
        }
      }

/*
        *(model->siesta.eden+k++) = *(f+j);
*/

/* DEBUG */
/*
if (i == n-1)
  {
  for (j=0 ; j<grid[0] ; j++)
    printf("%f\n", *(model->siesta.eden+j));
  }
*/
    }
  g_free(f);

/* this shouldn't happen (file moved while reading?) */
  if (k != n*grid[0])
    {
    gui_text_show(ERROR, "Fatal error reading binary file.\n");
/*
  g_free(model->siesta.eden);
  model->siesta.eden = NULL;
*/
    }
  }
}

/***********************************************************************/
/* process lattice constant line to produce a cartesian scaling factor */
/***********************************************************************/
#define DEBUG_SIESTA_SCALE 0
gdouble siesta_scale_parse(const gchar *text)
{
gint num_tokens;
gchar **buff;
gdouble scale = 1.0;

g_assert(text != NULL);

#if DEBUG_SIESTA_SCALE
printf("input line: %s", text);
#endif

/* process scale */
buff = tokenize(text, &num_tokens);
if (num_tokens > 1)
  scale *= str_to_float(*(buff+1));

/* process units */
if (num_tokens > 2)
  if ((g_ascii_strncasecmp(*(buff+2), "au", 2) == 0) ||
      (g_ascii_strncasecmp(*(buff+2), "Bohr", 4) == 0))
    scale *= AU2ANG;  

g_strfreev(buff);

#if DEBUG_SIESTA_SCALE
printf("output scale: %f\n", scale);
#endif

return(scale);
}

/********************************/
/* process a species block line */
/********************************/
void siesta_parse_species(const gchar *line, GSList **list)
{
gint num_tokens;
gchar **buff;
struct species_pak *species;

g_assert(line != NULL);

buff = tokenize(line, &num_tokens);

if (num_tokens > 2)
  {
  species = g_malloc(sizeof(struct species_pak));

  species->number = str_to_float(*(buff+1));
  species->label = g_strdup(*(buff+2));
  *list = g_slist_append(*list, species);
  }

g_strfreev(buff);
}

/************************************/
/* process a lattice parameter line */
/************************************/
void siesta_parse_cell(const gchar *line, gdouble scale, struct model_pak *model)
{
gint num_tokens;
gchar **buff;

g_assert(line != NULL);

buff = tokenize(line, &num_tokens);
if (num_tokens > 5)
  {
  model->pbc[0] = scale * str_to_float(*(buff+0));
  model->pbc[1] = scale * str_to_float(*(buff+1));
  model->pbc[2] = scale * str_to_float(*(buff+2));
  model->pbc[3] = D2R*str_to_float(*(buff+3));
  model->pbc[4] = D2R*str_to_float(*(buff+4));
  model->pbc[5] = D2R*str_to_float(*(buff+5));
  }

model->construct_pbc = FALSE;

/* hack for determining if it's a (GDIS created) surface */
if (fabs(model->pbc[2]) < FRACTION_TOLERANCE)
  model->periodic = 2;
else
  model->periodic = 3;

g_strfreev(buff);
}

/*********************************/
/* process a lattice vector line */
/*********************************/
void siesta_parse_lattice(const gchar *line, gint n, gdouble scale, struct model_pak *model)
{
gint num_tokens;
gchar **buff;

g_assert(line != NULL);

buff = tokenize(line, &num_tokens);
if (num_tokens > 2)
  {
  model->latmat[0+n] = scale * str_to_float(*(buff+0));
  model->latmat[3+n] = scale * str_to_float(*(buff+1));
  model->latmat[6+n] = scale * str_to_float(*(buff+2));
  }

model->construct_pbc = TRUE;
model->periodic = 3;

g_strfreev(buff);
}

/***********************************/
/* process a coordinate block line */
/***********************************/
void siesta_parse_coords(const gchar *line, GSList *species_list, gdouble scale, struct model_pak *model)
{
gint i, num_tokens;
gchar **buff;
struct species_pak *species;
struct core_pak *core;

g_assert(line != NULL);

buff = tokenize(line, &num_tokens);

if (num_tokens > 3)
  {
/* find corresponding number of this element in the species block */
/* NB: list counts from 0, fdf counts from 1 (hence the minus 1) */
  if (species_list)
    {
    i = abs((gint) str_to_float(*(buff+3)) - 1);
    species = g_slist_nth_data(species_list, i);

    core = new_core(species->label, model);
    
/* translucent ghost atom */
    if (species->number < 0)
      {
      core->ghost = TRUE;
      core->colour[3] = 0.5;
      }
    }
  else
    core = new_core("X", model);
 
  model->cores = g_slist_prepend(model->cores, core);
    
  core->x[0] = str_to_float(*(buff+0));
  core->x[1] = str_to_float(*(buff+1));
  core->x[2] = str_to_float(*(buff+2));

  VEC3MUL(core->x, scale);
  }

g_strfreev(buff);
}

/***************************/
/* zmatrix line processing */
/***************************/
void siesta_parse_zmatrix(const gchar *text,
                          GSList *species_list,
                          gdouble scale,
                          struct model_pak *model)
{
gint i, num_tokens;
gchar **buff;
struct core_pak *core;
struct species_pak *species_data;
static gint state=0;

/* keyword processing within zmatrix block */
if (g_ascii_strncasecmp("molecule", text, 8) == 0)
  {
  state = 1;
  if (strstr(text, "frac"))
    zmat_fractional_set(model->zmatrix);
  else
    zmat_cartesian_set(model->zmatrix);
  return;
  }

if (g_ascii_strncasecmp("variable", text, 8) == 0)
  {
  state = 2;
  return;
  }

if (g_ascii_strncasecmp("constant", text, 8) == 0)
  {
  state = 3;
  return;
  }

if (g_ascii_strncasecmp("fractional", text, 10) == 0)
  {
  state = 4;
  model->fractional = TRUE;
  return;
  }

if (g_ascii_strncasecmp("cartesian", text, 9) == 0)
  {
  state = 4;
  model->fractional = FALSE;
  return;
  }

/* data processing within zmatrix block */
switch (state)
  {
  case 1:
/* process zmatrix coords */
    zmat_core_add(text, model->zmatrix);
    break;

  case 2:
/* process zmatrix variables */
    zmat_var_add(text, model->zmatrix);
    break;

  case 3:
/* process zmatrix constants */
    zmat_const_add(text, model->zmatrix);
    break;

  case 4:
/* fractional input coords */
    buff = tokenize(text, &num_tokens);
    if (num_tokens > 3)
      {
/* find corresponding number of this element in the species block */
/* NB: list counts from 0, fdf counts from 1 (hence the minus 1) */
      if (species_list)
        {
        i = abs((gint) str_to_float(*(buff)) - 1);
        species_data = g_slist_nth_data(species_list, i);

        core = new_core(species_data->label, model);

/* translucent ghost atom */
        if (species_data->number < 0)
          {
          core->ghost = TRUE;
          core->colour[3] = 0.5;
          }
        }
      else
        core = new_core("X", model);

      model->cores = g_slist_prepend(model->cores, core);

      core->x[0] = scale * str_to_float(*(buff+1));
      core->x[1] = scale * str_to_float(*(buff+2));
      core->x[2] = scale * str_to_float(*(buff+3));
      }
    g_strfreev(buff);
    break;
  }
}

/******************************/
/* main fdf data block reader */
/******************************/
#define DEBUG_READ_FDF_BLOCK 0
gint read_fdf_block(FILE *fp, struct model_pak *model, gboolean read_inital_cooords)
{
gint count, num_tokens, block_type, siesta_sucks = FALSE;
gchar **buff, *line;
gdouble scale_lattice = 1.0, scale_coords = 1.0;
GSList *species_list=NULL;

g_assert(fp != NULL);
g_assert(model != NULL);

for (;;)
  {
  line = file_read_line(fp);

/* terminate if NULL returned */
  if (!line)
    break;

  buff = tokenize(line, &num_tokens);

/* get next line if blank */
  if (!buff)
    continue;

/* NEW - block/enblock processing */
  if (g_ascii_strncasecmp("%block", line, 6) == 0)
    {
    block_type = SIESTA_DEFAULT;

    if (num_tokens > 1)
      {
      if (g_ascii_strncasecmp("ChemicalSpecies", *(buff+1), 15) == 0)
        block_type = SIESTA_SPECIES;
      if (g_ascii_strncasecmp("LatticeParameters", *(buff+1), 17) == 0)
        block_type = SIESTA_LATTICE_PARAMETERS;
      if (g_ascii_strncasecmp("LatticeVectors", *(buff+1), 14) == 0)
        block_type = SIESTA_LATTICE_VECTORS;
      if (g_ascii_strncasecmp("AtomicCoordinates", *(buff+1), 17) == 0)
        block_type = SIESTA_COORDS;
      if (g_ascii_strncasecmp("zmatrix", *(buff+1), 7) == 0)
        {
        block_type = SIESTA_ZMATRIX;
        }
      }

#if DEBUG_READ_FDF_BLOCK
printf("processing block [type %d]\n", block_type);
#endif

    count = 0;
    for (;;)
      {
      line = file_read_line(fp);
      if (!line)
        goto siesta_done_file;
      if (g_ascii_strncasecmp("%endblock", line, 9) == 0)
        {
#if DEBUG_READ_FDF_BLOCK
printf("end of block.\n");
#endif
        goto siesta_done_line;
        }

      switch (block_type)
        {
        case SIESTA_SPECIES:
          siesta_parse_species(line, &species_list);
          break;

        case SIESTA_LATTICE_PARAMETERS:
          siesta_parse_cell(line, scale_lattice, model);
          break;

        case SIESTA_LATTICE_VECTORS:
          siesta_parse_lattice(line, count, scale_lattice, model);
          break;

        case SIESTA_COORDS:
          siesta_parse_coords(line, species_list, scale_coords, model);
          break;

        case SIESTA_ZMATRIX:
          siesta_parse_zmatrix(line, species_list, scale_coords, model);
          break;

        default:
          error_table_entry("Unrecognized block encountered.\n");
          break;
        }
      count++;
      }
/* done block processing - get next line */
    continue;
    }

/* cartesian/fractional */
  if (read_inital_cooords)
    {
    if (g_ascii_strncasecmp("AtomicCoordinatesFormat", line, 23) == 0)
      {
      if (g_strrstr(line, "ractional"))
        model->fractional = TRUE;
      else
        {
        model->fractional = FALSE;

/* cope with cartesians that arent really cartesian */
        if (g_strrstr(line, "ScaledCartesian"))
          siesta_sucks = TRUE;

        if (g_strrstr(line, "Bohr"))
          scale_coords = AU2ANG;
        }
      }
    }

/* scaling */
  if (g_ascii_strncasecmp("LatticeConstant", *buff, 15) == 0)
    scale_lattice = siesta_scale_parse(line);

  if (g_ascii_strncasecmp("NumberOfAtoms", *buff, 13) == 0)
    {
    if (num_tokens > 1)
      model->siesta.num_atoms = (gint) str_to_float(*(buff+1));
    }

  if (g_ascii_strncasecmp("SystemLabel", *buff, 11) == 0)
    {
    if (num_tokens > 1)
      {
      g_free(model->basename);
      model->basename = g_strdup(*(buff+1));
      }
    }

//Terry Siesta options - arrrgh!

  if (g_ascii_strncasecmp("kgrid_cutoff", *buff, 12) == 0)
    {
    if (num_tokens >= 3)             //kgrid_cutoff real length_type #othercrap
      {
      model->siesta.kgrid_cutoff = str_to_float(*(buff+1));
      if (set_length_units(*(buff+2), &model->siesta.kgrid_cutoff_units) !=0)
        error_table_entry("Units missing kgrid_cutoff.\n");
      }
    else
      error_table_entry("kgrid_cutoff line corrupt.\n");
    }

  if (g_ascii_strncasecmp("WriteCoorInitial", *buff, 16) == 0)
    {
    if (num_tokens >= 2)
      {
      if (g_ascii_strncasecmp("true", *(buff+1), 4) == 0)
          model->siesta.file_output_write_coor_init = TRUE;
      else if (g_ascii_strncasecmp(".true", *(buff+1), 5) == 0)
          model->siesta.file_output_write_coor_init = TRUE;
      else
          model->siesta.file_output_write_coor_init = FALSE;
      }
    else
      {
      error_table_entry("WriteCoorInitial line corrupt\n");
      }
    }

  if (g_ascii_strncasecmp("WriteCoorStep", *buff, 13) == 0)
    {
    if (num_tokens >= 2)
      {
      if (g_ascii_strncasecmp("true", *(buff+1), 4) == 0)
        model->siesta.file_output_write_coor_step = TRUE;
      else if (g_ascii_strncasecmp(".true", *(buff+1), 5) == 0)
        model->siesta.file_output_write_coor_step = TRUE;
      else
        model->siesta.file_output_write_coor_step = FALSE;
      }
    else
      error_table_entry("WriteCoorStep line corrupt\n");
    }

  if (g_ascii_strncasecmp("WriteForces", *buff, 11) == 0)
    {
    if (num_tokens >= 2)
      {
      if (g_ascii_strncasecmp("true", *(buff+1), 4) == 0)
        model->siesta.file_output_write_forces = TRUE;
      else if (g_ascii_strncasecmp(".true", *(buff+1), 5) == 0)
        model->siesta.file_output_write_forces = TRUE;
      else
        model->siesta.file_output_write_forces = FALSE;
      }
    else
      error_table_entry("WriteForcesStep line corrupt\n");
    }

  if (g_ascii_strncasecmp("WriteKpoints", *buff, 12) == 0)
    {
    if (num_tokens >= 2)
      {
      if (g_ascii_strncasecmp("true", *(buff+1), 4) == 0)
        model->siesta.file_output_write_kpoints = TRUE;
      else if (g_ascii_strncasecmp(".true", *(buff+1), 5) == 0)
        model->siesta.file_output_write_kpoints = TRUE;
      else
        model->siesta.file_output_write_kpoints = FALSE;
      }
    else
      error_table_entry("WriteKpoints line corrupt\n");
    }

  if (g_ascii_strncasecmp("WriteEigenvalues", *buff, 16) == 0)
    {
    if (num_tokens >= 2)
      {
      if (g_ascii_strncasecmp("true", *(buff+1), 4) == 0)
        model->siesta.file_output_write_eigenvalues = TRUE;
      else if (g_ascii_strncasecmp(".true", *(buff+1), 5) == 0)
        model->siesta.file_output_write_eigenvalues = TRUE;
      else
        model->siesta.file_output_write_eigenvalues = FALSE;
      }
    else
      error_table_entry("WriteEigenvalues line corrupt\n");
    }

  if (g_ascii_strncasecmp("WriteKbands", *buff, 11) == 0)
    {
    if (num_tokens >= 2)
      {
      if (g_ascii_strncasecmp("true", *(buff+1), 4) == 0)
        model->siesta.file_output_write_kbands = TRUE;
      else if (g_ascii_strncasecmp(".true", *(buff+1), 5) == 0)
        model->siesta.file_output_write_kbands = TRUE;
      else
        model->siesta.file_output_write_kbands = FALSE;
      }
    else
      error_table_entry("WriteKbands line corrupt\n");
    }

  if (g_ascii_strncasecmp("WriteBands", *buff, 10) == 0)
    {
    if (num_tokens >= 2)
      {
      if (g_ascii_strncasecmp("true", *(buff+1), 4) == 0)
        model->siesta.file_output_write_bands = TRUE;
      else if (g_ascii_strncasecmp(".true", *(buff+1), 5) == 0)
        model->siesta.file_output_write_bands = TRUE;
      else
        model->siesta.file_output_write_bands = FALSE;
      }
    else
      error_table_entry("WriteBands line corrupt\n");
    }

  if (g_ascii_strncasecmp("WriteWaveFunctions", *buff, 18) == 0)
    {
    if (num_tokens >= 2)
      {
      if (g_ascii_strncasecmp("true", *(buff+1), 4) == 0)
        model->siesta.file_output_write_wavefunctions = TRUE;
      else if (g_ascii_strncasecmp(".true", *(buff+1), 5) == 0)
        model->siesta.file_output_write_wavefunctions = TRUE;
      else
        model->siesta.file_output_write_wavefunctions = FALSE;
      }
    else
      error_table_entry("WriteWaveFunctions line corrupt\n");
    }

  if (g_ascii_strncasecmp("WriteDM", *buff, 7) == 0)
    {
    if (num_tokens >= 2)
      {
      if (g_ascii_strncasecmp("true", *(buff+1), 4) == 0)
        model->siesta.file_output_write_dm = TRUE;
      else if (g_ascii_strncasecmp(".true", *(buff+1), 5) == 0)
        model->siesta.file_output_write_dm = TRUE;
      else
        model->siesta.file_output_write_dm = FALSE;
      }
    else
      error_table_entry("WriteDM line corrupt\n");
    }

  if (g_ascii_strncasecmp("WriteCoorXmol", *buff, 13) == 0)
    {
    if (num_tokens >= 2)
      {
      if (g_ascii_strncasecmp("true", *(buff+1), 4) == 0)
        model->siesta.file_output_write_coor_xmol = TRUE;
      else if (g_ascii_strncasecmp(".true", *(buff+1), 5) == 0)
        model->siesta.file_output_write_coor_xmol = TRUE;
      else
        model->siesta.file_output_write_coor_xmol = FALSE;
      }
    else
      error_table_entry("WriteCoorXmol line corrupt\n");
    }

  if (g_ascii_strncasecmp("WriteCoorCerius", *buff, 15) == 0)
    {
    if (num_tokens >= 2)
      {
      if (g_ascii_strncasecmp("true", *(buff+1), 4) == 0)
        model->siesta.file_output_write_coor_cerius = TRUE;
      else if (g_ascii_strncasecmp(".true", *(buff+1), 5) == 0)
        model->siesta.file_output_write_coor_cerius = TRUE;
      else
        model->siesta.file_output_write_coor_cerius = FALSE;
      }
    else
      error_table_entry("WriteCoorCerius line corrupt\n");
    }

  if (g_ascii_strncasecmp("WriteMDXmol", *buff, 11) == 0)
    {
    if (num_tokens >= 2)
      {
      if (g_ascii_strncasecmp("true", *(buff+1), 4) == 0)
        model->siesta.file_output_write_md_xmol = TRUE;
      else if (g_ascii_strncasecmp(".true", *(buff+1), 5) == 0)
        model->siesta.file_output_write_md_xmol = TRUE;
      else
        model->siesta.file_output_write_md_xmol = FALSE;
      }
    else
      error_table_entry("WriteMDXmol line corrupt\n");
    }

  if (g_ascii_strncasecmp("WriteMDhistory", *buff, 14) == 0)
    {
    if (num_tokens >= 2)
      {
      if (g_ascii_strncasecmp("true", *(buff+1), 4) == 0)
        model->siesta.file_output_write_md_history = TRUE;
      else if (g_ascii_strncasecmp(".true", *(buff+1), 5) == 0)
        model->siesta.file_output_write_md_history = TRUE;
      else
        model->siesta.file_output_write_md_history = FALSE;
      }
    else
      error_table_entry("WriteMDhistory line corrupt\n");
    }
  
  if (g_ascii_strncasecmp("LongOutput", *buff, 10) == 0)
    {
    if (num_tokens >= 2)
      {
      if ((g_ascii_strncasecmp("true", *buff, 4) == 0) ||
         (g_ascii_strncasecmp(".true", *buff, 5) == 0) )
        {
        model->siesta.long_output = TRUE;
        model->siesta.file_output_write_coor_step = TRUE;
        model->siesta.file_output_write_forces = TRUE;
        model->siesta.file_output_write_kpoints = TRUE;
        model->siesta.file_output_write_eigenvalues = TRUE;
        model->siesta.file_output_write_bands = TRUE;
        model->siesta.file_output_write_kbands = TRUE;
        model->siesta.file_output_write_wavefunctions = TRUE;
        model->siesta.file_output_write_mullikenpop = 1.0;
        }
      else
        {
        model->siesta.long_output = FALSE;
        model->siesta.file_output_write_coor_step = FALSE;
        model->siesta.file_output_write_forces = FALSE;
        model->siesta.file_output_write_kpoints = FALSE;
        model->siesta.file_output_write_eigenvalues = FALSE;
        model->siesta.file_output_write_bands = FALSE;
        model->siesta.file_output_write_kbands = FALSE;
        model->siesta.file_output_write_wavefunctions = FALSE;
        model->siesta.file_output_write_mullikenpop = 0.0;
        }
      }
    else
      {
      model->siesta.long_output = FALSE;
      model->siesta.file_output_write_coor_step = FALSE;
      model->siesta.file_output_write_forces = FALSE;
      model->siesta.file_output_write_kpoints = FALSE;
      model->siesta.file_output_write_eigenvalues = FALSE;
      model->siesta.file_output_write_bands = FALSE;
      model->siesta.file_output_write_kbands = FALSE;
      model->siesta.file_output_write_wavefunctions = FALSE;
      model->siesta.file_output_write_mullikenpop = 0.0;

      error_table_entry("LongOutput line corrupt\n");
      }
    }

  if (g_ascii_strncasecmp("PAO.SplitNorm", *buff, 13) == 0)
    {
    if (num_tokens >= 2)
      model->siesta.split_zeta_norm = str_to_float(*(buff+1));
    else
      error_table_entry("PAO.SplitNorm line corrupt : using defaults.\n");
    }

  if (g_ascii_strncasecmp("PAO.BasisSize", *buff, 13) == 0)
    {
    if (num_tokens == 2)
      {
      if (g_ascii_strncasecmp("SZ", *(buff+1), 2) == 0)
          model->siesta.basis_set = SZ_ZETA;

      if (g_ascii_strncasecmp("SZP", *(buff+1), 3) == 0)
        model->siesta.basis_set = SZP_ZETA;

      if (g_ascii_strncasecmp("DZ", *(buff+1), 2) == 0)
        model->siesta.basis_set = DZ_ZETA;

      if (g_ascii_strncasecmp("DZP", *(buff+1), 3) == 0)
        model->siesta.basis_set = DZP_ZETA;
      }
    else
      {
      model->siesta.basis_set = DZP_ZETA;          
      error_table_entry("PAO.BasisSize line corrupt : using defaults\n");
      }
    }
  
  if (g_ascii_strncasecmp("PAO.EnergyShift", *buff, 15) == 0)
    {
    if (num_tokens >= 3)
      {
      model->siesta.energy_shift = str_to_float(*(buff+1));
      if (set_energy_units(*(buff+2), &model->siesta.energy_shift_units) != 0)
        error_table_entry("PAO.EnergyShift line corrupt : using defaults\n");
      }
    else
      error_table_entry("PAO.EnergyShift line corrupt : using defaults\n");
    }
  
  if (g_ascii_strncasecmp("SolutionMethod", *buff, 14) == 0)
    {
    if (g_ascii_strncasecmp("OrderN", *(buff+1), 6) == 0)
      model->siesta.solution_method_type = ORDERN_SOLMETH;
    else if (g_ascii_strncasecmp("diagon", *(buff+1), 6) == 0)
      model->siesta.solution_method_type = DIAGON_SOLMETH;
    else
      error_table_entry("SolutionMethod line corrupt\n");
    }

  if (g_ascii_strncasecmp("XC.functional", *buff, 13) == 0)
    {
    if (g_ascii_strncasecmp("LDA", *(buff+1), 3) == 0)
        model->siesta.xc_functional_type = LDA_XCFUNC;
    else if (g_ascii_strncasecmp("LSD", *(buff+1), 3) == 0)
        model->siesta.xc_functional_type = LSD_XCFUNC;
    else if (g_ascii_strncasecmp("GGA", *(buff+1), 3) == 0)
        model->siesta.xc_functional_type = GGA_XCFUNC;
    else
      error_table_entry("XC.functional line corrupt\n");
    }

  if (g_ascii_strncasecmp("XC.authors", *buff, 10) == 0)
    {
    if (g_ascii_strncasecmp("PZ", *(buff+1), 2) == 0)
        model->siesta.xc_author_type = PZ_XCAUTH;
    else if (g_ascii_strncasecmp("CA", *(buff+1), 2) == 0)
        model->siesta.xc_author_type = CA_XCAUTH;
    else if (g_ascii_strncasecmp("PW92", *(buff+1), 4) == 0)
        model->siesta.xc_author_type = PW_XCAUTH;
    else if (g_ascii_strncasecmp("PBE", *(buff+1), 3) == 0)
        model->siesta.xc_author_type = PBE_XCAUTH;
    else if (g_ascii_strncasecmp("RPBE", *(buff+1), 4) == 0)
        model->siesta.xc_author_type = RPBE_XCAUTH;
    else if (g_ascii_strncasecmp("LYP", *(buff+1), 4) == 0)
        model->siesta.xc_author_type = LYP_XCAUTH;
    else
      error_table_entry("XC.authors line corrupt\n");
    }

//mesh_cutoff -> MeshCutoff real siesta_unit_type
  if (g_ascii_strncasecmp("MeshCutoff", *buff, 10) == 0)
    {
    if (num_tokens >= 3)
      {
      model->siesta.mesh_cutoff = str_to_float(*(buff+1));
      if (set_energy_units(*(buff+2), &model->siesta.mesh_cutoff_units) != 0)
        error_table_entry("MeshCutoff UNITS corrupt : using default\n");
      }
    else
      error_table_entry("MeshCutoff line corrupt : using defaults\n");
    }

  if (g_ascii_strncasecmp("Harris_functional", *buff, 17) == 0)
    {
    if (num_tokens >=2)         //Harris_functional boolean #crap
      {
      if (g_ascii_strncasecmp("true", *(buff+1), 4) == 0)
        model->siesta.harris_functional = TRUE;
      else if (g_ascii_strncasecmp(".true", *(buff+1), 5) == 0)
        model->siesta.harris_functional = TRUE;
      else
        model->siesta.harris_functional = FALSE;
      }
    else
      error_table_entry("Harris_functional line corrupt\n");
    }

  if (g_ascii_strncasecmp("ElectronicTemperature", *buff, 21) == 0)
    {
    if (num_tokens > 2)             //temp real units #othercrap
      {
      model->siesta.electronic_temperature = str_to_float(*(buff+1));
      if (set_temperature_units(*(buff+2), &model->siesta.electronic_temperature_units) != 0)
        error_table_entry("ElectronicTemperature UNITS corrupt : using default\n");
      }
    else if (num_tokens == 2)
      {
      model->siesta.electronic_temperature = str_to_float(*(buff+1));
//default temp units? KELVIN FOR NOW!
      model->siesta.electronic_temperature_units = K_TEMP;
      error_table_entry("Units missing ElectronicTemperature : using Kelvin\n");
      }
    else
      error_table_entry("ElectronicTemperature line corrupt : using defaults\n");
    }

//spin_polarised -> SpinPolarized    true||.true.
  if (g_ascii_strncasecmp("SpinPolarized", *buff, 13) == 0)
    {
    if (num_tokens >= 2)             //SpinPolarized boolean #othercrap
      {
      if (g_ascii_strncasecmp("true", *(buff+1), 4) == 0)
        model->siesta.spin_polarised = TRUE;
      else if (g_ascii_strncasecmp(".true", *(buff+1), 5) == 0)
        model->siesta.spin_polarised = TRUE;
      else
        model->siesta.spin_polarised = FALSE;
      }
    else
      error_table_entry("SpinPolarized line corrupt : using defaults\n");
    }

//SIESTA - SCF OPTIONS
// -----> right no_of_cycles -> MaxSCFIterations

  if (g_ascii_strncasecmp("MaxSCFIterations", *buff, 16) == 0)
    {
    if (num_tokens >= 2)
      model->siesta.no_of_cycles = str_to_float(*(buff+1));
    else
      error_table_entry("MaxSCFIterations line corrupt : using defaults\n");
    }

  if (g_ascii_strncasecmp("DM.NumberPulay", *buff, 14) == 0)
    {
    if (num_tokens >= 2)
      model->siesta.no_of_pulay_matrices = str_to_float(*(buff+1));
    else
      error_table_entry("DM.NumberPulay line corrupt : using defaults\n");
    }

//mixing_weight -> DM.MixingWeight float #comments
  if (g_ascii_strncasecmp("DM.MixingWeight", *buff, 15) == 0)
    {
    if (num_tokens >= 2)
      model->siesta.mixing_weight = str_to_float(*(buff+1));
    else
      error_table_entry("DM.MixingWeight line corrupt : using defaults\n");
    }

//run_type - enum { SINGLE_POINT, OPTIMISATION, MOLECULAR_DYNAMICS, PHONON_CALCULATION };
  if (g_ascii_strncasecmp("MD.TypeOfRun", *buff, 12) == 0)
    {
    if (num_tokens >= 2)
      {
      if (g_ascii_strncasecmp("CG", *(buff+1), 2) == 0)
        {
//soo is it a single point or an optimisation?
//have to check for number of steps
        if (model->siesta.number_of_steps == 0)
          model->siesta.run_type = SINGLE_POINT;
        else
          model->siesta.run_type = OPTIMISATION;
        }
      else if (g_ascii_strncasecmp("Verlet", *(buff+1), 6) == 0)
        {
        model->siesta.run_type = MOLECULAR_DYNAMICS;
        model->siesta.md_type_of_run = VERLET_MDRUN;
        }
      else if (g_ascii_strncasecmp("Nose", *(buff+1), 4) == 0)
        {
        model->siesta.run_type = MOLECULAR_DYNAMICS;
        model->siesta.md_type_of_run = NOSE_MDRUN;
        }
      else if (g_ascii_strncasecmp("ParrinelloRahman", *(buff+1), 16) == 0)
        {
        model->siesta.run_type = MOLECULAR_DYNAMICS;
        model->siesta.md_type_of_run = PARRINELLOPAHMAN_MDRUN;
        }
      else if (g_ascii_strncasecmp("NoseParrinelloRahman", *(buff+1), 20) == 0)
        {
        model->siesta.run_type = MOLECULAR_DYNAMICS;
        model->siesta.md_type_of_run = NOSEPARRINELLOPAHMAN_MDRUN;
        }
      else if (g_ascii_strncasecmp("Anneal", *(buff+1), 6) == 0)
        {
        model->siesta.run_type = MOLECULAR_DYNAMICS;
        model->siesta.md_type_of_run = ANNEAL_MDRUN;
        }
      else if (g_ascii_strncasecmp("FC", *(buff+1), 2) == 0)
        {
        model->siesta.run_type = PHONON_CALCULATION;
        model->siesta.md_type_of_run = FC_MDRUN;
        }
      else if (g_ascii_strncasecmp("Phonon", *(buff+1), 6) == 0)
        {
        model->siesta.run_type = PHONON_CALCULATION;
        model->siesta.md_type_of_run = PHONON_MDRUN;
        }
      }
    else
      error_table_entry("MD.TypeOfRun line corrupt");
    }

  //model->siesta.md_variable_cell
  if (g_ascii_strncasecmp("MD.VariableCell", *buff, 15) == 0)
    {
    if (num_tokens >=2)
      {
      if (set_true_false(*(buff+1), &model->siesta.md_variable_cell) != 0)
        error_table_entry("Failed to set MD.VariableCell\n");
      }
    else
      error_table_entry("MD.VariableCell line corrupt : using defaults\n");
    }

  //model->siesta.md_num_cg_steps
  if (g_ascii_strncasecmp("MD.NumCGsteps", *buff, 13) == 0)
    {
    if (num_tokens >= 2)
      model->siesta.md_num_cg_steps = str_to_float(*(buff+1));
    else
      error_table_entry("MD.NumCGsteps line corrupt\n");
    }

  //model->siesta.md_max_cg_displacement
  if (g_ascii_strncasecmp("MD.MaxCGDispl", *buff, 13) == 0)
    {
    if (num_tokens >= 3)
      {
      model->siesta.md_max_cg_displacement = str_to_float(*(buff+1));
      if (set_length_units(*(buff+2), &model->siesta.md_max_cg_displacement_units) != 0)
        error_table_entry("MD.MaxCGDispl units unable to be set\n");
      }
    else
      error_table_entry("MD.MaxCGDispl line corrupt\n");
    }

  //model->siesta.md_precondition_variable_cell
  if (g_ascii_strncasecmp("MD.PreconditionVariableCell", *buff, 27) == 0)
    {
    if (num_tokens >= 3)
      {
      model->siesta.md_precondition_variable_cell = str_to_float(*(buff+1));
      if (set_length_units(*(buff+2), &model->siesta.md_precondition_variable_cell_units) != 0)
        error_table_entry("MD.PreconditionVariableCell units unable to be set\n");
      }
    else
        error_table_entry("MD.PreconditionVariableCell line corrupt\n");
    }

  //model->siesta.md_max_force_tol
  if (g_ascii_strncasecmp("MD.MaxForceTol", *buff, 14) == 0)
    {
    if (num_tokens >= 3)
      {
      model->siesta.md_max_force_tol = str_to_float(*(buff+1));
      if (set_force_units(*(buff+2), &model->siesta.md_max_force_tol_units) != 0)
        error_table_entry("MD.MaxForceTol units unable to be set\n");
      }
    else
      error_table_entry("MD.MaxForceTol line corrupt\n");
    }

  //model->siesta.md_max_stress_tol
  if (g_ascii_strncasecmp("MD.MaxStressTol", *buff, 14) == 0)
    {
    if (num_tokens >= 3)
      {
      model->siesta.md_max_stress_tol = str_to_float(*(buff+1));
      if (set_pressure_units(*(buff+2), &model->siesta.md_max_stress_tol_units) != 0)
        error_table_entry("MD.MaxStressTol units unable to be set\n");
      }
    else
      error_table_entry("MD.MaxStressTol line corrupt\n");
    }

  //model->siesta.md_inital_time_step
  if (g_ascii_strncasecmp("MD.InitialTimeStep", *buff, 18) == 0)
    {
    if (num_tokens >= 2)
      model->siesta.md_inital_time_step = str_to_float(*(buff+1));
    else
      error_table_entry("MD.InitialTimeStep line corrupt\n");
    }

  //model->siesta.md_final_time_step
  if (g_ascii_strncasecmp("MD.FinalTimeStep", *buff, 16) == 0)
    {
    if (num_tokens >= 2)
      model->siesta.md_final_time_step = str_to_float(*(buff+1));
    else
      error_table_entry("MD.FinalTimeStep line corrupt\n");
    }

  //model->siesta.timestep
  if (g_ascii_strncasecmp("MD.LengthTimeStep", *buff, 17) == 0)
    {
    if (num_tokens >= 3)
      {
      model->siesta.timestep = str_to_float(*(buff+1));
      if (set_time_units(*(buff+2), &model->siesta.timestep_units) != 0)
        error_table_entry("*ERROR* - MD.LengthTimeStep units corrupt\n");
      }
    else if (num_tokens == 2)
      {
      error_table_entry("MD.LengthTimeStep units MISSING\n");
      model->siesta.timestep = str_to_float(*(buff+1));
      set_time_units("", &model->siesta.timestep_units);
      }
    }

  //model->siesta.md_inital_temperature
  if (g_ascii_strncasecmp("MD.InitialTemperature", *buff, 21) == 0)
    {
    if (num_tokens >= 3)
      {
      model->siesta.md_inital_temperature = str_to_float(*(buff+1));
      if (set_temperature_units(*(buff+2), &model->siesta.md_inital_temperature_units) != 0)
        error_table_entry("MD.InitialTemperature units corrupt\n");
      }
    else
      error_table_entry("MD.InitialTemperature line corrupt\n");
    }

  //model->siesta.md_target_temperature
  if (g_ascii_strncasecmp("MD.TargetTemperature", *buff, 20) == 0)
    {
    if (num_tokens >= 3)
      {
      model->siesta.md_target_temperature = str_to_float(*(buff+1));
      if (set_temperature_units(*(buff+2), &model->siesta.md_target_temperature_units) != 0)
        error_table_entry("MD.TargetTemperature units corrupt\n");
      }
    else
      error_table_entry("MD.TargetTemperature line corrupt\n");
    }

  //model->siesta.md_quench
  if (g_ascii_strncasecmp("MD.Quench", *buff, 9) == 0)
    {
    if (num_tokens >= 2)
      {
      if (set_true_false(*(buff+1), &model->siesta.md_quench) != 0)
        error_table_entry("MD.Quench boolean corrupt\n");
      }
    else
     error_table_entry("MD.Quench line corrupt\n");
    }

  if (g_ascii_strncasecmp("MD.TargetPressure", *buff, 20) == 0)
    {
    if (num_tokens >= 3)
      {
      model->siesta.md_target_pressure = str_to_float(*(buff+1));
      if (set_pressure_units(*(buff+2), &model->siesta.md_target_pressure_units) != 0)
        error_table_entry("MD.TargetPressure units corrupt\n");
      }
    else
      error_table_entry("MD.TargetPressure line corrupt\n");
    }

#if OLD_STYLE
  if (g_ascii_strncasecmp("%block", *buff, 6) == 0 && num_tokens > 1)
  {
    if  (g_ascii_strncasecmp("MD.TargetStress", *(buff+1), 15) == 0)
    {
      for (;;)
      {
        g_strfreev(buff);
        buff = get_tokenized_line(fp, &num_tokens);
        if (!buff)
           break;
        if (g_ascii_strncasecmp("%endblock", *buff, 9) == 0)
          break;

        if (num_tokens == 6)
        {
            model->siesta.md_target_stress_xx = str_to_float(*(buff));
            model->siesta.md_target_stress_yy = str_to_float(*(buff+1));
            model->siesta.md_target_stress_zz = str_to_float(*(buff+2));
            model->siesta.md_target_stress_xy = str_to_float(*(buff+3));
            model->siesta.md_target_stress_xz = str_to_float(*(buff+4));
            model->siesta.md_target_stress_yz = str_to_float(*(buff+5));
        }
        else
        {
            //invalid stress block
            text = error_table_entry("*ERROR* - MD.TargetStress block corrupt\n");
            gui_text_show(ERROR, text);
            g_free(text);
        }

      }
    }
  }
#endif

  if (g_ascii_strncasecmp("MD.NoseMass", *buff, 11) == 0)
    {
    if (num_tokens >= 3)
      {
      model->siesta.md_nose_mass = str_to_float(*(buff+1));
      if (set_mominert_units(*(buff+2), &model->siesta.md_nose_mass_units) != 0)
        error_table_entry("*ERROR* - MD.NoseMass units corrupt\n");
      }
    else
      error_table_entry("*ERROR* - MD.NoseMass line corrupt\n");
    }

  if (g_ascii_strncasecmp("MD.ParrinelloRahmanMass", *buff, 23) == 0)
    {
    if (num_tokens >= 3)
      {
      model->siesta.md_parrinello_rahman_mass = str_to_float(*(buff+1));
      if (set_mominert_units(*(buff+2), &model->siesta.md_parrinello_rahman_mass_units) != 0)
        error_table_entry("MD.ParrinelloRahmanMass units corrupt\n");
      }
    else
      error_table_entry("MD.ParrinelloRahmanMass line corrupt\n");
    }

  if (g_ascii_strncasecmp("MD.AnnealOption", *buff, 23) == 0)
    {
    if (num_tokens >= 2)
      {
      if (g_ascii_strncasecmp("Temperature", *(buff+1), 11) == 0)
        model->siesta.md_anneal_option = TEMPERATURE_ANNEAL;
      else if (g_ascii_strncasecmp("Pressure", *(buff+1), 8) == 0)
        model->siesta.md_anneal_option = PRESSURE_ANNEAL;
      else if (g_ascii_strncasecmp("TemperatureandPressure", *(buff+1), 22) == 0)
        model->siesta.md_anneal_option = TEMP_AND_PRESS_ANNEAL;
      else
       error_table_entry("MD.AnnealOption string corrupt\n");
      }
    else
      error_table_entry("MD.AnnealOption line corrupt\n");
    }

  if (g_ascii_strncasecmp("MD.TauRelax", *buff, 11) == 0)
    {
    if (num_tokens >= 3)
      {
      model->siesta.md_tau_relax = str_to_float(*(buff+1));
      if (set_time_units(*(buff+2), &model->siesta.md_tau_relax_units) != 0)
        error_table_entry("MD.TauRelax units corrupt\n");
      }
    else
      error_table_entry("MD.TauRelax line corrupt\n");
    }

  if (g_ascii_strncasecmp("MD.BulkModulus", *buff, 14) == 0)
    {
    if (num_tokens >= 3)
      {
      model->siesta.md_bulk_modulus = str_to_float(*(buff+1));
      if (set_pressure_units(*(buff+2), &model->siesta.md_bulk_modulus_units) != 0)
        error_table_entry("MD.BulkModulus units corrupt\n");
      }
    else
      error_table_entry("MD.BulkModulus line corrupt\n");
    }

  if (g_ascii_strncasecmp("MD.FCDispl", *buff, 10) == 0)
    {
    if (num_tokens >= 3)
      {
      model->siesta.md_fc_displ = str_to_float(*(buff+1));
      if (set_length_units(*(buff+2), &model->siesta.md_fc_displ_units) != 0)
        error_table_entry("MD.FCDispl units corrupt\n");
      }
    else
      error_table_entry("MD.FCDispl line corrupt\n");
    }

  if (g_ascii_strncasecmp("MD.FCfirst", *buff, 10) == 0)
    {
    if (num_tokens >= 2)
      model->siesta.md_fc_first = str_to_float(*(buff+1));
    else
      error_table_entry("MD.FCfirst line corrupt\n");
    }

  if (g_ascii_strncasecmp("MD.FClast", *buff, 9) == 0)
    {
    if (num_tokens >= 2)
      model->siesta.md_fc_last = str_to_float(*(buff+1));
    else
      error_table_entry("MD.FClast line corrupt\n");
    }

/* zmatrix */
  if (g_ascii_strncasecmp("zm.unitslength", *buff, 14) == 0)
    {
    if (num_tokens)
      {
      if (g_ascii_strncasecmp("ang", *(buff+1), 3) == 0)
        zmat_distance_units_set(model->zmatrix, ANGSTROM);
      else
        zmat_distance_units_set(model->zmatrix, BOHR);
      }
    }
  if (g_ascii_strncasecmp("zm.unitsangle", *buff, 13) == 0)
    {
    if (num_tokens)
      {
      if (g_ascii_strncasecmp("deg", *(buff+1), 3) == 0)
        zmat_angle_units_set(model->zmatrix, DEGREES); 
      else
        zmat_angle_units_set(model->zmatrix, RADIANS); 
      }
    }
    
    // exit condition - "************************** End of input data file"
    if (g_ascii_strncasecmp("**************************", *buff, 26) == 0)
    {
        if (g_ascii_strncasecmp("End", *(buff+1), 3) == 0)
        {
            //leave routine -> break the for loop
            g_strfreev(buff);
            free_slist(species_list);        
            return 1;
        }
    }

/* loop cleanup */
  siesta_done_line:;

  g_strfreev(buff);
  g_free(line);
  }

siesta_done_file:;


/* convert scaled cartesians to proper cartesians */
if (siesta_sucks)
  {
  GSList *list;
  struct core_pak *core;

  for (list=model->cores ; list ; list=g_slist_next(list))
    {
    core = list->data;
    VEC3MUL(core->x, scale_lattice);
    }
  }


/* NEW - process zmatrix cores */
zmat_type(model->zmatrix, species_list);
zmat_process(model->zmatrix, model);

/* free the species list */
free_slist(species_list);
  
return 0;
}

/**************************/
/* fdf input file reading */
/**************************/
#define DEBUG_READ_FDF 0
gint read_fdf(gchar *filename, struct model_pak *model)
{
gint i;
gchar *text;

FILE *fp;

/* checks */
g_return_val_if_fail(model != NULL, 1);
g_return_val_if_fail(filename != NULL, 2);

fp = fopen(filename, "rt");
if (!fp)
  return(3);

error_table_clear();

/* terry's mod */
read_fdf_block(fp, model, TRUE);

/* check cores */
model->cores = g_slist_reverse(model->cores);
i = g_slist_length(model->cores);
if (model->siesta.num_atoms != i)
  {
  text = g_strdup_printf("Inconsistancy reading %s: expected %d cores, but found %d.\n",
                          filename, model->siesta.num_atoms, i);
  gui_text_show(ERROR, text);
  g_free(text);
  }

/* model setup */
strcpy(model->filename, filename);

model_prep(model);

/* NB: this relies on having a valid latmat */
siesta_density_files(filename, model);

error_table_print_all();

return(0);
}

/*******************************************/
/* read single SIESTA output configuration */
/*******************************************/
#define DEBUG_READ_SOUT 0
gint read_sout_block(FILE *fp, struct model_pak *model)
{
gint i, num_tokens;
gchar **temp, **buff, line[LINELEN];
GString *title, *energy_string, *grad_string, *pressure_string;
GSList *clist;
struct core_pak *core;

clist = model->cores;

/* get 1st line of coords */
if (fgetline(fp, line))
  return(1);
buff = tokenize(line, &num_tokens);

while (num_tokens > 4)
  {
  core = NULL;

  if (clist)
    {
    core = clist->data;
    clist = g_slist_next(clist);
    }
  else
    {
/* NEW - siesta can change the column the element symbol appears in */
    for (i=3 ; i<num_tokens ; i++)
      {
      if (elem_symbol_test(*(buff+i)))
        {
        core = new_core(*(buff+i), model);
        model->cores = g_slist_append(model->cores, core);
        break;
        }
      }
    }

  if (core)
    {
    core->x[0] = str_to_float(*(buff+0));
    core->x[1] = str_to_float(*(buff+1));
    core->x[2] = str_to_float(*(buff+2));
    }

/* get next line */
  g_strfreev(buff);
  if (fgetline(fp, line))
    return(2);
  buff = tokenize(line, &num_tokens);
  }
g_strfreev(buff);

/* attempt to get the lattice matrix */
/* if not found, use the initial values */
/* removed the following line as if initial structure used LatticeVectors, this stopped the pbcs being constructed */
/* model->construct_pbc = FALSE;*/
while (!fgetline(fp, line))
  {
/* if next frame encoutered - assume there is no outcell data & exit */
  if (g_ascii_strncasecmp(line, "siesta:        Begin CG move", 28) == 0)
    break;
  /* new version of siesta doesn't have the siesta: */
  if (g_ascii_strncasecmp(line, "                            Begin CG move", 41) == 0)
    break;
  /* Read in MD steps from siesta: */
  if (g_ascii_strncasecmp(line, "                            Begin MD step", 41) == 0)
    break;

/* acquire cell lengths */
  if (g_ascii_strncasecmp(line, "outcell: Cell vector modules", 28) == 0)
    {
    model->construct_pbc = FALSE;
    temp = g_strsplit(line, ":", 3);
    buff = tokenize(*(temp+2), &num_tokens);
    g_strfreev(temp);

    if (num_tokens > 2)
      {
      model->pbc[0] = str_to_float(*(buff+0));
      model->pbc[1] = str_to_float(*(buff+1));
      model->pbc[2] = str_to_float(*(buff+2));
      }
    else
      printf("Unexpected data format reading cell lengths.\n");
    g_strfreev(buff);
    }

/* acquire cell angles */
  if (g_ascii_strncasecmp(line, "outcell: Cell angles", 20) == 0)
    {
    temp = g_strsplit(line, ":", 3);
    buff = tokenize(*(temp+2), &num_tokens);
    g_strfreev(temp);

    if (num_tokens > 2)
      {
/* get the angles (in radians) */
      model->pbc[3] = D2R*str_to_float(*(buff+0));
      model->pbc[4] = D2R*str_to_float(*(buff+1));
      model->pbc[5] = D2R*str_to_float(*(buff+2));
      }
    else
      printf("Unexpected data format reading cell angles.\n");
    g_strfreev(buff);

/* no more data - break out */
/*
    break;
*/
    }

  if (g_ascii_strncasecmp(line, "siesta: E_KS(eV)", 16) == 0)
    {
    buff = tokenize(line, &num_tokens);
    model->siesta.energy = str_to_float(*(buff+3));
    model->siesta.have_energy = TRUE;
    g_strfreev(buff);
    }

  if (g_strrstr(line, "constrained") != NULL)
    {
    buff = tokenize(line, &num_tokens);
    model->siesta.max_grad = str_to_float(*(buff+1));
    model->siesta.have_max_grad = TRUE;
    g_strfreev(buff);
    }

  /* FIXME This only works at the end of an optimization. Need to fix for all frames */
  if (g_ascii_strncasecmp(line, "siesta: Pressure (static)", 25) == 0)
    {
    for (i=0; i<4; i++)
      {
      if (fgetline(fp, line))
        return(2);
      }
    buff = tokenize(line, &num_tokens);
    pressure_string = g_string_new("");
    g_string_append_printf(pressure_string, "%.4f %s", str_to_float(*(buff+1)), *(buff+3));
    property_add_ranked(5, "Final Pressure", pressure_string->str, model);
    g_string_free(pressure_string, TRUE);
    }
  }

g_free(model->title);

/* TODO read in units from sot file */
title = g_string_new("");
if (model->siesta.have_energy)
  {
  energy_string = g_string_new("");
  g_string_append_printf(energy_string, "%.5f eV", model->siesta.energy);
  property_add_ranked(3, "Energy", energy_string->str, model);
  g_string_free(energy_string, TRUE);
  g_string_append_printf(title, "E");
  g_string_append_printf(title, " = %.4f eV, ", model->siesta.energy);
  }
if (model->siesta.have_max_grad)
  {
  grad_string = g_string_new("");
  g_string_append_printf(grad_string, "%.4f eV/A", model->siesta.max_grad);
  property_add_ranked(4, "Maximum Gradient", grad_string->str, model);
  g_string_free(grad_string, TRUE);
  g_string_append_printf(title, "max grad = %.5f", model->siesta.max_grad);
  }
model->title = g_strdup(title->str);
g_string_free(title, TRUE);

return(0);
}

/*******************************/
/* SIESTA output frame reading */
/*******************************/
gint read_sout_frame(FILE *fp, struct model_pak *model)
{
    return(read_sout_block(fp, model));
}

/********************************/
/* Read in a SIESTA output file */
/********************************/
/* could almost use the read_fdf() code, as SIESTA spits out a copy */
/* however, we do need to know the number of frames for animation... */
gint read_sout(gchar *filename, struct model_pak *model)
{
gint frame=0;
gchar **buff, line[LINELEN], *text;
FILE *fp;
    
gint num_tokens, i;
gdouble scale = 1.0;

fp = fopen(filename, "rt");
if (!fp)
  return(1);

fgetline(fp, line);
/*
while (g_ascii_strncasecmp("************************** Dump of input data file", line, 50) != 0)
    fgetline(fp, line);
*/

while (!strstr(line, "Dump of input data file"))
    fgetline(fp, line);

//should be at the input file params -- dont want inital co-ords.
//read_fdf_block(fp,model, TRUE);
//should magicly return at the ======


while (!fgetline(fp, line))
  {
//---------------START TERRY COMMENT DUE TO read_fdf_pak_filler()----------      

/* lattice constant scaling */
  if (g_ascii_strncasecmp("LatticeConstant", line, 15) == 0)
    {
    scale = siesta_scale_parse(line);
    continue;
    }

/* default cell dimensions */
  if (g_ascii_strncasecmp("%block LatticeParameters", line, 24) == 0)
  {
     if (fgetline(fp, line))
     return(2);
     buff = tokenize(line, &num_tokens);
     if (num_tokens > 5)
      {
      model->pbc[0] = scale*str_to_float(*(buff+0));
      model->pbc[1] = scale*str_to_float(*(buff+1));
      model->pbc[2] = scale*str_to_float(*(buff+2));
      model->pbc[3] = D2R*str_to_float(*(buff+3));
      model->pbc[4] = D2R*str_to_float(*(buff+4));
      model->pbc[5] = D2R*str_to_float(*(buff+5));
      model->construct_pbc = FALSE;
      model->periodic = 3;
      }
    g_strfreev(buff);
  }
  else if (g_ascii_strncasecmp("%block LatticeVectors", line, 20) == 0)
  {
      for (i=0; i<3; i++)
      {
        if (fgetline(fp, line))
          return(2);
        buff = tokenize(line, &num_tokens);
        if (num_tokens == 3)
        {
          model->latmat[0+i] = scale*str_to_float(*(buff+0));
          model->latmat[3+i] = scale*str_to_float(*(buff+1));
          model->latmat[6+i] = scale*str_to_float(*(buff+2));
        }
        g_strfreev(buff);
        model->construct_pbc = TRUE;
        model->periodic = 3;
      }
  }

//---------------END TERRY COMMENT DUE TO read_fdf_pak_filler()----------     
 

/* functional */
  if (g_ascii_strncasecmp(line, "xc_check: Exchange-correlation functional:", 42) == 0)
    {
    if (fgetline(fp, line))
      return(2);
    if (g_ascii_strncasecmp(line, "xc_check: GGA Perdew, Burke & Ernzerhof 1996", 44) == 0)
      property_add_ranked(7, "Functional", "PBE", model);
    else if (g_ascii_strncasecmp(line, "xc_check: GGA Becke Lee Yang Parr", 33) == 0)
      property_add_ranked(7, "Functional", "BLYP", model);
    else
      property_add_ranked(7, "Functional", "Unknown", model);
    }

/* k-points */
  if (g_ascii_strncasecmp(line, "siesta: k-grid: Number of k-points", 34) == 0)
    {
    buff = g_strsplit(line, "=", 2);
    g_strchug(g_strchomp(*(buff+1)));
    property_add_ranked(10, "K-points", *(buff+1), model);
    g_strfreev(buff);
    }

/* Mesh cutoff */
  if (g_ascii_strncasecmp(line, "redata: Mesh Cutoff", 19) == 0)
  {
    buff = g_strsplit(line, "=", 2);
    text = format_value_and_units(*(buff+1), 2);
    property_add_ranked(6, "Mesh Cutoff", text, model);
    g_free(text);
    g_strfreev(buff);
  }

/* Energy Shift */
  if (g_ascii_strncasecmp(line, "SPLIT: energy shift",  19) == 0)
  {
    buff = g_strsplit(line, "=", 2);
    text = format_value_and_units(*(buff+1), 6);
    property_add_ranked(6, "Energy Shift", text, model);
    g_free(text);
    g_strfreev(buff);
  }

/* coordinates */
  if (model->siesta.md_type_of_run != FC_MDRUN)
  {
    if (g_ascii_strncasecmp(line, "outcoor: ", 9) == 0)
      {
        if (g_strrstr(line, "Ang") != NULL)
          model->fractional = FALSE;
        else if (g_strrstr(line, "fractional") != NULL)
          model->fractional = TRUE;
        else if (g_strrstr(line, "Bohr") != NULL)
          model->coord_units = BOHR;
        else
          {
          gui_text_show(ERROR, "unexpected coordinate type\n");
          return(2);
          }
        add_frame_offset(fp, model);
        read_sout_block(fp, model);
        frame++;
      }
  }


/* seems to increment the number of frames but doesn't set the animation pointer! */
/* line only occurs at end of a run that runs out of cycles so commenting out! */
/*  if (g_ascii_strncasecmp(line, "outcoor: Final", 14) == 0)
    frame++;*/
  }
  
  
  if (model->siesta.md_type_of_run == FC_MDRUN)
  {
    //something went wrong
    error_table_entry("MD.TypeOfRun == FC -> so inital co-ords from input shown\n");
  }
  
  

/* done */
strcpy(model->filename, filename);
g_free(model->basename);
model->basename = parse_strip(filename);
model->num_frames = model->cur_frame = frame;
model->cur_frame--;

model_prep(model);

return(0);
}

gint set_energy_units(gchar * buff, enum siesta_energy_measurement_type * siesta_unit_type)
{
    gint allok = 0;
    if (g_ascii_strncasecmp("J", buff, 1) == 0)
        *siesta_unit_type = J_ENRG;
    else if (g_ascii_strncasecmp("erg", buff, 3) == 0)
        *siesta_unit_type = ERG_ENRG;
    else if (g_ascii_strncasecmp("eV", buff, 2) == 0)
        *siesta_unit_type = EV_ENRG;
    else if (g_ascii_strncasecmp("meV", buff, 3) == 0)
        *siesta_unit_type = MEV_ENRG;
    else if (g_ascii_strncasecmp("Ry", buff, 2) == 0)
        *siesta_unit_type = RY_ENRG;
    else if (g_ascii_strncasecmp("mRy", buff, 3) == 0)
        *siesta_unit_type = MRY_ENRG;
    else if (g_ascii_strncasecmp("Hartree", buff, 7) == 0)
        *siesta_unit_type = HARTREE_ENRG;
    else if (g_ascii_strncasecmp("K", buff, 1) == 0)
        *siesta_unit_type = KO_ENRG;
    else if (g_ascii_strncasecmp("kcal/mol", buff, 8) == 0)
        *siesta_unit_type = KCALMOLE_ENRG;
    else if (g_ascii_strncasecmp("mHartree", buff, 8) == 0)
        *siesta_unit_type = MHARTREE_ENRG;
    else if (g_ascii_strncasecmp("kJ/mol", buff, 6) == 0)
        *siesta_unit_type = KJMOL_ENRG;
    else if (g_ascii_strncasecmp("Hz", buff, 3) == 0)
        *siesta_unit_type = HZ_ENRG;
    else if (g_ascii_strncasecmp("THz", buff, 3) == 0)
        *siesta_unit_type = THZ_ENRG;
    else if (g_ascii_strncasecmp("cm-1", buff, 4) == 0)
        *siesta_unit_type = CM_ENRG;
    else if (g_ascii_strncasecmp("cm**-1", buff, 6) == 0)
        *siesta_unit_type = CMCM_ENRG;
    else
    {
        //no match - set default
        *siesta_unit_type = RY_ENRG;
        allok = 1;
    }
    return allok;
}
gint set_temperature_units(gchar * buff, enum siesta_temperature_measurement_type * siesta_unit_type)
{
    gint allok = 0;
    if (g_ascii_strncasecmp("K", buff, 1) == 0)
        *siesta_unit_type = K_TEMP;
    else if (g_ascii_strncasecmp("C", buff, 1) == 0)
        *siesta_unit_type = C_TEMP;
    else
    {
        //no match - set default
        *siesta_unit_type = K_TEMP;
        allok = 1;
    }
    return allok;
}

gint set_mass_units(gchar * buff, enum siesta_mass_measurement_type * siesta_unit_type)
{
    //KG_MASS, G_MASS, AMU_MASS
    //
    gint allok = 0;
    if (g_ascii_strncasecmp("Kg", buff, 2) == 0)
        *siesta_unit_type = KG_MASS;
    else if (g_ascii_strncasecmp("g", buff, 1) == 0)
        *siesta_unit_type = G_MASS;
    else if (g_ascii_strncasecmp("amu", buff, 3) == 0)
        *siesta_unit_type = AMU_MASS;
    else
    {
        //no match - set default
        *siesta_unit_type = AMU_MASS;
        allok = 1;
    }
    return allok;
}

gint set_length_units(gchar * buff, enum siesta_length_measurement_type * siesta_unit_type)
{
    //M_LEN, CM_LEN, NM_LEN, ANG_LEN, BOHR_LEN
    gint allok = 0;
    if (g_ascii_strncasecmp("m", buff, 1) == 0)
        *siesta_unit_type = M_LEN;
    else if (g_ascii_strncasecmp("cm", buff, 2) == 0)
        *siesta_unit_type = CM_LEN;
    else if (g_ascii_strncasecmp("nm", buff, 2) == 0)
        *siesta_unit_type = NM_LEN;
    else if (g_ascii_strncasecmp("Ang", buff, 3) == 0)
        *siesta_unit_type = ANG_LEN;
    else if (g_ascii_strncasecmp("Bohr", buff, 4) == 0)
        *siesta_unit_type = BOHR_LEN;
    else
    {
        //no match - set default
        *siesta_unit_type = BOHR_LEN;
        allok = 1;
    }
    return allok;
}

gint set_time_units(gchar * buff, enum siesta_time_measurement_type * siesta_unit_type)
{
    //S_TIME, FS_TIME, PS_TIME, NS_TIME
    gint allok = 0;
    if (g_ascii_strncasecmp("S", buff, 1) == 0)
        *siesta_unit_type = S_TIME;
    else if (g_ascii_strncasecmp("fs", buff, 2) == 0)
        *siesta_unit_type = FS_TIME;
    else if (g_ascii_strncasecmp("ps", buff, 2) == 0)
        *siesta_unit_type = PS_TIME;
    else if (g_ascii_strncasecmp("ns", buff, 2) == 0)
        *siesta_unit_type = NS_TIME;
    else
    {
        //no match - set default
        *siesta_unit_type = FS_TIME;
        allok = 1;
    }
    return allok;
}

gint set_pressure_units(gchar * buff, enum siesta_pressure_measurement_type * siesta_unit_type)
{
    // P_PRES, MPA_PRES, GPA_PRES, ATM_PRES, BAR_PRES, MBAR_PRES, RYBOHR3_PRES, EVANG3_PRES
    gint allok = 0;
    if (g_ascii_strncasecmp("Pa", buff, 2) == 0)
        *siesta_unit_type = P_PRES;
    else if (g_ascii_strncasecmp("MPa", buff, 3) == 0)
        *siesta_unit_type = MPA_PRES;
    else if (g_ascii_strncasecmp("GPa", buff, 3) == 0)
        *siesta_unit_type = GPA_PRES;
    else if (g_ascii_strncasecmp("atm", buff, 3) == 0)
        *siesta_unit_type = ATM_PRES;
    else if (g_ascii_strncasecmp("bar", buff, 3) == 0)
        *siesta_unit_type = BAR_PRES;
    else if (g_ascii_strncasecmp("Mbar", buff, 4) == 0)
        *siesta_unit_type = MBAR_PRES;
    else if (g_ascii_strncasecmp("Ry/Bohr**3", buff, 10) == 0)
        *siesta_unit_type = RYBOHR3_PRES;
    else if (g_ascii_strncasecmp("eV/Ang**3", buff, 9) == 0)
        *siesta_unit_type = EVANG3_PRES;
    else
    {
        //no match - set default
        *siesta_unit_type = P_PRES;
        allok = 1;
    }
    return allok;
}

gint set_force_units(gchar * buff, enum siesta_force_measurement_type * siesta_unit_type)
{
    //N_FORCE, EVANG_FORCE, RYBOHR_FORCE
    gint allok = 0;
    if (g_ascii_strncasecmp("N", buff, 1) == 0)
        *siesta_unit_type = N_FORCE;
    else if (g_ascii_strncasecmp("eV/Ang", buff, 6) == 0)
        *siesta_unit_type = EVANG_FORCE;
    else if (g_ascii_strncasecmp("Ry/Bohr", buff, 7) == 0)
        *siesta_unit_type = RYBOHR_FORCE;
    else
    {
        //no match - set default
        *siesta_unit_type = N_FORCE;
        allok = 1;
    }
    return allok;
}

gint set_mominert_units(gchar * buff, enum siesta_mominert_measurement_type * siesta_unit_type)
{
    //KGM_MOMINERT, RYFS_MOMINERT
    gint allok = 0;
    if (g_ascii_strncasecmp("Kg*m**2", buff, 7) == 0)
        *siesta_unit_type = KGM_MOMINERT;
    else if (g_ascii_strncasecmp("Ry*fs**2", buff, 8) == 0)
        *siesta_unit_type = RYFS_MOMINERT;
    else
    {
        //no match - set default
        *siesta_unit_type = KGM_MOMINERT;
        allok = 1;
    }
    return allok;
}


gint set_true_false(gchar * buff, gboolean * siesta_boolean)
{
    gint allok = 0;
    if (g_ascii_strncasecmp("true", buff, 4) == 0)
        *siesta_boolean = TRUE;
    else if (g_ascii_strncasecmp(".true", buff, 5) == 0)
        *siesta_boolean = TRUE;
    else if (g_ascii_strncasecmp(".false", buff, 6) == 0)
        *siesta_boolean = FALSE;
    else if (g_ascii_strncasecmp("false", buff, 5) == 0)
        *siesta_boolean = FALSE;
    else {
        //no match - set default
        *siesta_boolean = FALSE;
        allok = 1;
    }
    return allok;
}



gint print_file_energy_units(FILE * fp, enum siesta_energy_measurement_type units)
{
    gint allok = 0;
    switch (units)
    {
        case J_ENRG:            fprintf(fp, " J\n"); break;
        case ERG_ENRG:          fprintf(fp, " erg\n"); break;
        case EV_ENRG:           fprintf(fp, " eV\n"); break;
        case MEV_ENRG:          fprintf(fp, " meV\n"); break;
        case RY_ENRG:           fprintf(fp, " Ry\n"); break;
        case MRY_ENRG:          fprintf(fp, " mRy\n"); break;
        case HARTREE_ENRG:      fprintf(fp, " Hartree\n"); break;
        case KO_ENRG:           fprintf(fp, " K\n"); break;
        case KCALMOLE_ENRG:     fprintf(fp, " kcal/mol\n"); break;
        case MHARTREE_ENRG:     fprintf(fp, " mHartree\n"); break;
        case KJMOL_ENRG:        fprintf(fp, " kJ/mol\n"); break;
        case HZ_ENRG:           fprintf(fp, " Hz\n"); break;
        case THZ_ENRG:          fprintf(fp, " THz\n"); break;
        case CM_ENRG:           fprintf(fp, " cm-1\n"); break;
        case CMCM_ENRG:         fprintf(fp, " cm**-1\n"); break;
        default:
            //no units given?
            //print new line
            fprintf(fp," \n");
            allok = 1;
            break;
    }
    return allok;
}

gint print_file_mass_units(FILE * fp, enum siesta_mass_measurement_type units)
{
    gint allok = 0;
    switch (units)
    {
        case KG_MASS:            fprintf(fp, " Kg\n"); break;
        case G_MASS:             fprintf(fp, " g\n"); break;
        case AMU_MASS:           fprintf(fp, " amu\n"); break;
        default:
            //no units given? print new line
            fprintf(fp," \n");
            allok = 1;
            break;
    }
    return allok;
}

gint print_file_length_units(FILE * fp, enum siesta_length_measurement_type units)
{
    gint allok = 0;
    switch (units)  //M_LEN, CM_LEN, NM_LEN, ANG_LEN, BOHR_LEN
    {
        case M_LEN:             fprintf(fp, " m\n"); break;
        case CM_LEN:            fprintf(fp, " cm\n"); break;
        case NM_LEN:            fprintf(fp, " nm\n"); break;
        case ANG_LEN:           fprintf(fp, " Ang\n"); break;
        case BOHR_LEN:          fprintf(fp, " Bohr\n"); break;
        default:
            //no units given? print new line
            fprintf(fp," \n");
            allok = 1;
            break;
    }
    return allok;
}

gint print_file_time_units(FILE * fp, enum siesta_time_measurement_type units)
{
    gint allok = 0;
    switch (units)  //S_TIME, FS_TIME, PS_TIME, NS_TIME
    {
        case S_TIME:            fprintf(fp, " s\n"); break;
        case FS_TIME:           fprintf(fp, " fs\n"); break;
        case PS_TIME:           fprintf(fp, " ps\n"); break;
        case NS_TIME:           fprintf(fp, " ns\n"); break;
        default:
            //no units given? print new line or force units?
            fprintf(fp," \n");
            allok = 1;
            break;
    }
    return allok;
}

//Pressure
gint print_file_pressure_units(FILE * fp, enum siesta_pressure_measurement_type units)
{
    gint allok = 0;
    switch (units)
    {
        case P_PRES:                fprintf(fp, " Pa\n"); break;
        case MPA_PRES:              fprintf(fp, " MPa\n"); break;
        case GPA_PRES:              fprintf(fp, " GPa\n"); break;
        case ATM_PRES:              fprintf(fp, " atm\n"); break;
        case BAR_PRES:              fprintf(fp, " bar\n"); break;
        case MBAR_PRES:             fprintf(fp, " Mbar\n"); break;
        case RYBOHR3_PRES:          fprintf(fp, " Ry/Bohr**3\n"); break;
        case EVANG3_PRES:           fprintf(fp, " eV/Ang**3\n"); break;
        default:
            //no units given? print new line or force units?
            fprintf(fp," \n");
            allok = 1;
            break;
    }
    return allok;
}

gint print_file_temperature_units(FILE * fp, enum siesta_temperature_measurement_type units)
{
    gint allok = 0;
    switch (units)
    {
        case C_TEMP:                fprintf(fp, " C\n"); break;
        case K_TEMP:                fprintf(fp, " K\n"); break;
        default:
            //no units given? print new line or force units?
            fprintf(fp," \n");
            allok = 1;
            break;
    }
    return allok;
}

gint print_file_force_units(FILE * fp, enum siesta_force_measurement_type units)
{
    gint allok = 0;
    switch (units) //N_FORCE, EVANG_FORCE, RYBOHR_FORCE
    {
        case N_FORCE:               fprintf(fp, " N\n"); break;
        case EVANG_FORCE:           fprintf(fp, " eV/Ang\n"); break;
        case RYBOHR_FORCE:          fprintf(fp, " Ry/Bohr\n"); break;
        default:
            //no units given? print new line or force units?
            fprintf(fp," \n");
            allok = 1;
            break;
    }
    return allok;
}

gint print_file_mominert_units(FILE * fp, enum siesta_mominert_measurement_type units)
{
    gint allok = 0;
    switch (units) //KGM_MOMINERT, RYFS_MOMINERT
    {
        case KGM_MOMINERT:            fprintf(fp, " Kg*m**2\n"); break;
        case RYFS_MOMINERT:           fprintf(fp, " Ry*fs**2\n"); break;
        default:
            //no units given? print new line or force units?
            fprintf(fp," \n");
            allok = 1;
            break;
    }
    return allok;
}
