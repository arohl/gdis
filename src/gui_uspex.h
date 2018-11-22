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
#define USPEX_PAGE_SYSTEM 0
#define USPEX_PAGE_STRUCTURES 1
#define USPEX_PAGE_CALCULATION 2
#define USPEX_PAGE_ADVANCED 3
#define USPEX_PAGE_SPECIFIC 4
/*fixed parameters*/
#define USPEX_MAX_NUM_OPT_STEPS 16
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
	GUI_OBJ *specific_page;/*because this page can be locked if nothing relevant*/
/*4.1 Type of run & System*/
	GUI_OBJ *calculationMethod;
	GUI_OBJ *calculationType;
	GUI_OBJ *_calctype_dim;
	gdouble _dim;/*absurd spin on double*/
	GUI_OBJ *_calctype_mol;
	GUI_OBJ *_calctype_var;
	GUI_OBJ *_calctype_mag;			/*VER 10.1*/
	GUI_OBJ *_calctype_mag_2;		/*VER 10.1*/
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
	/*numSpecies block*/
	GUI_OBJ *numSpecies;
	gchar *_tmp_blockSpecies;
	GUI_OBJ *blockSpecies;
	GUI_OBJ *Species_apply_button;
	GUI_OBJ *Species_delete_button;
	/*bond definition: apply/auto*/
	GUI_OBJ *goodBonds;
	GUI_OBJ *_bond_d;
	gchar *_tmp_bond_d;
	gboolean auto_bonds;
	//checkMolecules is auto-sync
	//checkConnectivity is auto-sync
	GUI_OBJ *fitLimit;			/*VER 10.1*/
	gchar *_tmp_ldaU;			/*VER 10.1*/
	GUI_OBJ *ldaU;				/*VER 10.1*/
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
	GUI_OBJ *minVectorLength;
	GUI_OBJ *MolCenters;
	GUI_OBJ *_centers;
	GUI_OBJ *_centers_button;
	GUI_OBJ *constraint_enhancement;
	gboolean auto_C_ion;
/*4.6 Cell*/
	GUI_OBJ *Latticevalues;
	GUI_OBJ *_latticeformat;
	GUI_OBJ *_latticevalue;
	gchar *_tmp_latticevalue;
	GUI_OBJ *splitInto;
	gboolean auto_C_lat;
/*4.7 restart*/
	GUI_OBJ *pickUpYN;
	GUI_OBJ *pickUpGen;
	GUI_OBJ *pickUpFolder;
	gboolean restart_cleanup;
/*4.8 Ab initio*/
	gdouble _tmp_num_opt_steps;/*because gtk_spin type is double*/
	GUI_OBJ *_num_opt_steps;/*usually determined indirectly*/
	gdouble _tmp_curr_step;/*because gtk_spin type is double*/
	GUI_OBJ *_curr_step;
	gboolean _tmp_isfixed;
	GUI_OBJ *_isfixed;
	GUI_OBJ *abinitioCode;
	GUI_OBJ *KresolStart;
	GUI_OBJ *vacuumSize;
	gboolean auto_step;
	gchar **_tmp_ai_input;
	GUI_OBJ *ai_input;
	GUI_OBJ *ai_input_button;
	gchar **_tmp_ai_opt;
	GUI_OBJ *ai_opt;
	GUI_OBJ *ai_opt_button;
	GUI_OBJ *ai_generate;
	gchar *_tmp_ai_spe;
	GUI_OBJ *ai_spe;
	gchar **_potentials;
	gchar *_tmp_ai_pot;
	GUI_OBJ *ai_pot;
	GUI_OBJ *ai_pot_button;
	GUI_OBJ *numProcessors;
	GUI_OBJ *numParallelCalcs;
	gchar **_tmp_commandExecutable;
	GUI_OBJ *commandExecutable;
	GUI_OBJ *whichCluster;
	GUI_OBJ *remoteFolder;
	//PhaseDiagram is auto-sync
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
/*4.12 developers*/
	GUI_OBJ *repeatForStatistics;
	GUI_OBJ *stopFitness;
	GUI_OBJ *fixRndSeed;
	//collectForces is auto-sync
/*4.13 seldom used*/
	//ordering_active is auto-sync
	//symmetrize is auto-sync
	GUI_OBJ *valenceElectr;
	GUI_OBJ *percSliceShift;
	GUI_OBJ *dynamicalBestHM;
	GUI_OBJ *softMutOnly;
	GUI_OBJ *maxDistHeredity;
	GUI_OBJ *manyParents;
	GUI_OBJ *minSlice;
	GUI_OBJ *maxSlice;
	GUI_OBJ *numberparents;
/*5.1 molecular: ADDITIONAL*/
	GUI_OBJ *mol_model;
	GUI_OBJ *mol_model_button;
	gdouble _tmp_num_mol;
	GUI_OBJ *num_mol;
	GUI_OBJ *mol_gdis;
	gdouble _tmp_curr_mol;
	GUI_OBJ *curr_mol;
	gboolean mol_as_gulp;
	GUI_OBJ *mol_gulp;
	gint *_tmp_mols_gdis;
	gboolean *_tmp_mols_gulp;
	GUI_OBJ *mol_apply_button;
/*5.2 BoltzTraP*/
	gboolean have_ZT;
	GUI_OBJ *BoltzTraP_T_max;
	GUI_OBJ *BoltzTraP_T_delta;
	GUI_OBJ *BoltzTraP_T_efcut;
	GUI_OBJ *TE_T_interest;
	GUI_OBJ *TE_threshold;
	GUI_OBJ *TE_goal;
	gchar *_tmp_cmd_BoltzTraP;/*optional, and unused (for now)*/
	GUI_OBJ *cmd_BoltzTraP;/*optional, and unused (for now)*/
	GUI_OBJ *cmd_BoltzTraP_button;
/*5.3 Surfaces*/
	GUI_OBJ *thicknessS;
	GUI_OBJ *thicknessB;
	GUI_OBJ *reconstruct;
	GUI_OBJ *StoichiometryStart;/*almost undocumented*/
	GUI_OBJ *substrate_model;/*additional*/
	GUI_OBJ *substrate_model_button;/*additional*/
/*5.4 Clusters*/
/*5.5 variable composition*/
	GUI_OBJ *firstGeneMax;
	GUI_OBJ *minAt;
	GUI_OBJ *maxAt;
	GUI_OBJ *fracTrans;
	GUI_OBJ *howManyTrans;
	GUI_OBJ *specificTrans;
/*5.6 metadynamics*/
	GUI_OBJ *ExternalPressure;
	GUI_OBJ *GaussianWidth;
	GUI_OBJ *GaussianHeight;
	GUI_OBJ *FullRelax;
	GUI_OBJ *maxVectorLength;
	GUI_OBJ *meta_model;/*additional*/
	GUI_OBJ *meta_model_button;/*additional*/
/*5.7 Particle swarn optimization*/
	GUI_OBJ *PSO_softMut;
	GUI_OBJ *PSO_BestStruc;
	GUI_OBJ *PSO_BestEver;
/*6.1 Variable-Cell Nudged Elastic Band */
	GUI_OBJ *vcnebType;
	GUI_OBJ *_vcnebtype_method;
	GUI_OBJ *_vcnebtype_img_num;
	GUI_OBJ *_vcnebtype_spring;
	GUI_OBJ *numImages;
	GUI_OBJ *numSteps;
	GUI_OBJ *optReadImages;
	GUI_OBJ *optimizerType;
	GUI_OBJ *optRelaxType;
	GUI_OBJ *dt;
	GUI_OBJ *ConvThreshold;
	GUI_OBJ *VarPathLength;
	GUI_OBJ *K_min;
	GUI_OBJ *K_max;
	GUI_OBJ *Kconstant;
	GUI_OBJ *optFreezing;
	GUI_OBJ *optMethodCIDI;
	GUI_OBJ *startCIDIStep;
	GUI_OBJ *pickupImages;
	GUI_OBJ *FormatType;
	GUI_OBJ *PrintStep;
	GUI_OBJ *img_model;/*additional*/
	GUI_OBJ *img_model_button;/*additional*/
/*6.2 Transition Path Sampling*/
	GUI_OBJ *numIterations;
	GUI_OBJ *speciesSymbol;
	GUI_OBJ *mass;
	GUI_OBJ *amplitudeShoot_AB;
	GUI_OBJ *amplitudeShoot_BA;
	GUI_OBJ *magnitudeShoot_success;
	GUI_OBJ *magnitudeShoot_failure;
	GUI_OBJ *shiftRatio;
	GUI_OBJ *orderParaType;
	GUI_OBJ *opCriteria_start;
	GUI_OBJ *opCriteria_end;
	GUI_OBJ *cmdOrderParameter;
	GUI_OBJ *cmdOrderParameter_button;
	GUI_OBJ *cmdEnthalpyTemperature;
	GUI_OBJ *cmdEnthalpyTemperature_button;
	GUI_OBJ *orderParameterFile;
	GUI_OBJ *orderParameterFile_button;
	GUI_OBJ *enthalpyTemperatureFile;
	GUI_OBJ *enthalpyTemperatureFile_button;
	GUI_OBJ *trajectoryFile;
	GUI_OBJ *trajectoryFile_button;
	GUI_OBJ *MDrestartFile;
	GUI_OBJ *MDrestartFile_button;
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
