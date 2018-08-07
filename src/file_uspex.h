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

/* header for USPEX file */

/*defines*/
#define _UO (*uspex_output)
#define _UC (*uspex_calc)



/*structures*/
typedef enum {/*released USPEX calculation methods as of 2018 (5)*/
	US_CM_USPEX,	/*USPEX*/
	US_CM_META,	/*metadynamics*/
	US_CM_VCNEB,	/*VC-NEB*/
	US_CM_PSO,	/*PSO, unsupported*/
	US_CM_UNKNOWN,	/*reserved for future use*/
} uspex_method;

typedef enum {
	US_OT_ENTHALPY=1,			/*most stable structure*/
	US_OT_VOLUME=2,				/*min volume*/
	US_OT_HARDNESS=3,			/*max hardness*/
	US_OT_ORDER=4,				/*most ordered structure*/
	US_OT_DISTANCE=5,			/*max structure <differences>*/
	US_OT_DIELEC_S=6,			/*max dielectric suseptibility*/
	US_OT_GAP=7,				/*max band gap*/
	US_OT_DIELEC_GAP=8,			/*max energy storage capacity*/
	US_OT_MAG=9,				/*max magnetization*/
	US_OT_QE=10,				/*max stucture quasientropy*/
/*elasticity related*/
	US_OT_BULK_M=1101,			/*max bulk modulus*/
	US_OT_SHEAR_M=1102,			/*max shear modulus*/
	US_OT_YOUNG_M=1103,			/*max Young modulus*/
	US_OT_POISSON=1104,			/*max Poisson modulus*/
	US_OT_PUGH_R=1105,			/*max Pugh modulus ratio*/
	US_OT_VICKERS_H=1106,			/*max Vickers hardness*/
	US_OT_FRACTURE=1107,			/*max fracture toughness*/
	US_OT_DEBYE_T=1108,			/*max Debye temperature*/
	US_OT_SOUND_V=1109,			/*max sound velocity*/
	US_OT_SWAVE_V=1110,			/*max S-wave velocity*/
	US_OT_PWAVE_V=1111,			/*max P-wave velocity*/
	US_OT_UNKNOWN=0,			/*reserved for future use*/
} uspex_opt;

typedef struct {
        gint gen;				/*generation of this individual*/
	gint natoms;				/*number of atoms*/
	gint *atoms;				/*number of atoms per species*/
        gdouble energy;				/*total energy*/
	gdouble E;				/*reduce energy / atoms*/
	gdouble fitness;			/*NEW: fitness*/
        gdouble volume;				/*total volume*/
} uspex_individual;

/* USPEX calculation structure */
typedef struct {
/*additional controls*/
	gchar *name;				/*the job name*/
	gint special;				/*special tag*/
/*4.1 Type of run & System*/
        uspex_method    calculationMethod;	/*supported calculation methods*/
        gint            calculationType;	/*supported calculation types*/
        gint            _calctype_dim;		/*HIDDEN: X.. digit in calculation type: dimension*/
        gboolean        _calctype_mol;		/*HIDDEN: .X. digit in calculation type: molecular*/
        gboolean        _calctype_var;		/*HIDDEN: ..X digit in calculation type: variable composition*/
        uspex_opt	optType;		/* optimization property*/
	gint		_nspecies;		/*HIDDEN: number of different species*/
	gint 		*atomType;		/*atomic number of each species*/
	gint 		_var_nspecies;		/*HIDDEN: numer of numSpecies blocks*/
	gint 		*numSpecies;		/*number of species for each type*/
	gdouble		ExternalPressure;	/*external pressure*/
	gint 		*valences;		/*valence for each species*/
	gdouble		*goodBonds;		/*minimum bond valences*/
	gboolean	checkMolecules;		/*check molecules*/
	gboolean	checkConnectivity;	/*check connectivity*/
/*4.2 Population*/
	gint 		populationSize;		/*number of structure in each generation*/
	gint 		initialPopSize;		/*number of structure in initial generation*/
	gint 		numGenerations;		/*maximum number of generation*/
	gint 		stopCrit;		/*stop criterion (max generation with same lowest structure)*/
/*4.3 Survival of the fittest & Selection*/
	gdouble 	bestFrac;		/*fraction of structures in current generation used for producing the next*/
	gint 		keepBestHM;		/*number of surviving best structures*/
	gboolean	reoptOld;		/*reoptimization of survivors*/
/*4.4 Structure generation and variation operators*/
	gchar		*symmetries;		/*possible symmetries*/
	gdouble		fracGene;		/*fraction of structures generated by heredity*/
	gdouble		fracRand;		/*fraction of structures generated randomly*/
	gdouble		fracPerm;		/*fraction of structures generated by permutation*/
	gdouble		fracAtomsMut;		/*fraction of structures generated by softmutation*/
	gdouble		fracRotMut;		/*fraction of structures generated by orientation mutation*/
	gdouble		fracLatMut;		/*fraction of structures generated by lattice mutation*/
	gint		howManySwaps;		/*number of pairwise swaps*/
	gchar 		*specificSwaps;		/*which atoms are swapped*/
	gdouble		mutationDegree;		/*max displacement (Ang) in softmutation*/
	gdouble		mutationRate;		/*standard deviation in the strain matrix for lattice mutation*/
	gdouble		DisplaceInLatmutation;	/*max displacement in sofmutation of lattice mutation*/
	gboolean	AutoFrac;		/*automatic evolution of variation operators*/
/*4.5 Constrains*/
	gdouble		minVectorLength;	/*minimum length of a cell parameter*/
	gdouble 	*IonDistances;		/*triangular-matrix of the minimum allowed interatomic distances*/
	gint		constraint_enhancement;	/*allow c_e times stricter constraints of IonDistances*/
	gint 		_nmolecules;		/*HIDDEN: number of molecules*/
	gdouble		*MolCenters;		/*triangular-matrix of minimal distances between centers of molecules*/
/*4.6 Cell*/
	gint 		_nlatticevalues;	/*HIDDEN: number of lines in Latticevalues*/
	gdouble 	*Latticevalues;		/*initial volume/parameter of the unit cell*/
	gint 		_nsplits;		/*HIDDEN: number of splits*/
	gint 		*splitInto;		/*number of subcells into which unit cell is split*/
/*4.7 Restart*/
			/*Restart is unsupported*/
/*4.8 Ab initio calculations*/
	gint		_num_opt_steps;		/*HIDDEN: number of optimization steps*/
	gint 		*abinitioCode;		/*calculation code used for each optimization steps*/
	gdouble 	*KresolStart;		/*reciprocal-space resolution (2*pi*ang-1) for k-points generation*/
	gdouble		*vacuumSize;		/*amount of vacuum added at each optimization step*/
	gint 		numParallelCalcs;	/*number of structure relaxations at a time (in parallel)*/
	gchar 		*commandExecutable;	/*executable file/command for each optimization step*/
	gint		whichCluster;		/*type of job-submission*/
	gchar 		*remoteFolder;		/*remote folder where calculation is performed (whichCluster=2)*/
	gboolean	PhaseDiagram;		/*calculate a crude phase diagram w/ (rough) transition pressure*/
/*4.9 fingerprint*/
	gdouble		RmaxFing;		/*the distance cutoff (Ang)*/
	gdouble 	deltaFing;		/*discretization (Ang) of fingerprint function*/
	gdouble		sigmaFing;		/*gaussian broadening of interatomic distances*/
/*4.10 Antiseed*/
	gint 		antiSeedsActivation;	/*generation at which antiseeds will be active*/
	gdouble		antiSeedsMax;		/*gaussian height in msd of enthalpy in the generation among potential parents*/
	gdouble		antiSeedsSigma;		/*gaussian width in average distance in the generation among potential parents*/
/*4.11 Space group*/
	gboolean	doSpaceGroup;		/*determine space group*/
	gdouble		SymTolerance;		/*precision (ang) in the symmetry determination*/
/*4.12 Developers*/
	gint 		repeatForStatistics;	/*repeat USPEX calculation rFS times*/
	gdouble		stopFitness;		/*set halting fitness condition*/
	gboolean	collectForces;		/*collect all relaxation information (from VASP)*/
/*4.13 Seldom*/
	gboolean	ordering_active;	/*biasing of variation parameters*/
	gboolean	symmetrize;		/*all structures into std symmetrized crystallographic settings*/
	gint 		*valenceElectr;		/*number of valence electrons for each atom type*/
	gdouble		percSliceShift;		/*probability of shifting slabs*/
	gint		dynamicalBestHM;	/*number of survivors will vary: 0:no 1:yes 2:yes, with max diversity*/
	gchar 		*softMutOnly;		/*generations that are produced by softmutation*/
	gdouble		maxDistHeredity;	/*max cos distances between structure involved in heredity*/
	gint		manyParents;		/*more than 2 slices (or parents) are used for heredity*/
	gdouble		minSlice;		/*minimal thickness (Ang) of a slice*/
	gdouble		maxSlice;		/*maximal thickness (Ang) of a slice*/
	gint 		numberparents;		/*number of parents for cluster*/
/*5 Additional*/
/*5.1 molecular crystal*/
			/*no keyword in this section*/
/*5.2 Surfaces*/
	gdouble		thicknessS;		/*thickness (Ang) of the surface (ie adatoms) surface region*/
	gdouble		thicknessB;		/*thickness (Ang) of the buffer (no adatoms) surface region*/
	gint		reconstruct;		/*max number of multiplication of the surface cell (for reconstruction)*/
/*5.3 Variable composition*/
	gint		firstGeneMax;		/*how many different composition are sampled in the 1st generation*/
	gint		minAt;			/*minimum number of atoms/unit cell in the 1st generation*/
	gint 		maxAt;			/*maximum number of atoms/unit cell in the 1st generation*/
	gdouble		fracTrans;		/*fraction of structures generated by transmutation*/
	gdouble		howManyTrans;		/*maximum fraction of atoms in the structure that are being transmuted*/
	gint 		_nspetrans;		/*HIDDEN: number of specific transmutations*/
	gint 		*specificTrans;		/*specific transmutations*/
/*5.4 metadynamics*/
	gdouble		GaussianWidth;		/*width of each gaussian added to surface energy*/
	gdouble		GaussianHeight;		/*height of each gaussian added to surface energy*/
	gint 		FullRelax;		/*whether relax within a fixed cell, fully relax best, or all structures*/
	gdouble		maxVectorLength;	/*maximum (Ang) unit cell length (if larger then steep force correction is added)*/
/*5.5 Particle swarm optimization*/
	gdouble		PSO_softMut;		/*weight of softmutation*/
	gdouble		PSO_BestStruc;		/*weight of heredity with best position, within same PSO particle*/
	gdouble		PSO_BestEver;		/*weight of heredity with globally best PSO particle*/
/*6 Phase Transition Pathways*/
/*6.1 VC-NEB*/
		/*no keyword in this section*/
/*6.2 VC-NEB options*/
	gint		vcnebType;		/*type of VC-NEB*/
	gint 		_vcnebtype_method;	/*HIDDEN: X.. digit in vcnebType: use VCNEB or just relax*/
	gboolean	_vcnebtype_img_num;	/*HIDDEN: .X. digit in vcnebType: fixed or variable image number*/
	gboolean	_vcnebtype_spring;	/*HIDDEN: ..X digit in vcnebType: fixed or variable spring constant*/
	gint		numImages;		/*initial number of images*/
	gint 		numSteps;		/*maximum VC-NEB steps*/
	gint 		optReadImages;		/*reading the Image file (all image, initial and final, or +intermediate)*/
	gboolean	optimizerType;		/*optimization algorithm*/
	gint		optRelaxType;		/*relaxation mode (fixed cell, only cell, full relaxation)*/
	gdouble		dt;			/*time step*/
	gdouble		ConvThreshold;		/*halting RMS (eV/Ang) condition*/
	gdouble		VarPathLength;		/*min distance (* 1.5) between 2 image to generate a new image*/
	gdouble		K_min;			/*min spring constant, for variable spring cte calculation*/
	gdouble		K_max;			/*max spring constant, for variable spring cte calculation*/
	gdouble		Kconstant;		/*spring constant, for variable spring cte calculation*/
	gboolean	optFreezing;		/*freeze image structure (when ConvThreshold is achieve)*/
	gint		optMethodCIDI;		/*use climbing image / downing image method*/
	gint		startCIDIStep;		/*CI/DI method starting step*/
	gint 		_npickimg;		/*HIDDEN: number of picked-up images*/
	gint 		*pickupImages;		/*# of images chosen for CI/DI*/
	gint		FormatType;		/*structure format in the output pathway (ONLY VASP5 is supported!)*/
	gint		PrintStep;		/*save VC-NEB restart files in STEP directory every PS steps*/
/*6.3 Set initial pathway in VC-NEB*/
		/*no keyword in this section*/
} uspex_calc_struct;
/* global output structure */
typedef struct {
/*name*/
        gchar *name;
        gint version;
	uspex_calc_struct calc;
/*system parameters*/
        uspex_method method;
        gint type;
        /* in detail */
        gint dim;
        gboolean mol;
        gboolean var;
        /* optimization */
        gint opt_type;
        gboolean have_fitness;
/*structures*/
        gint num_gen;
        gint num_struct;
        gint nspecies;          /*number of species*/
        gint *spe_Z;            /*atomic number of each species*/
        gint *red_index;        /*in case not all structure are displayed (ie. META)*/
        gdouble min_E;
        gdouble max_E;
        uspex_individual *ind;
        gint num_best;
        gint *best_ind;
/*interpretation*/
        gpointer graph;
        gpointer graph_best;
        gpointer graph_comp;
} uspex_output_struct;
/* execution structure */
typedef struct {
	int job_id;
	gboolean have_result;
	/*job related*/
	gchar *job_uspex_exe;
	gchar *job_path;
} uspex_exec_struct;


/*methods of interest*/


