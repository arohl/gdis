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

struct keyword_pak gulp_keyword[] = 
{
{"sing", E_SINGLE},
{"opti", E_OPTIMIZE},
{"md",   MD},
{"bfgs", BFGS_OPT},
{"rfo",  RFO_OPT},
{"conj", CONJ_OPT},
{"conp", CONP},
{"conv", CONV},
/* technically, cellonly */
{"cell", CELL},
{"noflag", NOFLAGS},
{"mole", MOLE},
{"molm", MOLMEC},
{"molq", MOLQ},
{"qeq", QEQ},
{"free", FREE_ENERGY},
{"zsisa", ZSISA},
{"comp", COMPARE},
{"nosym", NOSYM},
{"fix", FIX},
{"noauto", GULP_NOAUTOBOND},
{"cosmi", GULP_SOLVATION_COSMIC},
{"cosmo", GULP_SOLVATION_COSMO},
{"phon", PHONON},
{"eigen", EIGEN},
{"unit", UNIT_HESSIAN},
{NULL, -1}
};

struct keyword_pak gulp_option[] = 
{
{"kpoin", KPOINTS},
{"title", TITLE},
{"endf", ENDFORCE},
{"end",  END},
{"name", NAME},
{"pcel", POLYMER_CELL},
{"pvec", POLYMER_VECTOR},
{"scel", SURFACE_CELL},
{"svec", SURFACE_VECTORS},
{"vect", LATTICE_VECTORS},
{"cell", CELL},
{"pfrac", PFRAC},
{"sfrac", SFRAC},
{"frac", FRAC},
{"cart", CART},
{"spac", SPACE},
{"orig", ORIGIN},
{"spec", SPECIES},
{"maxc", MAXCYC},
{"tota", TOTAL_ENERGY},
{"dhkl", D_HKL},
{"sbulk", SBULK_ENERGY},
{"dump", DUMP_FILE},
{"obse", GULP_OBSERVABLES},
{"ensem", ENSEMBLE},
{"pres", PRESSURE},
{"temp", TEMPERATURE},
{"timestep", TIMESTEP},
{"equi", EQUILIBRATION},
{"prod", PRODUCTION},
{"samp", SAMPLE},
{"conn", GULP_CONNECT},

/* NEW */
{"cosmoshape", GULP_COSMO_SHAPE},
{"points", GULP_COSMO_POINTS},
{"segments", GULP_COSMO_SEGMENTS},
{"solventeps", GULP_COSMO_SOLVENT_EPSILON},
{"solventrad", GULP_COSMO_SOLVENT_RADIUS},
{"solventrmax", GULP_COSMO_SOLVENT_RMAX},
{"rangeforsmooth", GULP_COSMO_SMOOTHING},

{"writ", WRITE},
{"prin", PRINT},
{"super", SUPER_CELL},
{"outp", OUTPUT},
{"libr", GULP_LIBRARY},

/* special case potentials group */
{"qwolf", GULP_SINGLE_LINE_POTENTIAL},
{"bren", GULP_SINGLE_LINE_POTENTIAL},

/* NEW - don't explicitly name FF type */
{"atom", GULP_POTENTIAL},
{"axil", GULP_POTENTIAL},
{"bacross", GULP_POTENTIAL},
{"bcos", GULP_POTENTIAL},
{"bcross", GULP_POTENTIAL},
{"botwo", GULP_POTENTIAL},
{"borep", GULP_POTENTIAL},
{"boatt", GULP_POTENTIAL},
{"bsm", GULP_POTENTIAL},
{"buck", GULP_POTENTIAL},
{"coul", GULP_POTENTIAL},
{"covexp", GULP_POTENTIAL},
{"damp", GULP_POTENTIAL},
{"eam_", GULP_POTENTIAL},
{"epsi", GULP_POTENTIAL},
{"expon", GULP_POTENTIAL},
{"fermi", GULP_POTENTIAL},
{"four", GULP_POTENTIAL},
{"gene", GULP_POTENTIAL},
{"harm", GULP_POTENTIAL},
{"igauss", GULP_POTENTIAL},
{"lenn", GULP_POTENTIAL},
{"lin3", GULP_POTENTIAL},
{"ljbuff", GULP_POTENTIAL},
{"many", GULP_POTENTIAL},
{"mors", GULP_POTENTIAL},
{"murrell", GULP_POTENTIAL},
{"outofp", GULP_POTENTIAL},
{"poly", GULP_POTENTIAL},
{"qtaper", GULP_POTENTIAL},
{"qerfc", GULP_POTENTIAL},
{"ryck", GULP_POTENTIAL},
{"rydb", GULP_POTENTIAL},
{"spline", GULP_POTENTIAL},
{"spri", GULP_POTENTIAL},
{"squared", GULP_POTENTIAL},
{"stillinger", GULP_POTENTIAL},
{"sw2", GULP_POTENTIAL},
{"sw3", GULP_POTENTIAL},
{"thre", GULP_POTENTIAL},
{"tors", GULP_POTENTIAL},
{"tortaper", GULP_POTENTIAL},
{"torexp", GULP_POTENTIAL},
{"torharm", GULP_POTENTIAL},
{"tsune", GULP_POTENTIAL},
{"urey", GULP_POTENTIAL},

{"switch", SWITCH_ON},
{"bfgs", BFGS_OPT},
{"rfo",  RFO_OPT},
{"conj", CONJ_OPT},
{"gnorm", GNORM},
{"cycle", CYCLE},
{"keyword", KEYWORD},
{"elem", ELEMENT},
{"symbol", SYMBOL},
{"number", NUMBER},
{"weight", WEIGHT},
{"variable", GULP_VARIABLES},
{"cova", COVALENT},
{"ionic", IONIC},
{"vdw", VDW},
{"ignore", IGNORE},
{"erongi", ERONGI},
{NULL, -1}
};

