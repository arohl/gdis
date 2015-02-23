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

along with this program; if not, write to the Free Software
Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.

The GNU GPL can also be found at http://www.gnu.org
*/

/*******************************/
/* model sub-structure (atoms) */
/*******************************/
struct core_pak
{
/* identifiers */
gint atom_code;
guint atom_order;
gchar *atom_label; /* NEW - replacement for label[] */
gchar *atom_type;  /* NEW - FF type */

/*VZ*/
gdouble atom_nmr_shift;//isotropic shift
gdouble atom_nmr_aniso;//CSA anisotropy
gdouble atom_nmr_asym;//CSA assymetry
gdouble atom_nmr_cq;//EFG Cq
gdouble atom_nmr_efgasym;//EFG assymetry

/* TODO - data structure that contains res_name and res_no and each atom points to this */
gchar *res_name;   /* NEW - residue name */
gint res_no;       /* NEW - residue number */
gchar *chain;       /* NEW - chaincw name */

/* flags */
gint status;
gint primary;          /* one of the initial atoms? */
gint orig;             /* original atom? (ie part of whole unit cell) */
gint ghost;            /* QM ghost atom */
gint breathe;
gint growth;           /* growth slice */
gint translate;        /* translation marker */
gint render_mode;      /* render type */
gint render_wire;      /* render as wire-frame */

/* connectivity data */
gint molecule;
GSList *bonds;
struct shel_pak *shell;
struct mol_pak *mol;

/* symmetry related core (if non primary) */
struct core_pak *primary_core;

/* region type */
gint region;
/* coordinate data, fractional or cartesian (inhomogeneous) */
gdouble x[4];
/* rotated cartesian (inhomogeneous) */
gdouble rx[4];
/* velocities */
gdouble v[3];

/* vibration vector lists */
GSList *vibx_list;
GSList *viby_list;
GSList *vibz_list;

/* NEW: is a hydrogen capable of forming a h-bond
   instead of connecting all hydrogens within a cut-off
   distance to acceptor atoms, only connect them if this 
   flag is true */
gboolean hydrogen_bond;

/* bonding cutoff */
gdouble bond_cutoff;
/* site occupancy factor */
gint has_sof;
gdouble sof;
/* breathing radius */
gdouble radius;
/* charge info */
gint lookup_charge;
gdouble charge;
/* fitting flags */
gchar *flags;

/* display data */
gdouble colour[4];
gdouble offset[3];
};

/********************************/
/* model sub-structure (shells) */
/********************************/
struct shel_pak
{
/* identifiers */
gchar *shell_label;
/*
gchar element[ELEM_LABEL_SIZE];
*/
/* flags */
gint atom_code;    /* atom type of the shell */
gint status;       /* normal, deleted, highlighted etc. */
gint primary;      /* one of the initial atoms? */
gint orig;         /* original atom? (ie part of whole unit cell) */
gint breathe;
gint translate;    /* translation marker */
/* associated core (if any) */
struct core_pak *core;
/* symmetry related shell (if non primary) */
struct shel_pak *primary_shell;
/* breathing radius */
gdouble radius;
/* region type */
gint region;
/* coord data */
gdouble x[4];
gdouble rx[4];
/* velocity */
gdouble v[3];
/* periodic image coordinate */
gint pic[3];
/* charge info */
gint lookup_charge;
gdouble charge;
/* site occupancy factor */
gint has_sof;
gdouble sof;
/* display data */
gdouble colour[3];
/* fitting flags */
gchar *flags;
/* offset (shell only, NOT global) */
gdouble offset[3];
};

/*******************************/
/* model sub-structure (bonds) */
/*******************************/
struct bond_pak
{
/* deleted/normal/hidden */
gint status;
/* single/double etc. */
gint type;
/* relative (to atom1) fractional offset for bond midpoint */
gdouble offset[3];

/* constituent atom indices */
struct core_pak *atom1;
struct core_pak *atom2;
};

/***********************************/
/* model sub-structure (molecules) */
/***********************************/
struct mol_pak
{
GSList *cores;
gdouble centroid[3];
};

/****************************/
/* periodic image structure */
/****************************/
struct image_pak
{
/* image coordinates */
gint pic[3];
/* cartesian coordinates */
gdouble rx[3];
};

/**********************************/
/* coords/connectivity prototypes */
/**********************************/

/* debugging */
void print_core(struct core_pak *);
void print_core_cart(struct core_pak *);
void print_shell(struct shel_pak *);
void print_cores(struct model_pak *);
void print_cores_cart(struct model_pak *);
void print_shells(struct model_pak *);
void print_core_shell(struct model_pak *);
void print_core_list(GSList *);
void print_core_list_cart(GSList *);
void print_shell_list(GSList *);

/* main */
void coords_init(gint, struct model_pak *);
void coords_init_units(struct model_pak *);
void coords_compute(struct model_pak *);
gint coords_center(struct model_pak *);
void coords_confine_cores(GSList *, struct model_pak *);
void coords_confine_centroid(struct mol_pak *, struct model_pak *);
void coords_make_fractional(struct model_pak *);
void coords_make_cartesian(struct model_pak *);

void fractional_clamp(gdouble *, gint *, gint);
void fractional_min(gdouble *, gint);

void core_free(gpointer);
void free_core_list(struct model_pak *);
void free_mol_list(struct model_pak *);

void core_init(gchar *, struct core_pak *, struct model_pak *);
gpointer core_new(gchar *, gchar *, struct model_pak *);
gpointer shell_new(gchar *, gchar *, struct model_pak *);

/* deprec */
struct core_pak *new_core(gchar *, struct model_pak *);
struct shel_pak *new_shell(gchar *, struct model_pak *);

struct core_pak *dup_core(struct core_pak *);
struct shel_pak *dup_shel(struct shel_pak *);
struct bond_pak *dup_bond(struct bond_pak *);

GSList *dup_core_list(GSList *);
GSList *dup_shell_list(GSList *);

struct core_pak *copy_core(struct core_pak *, struct model_pak *, struct model_pak *);

void delete_commit(struct model_pak *);
void delete_core(struct core_pak *);
void delete_duplicate_cores(struct model_pak *);
void add_atom(gint, gint, struct model_pak *);

void elem_init(struct core_pak *, struct model_pak *);

void atom_colour_scheme(gint, struct core_pak *, struct model_pak *);
void model_colour_scheme(gint, struct model_pak *);

void init_atom_colour(struct core_pak *, struct model_pak *);
void init_atom_charge(struct core_pak *, struct model_pak *);
void init_model_charges(struct model_pak *);
gdouble atom_charge(struct core_pak *);
void calc_emp(struct model_pak *);

GSList *find_unique(gint, struct model_pak *);

void cor_calc_xlimits(gdouble *, gdouble *, GSList *);

void shell_make_links(struct model_pak *);

void atom_numbers_update(struct model_pak *);

void connect_molecules(struct model_pak *);
void connect_bonds(struct model_pak *);
void connect_atom_clear(struct core_pak *, struct model_pak *);
void connect_atom_compute(struct core_pak *, struct model_pak *);
void connect_atom_refresh(struct core_pak *, struct model_pak *);
void connect_refresh(struct model_pak *);
void connect_refresh_global(void);
void connect_centroid_compute(struct mol_pak *);
void wipe_bonds(struct model_pak *);
gint connect_split(struct bond_pak *);

void connect_image_vector(gdouble *, struct core_pak *, struct core_pak *, gint);

void connect_user_bond(struct core_pak *, struct core_pak *, gint, struct model_pak *);
void connect_make_bond(struct core_pak *, gint, struct model_pak *);
void connect_merge_user(struct model_pak *);

void connect_fragment_init(struct model_pak *);
GSList *connect_fragment_get(struct core_pak *,
                             struct core_pak *,
                             struct model_pak *);
GSList *connect_neighbours(struct core_pak *) ;

