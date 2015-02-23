
/* surface generation */
#define EPSILON 1e-6
#define FILL_EPS 0.00025

/* surface convergence */
#define MIN_THICKNESS 8.0
#define MAX_CONV_CYCLES 20
/* result (almost always) good to 4dp */
#define MAX_DESURF 0.0001
#define MAX_DEEATT 0.001

/* mode control masks */
enum
{
SINGLE_SHIFT, SINGLE_PLANE, SINGLE_SURFACE,
ALL_SHIFTS, ALL_PLANES, INVALID_SHIFTS
};

/* 2D vdw surface parameter specifiers */
enum 
{
MS_ACCESSIBLE, MS_MOLECULAR, MS_EDEN,
MS_DSIZE, MS_DEPTH, MS_OFFSET, MS_PRAD, MS_ACCURACY,
MS_TOUCH, MS_AFM, MS_EPOT, MS_HIRSHFELD, MS_SOLVENT,
MS_DE, MS_CURVEDNESS, MS_SHAPE_INDEX, MS_SSATOMS
};

/* surface specifiers */
enum 
{
CALC_SHIFTS, CALC_ENERGY, CONV_REGIONS, RANK_FACES, MAKE_FACES,
ADD_SHIFT, DELETE_SHIFT, GENERATE 
};

/* layer structure */

struct layer_pak
{
gdouble width;
gdouble centroid[3];
GSList *cores;
};

/* prototypes */
gint generate_surface(struct model_pak *, struct model_pak *);

gpointer plane_new(gdouble *, struct model_pak *);
gpointer shift_new(gdouble);
void plane_free(gpointer);
void shift_free(gpointer);
void plane_data_free(GSList *);
void shift_data_free(GSList *);
gpointer plane_find(gdouble *, struct model_pak *);

void update_plane_energy(struct plane_pak *, struct model_pak *data);

GSList *get_ranked_faces(gint, gdouble, struct model_pak *);
gint rank_faces(void);
gint calc_shifts(void);
void calc_emp(struct model_pak *);
gint region_max(struct model_pak *);
gint surf_sysabs(struct model_pak *, gint, gint, gint);
void surf_symmetry_generate(struct model_pak *);

GSList *get_facet_equiv(struct model_pak *, gint *);

void update_surface_dialog(struct model_pak *);

void select_shift(GtkWidget *, gint, gint);

void sort_coords(struct model_pak *);

int GCD(int, int);

void dock_selection(gchar *, struct model_pak *);

void diffract_layer_setup(struct model_pak *);

gint facet_equiv(struct model_pak *, gint *, gint *);

gint dhkl_compare(gpointer, gpointer);

void free_vertices(struct model_pak *);

void morph_sculpt(GtkWidget *, gpointer);

gpointer plane_dup(struct plane_pak *);

void surf_shift_explore(struct model_pak *, struct surface_pak *);

gint region_move_atom(struct core_pak *, gint, struct model_pak *);

gpointer make_surface(struct model_pak *,
                      struct plane_pak *,
                      struct shift_pak *);

