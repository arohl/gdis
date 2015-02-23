
/**************/
/* prototypes */
/**************/

gint core_match(const gchar *, struct core_pak *);
gint shel_match(const gchar *, struct shel_pak *);
gint pair_match(const gchar *, const gchar *,
                struct core_pak *, struct core_pak *);

void info_bond(GtkWidget *, gint, gint, struct model_pak *);
void info_dist(gint, gint, struct model_pak *);
void info_angle(gint, gint, struct model_pak *);
void info_torsion(gint, gint, struct model_pak *);

void elem_change_colour(GtkWidget *, gpointer *);

void gulp_files_init(struct model_pak *);
void gulp_data_free(struct model_pak *);
void gulp_files_free(struct model_pak *);
void gulp_data_copy(struct model_pak *, struct model_pak *);
void gulp_extra_copy(struct model_pak *, struct model_pak *);

gint gulp_cosmo_points(gint, gint, gint);

gint free_moldy_data(struct model_pak *);

gint search_basename(const gchar *);
gint dialog_active(gint);

void unhide_atoms(void);

void type_model(const gchar *, struct model_pak *);

void periodicity_page(GtkWidget *);

void atom_numbers_box(GtkWidget *);

void delete_commit(struct model_pak *);

void core_delete_single(struct core_pak *, struct model_pak *);
void core_delete_all(struct model_pak *);


