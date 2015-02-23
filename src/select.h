
void set_selection_mode(GtkWidget *, gpointer *);

void select_clear(struct model_pak *);
void select_copy(void);
void select_paste(void);

void select_colour(void);
void select_delete(void);
void select_hide(void);
void select_invert(void);
void select_all(void);

void select_flag_ghost(void);
void select_flag_normal(void);

void select_core(struct core_pak *, gint, struct model_pak *);

gint select_add_core(struct core_pak *, struct model_pak *);
void select_add_mol(struct mol_pak *, struct model_pak *);
void select_add_fragment(struct core_pak *, struct model_pak *);
void select_add_region(struct core_pak *, struct model_pak *);
void select_del_core(struct core_pak *, struct model_pak *);
void select_all_labels(struct core_pak *, struct model_pak *);
void select_all_types(struct core_pak *, struct model_pak *);
void select_all_elem(gint, struct model_pak *);
void select_all_elem_in_mol(struct core_pak *, struct model_pak *);
void select_model_labels(void);

void core_render_mode_set(gint, GSList *);
void core_render_wire_set(gint, GSList *);

void rotate_select(struct model_pak *, gdouble *);
void select_translate(gint, gint, struct model_pak *);

void construct_backbone(struct core_pak *, struct model_pak *);
void construct_ribbon(GList *, struct model_pak *);

GSList *exit_branches(struct core_pak *, struct core_pak *);

/* NEW - structures for exploring the connectivity/ribbon */
struct link_pak
{
gint id;
gint type;
gint size;

/* NB: first & last atoms of the link are nodes */
struct core_pak *first;
struct core_pak *last;
GList *chain;
};

struct node_pak
{
struct core_pak *core;
gdouble x[3];

/* FIXME - need better way of indicating explored branches */
gint num_branches;
gint num_explored;

/* TODO - a list dedicated to the 1st atom in the explored pathways */
/* ie makes it easier to check for unexplored exit pathways */
GSList *exit_list;

/* atom chains (and status) */
/* each element in list is a pointer to a link_pak */
GList *link_list;
};


/****************************/
/* ribbon segment structure */
/****************************/
struct ribbon_pak
{
gdouble colour[4];
/* cyclic group id's */
gint id1;
gint id2;
/* two connected cyclic groups (orig position) */
gdouble x1[3];
gdouble x2[3];
/* two connected cyclic groups (transformed position) */
gdouble r1[3];
gdouble r2[3];
/* normal at each group (orig position) */
gdouble u1[3];
gdouble u2[3];
/* normal at each group (transformed position) */
gdouble n1[3];
gdouble n2[3];
/* ribbon orientation vectors at each group (orig position) */
gdouble v1[3];
gdouble v2[3];
/* ribbon orientation vectors at each group (transformed position) */
gdouble o1[3];
gdouble o2[3];
};

