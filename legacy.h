#ifndef _LEGACY_H
#define _LEGACY_H 1

/* Frame data pak... */
struct rdf_frame_pak
{
  gint atom_count[2];
  struct cartesian *firstatom;
  struct cartesian *secondatom;
};

/* RDF data pak... */
struct rdf_pak
{
  gint frames_count;
  gint frame_atoms[2];
  gint atom_element[2];
  gint atom_total[2];
  struct rdf_frame_pak *frame;
};

/* Interface functions... */
void create_rdfwindow(void);
void on_atom1select_toggled (GtkToggleButton *togglebutton, gpointer user_data);
void on_atom2select_toggled (GtkToggleButton *togglebutton, gpointer user_data);
void on_atom3select_toggled (GtkToggleButton *togglebutton, gpointer user_data);

/* Queueing system interface function... */
void init_rdf_task(void);

/* Main process loop */
void exec_rdf_task(gpointer *ptr);

/* Memory allocation / deallocation / scrubbing functions... */
struct rdf_pak *rdf_init(gint num_frames);
gint rdf_cleanup(struct rdf_pak *rdf);
gint rdf_add_frame(struct rdf_pak *rdf, gint frames_total);
gint rdf_add_coords(struct rdf_frame_pak *frame, gint atom);

/* Data sorting / retrieval functions... */
gint rdf_frame_atoms(struct model_pak *model, gint atom_select, gint atom_element);
gint get_rdf_data(struct rdf_pak *rdf, struct model_pak *model, gint atom_select);

/* Calculation functions... */
gint distances_atoms(struct rdf_pak *rdf, struct model_pak *model, gint *distances, gint hist_points, gdouble scale, gint count_method);
gint distances_elements(struct rdf_pak *rdf, struct model_pak *model, gint *distances, gint hist_points, gdouble scale, gint count_method);

#endif
