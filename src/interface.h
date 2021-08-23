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

/* spacing */
#define PANEL_SPACING 4
/* sub windows */
#define W_BORDER 5
#define H_BORDER 5
/* graphics window defs */
#define START_WIDTH 640
#define START_HEIGHT 600
#define MIN_WIDTH 200
#define MIN_HEIGHT 200
/* use GTKGL_ as GL_ are often reserved by OpenGL */
#define GTKGL_DIM 600
#define GTKGL_LINE_WIDTH 1.5
/* angle drawing */
#define MIN_ARC_RAD 15

#define PIX2SCALE 0.005

/* main operation modes */
enum 
{
FREE, XLAT, ROLL, YAW, PITCH, LOCKED, QUIT,
SELECT_FRAGMENT,
CREATE_RIBBON,
ATOM_ADD, ATOM_CHANGE_TYPE,
BOND_NORMAL, BOND_SPLIT,
BOND_DELETE, BOND_SINGLE, BOND_DOUBLE, BOND_TRIPLE, BOND_HBOND, BOND_ZEOLITE,
BOND_INFO, DIST_INFO, ANGLE_INFO, DIHEDRAL_INFO,
DEFINE_RIBBON, DEFINE_VECTOR, DEFINE_PLANE,
ANIMATING, RECORD
};

/* widget & file/package types */
/* including dialog request types */
enum 
{
GDIS, SGINFO, SYMMETRY, GPERIODIC, ELEM_EDIT, POVRAY, GENSURF,
ABOUT, ANIM, PLOTS, SURF, DISPLAY, GEOMETRY, SPATIAL, TASKMAN, MANUAL, SETUP,
FILE_SELECT, FILE_LOAD, FILE_SAVE, FILE_SAVE_AS,
NODATA, DATA, BIOSYM, CIF, FDF, GULP, MONTY, MARVIN, MORPH,
META_DATA, XML, XTL, XYZ, MOL2,
AUTO, CSSR, GAMOUT, PDB, GEOMVIEW_OFF, PROJECT, DOCKING,
MDI, CREATOR, MVNOUT, GULPOUT, GULP_TRJ, SIESTA_OUT, VASP, CVASP,
USPEX, CUSPEX, REPLACE_ATOMS, DYNAMICS, ZMATRIX,
OPENGL, OPENGL_OPTIONS, GAMESS, GAMESS_OUT, DIFFAX, DIFFAX_INP, DMOL_INPUT,
ABINIT, ABINIT_OUT, NWCHEM, NWCHEM_OUT, CASTEP, CASTEP_OUT, GAUSS, GAUSS_OUT,
QE, QE_OUT,
HIRSHFELD, CONNECT, RDF, LIQUIDS, VACF, MD_ANALYSIS,
PICTURE, TEXT, RIETICA, CEL, DLPOLY, CRYSTAL_GRAPH, BIOGRAF, DLP, GROMACS,
LAST
};

/* model display features */
enum 
{
PLANE_LABELS, ASYM_TOGGLE,
ATOM_LABELS, FRAME, MESH, SHELLS, AXES_TYPE, 
PBC_CONFINE_NONE, PBC_CONFINE_ATOMS, PBC_CONFINE_MOLS 
};

/* switch_view() call modes */
enum {CANVAS_INIT, CANVAS_SINGLE, CANVAS_VSPLIT, CANVAS_HSPLIT, CANVAS_HVSPLIT,
      PREV_MODEL, NEXT_MODEL};

/* view actions */
enum {ROTATION, UPDATE_X, UPDATE_Y, UPDATE_Z};

enum {GUI_CANVAS, GUI_MODEL_TREE, GUI_MODEL_PROPERTIES, GUI_TEXT_BUFFER};

/* prototypes */

/* main */
void gui_init(int, char **);
gint gui_motion_event(GtkWidget *, GdkEventMotion *);
gint gui_press_event(GtkWidget *, GdkEventButton *);
gint gui_scroll_event(GtkWidget *, GdkEventScroll *);
gint gui_release_event(GtkWidget *, GdkEventButton *);

/* display control */
void gui_mode_switch(gint);
void gtk_mode_switch(GtkWidget *, gint);
void gui_text_show(gint, gchar *);
void gui_model_select(struct model_pak *);
void gui_angles_refresh(void);
void gui_angles_reset(void);
void gui_camera_refresh(void);

/* model tree */
void tree_select_delete(void);
void tree_select_active(void);
void tree_select_model(struct model_pak *);
void tree_model_add(struct model_pak *);
void tree_model_refresh(struct model_pak *);
void tree_init(GtkWidget *);

/* NEW: export edit_model_create for use in button toolbar */
void edit_model_create(void);

void diffract_select_peak(gint, gint, struct model_pak *);
void analysis_export_dialog(void);
void symmetry_widget_redraw(void);
void image_export_dialog(void);
void image_import_dialog(void);
void image_write(gpointer);

void gui_grid_dialog(void);
void gui_job_dialog(void);
void gui_animate_dialog(void);
void gui_diffract_dialog(void);
void gui_gamess_widget(GtkWidget *, struct model_pak *);
void gui_gulp_task(GtkWidget *, struct model_pak *);
void gui_gamess_task(GtkWidget *, struct model_pak *);
void gui_measure_dialog(void);
void gui_edit_dialog(void);
void gui_analysis_dialog(void);
void gui_isosurf_dialog(void);
void gui_gperiodic_dialog(void);
void gui_plots_dialog(void);
void gui_render_dialog(void);
void gui_mdi_dialog(void);
void gui_help_dialog(void);
void gui_about_dialog(void);
void gui_dock_dialog(void);
void gui_defect_dialog(void);
void gui_setup_dialog(void);
void gui_siesta_dialog(void);
void gui_vasp_dialog(void);
void gui_uspex_dialog(void);

void gui_surface_widget(GtkWidget *);
void gui_surface_setup(GtkWidget *);
void gui_surface_refresh(GtkWidget *);

void gui_edit_widget(GtkWidget *);
void gui_display_widget(GtkWidget *);
void gui_symmetry_refresh(GtkWidget *);

void gui_active_refresh(void);

void gui_view_x(void);
void gui_view_y(void);
void gui_view_z(void);
void gui_view_a(void);
void gui_view_b(void);
void gui_view_c(void);

void canvas_new(gint, gint, gint, gint);
void canvas_init(GtkWidget *);
void canvas_shuffle(void);
void canvas_single(void);
void canvas_create(void);
void canvas_delete(void);
void canvas_resize(void);
void canvas_select(gint, gint);
gpointer canvas_find(struct model_pak *);

void help_init(void);

void redraw_canvas(gint);

gint gui_canvas_handler(gpointer *);
gint gui_widget_handler(void);

void gui_refresh(gint);

void gui_refresh_selection(void);

gpointer camera_new(void);
gpointer camera_dup(gpointer);
void camera_copy(gpointer, gpointer);
void camera_init(struct model_pak *);
void camera_reset(struct model_pak *);
void camera_dump(gpointer);
void camera_view(gdouble *, gpointer);
void camera_waypoint_animate(gint, gint, struct model_pak *);
void camera_rotate_animate(gint, gdouble *, gint, struct model_pak *);
gdouble *camera_q(gpointer);
void camera_rescale(gdouble, gpointer);

void gui_track_output(void);
