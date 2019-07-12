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

#include "gui_defs.h"

typedef enum {
        GRAPH_REGULAR=0,	/*previous plain graph*/
	/*  NEW - graph: generic data types*/
	GRAPH_IY_TYPE,		/*nx integer {x} + ny sets of {nx real y for all x}*/
	GRAPH_XY_TYPE,		/*nx real    {x} + ny sets of {nx real y for all x}*/
	GRAPH_YX_TYPE,		/*interleave {y,x} set (for BANDOS)*/
	GRAPH_IX_TYPE,		/*nx integer {x} + nx sets of {ny real y for one x}*/
	GRAPH_XX_TYPE,		/*nx real    {x} + nx sets of {ny real y for one x}*/
} graph_type;
typedef enum {
	/* NEW - graph: symbol shapes */
	GRAPH_SYMB_NONE=0,
	GRAPH_SYMB_CROSS,
	GRAPH_SYMB_SQUARE,
	GRAPH_SYMB_TRI_DN,
	GRAPH_SYMB_TRI_UP,
	GRAPH_SYMB_DIAM,
} graph_symbol;
typedef enum {
	/* NEW - graph: line shapes */
	GRAPH_LINE_NONE=0,	/*no lines*/
	GRAPH_LINE_SINGLE,	/*single line*/
	GRAPH_LINE_DASH,	/*dashed line*/
	GRAPH_LINE_DOT,		/*dotted line*/
	GRAPH_LINE_THICK,	/*thick line*/
} graph_line;
typedef enum {
	/* NEW - graph: colors */
	GRAPH_COLOR_DEFAULT=-1, /*whatever title color*/
	GRAPH_COLOR_BLACK=0,	/*0.0,	0.0,	0.0*/
	GRAPH_COLOR_WHITE,	/*1.0,	1.0,	1.0*/
	GRAPH_COLOR_BLUE,	/*0.0,	0.0,	1.0*/
	GRAPH_COLOR_GREEN,	/*0.0,	0.5,	0.0*/
	GRAPH_COLOR_RED,	/*1.0,	0.0,	0.0*/
	GRAPH_COLOR_YELLOW,	/*1.0,	1.0,	0.0*/
	GRAPH_COLOR_GRAY,	/*0.5,	0.5,	0.5*/
	GRAPH_COLOR_NAVY,	/*0.0,	0.0,	0.5*/
	GRAPH_COLOR_LIME,	/*0.0,	1.0,	0.0*/
	GRAPH_COLOR_TEAL,	/*0.0,	0.5,	0.5*/
	GRAPH_COLOR_AQUA,	/*0.0,	1.0,	1.0*/
	GRAPH_COLOR_MAROON,	/*0.5,	0.0,	0.0*/
	GRAPH_COLOR_PURPLE,	/*0.5,	0.0,	0.5*/
	GRAPH_COLOR_OLIVE,	/*0.5,	0.5,	0.0*/
	GRAPH_COLOR_SILVER,	/*0.75,	0.75,	0.75*/
	GRAPH_COLOR_FUSHIA,	/*1.0,	0.0,	1.0*/
} graph_color;
typedef struct {
	gint x_size;            /*size of the x array*/
	gdouble *x;             /* -> x array*/
} g_data_x;
typedef struct {
	gint y_size;            /*size of a y array*/
	gdouble *y;             /* -> y array*/
	gint32 *idx;            /*index up to 2,147,483,647 (<0 = missing structure)*/
	graph_type type;	/* NEW - type can change for each set! */
	graph_symbol *symbol;	/*symbol for each data*/
	graph_color *sym_color; /*a color for each symbol*/
	gboolean mixed_symbol;  /*preserve symbols from being changed all at once*/
	graph_line line;	/*one line type per set*/
	graph_color color;	/*TODO graph set color*/
} g_data_y;
/************************************/
/* UI definition for graph_controls */
/************************************/
struct graph_ui{
	struct model_pak *ref_model;
	gpointer graph_active;
	GUI_OBJ *title;
	GUI_OBJ *sub_title;
	GUI_OBJ *x_title;
	GUI_OBJ *y_title;
	GUI_OBJ *xmin;
	GUI_OBJ *xmax;
	GUI_OBJ *ymin;
	GUI_OBJ *ymax;
	gdouble xticks;/*GTK_SPIN*/
	GUI_OBJ *xtics;
	gdouble yticks;/*GTK_SPIN*/
	GUI_OBJ *ytics;
	gboolean auto_x;
	gboolean auto_y;
	GUI_OBJ *set;
	gdouble set_number;/*GTK SPIN*/
	GUI_OBJ *size;
	GUI_OBJ *by_val;
	gboolean by_value;
	GUI_OBJ *num;
	gdouble num_number;/*GTK SPIN*/
	GUI_OBJ *type;
	GUI_OBJ *symbol;
	GUI_OBJ *idx;
	GUI_OBJ *x_val;
	GUI_OBJ *y_val;
	GUI_OBJ *line;
	GUI_OBJ *color;
};
/*******************/
/* data structures */
/*******************/
struct graph_pak
{
	gint grafted;
	gchar *treename;
	gint treenumber;
	/* graph generation parameters */
	gdouble wavelength;
	/* flags */
	gint xlabel;
	gint ylabel;
	graph_type type;	/*TODO eliminate this*/
	gboolean require_xaxis;
	gboolean require_yaxis;
	/* graph layout */
	gint xticks;
	gint yticks;
	gdouble xmin;
	gdouble xmax;
	gdouble ymin;
	gdouble ymax;
	/* NB: all sets are required to be <= size */
	gint size;
	GSList *set_list;/*in case of 1D graph, set_list={x[]} for other set_list={g_data_x,g_data_y[]}*/
	/* peak selection */
	gint select;
	gdouble select_2;/*peak selections for 2D graph*/
	gchar *select_label;
	/* for graph_controls */
	gchar *title;
	gchar *sub_title;
	gchar *x_title;
	gchar *y_title;
};
gpointer graph_new(const gchar *, struct model_pak *);
void graph_reset_data(struct graph_pak *graph);
void graph_reset(struct graph_pak *graph);
void graph_free_list(struct model_pak *);
void graph_free(gpointer, struct model_pak *);
void graph_set_color(gint size,gchar *color,gpointer data);
void graph_add_data(gint, gdouble *, gdouble, gdouble, gpointer);
void graph_add_borned_data(gint size,gdouble *x,gdouble x_min,gdouble x_max,gdouble y_min,gdouble y_max,gint type,gpointer data);
void graph_frequency_select(gint x, gint y, struct model_pak *model);
void graph_uspex_select(gint x, gint y, struct model_pak *model);
void graph_uspex_2d_select(gint x, gint y, struct model_pak *model);
void graph_set_grafted(gint, gpointer);
void graph_set_xticks(gint, gint, gpointer);
void graph_set_yticks(gint, gint, gpointer);
void graph_set_wavelength(gdouble, gpointer);
void graph_set_select(gdouble, gchar *, gpointer);
void graph_write(gchar *, gpointer);
void graph_read(gchar *);
gchar *graph_treename(gpointer);

gint graph_grafted(gpointer);
gint graph_xlabel(gpointer);
gint graph_ylabel(gpointer);
gdouble graph_xmin(gpointer);
gdouble graph_xmax(gpointer);
gdouble graph_ymin(gpointer);
gdouble graph_ymax(gpointer);
/*for diffraction*/
gdouble graph_wavelength(gpointer);
GSList *diff_get_ranked_faces(gdouble min, struct model_pak *model);
/* dat_graph interface */
void dat_graph_toggle_xaxis(gpointer pgraph);
void dat_graph_toggle_yaxis(gpointer pgraph);
void dat_graph_set_title(const gchar *title,gpointer pgraph);
void dat_graph_set_sub_title(const gchar *title,gpointer pgraph);
void dat_graph_set_x_title(const gchar *x_title,gpointer pgraph);
void dat_graph_set_y_title(const gchar *y_title,gpointer pgraph);
void dat_graph_set_x(g_data_x dx,gpointer pgraph);
void dat_graph_add_y(g_data_y dy,gpointer pgraph);
void dat_graph_set_limits(gdouble x_min,gdouble x_max,gdouble y_min,gdouble y_max,gpointer pgraph);
void dat_graph_set_type(graph_type type,gpointer pgraph);
void dat_graph_select(gint x, gint y, struct model_pak *model);
/* control interface */
void gui_graph_controls(void);


/* TODO - relocate */
gint anim_next_frame(struct model_pak *);

