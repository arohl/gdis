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

/* header for plotting functions */

#include "task.h"

typedef enum {
        PLOT_NONE=0,
        PLOT_ENERGY=1,
        PLOT_FORCE=2,
        PLOT_VOLUME=4,
        PLOT_PRESSURE=8,
        PLOT_BAND=16,
        PLOT_DOS=32,
        PLOT_BANDOS=48,
        PLOT_FREQUENCY=64,
	PLOT_RAMAN=128
} plot_type;
typedef struct {
        gint size;
        gdouble *data;
        gdouble xmin;
        gdouble xmax;
        gdouble ymin;
        gdouble ymax;
} preload_plot_data;
struct plot_pak {
        /* store things differently than graph */
	struct model_pak *model;/*connection to model*/
        plot_type type;/*What we want to plot*/
        gint plot_mask;/*What we can plot*/
	gint plot_sel;/*what user selected in all pages*/
        gpointer graph;/*connection to graph*/
        preload_plot_data energy;
        preload_plot_data force;
        preload_plot_data volume;
        preload_plot_data pressure;
        preload_plot_data band;
        preload_plot_data dos;
        preload_plot_data frequency;
	preload_plot_data raman;
        /* plot parameter */
	gint ndos;
	gint nbands;
	gint nfreq;
	gint nraman;
        gdouble xmin;
        gdouble xmax;
        gdouble ymin;
        gdouble ymax;
        gdouble xtics;
        gdouble ytics;
        gboolean auto_x;
        gboolean auto_y;
        gboolean data_changed;
};

/* prototypes */
void plot_prepare_data(struct plot_pak *plot);
void plot_load_data(struct plot_pak *plot,struct task_pak *task);

void plot_draw_graph(struct plot_pak *plot,struct task_pak *task);
void plot_show_graph(struct plot_pak *plot);


