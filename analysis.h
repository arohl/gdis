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

#include <glib.h>
#include "task.h"

/* position storage */
struct gd3_pak
{
gdouble x[3];
};

/* matrix storage */
struct gd9_pak
{
gdouble m[9];
};

struct analysis_pak
{
/* analysis setup */
gint num_bins;
gdouble start;
gdouble stop;
gdouble step;
gchar *atom1;
gchar *atom2;
gint rdf_normalize;
/* model setup */
gint num_atoms;
gint num_frames;
gdouble time_start;
gdouble time_stop;
/* packed frame data */
struct gd9_pak *latmat;
struct gd3_pak *position;
struct gd3_pak *velocity;
gdouble *time;
gdouble *ke;
gdouble *pe;
gdouble *temp;
};

/* prototypes */
gint analysis_init(struct model_pak *);
void analysis_free(struct analysis_pak *);
gpointer analysis_new(void);
gpointer analysis_dup(gpointer);
void analysis_show(struct model_pak *);

void analysis_plot_rdf(struct analysis_pak *, struct task_pak *);
void analysis_plot_meas(struct analysis_pak *, struct task_pak *);
void analysis_plot_vacf(struct analysis_pak *, struct task_pak *);
void analysis_plot_temp(struct analysis_pak *, struct task_pak *);
void analysis_plot_ke(struct analysis_pak *, struct task_pak *);
void analysis_plot_pe(struct analysis_pak *, struct task_pak *);

