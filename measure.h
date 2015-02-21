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

/* constants */

enum {MEASURE_BOND, MEASURE_INTRA, MEASURE_INTER, MEASURE_DISTANCE,
      MEASURE_ANGLE, MEASURE_TORSION};

/* prototypes */

gpointer measure_bond_test(struct core_pak **,
                           gdouble,
                           gdouble,
                           struct model_pak *);
gpointer measure_distance_test(gint,
                               struct core_pak **,
                               gdouble,
                               gdouble,
                               struct model_pak *);
gpointer measure_angle_test(struct core_pak **,
                            gdouble,
                            gdouble,
                            struct model_pak *);
gpointer measure_torsion_test(struct core_pak **,
                              gdouble,
                              gdouble,
                              struct model_pak *);

void measure_bond_search(const gchar **, gdouble, gdouble, struct model_pak *);
void measure_distance_search(const gchar **, gint, gdouble, gdouble, struct model_pak *);
void measure_bangle_search(const gchar **, gdouble, gdouble, struct model_pak *);
void measure_angle_search(const gchar **, gdouble *, struct model_pak *);

gint measure_type_get(gpointer);
gchar *measure_value_get(gpointer);
GSList *measure_cores_get(gpointer);
void measure_colour_get(gdouble *, gpointer);
void measure_coord_get(gdouble *, gint, gpointer, struct model_pak *);
gboolean measure_has_core(struct core_pak *, gpointer);

gchar *measure_type_label_create(gpointer);
gchar *measure_constituents_create(gpointer);

void measure_colour_set(gdouble, gdouble, gdouble, gpointer);

void measure_free(gpointer, struct model_pak *);
void measure_free_all(struct model_pak *);
void measure_dump_all(struct model_pak *);
void measure_select_all(void);

void meas_prune_model(struct model_pak *);
void meas_graft_model(struct model_pak *);

gdouble measure_update_single(gpointer, struct model_pak *);
void measure_update_global(struct model_pak *);

gdouble measure_distance(gdouble *, gdouble *);
gdouble measure_angle(gdouble *, gdouble *, gdouble *);
gdouble measure_dihedral(gdouble *, gdouble *, gdouble *, gdouble *);

/* bonded case of dihedral */
gdouble measure_torsion(struct core_pak **);

