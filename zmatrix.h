/*
Copyright (C) 2000 by Sean David Fleming

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

void zmat_debug(gpointer);
void zmat_mol_print(FILE *, gpointer);
void zmat_var_print(FILE *, gpointer);
void zmat_const_print(FILE *, gpointer);
void zmat_coord_print(FILE *, GSList *, gpointer);

gpointer zmat_new(void);

void zmat_build(void);

void zmat_free(gpointer);

gpointer zval_new(void);

gint zmat_entries_get(gpointer);

void zmat_cartesian_set(gpointer);
void zmat_fractional_set(gpointer);

void zmat_distance_units_set(gpointer, gint);
void zmat_angle_units_set(gpointer, gint);

const gchar *zmat_distance_units_get(gpointer);
const gchar *zmat_angle_units_get(gpointer);

void zmat_core_add(const gchar *, gpointer);
void zmat_var_add(const gchar *, gpointer);
void zmat_const_add(const gchar *, gpointer);

void zmat_type(gpointer, GSList *);

void zmat_process(gpointer, struct model_pak *);

void gui_zmat_dialog(void);

/* a SIESTA special... */
struct species_pak
{
gint number;   /* can be -ve -> ghost */
gchar *label;  /* same as core->label */
};

