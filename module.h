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

/* public */
void module_setup(void);
void module_free(void);
void module_label_set(gchar *, gpointer);
void module_symbol_add(gchar *s, gpointer);
void module_resident_set(gint, gpointer);

gchar *module_label(gpointer);
GSList *module_symbols(gpointer);
gint module_invoke(gchar *, gpointer, gpointer);

/* Hirshfeld surface prototypes */
void hfs_init(void);
gdouble hfs_calc_wf(gdouble *, struct model_pak *, GSList *);
GSList *hfs_calc_env(GSList *, struct model_pak *);
int hfs_calc_normals(GSList *, struct model_pak *, GSList *);
int hfs_calulation_limits(GSList *selection, struct model_pak *, gdouble *min, gdouble *max);
