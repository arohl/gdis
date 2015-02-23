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

/* types */
enum
{
FF_HARMONIC, FF_MORSE, FF_BUCKINGHAM, FF_LENNARD,
FF_3B_HARMONIC, FF_DIHEDRAL, FF_IMPROPER, FF_DIHEDRAL_RB
};

/* units */
enum
{
FF_UNKNOWN,
FF_ANG, FF_AU, FF_DEG, FF_RAD,
FF_EV, FF_KJ, FF_KCAL
};

/* prototypes */

gpointer ff_type_new(gint);
gpointer ff_dup(gpointer);

gint ff_type_get(gpointer);
gint ff_data_size_get(gpointer);
gint ff_data_add(gpointer, gchar *);
gdouble *ff_data_get(gpointer);
gpointer ff_search(struct core_pak **, gint, GSList *);
GSList *ff_filter_list(gint, GSList *);

void ff_dump_all(GSList *);
void ff_dump_type(gint, GSList *);

GSList *ff_gulp_parse(const gchar *);
gpointer ff_gulp_new(const gchar *);
gchar *ff_gulp_string(gpointer);

